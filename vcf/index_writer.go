package vcf

import (
	"errors"
	"fmt"
	"io"
	"os"

	"github.com/inful/htsfile/csi"
	"github.com/inful/htsfile/internal/bgzf"
	"github.com/inful/htsfile/tabix"
)

// IndexOpts controls how BuildIndex emits an index.
//
// Format selects the file format (.tbi or .csi).
// MinShift and Depth only matter when Format is FormatCSI:
// they configure the CSI bin geometry. For tabix the
// binning uses the canonical 16 KiB / 5-level scheme.
type IndexOpts struct {
	Format   IndexFormat
	MinShift int32
	Depth    int32
}

// IndexFormat selects the on-disk index format.
type IndexFormat int

const (
	// FormatTabix emits a .tbi index (samtools tabix).
	// This is the default if Format is unspecified.
	FormatTabix IndexFormat = iota
	// FormatCSI emits a .csi index with the given
	// min_shift and depth.
	FormatCSI
)

// BuildIndex reads a BGZF-compressed VCF file at path
// and writes a tabix .tbi or CSI .csi index next to it.
// The VCF must be coordinate-sorted (CHROM, then POS) per
// the tabix / CSI specification.
//
// Algorithm: the vcf.Decoder buffers blocks ahead, so
// per-record offsets inside a block are unreliable. We
// instead use the BGZF block boundary as the natural
// chunk boundary — every block produces at most one
// chunk per bin. This matches the per-block chunk
// layout that samtools `tabix -p vcf` and
// `bcftools index` emit by default.
//
//   - When the reader crosses into a new BGZF block,
//     close any open chunks with End = start of the new
//     block.
//   - On each record, for every bin in
//     reg2bins([POS-1, POS)): if the bin already has an
//     open chunk in this block, reuse it; otherwise
//     open a new chunk with Beg = this block's start.
//   - The linear (interval) index records the BGZF
//     offset of the first record seen in each
//     (2^min_shift)-bp window.
//
// At EOF any chunks still open get End = 0xFFFF...
// (the canonical "read to EOF" sentinel that tabix and
// CSI readers both treat as unbounded).
//
// Caveat: a single BGZF block can hold records from
// multiple CHROMs. The resulting index chunks cover all
// records in the block regardless of CHROM; the caller
// (OpenTabix / OpenCSI) must filter by CHROM at read
// time. TabixReader.Query does this in its inner loop.
func BuildIndex(path string, opts IndexOpts) error {
	f, err := os.Open(path) //nolint:gosec // caller-supplied path
	if err != nil {
		return fmt.Errorf("vcf: open %s: %w", path, err)
	}
	defer func() { _ = f.Close() }()

	seek, err := bgzf.NewSeekReader(f)
	if err != nil {
		return fmt.Errorf("vcf: open bgzf %s: %w", path, err)
	}
	defer func() { _ = seek.Close() }()

	// Parse the header first; we need its contigs to
	// initialise per-reference index state.
	if err := seek.Seek(0); err != nil {
		return fmt.Errorf("vcf: seek header: %w", err)
	}
	dec := NewDecoder(seek)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("vcf: parse header: %w", err)
	}
	chromOrder := make(map[string]int, len(hdr.Contigs))
	for i, c := range hdr.Contigs {
		chromOrder[c.ID] = i
	}

	b := &indexBuilder{
		chromOrder: chromOrder,
		refs:       make([]indexRef, len(hdr.Contigs)),
		openChunks: make(map[[2]int]*ChunkBufIndex),
	}
	switch opts.Format {
	case FormatCSI:
		b.minShift = opts.MinShift
		b.depth = opts.Depth
	default:
		// Tabix: 16 KiB / 5-level binning.
		b.minShift = 14
		b.depth = 5
	}
	if b.minShift < 0 || b.minShift > 30 {
		return fmt.Errorf("vcf: invalid min_shift %d (must be 0..30)", b.minShift)
	}
	if b.depth < 1 || b.depth > 10 {
		return fmt.Errorf("vcf: invalid depth %d (must be 1..10)", b.depth)
	}
	b.windowSize = int64(1) << uint(b.minShift)

	// eofVirtOff is the BGZF virtual offset of the
	// 28-byte empty EOF marker at the end of the
	// file. Set as the End of every chunk in the last
	// block — what samtools / bcftools emit by default
	// — and what makes the chunk-End stop actually
	// fire for single-block files.
	stat, err := f.Stat()
	if err != nil {
		return fmt.Errorf("vcf: stat %s: %w", path, err)
	}
	const eofMarkerSize = 28
	fileSize := stat.Size()
	if fileSize < eofMarkerSize {
		return fmt.Errorf("vcf: %s is too small to be BGZF", path)
	}
	b.eofVirtOff = uint64(fileSize-eofMarkerSize) << 16

	prevBlockOff := uint64(0)
	curBlockStartOff := uint64(0)
	for {
		off := seek.VirtualOffset()
		curBlockOff := off >> 16
		if curBlockOff != prevBlockOff {
			// New block: close any open chunks
			// with End = start of the new block.
			for _, oc := range b.openChunks {
				b.refs[oc.refIdx].chunks[oc.bin][oc.idx].End = tabix.VirtualOffset(off)
			}
			b.openChunks = make(map[[2]int]*ChunkBufIndex)
			curBlockStartOff = off
			prevBlockOff = curBlockOff
		}
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("vcf: read record: %w", err)
		}
		if rec == nil {
			break
		}
		refIdx, ok := chromOrder[rec.Chrom]
		if !ok {
			// Unknown contig: skip. A contig in
			// the data but not the header is a
			// VCF spec violation; dropping it
			// is safer than mis-numbering the
			// bin index.
			continue
		}
		bins := csi.BinFromRange(rec.Pos-1, rec.Pos)
		if len(bins) == 0 {
			continue
		}
		ref := &b.refs[refIdx]
		if ref.chunks == nil {
			ref.chunks = make(map[uint32][]tabix.Chunk)
		}
		for _, bin := range bins {
			key := [2]int{refIdx, int(bin)}
			if _, ok := b.openChunks[key]; !ok {
				// First record in this block
				// for this bin: open a chunk.
				idx := len(ref.chunks[bin])
				ref.chunks[bin] = append(ref.chunks[bin], tabix.Chunk{
					Beg: tabix.VirtualOffset(curBlockStartOff),
				})
				b.openChunks[key] = &ChunkBufIndex{
					refIdx: refIdx,
					bin:    bin,
					idx:    idx,
				}
			}
			// Otherwise: the existing chunk for
			// this bin in this block is reused.
			// Multiple records in the same block
			// share one chunk per bin.
		}
		// Linear (interval) index: first record per
		// (2^min_shift)-bp window sets the offset.
		winIdx := (rec.Pos - 1) >> b.minShift
		for int64(len(ref.intervals)) <= winIdx {
			ref.intervals = append(ref.intervals, 0)
		}
		if ref.intervals[winIdx] == 0 {
			ref.intervals[winIdx] = off
		}
	}
	// Close any chunks still open at EOF with the
	// BGZF virtual offset of the EOF marker.
	for _, oc := range b.openChunks {
		b.refs[oc.refIdx].chunks[oc.bin][oc.idx].End = tabix.VirtualOffset(b.eofVirtOff)
	}

	switch opts.Format {
	case FormatCSI:
		return b.writeCSI(path)
	default:
		return b.writeTabix(path)
	}
}

// ChunkBufIndex tracks an open chunk in the index
// builder. Used to set the chunk's End when a block
// boundary is crossed.
type ChunkBufIndex struct {
	refIdx int
	bin    uint32
	idx    int
}

// indexRef is the per-reference (per-contig) state
// during a BuildIndex run.
type indexRef struct {
	chunks    map[uint32][]tabix.Chunk
	intervals []uint64
}

// indexBuilder accumulates index state across the
// whole file.
type indexBuilder struct {
	chromOrder map[string]int
	refs       []indexRef
	openChunks map[[2]int]*ChunkBufIndex
	minShift   int32
	depth      int32
	windowSize int64
	// eofVirtOff is the BGZF virtual offset of the
	// 28-byte EOF marker at the end of the file.
	// Set once at construction from f.Stat().
	eofVirtOff uint64
}

// writeTabix serialises b's state to a .tbi file
// next to the VCF. The .tbi file's references mirror
// the input header's contig order.
func (b *indexBuilder) writeTabix(path string) error {
	refs := make([]*tabix.Reference, 0, len(b.refs))
	for i := range b.refs {
		r := &b.refs[i]
		if r.chunks == nil {
			r.chunks = map[uint32][]tabix.Chunk{}
		}
		intervals := make([]tabix.VirtualOffset, len(r.intervals))
		for j, iv := range r.intervals {
			intervals[j] = tabix.VirtualOffset(iv)
		}
		refs = append(refs, &tabix.Reference{
			Bins:      r.chunks,
			Intervals: intervals,
		})
	}
	tbiPath := path + ".tbi"
	tbiFile, err := os.Create(tbiPath) //nolint:gosec // sibling of caller-supplied VCF
	if err != nil {
		return fmt.Errorf("vcf: create %s: %w", tbiPath, err)
	}
	defer func() { _ = tbiFile.Close() }()
	if err := tabix.NewWriter(tbiFile).Write(&tabix.Index{References: refs}); err != nil {
		return fmt.Errorf("vcf: write %s: %w", tbiPath, err)
	}
	return nil
}

// writeCSI serialises b's state to a .csi file next
// to the VCF.
func (b *indexBuilder) writeCSI(path string) error {
	refs := make([]*csi.Reference, 0, len(b.refs))
	for i := range b.refs {
		r := &b.refs[i]
		if r.chunks == nil {
			r.chunks = map[uint32][]tabix.Chunk{}
		}
		chunks := make(map[uint32][]csi.Chunk, len(r.chunks))
		for bin, cs := range r.chunks {
			out := make([]csi.Chunk, len(cs))
			for j, c := range cs {
				out[j] = csi.Chunk{Beg: uint64(c.Beg), End: uint64(c.End)} //nolint:gosec // csi.Chunk fields are uint64 by definition
			}
			chunks[bin] = out
		}
		intervals := make([]csi.VirtualOffset, len(r.intervals))
	copy(intervals, r.intervals)
		refs = append(refs, &csi.Reference{
			Bins:      chunks,
			Intervals: intervals,
		})
	}
	csiPath := path + ".csi"
	csiFile, err := os.Create(csiPath) //nolint:gosec // sibling of caller-supplied VCF
	if err != nil {
		return fmt.Errorf("vcf: create %s: %w", csiPath, err)
	}
	defer func() { _ = csiFile.Close() }()
	if err := csi.NewWriter(csiFile).Write(&csi.Index{
		MinShift:   b.minShift,
		Depth:      b.depth,
		References: refs,
	}); err != nil {
		return fmt.Errorf("vcf: write %s: %w", csiPath, err)
	}
	return nil
}