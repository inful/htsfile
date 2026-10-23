package vcf

import (
	"bytes"
	"encoding/binary"
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
	// Exact reads: the index needs the virtual offset of every
	// record, and the scanner's read-ahead would put the
	// reader's position past the start of the next record.
	dec.SetExactReads(true)
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
	}
	// The reference names go in the index: the .tbi carries
	// them in its name table and the .csi in its aux field.
	for i, c := range hdr.Contigs {
		b.refs[i].name = c.ID
		b.names = append(b.names, c.ID)
	}
	b.format = opts.Format
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

	// htslib's chunking convention (hts.c hts_idx_push and
	// insert_to_l), which the tabix format depends on:
	//
	//   - a chunk runs from the start of the first record in
	//     its bin to the start of the first record in the next
	//     bin, so a reader that starts at chunk.beg and scans
	//     forward sees every record the bin covers;
	//   - a chunk is closed when the bin changes;
	//   - the linear index entry for a window is the start of
	//     the first record at or after that window, and a
	//     window with no record of its own inherits the entry
	//     of the window after it.
	//
	// saveOff is the start of the chunk being built, saveBin
	// its bin, and the unset marker for the linear index is
	// all-ones, as htslib's is.
	const unsetOffset = ^uint64(0)
	var saveOff uint64
	saveBin := uint32(0xffffffff)
	saveRef := -1
	// offBeg is the start of the first record and nRecords the
	// number of records seen; together they form the tabix
	// meta bin (htslib's META_BIN).
	var offBeg uint64
	nRecords := 0
	for {
		off := seek.VirtualOffset()
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
		// The bin a record is filed under: the finest bin
		// that covers its position. An index holds only
		// those leaf bins - htslib's tabix output has one
		// bin per record position group, not the whole
		// ancestor chain - and a CSI has one more level
		// than a tabix index, so its leaf bins are numbered
		// from 37449 rather than 4681.
		levels := b.depth
		if b.format == FormatCSI {
			levels = b.depth + 1
		}
		bins := csi.BinFromRangeWith(rec.Pos-1, rec.Pos, b.minShift, levels)
		if len(bins) == 0 {
			continue
		}
		ref := &b.refs[refIdx]
		if ref.chunks == nil {
			ref.chunks = make(map[uint32][]tabix.Chunk)
		}
		if nRecords == 0 {
			offBeg = off
		}
		nRecords++
		bin := bins[len(bins)-1]
		if bin != saveBin || refIdx != saveRef {
			// The bin changed: close the previous chunk
			// and start a new one. A chunk that spans
			// several bins is not an error, but htslib
			// closes on a bin change and so do we, so
			// the two indexes are comparable.
			if saveBin != 0xffffffff && saveRef >= 0 {
				prev := &b.refs[saveRef]
				prev.chunks[saveBin] = append(prev.chunks[saveBin], tabix.Chunk{
					Beg: tabix.VirtualOffset(saveOff),
					End: tabix.VirtualOffset(off),
				})
			}
			saveOff = off
			saveBin = bin
			saveRef = refIdx
		}
		// Linear (interval) index: the window's entry is
		// the start of the first record in it. A window
		// with no record of its own inherits the entry of
		// the window after it (filled in below).
		winIdx := (rec.Pos - 1) >> b.minShift
		for int64(len(ref.intervals)) <= winIdx {
			ref.intervals = append(ref.intervals, unsetOffset)
		}
		if ref.intervals[winIdx] == unsetOffset {
			ref.intervals[winIdx] = off
		}
	}
	// Close the last chunk at the EOF marker.
	if saveBin != 0xffffffff && saveRef >= 0 {
		last := &b.refs[saveRef]
		last.chunks[saveBin] = append(last.chunks[saveBin], tabix.Chunk{
			Beg: tabix.VirtualOffset(saveOff),
			End: tabix.VirtualOffset(b.eofVirtOff),
		})
	}
	// The meta bin holds the file's byte range and the mapped /
	// unmapped record counts (htslib's META_BIN). A CSI has one
	// more level than a tabix index, so its meta bin is numbered
	// from 299594 rather than 37450.
	if saveRef >= 0 {
		levels := b.depth
		if b.format == FormatCSI {
			levels = b.depth + 1
		}
		meta := csi.BinLimit(levels) + 1
		ref := &b.refs[saveRef]
		if ref.chunks == nil {
			ref.chunks = make(map[uint32][]tabix.Chunk)
		}
		ref.chunks[uint32(meta)] = []tabix.Chunk{ //nolint:gosec // the meta bin number is far below 2^32
			{Beg: tabix.VirtualOffset(offBeg), End: tabix.VirtualOffset(b.eofVirtOff)},
			{Beg: tabix.VirtualOffset(nRecords), End: 0},
		}
	}
	// Fill the gaps in the linear index: a window with no
	// record of its own inherits the entry of the window
	// after it, so a reader that starts there and scans
	// forward finds the records that overlap the window.
	for i := range b.refs {
		iv := b.refs[i].intervals
		for j := len(iv) - 2; j >= 0; j-- {
			if iv[j] == unsetOffset {
				iv[j] = iv[j+1]
			}
		}
	}

	switch opts.Format {
	case FormatCSI:
		return b.writeCSI(path)
	default:
		return b.writeTabix(path)
	}
}

// indexRef is the per-reference (per-contig) state
// during a BuildIndex run.
type indexRef struct {
	chunks    map[uint32][]tabix.Chunk
	intervals []uint64
	// name is the reference sequence name. The .tbi carries
	// the names in its name table, and the .csi carries them
	// in its aux field; without them a region query cannot
	// resolve a reference.
	name string
}

// indexBuilder accumulates index state across the
// whole file.
type indexBuilder struct {
	chromOrder map[string]int
	refs       []indexRef
	minShift   int32
	depth      int32
	windowSize int64
	// names are the reference names in index order.
	names []string
	// format selects the tabix or CSI layout.
	format IndexFormat
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
			Name:      r.name,
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
	idx := &tabix.Index{
		// The VCF preset, as htslib's tabix writes it: the
		// sequence is column 1, the position is column 2,
		// and the end column is inferred.
		Format:     2, // VCF
		ColSeq:     1,
		ColBeg:     2,
		ColEnd:     0,
		MetaChar:   '#',
		SkipLines:  0,
		Names:      b.names,
		References: refs,
	}
	if err := tabix.NewWriter(tbiFile).Write(idx); err != nil {
		return fmt.Errorf("vcf: write %s: %w", tbiPath, err)
	}
	return nil
}

// writeCSI serialises b's state to a .csi file next
// to the VCF.
func (b *indexBuilder) writeCSI(path string) error {
	// The meta bin, whose loffset htslib leaves at zero.
	metaBin := uint32(csi.BinLimit(b.depth+1) + 1) //nolint:gosec // the meta bin number is far below 2^32
	refs := make([]*csi.Reference, 0, len(b.refs))
	for i := range b.refs {
		r := &b.refs[i]
		if r.chunks == nil {
			r.chunks = map[uint32][]tabix.Chunk{}
		}
		chunks := make(map[uint32][]csi.Chunk, len(r.chunks))
		loffsets := make(map[uint32]csi.VirtualOffset, len(r.chunks))
		for bin, cs := range r.chunks {
			out := make([]csi.Chunk, len(cs))
			for j, c := range cs {
				out[j] = csi.Chunk{Beg: uint64(c.Beg), End: uint64(c.End)} //nolint:gosec // csi.Chunk fields are uint64 by definition
			}
			chunks[bin] = out
			// The loffset is the virtual offset of the first
			// record in the bin. htslib leaves the meta
			// bin's loffset at zero.
			if len(out) == 0 {
				continue
			}
			if bin == metaBin {
				loffsets[bin] = 0
				continue
			}
			loffsets[bin] = out[0].Beg
		}
		refs = append(refs, &csi.Reference{
			Bins:     chunks,
			Loffsets: loffsets,
		})
	}
	csiPath := path + ".csi"
	csiFile, err := os.Create(csiPath) //nolint:gosec // sibling of caller-supplied VCF
	if err != nil {
		return fmt.Errorf("vcf: create %s: %w", csiPath, err)
	}
	defer func() { _ = csiFile.Close() }()
	if err := csi.NewWriter(csiFile).Write(&csi.Index{
		MinShift: b.minShift,
		// A CSI records the number of bin *levels*, which is
		// one more than the tabix depth. htslib writes 6 for
		// the tabix depth of 5.
		Depth: b.depth + 1,
		// The aux field carries the tabix preset and the
		// reference names, the same way htslib does.
		Aux:        b.csiAux(),
		References: refs,
	}); err != nil {
		return fmt.Errorf("vcf: write %s: %w", csiPath, err)
	}
	return nil
}

// csiAux returns the CSI aux field: the tabix preset fields
// followed by the reference names, exactly as htslib writes
// them (format, col_seq, col_beg, col_end, meta, skip,
// l_nm, then the NUL-separated names).
func (b *indexBuilder) csiAux() []byte {
	var buf bytes.Buffer
	put := func(v int32) {
		var x [4]byte
		binary.LittleEndian.PutUint32(x[:], uint32(v)) //nolint:gosec // the preset fields are small positive numbers
		buf.Write(x[:])
	}
	put(2) // VCF
	put(1) // col_seq
	put(2) // col_beg
	put(0) // col_end
	put('#')
	put(0) // skip
	// The name table is each name followed by a NUL,
	// including the last one.
	var nameBuf []byte
	for _, n := range b.names {
		nameBuf = append(nameBuf, n...)
		nameBuf = append(nameBuf, 0)
	}
	put(int32(len(nameBuf))) //nolint:gosec // the name table is far below 2^31 bytes
	buf.Write(nameBuf)
	return buf.Bytes()
}
