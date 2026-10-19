package vcf

import (
	"errors"
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/inful/htsfile/csi"
	"github.com/inful/htsfile/internal/bgzf"
)

// CSIReader is the CSI (Coordinate Sorting Index) sibling
// of [TabixReader]. It supports random access by
// genomic region using a CSI index (.csi) alongside the
// BGZF-compressed VCF file.
//
// A CSIReader is not safe for concurrent use.
type CSIReader struct {
	seek     *bgzf.SeekReader
	idx      *csi.Index
	hdr      *Header
	path     string
	chromIdx map[string]int
}

// OpenCSI opens a CSI-indexed VCF file at path. The
// sibling .csi file at "<path>.csi" must exist. The VCF
// file is expected to be BGZF-compressed.
func OpenCSI(path string) (*CSIReader, error) {
	f, err := os.Open(path) //nolint:gosec // path is caller-supplied
	if err != nil {
		return nil, fmt.Errorf("vcf: open %s: %w", path, err)
	}
	seek, err := bgzf.NewSeekReader(f)
	if err != nil {
		_ = f.Close()
		return nil, fmt.Errorf("vcf: open bgzf %s: %w", path, err)
	}
	csiPath := path + ".csi"
	csiFile, err := os.Open(csiPath) //nolint:gosec // sibling of caller-supplied
	if err != nil {
		_ = seek.Close()
		return nil, fmt.Errorf("vcf: open index %s: %w", csiPath, err)
	}
	idx, err := csi.NewReader(csiFile).Read()
	_ = csiFile.Close()
	if err != nil {
		_ = seek.Close()
		return nil, fmt.Errorf("vcf: parse index %s: %w", csiPath, err)
	}
	return &CSIReader{seek: seek, idx: idx, path: path, chromIdx: make(map[string]int)}, nil
}

// Close releases the underlying BGZF reader.
func (r *CSIReader) Close() error {
	if r.seek == nil {
		return nil
	}
	err := r.seek.Close()
	r.seek = nil
	return err
}

// Header returns the VCF header. The first call seeks to
// the start of the BGZF stream and parses the header.
// Subsequent calls return the cached value.
func (r *CSIReader) Header() (*Header, error) {
	if r.hdr != nil {
		return r.hdr, nil
	}
	if err := r.seek.Seek(0); err != nil {
		return nil, fmt.Errorf("vcf: seek to header: %w", err)
	}
	dec := NewDecoder(r.seek)
	hdr, err := dec.Header()
	if err != nil {
		return nil, fmt.Errorf("vcf: parse header: %w", err)
	}
	r.hdr = hdr
	return hdr, nil
}

// refByHeaderIndex returns the CSI reference whose
// position in the index matches the given header contig
// index. CSI references are ordered the same as the
// header contigs.
func (r *CSIReader) refByHeaderIndex(idx int) *csi.Reference {
	if idx < 0 || idx >= len(r.idx.References) {
		return nil
	}
	return r.idx.References[idx]
}

// Query returns the data records that overlap the
// half-open 1-based region [start, end) on chrom. The
// returned records are in file order.
//
// Unknown chroms return an empty slice (no error). If
// the region is malformed (start < 1 or end <= start),
// Query returns an error.
func (r *CSIReader) Query(chrom string, start, end int64) ([]*Record, error) {
	if start < 1 {
		return nil, fmt.Errorf("vcf: query start %d is 0-based; use 1-based", start)
	}
	if end <= start {
		return nil, fmt.Errorf("vcf: query end %d must be > start %d", end, start)
	}
	hdr, err := r.Header()
	if err != nil {
		return nil, err
	}
	idx := -1
	for i, c := range hdr.Contigs {
		r.chromIdx[c.ID] = i
		if c.ID == chrom {
			idx = i
			break
		}
	}
	if idx < 0 {
		return nil, nil
	}
	ref := r.refByHeaderIndex(idx)
	if ref == nil {
		return nil, nil
	}
	refBins := csi.BinFromRange(start, end)
	chunks := ref.LookupChunks(refBins)
	if len(chunks) == 0 {
		return nil, nil
	}
	type pending struct {
		end  uint64
		virt uint64
	}
	pendingChunks := make([]pending, len(chunks))
	for i, c := range chunks {
		pendingChunks[i] = pending{
			end:  c.End,
			virt: c.Beg,
		}
	}
	sort.Slice(pendingChunks, func(i, j int) bool {
		return pendingChunks[i].virt < pendingChunks[j].virt
	})
	var recs []*Record
	dedup := make(map[uint64]struct{})
	for _, c := range pendingChunks {
		chunkRecs, err := r.decodeChunk(c.virt, c.end)
		if err != nil {
			return nil, err
		}
		for _, rec := range chunkRecs {
			if rec.Chrom != chrom {
				continue
			}
			if rec.Pos < start || rec.Pos >= end {
				continue
			}
			key := dedupKey(r.chromIdx, rec.Chrom, rec.Pos)
			if _, ok := dedup[key]; ok {
				continue
			}
			dedup[key] = struct{}{}
			recs = append(recs, rec)
		}
	}
	return recs, nil
}

// decodeChunk seeks to virtualOff and decodes VCF records
// until the BGZF virtual offset reaches chunkEndOff (or
// EOF, whichever comes first). Records are filtered by
// region by the caller.
func (r *CSIReader) decodeChunk(virtualOff, chunkEndOff uint64) ([]*Record, error) {
	if err := r.seek.Seek(virtualOff); err != nil {
		return nil, fmt.Errorf("vcf: seek chunk: %w", err)
	}
	dec := NewDecoder(r.seek)
	var recs []*Record
	for {
		off := r.seek.VirtualOffset()
		if chunkEndOff > 0 && off >= chunkEndOff {
			break
		}
		rec, err := dec.Record()
		if err != nil {
			if errors.Is(err, io.EOF) {
				break
			}
			return nil, fmt.Errorf("vcf: read record: %w", err)
		}
		if rec == nil {
			break
		}
		recs = append(recs, rec)
	}
	return recs, nil
}