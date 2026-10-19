package vcf

import (
	"errors"
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/inful/htsfile/internal/bgzf"
	"github.com/inful/htsfile/tabix"
)

// TabixReader is a read-only handle on a tabix-indexed
// VCF file. It supports random access by genomic region:
// callers compute the bins for the region via the
// tabix index, look up the chunks, and decode only those
// overlapping the region.
//
// A TabixReader is not safe for concurrent use.
type TabixReader struct {
	seek *bgzf.SeekReader
	idx  *tabix.Index
	hdr  *Header
	path string
}

// OpenTabix opens a tabix-indexed VCF file at path. The
// sibling .tbi file at "<path>.tbi" must exist. The VCF
// file is expected to be BGZF-compressed (as is
// conventional for indexed VCF).
func OpenTabix(path string) (*TabixReader, error) {
	f, err := os.Open(path) //nolint:gosec // path is caller-supplied
	if err != nil {
		return nil, fmt.Errorf("vcf: open %s: %w", path, err)
	}
	seek, err := bgzf.NewSeekReader(f)
	if err != nil {
		_ = f.Close()
		return nil, fmt.Errorf("vcf: open bgzf %s: %w", path, err)
	}
	tbiPath := path + ".tbi"
	tbiFile, err := os.Open(tbiPath) //nolint:gosec // sibling of caller-supplied
	if err != nil {
		_ = seek.Close()
		return nil, fmt.Errorf("vcf: open index %s: %w", tbiPath, err)
	}
	idx, err := tabix.NewReader(tbiFile).Read()
	_ = tbiFile.Close()
	if err != nil {
		_ = seek.Close()
		return nil, fmt.Errorf("vcf: parse index %s: %w", tbiPath, err)
	}
	return &TabixReader{seek: seek, idx: idx, path: path}, nil
}

// Close releases the underlying BGZF reader.
func (r *TabixReader) Close() error {
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
func (r *TabixReader) Header() (*Header, error) {
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

// refByHeaderIndex returns the index reference whose
// position in the tabix index matches the given header
// contig index. The tabix index stores references in the
// same order as the header, so we look up by position.
func (r *TabixReader) refByHeaderIndex(idx int) *tabix.Reference {
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
func (r *TabixReader) Query(chrom string, start, end int64) ([]*Record, error) {
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
	// Find the header contig for this chrom. tabix
	// references are ordered the same as the header
	// contigs, so we match by index.
	idx := -1
	for i, c := range hdr.Contigs {
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
	refBins := tabix.BinFromRange(start, end)
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
			end:  uint64(c.End),
			virt: uint64(c.Beg),
		}
	}
	sort.Slice(pendingChunks, func(i, j int) bool {
		return pendingChunks[i].virt < pendingChunks[j].virt
	})
	_ = pendingChunks // silence unused warning if needed
	var recs []*Record
	for _, c := range pendingChunks {
		chunkRecs, err := r.decodeChunk(c.virt)
		if err != nil {
			return nil, err
		}
		for _, rec := range chunkRecs {
			if rec.Pos >= start && rec.Pos < end {
				recs = append(recs, rec)
			}
		}
	}
	return recs, nil
}

// decodeChunk seeks to virtualOff and decodes VCF records
// until EOF. Records are filtered by region by the caller.
// Stopping at the chunk's End offset would require the
// Decoder to expose per-record byte offsets; for now
// callers filter overlapping records on POS.
func (r *TabixReader) decodeChunk(virtualOff uint64) ([]*Record, error) {
	if err := r.seek.Seek(virtualOff); err != nil {
		return nil, fmt.Errorf("vcf: seek chunk: %w", err)
	}
	dec := NewDecoder(r.seek)
	var recs []*Record
	for {
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