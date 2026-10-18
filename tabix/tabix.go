package tabix

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
)

// Magic is the four-byte Tabix index signature: "TBI\x01".
var Magic = [4]byte{'T', 'B', 'I', 0x01}

// VirtualOffset is a BGZF virtual address: the high 48 bits
// are the compressed-block file offset; the low 16 bits are
// the offset within the decompressed block.
type VirtualOffset uint64

// Chunk is one contiguous block of alignments assigned to a
// bin. Reads at cnk_beg should continue until cnk_end (which
// may be the start of the next chunk or the EOF marker).
type Chunk struct {
	Beg VirtualOffset
	End VirtualOffset
}

// Bin is one bin entry: a bin number (samtools 5-level scheme)
// plus the chunks that cover it.
type Bin struct {
	Number uint32
	Chunks []Chunk
}

// Reference is the index for one reference sequence.
type Reference struct {
	// Bins maps bin number → chunks. Samtools bins go up to
	// 4680 (level 4); level 5 (4681-37448) is the smallest
	// bins used by Tabix in some encoders.
	Bins map[uint32][]Chunk

	// Intervals is the linear index: one BGZF virtual offset
	// per 16 KiB window of the reference. Intervals[i] is the
	// offset of the first alignment whose start lies in
	// window i.
	Intervals []VirtualOffset

	// NoCoor is the count of records with no coordinate
	// (optional in older files; zero if absent).
	NoCoor uint64
}

// Index is the parsed Tabix file.
type Index struct {
	References []*Reference
}

// Reader parses a Tabix index from an underlying [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("sample.bcf.tbi")
//	defer f.Close()
//	idx, err := tabix.NewReader(f).Read()
//	if err != nil { ... }
//	for _, ref := range idx.References {
//	    // ref.Name, ref.Bins, ref.Intervals
//	}
//
// A Reader is not safe for concurrent use.
type Reader struct {
	r io.Reader
}

// NewReader wraps r in a [Reader].
func NewReader(r io.Reader) *Reader {
	return &Reader{r: r}
}

// Read parses the entire Tabix file from r.
func (r *Reader) Read() (*Index, error) {
	var magic [4]byte
	if _, err := io.ReadFull(r.r, magic[:]); err != nil {
		if errors.Is(err, io.EOF) {
			return nil, &ReadError{Err: "empty stream"}
		}
		return nil, &ReadError{Err: "reading magic: " + err.Error()}
	}
	if magic != Magic {
		return nil, &ReadError{Err: fmt.Sprintf("not a Tabix file: got %q", magic[:])}
	}
	var nRef uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nRef); err != nil {
		return nil, &ReadError{Err: "reading n_ref: " + err.Error()}
	}
	idx := &Index{References: make([]*Reference, 0, nRef)}
	for i := uint32(0); i < nRef; i++ {
		ref, err := r.readReference()
		if err != nil {
			return nil, fmt.Errorf("reference %d: %w", i, err)
		}
		idx.References = append(idx.References, ref)
	}
	return idx, nil
}

// readReference parses one reference's index.
func (r *Reader) readReference() (*Reference, error) {
	// Bin index.
	var nBin uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nBin); err != nil {
		return nil, &ReadError{Err: "reading n_bin: " + err.Error()}
	}
	bins := make(map[uint32][]Chunk)
	for i := uint32(0); i < nBin; i++ {
		var bin uint32
		if err := binary.Read(r.r, binary.LittleEndian, &bin); err != nil {
			return nil, fmt.Errorf("reading bin[%d].number: %w", i, err)
		}
		var nChunk uint32
		if err := binary.Read(r.r, binary.LittleEndian, &nChunk); err != nil {
			return nil, fmt.Errorf("reading bin[%d].n_chunk: %w", i, err)
		}
		chunks := make([]Chunk, nChunk)
		for j := uint32(0); j < nChunk; j++ {
			var beg, end uint64
			if err := binary.Read(r.r, binary.LittleEndian, &beg); err != nil {
				return nil, fmt.Errorf("bin[%d].chunk[%d].beg: %w", i, j, err)
			}
			if err := binary.Read(r.r, binary.LittleEndian, &end); err != nil {
				return nil, fmt.Errorf("bin[%d].chunk[%d].end: %w", i, j, err)
			}
			chunks[j] = Chunk{Beg: VirtualOffset(beg), End: VirtualOffset(end)}
		}
		bins[bin] = chunks
	}
	// Linear (interval) index.
	var nIntv uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nIntv); err != nil {
		return nil, &ReadError{Err: "reading n_intv: " + err.Error()}
	}
	intervals := make([]VirtualOffset, nIntv)
	for i := uint32(0); i < nIntv; i++ {
		var ioff uint64
		if err := binary.Read(r.r, binary.LittleEndian, &ioff); err != nil {
			return nil, fmt.Errorf("reading interval[%d]: %w", i, err)
		}
		intervals[i] = VirtualOffset(ioff)
	}
	// Optional n_no_coor field. The Tabix format does not
	// include a version byte, so there is no in-band way to
	// tell whether n_no_coor is present. In practice modern
	// Tabix files (since samtools 0.1.18) always include it.
	//
	// v1 reads n_no_coor unconditionally. For multi-reference
	// files written by an older encoder without the field, this
	// will misalign the next reference's data and produce
	// invalid results — there's no clean way to detect the
	// absence without a version marker. We document the
	// limitation: if you need to read such files, use CSI
	// (which has an explicit version byte).
	ref := &Reference{Bins: bins, Intervals: intervals}
	var noCoor uint64
	if err := binary.Read(r.r, binary.LittleEndian, &noCoor); err != nil {
		if errors.Is(err, io.EOF) || errors.Is(err, io.ErrUnexpectedEOF) {
			// EOF: this is the last reference in the file
			// and there's no n_no_coor (legacy format with
			// a single ref, or already-truncated stream).
		} else {
			return nil, fmt.Errorf("reading n_no_coor: %w", err)
		}
	} else {
		ref.NoCoor = noCoor
	}
	return ref, nil
}

// ReadError describes a failure to parse a Tabix index.
type ReadError struct {
	Err string
}

// Error implements the error interface.
func (e *ReadError) Error() string { return "tabix: " + e.Err }

// Unwrap returns nil: ReadError is a leaf error.
func (e *ReadError) Unwrap() error { return nil }

// BinFromRange returns the bin numbers that may overlap the
// half-open genomic range [beg, end). The bin numbers cover
// the samtools 5-level scheme used by Tabix.
//
// beg and end are 1-based inclusive (GFF/VCF convention) at
// the API level; internally they are converted to 0-based
// half-open for the bin math.
func BinFromRange(beg, end int64) []uint32 {
	return reg2bins(beg, end)
}

// reg2bins implements the samtools reg2bins algorithm: return
// the set of bin numbers that overlap the half-open range
// [beg, end). Used by callers who want to ask "which bins
// might contain my region?" before scanning the chunks.
//
// beg and end are 0-based half-open. The function is
// inclusive of both endpoints (end-1 is included).
func reg2bins(beg, end int64) []uint32 {
	if end <= beg {
		return nil
	}
	end--
	var bins []uint32
	for k := 0; k < 5; k++ {
		switch k {
		case 0:
			bins = append(bins, 0)
		case 1:
			for b := 1 + (beg >> 26); b <= 1+(end>>26); b++ {
				bins = append(bins, uint32(b))
			}
		case 2:
			for b := 9 + (beg >> 23); b <= 9+(end>>23); b++ {
				bins = append(bins, uint32(b))
			}
		case 3:
			for b := 73 + (beg >> 20); b <= 73+(end>>20); b++ {
				bins = append(bins, uint32(b))
			}
		case 4:
			for b := 585 + (beg >> 17); b <= 585+(end>>17); b++ {
				bins = append(bins, uint32(b))
			}
		}
	}
	return bins
}

// LinearIndexOffset returns the linear-index window index for
// a 0-based half-open position. The Tabix linear index uses
// 16 KiB windows (1 << 14); for a region [beg, end), use
// the minimum window index covered by the region as the
// starting point.
func LinearIndexOffset(pos int64) int64 {
	return pos >> 14
}

// LookupChunks returns the chunks whose [Beg, End) range
// overlaps [refBins], or all chunks if refBins is empty.
// This is a thin convenience for callers that don't care
// about which specific bins overlap their region.
func (ref *Reference) LookupChunks(refBins []uint32) []Chunk {
	if len(refBins) == 0 {
		// Return all chunks across all bins.
		var out []Chunk
		for _, chunks := range ref.Bins {
			out = append(out, chunks...)
		}
		return out
	}
	set := make(map[uint32]bool, len(refBins))
	for _, b := range refBins {
		set[b] = true
	}
	var out []Chunk
	for b, chunks := range ref.Bins {
		if set[b] {
			out = append(out, chunks...)
		}
	}
	return out
}
