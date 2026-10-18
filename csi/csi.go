package csi

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
)

// Magic is the four-byte CSI file signature: "CSI\x01".
var Magic = [4]byte{'C', 'S', 'I', 0x01}

// VirtualOffset is a BGZF virtual address: the high 48 bits
// are the compressed-block file offset; the low 16 bits are
// the offset within the decompressed block.
type VirtualOffset = uint64

// Chunk is one contiguous block of alignments assigned to a
// bin.
type Chunk struct {
	Beg VirtualOffset
	End VirtualOffset
}

// Bin is one bin entry: a bin number plus its chunks.
//
// CSI bin numbering uses 5 levels (depth=5) by default; the
// smallest level covers 2^min_shift bp; each successive level
// covers 8× more. The exact bin numbers match the Tabix
// scheme when (min_shift=14, depth=5) is used.
type Bin struct {
	Number uint32
	Chunks []Chunk
}

// Reference is the index for one reference sequence.
type Reference struct {
	// Bins maps bin number → chunks. Same numbering as Tabix
	// when min_shift=14.
	Bins map[uint32][]Chunk

	// Intervals is the linear index: one BGZF virtual offset
	// per (2^MinShift)-bp window.
	Intervals []VirtualOffset
}

// Index is the parsed CSI file.
type Index struct {
	// MinShift is the binning window size exponent
	// (window = 2^MinShift bp).
	MinShift int32

	// Depth is the number of bin levels (default 5).
	Depth int32

	// References is one entry per indexed reference.
	References []*Reference
}

// WindowSize returns the linear-index window size in bp
// (i.e. 2^MinShift).
func (idx *Index) WindowSize() int64 {
	return int64(1) << uint(idx.MinShift)
}

// Reader parses a CSI index from an underlying [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("sample.bcf.csi")
//	defer f.Close()
//	idx, err := csi.NewReader(f).Read()
//	if err != nil { ... }
//	for _, ref := range idx.References {
//	    // ref.Bins, ref.Intervals
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

// Read parses the entire CSI file from r.
func (r *Reader) Read() (*Index, error) {
	var magic [4]byte
	if _, err := io.ReadFull(r.r, magic[:]); err != nil {
		if errors.Is(err, io.EOF) {
			return nil, &ReadError{Err: "empty stream"}
		}
		return nil, &ReadError{Err: "reading magic: " + err.Error()}
	}
	if magic != Magic {
		return nil, &ReadError{Err: fmt.Sprintf("not a CSI file: got %q", magic[:])}
	}
	idx := &Index{}
	if err := binary.Read(r.r, binary.LittleEndian, &idx.MinShift); err != nil {
		return nil, &ReadError{Err: "reading min_shift: " + err.Error()}
	}
	if err := binary.Read(r.r, binary.LittleEndian, &idx.Depth); err != nil {
		return nil, &ReadError{Err: "reading depth: " + err.Error()}
	}
	var nRef uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nRef); err != nil {
		return nil, &ReadError{Err: "reading n_ref: " + err.Error()}
	}
	idx.References = make([]*Reference, 0, nRef)
	for i := uint32(0); i < nRef; i++ {
		ref, err := r.readReference()
		if err != nil {
			return nil, fmt.Errorf("reference %d: %w", i, err)
		}
		idx.References = append(idx.References, ref)
	}
	return idx, nil
}

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
			chunks[j] = Chunk{Beg: beg, End: end}
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
		intervals[i] = ioff
	}
	return &Reference{Bins: bins, Intervals: intervals}, nil
}

// ReadError describes a failure to parse a CSI index.
type ReadError struct {
	Err string
}

// Error implements the error interface.
func (e *ReadError) Error() string { return "csi: " + e.Err }

// Unwrap returns nil: ReadError is a leaf error.
func (e *ReadError) Unwrap() error { return nil }

// BinFromRange returns the bin numbers that may overlap the
// half-open genomic range [beg, end). CSI uses the same
// 5-level binning scheme as Tabix when (min_shift=14,
// depth=5), so the output is the same set of bin numbers
// (but the linear-index window size differs based on the
// file's MinShift).
//
// beg and end are 1-based inclusive (GFF/VCF convention);
// they are converted to 0-based half-open for the bin math.
func BinFromRange(beg, end int64) []uint32 {
	return reg2bins(beg, end)
}

// reg2bins implements the samtools reg2bins algorithm: return
// the set of bin numbers that overlap the half-open range
// [beg, end).
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

// LookupChunks returns the chunks whose [Beg, End) range
// overlaps [refBins], or all chunks if refBins is empty.
func (ref *Reference) LookupChunks(refBins []uint32) []Chunk {
	if len(refBins) == 0 {
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
