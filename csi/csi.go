package csi

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"sort"
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

// Bin is one bin entry: a bin number plus its chunks. CSI
// bin numbering uses a configurable number of levels (the
// file's Depth, which counts the level number of the smallest
// bin; depth=5 means six levels — 0..5). For (min_shift=14,
// depth=5) the bin numbers match the Tabix scheme up to
// bin 37448.
type Bin struct {
	Number uint32
	Chunks []Chunk
}

// Reference is the index for one reference sequence.
type Reference struct {
	// Name carries the NUL-stripped name from the index
	// header when the file has an aux blob.
	Name string

	// Bins maps bin number to chunks.
	Bins map[uint32][]Chunk

	// Intervals is the linear index: one BGZF virtual offset
	// per (2^MinShift)-bp window.
	Intervals []VirtualOffset

	// NoCoor is the count of unmapped unplaced reads (RNAME *).
	// Zero when the field is absent from the file.
	NoCoor uint64
}

// Index is the parsed CSI file.
type Index struct {
	MinShift int32
	Depth    int32
	Aux      []byte

	References []*Reference
}

// WindowSize returns the linear-index window size in bp
// (i.e. 2^MinShift).
func (idx *Index) WindowSize() int64 {
	return int64(1) << uint(idx.MinShift)
}

// BinLimit returns the maximum valid bin number for this
// index ((2^((depth+1)*3) - 1) / 7). For the canonical
// (min_shift=14, depth=5) configuration this is 37448.
func (idx *Index) BinLimit() int {
	return (1<<((idx.Depth+1)*3) - 1) / 7
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
			return nil, &ReadError{Err: errors.New("empty stream")}
		}
		return nil, &ReadError{Err: fmt.Errorf("reading magic: %w", err)}
	}
	if magic != Magic {
		return nil, &ReadError{Err: fmt.Errorf("not a CSI file: got %q", magic[:])}
	}
	idx := &Index{}
	if err := binary.Read(r.r, binary.LittleEndian, &idx.MinShift); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading min_shift: %w", err)}
	}
	if err := binary.Read(r.r, binary.LittleEndian, &idx.Depth); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading depth: %w", err)}
	}
	var lAux int32
	if err := binary.Read(r.r, binary.LittleEndian, &lAux); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading l_aux: %w", err)}
	}
	if lAux < 0 {
		return nil, &ReadError{Err: fmt.Errorf("negative l_aux: %d", lAux)}
	}
	if lAux > 0 {
		idx.Aux = make([]byte, lAux)
		if _, err := io.ReadFull(r.r, idx.Aux); err != nil {
			return nil, &ReadError{Err: fmt.Errorf("reading aux: %w", err)}
		}
	}
	var nRef uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nRef); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading n_ref: %w", err)}
	}
	idx.References = make([]*Reference, 0, nRef)
	for i := uint32(0); i < nRef; i++ {
		ref, err := r.readReference(idx.MinShift, idx.Depth)
		if err != nil {
			return nil, fmt.Errorf("reference %d: %w", i, err)
		}
		idx.References = append(idx.References, ref)
	}
	return idx, nil
}

// readReference parses one reference's index. The CSI spec
// requires a per-bin loffset (uint64) between bin and n_chunk;
// the optional n_no_coor (uint64) at the end of the
// reference may be absent in legacy files.
func (r *Reader) readReference(minShift, depth int32) (*Reference, error) {
	var nBin uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nBin); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading n_bin: %w", err)}
	}
	bins := make(map[uint32][]Chunk, nBin)
	for i := uint32(0); i < nBin; i++ {
		var bin uint32
		if err := binary.Read(r.r, binary.LittleEndian, &bin); err != nil {
			return nil, fmt.Errorf("bin[%d].number: %w", i, err)
		}
		// Per-bin loffset: virtual file offset of the first
		// overlapping record. The spec requires it; this reader
		// discards the value (the linear index already covers
		// the lookup path) but must consume the bytes or the
		// stream misaligns.
		var loffset uint64
		if err := binary.Read(r.r, binary.LittleEndian, &loffset); err != nil {
			return nil, fmt.Errorf("bin[%d].loffset: %w", i, err)
		}
		var nChunk uint32
		if err := binary.Read(r.r, binary.LittleEndian, &nChunk); err != nil {
			return nil, fmt.Errorf("bin[%d].n_chunk: %w", i, err)
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
		_ = loffset
	}
	var nIntv uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nIntv); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading n_intv: %w", err)}
	}
	intervals := make([]VirtualOffset, nIntv)
	for i := uint32(0); i < nIntv; i++ {
		var ioff uint64
		if err := binary.Read(r.r, binary.LittleEndian, &ioff); err != nil {
			return nil, fmt.Errorf("interval[%d]: %w", i, err)
		}
		intervals[i] = ioff
	}
	ref := &Reference{Bins: bins, Intervals: intervals}
	// Optional n_no_coor: EOF on the trailing uint64 is treated
	// as zero in legacy files. This branch only matters on the
	// very last reference; the next reference's loop will see
	// a short read and error out cleanly.
	var noCoor uint64
	err := binary.Read(r.r, binary.LittleEndian, &noCoor)
	switch {
	case err == nil:
		ref.NoCoor = noCoor
	case errors.Is(err, io.EOF), errors.Is(err, io.ErrUnexpectedEOF):
	default:
		return nil, fmt.Errorf("reading n_no_coor: %w", err)
	}
	_ = minShift
	_ = depth
	return ref, nil
}

// ReadError describes a failure to parse a CSI index.
type ReadError struct {
	Err error
}

// Error implements the error interface.
func (e *ReadError) Error() string {
	return "csi: " + e.Err.Error()
}

// Unwrap returns the underlying error.
func (e *ReadError) Unwrap() error { return e.Err }

// BinFromRange returns the bin numbers that may overlap the
// half-open genomic range [beg, end). CSI uses the same
// 6-level scheme as Tabix when (min_shift=14, depth=5), so
// the output is the same set of bin numbers. The shift/offset
// tables for level 5 (4681-37448) are computed at the standard
// min_shift=14; callers using a different min_shift can call
// reg2bins directly with the file's parameters.
//
// beg and end are 0-based half-open.
func BinFromRange(beg, end int64) []uint32 {
	return reg2bins(beg, end, 14, 5)
}

// reg2bins implements the samtools CSI reg2bins algorithm
// generalized for arbitrary (min_shift, depth). The number of
// levels emitted is depth+1, matching the spec's "level 0
// through level depth" (inclusive) scheme.
//
// beg and end are 0-based half-open. The function is
// inclusive of both endpoints (end-1 is included).
func reg2bins(beg, end int64, minShift, depth int32) []uint32 { //nolint:unparam // callers may pass other min_shift/depth
	if end <= beg {
		return nil
	}
	end--
	bins := []uint32{}
	s := uint64(minShift) + uint64(depth)*3
	t := uint64(0)
	for l := int32(0); l <= depth; l++ {
		lo := uint32(t + (uint64(beg) >> s))
		hi := uint32(t + (uint64(end) >> s))
		for b := lo; b <= hi; b++ {
			bins = append(bins, b)
		}
		s -= 3
		t += 1 << (uint64(l) * 3) // update for the NEXT iteration
	}
	sort.Slice(bins, func(i, j int) bool { return bins[i] < bins[j] })
	return bins
}

// LookupChunks returns the chunks whose bin number is in
// refBins, or all chunks if refBins is empty.
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
