// Package tabix provides a reader and writer for the Tabix
// index format (file extension .tbi), the general-purpose
// index for sorted, bgzipped, position-sorted text files
// (VCF, GFF, BED, SAM, ...).
//
// The Tabix format is line-oriented and per-reference (per
// hts-specs/tabix.tex):
//
//	magic      : "TBI\x01"            (4 bytes)
//	n_ref      : int32 little-endian
//	format     : int32   (0 generic; 1 SAM; 2 VCF; 0x10000 = BED)
//	col_seq    : int32   1-based column of the sequence name
//	col_beg    : int32   1-based column of the start
//	col_end    : int32   1-based column of the end (0 if implicit)
//	meta       : int32   first character of comment lines
//	skip       : int32   # of lines to skip at the top
//	l_nm       : int32   length of the concatenated reference names
//	names      : char[l_nm]   NUL-terminated sequence names
//	For each reference:
//	    n_bin       : int32
//	    For each bin:
//	        bin     : uint32
//	        n_chunk : int32
//	        For each chunk:
//	            cnk_beg : uint64 BGZF virtual offset
//	            cnk_end : uint64 BGZF virtual offset
//	    n_intv      : int32
//	    For each interval:
//	        ioff       : uint64 BGZF virtual offset
//	    n_no_coor    : uint64   (count of unmapped records;
//	                              may be absent in very old files
//	                              with a single reference; this
//	                              reader tolerates that, the writer
//	                              always emits it)
//
// reg2bins returns all six levels of the samtools scheme
// (bins 0, 1-8, 9-72, 73-584, 585-4680, 4681-37448).
package tabix

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"sort"
	"strings"
)

// Magic is the four-byte Tabix index signature: "TBI\x01".
var Magic = [4]byte{'T', 'B', 'I', 0x01}

// Format constants for the per-file format field.
const (
	FormatGeneric int32 = 0
	FormatSAM     int32 = 1
	FormatVCF     int32 = 2
	FormatBEDBit  int32 = 0x10000 // bit set on FormatGeneric for BED
)

// VirtualOffset is a BGZF virtual address: the high 48 bits
// are the compressed-block file offset; the low 16 bits are
// the offset within the decompressed block.
type VirtualOffset uint64

// Chunk is one contiguous block of alignments assigned to a
// bin.
type Chunk struct {
	Beg VirtualOffset
	End VirtualOffset
}

// Reference is the index for one reference sequence.
type Reference struct {
	Name      string
	Bins      map[uint32][]Chunk
	Intervals []VirtualOffset
	NoCoor    uint64
}

// Index is the parsed Tabix file.
type Index struct {
	Format     int32
	ColSeq     int32
	ColBeg     int32
	ColEnd     int32
	MetaChar   int32
	SkipLines  int32
	Names      []string
	References []*Reference
}

// Reader parses a Tabix index from an underlying io.Reader.
//
// A Reader is not safe for concurrent use.
type Reader struct {
	r io.Reader
}

// NewReader returns a Reader reading from r.
func NewReader(r io.Reader) *Reader {
	return &Reader{r: r}
}

// Read parses the entire Tabix file from r.
func (r *Reader) Read() (*Index, error) {
	var magic [4]byte
	if _, err := io.ReadFull(r.r, magic[:]); err != nil {
		if errors.Is(err, io.EOF) {
			return nil, &ReadError{Err: errors.New("empty stream")}
		}
		return nil, &ReadError{Err: fmt.Errorf("reading magic: %w", err)}
	}
	if magic != Magic {
		return nil, &ReadError{Err: fmt.Errorf("not a Tabix file: got %q", magic[:])}
	}
	var nRef uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nRef); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading n_ref: %w", err)}
	}
	idx := &Index{References: make([]*Reference, 0, nRef)}
	// Eight header fields between n_ref and the name table.
	var format, colSeq, colBeg, colEnd, meta, skip int32
	if err := binary.Read(r.r, binary.LittleEndian, &format); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading format: %w", err)}
	}
	if err := binary.Read(r.r, binary.LittleEndian, &colSeq); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading col_seq: %w", err)}
	}
	if err := binary.Read(r.r, binary.LittleEndian, &colBeg); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading col_beg: %w", err)}
	}
	if err := binary.Read(r.r, binary.LittleEndian, &colEnd); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading col_end: %w", err)}
	}
	if err := binary.Read(r.r, binary.LittleEndian, &meta); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading meta: %w", err)}
	}
	if err := binary.Read(r.r, binary.LittleEndian, &skip); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading skip: %w", err)}
	}
	idx.Format = format
	idx.ColSeq = colSeq
	idx.ColBeg = colBeg
	idx.ColEnd = colEnd
	idx.MetaChar = meta
	idx.SkipLines = skip

	var lNm int32
	if err := binary.Read(r.r, binary.LittleEndian, &lNm); err != nil {
		return nil, &ReadError{Err: fmt.Errorf("reading l_nm: %w", err)}
	}
	if lNm < 0 {
		return nil, &ReadError{Err: fmt.Errorf("negative l_nm: %d", lNm)}
	}
	if lNm > 0 {
		nameBuf := make([]byte, lNm)
		if _, err := io.ReadFull(r.r, nameBuf); err != nil {
			return nil, &ReadError{Err: fmt.Errorf("reading names: %w", err)}
		}
		idx.Names = splitNames(nameBuf)
	}

	for i := uint32(0); i < nRef; i++ {
		ref, err := r.readReference()
		if err != nil {
			return nil, fmt.Errorf("reference %d: %w", i, err)
		}
		idx.References = append(idx.References, ref)
	}
	for i, ref := range idx.References {
		if i < len(idx.Names) {
			ref.Name = idx.Names[i]
		}
	}
	return idx, nil
}

func splitNames(buf []byte) []string {
	var out []string
	for len(buf) > 0 {
		i := strings.IndexByte(string(buf), 0)
		if i < 0 {
			out = append(out, string(buf))
			break
		}
		out = append(out, string(buf[:i]))
		buf = buf[i+1:]
	}
	return out
}

// readReference parses one reference's index. The trailing
// n_no_coor (uint64) is present in every modern file but
// may be absent in single-reference legacy files; the reader
// distinguishes the two by treating EOF on the trailing
// field as zero.
func (r *Reader) readReference() (*Reference, error) {
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
			chunks[j] = Chunk{Beg: VirtualOffset(beg), End: VirtualOffset(end)}
		}
		bins[bin] = chunks
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
		intervals[i] = VirtualOffset(ioff)
	}
	ref := &Reference{Bins: bins, Intervals: intervals}
	var noCoor uint64
	err := binary.Read(r.r, binary.LittleEndian, &noCoor)
	switch {
	case err == nil:
		ref.NoCoor = noCoor
	case errors.Is(err, io.EOF), errors.Is(err, io.ErrUnexpectedEOF):
	default:
		return nil, fmt.Errorf("reading n_no_coor: %w", err)
	}
	return ref, nil
}

// ReadError describes a failure to parse a Tabix index. It
// implements errors.Unwrap so callers can match the
// underlying io.EOF / io.ErrUnexpectedEOF if they need to
// distinguish the legacy "no n_no_coor" case from a true
// truncation.
type ReadError struct {
	Err error
}

// Error implements error.
func (e *ReadError) Error() string {
	return "tabix: " + e.Err.Error()
}

// Unwrap returns the underlying error.
func (e *ReadError) Unwrap() error { return e.Err }

// BinFromRange returns the set of bin numbers that overlap
// the half-open range [beg, end). It is the public form of
// the samtools reg2bins algorithm. beg and end are 0-based
// half-open; the function is inclusive of both endpoints
// (end-1 is included). It produces the full 6-level scheme:
// bins 0, 1-8, 9-72, 73-584, 585-4680, and 4681-37448.
//
// Old callers asking for "5 levels" (bins 0..4680) get the
// correct 6-level set; a typical 4 KiB region maps to a
// single bin in the smallest (16 KiB) level, so the new
// behaviour is a strict refinement.
func BinFromRange(beg, end int64) []uint32 {
	if end <= beg {
		return nil
	}
	end--
	bins := []uint32{0}
	shifts := []uint32{26, 23, 20, 17, 14}
	offsets := []uint32{1, 9, 73, 585, 4681}
	for i, shift := range shifts {
		lo := offsets[i] + uint32(beg>>shift)
		hi := offsets[i] + uint32(end>>shift)
		for b := lo; b <= hi; b++ {
			bins = append(bins, b)
		}
	}
	sort.Slice(bins, func(i, j int) bool { return bins[i] < bins[j] })
	return bins
}

// LinearIndexOffset returns the linear-index window index for
// a 0-based half-open position. The Tabix linear index uses
// 16 KiB windows (1 << 14).
func LinearIndexOffset(pos int64) int64 {
	return pos >> 14
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
