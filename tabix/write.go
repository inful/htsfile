package tabix

import (
	"encoding/binary"
	"fmt"
	"io"
	"sort"
)

// Writer serializes a Tabix index to an underlying [io.Writer].
//
// Typical usage:
//
//	var buf bytes.Buffer
//	if err := tabix.NewWriter(&buf).Write(idx); err != nil { ... }
//	idx, err := tabix.NewReader(&buf).Read()
//
// The Writer is the symmetric counterpart of [Reader]: the bytes
// it produces round-trip cleanly through Reader.Read, with the
// caveat documented there about the optional n_no_coor field
// (always emitted by this Writer).
//
// A Writer is not safe for concurrent use.
type Writer struct {
	w io.Writer
}

// NewWriter wraps w in a [Writer].
func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w}
}

// Write serializes idx to the underlying writer. The complete
// file is emitted in one call: magic + n_ref + per-reference
// bin index + linear index + n_no_coor. The reference's bin
// chunks are sorted by Beg to produce deterministic output
// (the Reader also sorts on read; matching the sort keeps
// byte-for-byte round-trip stable).
//
// n_no_coor is always emitted (matching modern samtools Tabix
// output since 0.1.18). Callers that need legacy single-ref
// files without the trailer can post-process the output.
func (wr *Writer) Write(idx *Index) error {
	if idx == nil {
		return fmt.Errorf("tabix: Write: nil index")
	}
	// Magic.
	if _, err := wr.w.Write(Magic[:]); err != nil {
		return fmt.Errorf("writing magic: %w", err)
	}
	// n_ref.
	var nRef [4]byte
	binary.LittleEndian.PutUint32(nRef[:], uint32(len(idx.References)))
	if _, err := wr.w.Write(nRef[:]); err != nil {
		return fmt.Errorf("writing n_ref: %w", err)
	}
	// Per-reference data.
	for i, ref := range idx.References {
		if err := wr.writeReference(ref); err != nil {
			return fmt.Errorf("reference %d: %w", i, err)
		}
	}
	return nil
}

// writeReference serializes one reference's bin index + linear
// index + n_no_coor. Bins are emitted in ascending bin-number
// order for deterministic output; chunks within each bin are
// sorted by Beg.
func (wr *Writer) writeReference(ref *Reference) error {
	// Bin index.
	binNumbers := make([]uint32, 0, len(ref.Bins))
	for bin := range ref.Bins {
		binNumbers = append(binNumbers, bin)
	}
	sort.Slice(binNumbers, func(i, j int) bool {
		return binNumbers[i] < binNumbers[j]
	})
	var nBin [4]byte
	binary.LittleEndian.PutUint32(nBin[:], uint32(len(binNumbers)))
	if _, err := wr.w.Write(nBin[:]); err != nil {
		return fmt.Errorf("writing n_bin: %w", err)
	}
	for _, bin := range binNumbers {
		chunks := ref.Bins[bin]
		// Sort chunks by Beg for deterministic output.
		sorted := make([]Chunk, len(chunks))
		copy(sorted, chunks)
		sort.Slice(sorted, func(i, j int) bool {
			return sorted[i].Beg < sorted[j].Beg
		})
		var binBytes [4]byte
		binary.LittleEndian.PutUint32(binBytes[:], bin)
		if _, err := wr.w.Write(binBytes[:]); err != nil {
			return fmt.Errorf("writing bin %d number: %w", bin, err)
		}
		var nChunk [4]byte
		binary.LittleEndian.PutUint32(nChunk[:], uint32(len(sorted)))
		if _, err := wr.w.Write(nChunk[:]); err != nil {
			return fmt.Errorf("writing bin %d n_chunk: %w", bin, err)
		}
		for _, c := range sorted {
			var begBytes [8]byte
			binary.LittleEndian.PutUint64(begBytes[:], uint64(c.Beg))
			if _, err := wr.w.Write(begBytes[:]); err != nil {
				return fmt.Errorf("writing chunk beg: %w", err)
			}
			var endBytes [8]byte
			binary.LittleEndian.PutUint64(endBytes[:], uint64(c.End))
			if _, err := wr.w.Write(endBytes[:]); err != nil {
				return fmt.Errorf("writing chunk end: %w", err)
			}
		}
	}
	// Linear (interval) index.
	var nIntv [4]byte
	binary.LittleEndian.PutUint32(nIntv[:], uint32(len(ref.Intervals)))
	if _, err := wr.w.Write(nIntv[:]); err != nil {
		return fmt.Errorf("writing n_intv: %w", err)
	}
	for i, iv := range ref.Intervals {
		var ivBytes [8]byte
		binary.LittleEndian.PutUint64(ivBytes[:], uint64(iv))
		if _, err := wr.w.Write(ivBytes[:]); err != nil {
			return fmt.Errorf("writing interval[%d]: %w", i, err)
		}
	}
	// n_no_coor (always emitted — see Reader doc for the
	// multi-ref legacy limitation).
	var noCoorBytes [8]byte
	binary.LittleEndian.PutUint64(noCoorBytes[:], ref.NoCoor)
	if _, err := wr.w.Write(noCoorBytes[:]); err != nil {
		return fmt.Errorf("writing n_no_coor: %w", err)
	}
	return nil
}
