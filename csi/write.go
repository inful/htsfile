package csi

import (
	"encoding/binary"
	"fmt"
	"io"
	"sort"
)

// Writer serializes a CSI index to an underlying [io.Writer].
//
// Typical usage:
//
//	var buf bytes.Buffer
//	if err := csi.NewWriter(&buf).Write(idx); err != nil { ... }
//	idx, err := csi.NewReader(&buf).Read()
//
// The Writer is the symmetric counterpart of [Reader]. CSI has
// an explicit version byte (the "CSI\x01" magic) so there is
// no ambiguity in the format — bytes produced by this Writer
// always round-trip through Reader.Read.
//
// A Writer is not safe for concurrent use.
type Writer struct {
	w io.Writer
}

// NewWriter wraps w in a [Writer].
func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w}
}

// Write serializes idx to the underlying writer: magic +
// min_shift + depth + n_ref + per-reference bin index + linear
// index. CSI does not have an n_no_coor field (it uses the
// version byte to disambiguate format variants, so EOF is the
// sole end-of-stream signal). Bin numbers are emitted in
// ascending order and chunks within each bin are sorted by Beg
// to produce deterministic output matching the Reader's
// post-sort behavior.
func (wr *Writer) Write(idx *Index) error {
	if idx == nil {
		return fmt.Errorf("csi: Write: nil index")
	}
	if _, err := wr.w.Write(Magic[:]); err != nil {
		return fmt.Errorf("writing magic: %w", err)
	}
	// min_shift and depth are int32 LE per the CSI spec.
	var minShiftBytes [4]byte
	binary.LittleEndian.PutUint32(minShiftBytes[:], uint32(idx.MinShift))
	if _, err := wr.w.Write(minShiftBytes[:]); err != nil {
		return fmt.Errorf("writing min_shift: %w", err)
	}
	var depthBytes [4]byte
	binary.LittleEndian.PutUint32(depthBytes[:], uint32(idx.Depth))
	if _, err := wr.w.Write(depthBytes[:]); err != nil {
		return fmt.Errorf("writing depth: %w", err)
	}
	var nRef [4]byte
	binary.LittleEndian.PutUint32(nRef[:], uint32(len(idx.References)))
	if _, err := wr.w.Write(nRef[:]); err != nil {
		return fmt.Errorf("writing n_ref: %w", err)
	}
	for i, ref := range idx.References {
		if err := wr.writeReference(ref); err != nil {
			return fmt.Errorf("reference %d: %w", i, err)
		}
	}
	return nil
}

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
			binary.LittleEndian.PutUint64(begBytes[:], c.Beg)
			if _, err := wr.w.Write(begBytes[:]); err != nil {
				return fmt.Errorf("writing chunk beg: %w", err)
			}
			var endBytes [8]byte
			binary.LittleEndian.PutUint64(endBytes[:], c.End)
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
		binary.LittleEndian.PutUint64(ivBytes[:], iv)
		if _, err := wr.w.Write(ivBytes[:]); err != nil {
			return fmt.Errorf("writing interval[%d]: %w", i, err)
		}
	}
	return nil
}
