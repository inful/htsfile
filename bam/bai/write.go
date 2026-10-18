package bai

import (
	"encoding/binary"
	"fmt"
	"io"
	"sort"
)

// Writer serializes a BAI index to an underlying [io.Writer].
//
// Typical usage:
//
//	var buf bytes.Buffer
//	if err := bai.NewWriter(&buf).Write(idx); err != nil { ... }
//	idx, err := bai.NewReader(&buf).Read()
//
// The Writer is the symmetric counterpart of [Reader]: bytes
// it produces round-trip cleanly through Reader.Read. Bins are
// emitted in ascending Number order to match the linear scan
// the Reader uses.
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
// n_ref + per-reference (n_bin + bins + n_intv + intervals
// + n_no_coor).
func (wr *Writer) Write(idx *Index) error {
	if idx == nil {
		return fmt.Errorf("bai: Write: nil index")
	}
	if _, err := wr.w.Write(Magic[:]); err != nil {
		return fmt.Errorf("writing magic: %w", err)
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
	// Sort bins by Number ascending; the Reader scans bins in
	// order and assumes monotonically increasing numbers.
	bins := make([]*Bin, len(ref.Bins))
	copy(bins, ref.Bins)
	sort.Slice(bins, func(i, j int) bool {
		return bins[i].Number < bins[j].Number
	})
	var nBin [4]byte
	binary.LittleEndian.PutUint32(nBin[:], uint32(len(bins)))
	if _, err := wr.w.Write(nBin[:]); err != nil {
		return fmt.Errorf("writing n_bin: %w", err)
	}
	for _, b := range bins {
		var bnum [4]byte
		binary.LittleEndian.PutUint32(bnum[:], b.Number)
		if _, err := wr.w.Write(bnum[:]); err != nil {
			return fmt.Errorf("writing bin.number: %w", err)
		}
		var lo [8]byte
		binary.LittleEndian.PutUint64(lo[:], b.Lo)
		if _, err := wr.w.Write(lo[:]); err != nil {
			return fmt.Errorf("writing bin.lo: %w", err)
		}
		var hi [8]byte
		binary.LittleEndian.PutUint64(hi[:], b.Hi)
		if _, err := wr.w.Write(hi[:]); err != nil {
			return fmt.Errorf("writing bin.hi: %w", err)
		}
	}
	// Intervals.
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
	// n_no_coor (always emitted — the Reader probes for it
	// with EOF tolerance so legacy files without it still
	// parse, but the Writer always writes the modern form).
	var noCoor [8]byte
	binary.LittleEndian.PutUint64(noCoor[:], ref.NoCoor)
	if _, err := wr.w.Write(noCoor[:]); err != nil {
		return fmt.Errorf("writing n_no_coor: %w", err)
	}
	return nil
}
