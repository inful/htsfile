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
// The Writer is the symmetric counterpart of [Reader]: it
// emits the full CSI v1 format with l_aux, per-bin loffset,
// and per-ref n_no_coor. Bytes always round-trip through
// Reader.Read.
//
// A Writer is not safe for concurrent use.
type Writer struct {
	w io.Writer
}

// NewWriter wraps w in a [Writer].
func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w}
}

// Write serializes idx to the underlying writer: magic,
// min_shift, depth, l_aux, aux, n_ref, then per-reference
// bin index (with per-bin loffset), linear index, and
// n_no_coor. Bin numbers are emitted in ascending order and
// chunks within each bin are sorted by Beg to produce
// deterministic output matching the Reader's post-sort
// behavior.
func (wr *Writer) Write(idx *Index) error {
	if idx == nil {
		return fmt.Errorf("csi: Write: nil index")
	}
	if _, err := wr.w.Write(Magic[:]); err != nil {
		return fmt.Errorf("writing magic: %w", err)
	}
	var i32 [4]byte
	binary.LittleEndian.PutUint32(i32[:], uint32(idx.MinShift))
	if _, err := wr.w.Write(i32[:]); err != nil {
		return fmt.Errorf("writing min_shift: %w", err)
	}
	binary.LittleEndian.PutUint32(i32[:], uint32(idx.Depth))
	if _, err := wr.w.Write(i32[:]); err != nil {
		return fmt.Errorf("writing depth: %w", err)
	}
	binary.LittleEndian.PutUint32(i32[:], uint32(len(idx.Aux)))
	if _, err := wr.w.Write(i32[:]); err != nil {
		return fmt.Errorf("writing l_aux: %w", err)
	}
	if len(idx.Aux) > 0 {
		if _, err := wr.w.Write(idx.Aux); err != nil {
			return fmt.Errorf("writing aux: %w", err)
		}
	}
	binary.LittleEndian.PutUint32(i32[:], uint32(len(idx.References)))
	if _, err := wr.w.Write(i32[:]); err != nil {
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
	binNumbers := make([]uint32, 0, len(ref.Bins))
	for bin := range ref.Bins {
		binNumbers = append(binNumbers, bin)
	}
	sort.Slice(binNumbers, func(i, j int) bool {
		return binNumbers[i] < binNumbers[j]
	})
	var i32 [4]byte
	binary.LittleEndian.PutUint32(i32[:], uint32(len(binNumbers)))
	if _, err := wr.w.Write(i32[:]); err != nil {
		return fmt.Errorf("writing n_bin: %w", err)
	}
	for _, bin := range binNumbers {
		chunks := ref.Bins[bin]
		sorted := make([]Chunk, len(chunks))
		copy(sorted, chunks)
		sort.Slice(sorted, func(i, j int) bool {
			return sorted[i].Beg < sorted[j].Beg
		})
		binary.LittleEndian.PutUint32(i32[:], bin)
		if _, err := wr.w.Write(i32[:]); err != nil {
			return fmt.Errorf("writing bin %d number: %w", bin, err)
		}
		// Per-bin loffset: virtual file offset of the first
		// overlapping record. The Reader discards this; we
		// emit zero, which is the value samtools uses when
		// the linear index covers the lookup path.
		var u64 [8]byte
		binary.LittleEndian.PutUint64(u64[:], 0)
		if _, err := wr.w.Write(u64[:]); err != nil {
			return fmt.Errorf("writing bin %d loffset: %w", bin, err)
		}
		binary.LittleEndian.PutUint32(i32[:], uint32(len(sorted)))
		if _, err := wr.w.Write(i32[:]); err != nil {
			return fmt.Errorf("writing bin %d n_chunk: %w", bin, err)
		}
		for _, c := range sorted {
			binary.LittleEndian.PutUint64(u64[:], c.Beg)
			if _, err := wr.w.Write(u64[:]); err != nil {
				return fmt.Errorf("writing chunk beg: %w", err)
			}
			binary.LittleEndian.PutUint64(u64[:], c.End)
			if _, err := wr.w.Write(u64[:]); err != nil {
				return fmt.Errorf("writing chunk end: %w", err)
			}
		}
	}
	binary.LittleEndian.PutUint32(i32[:], uint32(len(ref.Intervals)))
	if _, err := wr.w.Write(i32[:]); err != nil {
		return fmt.Errorf("writing n_intv: %w", err)
	}
	for i, iv := range ref.Intervals {
		var u64 [8]byte
		binary.LittleEndian.PutUint64(u64[:], iv)
		if _, err := wr.w.Write(u64[:]); err != nil {
			return fmt.Errorf("writing interval[%d]: %w", i, err)
		}
	}
	// n_no_coor: per-spec, optional. The reader tolerates
	// absence (treats as zero); we always emit for forward
	// compatibility.
	var u64 [8]byte
	binary.LittleEndian.PutUint64(u64[:], ref.NoCoor)
	if _, err := wr.w.Write(u64[:]); err != nil {
		return fmt.Errorf("writing n_no_coor: %w", err)
	}
	return nil
}
