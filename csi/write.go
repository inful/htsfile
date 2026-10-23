package csi

import (
	"encoding/binary"
	"fmt"
	"github.com/inful/htsfile/internal/bgzf"
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
// bin index (with per-bin loffset), and n_no_coor.
//
// A CSI has no linear index: the per-bin loffset replaces it.
// Bin numbers are emitted in ascending order and chunks
// within each bin are sorted by Beg to produce deterministic
// output matching the Reader's post-sort behavior.
func (wr *Writer) Write(idx *Index) error {
	if idx == nil {
		return fmt.Errorf("csi: Write: nil index")
	}
	// htslib writes the index BGZF-compressed; an
	// uncompressed .csi is rejected by htslib outright.
	zw := bgzf.NewWriter(wr.w)
	orig := wr.w
	wr.w = zw
	defer func() { wr.w = orig }()
	if err := wr.writeIndex(idx); err != nil {
		_ = zw.Close()
		return err
	}
	if err := zw.Close(); err != nil {
		return fmt.Errorf("csi: closing bgzf stream: %w", err)
	}
	return nil
}

// writeIndex serialises idx to the writer's current target.
func (wr *Writer) writeIndex(idx *Index) error {
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
		if err := wr.writeReference(ref, idx.Depth); err != nil {
			return fmt.Errorf("reference %d: %w", i, err)
		}
	}
	// n_no_coor: written once, at the end of the file.
	var u64 [8]byte
	binary.LittleEndian.PutUint64(u64[:], idx.NoCoor)
	if _, err := wr.w.Write(u64[:]); err != nil {
		return fmt.Errorf("writing n_no_coor: %w", err)
	}
	return nil
}

func (wr *Writer) writeReference(ref *Reference, depth int32) error {
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
	metaBin := uint32(BinLimit(depth) + 1)
	for _, bin := range binNumbers {
		chunks := ref.Bins[bin]
		sorted := make([]Chunk, len(chunks))
		copy(sorted, chunks)
		if bin != metaBin {
			// The meta bin holds two pseudo-chunks - the
			// file's byte range and the record counts -
			// whose order is fixed by the format.
			sort.Slice(sorted, func(i, j int) bool {
				return sorted[i].Beg < sorted[j].Beg
			})
		}
		binary.LittleEndian.PutUint32(i32[:], bin)
		if _, err := wr.w.Write(i32[:]); err != nil {
			return fmt.Errorf("writing bin %d number: %w", bin, err)
		}
		// Per-bin loffset: the virtual file offset of the
		// first record in the bin. A CSI has no linear index,
		// so this is what tells a reader where the bin's
		// records start. When the index does not carry one,
		// the first chunk's start is the same value.
		loffset, ok := ref.Loffsets[bin]
		if !ok && len(sorted) > 0 {
			loffset = sorted[0].Beg
		}
		var u64 [8]byte
		binary.LittleEndian.PutUint64(u64[:], loffset)
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
	// A CSI has no linear index: the per-bin loffset replaces
	// it, and htslib writes none.
	return nil
}
