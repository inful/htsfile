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
// it produces round-trip cleanly through Reader.Read on the
// six-level scheme that this package supports.
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
// file is emitted in one call: magic, n_ref, the eight header
// fields, the concatenated name table, then per-reference
// bin index, linear index, and n_no_coor. The reference's bin
// chunks are sorted by Beg to produce deterministic output.
func (wr *Writer) Write(idx *Index) error {
	if idx == nil {
		return fmt.Errorf("tabix: Write: nil index")
	}
	if _, err := wr.w.Write(Magic[:]); err != nil {
		return fmt.Errorf("writing magic: %w", err)
	}
	var nRef [4]byte
	binary.LittleEndian.PutUint32(nRef[:], uint32(len(idx.References)))
	if _, err := wr.w.Write(nRef[:]); err != nil {
		return fmt.Errorf("writing n_ref: %w", err)
	}
	if err := wr.writeInt32(idx.Format); err != nil {
		return fmt.Errorf("writing format: %w", err)
	}
	if err := wr.writeInt32(idx.ColSeq); err != nil {
		return fmt.Errorf("writing col_seq: %w", err)
	}
	if err := wr.writeInt32(idx.ColBeg); err != nil {
		return fmt.Errorf("writing col_beg: %w", err)
	}
	if err := wr.writeInt32(idx.ColEnd); err != nil {
		return fmt.Errorf("writing col_end: %w", err)
	}
	if err := wr.writeInt32(idx.MetaChar); err != nil {
		return fmt.Errorf("writing meta: %w", err)
	}
	if err := wr.writeInt32(idx.SkipLines); err != nil {
		return fmt.Errorf("writing skip: %w", err)
	}
	// Concatenated reference names, NUL-terminated.
	var nameBuf []byte
	if len(idx.Names) == 0 && len(idx.References) > 0 {
		// Fall back to Reference.Name if Names was not
		// populated. This keeps the writer working for
		// callers that build an Index by setting per-ref
		// Name only.
		for _, ref := range idx.References {
			nameBuf = append(nameBuf, []byte(ref.Name)...)
			nameBuf = append(nameBuf, 0)
		}
	} else {
		for _, name := range idx.Names {
			nameBuf = append(nameBuf, []byte(name)...)
			nameBuf = append(nameBuf, 0)
		}
	}
	if err := wr.writeInt32(int32(len(nameBuf))); err != nil {
		return fmt.Errorf("writing l_nm: %w", err)
	}
	if len(nameBuf) > 0 {
		if _, err := wr.w.Write(nameBuf); err != nil {
			return fmt.Errorf("writing names: %w", err)
		}
	}

	for i, ref := range idx.References {
		if err := wr.writeReference(ref); err != nil {
			return fmt.Errorf("reference %d: %w", i, err)
		}
	}
	return nil
}

func (wr *Writer) writeInt32(v int32) error {
	var b [4]byte
	binary.LittleEndian.PutUint32(b[:], uint32(v))
	_, err := wr.w.Write(b[:])
	return err
}

func (wr *Writer) writeUint64(v uint64) error {
	var b [8]byte
	binary.LittleEndian.PutUint64(b[:], v)
	_, err := wr.w.Write(b[:])
	return err
}

// writeReference serializes one reference's bin index, linear
// index, and n_no_coor. Bins are emitted in ascending bin-number
// order for deterministic output; chunks within each bin are
// sorted by Beg.
func (wr *Writer) writeReference(ref *Reference) error {
	binNumbers := make([]uint32, 0, len(ref.Bins))
	for bin := range ref.Bins {
		binNumbers = append(binNumbers, bin)
	}
	sort.Slice(binNumbers, func(i, j int) bool {
		return binNumbers[i] < binNumbers[j]
	})
	if err := wr.writeInt32(int32(len(binNumbers))); err != nil {
		return fmt.Errorf("writing n_bin: %w", err)
	}
	for _, bin := range binNumbers {
		chunks := ref.Bins[bin]
		sorted := make([]Chunk, len(chunks))
		copy(sorted, chunks)
		sort.Slice(sorted, func(i, j int) bool {
			return sorted[i].Beg < sorted[j].Beg
		})
		if err := wr.writeInt32(int32(bin)); err != nil {
			return fmt.Errorf("writing bin %d number: %w", bin, err)
		}
		if err := wr.writeInt32(int32(len(sorted))); err != nil {
			return fmt.Errorf("writing bin %d n_chunk: %w", bin, err)
		}
		for _, c := range sorted {
			if err := wr.writeUint64(uint64(c.Beg)); err != nil {
				return fmt.Errorf("writing chunk beg: %w", err)
			}
			if err := wr.writeUint64(uint64(c.End)); err != nil {
				return fmt.Errorf("writing chunk end: %w", err)
			}
		}
	}
	if err := wr.writeInt32(int32(len(ref.Intervals))); err != nil {
		return fmt.Errorf("writing n_intv: %w", err)
	}
	for i, iv := range ref.Intervals {
		if err := wr.writeUint64(uint64(iv)); err != nil {
			return fmt.Errorf("writing interval[%d]: %w", i, err)
		}
	}
	if err := wr.writeUint64(ref.NoCoor); err != nil {
		return fmt.Errorf("writing n_no_coor: %w", err)
	}
	return nil
}
