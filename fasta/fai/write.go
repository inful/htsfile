package fai

import (
	"errors"
	"fmt"
	"io"
)

// Writer serialises an [*Index] to an underlying [io.Writer]
// in the canonical .fai format (one tab-separated record per
// line, LF line endings, UTF-8).
//
// Typical usage:
//
//	idx, _ := fai.Build(f)
//	fai.NewWriter(out).Write(idx)
//
// A Writer is not safe for concurrent use.
type Writer struct {
	w io.Writer
}

// NewWriter wraps w in a [Writer] that emits .fai-formatted
// output.
func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w}
}

// Write emits idx in .fai format. The empty index is a no-op
// (no bytes written). A nil index is also a no-op.
//
// Output format per reference:
//
//	name<TAB>length<TAB>offset<TAB>lineBases<TAB>lineWidth<LF>
//
// Length, offset, lineBases, and lineWidth are written as
// base-10 integers. Names are written verbatim — no escaping
// or quoting. Tab and newline characters in Name are rejected
// because they would corrupt the line-oriented format.
func (w *Writer) Write(idx *Index) error {
	if idx == nil {
		return nil
	}
	for _, r := range idx.References {
		if r == nil {
			continue
		}
		if err := validateName(r.Name); err != nil {
			return err
		}
		if _, err := fmt.Fprintf(w.w, "%s\t%d\t%d\t%d\t%d\n",
			r.Name, r.Length, r.Offset, r.LineBases, r.LineWidth); err != nil {
			return err
		}
	}
	return nil
}

// validateName rejects any character that would break the
// tab/line-oriented .fai format.
func validateName(name string) error {
	for i := 0; i < len(name); i++ {
		c := name[i]
		if c == '\t' || c == '\n' || c == '\r' {
			return fmt.Errorf("fai: name %q contains disallowed byte 0x%02x at offset %d",
				name, c, i)
		}
	}
	return nil
}

// errMissingGT is returned by Build when the input does not
// start with a '>' header line.
var errMissingGT = errors.New("fai: FASTA input missing '>' header on first non-empty line")