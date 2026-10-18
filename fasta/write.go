package fasta

import (
	"bufio"
	"bytes"
	"io"
)

// Encoder writes a FASTA document to an underlying [io.Writer].
//
// Typical usage:
//
//	enc := fasta.NewEncoder(w)
//	for _, r := range recs { _ = enc.EncodeRecord(r) }
//	_ = enc.Flush()
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w          *bufio.Writer
	lineEnding []byte
	lineWidth  int
}

// DefaultLineWidth is the per-line sequence byte cap used by
// [NewEncoder]. 80 chars is the de facto FASTA convention.
const DefaultLineWidth = 80

// NewEncoder wraps w in an [Encoder] with default LF line
// endings and [DefaultLineWidth] sequence-line wrapping.
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{
		w:          bufio.NewWriter(w),
		lineEnding: []byte{'\n'},
		lineWidth:  DefaultLineWidth,
	}
}

// NewEncoderCRLF wraps w in an [Encoder] that uses CR+LF line
// endings.
func NewEncoderCRLF(w io.Writer) *Encoder {
	return &Encoder{
		w:          bufio.NewWriter(w),
		lineEnding: []byte{'\r', '\n'},
		lineWidth:  DefaultLineWidth,
	}
}

// SetLineWidth overrides the per-line sequence byte cap. Pass
// <= 0 to keep the existing width. The wrap is applied to
// Sequence bytes verbatim (no whitespace insertion — pure
// byte-count line breaks), so nucleotide and protein sequences
// line up the same way regardless of alphabet.
func (e *Encoder) SetLineWidth(n int) {
	if n > 0 {
		e.lineWidth = n
	}
}

// EncodeRecord writes one FASTA record. If r.Raw is set it
// is emitted verbatim; otherwise the header is synthesised as
// `>ID Description` (with a single space between ID and
// Description when Description is non-empty). Sequence bytes
// are emitted with lineWidth-byte line breaks.
func (e *Encoder) EncodeRecord(r *Record) error {
	if r == nil {
		return nil
	}
	// Header
	if err := e.w.WriteByte('>'); err != nil {
		return err
	}
	if r.Raw != "" {
		if _, err := e.w.WriteString(r.Raw); err != nil {
			return err
		}
	} else {
		if _, err := e.w.WriteString(r.ID); err != nil {
			return err
		}
		if r.Description != "" {
			if err := e.w.WriteByte(' '); err != nil {
				return err
			}
			if _, err := e.w.WriteString(r.Description); err != nil {
				return err
			}
		}
	}
	if _, err := e.w.Write(e.lineEnding); err != nil {
		return err
	}
	// Sequence with lineWidth-byte wrapping.
	if err := e.writeWrapped(r.Sequence); err != nil {
		return err
	}
	// Trailing line ending after the last sequence line.
	_, err := e.w.Write(e.lineEnding)
	return err
}

func (e *Encoder) writeWrapped(seq []byte) error {
	if len(seq) == 0 {
		return nil
	}
	w := e.lineWidth
	var b bytes.Buffer
	for i := 0; i < len(seq); i += w {
		end := i + w
		if end > len(seq) {
			end = len(seq)
		}
		b.Reset()
		b.Write(seq[i:end])
		b.Write(e.lineEnding)
		if _, err := e.w.Write(b.Bytes()); err != nil {
			return err
		}
	}
	return nil
}

// Flush writes any buffered data to the underlying writer.
func (e *Encoder) Flush() error { return e.w.Flush() }