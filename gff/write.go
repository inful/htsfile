package gff

import (
	"bufio"
	"bytes"
	"io"
	"strconv"
)

// Encoder writes a GFF3 document to an underlying [io.Writer].
//
// Typical usage:
//
//	enc := gff.NewEncoder(w)
//	if hdr != nil { _ = enc.EncodeHeader(hdr) }
//	for _, r := range recs { _ = enc.EncodeRecord(r) }
//	_ = enc.Flush()
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w            *bufio.Writer
	lineEnding   []byte
	fastaWritten bool
}

// NewEncoder wraps w in an [Encoder] with default LF line endings.
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{w: bufio.NewWriter(w), lineEnding: []byte{'\n'}}
}

// NewEncoderCRLF wraps w in an [Encoder] that uses CR+LF line endings.
func NewEncoderCRLF(w io.Writer) *Encoder {
	return &Encoder{w: bufio.NewWriter(w), lineEnding: []byte{'\r', '\n'}}
}

// EncodeHeader writes the GFF3 header. The output is the verbatim
// sequence of `##` pragma lines from hdr.MetaLines (in source
// order), preserving unknown pragmas, followed by the bare `#`
// comments (if any). Each call to EncodeHeader writes the
// header once; subsequent calls are no-ops.
func (e *Encoder) EncodeHeader(h *Header) error {
	if h == nil {
		return nil
	}
	for _, m := range h.MetaLines {
		if _, err := e.w.WriteString("##"); err != nil {
			return err
		}
		if _, err := e.w.WriteString(m); err != nil {
			return err
		}
		if _, err := e.w.Write(e.lineEnding); err != nil {
			return err
		}
	}
	for _, c := range h.Comments {
		if _, err := e.w.WriteString(c); err != nil {
			return err
		}
		if _, err := e.w.Write(e.lineEnding); err != nil {
			return err
		}
	}
	return nil
}

// EncodeRecord writes one feature record as a tab-separated
// line followed by the line ending.
func (e *Encoder) EncodeRecord(r *Record) error {
	if r == nil {
		return nil
	}
	var b bytes.Buffer
	b.WriteString(orDot(r.SeqID))
	b.WriteByte('\t')
	b.WriteString(orDot(r.Source))
	b.WriteByte('\t')
	b.WriteString(orDot(r.Type))
	b.WriteByte('\t')
	b.WriteString(strconv.FormatInt(r.Start, 10))
	b.WriteByte('\t')
	b.WriteString(strconv.FormatInt(r.End, 10))
	b.WriteByte('\t')
	b.WriteString(formatScore(r.Score, r.ScoreKnown))
	b.WriteByte('\t')
	b.WriteString(strandByte(r.Strand))
	b.WriteByte('\t')
	b.WriteString(phaseByte(r.Phase))
	b.WriteByte('\t')
	b.WriteString(formatAttributes(r.Attributes))
	_, _ = b.Write(e.lineEnding)
	_, err := e.w.Write(b.Bytes())
	return err
}

// strandByte returns the strand character as a one-rune string,
// or "." if the byte is zero (unset).
func strandByte(s byte) string {
	if s == 0 {
		return "."
	}
	return string(s)
}

// phaseByte returns the phase as a one- or two-character
// string, or "." if Phase is PhaseUnknown (-1).
func phaseByte(p int) string {
	if p < 0 {
		return "."
	}
	return strconv.Itoa(p)
}

// Flush writes any buffered data to the underlying writer.
func (e *Encoder) Flush() error { return e.w.Flush() }

// orDot returns s if non-empty, else ".". The GFF3 spec uses
// "." as the canonical sentinel for "field unavailable" across
// the seqid, source, and type columns.
func orDot(s string) string {
	if s == "" {
		return "."
	}
	return s
}
