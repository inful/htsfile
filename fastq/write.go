package fastq

import (
	"bufio"
	"io"
)

// Encoder writes a FASTQ document to an underlying [io.Writer].
//
// Typical usage:
//
//	enc := fastq.NewEncoder(w)
//	for _, r := range recs {
//	    _ = enc.EncodeRecord(r)
//	}
//	_ = enc.Flush()
//
// Each record is written as exactly four lines:
//
//	@<Header>
//	<Seq>
//	+[<PlusLine>]
//	<Qual>
//
// The Encoder does not modify Seq or Qual. The PhredOffset
// field is informational and is not emitted.
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w          *bufio.Writer
	lineEnding []byte
}

// NewEncoder wraps w in an [Encoder] that emits LF line
// endings.
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{
		w:          bufio.NewWriter(w),
		lineEnding: []byte{'\n'},
	}
}

// NewEncoderCRLF wraps w in an [Encoder] that emits CR+LF
// line endings.
func NewEncoderCRLF(w io.Writer) *Encoder {
	return &Encoder{
		w:          bufio.NewWriter(w),
		lineEnding: []byte{'\r', '\n'},
	}
}

// EncodeRecord writes one FASTQ record. A nil record is a
// no-op. Seq and Qual are emitted verbatim — callers must
// have already applied any case normalisation, whitespace
// stripping, or quality offset arithmetic they require.
func (e *Encoder) EncodeRecord(r *Record) error {
	if r == nil {
		return nil
	}
	// Header.
	if err := e.w.WriteByte('@'); err != nil {
		return err
	}
	if _, err := e.w.WriteString(r.Header); err != nil {
		return err
	}
	if _, err := e.w.Write(e.lineEnding); err != nil {
		return err
	}
	// Sequence.
	if _, err := e.w.Write(r.Seq); err != nil {
		return err
	}
	if _, err := e.w.Write(e.lineEnding); err != nil {
		return err
	}
	// Plus line.
	if err := e.w.WriteByte('+'); err != nil {
		return err
	}
	if r.PlusLine != "" {
		if _, err := e.w.WriteString(r.PlusLine); err != nil {
			return err
		}
	}
	if _, err := e.w.Write(e.lineEnding); err != nil {
		return err
	}
	// Qualities.
	if _, err := e.w.Write(r.Qual); err != nil {
		return err
	}
	_, err := e.w.Write(e.lineEnding)
	return err
}

// Flush writes any buffered data to the underlying writer.
func (e *Encoder) Flush() error { return e.w.Flush() }
