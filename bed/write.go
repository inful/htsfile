package bed

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strconv"
)

// Encoder writes BED records to an underlying [io.Writer].
//
// Typical usage:
//
//	enc := bed.NewEncoder(w)
//	if track != nil { _ = enc.EncodeTrack(track) }
//	for _, r := range recs { _ = enc.EncodeRecord(r) }
//	_ = enc.Flush()
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w          *bufio.Writer
	lineEnding []byte
}

// NewEncoder wraps w in an [Encoder] with default LF line endings.
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{
		w:          bufio.NewWriter(w),
		lineEnding: []byte{'\n'},
	}
}

// NewEncoderCRLF wraps w in an [Encoder] that uses CR+LF line endings.
func NewEncoderCRLF(w io.Writer) *Encoder {
	return &Encoder{
		w:          bufio.NewWriter(w),
		lineEnding: []byte{'\r', '\n'},
	}
}

// EncodeTrack writes a "track ..." header line followed by the line
// ending. If t.Raw is set, that exact line is emitted verbatim;
// otherwise a line is synthesized from the structured fields. The
// Raw path is preferred for round-trip fidelity (it preserves
// unknown attributes this package does not model). EncodeTrack is
// optional — the encoder will work without it.
func (e *Encoder) EncodeTrack(t *Track) error {
	if t == nil {
		return nil
	}
	if t.Raw != "" {
		// Emit the original line verbatim. It already starts with
		// "track" by construction (the parser stores the full line).
		if _, err := e.w.WriteString(t.Raw); err != nil {
			return err
		}
		_, err := e.w.Write(e.lineEnding)
		return err
	}
	// Synthesize from the structured fields.
	var b bytes.Buffer
	b.WriteString("track")
	if t.Name != "" {
		fmt.Fprintf(&b, " name=%q", t.Name)
	}
	if t.Description != "" {
		fmt.Fprintf(&b, " description=%q", t.Description)
	}
	if t.UseScore {
		b.WriteString(" useScore=1")
	}
	if t.ItemRgbOn {
		b.WriteString(` itemRgb="On"`)
	}
	if t.Visibility != "" {
		fmt.Fprintf(&b, " visibility=%q", t.Visibility)
	}
	_, _ = b.Write(e.lineEnding)
	_, err := e.w.Write(b.Bytes())
	return err
}

// EncodeRecord writes one data record, emitting the number of
// columns implied by r.ColumnCount(). Per the BED spec, a file must
// have a consistent column count across all records; this method
// does not enforce that — the validator does.
func (e *Encoder) EncodeRecord(r *Record) error {
	var b bytes.Buffer
	b.WriteString(r.Chrom)
	b.WriteByte('\t')
	b.WriteString(strconv.FormatInt(r.ChromStart, 10))
	b.WriteByte('\t')
	b.WriteString(strconv.FormatInt(r.ChromEnd, 10))

	if r.Name != nil {
		b.WriteByte('\t')
		b.WriteString(*r.Name)
	}
	if r.Score != nil {
		b.WriteByte('\t')
		b.WriteString(strconv.Itoa(*r.Score))
	}
	if r.Strand != nil {
		b.WriteByte('\t')
		b.WriteString(r.Strand.String())
	}
	if r.ThickStart != nil {
		b.WriteByte('\t')
		b.WriteString(strconv.FormatInt(*r.ThickStart, 10))
	}
	if r.ThickEnd != nil {
		b.WriteByte('\t')
		b.WriteString(strconv.FormatInt(*r.ThickEnd, 10))
	}
	if r.ItemRgb != nil {
		b.WriteByte('\t')
		fmt.Fprintf(&b, "%d,%d,%d", r.ItemRgb.R, r.ItemRgb.G, r.ItemRgb.B)
	}
	if r.BlockCount != nil {
		b.WriteByte('\t')
		b.WriteString(strconv.Itoa(*r.BlockCount))
	}
	if r.BlockSizes != nil {
		b.WriteByte('\t')
		writeIntList(&b, r.BlockSizes)
	}
	if r.BlockStarts != nil {
		b.WriteByte('\t')
		writeIntList(&b, r.BlockStarts)
	}
	_, _ = b.Write(e.lineEnding)
	_, err := e.w.Write(b.Bytes())
	return err
}

// Flush writes any buffered data to the underlying writer.
func (e *Encoder) Flush() error { return e.w.Flush() }

// writeIntList writes a comma-separated list of integers to b.
func writeIntList(b *bytes.Buffer, vs []int64) {
	for i, n := range vs {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.FormatInt(n, 10))
	}
}
