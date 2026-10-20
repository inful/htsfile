package fastq

import (
	"bytes"
	"strings"
	"testing"
)

func TestEncoder_BasicRecord(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	rec := &Record{
		Header:      "read1",
		Seq:         []byte("ACGT"),
		PlusLine:    "",
		Qual:        []byte("!!!!"),
		PhredOffset: 33,
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "@read1\nACGT\n+\n!!!!\n"
	if buf.String() != want {
		t.Errorf("got %q want %q", buf.String(), want)
	}
}

func TestEncoder_WithPlusLine(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	rec := &Record{
		Header:   "r1 desc",
		Seq:      []byte("ACGT"),
		PlusLine: "r1 desc",
		Qual:     []byte("!!!!"),
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "@r1 desc\nACGT\n+r1 desc\n!!!!\n"
	if buf.String() != want {
		t.Errorf("got %q want %q", buf.String(), want)
	}
}

func TestEncoder_NilRecord(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeRecord(nil); err != nil {
		t.Errorf("nil record should be no-op, got %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if buf.Len() != 0 {
		t.Errorf("expected empty output, got %q", buf.String())
	}
}

func TestEncoder_EmptySeqAndQual(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	rec := &Record{Header: "empty", Seq: []byte{}, Qual: []byte{}}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "@empty\n\n+\n\n"
	if buf.String() != want {
		t.Errorf("got %q want %q", buf.String(), want)
	}
}

func TestEncoder_CRLF(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoderCRLF(&buf)
	rec := &Record{Header: "r", Seq: []byte("AC"), Qual: []byte("II")}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "@r\r\nAC\r\n+\r\nII\r\n"
	if buf.String() != want {
		t.Errorf("got %q want %q", buf.String(), want)
	}
}

func TestEncoder_RoundTrip(t *testing.T) {
	// Write a few records, read them back, verify equality.
	// Use Phred+33 throughout so the Decoder's auto-detect
	// stays consistent across records (real FASTQ streams
	// never mix offsets; the round-trip is bytes-faithful
	// but the PhredOffset field is detected once per stream).
	originals := []*Record{
		{Header: "r1", Seq: []byte("ACGT"), Qual: []byte("!!!!"), PhredOffset: 33},
		{Header: "r2 desc", Seq: []byte("GGGG"), PlusLine: "r2 desc",
			Qual: []byte("JJJJ"), PhredOffset: 33},
		{Header: "r3", Seq: []byte("ACGTACGT"), Qual: []byte("87654321"), PhredOffset: 33},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	for _, r := range originals {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatalf("EncodeRecord: %v", err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	dec := NewDecoder(strings.NewReader(buf.String()))
	for i, want := range originals {
		got, err := dec.Record()
		if err != nil {
			t.Fatalf("record %d: %v", i, err)
		}
		if got.Header != want.Header {
			t.Errorf("rec %d Header: got %q want %q", i, got.Header, want.Header)
		}
		if string(got.Seq) != string(want.Seq) {
			t.Errorf("rec %d Seq: got %q want %q", i, got.Seq, want.Seq)
		}
		if got.PlusLine != want.PlusLine {
			t.Errorf("rec %d PlusLine: got %q want %q", i, got.PlusLine, want.PlusLine)
		}
		if string(got.Qual) != string(want.Qual) {
			t.Errorf("rec %d Qual: got %q want %q", i, got.Qual, want.Qual)
		}
		if got.PhredOffset != 33 {
			t.Errorf("rec %d PhredOffset: got %d want 33", i, got.PhredOffset)
		}
	}
	if _, err := dec.Record(); !isEOF(err) {
		t.Errorf("expected EOF after round-trip, got %v", err)
	}
}

// isEOF returns true if err is io.EOF or wraps it.
func isEOF(err error) bool {
	if err == nil {
		return false
	}
	return err.Error() == "EOF" ||
		strings.Contains(err.Error(), "EOF")
}

func TestEncoder_MultipleRecords(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	for i := 0; i < 3; i++ {
		rec := &Record{
			Header: "r" + string(rune('0'+i)),
			Seq:    []byte("ACGT"),
			Qual:   []byte("!!!!"),
		}
		if err := enc.EncodeRecord(rec); err != nil {
			t.Fatalf("EncodeRecord: %v", err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "@r0\nACGT\n+\n!!!!\n" +
		"@r1\nACGT\n+\n!!!!\n" +
		"@r2\nACGT\n+\n!!!!\n"
	if buf.String() != want {
		t.Errorf("got %q want %q", buf.String(), want)
	}
}
