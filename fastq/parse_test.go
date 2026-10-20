package fastq

import (
	"errors"
	"io"
	"strings"
	"testing"
)

func TestDecoder_BasicRecord(t *testing.T) {
	src := "@read1\nACGT\n+\n!!!!\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.Header != "read1" {
		t.Errorf("Header: got %q want %q", rec.Header, "read1")
	}
	if string(rec.Seq) != "ACGT" {
		t.Errorf("Seq: got %q want %q", rec.Seq, "ACGT")
	}
	if rec.PlusLine != "" {
		t.Errorf("PlusLine: got %q want empty", rec.PlusLine)
	}
	if string(rec.Qual) != "!!!!" {
		t.Errorf("Qual: got %q want %q", rec.Qual, "!!!!")
	}
	if rec.PhredOffset != 33 {
		t.Errorf("PhredOffset: got %d want 33", rec.PhredOffset)
	}
	if rec.Line != 1 {
		t.Errorf("Line: got %d want 1", rec.Line)
	}
}

func TestDecoder_TwoRecords(t *testing.T) {
	src := "@r1\nACGT\n+\n!!!!\n@r2 description here\nGGGG\n+ r2 description here\nJJJJ\n"
	dec := NewDecoder(strings.NewReader(src))
	r1, err := dec.Record()
	if err != nil {
		t.Fatalf("r1: %v", err)
	}
	if r1.Header != "r1" || string(r1.Seq) != "ACGT" {
		t.Errorf("r1: %+v", r1)
	}
	r2, err := dec.Record()
	if err != nil {
		t.Fatalf("r2: %v", err)
	}
	if r2.Header != "r2 description here" {
		t.Errorf("r2 Header: got %q", r2.Header)
	}
	if r2.PlusLine != "r2 description here" {
		t.Errorf("r2 PlusLine: got %q", r2.PlusLine)
	}
	if string(r2.Seq) != "GGGG" || string(r2.Qual) != "JJJJ" {
		t.Errorf("r2 Seq/Qual: %q / %q", r2.Seq, r2.Qual)
	}
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Errorf("expected EOF, got %v", err)
	}
}

func TestDecoder_BlankLinesBetween(t *testing.T) {
	// Blank lines between records should be tolerated; the
	// parser resyncs on the next '@'.
	src := "@r1\nAC\n+\nII\n\n\n@r2\nGT\n+\nJJ\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Record(); err != nil {
		t.Fatalf("r1: %v", err)
	}
	if _, err := dec.Record(); err != nil {
		t.Fatalf("r2: %v", err)
	}
}

func TestDecoder_EmptyQualities(t *testing.T) {
	// Sequence and quality both empty — allowed in some
	// pipelines (zero-length reads). Should not error.
	src := "@empty\n\n+\n\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if len(rec.Seq) != 0 || len(rec.Qual) != 0 {
		t.Errorf("expected empty Seq/Qual, got Seq=%q Qual=%q", rec.Seq, rec.Qual)
	}
}

func TestDecoder_AutoDetectPhred64(t *testing.T) {
	// Phred+64 data: quality bytes in the [64, 104] range.
	// "hhhh" decodes to Phred 40 in Phred+64 (104-64=40).
	src := "@r\nACGT\n+\nhhhh\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.PhredOffset != 64 {
		t.Errorf("PhredOffset: got %d want 64", rec.PhredOffset)
	}
}

func TestDecoder_SetPhredOffsetOverride(t *testing.T) {
	// Caller forces Phred+64 even though the data looks like
	// Phred+33 — SetPhredOffset disables auto-detect.
	src := "@r\nACGT\n+\n!!!!\n"
	dec := NewDecoder(strings.NewReader(src))
	dec.SetPhredOffset(64)
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.PhredOffset != 64 {
		t.Errorf("PhredOffset: got %d want 64", rec.PhredOffset)
	}
}

func TestDecoder_LineTracking(t *testing.T) {
	// Each record's Line field points at the @header line.
	src := "\n\n@r1\nAC\n+\nII\n\n@r2\nGT\n+\nJJ\n"
	dec := NewDecoder(strings.NewReader(src))
	r1, err := dec.Record()
	if err != nil {
		t.Fatalf("r1: %v", err)
	}
	if r1.Line != 3 {
		t.Errorf("r1.Line: got %d want 3", r1.Line)
	}
	r2, err := dec.Record()
	if err != nil {
		t.Fatalf("r2: %v", err)
	}
	if r2.Line != 8 {
		t.Errorf("r2.Line: got %d want 8", r2.Line)
	}
}

func TestDecoder_CRLFLineEndings(t *testing.T) {
	// CR+LF input should parse identically to LF.
	src := "@r1\r\nACGT\r\n+\r\n!!!!\r\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if string(rec.Seq) != "ACGT" {
		t.Errorf("Seq: got %q", rec.Seq)
	}
	if string(rec.Qual) != "!!!!" {
		t.Errorf("Qual: got %q", rec.Qual)
	}
}

func TestDecoder_MalformedHeader(t *testing.T) {
	// First non-empty line doesn't start with '@'.
	src := "ACGT\n+\n!!!!\n"
	dec := NewDecoder(strings.NewReader(src))
	_, err := dec.Record()
	if !IsParseError(err) {
		t.Errorf("expected ParseError, got %v", err)
	}
}

func TestDecoder_MalformedPlus(t *testing.T) {
	// Third line doesn't start with '+'.
	src := "@r\nACGT\nnotaplus\n!!!!\n"
	dec := NewDecoder(strings.NewReader(src))
	_, err := dec.Record()
	if !IsParseError(err) {
		t.Errorf("expected ParseError, got %v", err)
	}
}

func TestDecoder_TruncatedAfterPlus(t *testing.T) {
	// EOF before the qualities line.
	src := "@r\nACGT\n+\n"
	dec := NewDecoder(strings.NewReader(src))
	_, err := dec.Record()
	if !IsParseError(err) {
		t.Errorf("expected ParseError on truncation, got %v", err)
	}
}

func TestDecoder_TruncatedMidFile(t *testing.T) {
	// One valid record, then EOF in the middle of the next.
	src := "@r1\nAC\n+\nII\n@r2\nGT\n+\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Record(); err != nil {
		t.Fatalf("r1: %v", err)
	}
	_, err := dec.Record()
	if !IsParseError(err) {
		t.Errorf("expected ParseError on truncated r2, got %v", err)
	}
}

func TestDecoder_EOFOnEmpty(t *testing.T) {
	dec := NewDecoder(strings.NewReader(""))
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Errorf("expected EOF on empty input, got %v", err)
	}
}

func TestDecoder_EOFAfterFinalNewline(t *testing.T) {
	// Trailing newline after the last record's qualities.
	src := "@r\nAC\n+\nII\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Record(); err != nil {
		t.Fatalf("r: %v", err)
	}
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Errorf("expected EOF, got %v", err)
	}
}

func TestDecoder_DoSBoundedRead(t *testing.T) {
	// A single line longer than the configured max produces a
	// ParseError rather than an unbounded allocation.
	long := strings.Repeat("A", 1000) + "\n@r\nAC\n+\nII\n"
	dec := NewDecoderWithMaxLine(strings.NewReader(long), 100)
	_, err := dec.Record()
	if !IsParseError(err) {
		t.Errorf("expected ParseError on over-long line, got %v", err)
	}
}

func TestDecoder_LongSequence(t *testing.T) {
	// Sequences longer than the default scan buffer should still
	// work (DefaultMaxLineLength is 10 MB; a 1 MB sequence is
	// well within that).
	seq := strings.Repeat("A", 1<<20) // 1 MiB
	qual := strings.Repeat("I", 1<<20)
	src := "@big\n" + seq + "\n+\n" + qual + "\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if len(rec.Seq) != 1<<20 || len(rec.Qual) != 1<<20 {
		t.Errorf("len: Seq=%d Qual=%d", len(rec.Seq), len(rec.Qual))
	}
}
