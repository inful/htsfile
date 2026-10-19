package fastq

import (
	"strings"
	"testing"
)

func TestPhredEncoding_String(t *testing.T) {
	if got := PhredEncodingPhred.String(); got != "phred" {
		t.Errorf("Phred: %q", got)
	}
	if got := PhredEncodingSolexa.String(); got != "solexa" {
		t.Errorf("Solexa: %q", got)
	}
	if got := PhredEncoding(99).String(); got != "phred" {
		t.Errorf("invalid: %q", got)
	}
}

func TestPhredEncoding_IsValid(t *testing.T) {
	if !PhredEncodingPhred.IsValid() {
		t.Error("Phred should be valid")
	}
	if !PhredEncodingSolexa.IsValid() {
		t.Error("Solexa should be valid")
	}
	if PhredEncoding(99).IsValid() {
		t.Error("99 should be invalid")
	}
}

func TestDecoder_DefaultEncodingIsPhred(t *testing.T) {
	dec := NewDecoder(strings.NewReader(""))
	if got := dec.PhredEncoding(); got != PhredEncodingPhred {
		t.Errorf("default encoding: got %d, want Phred", got)
	}
}

func TestDecoder_SetSolexaEncoding(t *testing.T) {
	// Solexa+64 has the same byte range as Phred+64; the
	// test verifies the flag propagates from Decoder to
	// each parsed Record.
	src := "@r1\nACGT\n+\nIIII\n"
	dec := NewDecoder(strings.NewReader(src))
	dec.SetPhredOffset(64)
	dec.SetPhredEncoding(PhredEncodingSolexa)
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.PhredEncoding != PhredEncodingSolexa {
		t.Errorf("record encoding: got %d, want Solexa", rec.PhredEncoding)
	}
	if rec.PhredOffset != 64 {
		t.Errorf("record offset: got %d, want 64", rec.PhredOffset)
	}
}

func TestRecord_EncodingDefaultIsPhred(t *testing.T) {
	src := "@r1\nACGT\n+\nIIII\n"
	rec, err := NewDecoder(strings.NewReader(src)).Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.PhredEncoding != PhredEncodingPhred {
		t.Errorf("default record encoding: got %d, want Phred", rec.PhredEncoding)
	}
}
