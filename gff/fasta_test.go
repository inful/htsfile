package gff

import (
	"errors"
	"strings"
	"testing"
)

func TestEncodeFastaSection_EmitsTerminatorAndContent(t *testing.T) {
	recs := []*Record{
		{SeqID: "chr1", Source: "test", Type: "gene", Start: 100, End: 200},
	}
	var out strings.Builder
	enc := NewEncoder(&out)
	if err := enc.EncodeRecord(recs[0]); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	fasta := ">chr1\nACGTACGT\nACGT\n"
	if err := enc.EncodeFastaSection(strings.NewReader(fasta)); err != nil {
		t.Fatalf("EncodeFastaSection: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	got := out.String()
	// Note: phase column emits "0" by default (a valid
	// GFF3 value), so the per-record line is
	//   chr1 test gene 100 200 . . 0 .
	want := "chr1\ttest\tgene\t100\t200\t.\t.\t0\t.\n###\n>chr1\nACGTACGT\nACGT\n"
	if got != want {
		t.Errorf("output mismatch:\ngot:  %q\nwant: %q", got, want)
	}
}

func TestEncodeFastaSection_NilReader(t *testing.T) {
	var out strings.Builder
	enc := NewEncoder(&out)
	if err := enc.EncodeFastaSection(nil); err != nil {
		t.Fatalf("EncodeFastaSection(nil): %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if got := out.String(); got != "###\n" {
		t.Errorf("nil reader: got %q, want %q", got, "###\n")
	}
}

func TestEncodeFastaSection_DoubleCallRejected(t *testing.T) {
	var out strings.Builder
	enc := NewEncoder(&out)
	if err := enc.EncodeFastaSection(strings.NewReader(">x\nA\n")); err != nil {
		t.Fatalf("first call: %v", err)
	}
	err := enc.EncodeFastaSection(strings.NewReader(">y\nC\n"))
	if !errors.Is(err, errFastaAlreadyWritten) {
		t.Errorf("second call: got err %v, want errFastaAlreadyWritten", err)
	}
}

func TestEncodeFastaSection_CRLF(t *testing.T) {
	var out strings.Builder
	enc := NewEncoderCRLF(&out)
	if err := enc.EncodeFastaSection(strings.NewReader(">x\nA\n")); err != nil {
		t.Fatalf("EncodeFastaSection: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "###\r\n>x\nA\n"
	if got := out.String(); got != want {
		t.Errorf("CRLF: got %q, want %q", got, want)
	}
}