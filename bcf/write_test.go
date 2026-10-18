package bcf

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

func TestEncoder_RoundTrip_Basic(t *testing.T) {
	hdr := ParseHeader(
		"##fileformat=VCFv4.2\n" +
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
			"##contig=<ID=chr1,length=1000>\n")
	records := []*Record{
		{
			CHROM: "chr1", Pos: 100,
			ID:   "rs1",
			Ref:  "A",
			Alt:  []string{"G", "C"},
			Qual: 99.5, QualKnown: true,
			Filter: []string{"PASS"},
			Info: map[string]InfoValue{
				"DP": {Type: "Integer", Ints: []int64{42}},
			},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, r := range records {
		if err := enc.WriteRecord(r); err != nil {
			t.Fatalf("WriteRecord: %v", err)
		}
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	dec := NewDecoder(&buf)
	gotHdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if gotHdr.FileFormat != hdr.FileFormat {
		t.Errorf("FileFormat: got %s want %s", gotHdr.FileFormat, hdr.FileFormat)
	}
	var got []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		got = append(got, rec)
	}
	if len(got) != len(records) {
		t.Fatalf("records: got %d want %d", len(got), len(records))
	}
	for i, want := range records {
		g := got[i]
		if g.CHROM != want.CHROM {
			t.Errorf("[%d] CHROM: got %q want %q", i, g.CHROM, want.CHROM)
		}
		if g.Pos != want.Pos {
			t.Errorf("[%d] Pos: got %d want %d", i, g.Pos, want.Pos)
		}
		if g.ID != want.ID {
			t.Errorf("[%d] ID: got %q want %q", i, g.ID, want.ID)
		}
		if g.Ref != want.Ref {
			t.Errorf("[%d] Ref: got %q want %q", i, g.Ref, want.Ref)
		}
		if !equalStrings(g.Alt, want.Alt) {
			t.Errorf("[%d] Alt: got %v want %v", i, g.Alt, want.Alt)
		}
		if g.QualKnown != want.QualKnown || g.Qual != want.Qual {
			t.Errorf("[%d] Qual: got %v/%v want %v/%v",
				i, g.Qual, g.QualKnown, want.Qual, want.QualKnown)
		}
		if !equalStrings(g.Filter, want.Filter) {
			t.Errorf("[%d] Filter: got %v want %v", i, g.Filter, want.Filter)
		}
		if v, ok := g.Info["DP"]; !ok || v.Type != "Integer" || len(v.Ints) != 1 || v.Ints[0] != 42 {
			t.Errorf("[%d] Info[DP]: ok=%v val=%+v", i, ok, v)
		}
	}
}

func TestEncoder_EmptyStream(t *testing.T) {
	hdr := ParseHeader("##fileformat=VCFv4.2\n")
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	dec := NewDecoder(&buf)
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Errorf("expected EOF, got %v", err)
	}
}

func TestEncoder_NoHeader(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	err := enc.WriteRecord(&Record{CHROM: "chr1"})
	if err == nil {
		t.Error("WriteRecord before WriteHeader: expected error, got nil")
	}
}

func TestEncoder_RecordAfterClose(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if err := enc.WriteRecord(&Record{CHROM: "chr1"}); err == nil {
		t.Error("WriteRecord after Close: expected error, got nil")
	}
}

func TestEncoder_NilHeader(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(nil); err == nil {
		t.Error("WriteHeader(nil): expected error, got nil")
	}
}

func TestEncoder_NilRecord(t *testing.T) {
	hdr := ParseHeader("##fileformat=VCFv4.2\n")
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(nil); err != nil {
		t.Errorf("WriteRecord(nil): expected no-op, got %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

// TestEncoder_NoContig verifies that a record referencing a
// contig not in the header produces an error.
func TestEncoder_NoContig(t *testing.T) {
	hdr := ParseHeader("##fileformat=VCFv4.2\n")
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	err := enc.WriteRecord(&Record{CHROM: "chr1", Pos: 1, Ref: "A", Alt: []string{"G"}})
	if err == nil {
		t.Error("WriteRecord with unknown contig: expected error, got nil")
	}
}

// equalStrings compares two []string for equality.
func equalStrings(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// ensure strings package is used
var _ = strings.HasPrefix
