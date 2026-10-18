package fasta

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

func TestRecord_BasicFields(t *testing.T) {
	src := ">seq1 E. coli K-12\nACGT\nACGT\n>seq2 some other desc\nGATTACA\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.ID != "seq1" {
		t.Errorf("ID: %q", rec.ID)
	}
	if rec.Description != "E. coli K-12" {
		t.Errorf("Description: %q", rec.Description)
	}
	if string(rec.Sequence) != "ACGTACGT" {
		t.Errorf("Sequence: %q", rec.Sequence)
	}
	if rec.Raw != "seq1 E. coli K-12" {
		t.Errorf("Raw: %q", rec.Raw)
	}
	rec2, err := dec.Record()
	if err != nil {
		t.Fatalf("Record 2: %v", err)
	}
	if rec2.ID != "seq2" || string(rec2.Sequence) != "GATTACA" {
		t.Errorf("Record 2: %+v", rec2)
	}
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Errorf("expected EOF, got %v", err)
	}
}

func TestRecord_NoDescription(t *testing.T) {
	src := ">seq1\nACGT\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatal(err)
	}
	if rec.ID != "seq1" || rec.Description != "" {
		t.Errorf("ID/Desc: %q / %q", rec.ID, rec.Description)
	}
}

func TestRecord_NoWhitespaceInHeader(t *testing.T) {
	// Header is just "seq1" with no whitespace.
	src := ">seq1\nACGT\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, _ := dec.Record()
	if rec.Description != "" {
		t.Errorf("Description: %q", rec.Description)
	}
}

func TestRecord_TabSeparatorInHeader(t *testing.T) {
	src := ">seq1\torganism=E.coli\nACGT\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, _ := dec.Record()
	if rec.ID != "seq1" {
		t.Errorf("ID: %q", rec.ID)
	}
	if rec.Description != "organism=E.coli" {
		t.Errorf("Description: %q", rec.Description)
	}
}

func TestRecord_SingleLineSequence(t *testing.T) {
	src := ">seq1\nACGTACGTACGTACGT\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, _ := dec.Record()
	if string(rec.Sequence) != "ACGTACGTACGTACGT" {
		t.Errorf("Sequence: %q", rec.Sequence)
	}
}

func TestRecord_EmptySequence(t *testing.T) {
	src := ">empty\n>next\nACGT\n"
	dec := NewDecoder(strings.NewReader(src))
	r1, _ := dec.Record()
	if len(r1.Sequence) != 0 {
		t.Errorf("empty sequence: got %q", r1.Sequence)
	}
	r2, _ := dec.Record()
	if r2.ID != "next" {
		t.Errorf("next ID: %q", r2.ID)
	}
}

func TestRecord_BlankLinesBetween(t *testing.T) {
	src := ">seq1\nACGT\n\nACGT\n\n>seq2\nGATTACA\n"
	dec := NewDecoder(strings.NewReader(src))
	r1, _ := dec.Record()
	if string(r1.Sequence) != "ACGTACGT" {
		t.Errorf("r1: %q", r1.Sequence)
	}
	r2, _ := dec.Record()
	if string(r2.Sequence) != "GATTACA" {
		t.Errorf("r2: %q", r2.Sequence)
	}
}

func TestRecord_MissingHeaderAtStart(t *testing.T) {
	dec := NewDecoder(strings.NewReader("ACGT\n>seq1\nGATTACA\n"))
	if _, err := dec.Record(); err == nil {
		t.Fatal("expected error on missing header at start")
	}
}

func TestRoundTrip(t *testing.T) {
	src := ">seq1 E. coli K-12\n" +
		"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
		"GATTACA\n"
	dec := NewDecoder(strings.NewReader(src))
	var recs []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		recs = append(recs, rec)
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	enc.SetLineWidth(40) // override default 80
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatal(err)
	}
	out := buf.String()
	// Verify the output has the right structure.
	for _, want := range []string{">seq1 E. coli K-12", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "GATTACA"} {
		if !strings.Contains(out, want) {
			t.Errorf("round-trip output missing %q\noutput:\n%s", want, out)
		}
	}
}

func TestRoundTrip_DefaultLineWidth(t *testing.T) {
	// 200-base sequence; default LineWidth=80 should produce
	// 80 + 80 + 40 wrapping.
	src := ">seq\n" + strings.Repeat("A", 200) + "\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, _ := dec.Record()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatal(err)
	}
	_ = enc.Flush()
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	// Expect: header + 80 + 80 + 40 = 4 lines.
	if len(lines) != 4 {
		t.Errorf("lines: %d (want 4)", len(lines))
	}
	if len(lines[1]) != 80 {
		t.Errorf("first seq line length: %d", len(lines[1]))
	}
}

func TestValidator_Default(t *testing.T) {
	recs := []*Record{
		{Line: 1, ID: "seq1", Sequence: []byte("ACGT")},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(recs)
	if len(issues) != 0 {
		t.Errorf("clean record issues: %v", issues)
	}
}

func TestValidator_MissingID(t *testing.T) {
	recs := []*Record{
		{Line: 1, Sequence: []byte("ACGT")},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(recs)
	if len(issues) == 0 {
		t.Fatal("expected missing-ID error")
	}
}

func TestValidator_EmptySequence(t *testing.T) {
	recs := []*Record{
		{Line: 1, ID: "seq1"},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(recs)
	if len(issues) == 0 {
		t.Fatal("expected empty-sequence error")
	}
}

func TestValidator_StrictAlphabet(t *testing.T) {
	recs := []*Record{
		{Line: 1, ID: "seq1", Sequence: []byte("ACGT")},
		{Line: 5, ID: "seq2", Sequence: []byte("ACGTX")}, // X is non-IUPAC
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(recs)
	if len(issues) == 0 {
		t.Fatal("expected strict alphabet error")
	}
}

func TestValidator_DuplicateID(t *testing.T) {
	recs := []*Record{
		{Line: 1, ID: "seq1", Sequence: []byte("ACGT")},
		{Line: 5, ID: "seq1", Sequence: []byte("GATT")},
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(recs)
	var dup bool
	for _, iss := range issues {
		if iss.Field == "ID" {
			dup = true
		}
	}
	if !dup {
		t.Errorf("expected duplicate-ID warning, got %v", issues)
	}
}

func TestParseError_Helper(t *testing.T) {
	e := &ParseError{Line: 5, Msg: "bad"}
	if !strings.Contains(e.Error(), "line 5") {
		t.Errorf("Error: %q", e.Error())
	}
	if !IsParseError(e) {
		t.Error("IsParseError(*ParseError)")
	}
	if IsParseError(errors.New("plain")) {
		t.Error("IsParseError(plain)")
	}
}
