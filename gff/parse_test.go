package gff

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

func TestHeader_Basic(t *testing.T) {
	src := "##gff-version 3\n" +
		"##sequence-region ctgA 1 1497228\n" +
		"##date 2024-01-15\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr.Version != 3 {
		t.Errorf("Version: %d", hdr.Version)
	}
	if len(hdr.SequenceRegions) != 1 || hdr.SequenceRegions[0].SeqID != "ctgA" {
		t.Errorf("SequenceRegions: %+v", hdr.SequenceRegions)
	}
	if hdr.Date != "2024-01-15" {
		t.Errorf("Date: %q", hdr.Date)
	}
	if len(hdr.MetaLines) < 3 {
		t.Errorf("MetaLines: %d", len(hdr.MetaLines))
	}
}

func TestHeader_NoPragmas(t *testing.T) {
	// File starts directly with a record.
	src := "ctgA\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene1\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if len(hdr.MetaLines) != 0 {
		t.Errorf("MetaLines: %d", len(hdr.MetaLines))
	}
}

func TestRecord_BasicFields(t *testing.T) {
	src := "##gff-version 3\n" +
		"ctgA\tEMBL\tgene\t1000\t9000\t.\t+\t.\tID=gene1;Name=ABC1\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.SeqID != "ctgA" {
		t.Errorf("SeqID: %q", rec.SeqID)
	}
	if rec.Source != "EMBL" {
		t.Errorf("Source: %q", rec.Source)
	}
	if rec.Type != "gene" {
		t.Errorf("Type: %q", rec.Type)
	}
	if rec.Start != 1000 || rec.End != 9000 {
		t.Errorf("Start/End: %d/%d", rec.Start, rec.End)
	}
	if rec.ScoreKnown {
		t.Errorf("ScoreKnown: got true (source had '.')")
	}
	if rec.Strand != '+' {
		t.Errorf("Strand: %q", string(rec.Strand))
	}
	if rec.Phase != PhaseUnknown {
		t.Errorf("Phase: %d", rec.Phase)
	}
	if rec.Attribute("ID") != "gene1" {
		t.Errorf("Attribute(ID): %q", rec.Attribute("ID"))
	}
	if rec.Attribute("Name") != "ABC1" {
		t.Errorf("Attribute(Name): %q", rec.Attribute("Name"))
	}
}

func TestRecord_MultivalueAttributes(t *testing.T) {
	src := "##gff-version 3\n" +
		"ctgA\t.\tmRNA\t1000\t9000\t.\t+\t.\tID=mrna1;Parent=gene1,transcript1\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatal(err)
	}
	parents := rec.AttributeList("Parent")
	if len(parents) != 2 {
		t.Fatalf("Parents: %d", len(parents))
	}
	if parents[0] != "gene1" || parents[1] != "transcript1" {
		t.Errorf("Parents: %v", parents)
	}
	if rec.Attribute("Parent") != "gene1" {
		t.Errorf("Attribute(Parent) first: %q", rec.Attribute("Parent"))
	}
}

func TestRecord_ScoreParsing(t *testing.T) {
	src := "##gff-version 3\n" +
		"ctgA\t.\tgene\t1000\t9000\t45.5\t+\t.\tID=g1\n" +
		"ctgA\t.\tgene\t2000\t3000\t.\t+\t.\tID=g2\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	r1, err := dec.Record()
	if err != nil {
		t.Fatal(err)
	}
	if !r1.ScoreKnown || r1.Score != 45.5 {
		t.Errorf("r1 Score: known=%v val=%v", r1.ScoreKnown, r1.Score)
	}
	r2, err := dec.Record()
	if err != nil {
		t.Fatal(err)
	}
	if r2.ScoreKnown {
		t.Errorf("r2 ScoreKnown: got true (source had '.')")
	}
}

func TestRecord_AllStrands(t *testing.T) {
	for _, want := range []byte{'+', '-', '.', '?'} {
		src := "##gff-version 3\n" +
			"ctgA\t.\tgene\t1\t100\t.\t" + string(want) + "\t.\tID=g\n"
		dec := NewDecoder(strings.NewReader(src))
		if _, err := dec.Header(); err != nil {
			t.Fatal(err)
		}
		rec, err := dec.Record()
		if err != nil {
			t.Fatalf("strand=%q: %v", string(want), err)
		}
		if rec.Strand != want {
			t.Errorf("strand: got %q want %q", string(rec.Strand), string(want))
		}
	}
}

func TestRecord_BadStrand(t *testing.T) {
	src := "##gff-version 3\nctgA\t.\tgene\t1\t100\t.\tX\t.\tID=g\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	if _, err := dec.Record(); err == nil {
		t.Fatal("expected strand error")
	}
}

func TestRecord_BadPhase(t *testing.T) {
	src := "##gff-version 3\nctgA\t.\tCDS\t1\t100\t.\t+\t9\tID=c\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	if _, err := dec.Record(); err == nil {
		t.Fatal("expected phase error")
	}
}

func TestRecord_BadStart(t *testing.T) {
	src := "##gff-version 3\nctgA\t.\tgene\t0\t100\t.\t+\t.\tID=g\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	if _, err := dec.Record(); err == nil {
		t.Fatal("expected start error (GFF3 is 1-based)")
	}
}

func TestRecord_EndLessThanStart(t *testing.T) {
	src := "##gff-version 3\nctgA\t.\tgene\t100\t50\t.\t+\t.\tID=g\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	if _, err := dec.Record(); err == nil {
		t.Fatal("expected end < start error")
	}
}

func TestRecord_WrongColumnCount(t *testing.T) {
	src := "##gff-version 3\nctgA\t.\tgene\t100\t9000\t.\t+\t.\n" // missing attributes
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	if _, err := dec.Record(); err == nil {
		t.Fatal("expected column count error")
	}
}

func TestRecord_Terminator(t *testing.T) {
	src := "##gff-version 3\n" +
		"ctgA\t.\tgene\t1000\t9000\t.\t+\t.\tID=g1\n" +
		"###\n" +
		"ctgA\t.\tgene\t2000\t3000\t.\t+\t.\tID=g2\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	r1, err := dec.Record()
	if err != nil || r1.Attribute("ID") != "g1" {
		t.Fatalf("first record: %v %v", r1, err)
	}
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Fatalf("after ###: want EOF, got %v", err)
	}
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Fatalf("after EOF: want EOF, got %v", err)
	}
}

func TestRecord_NonNumericStart(t *testing.T) {
	src := "##gff-version 3\nctgA\t.\tgene\tabc\t100\t.\t+\t.\tID=g\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	if _, err := dec.Record(); err == nil {
		t.Fatal("expected non-numeric start error")
	}
}

func TestRecord_MidFilePragmaIgnored(t *testing.T) {
	// Some producers emit tool-specific ## pragmas in the
	// data section; the decoder should skip them, not error.
	src := "##gff-version 3\n" +
		"##feature-ontology http://example.org/so.obo\n" +
		"ctgA\t.\tgene\t1000\t9000\t.\t+\t.\tID=g1\n" +
		"##source-ontology http://example.org/so.obo\n" +
		"ctgA\t.\texon\t1000\t2000\t.\t+\t.\tParent=g1\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	if _, err := dec.Record(); err != nil {
		t.Fatalf("record 1: %v", err)
	}
	r2, err := dec.Record()
	if err != nil {
		t.Fatalf("record 2 (after mid-file pragma): %v", err)
	}
	if r2.Attribute("ID") != "" || r2.Attribute("Parent") != "g1" {
		t.Errorf("r2 attrs: ID=%q Parent=%q", r2.Attribute("ID"), r2.Attribute("Parent"))
	}
}

func TestRoundTrip(t *testing.T) {
	src := "##gff-version 3\n" +
		"##sequence-region ctgA 1 1497228\n" +
		"ctgA\tEMBL\tgene\t1000\t9000\t45.5\t+\t.\tID=gene1;Name=ABC1\n" +
		"ctgA\tEMBL\texon\t1000\t2000\t.\t+\t.\tParent=gene1\n" +
		"ctgA\tEMBL\tCDS\t1200\t1500\t.\t+\t0\tParent=transcript1\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatal(err)
	}
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
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatal(err)
	}
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatal(err)
	}
	out := buf.String()
	// Sanity: contains version, all record IDs.
	for _, want := range []string{"##gff-version 3", "ID=gene1", "Parent=gene1", "ctgA\tEMBL"} {
		if !strings.Contains(out, want) {
			t.Errorf("round-trip missing %q\noutput:\n%s", want, out)
		}
	}
}

func TestValidator_CleanRecords(t *testing.T) {
	src := "##gff-version 3\n" +
		"ctgA\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene1\n" +
		"ctgA\t.\tmRNA\t1000\t9000\t.\t+\t.\tID=mrna1;Parent=gene1\n" +
		"ctgA\t.\texon\t1000\t2000\t.\t+\t.\tParent=mrna1\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatal(err)
	}
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
	val := NewValidator(ModeStrict)
	issues := val.Validate(hdr, recs)
	if len(issues) != 0 {
		t.Errorf("expected no issues, got %v", issues)
	}
}

func TestValidator_BadParentReference(t *testing.T) {
	src := "##gff-version 3\n" +
		"ctgA\t.\tmRNA\t1000\t9000\t.\t+\t.\tID=mrna1;Parent=gene_missing\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatal(err)
	}
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
	val := NewValidator(ModeStrict)
	issues := val.Validate(hdr, recs)
	if len(issues) == 0 {
		t.Fatal("expected parent reference error")
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
