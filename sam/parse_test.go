package sam

import (
	"errors"
	"io"
	"os"
	"strings"
	"testing"
)

// readAll reads through the decoder, returning (header, records, error).
func readAll(t *testing.T, path string) (*Header, []*Record) {
	t.Helper()
	f, err := os.Open(path) //nolint:gosec // test code, not user input
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close() //nolint:errcheck

	dec := NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	var recs []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, ioEOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		recs = append(recs, rec)
	}
	return hdr, recs
}

var ioEOF = io.EOF

func TestDecoder_Example(t *testing.T) {
	t.Parallel()
	hdr, recs := readAll(t, "testdata/example.sam")

	if hdr == nil {
		t.Fatal("Header = nil; want populated")
	}
	if hdr.HD == nil {
		t.Fatal("HD = nil; want populated")
	}
	if hdr.HD.VN != "1.6" {
		t.Errorf("HD.VN = %q; want 1.6", hdr.HD.VN)
	}
	if hdr.HD.SO != "coordinate" {
		t.Errorf("HD.SO = %q; want coordinate", hdr.HD.SO)
	}
	if len(hdr.SQs) != 2 {
		t.Errorf("len(SQs) = %d; want 2", len(hdr.SQs))
	}
	if hdr.SQs[0].SN != "chr1" || hdr.SQs[0].LN != 248956422 {
		t.Errorf("SQ[0] = %+v", hdr.SQs[0])
	}
	if len(hdr.RGs) != 2 {
		t.Errorf("len(RGs) = %d; want 2", len(hdr.RGs))
	}
	if hdr.RGs[0].ID != "RG1" || hdr.RGs[0].SM != "sample1" {
		t.Errorf("RG[0] = %+v", hdr.RGs[0])
	}
	if len(hdr.PGs) != 1 {
		t.Errorf("len(PGs) = %d; want 1", len(hdr.PGs))
	}
	if hdr.PGs[0].ID != "bwa" {
		t.Errorf("PG[0].ID = %q; want bwa", hdr.PGs[0].ID)
	}
	if len(hdr.COs) == 0 {
		t.Error("len(COs) = 0; want at least 1")
	}

	if len(recs) != 5 {
		t.Fatalf("len(recs) = %d; want 5", len(recs))
	}
	// r1: 0 flag, no mate, NM:i:1, MD:Z:99, RG:Z:RG1
	r0 := recs[0]
	if r0.QName != "r1" || r0.Flag != 0 || r0.RName != "chr1" {
		t.Errorf("r0 = %+v", r0)
	}
	if r0.Pos != 100 || r0.MapQ != 60 {
		t.Errorf("r0.Pos/MapQ = %d/%d", r0.Pos, r0.MapQ)
	}
	if len(r0.Tags) != 3 {
		t.Fatalf("r0.Tags = %v; want 3 tags", r0.Tags)
	}
	if r0.Tags[0].Key != "NM" || r0.Tags[0].Type != Tagi || r0.Tags[0].Value != "1" {
		t.Errorf("r0.Tags[0] = %+v", r0.Tags[0])
	}
	if r0.Tags[1].Key != "MD" || r0.Tags[1].Type != TagZ || r0.Tags[1].Value != "99" {
		t.Errorf("r0.Tags[1] = %+v", r0.Tags[1])
	}
	// r3: 83 flag (paired+proper_pair+first+reverse), RG2
	r2 := recs[2]
	if r2.Flag != 83 {
		t.Errorf("r2.Flag = %d; want 83", r2.Flag)
	}
	if r2.Flag&FlagPaired == 0 || r2.Flag&FlagProperPair == 0 ||
		r2.Flag&FlagFirst == 0 || r2.Flag&FlagReverse == 0 {
		t.Errorf("r2.Flag = %d; expected paired+proper_pair+first+reverse", r2.Flag)
	}
}

func TestDecoder_HeaderOnly(t *testing.T) {
	t.Parallel()
	src := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr == nil {
		t.Fatal("Header = nil")
	}
	if hdr.HD == nil || hdr.HD.VN != "1.6" {
		t.Errorf("HD = %+v", hdr.HD)
	}
	if len(hdr.SQs) != 1 || hdr.SQs[0].SN != "chr1" {
		t.Errorf("SQs = %+v", hdr.SQs)
	}
	// Reading a record should give io.EOF since there are no
	// records after the header.
	rec, err := dec.Record()
	if !errors.Is(err, ioEOF) {
		t.Errorf("Record err = %v; want io.EOF", err)
	}
	if rec != nil {
		t.Errorf("Record = %+v; want nil", rec)
	}
}

func TestDecoder_NoHeader(t *testing.T) {
	t.Parallel()
	// Records-only file (no @ lines).
	src := "r1\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGT\tIIII\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr != nil {
		t.Errorf("Header = %+v; want nil for records-only file", hdr)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.QName != "r1" {
		t.Errorf("QName = %q; want r1", rec.QName)
	}
}

func TestDecoder_TagWithColonInValue(t *testing.T) {
	t.Parallel()
	// Z tag values may contain colons (e.g. a SAM header reference).
	src := "r1\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGT\tIIII\tCO:Z:note:with:colons\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if len(rec.Tags) != 1 {
		t.Fatalf("len(Tags) = %d; want 1", len(rec.Tags))
	}
	tag := rec.Tags[0]
	if tag.Key != "CO" || tag.Type != TagZ {
		t.Errorf("tag = %+v", tag)
	}
	if tag.Value != "note:with:colons" {
		t.Errorf("Value = %q; want %q", tag.Value, "note:with:colons")
	}
}

func TestDecoder_MalformedHeaderTypeTag(t *testing.T) {
	t.Parallel()
	// Header line with only one type character: should error,
	// not be silently stored in COs.
	src := "@H\trandom\nr1\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGT\tIIII\n"
	dec := NewDecoder(strings.NewReader(src))
	_, err := dec.Header()
	if err == nil {
		t.Fatal("expected error for malformed header type tag '@H'")
	}
	if !IsParseError(err) {
		t.Errorf("error type = %T; want *ParseError", err)
	}
}

func TestDecoder_DuplicateHD(t *testing.T) {
	t.Parallel()
	// Two @HD lines: second should be rejected (spec allows
	// exactly one @HD). The first is preserved.
	src := "@HD\tVN:1.6\n@HD\tVN:1.0\nr1\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGT\tIIII\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err == nil {
		t.Fatal("expected error for duplicate @HD")
	}
	if hdr == nil || hdr.HD == nil {
		t.Fatal("expected first @HD to be preserved before error")
	}
	if hdr.HD.VN != "1.6" {
		t.Errorf("first @HD.VN = %q; want 1.6", hdr.HD.VN)
	}
}

func TestDecoder_EmptyFile(t *testing.T) {
	t.Parallel()
	dec := NewDecoder(strings.NewReader(""))
	hdr, err := dec.Header()
	if err != nil {
		t.Errorf("Header err = %v; want nil", err)
	}
	if hdr != nil {
		t.Errorf("Header = %+v; want nil on empty file", hdr)
	}
	rec, err := dec.Record()
	if !errors.Is(err, ioEOF) {
		t.Errorf("Record err = %v; want io.EOF", err)
	}
	if rec != nil {
		t.Errorf("Record = %+v; want nil", rec)
	}
}

func TestDecoder_TooFewFields(t *testing.T) {
	t.Parallel()
	src := "r1\t0\tchr1\n" // 3 columns; need 11
	dec := NewDecoder(strings.NewReader(src))
	_, _ = dec.Header()
	_, err := dec.Record()
	if err == nil {
		t.Fatal("expected error for line with < 11 fields")
	}
	if !IsParseError(err) {
		t.Errorf("error type = %T; want *ParseError", err)
	}
}

func TestDecoder_AutoCallsHeaderOnFirstRecord(t *testing.T) {
	t.Parallel()
	src := "@HD\tVN:1.6\nr1\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGT\tIIII\n"
	dec := NewDecoder(strings.NewReader(src))
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.QName != "r1" {
		t.Errorf("QName = %q", rec.QName)
	}
	hdr, err := dec.Header()
	if err != nil {
		t.Errorf("Header: %v", err)
	}
	if hdr == nil || hdr.HD == nil || hdr.HD.VN != "1.6" {
		t.Errorf("Header = %+v", hdr)
	}
}

func TestParseTagPairs(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in   string
		want map[string]string
	}{
		{
			"VN:1.6\tSO:coordinate",
			map[string]string{"VN": "1.6", "SO": "coordinate"},
		},
		{
			"",
			map[string]string{},
		},
		{
			"SN:chr1 LN:1000 M5:abc",
			map[string]string{"SN": "chr1", "LN": "1000", "M5": "abc"},
		},
	}
	for _, c := range cases {
		got := parseTagPairs(c.in)
		if len(got) != len(c.want) {
			t.Errorf("parseTagPairs(%q) = %v; want %v", c.in, got, c.want)
			continue
		}
		for k, v := range c.want {
			if got[k] != v {
				t.Errorf("parseTagPairs(%q)[%q] = %q; want %q", c.in, k, got[k], v)
			}
		}
	}
}

// TestDecoder_MaxLineLength pins the bounded-read behavior: a
// line longer than MaxLineLength must error, not silently allocate
// gigabytes.
func TestDecoder_MaxLineLength(t *testing.T) {
	t.Parallel()
	// 200 'A' characters in @SQ URL field, with a 100-byte limit.
	longURL := strings.Repeat("a", 200)
	src := "@HD\tVN:1.6\n" +
		"@SQ\tSN:chr1\tLN:1000\tUR:" + longURL + "\n"

	dec := NewDecoderWithMaxLine(strings.NewReader(src), 100)
	_, err := dec.Header()
	if err == nil {
		t.Fatal("expected error for header line exceeding MaxLineLength")
	}
	if !IsParseError(err) {
		t.Errorf("error type = %T; want *ParseError", err)
	}
	if !strings.Contains(err.Error(), "exceeds") {
		t.Errorf("error msg = %q; want it to mention 'exceeds'", err.Error())
	}
}
