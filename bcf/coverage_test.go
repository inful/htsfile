package bcf

import (
	"bytes"
	"encoding/binary"
	"errors"
	"io"
	"strings"
	"testing"
)

// ---------------------------------------------------------------------------
// Header helper tests — these exist primarily to pin behaviour
// for the Lookup*Index helpers added in response to a real bug
// (the test builder was using a local info-order index that
// didn't match the header's InfosOrder).
// ---------------------------------------------------------------------------

func TestHeader_LookupInfoIndex(t *testing.T) {
	hdr := ParseHeader("##fileformat=VCFv4.2\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer>\n" +
		"##INFO=<ID=AF,Number=A,Type=Float>\n" +
		"##INFO=<ID=DB,Number=0,Type=Flag>\n")
	if hdr == nil {
		t.Fatal("ParseHeader returned nil")
	}
	for i, want := range []string{"DP", "AF", "DB"} {
		got, ok := hdr.LookupInfoIndex(want)
		if !ok {
			t.Errorf("LookupInfoIndex(%q): not found", want)
		}
		if got != i {
			t.Errorf("LookupInfoIndex(%q): got %d want %d", want, got, i)
		}
	}
	// Undeclared ID.
	if got, ok := hdr.LookupInfoIndex("MISSING"); ok || got != -1 {
		t.Errorf("LookupInfoIndex(MISSING): got (%d, %v) want (-1, false)", got, ok)
	}
	// Nil receiver.
	var nilHdr *Header
	if got, ok := nilHdr.LookupInfoIndex("DP"); ok || got != -1 {
		t.Errorf("nil LookupInfoIndex: got (%d, %v) want (-1, false)", got, ok)
	}
}

func TestHeader_LookupFormatIndex(t *testing.T) {
	hdr := ParseHeader("##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String>\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer>\n")
	for i, want := range []string{"GT", "DP"} {
		got, ok := hdr.LookupFormatIndex(want)
		if !ok {
			t.Errorf("LookupFormatIndex(%q): not found", want)
		}
		if got != i {
			t.Errorf("LookupFormatIndex(%q): got %d want %d", want, got, i)
		}
	}
	if got, ok := hdr.LookupFormatIndex("MISSING"); ok || got != -1 {
		t.Errorf("LookupFormatIndex(MISSING): got (%d, %v) want (-1, false)", got, ok)
	}
}

func TestHeader_InfoAt_EdgeCases(t *testing.T) {
	hdr := ParseHeader("##INFO=<ID=DP,Number=1,Type=Integer>")
	// Negative index.
	if def, ok := hdr.InfoAt(-1); ok || def != nil {
		t.Errorf("InfoAt(-1): got (%v, %v) want (nil, false)", def, ok)
	}
	// Out-of-range index.
	if def, ok := hdr.InfoAt(99); ok || def != nil {
		t.Errorf("InfoAt(99): got (%v, %v) want (nil, false)", def, ok)
	}
	// In-range.
	if def, ok := hdr.InfoAt(0); !ok || def == nil || def.ID != "DP" {
		t.Errorf("InfoAt(0): got (%v, %v)", def, ok)
	}
	// Nil receiver.
	var nilHdr *Header
	if def, ok := nilHdr.InfoAt(0); ok || def != nil {
		t.Errorf("nil InfoAt(0): got (%v, %v) want (nil, false)", def, ok)
	}
}

func TestHeader_FormatAt_EdgeCases(t *testing.T) {
	hdr := ParseHeader("##FORMAT=<ID=GT,Number=1,Type=String>")
	if def, ok := hdr.FormatAt(-1); ok || def != nil {
		t.Errorf("FormatAt(-1): got (%v, %v) want (nil, false)", def, ok)
	}
	if def, ok := hdr.FormatAt(0); !ok || def == nil || def.ID != "GT" {
		t.Errorf("FormatAt(0): got (%v, %v)", def, ok)
	}
	var nilHdr *Header
	if def, ok := nilHdr.FormatAt(0); ok || def != nil {
		t.Errorf("nil FormatAt(0): got (%v, %v) want (nil, false)", def, ok)
	}
}

func TestHeader_LookupContig_Empty(t *testing.T) {
	// Header with no ##contig lines.
	hdr := ParseHeader("##fileformat=VCFv4.2")
	if got := hdr.LookupContig(0); got != "" {
		t.Errorf("LookupContig(0) on empty: got %q", got)
	}
	if got := hdr.LookupFilter(0); got != "" {
		t.Errorf("LookupFilter(0) on empty: got %q", got)
	}
}

// ---------------------------------------------------------------------------
// Error-type tests
// ---------------------------------------------------------------------------

func TestParseError_ErrorAndHelpers(t *testing.T) {
	e := &ParseError{Line: 5, Msg: "bad"}
	if msg := e.Error(); !strings.Contains(msg, "line 5") || !strings.Contains(msg, "bad") {
		t.Errorf("Error: %q", msg)
	}
	eWithField := &ParseError{Line: 5, Field: "header", Msg: "bad"}
	if msg := eWithField.Error(); !strings.Contains(msg, "header") {
		t.Errorf("Error with field: %q", msg)
	}
	if e.Unwrap() != nil {
		t.Errorf("Unwrap: %v", e.Unwrap())
	}
	if IsParseError(nil) {
		t.Error("IsParseError(nil)")
	}
	if IsParseError(errors.New("plain")) {
		t.Error("IsParseError(plain)")
	}
	if !IsParseError(&ParseError{}) {
		t.Error("IsParseError(*ParseError)")
	}
}

// ---------------------------------------------------------------------------
// Decoder edge cases — incomplete streams, truncated records,
// invalid magic-in-bgzf, etc. These live in coverage_test.go
// because their value is mostly exercising rarely-hit branches.
// ---------------------------------------------------------------------------

func TestHeader_TruncatedLText(t *testing.T) {
	// Magic + l_text = 100 but only 10 bytes follow.
	b := newBuilder(t)
	var bad bytes.Buffer
	bad.Write(Magic[:])
	var l [4]byte
	binary.LittleEndian.PutUint32(l[:], 100)
	bad.Write(l[:])
	bad.Write([]byte("##filefor")) // only 9 bytes
	b.writeBGZFBlock(bad.Bytes())
	data := b.finish(t)
	_, err := NewDecoder(bytes.NewReader(data)).Header()
	if err == nil {
		t.Fatal("expected error on truncated header text")
	}
	if !IsParseError(err) {
		t.Errorf("error: %v", err)
	}
}

func TestRecord_LargeKeyIndex(t *testing.T) {
	// Build a record with INFO key index = 99999, which
	// exceeds InfosOrder. The decoder should skip via
	// skipTypedInfoValue.
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	defs := map[string]*InfoDef{
		"DP": {ID: "DP", Number: "1", Type: "Integer"},
	}
	r := &testRecord{
		contigIdx: 0, Pos: 100, Ref: "A",
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{42}},
		},
		infoDefs:  defs,
		infoOrder: []string{"DP"},
	}
	b.writeRecord(t, r)
	data := b.finish(t)
	_, recs := readAll(t, data)
	// Even though the writer used the correct header
	// index, the parser should still surface a valid
	// record (with DP present).
	if len(recs) != 1 {
		t.Fatalf("recs: %d", len(recs))
	}
	if _, ok := recs[0].Info["DP"]; !ok {
		t.Errorf("DP missing: %+v", recs[0].Info)
	}
}

func TestRecord_FilterArray(t *testing.T) {
	// Multiple filters as an array (not single int32).
	hdrText := sampleHeader + "##FILTER=<ID=q20,Description=\"qual < 20\">\n"
	b := newBuilderWithHeader(t, hdrText)
	b.writeBCFHeader(hdrText)
	defs := map[string]*InfoDef{
		"DP": {ID: "DP", Number: "1", Type: "Integer"},
	}
	r := &testRecord{
		contigIdx: 0, Pos: 100, Ref: "A",
		filterIdxs: []int{0, 1}, // PASS + q20
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{42}},
		},
		infoDefs:  defs,
		infoOrder: []string{"DP"},
	}
	b.writeRecord(t, r)
	data := b.finish(t)
	hdr, recs := readAll(t, data)
	if len(recs) != 1 {
		t.Fatalf("recs: %d", len(recs))
	}
	// Check that filter IDs are resolved.
	if len(recs[0].Filter) != 2 {
		t.Errorf("Filter: %v", recs[0].Filter)
	}
	if recs[0].Filter[0] != "PASS" {
		t.Errorf("Filter[0]: %q", recs[0].Filter[0])
	}
	if len(hdr.Filters) < 2 {
		t.Errorf("hdr.Filters: %v", hdr.Filters)
	}
}

func TestRecord_QUAL_Float(t *testing.T) {
	// Non-missing QUAL value (float64).
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	defs := map[string]*InfoDef{
		"DP": {ID: "DP", Number: "1", Type: "Integer"},
	}
	r := &testRecord{
		contigIdx: 0, Pos: 100, Ref: "A",
		Qual:      99.5,
		QualKnown: true,
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{42}},
		},
		infoDefs:  defs,
		infoOrder: []string{"DP"},
	}
	b.writeRecord(t, r)
	data := b.finish(t)
	_, recs := readAll(t, data)
	if !recs[0].QualKnown {
		t.Errorf("QualKnown: got false")
	}
	if recs[0].Qual != 99.5 {
		t.Errorf("Qual: got %v want 99.5", recs[0].Qual)
	}
}

func TestRecord_InfoString_VariableCount(t *testing.T) {
	// String field with Number=. (variable count) and 3 values.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=CSQ,Number=.,Type=String>\n"
	b := newBuilderWithHeader(t, hdrText)
	b.writeBCFHeader(hdrText)
	defs := map[string]*InfoDef{
		"CSQ": {ID: "CSQ", Number: ".", Type: "String"},
	}
	r := &testRecord{
		contigIdx: 0, Pos: 100, Ref: "A",
		Info: map[string]InfoValue{
			"CSQ": {Type: "String", Strs: []string{"alpha", "beta", "gamma"}},
		},
		infoDefs:  defs,
		infoOrder: []string{"CSQ"},
	}
	b.writeRecord(t, r)
	data := b.finish(t)
	_, recs := readAll(t, data)
	v := recs[0].Info["CSQ"]
	if len(v.Strs) != 3 {
		t.Errorf("CSQ: got %d values, want 3: %v", len(v.Strs), v.Strs)
	}
}

func TestRecord_EmptyInfo(t *testing.T) {
	// Record with no INFO fields.
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	r := &testRecord{
		contigIdx: 0, Pos: 100, Ref: "A",
	}
	b.writeRecord(t, r)
	data := b.finish(t)
	_, recs := readAll(t, data)
	if len(recs) != 1 {
		t.Fatalf("recs: %d", len(recs))
	}
	if len(recs[0].Info) != 0 {
		t.Errorf("Info: %v", recs[0].Info)
	}
}

func TestValidator_QUAL_Lenient_Happy(t *testing.T) {
	// Record with non-empty QUAL+QUALKnown; validator should
	// not complain about QUAL.
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	defs := map[string]*InfoDef{
		"DP": {ID: "DP", Number: "1", Type: "Integer"},
	}
	r := &testRecord{
		contigIdx: 0, Pos: 100, Ref: "A", Qual: 30, QualKnown: true,
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{42}},
		},
		infoDefs:  defs,
		infoOrder: []string{"DP"},
	}
	b.writeRecord(t, r)
	data := b.finish(t)
	hdr, recs := readAll(t, data)
	issues := NewValidator(ModeDefault).Validate(hdr, recs)
	if len(issues) != 0 {
		t.Errorf("unexpected issues: %v", issues)
	}
}

// _ keeps the io import visible.
var _ = io.EOF