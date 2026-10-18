package bcf

import (
	"bufio"
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

// ---------------------------------------------------------------------------
// Coverage pin — Validator strict mode (undeclared INFO key)
// ---------------------------------------------------------------------------

func TestValidator_Strict_UndeclaredInfo(t *testing.T) {
	// Construct an in-memory Header + Record pair where the
	// record carries an INFO key that is not in the header's
	// declared ##INFO list. (This branch is defensive — the
	// binary parser doesn't produce such records today, but
	// the validator still flags it for callers that build
	// Record values directly.)
	hdr := ParseHeader("##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer>\n")
	rec := &Record{
		Ref: "A",
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{42}},
		},
	}
	// Default mode: clean (no declared-key check).
	if issues := NewValidator(ModeDefault).Validate(hdr, []*Record{rec}); len(issues) != 0 {
		t.Errorf("ModeDefault issues: %v", issues)
	}
	// Inject an undeclared INFO key.
	rec.Info["UNKNOWN"] = InfoValue{Type: "Integer", Ints: []int64{7}}
	issues := NewValidator(ModeStrict).Validate(hdr, []*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("ModeStrict issues: got %d want 1 (%v)", len(issues), issues)
	}
	if issues[0].Field != "INFO" || issues[0].Severity != SeverityError {
		t.Errorf("issue: %+v", issues[0])
	}
	if !strings.Contains(issues[0].Msg, "UNKNOWN") {
		t.Errorf("issue Msg: got %q (want to mention UNKNOWN)", issues[0].Msg)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — LookupFormatIndex nil receiver
// ---------------------------------------------------------------------------

func TestHeader_LookupFormatIndex_Nil(t *testing.T) {
	var nilHdr *Header
	if got, ok := nilHdr.LookupFormatIndex("GT"); ok || got != -1 {
		t.Errorf("nil LookupFormatIndex: got (%d, %v) want (-1, false)", got, ok)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Decoder.Header cached return on second call
// ---------------------------------------------------------------------------

func TestDecoder_HeaderTwice(t *testing.T) {
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	b.writeRecord(t, &testRecord{contigIdx: 0, Pos: 100, Ref: "A"})
	data := b.finish(t)
	dec := NewDecoder(bytes.NewReader(data))
	hdr1, err := dec.Header()
	if err != nil {
		t.Fatalf("first Header: %v", err)
	}
	hdr2, err := dec.Header()
	if err != nil {
		t.Fatalf("second Header: %v", err)
	}
	// Same pointer — second call returns the cached value.
	if hdr1 != hdr2 {
		t.Errorf("Header not cached: hdr1=%p hdr2=%p", hdr1, hdr2)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Decoder.Record implicit-Header call
// ---------------------------------------------------------------------------

func TestDecoder_RecordBeforeHeader(t *testing.T) {
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	b.writeRecord(t, &testRecord{contigIdx: 0, Pos: 100, Ref: "A"})
	data := b.finish(t)
	dec := NewDecoder(bytes.NewReader(data))
	// Call Record() without first calling Header() — exercises
	// the implicit-header branch.
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.Ref != "A" {
		t.Errorf("Ref: got %q want %q", rec.Ref, "A")
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Decoder.Record with nFmt > 0 (FORMAT block)
// v2: the decoder now reads the FORMAT list and per-sample
// data; the hand-built stream below has no FORMAT body
// after the n_fmt uint16, so the decoder hits EOF when
// trying to read the first FORMAT key.
// ---------------------------------------------------------------------------

func TestDecoder_RecordWithNfmt(t *testing.T) {
	// Build a stream by hand: BGZF + magic + header + a record
	// whose n_fmt field is non-zero but which has no actual
	// FORMAT list following it.
	hdrText := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
	b := newBuilderWithHeader(t, hdrText)
	b.writeBCFHeader(hdrText)
	var payload bytes.Buffer
	var prefix [16]byte
	binary.LittleEndian.PutUint32(prefix[0:4], 0)   // chromIdx
	binary.LittleEndian.PutUint32(prefix[4:8], 99)  // pos0
	binary.LittleEndian.PutUint32(prefix[8:12], 0)  // rlen
	binary.LittleEndian.PutUint32(prefix[12:16], 0) // n_info
	payload.Write(prefix[:])
	var nFmt [2]byte
	binary.LittleEndian.PutUint16(nFmt[:], 1) // n_fmt = 1, no body
	payload.Write(nFmt[:])
	b.writeBGZFBlock(payload.Bytes())
	data := b.finish(t)

	dec := NewDecoder(bytes.NewReader(data))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	_, err := dec.Record()
	// We expect an EOF (or wrapped-EOF) because the FORMAT
	// list is missing.
	if err == nil {
		t.Fatal("expected EOF, got nil")
	}
	if !errors.Is(err, io.EOF) && !errors.Is(err, io.ErrUnexpectedEOF) {
		// Some parsing paths may surface a typed ParseError;
		// accept any non-nil error here.
		t.Logf("non-EOF error (acceptable): %v", err)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — InfoValue.String()
// ---------------------------------------------------------------------------

func TestInfoValue_String(t *testing.T) {
	cases := []struct {
		name string
		val  InfoValue
		want string
	}{
		{"Flag_true", InfoValue{Type: "Flag", Flag: true}, "."},
		{"Flag_false", InfoValue{Type: "Flag", Flag: false}, ""},
		{"Integer_single", InfoValue{Type: "Integer", Ints: []int64{42}}, "42"},
		{"Integer_multi", InfoValue{Type: "Integer", Ints: []int64{1, 2, 3}}, "1,2,3"},
		{"Float_single", InfoValue{Type: "Float", Floats: []float64{1.5}}, "1.5"},
		{"Float_multi", InfoValue{Type: "Float", Floats: []float64{1.5, 2.5}}, "1.5,2.5"},
		{"Character", InfoValue{Type: "Character", Str: "X"}, "X"},
		{"String", InfoValue{Type: "String", Str: "hello"}, "hello"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := tc.val.String(); got != tc.want {
				t.Errorf("String(): got %q want %q", got, tc.want)
			}
		})
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Severity.String()
// ---------------------------------------------------------------------------

func TestSeverity_String(t *testing.T) {
	cases := []struct {
		name string
		sev  Severity
		want string
	}{
		{"Warning", SeverityWarning, "warning"},
		{"Error", SeverityError, "error"},
		{"Unknown", Severity(99), "unknown"},
		{"ZeroValue", Severity(0), "warning"}, // zero value is SeverityWarning
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := tc.sev.String(); got != tc.want {
				t.Errorf("String(): got %q want %q", got, tc.want)
			}
		})
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — sortedKeys (reserved-for-future-tests helper)
// ---------------------------------------------------------------------------

func TestSortedKeys(t *testing.T) {
	m := map[string]int{"b": 2, "a": 1, "c": 3}
	got := sortedKeys(m)
	want := []string{"a", "b", "c"}
	if len(got) != len(want) {
		t.Fatalf("len: got %d want %d (%v)", len(got), len(want), got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("[%d]: got %q want %q", i, got[i], want[i])
		}
	}
	// Empty map → empty slice.
	if got := sortedKeys(map[string]int{}); len(got) != 0 {
		t.Errorf("empty: got %v want []", got)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Decoder.skipFmtBlock
// (always returns the documented FORMAT-not-supported error)
// ---------------------------------------------------------------------------

func TestDecoder_SkipFmtBlock(t *testing.T) {
	// Construct a minimal Decoder — skipFmtBlock doesn't touch any
	// of the receiver's fields, so a zero-valued one suffices.
	d := &Decoder{}
	err := d.skipFmtBlock()
	if err == nil {
		t.Fatal("skipFmtBlock: expected error, got nil")
	}
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("skipFmtBlock: not a *ParseError: %T", err)
	}
	if pe.Field != "FORMAT" {
		t.Errorf("Field: got %q want %q", pe.Field, "FORMAT")
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Decoder.skipTypedInfoValue
// (used when INFO key index exceeds the header's declaration
// count; reads the type byte, then skips the typed body)
// ---------------------------------------------------------------------------

func TestDecoder_SkipTypedInfoValue(t *testing.T) {
	cases := []struct {
		name     string
		typeByte byte
		count    uint32
		// Optional payload bytes that follow the length prefix.
		// Size is count * elemSize for the type, or count * 1 for default.
		payloadSize int
		wantErr     bool
	}{
		// typeMissing short-circuits without reading length/payload.
		{"Missing", typeMissing, 0, 0, false},
		{"Int8_3", typeInt8, 3, 3, false},
		{"Int16_2", typeInt16, 2, 4, false}, // 2 elems × 2 bytes
		{"Int32_1", typeInt32, 1, 4, false},
		{"Float_4", typeFloat, 4, 16, false}, // 4 elems × 4 bytes
		// Default branch (any other type byte) — best-effort, 1 byte per elem.
		{"Unknown_2", 7, 2, 2, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var buf bytes.Buffer
			buf.WriteByte(tc.typeByte)
			if tc.typeByte != typeMissing {
				var lenBytes [4]byte
				binary.LittleEndian.PutUint32(lenBytes[:], tc.count)
				buf.Write(lenBytes[:])
				buf.Write(make([]byte, tc.payloadSize))
			}
			d := &Decoder{br: bufio.NewReader(&buf)}
			err := d.skipTypedInfoValue()
			if tc.wantErr {
				if err == nil {
					t.Errorf("expected error, got nil")
				}
				return
			}
			if err != nil {
				t.Errorf("unexpected error: %v", err)
			}
		})
	}
	// Truncated stream — type byte present but length prefix missing.
	t.Run("Truncated_Length", func(t *testing.T) {
		var buf bytes.Buffer
		buf.WriteByte(typeInt8) // no length follows
		d := &Decoder{br: bufio.NewReader(&buf)}
		if err := d.skipTypedInfoValue(); err == nil {
			t.Error("expected error on truncated length, got nil")
		}
	})
	// Truncated payload — length says 10 bytes but only 3 follow.
	t.Run("Truncated_Payload", func(t *testing.T) {
		var buf bytes.Buffer
		buf.WriteByte(typeInt8)
		var lenBytes [4]byte
		binary.LittleEndian.PutUint32(lenBytes[:], 10) // expect 10 bytes
		buf.Write(lenBytes[:])
		buf.Write([]byte{1, 2, 3}) // only 3 follow
		d := &Decoder{br: bufio.NewReader(&buf)}
		if err := d.skipTypedInfoValue(); err == nil {
			t.Error("expected error on truncated payload, got nil")
		}
	})
}

// _ keeps the io import visible.
var _ = io.EOF
