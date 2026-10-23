package bam

import (
	"bytes"
	"encoding/binary"
	"errors"
	"strings"
	"testing"
)

// ---------------------------------------------------------------------------
// Header / Lookup tests
// ---------------------------------------------------------------------------

func TestHeader_Lookups(t *testing.T) {
	hdr := &Header{
		SQs: []*SQ{{Name: "chr1", Len: 1000}, {Name: "chr2", Len: 2000}},
		RGs: []*RG{{ID: "rg1"}, {ID: "rg2"}},
		PGs: []*PG{{ID: "prog1"}, {ID: "prog2"}},
	}
	if sq, ok := hdr.LookupSQ("chr1"); !ok || sq.Name != "chr1" {
		t.Errorf("LookupSQ(chr1): %v %v", sq, ok)
	}
	if sq, ok := hdr.LookupSQ("chrX"); ok || sq != nil {
		t.Errorf("LookupSQ(chrX): unexpected %v %v", sq, ok)
	}
	if rg, ok := hdr.LookupRG("rg2"); !ok || rg.ID != "rg2" {
		t.Errorf("LookupRG(rg2): %v %v", rg, ok)
	}
	if pg, ok := hdr.LookupPG("prog1"); !ok || pg.ID != "prog1" {
		t.Errorf("LookupPG(prog1): %v %v", pg, ok)
	}
	// Nil receiver should not panic.
	var nilHdr *Header
	if _, ok := nilHdr.LookupSQ("chr1"); ok {
		t.Error("nil LookupSQ returned ok=true")
	}
}

func TestSQ_String(t *testing.T) {
	if (&SQ{Name: "chr1", Len: 1000}).String() != "chr1:1000" {
		t.Errorf("SQ.String: got %q", (&SQ{Name: "chr1", Len: 1000}).String())
	}
	var nilSq *SQ
	if nilSq.String() != "" {
		t.Errorf("nil SQ.String: got %q", nilSq.String())
	}
}

// ---------------------------------------------------------------------------
// Error type tests
// ---------------------------------------------------------------------------

func TestParseError_Error(t *testing.T) {
	e := &ParseError{Line: 5, Msg: "bad"}
	if got := e.Error(); !strings.Contains(got, "line 5") || !strings.Contains(got, "bad") {
		t.Errorf("Error: %q", got)
	}
	e2 := &ParseError{Line: 5, Field: "header", Msg: "bad"}
	if got := e2.Error(); !strings.Contains(got, "header") {
		t.Errorf("Error with field: %q", got)
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
// Tag type tests
// ---------------------------------------------------------------------------

func TestTagType_StringAndParse(t *testing.T) {
	for _, tt := range []TagType{TagA, Tagc, TagC, Tags, TagS, Tagi, TagI, Tagf, Tagd, TagZ, TagH, TagB} {
		if !tt.IsValid() {
			t.Errorf("TagType %q should be valid", string(byte(tt)))
		}
		if tt.String() != string(byte(tt)) {
			t.Errorf("TagType.String: %q != %q", tt.String(), string(byte(tt)))
		}
		parsed, err := ParseTagType(byte(tt))
		if err != nil {
			t.Errorf("ParseTagType(%q): %v", string(byte(tt)), err)
		}
		if parsed != tt {
			t.Errorf("ParseTagType round-trip: %q != %q", string(byte(parsed)), string(byte(tt)))
		}
	}
	if _, err := ParseTagType('z'); err == nil {
		t.Error("expected error on unknown tag type")
	}
}

// ---------------------------------------------------------------------------
// Validator edge cases — strict-mode branches, QUAL length
// mismatch, and the strict-mode POS check. These live in
// coverage_test.go because they exist primarily to pin the
// otherwise-untested branches in checkOne; the round-trip
// path through checkOne is otherwise only hit by the
// happy-path TestValidator_CleanRecords /
// TestValidator_BadRefID tests in parse_test.go.
// ---------------------------------------------------------------------------

func TestValidator_BadQUALLength(t *testing.T) {
	// The BAM parser enforces SEQ/QUAL length agreement via
	// l_seq, so this test constructs the mismatch directly
	// on the exported Record struct (bypassing the binary
	// path) to exercise the validator branch.
	rec := &Record{
		Line:  1,
		Seq:   "ACGTA",
		Qual:  []byte("IIII"), // 4 vs 5
		Cigar: []uint32{cigar(CigarM, 5)},
	}
	issues := NewValidator(ModeDefault).Validate(nil, []*Record{rec})
	if len(issues) == 0 {
		t.Fatal("expected QUAL-length error")
	}
	var found bool
	for _, iss := range issues {
		if iss.Field == "QUAL" {
			found = true
		}
	}
	if !found {
		t.Errorf("no QUAL issue found: %v", issues)
	}
}

func TestValidator_StrictModeBadPOS(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000}}); err != nil {
		t.Fatal(err)
	}
	// POS = 0 with flag bit 0x4 (unmapped) clear — strict mode rejects.
	if err := b.writeRecord(&Record{
		QName: "bad", RefID: 0, Pos: 0, MapQ: 60, Bin: 4681,
		Flag:  0, // not unmapped
		Cigar: []uint32{cigar(CigarM, 5)},
		Seq:   "ACGTA", Qual: []byte("IIIII"),
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	hdr, recs := readAll(t, data)
	issues := NewValidator(ModeStrict).Validate(hdr, recs)
	if len(issues) == 0 {
		t.Fatal("expected strict-mode POS error")
	}
	var found bool
	for _, iss := range issues {
		if iss.Field == "POS" {
			found = true
		}
	}
	if !found {
		t.Errorf("no POS issue found: %v", issues)
	}
}

func TestValidator_StrictModeDuplicateTag(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000}}); err != nil {
		t.Fatal(err)
	}
	// Two tags with key NM — strict mode warns on duplicates.
	if err := b.writeRecord(&Record{
		QName: "dup", RefID: 0, Pos: 100, MapQ: 60, Bin: 4681,
		Cigar: []uint32{cigar(CigarM, 5)},
		Seq:   "ACGTA", Qual: []byte("IIIII"),
		Tags: []Tag{
			{Key: "NM", Type: Tagi, Value: "1"},
			{Key: "NM", Type: Tagi, Value: "2"},
		},
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	hdr, recs := readAll(t, data)
	issues := NewValidator(ModeStrict).Validate(hdr, recs)
	var found bool
	for _, iss := range issues {
		if iss.Field == "tag" {
			found = true
		}
	}
	if !found {
		t.Errorf("expected duplicate-tag warning, got %v", issues)
	}
}

func TestValidator_StrictModeBadNextRefID(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000}}); err != nil {
		t.Fatal(err)
	}
	// NextRefID = 5 with only 1 reference — strict mode rejects.
	if err := b.writeRecord(&Record{
		QName: "bad", RefID: 0, Pos: 100, MapQ: 60, Bin: 4681,
		Cigar: []uint32{cigar(CigarM, 5)},
		Seq:   "ACGTA", Qual: []byte("IIIII"),
		NextRefID: 5,
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	hdr, recs := readAll(t, data)
	issues := NewValidator(ModeStrict).Validate(hdr, recs)
	if len(issues) == 0 {
		t.Fatal("expected strict-mode next_refID error")
	}
	var found bool
	for _, iss := range issues {
		if iss.Field == "next_refID" {
			found = true
		}
	}
	if !found {
		t.Errorf("no next_refID issue found: %v", issues)
	}
}

func TestCigarOpChar_EdgeCases(t *testing.T) {
	tests := []struct {
		p    uint32
		want byte
	}{
		{cigar(CigarM, 0), 'M'},
		{cigar(CigarI, 100), 'I'},
		{cigar(CigarS, 50), 'S'},
		{cigar(CigarEq, 10), '='},
		{cigar(CigarX, 25), 'X'},
		// Unknown op: op code 15 is not one of the spec's 0-9.
		// Per the BAM layout op_len<<4 | op the op code is the LOW
		// 4 bits, so this is a zero-length op 15 rather than a huge M.
		{uint32(15), '?'},
	}
	for _, tc := range tests {
		if got := CigarOpChar(tc.p); got != tc.want {
			t.Errorf("CigarOpChar(%#x): got %q want %q", tc.p, got, tc.want)
		}
	}
}

func TestBaseChar_KnownCodes(t *testing.T) {
	for i := 0; i < 16; i++ {
		if baseChar(byte(i)) == 0 {
			t.Errorf("baseChar(%d) returned 0", i)
		}
	}
	if baseChar(42) != 'N' {
		t.Errorf("baseChar(42): got %c", baseChar(42))
	}
}

// ---------------------------------------------------------------------------
// Severity string tests
// ---------------------------------------------------------------------------

func TestSeverity_String(t *testing.T) {
	if SeverityWarning.String() != "warning" {
		t.Errorf("SeverityWarning: %q", SeverityWarning.String())
	}
	if SeverityError.String() != "error" {
		t.Errorf("SeverityError: %q", SeverityError.String())
	}
	if Severity(99).String() != "unknown" {
		t.Errorf("Severity(99): %q", Severity(99).String())
	}
}

// ---------------------------------------------------------------------------
// Parser edge cases
// ---------------------------------------------------------------------------

func TestParser_TruncatedHeader(t *testing.T) {
	// Magic + l_text=10 but no text.
	data := []byte{'B', 'A', 'M', 0x01, 10, 0, 0, 0}
	_, err := NewDecoder(bytes.NewReader(data)).Header()
	if err == nil {
		t.Fatal("expected error on truncated header")
	}
}

// ---------------------------------------------------------------------------
// B-array decode (hand-crafted aux bytes)
// ---------------------------------------------------------------------------

func TestRecord_BArray(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000}}); err != nil {
		t.Fatal(err)
	}
	if err := b.writeRecordBArray(&Record{
		QName: "bar", RefID: 0, Pos: 100, MapQ: 60, Bin: 4681,
		Cigar: []uint32{cigar(CigarM, 5)}, Seq: "ACGTA", Qual: []byte("IIIII"),
		NextRefID: -1,
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	_, recs := readAll(t, data)
	if len(recs) != 1 || len(recs[0].Tags) != 1 {
		t.Fatalf("expected 1 record with 1 tag, got %d records / %d tags",
			len(recs), len(recs[0].Tags))
	}
	tag := recs[0].Tags[0]
	if tag.Key != "XB" {
		t.Errorf("tag key: %q", tag.Key)
	}
	if tag.Type != TagB {
		t.Errorf("tag type: %q", string(byte(tag.Type)))
	}
	if tag.Value != "1,2,-3" {
		t.Errorf("tag value: %q want 1,2,-3", tag.Value)
	}
}

// writeRecordBArray writes a record whose only aux tag is a
// B-array with sub-type c (int8) and values [1, 2, -3].
func (b *builder) writeRecordBArray(r *Record) error {
	var rem bytes.Buffer
	rem.WriteString(r.QName)
	rem.WriteByte(0)
	for _, op := range r.Cigar {
		var c [4]byte
		binary.LittleEndian.PutUint32(c[:], op)
		rem.Write(c[:])
	}
	lSeq := int32(len(r.Seq))
	if lSeq > 0 {
		packed := make([]byte, (int(lSeq)+1)/2)
		for i := 0; i < int(lSeq); i++ {
			code := baseCode(r.Seq[i])
			if i%2 == 0 {
				packed[i/2] = code << 4
			} else {
				packed[i/2] |= code & 0x0F
			}
		}
		rem.Write(packed)
	}
	if lSeq > 0 {
		// Raw Phred bytes, not the ASCII "+33" form.
		for _, q := range r.Qual {
			rem.WriteByte(q - 33)
		}
	}
	// Aux: XB:B:c with [1, 2, -3].
	rem.WriteString("XB")
	rem.WriteByte(byte(TagB))
	rem.WriteByte(byte(Tagc))
	var cnt [4]byte
	binary.LittleEndian.PutUint32(cnt[:], 3)
	rem.Write(cnt[:])
	rem.WriteByte(1)
	rem.WriteByte(2)
	v := int8(-3)
	rem.WriteByte(byte(v))
	core := make([]byte, 32)
	binary.LittleEndian.PutUint32(core[0:4], uint32(r.RefID))
	binary.LittleEndian.PutUint32(core[4:8], uint32(int32(r.Pos)-1))
	core[8] = uint8(len(r.QName)) + 1
	core[9] = r.MapQ
	binary.LittleEndian.PutUint16(core[10:12], r.Bin)
	binary.LittleEndian.PutUint16(core[12:14], uint16(len(r.Cigar)))
	binary.LittleEndian.PutUint16(core[14:16], r.Flag)
	binary.LittleEndian.PutUint32(core[16:20], uint32(lSeq))
	binary.LittleEndian.PutUint32(core[20:24], uint32(r.NextRefID))
	binary.LittleEndian.PutUint32(core[24:28], 0xFFFFFFFF)
	binary.LittleEndian.PutUint32(core[28:32], 0)
	blockSize := uint32(32 + rem.Len())
	var bs [4]byte
	binary.LittleEndian.PutUint32(bs[:], blockSize)
	var payload bytes.Buffer
	payload.Write(bs[:])
	payload.Write(core)
	payload.Write(rem.Bytes())
	_, err := b.bw.Write(payload.Bytes())
	return err
}
