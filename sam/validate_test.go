package sam

import (
	"strings"
	"testing"
)

func TestValidator_Clean(t *testing.T) {
	t.Parallel()
	hdr := &Header{
		SQs: []*SQ{{SN: "chr1", LN: 1000}},
	}
	recs := []*Record{
		{Line: 2, QName: "r1", Flag: 0, RName: "chr1", Pos: 100,
			MapQ: 60, Cigar: "10M", RNext: "*", PNext: 0, TLen: 0,
			Seq: "ACGTACGTAC", Qual: "IIIIIIIIII"},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 0 {
		t.Errorf("issues = %+v; want none", issues)
	}
}

func TestValidator_MapQOutOfRange(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, MapQ: 300, Seq: "ACGT", Qual: "IIII"},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("issues = %+v; want 1", issues)
	}
	if !strings.Contains(issues[0].Msg, "MAPQ must be in") {
		t.Errorf("Msg = %q", issues[0].Msg)
	}
}

func TestValidator_RNameNotDeclared(t *testing.T) {
	t.Parallel()
	hdr := &Header{
		SQs: []*SQ{{SN: "chr1", LN: 1000}},
	}
	recs := []*Record{
		{Line: 1, RName: "chr2", Seq: "ACGT", Qual: "IIII"},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("issues = %+v; want 1", issues)
	}
	if !strings.Contains(issues[0].Msg, "not declared") {
		t.Errorf("Msg = %q", issues[0].Msg)
	}
}

func TestValidator_RNameStarIsOK(t *testing.T) {
	t.Parallel()
	// RNAME "*" is allowed even if @SQ doesn't contain it.
	hdr := &Header{
		SQs: []*SQ{{SN: "chr1", LN: 1000}},
	}
	recs := []*Record{
		{Line: 1, RName: "*", Seq: "ACGT", Qual: "IIII"},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 0 {
		t.Errorf("issues = %+v; want none (RNAME * is allowed)", issues)
	}
}

func TestValidator_SeqQualLengthMismatch(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, Seq: "ACGT", Qual: "III"}, // Qual too short
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("issues = %+v; want 1", issues)
	}
	if !strings.Contains(issues[0].Msg, "QUAL length") {
		t.Errorf("Msg = %q", issues[0].Msg)
	}
}

func TestValidator_SeqStarQualStar(t *testing.T) {
	t.Parallel()
	// Both "*" is the documented "not stored" case.
	recs := []*Record{
		{Line: 1, Seq: "*", Qual: "*"},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 0 {
		t.Errorf("issues = %+v; want none (both SEQ and QUAL are *)", issues)
	}
}

func TestValidator_SeqStarQualNotStar(t *testing.T) {
	t.Parallel()
	// One is "*" but not both: invalid.
	recs := []*Record{
		{Line: 1, Seq: "*", Qual: "IIII"},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("issues = %+v; want 1", issues)
	}
}

func TestValidator_StrictModePosForMapped(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, Flag: 0, Pos: 0, Seq: "ACGT", Qual: "IIII"}, // mapped but Pos=0
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("issues = %+v; want 1", issues)
	}
	if !strings.Contains(issues[0].Msg, "POS") {
		t.Errorf("Msg = %q", issues[0].Msg)
	}
}

func TestValidator_StrictModeUnmappedPosZeroOK(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, Flag: FlagUnmap, Pos: 0, Seq: "ACGT", Qual: "IIII"},
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(nil, recs)
	if len(issues) != 0 {
		t.Errorf("issues = %+v; want none (unmapped with Pos=0 is fine)", issues)
	}
}

func TestValidator_DuplicateTagKey(t *testing.T) {
	t.Parallel()
	// Flag=FlagUnmap so the strict-mode POS check doesn't also fire.
	recs := []*Record{
		{Line: 1, Flag: FlagUnmap, Seq: "ACGT", Qual: "IIII", Tags: []Tag{
			{Key: "NM", Type: Tagi, Value: "1"},
			{Key: "NM", Type: Tagi, Value: "2"}, // duplicate
		}},
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("issues = %+v; want 1", issues)
	}
	if issues[0].Severity != SeverityWarning {
		t.Errorf("Severity = %v; want SeverityWarning", issues[0].Severity)
	}
}

func TestSeverity_String(t *testing.T) {
	t.Parallel()
	if SeverityWarning.String() != "warning" {
		t.Errorf("Warning = %q; want warning", SeverityWarning.String())
	}
	if SeverityError.String() != "error" {
		t.Errorf("Error = %q; want error", SeverityError.String())
	}
	if Severity(99).String() != "unknown" {
		t.Errorf("Severity(99) = %q; want unknown", Severity(99).String())
	}
}
