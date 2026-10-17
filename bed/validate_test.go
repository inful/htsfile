package bed

import (
	"strings"
	"testing"
)

func TestValidator_CleanBED6(t *testing.T) {
	t.Parallel()
	plus := StrandPlus
	recs := []*Record{
		{Line: 2, Chrom: "chr1", ChromStart: 100, ChromEnd: 200,
			Name: strPtr("a"), Score: intPtr(500), Strand: &plus},
		{Line: 3, Chrom: "chr1", ChromStart: 300, ChromEnd: 400,
			Name: strPtr("b"), Score: intPtr(900), Strand: &plus},
		{Line: 4, Chrom: "chr2", ChromStart: 500, ChromEnd: 600,
			Name: strPtr("c"), Score: intPtr(750), Strand: &plus},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 0 {
		t.Errorf("issues = %+v; want none", issues)
	}
}

func TestValidator_InconsistentColumnCount(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, Chrom: "chr1", ChromStart: 100, ChromEnd: 200}, // BED3
		{Line: 2, Chrom: "chr1", ChromStart: 300, ChromEnd: 400,
			Name: strPtr("x")}, // BED4
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("len(issues) = %d; want 1", len(issues))
	}
	if !strings.Contains(issues[0].Msg, "column count") {
		t.Errorf("Msg = %q; want it to mention column count", issues[0].Msg)
	}
}

func TestValidator_NonContiguousPrefix(t *testing.T) {
	t.Parallel()
	// Score is set but Name is nil — non-contiguous prefix.
	recs := []*Record{
		{Line: 1, Chrom: "chr1", ChromStart: 100, ChromEnd: 200,
			Score: intPtr(500)},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) == 0 {
		t.Fatal("expected an issue for non-contiguous prefix")
	}
	if !strings.Contains(issues[0].Msg, "non-contiguous") {
		t.Errorf("Msg = %q; want 'non-contiguous'", issues[0].Msg)
	}
}

func TestValidator_NegativeStart(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, Chrom: "chr1", ChromStart: -1, ChromEnd: 200},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("len(issues) = %d; want 1", len(issues))
	}
	if !strings.Contains(issues[0].Msg, ">= 0") {
		t.Errorf("Msg = %q", issues[0].Msg)
	}
}

func TestValidator_StartAfterEnd(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, Chrom: "chr1", ChromStart: 200, ChromEnd: 100},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("len(issues) = %d; want 1", len(issues))
	}
	if !strings.Contains(issues[0].Msg, ">") {
		t.Errorf("Msg = %q", issues[0].Msg)
	}
}

func TestValidator_ScoreOutOfRange(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, Chrom: "chr1", ChromStart: 100, ChromEnd: 200,
			Name: strPtr("a"), Score: intPtr(1500), Strand: strPtrHelper(StrandPlus)},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 1 {
		t.Fatalf("len(issues) = %d; want 1", len(issues))
	}
	if !strings.Contains(issues[0].Msg, "[0, 1000]") {
		t.Errorf("Msg = %q", issues[0].Msg)
	}
}

func TestValidator_BED12BlockInvariants(t *testing.T) {
	t.Parallel()
	cases := []struct {
		name        string
		modify      func(r *Record)
		wantContain string
	}{
		{
			name: "blockCount mismatch",
			modify: func(r *Record) {
				n := 2
				r.BlockCount = &n
				// BlockSizes and BlockStarts have 3 entries
			},
			wantContain: "blockSizes has 3 entries",
		},
		{
			name: "first blockStart not 0",
			modify: func(r *Record) {
				r.BlockStarts = []int64{5, 100}
			},
			wantContain: "first blockStart must be 0",
		},
		{
			name: "blockSize 0",
			modify: func(r *Record) {
				r.BlockSizes = []int64{100, 0, 300}
			},
			wantContain: "blockSize[1] = 0",
		},
		{
			name: "overlapping blocks",
			modify: func(r *Record) {
				r.BlockStarts = []int64{0, 50}
			},
			wantContain: "overlaps",
		},
		{
			name: "sum invariant violated",
			modify: func(r *Record) {
				// region length is 4000; last block (start=3500, size=200)
				// ends at 3700, not 4000.
				r.BlockStarts = []int64{0, 3500}
				r.BlockSizes = []int64{100, 200}
			},
			wantContain: "ends at 3700",
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			t.Parallel()
			r := &Record{
				Line: 1, Chrom: "chr1", ChromStart: 1000, ChromEnd: 5000,
				Name: strPtr("g"), Score: intPtr(900), Strand: strPtrHelper(StrandPlus),
				ThickStart: int64Ptr(1200), ThickEnd: int64Ptr(4800),
				ItemRgb:     &RGB{0, 0, 255},
				BlockCount:  intPtr(3),
				BlockSizes:  []int64{100, 200, 300},
				BlockStarts: []int64{0, 1500, 3700},
			}
			c.modify(r)
			v := NewValidator(ModeDefault)
			issues := v.Validate(nil, []*Record{r})
			if len(issues) == 0 {
				t.Fatal("expected an issue")
			}
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, c.wantContain) {
					found = true
					break
				}
			}
			if !found {
				t.Errorf("issues did not contain %q: %+v", c.wantContain, issues)
			}
		})
	}
}

func TestValidator_ZeroLengthRegionIsValid(t *testing.T) {
	t.Parallel()
	// Per spec, chromStart == chromEnd is a valid insertion point.
	recs := []*Record{
		{Line: 1, Chrom: "chr1", ChromStart: 100, ChromEnd: 100},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(nil, recs)
	if len(issues) != 0 {
		t.Errorf("zero-length region should be valid; got %+v", issues)
	}
}

func TestValidator_StrictModeDetectsUnsorted(t *testing.T) {
	t.Parallel()
	recs := []*Record{
		{Line: 1, Chrom: "chr1", ChromStart: 200, ChromEnd: 300},
		{Line: 2, Chrom: "chr1", ChromStart: 100, ChromEnd: 200}, // out of order
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(nil, recs)
	if len(issues) == 0 {
		t.Fatal("expected sort-order warning in strict mode")
	}
	if issues[0].Severity != SeverityWarning {
		t.Errorf("Severity = %v; want SeverityWarning", issues[0].Severity)
	}
}

func TestSeverity_String(t *testing.T) {
	t.Parallel()
	if SeverityWarning.String() != "warning" {
		t.Errorf("Warning.String() = %q; want warning", SeverityWarning.String())
	}
	if SeverityError.String() != "error" {
		t.Errorf("Error.String() = %q; want error", SeverityError.String())
	}
	if Severity(99).String() != "unknown" {
		t.Errorf("Severity(99).String() = %q; want unknown", Severity(99).String())
	}
}

func strPtrHelper(s Strand) *Strand { return &s }
