package sam

import (
	"strings"
	"testing"
)

func TestParseCigarValid(t *testing.T) {
	cases := []struct {
		in   string
		want []CigarUnit
	}{
		// Empty / "*" sentinels.
		{"*", nil},

		// Single op.
		{"5M", []CigarUnit{{5, CigarM}}},
		{"1I", []CigarUnit{{1, CigarI}}},
		{"10=", []CigarUnit{{10, CigarEq}}},

		// Multi-op.
		{"5M2I3M", []CigarUnit{{5, CigarM}, {2, CigarI}, {3, CigarM}}},
		{"10M200N10M", []CigarUnit{{10, CigarM}, {200, CigarN}, {10, CigarM}}},
		{"5S20M5S", []CigarUnit{{5, CigarS}, {20, CigarM}, {5, CigarS}}},
		{"5H10M5H", []CigarUnit{{5, CigarH}, {10, CigarM}, {5, CigarH}}},
		{"2X3=2X", []CigarUnit{{2, CigarX}, {3, CigarEq}, {2, CigarX}}},

		// All nine ops.
		{"1M1I1D1N1S1H1P1=1X", []CigarUnit{
			{1, CigarM}, {1, CigarI}, {1, CigarD}, {1, CigarN},
			{1, CigarS}, {1, CigarH}, {1, CigarP},
			{1, CigarEq}, {1, CigarX},
		}},

		// Large lengths.
		{"1000000M", []CigarUnit{{1000000, CigarM}}},
	}
	for _, c := range cases {
		got, err := ParseCigar(c.in)
		if err != nil {
			t.Errorf("ParseCigar(%q): unexpected error %v", c.in, err)
			continue
		}
		if !cigarUnitsEqual(got, c.want) {
			t.Errorf("ParseCigar(%q): got %v, want %v", c.in, got, c.want)
		}
	}
}

func TestParseCigarLowercase(t *testing.T) {
	// Lower-case op chars are accepted (per SAM v1 the file is
	// case-sensitive, but real-world aligners emit lowercase 'm'
	// in extended-CIGAR contexts).
	got, err := ParseCigar("5m2i3m")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	want := []CigarUnit{{5, CigarM}, {2, CigarI}, {3, CigarM}}
	if !cigarUnitsEqual(got, want) {
		t.Errorf("ParseCigar(\"5m2i3m\"): got %v, want %v", got, want)
	}
}

func TestParseCigarErrors(t *testing.T) {
	cases := []string{
		"",    // empty
		"M",   // missing length
		"5",   // missing op
		"5Q",  // unknown op
		"5MM", // missing length before second M
		"abc", // non-numeric length
		"-1M", // negative length (the '-' trips the parser)
		"5M ", // trailing whitespace
		"5M2", // trailing length without op
	}
	for _, in := range cases {
		_, err := ParseCigar(in)
		if err == nil {
			t.Errorf("ParseCigar(%q): expected error, got nil", in)
			continue
		}
		if !strings.Contains(err.Error(), "cigar") {
			t.Errorf("ParseCigar(%q): error %v should mention \"cigar\"", in, err)
		}
	}
}

func TestFormatCigarRoundTrip(t *testing.T) {
	cases := []string{
		"5M",
		"3M2I3M",
		"10M200N10M",
		"5S20M5S",
		"5H10M5H",
		"2X3=2X",
		"1M1I1D1N1S1H1P1=1X",
	}
	for _, in := range cases {
		parsed, err := ParseCigar(in)
		if err != nil {
			t.Fatalf("ParseCigar(%q): %v", in, err)
		}
		got := FormatCigar(parsed)
		if got != in {
			t.Errorf("round-trip %q: got %q", in, got)
		}
	}
}

func TestFormatCigarEmpty(t *testing.T) {
	if got := FormatCigar(nil); got != "*" {
		t.Errorf("FormatCigar(nil): got %q, want %q", got, "*")
	}
	if got := FormatCigar([]CigarUnit{}); got != "*" {
		t.Errorf("FormatCigar([]): got %q, want %q", got, "*")
	}
}

func TestFormatCigarDropsZeroLength(t *testing.T) {
	ops := []CigarUnit{{5, CigarM}, {0, CigarI}, {3, CigarM}}
	got := FormatCigar(ops)
	if got != "5M3M" {
		t.Errorf("zero-length drop: got %q, want %q", got, "5M3M")
	}
	if got := FormatCigar([]CigarUnit{{0, CigarM}}); got != "*" {
		t.Errorf("all-zero: got %q, want %q", got, "*")
	}
}

func TestCigarOpStringAndChar(t *testing.T) {
	for op := CigarM; op <= CigarX; op++ {
		if !op.IsValid() {
			t.Errorf("op %d should be valid", op)
		}
		if op.Char() == 0 {
			t.Errorf("op %d: Char() returned 0", op)
		}
		if op.String() == "?" {
			t.Errorf("op %d: String() returned \"?\"", op)
		}
	}
	if CigarOp(99).IsValid() {
		t.Error("CigarOp(99) should be invalid")
	}
	if CigarOp(99).Char() != 0 {
		t.Error("CigarOp(99).Char() should be 0")
	}
	if CigarOp(99).String() != "?" {
		t.Errorf("CigarOp(99).String(): got %q, want \"?\"", CigarOp(99).String())
	}
}

func TestCigarOpConsumes(t *testing.T) {
	// Reference-consuming ops: M, D, N, =, X.
	for _, op := range []CigarOp{CigarM, CigarD, CigarN, CigarEq, CigarX} {
		if !op.ConsumesReference() {
			t.Errorf("op %s should consume reference", op)
		}
	}
	for _, op := range []CigarOp{CigarI, CigarS, CigarH, CigarP} {
		if op.ConsumesReference() {
			t.Errorf("op %s should NOT consume reference", op)
		}
	}
	// Query-consuming ops: M, I, S, =, X.
	for _, op := range []CigarOp{CigarM, CigarI, CigarS, CigarEq, CigarX} {
		if !op.ConsumesQuery() {
			t.Errorf("op %s should consume query", op)
		}
	}
	for _, op := range []CigarOp{CigarD, CigarN, CigarH, CigarP} {
		if op.ConsumesQuery() {
			t.Errorf("op %s should NOT consume query", op)
		}
	}
}

func TestReferenceAndQueryLength(t *testing.T) {
	ops, err := ParseCigar("5M2I3M2D4N5S5H3=2X")
	if err != nil {
		t.Fatalf("ParseCigar: %v", err)
	}
	// Reference: 5 (M) + 3 (M) + 2 (D) + 4 (N) + 3 (=) + 2 (X) = 19
	if got := ReferenceLength(ops); got != 19 {
		t.Errorf("ReferenceLength: got %d, want 19", got)
	}
	// Query: 5 (M) + 2 (I) + 3 (M) + 5 (S) + 3 (=) + 2 (X) = 20
	if got := QueryLength(ops); got != 20 {
		t.Errorf("QueryLength: got %d, want 20", got)
	}
}

func TestReferenceAndQueryLengthEmpty(t *testing.T) {
	if got := ReferenceLength(nil); got != 0 {
		t.Errorf("ReferenceLength(nil): got %d, want 0", got)
	}
	if got := QueryLength(nil); got != 0 {
		t.Errorf("QueryLength(nil): got %d, want 0", got)
	}
	if got := ReferenceLength([]CigarUnit{{5, CigarH}}); got != 0 {
		t.Errorf("ReferenceLength(hard-clip-only): got %d, want 0", got)
	}
}

func TestAlignmentEnd(t *testing.T) {
	ops, err := ParseCigar("5M2I3M")
	if err != nil {
		t.Fatalf("ParseCigar: %v", err)
	}
	// ReferenceLength(ops) = 8, so AlignmentEnd(100, ops) = 108.
	if got := AlignmentEnd(100, ops); got != 108 {
		t.Errorf("AlignmentEnd(100, 5M2I3M): got %d, want 108", got)
	}
	// Unmapped (all non-reference-consuming ops): alignment end
	// collapses to pos (typically 0 per SAM convention).
	if got := AlignmentEnd(0, []CigarUnit{{5, CigarS}, {5, CigarS}}); got != 0 {
		t.Errorf("AlignmentEnd(unmapped): got %d, want 0", got)
	}
}

func TestParseCigarBoundaryLength(t *testing.T) {
	// Length at int32 max must succeed.
	_, err := ParseCigar("2147483647M")
	if err != nil {
		t.Errorf("int32-max length: unexpected error %v", err)
	}
	// Length above int32 max must fail.
	_, err = ParseCigar("2147483648M")
	if err == nil {
		t.Error("int32-max+1 length: expected error")
	}
}

// cigarUnitsEqual compares two CigarUnit slices element-wise.
// Returns true when both are the same length and each pair
// matches in Len and Op.
func cigarUnitsEqual(a, b []CigarUnit) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i].Len != b[i].Len || a[i].Op != b[i].Op {
			return false
		}
	}
	return true
}
