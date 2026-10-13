package vcf

import "testing"

func TestIsBreakend(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in   string
		want bool
	}{
		// Real breakend forms from VCFv4.5 §1.6.5 examples
		{"G]17:198982]", true},
		{"C[2:321682[", true},
		{"T[chrA:5[", true},
		{"]chrA:2]A", true},
		{"CAGTNNNNNCA[2:321682[", true},
		{"]13:123456]AGTNNNNNCAT", true},
		{"C[<ctg1>:1[", true},
		{"]<ctg1>:329]A", true},
		// Suffix-only style (replacement after the bracket pair)
		{"[chr1:5[A", true},
		{"]chr1:5]T", true},
		// Prefix-only style (replacement before the bracket pair)
		{"N]17:198982]", true},
		{"N[chr1:5[", true},
		// Empty bracket pair
		{"[]", false},
		{"[", false},
		{"[[", false},
		// Not a breakend — no position colon
		{"t[p[", false},
		{"t]p]", false},
		// Not a breakend — mismatched brackets
		{"t[p]", false},
		{"t]p[", false},
		// Plain nucleotide
		{"A", false},
		// Symbolic alleles
		{"<DEL>", false},
		{"<*>", false},
	}
	for _, c := range cases {
		if got := isBreakend(c.in); got != c.want {
			t.Errorf("isBreakend(%q) = %v; want %v", c.in, got, c.want)
		}
	}
}

func TestValidAltAllele(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in   string
		want bool
	}{
		{"A", true},
		{"ACGT", true},
		{"ACGTN", true},
		{"acgtn", true},
		{"*", true},
		{"<*>", true},
		{"<DEL>", true},
		{"<INS:ME>", true},
		{"<CNV:TR>", true},
		{"<bad space>", false},
		{"<>", false},
		{"<DEL<>>", false},
		{"", false},
		{"Z", false},
		{"T[chrA:5[", true},
		{"]chrA:2]A", true},
		{"[chr1:1[A", true},
	}
	for _, c := range cases {
		if got := validAltAllele(c.in); got != c.want {
			t.Errorf("validAltAllele(%q) = %v; want %v", c.in, got, c.want)
		}
	}
}

func TestParseVCFFloat(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in      string
		want    float64
		wantNaN bool
	}{
		// The spec calls for IEEE-754 32-bit floats, so we test with
		// values that round-trip through float32 exactly.
		{"1.5", 1.5, false},
		{"0", 0, false},
		{"0.5", 0.5, false},
		{"-1", -1, false},
		{"INF", positiveInf, false},
		{"inf", positiveInf, false},
		{"Infinity", positiveInf, false},
		{"-INF", negativeInf, false},
		{"-INFINITY", negativeInf, false},
		{"NaN", 0, true},
		{"nan", 0, true},
		{"foo", 0, false}, // error; we just check structural
	}
	for _, c := range cases {
		got, err := parseVCFFloat(c.in)
		if c.in == "foo" {
			if err == nil {
				t.Errorf("parseVCFFloat(%q) should error", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("parseVCFFloat(%q): %v", c.in, err)
			continue
		}
		if c.wantNaN {
			if !isNaN(got) {
				t.Errorf("parseVCFFloat(%q) = %v; want NaN", c.in, got)
			}
			continue
		}
		if got != c.want {
			t.Errorf("parseVCFFloat(%q) = %v; want %v", c.in, got, c.want)
		}
	}
}

func isNaN(f float64) bool {
	// math.IsNaN is the proper check.
	return f != f
}
