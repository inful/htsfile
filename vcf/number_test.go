package vcf

import (
	"errors"
	"testing"
)

func TestParseNumber(t *testing.T) {
	t.Parallel()

	cases := []struct {
		in      string
		want    Number
		wantErr bool
	}{
		{"0", Number{Number: 0}, false},
		{"1", Number{Number: 1}, false},
		{"2", Number{Number: 2}, false},
		{"100", Number{Number: 100}, false},
		{"A", Number{Fixed: SymbolA}, false},
		{"R", Number{Fixed: SymbolR}, false},
		{"G", Number{Fixed: SymbolG}, false},
		{".", Number{Fixed: SymbolDot}, false},
		{"LR", Number{Fixed: SymbolLR}, false},
		{"LA", Number{Fixed: SymbolLA}, false},
		{"LG", Number{Fixed: SymbolLG}, false},
		// P and M map to SymbolDot (variable cardinality).
		{"P", Number{Fixed: SymbolDot}, false},
		{"M", Number{Fixed: SymbolDot}, false},
		// Errors
		{"-1", Number{}, true},
		{"foo", Number{}, true},
		{"", Number{}, true},
		{"1.5", Number{}, true},
	}
	for _, c := range cases {
		t.Run(c.in, func(t *testing.T) {
			t.Parallel()
			got, err := ParseNumber(c.in)
			if c.wantErr {
				if err == nil {
					t.Fatalf("expected error for %q, got %+v", c.in, got)
				}
				if !errors.Is(err, err) { // placeholder; just ensure err is non-nil
					t.Fatal("error is nil")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if got != c.want {
				t.Fatalf("got %+v, want %+v", got, c.want)
			}
		})
	}
}

func TestNumberString(t *testing.T) {
	t.Parallel()

	cases := []Number{
		{Number: 1},
		{Number: 0},
		{Fixed: SymbolA},
		{Fixed: SymbolR},
		{Fixed: SymbolG},
		{Fixed: SymbolDot},
		{Fixed: SymbolLA},
		{Fixed: SymbolLR},
		{Fixed: SymbolLG},
	}
	wants := []string{"1", "0", "A", "R", "G", ".", "LA", "LR", "LG"}
	for i, n := range cases {
		if got := n.String(); got != wants[i] {
			t.Errorf("Number%+v.String() = %q; want %q", n, got, wants[i])
		}
	}
}

func TestNumberCardinality(t *testing.T) {
	t.Parallel()

	cases := []struct {
		name   string
		n      Number
		alts   int
		ploidy int
		want   int
	}{
		{"literal-1 one-alt", Number{Number: 1}, 2, 0, 1},
		{"literal-0 zero", Number{Number: 0}, 5, 2, 0},
		{"literal-3 three", Number{Number: 3}, 0, 0, 3},
		{"A two-alts", Number{Fixed: SymbolA}, 2, 0, 2},
		{"A zero-alts", Number{Fixed: SymbolA}, 0, 0, 0},
		{"R two-alts", Number{Fixed: SymbolR}, 2, 0, 3},          // REF + 2 alts
		{"R zero-alts", Number{Fixed: SymbolR}, 0, 0, 1},         // REF only
		{"G diploid biallelic", Number{Fixed: SymbolG}, 1, 2, 3}, // 00,01,11
		{"G diploid 2 alts", Number{Fixed: SymbolG}, 2, 2, 6},    // 00,01,11,02,12,22
		{"G triploid 2 alts", Number{Fixed: SymbolG}, 2, 3, 10},  // 000,001,011,111,002,012,112,022,122,222
		{"G haploid biallelic", Number{Fixed: SymbolG}, 1, 1, 2}, // 0,1
		{"G missing ploidy", Number{Fixed: SymbolG}, 1, 0, -1},
		{"Dot variable", Number{Fixed: SymbolDot}, 1, 1, -1},
		{"LA unknown LAA", Number{Fixed: SymbolLA}, 3, 0, -1},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			t.Parallel()
			if got := c.n.Cardinality(c.alts, c.ploidy); got != c.want {
				t.Errorf("%+v.Cardinality(alts=%d, ploidy=%d) = %d; want %d",
					c.n, c.alts, c.ploidy, got, c.want)
			}
		})
	}
}

func TestGenotypeCount(t *testing.T) {
	t.Parallel()
	cases := []struct {
		alleles, ploidy, want int
	}{
		{2, 1, 2},  // 0,1
		{2, 2, 3},  // 00,01,11
		{3, 2, 6},  // 00,01,11,02,12,22
		{3, 3, 10}, // 000,...222
		{1, 5, 1},  // 00000
		{4, 2, 10}, // 00,01,11,02,12,22,03,13,23,33
		{5, 2, 15},
		{0, 2, -1}, // error
		{2, 0, -1}, // error
	}
	for _, c := range cases {
		got := genotypeCount(c.alleles, c.ploidy)
		if got != c.want {
			t.Errorf("genotypeCount(alleles=%d, ploidy=%d) = %d; want %d",
				c.alleles, c.ploidy, got, c.want)
		}
	}
}
