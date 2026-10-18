package vcf

import (
	"strings"
	"testing"
)

// iupacHeader is a minimal two-column header used by the
// IUPAC reduction tests.
const iupacHeader = "##fileformat=VCFv4.2\n" +
	"##contig=<ID=chr1,length=1000>\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

func TestReduceIUPACBases(t *testing.T) {
	// The helper is a pure string function.
	cases := []struct {
		in, want string
	}{
		// Plain ACGTN passes through.
		{"ACGTN", "ACGTN"},
		{"acgtn", "acgtn"},
		// Lower case IUPAC codes are reduced too.
		{"R", "A"},
		{"r", "a"},
		{"Y", "C"},
		{"y", "c"},
		{"S", "C"},
		{"W", "A"},
		{"K", "G"},
		{"k", "g"},
		{"M", "A"},
		{"m", "a"},
		{"B", "C"},
		{"D", "A"},
		{"H", "A"},
		{"V", "A"},
		{"v", "a"},
		// Mixtures: each IUPAC base is reduced in place.
		{"RYSWKMBDHV", "ACCAGACAAA"}, // RY->AC; SW->CA; KM->AG; BDHV->CAAA
		{"AR", "AA"},
		{"RC", "AC"},
		// Strings without IUPAC codes are returned unchanged
		// (and the helper does not allocate for them).
		{"", ""},
		{"ACGT", "ACGT"},
	}
	for _, c := range cases {
		got := reduceIUPACBases(c.in)
		if got != c.want {
			t.Errorf("reduceIUPACBases(%q) = %q, want %q", c.in, got, c.want)
		}
	}
}

func TestDecoder_REFIUPACReduction(t *testing.T) {
	// Each IUPAC ambiguity code in REF should be reduced
	// to one of its constituent bases by the parser.
	cases := []struct{ in, want string }{
		{"R", "A"},
		{"Y", "C"},
		{"S", "C"},
		{"W", "A"},
		{"K", "G"},
		{"M", "A"},
		{"B", "C"},
		{"D", "A"},
		{"H", "A"},
		{"V", "A"},
		{"r", "a"},
		{"ACRYGT", "ACACGT"},
	}
	for _, c := range cases {
		rec := parseSingleRef(t, c.in)
		if rec.Ref != c.want {
			t.Errorf("REF=%q parsed as %q, want %q", c.in, rec.Ref, c.want)
		}
		// After reduction, REF should be valid (only
		// A/C/G/T/N), so strict validation must accept
		// it.
		dec := NewDecoder(strings.NewReader(iupacHeader + "chr1\t100\t.\t" + c.in + "\tA\t.\t.\t.\n"))
		hdr, err := dec.Header()
		if err != nil {
			t.Fatalf("Header: %v", err)
		}
		rec2, err := dec.Record()
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		for _, iss := range NewValidator(ModeStrict).Validate(hdr, []*Record{rec2}) {
			if strings.Contains(iss.Msg, "invalid REF allele") {
				t.Errorf("REF=%q: strict validator still flags REF after reduction: %+v", c.in, iss)
			}
		}
		_ = rec
	}
}

func parseSingleRef(t *testing.T, ref string) *Record {
	t.Helper()
	dec := NewDecoder(strings.NewReader(iupacHeader + "chr1\t100\t.\t" + ref + "\tA\t.\t.\t.\n"))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	return rec
}