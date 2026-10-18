package vcf

import (
	"strings"
	"testing"
)

func TestIsReservedInfoKey(t *testing.T) {
	for _, key := range []string{"RN", "RUS", "RUL", "RUC", "RB", "CIRUC", "CIRB", "RUB", "LAA", "LAD"} {
		if !IsReservedInfoKey(key) {
			t.Errorf("IsReservedInfoKey(%q): should be true", key)
		}
	}
	for _, key := range []string{"DP", "AF", "AC", "NS", "MN", "RM", "FOO"} {
		if IsReservedInfoKey(key) {
			t.Errorf("IsReservedInfoKey(%q): should be false", key)
		}
	}
}

func TestIsReservedFormatKey(t *testing.T) {
	for _, key := range []string{"LAA", "LAD", "LADF", "LADR", "LEC", "LGL", "LGP", "LPL", "LPP"} {
		if !IsReservedFormatKey(key) {
			t.Errorf("IsReservedFormatKey(%q): should be true", key)
		}
	}
	for _, key := range []string{"GT", "DP", "PL", "GL", "FOO"} {
		if IsReservedFormatKey(key) {
			t.Errorf("IsReservedFormatKey(%q): should be false", key)
		}
	}
}

func TestIsBaseModKey(t *testing.T) {
	cases := []struct {
		in   string
		want bool
	}{
		// Positive cases.
		{"M0A", true},
		{"M5mC", true},
		{"M5hmC", true},
		{"M5fC", true},
		{"M6mA", true},
		{"M1C", true},
		{"M100T", true},
		{"M0", true}, // digit-only after M is also valid
		// Negative cases.
		{"M", false},  // no digit
		{"MC", false}, // no digit before code
		{"N5mC", false}, // wrong prefix
		{"", false},
		{"m5mC", false}, // case-sensitive: M is uppercase per spec
	}
	for _, c := range cases {
		if got := IsBaseModKey(c.in); got != c.want {
			t.Errorf("IsBaseModKey(%q): got %v, want %v", c.in, got, c.want)
		}
	}
}

func TestValidator_ReservedKeyMetadataMismatch(t *testing.T) {
	// RN is reserved Number=1 Type=Integer. A header that
	// declares RN as Type=String should trigger the warning.
	src := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=RN,Number=1,Type=String,Description=\"tandem repeat ref number\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(hdr, nil)
	found := false
	for _, iss := range issues {
		if iss.Field == "header/INFO/RN" {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected reserved-key warning for RN, got issues: %+v", issues)
	}
}

func TestValidator_ReservedKeyMetadataOK(t *testing.T) {
	// RN with the spec-correct Number=1, Type=Integer must
	// not produce a warning.
	src := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=RN,Number=1,Type=Integer,Description=\"tandem repeat ref number\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(hdr, nil)
	for _, iss := range issues {
		if iss.Field == "header/INFO/RN" {
			t.Errorf("unexpected reserved-key warning for correctly-declared RN: %+v", iss)
		}
	}
}