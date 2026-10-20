package vcf

import (
	"strings"
	"testing"
)

func TestIsReservedInfoKey(t *testing.T) {
	// The v4.5-new structural-variant subset.
	for _, key := range []string{"RN", "RUS", "RUL", "RUC", "RB", "CIRUC", "CIRB", "RUB", "LAA", "LAD"} {
		if !IsReservedInfoKey(key) {
			t.Errorf("IsReservedInfoKey(%q): should be true", key)
		}
	}
	// The Reserved INFO keys table (spec lines 460-480). T4.25
	// added these; before it DP/AF/AC/NS were absent and this
	// test asserted they were NOT reserved, which pinned the gap.
	for _, key := range []string{
		"AA", "AC", "AD", "ADF", "ADR", "AF", "AN", "BQ", "CIGAR",
		"DB", "DP", "END", "H2", "H3", "MQ", "MQ0", "NS", "SB",
		"SOMATIC", "VALIDATED", "1000G",
	} {
		if !IsReservedInfoKey(key) {
			t.Errorf("IsReservedInfoKey(%q): should be true (spec Reserved INFO keys table)", key)
		}
	}
	for _, key := range []string{"MN", "RM", "FOO", "ZZZ"} {
		if IsReservedInfoKey(key) {
			t.Errorf("IsReservedInfoKey(%q): should be false", key)
		}
	}
}

func TestIsReservedFormatKey(t *testing.T) {
	for _, key := range []string{
		"LAA", "LAD", "LADF", "LADR", "LEC", "LGL", "LGP", "LPL", "LPP",
		"PSL", "PSO", "PSQ",
		"PL", "GL", "GP",
	} {
		if !IsReservedFormatKey(key) {
			t.Errorf("IsReservedFormatKey(%q): should be true", key)
		}
	}
	// The Reserved genotype keys table (spec lines 534-596),
	// excluding the M/DPM/ADM base-modification families which
	// resolve through IsBaseModKey / MFamily instead.
	for _, key := range []string{
		"AD", "ADF", "ADR", "DP", "EC", "LEN", "FT", "GQ", "GT",
		"HQ", "LA", "MQ", "PP", "PQ", "PS",
	} {
		if !IsReservedFormatKey(key) {
			t.Errorf("IsReservedFormatKey(%q): should be true (spec Reserved genotype keys table)", key)
		}
	}
	for _, key := range []string{"FOO", "ZZZ"} {
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
		{"M", false},    // no digit
		{"MC", false},   // no digit before code
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
	// RN with the spec-correct Number=A, Type=Integer must
	// not produce a warning (VCFv4.5 §1.6.1.7).
	src := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=RN,Number=A,Type=Integer,Description=\"tandem repeat ref number\">\n" +
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

func TestMFamily(t *testing.T) {
	cases := []struct {
		in   string
		want string
	}{
		{"M0A", "M"},
		{"M5mC", "M"},
		{"M5hmC", "M"},
		{"M5fC", "M"},
		{"M5caC", "M"},
		{"M5hmU", "M"},
		{"M6mA", "M"},
		{"M8oxoG", "M"},
		{"M27551C", "M"},
		{"M1", "M"},
		// DPM and ADM families.
		{"DPM5mC", "DPM"},
		{"DPM0A", "DPM"},
		{"ADM5mC", "ADM"},
		{"ADM0A", "ADM"},
		// Non-matches.
		{"MN", ""},
		{"MR", ""},
		{"M", ""},
		{"MC", ""},    // no digit
		{"M5mC ", ""}, // trailing whitespace
		{"", ""},
		{"m5mC", ""},    // case-sensitive: M is uppercase
		{"X5mC", ""},    // wrong prefix
		{"DPM", ""},     // no digit
		{"ADPM5mC", ""}, // wrong prefix
	}
	for _, c := range cases {
		if got := MFamily(c.in); got != c.want {
			t.Errorf("MFamily(%q): got %q, want %q", c.in, got, c.want)
		}
	}
}

func TestValidator_MValueRange(t *testing.T) {
	// M field: float in [0, 1]. Out-of-range values trigger an
	// error in default mode.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=M5mC,Number=.,Type=Float,Description=\"5mC fraction\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	dec := NewDecoder(strings.NewReader(hdrText))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	// Out-of-range: -0.1.
	rec := &Record{
		Chrom: "chr1", Pos: 100, Ref: "C", Alt: []string{"."},
		Format: []string{"M5mC"},
		Samples: []Genotype{{Fields: map[string]FieldValue{
			"M5mC": {typedCore: typedCore{Type: TypeFloat, Floats: []float64{-0.1}}},
		}}},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "FORMAT/M5mC/sample[0]" {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected M5mC out-of-range issue, got: %+v", issues)
	}
	// In-range: 0.95 — no issue.
	rec.Samples[0].Fields["M5mC"] = FieldValue{typedCore: typedCore{Type: TypeFloat, Floats: []float64{0.95}}}
	issues = v.Validate(hdr, []*Record{rec})
	for _, iss := range issues {
		if iss.Field == "FORMAT/M5mC/sample[0]" {
			t.Errorf("unexpected M5mC issue for in-range value: %+v", iss)
		}
	}
}

func TestValidator_MValueRangeUpper(t *testing.T) {
	// M field: 1.5 is out of range.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=M5mC,Number=.,Type=Float,Description=\"5mC fraction\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	dec := NewDecoder(strings.NewReader(hdrText))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec := &Record{
		Chrom: "chr1", Pos: 100, Ref: "C", Alt: []string{"."},
		Format: []string{"M5mC"},
		Samples: []Genotype{{Fields: map[string]FieldValue{
			"M5mC": {typedCore: typedCore{Type: TypeFloat, Floats: []float64{1.5}}},
		}}},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "FORMAT/M5mC/sample[0]" {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected M5mC upper-bound issue, got: %+v", issues)
	}
}

func TestValidator_DPMNegative(t *testing.T) {
	// DPM field: negative depth is an error.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=DPM5mC,Number=.,Type=Integer,Description=\"total depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	dec := NewDecoder(strings.NewReader(hdrText))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec := &Record{
		Chrom: "chr1", Pos: 100, Ref: "C", Alt: []string{"."},
		Format: []string{"DPM5mC"},
		Samples: []Genotype{{Fields: map[string]FieldValue{
			"DPM5mC": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{-1}}},
		}}},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "FORMAT/DPM5mC/sample[0]" && strings.Contains(iss.Msg, "negative") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected DPM5mC negative-depth issue, got: %+v", issues)
	}
}

func TestValidator_PSLMetadataMismatch(t *testing.T) {
	// PSL is reserved Number=P, Type=String. A header that
	// declares PSL as Type=Integer triggers a warning.
	src := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=PSL,Number=P,Type=Integer,Description=\"phase set list\">\n" +
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
		if iss.Field == "header/FORMAT/PSL" {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected reserved-key warning for PSL, got issues: %+v", issues)
	}
}
