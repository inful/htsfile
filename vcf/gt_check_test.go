package vcf

import (
	"strings"
	"testing"
)

func TestCheckGT_Valid(t *testing.T) {
	// nAlts=2: indices 0 (REF), 1, 2 are valid; "." is the
	// missing-allele marker.
	for _, gt := range []string{
		"0/1", "0|1", "1/2", "0/0", "0/1/2", "0|1|2",
		"./.", "./.", "0/.", "0/.",
	} {
		if err := checkGT(gt, 2); err != nil {
			t.Errorf("checkGT(%q, 2): unexpected err: %v", gt, err)
		}
	}
}

func TestCheckGT_MixedSeparators(t *testing.T) {
	// Mixing / and | within a single GT is forbidden.
	for _, gt := range []string{
		"0/1|2", "0|1/2", "0/1|2/3", "0|.|1/2",
	} {
		err := checkGT(gt, 3)
		if err == nil {
			t.Errorf("checkGT(%q, 3): expected error, got nil", gt)
			continue
		}
		if !strings.Contains(err.Error(), "mixes") {
			t.Errorf("checkGT(%q, 3): error %q should mention mixed separators", gt, err)
		}
	}
}

func TestCheckGT_OutOfRange(t *testing.T) {
	// nAlts=2: index 3 or higher is out of range.
	for _, gt := range []string{
		"0/3", "3|0", "0/1/4", "9/0",
	} {
		err := checkGT(gt, 2)
		if err == nil {
			t.Errorf("checkGT(%q, 2): expected error, got nil", gt)
			continue
		}
		if !strings.Contains(err.Error(), "out of range") {
			t.Errorf("checkGT(%q, 2): error %q should mention out of range", gt, err)
		}
	}
}

func TestCheckGT_NonInteger(t *testing.T) {
	// Non-integer tokens are not valid allele indices.
	for _, gt := range []string{
		"X/Y", "0/a", "x|y", "0/1.5",
	} {
		err := checkGT(gt, 2)
		if err == nil {
			t.Errorf("checkGT(%q, 2): expected error, got nil", gt)
			continue
		}
		if !strings.Contains(err.Error(), "not an integer") {
			t.Errorf("checkGT(%q, 2): error %q should mention not an integer", gt, err)
		}
	}
}

func TestCheckGT_Empty(t *testing.T) {
	// An empty GT is treated as "no value" and passes.
	if err := checkGT("", 2); err != nil {
		t.Errorf("checkGT(\"\", 2): unexpected err: %v", err)
	}
}

func TestCheckGT_ZeroAlts(t *testing.T) {
	// With no ALTs, the only valid allele is 0 (REF).
	if err := checkGT("0/0", 0); err != nil {
		t.Errorf("checkGT(\"0/0\", 0): unexpected err: %v", err)
	}
	if err := checkGT("0/1", 0); err == nil {
		t.Error("checkGT(\"0/1\", 0): expected error (1 is out of range)")
	}
}

func TestValidator_GTIssuesInCheckSamples(t *testing.T) {
	// Wire check: a record with a mixed-separator GT surfaces
	// an error from Validator.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	dec := NewDecoder(strings.NewReader(hdrText))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	// Build a record by hand with a mixed-separator GT.
	rec := &Record{
		Line: 1, Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"C", "G"},
		Format: []string{"GT"},
		Samples: []Genotype{{Fields: map[string]FieldValue{
			"GT": {typedCore: typedCore{Type: TypeString, Strs: []string{"0/1|2"}}},
		}}},
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "FORMAT/GT/sample[0]" && strings.Contains(iss.Msg, "mixes") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected GT-mixed-separator issue, got: %+v", issues)
	}
}
