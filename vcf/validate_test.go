package vcf

import (
	"strings"
	"testing"
)

func TestValidator_CleanSpecExample(t *testing.T) {
	t.Parallel()
	hdr, recs := parseFile(t, "testdata/spec_example.vcf")
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) > 0 {
		t.Fatalf("expected no issues, got %d:\n%v", len(issues), issues)
	}
}

func TestValidator_CleanSpecExampleStrict(t *testing.T) {
	t.Parallel()
	hdr, recs := parseFile(t, "testdata/spec_example.vcf")
	v := NewValidator(ModeStrict)
	issues := v.Validate(hdr, recs)
	if len(issues) > 0 {
		t.Fatalf("expected no issues in strict mode, got %d:\n%v", len(issues), issues)
	}
}

func TestValidator_CleanSVExample(t *testing.T) {
	t.Parallel()
	hdr, recs := parseFile(t, "testdata/sv_example.vcf")
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) > 0 {
		t.Fatalf("expected no issues, got %d:\n%v", len(issues), issues)
	}
}

func TestValidator_UndeclaredInfoKey(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t2\t.\tA\tC\t.\t.\tZZ=42;DP=10\n"
	dec := NewDecoder(stringReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatal(err)
	}
	var recs []*Record
	for {
		r, err := dec.Record()
		if err != nil {
			break
		}
		recs = append(recs, r)
	}
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 issue, got %d:\n%v", len(issues), issues)
	}
	if issues[0].Severity != SeverityWarning {
		t.Errorf("default mode should warn, got %v", issues[0].Severity)
	}
	if issues[0].Field != "INFO/ZZ" {
		t.Errorf("Field = %q", issues[0].Field)
	}

	// Strict mode should error.
	v = NewValidator(ModeStrict)
	issues = v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 issue, got %d", len(issues))
	}
	if issues[0].Severity != SeverityError {
		t.Errorf("strict mode should error, got %v", issues[0].Severity)
	}
}

func TestValidator_NumberCardinalityMismatch(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"freq\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t2\t.\tA\tC,G\t.\t.\tAF=0.1\n" // AF=0.1 but ALT has 2 alleles
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 issue, got %d:\n%v", len(issues), issues)
	}
	if issues[0].Field != "INFO/AF" {
		t.Errorf("Field = %q", issues[0].Field)
	}
}

func TestValidator_InfoTypeMismatch(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t2\t.\tA\tC\t.\t.\tDP=oops\n"
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 issue, got %d:\n%v", len(issues), issues)
	}
	if issues[0].Field != "INFO/DP" {
		t.Errorf("Field = %q", issues[0].Field)
	}
}

func TestValidator_FilterNotDeclared(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##FILTER=<ID=PASS,Description=\"ok\">\n" +
		"##FILTER=<ID=q10,Description=\"q<10\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t2\t.\tA\tC\t.\tq20\t.\n"
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 issue, got %d", len(issues))
	}
	if issues[0].Field != "FILTER" {
		t.Errorf("Field = %q", issues[0].Field)
	}
	if issues[0].Severity != SeverityWarning {
		t.Errorf("default mode warning, got %v", issues[0].Severity)
	}
}

func TestValidator_PASSAccepted(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t2\t.\tA\tC\t.\tPASS\t.\n"
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 0 {
		t.Fatalf("PASS should not be flagged; got %d issues", len(issues))
	}
}

func TestValidator_GTNotFirst(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" +
		"1\t2\t.\tA\tC\t.\t.\t.\tDP:GT\t5:0/1\n"
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) == 0 {
		t.Fatal("expected GT-not-first issue")
	}
	found := false
	for _, i := range issues {
		if i.Field == "FORMAT" {
			found = true
		}
	}
	if !found {
		t.Errorf("expected FORMAT issue; got %v", issues)
	}
}

func TestValidator_UndeclaredContig(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr2\t2\t.\tA\tC\t.\t.\t.\n"
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 issue, got %d", len(issues))
	}
	if issues[0].Severity != SeverityWarning {
		t.Errorf("default mode warning; got %v", issues[0].Severity)
	}

	v = NewValidator(ModeStrict)
	issues = v.Validate(hdr, recs)
	if issues[0].Severity != SeverityError {
		t.Errorf("strict mode error; got %v", issues[0].Severity)
	}
}

func TestValidator_StrictRefBases(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t2\t.\tACGX\tC\t.\t.\t.\n"
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeStrict)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("strict mode should reject non-ACGTN REF, got %d issues:\n%v",
			len(issues), issues)
	}
}

func TestValidator_StrictPOSMonotonic(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t100\t.\tA\tC\t.\t.\t.\n" +
		"1\t50\t.\tA\tC\t.\t.\t.\n" +
		"2\t10\t.\tA\tC\t.\t.\t.\n" +
		"2\t20\t.\tA\tC\t.\t.\t.\n"
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeStrict)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 issue, got %d:\n%v", len(issues), issues)
	}
	if issues[0].Field != "POS" {
		t.Errorf("Field = %q", issues[0].Field)
	}
}

func TestValidator_StrictALTEqualsRef(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t2\t.\tA\tA\t.\t.\t.\n"
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeStrict)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 issue, got %d:\n%v", len(issues), issues)
	}
}

func TestValidator_FormatCardinalityG(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"pl\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" +
		"1\t2\t.\tA\tC,G\t.\t.\t.\tGT:PL\t0/1:1,2,3,4,5\n" // 5 values, expected 6 for diploid biallelic
	hdr, recs := parseFromString(t, src)
	v := NewValidator(ModeDefault)
	issues := v.Validate(hdr, recs)
	if len(issues) != 1 {
		t.Fatalf("expected 1 cardinality issue, got %d:\n%v", len(issues), issues)
	}
	if issues[0].Field != "FORMAT/PL/sample[0]" {
		t.Errorf("Field = %q", issues[0].Field)
	}
}

// helpers

func parseFromString(t *testing.T, src string) (*Header, []*Record) {
	t.Helper()
	dec := NewDecoder(stringReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	var recs []*Record
	for {
		r, err := dec.Record()
		if err != nil {
			break
		}
		recs = append(recs, r)
	}
	return hdr, recs
}

// stringReader returns an io.Reader over the given string.
func stringReader(s string) *strings.Reader { return strings.NewReader(s) }
