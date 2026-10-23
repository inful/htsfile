package vcf

import (
	"os"
	"strings"
	"testing"
)

// The spec's own gVCF example (VCFv4.5 §1.7, tex lines
// 1862-1874) has a typo: five rows separate the last two
// FORMAT values with a semicolon where a colon belongs, so
// "PL;LEN" arrives as a single value "900;14" that is not an
// Integer. The LEN value is then missing, and the END check
// fails with it, because a reference block's END is POS+LEN-1.
//
// The fixture testdata/gvcf_example.vcf is that example with
// the semicolons corrected; it must validate clean. The
// example as printed must produce exactly the diagnostics the
// typo causes - that is what pins the behaviour, so a future
// change cannot quietly accept the typo.

// gvcfExampleAsPrinted is the spec's example with its typo
// intact: the last FORMAT value of rows 1, 2, 4, 5 and 7 uses
// a semicolon instead of a colon.
var gvcfExampleAsPrinted = strings.ReplaceAll(
	strings.ReplaceAll(
		strings.ReplaceAll(
			strings.ReplaceAll(
				strings.ReplaceAll(gvcfExample, "900:14", "900;14"),
				"630:4", "630;4"),
			"315:1", "315;1"),
		"945:4", "945;4"),
	"593:19", "593;19")

// gvcfExample is the spec's example, read from the fixture.
var gvcfExample = func() string {
	data, err := os.ReadFile("testdata/gvcf_example.vcf")
	if err != nil {
		// The fixture is missing: the tests below fail with a
		// clear message rather than a panic.
		return ""
	}
	return string(data)
}()

// gvcfIssues validates a gVCF file and returns the issues.
func gvcfIssues(t *testing.T, content string) []string {
	t.Helper()
	dec := NewDecoder(strings.NewReader(content))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("header: %v", err)
	}
	v := NewValidator(ModeStrict)
	var issues []string
	for _, is := range v.HeaderChecks(hdr) {
		issues = append(issues, is.Msg)
	}
	for {
		rec, err := dec.Record()
		if err != nil {
			break
		}
		if rec == nil {
			break
		}
		for _, is := range v.Next(hdr, rec) {
			issues = append(issues, is.Msg)
		}
	}
	return issues
}

// TestGVCFSpecExampleIsClean checks that the spec's gVCF
// example, with its semicolon typos corrected, validates
// without a single issue.
func TestGVCFSpecExampleIsClean(t *testing.T) {
	if gvcfExample == "" {
		t.Fatal("testdata/gvcf_example.vcf is missing")
	}
	if issues := gvcfIssues(t, gvcfExample); len(issues) != 0 {
		t.Errorf("the spec's gVCF example should validate clean, got %d issue(s):", len(issues))
		for _, is := range issues {
			t.Errorf("  %s", is)
		}
	}
}

// TestGVCFSpecExampleAsPrinted checks that the example exactly
// as the spec prints it produces the diagnostics its typos
// cause: five FORMAT values that are not Integers, and the four
// END checks that fail because the LEN value they need never
// arrived.
func TestGVCFSpecExampleAsPrinted(t *testing.T) {
	if gvcfExample == "" {
		t.Fatal("testdata/gvcf_example.vcf is missing")
	}
	issues := gvcfIssues(t, gvcfExampleAsPrinted)
	if len(issues) != 9 {
		t.Fatalf("the spec's example as printed should produce 9 issues, got %d:", len(issues))
		for _, is := range issues {
			t.Errorf("  %s", is)
		}
	}
	var notInteger, endMismatch int
	for _, is := range issues {
		switch {
		case strings.Contains(is, "does not match declared Type=\"Integer\""):
			notInteger++
		case strings.Contains(is, "declared END"):
			endMismatch++
		}
	}
	if notInteger != 5 {
		t.Errorf("FORMAT values that are not Integers: got %d want 5", notInteger)
	}
	if endMismatch != 4 {
		t.Errorf("END mismatches: got %d want 4", endMismatch)
	}
}
