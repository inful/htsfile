package vcf

import (
	"strings"
	"testing"
)

// minimalVCFWithID is a small VCF body used by the streaming
// tests below. The first record carries an explicit ID
// ("rs1") so cross-record uniqueness can be exercised.
const minimalVCFWithID = "##fileformat=VCFv4.5\n" +
	"##contig=<ID=chr1,length=1000>\n" +
	"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
	"chr1\t100\trs1\tA\tC\t30\tPASS\tDP=42\n" +
	"chr1\t200\trs2\tA\tC\t30\tPASS\tDP=42\n"

// TestValidatorResetClearsState locks in the doc-comment
// promise from Validator.Reset: a Validator that runs two
// files through the streaming API (HeaderChecks / Next /
// Finish) without calling Reset() in between will leak the
// first file's seenIDs into the second file's first
// record. After Reset(), the second file's first record
// passes cleanly because the cross-record ID set has been
// cleared.
func TestValidatorResetClearsState(t *testing.T) {
	// Each "file" is the same minimal VCF, but the cross-
	// record ID set should still be reset between files
	// when the consumer says so.
	dec := NewDecoder(strings.NewReader(minimalVCFWithID))
	h1, err := dec.Header()
	if err != nil {
		t.Fatalf("first header: %v", err)
	}

	v := NewValidator(ModeDefault)
	if issues := v.HeaderChecks(h1); len(issues) != 0 {
		t.Fatalf("first HeaderChecks issues: %+v", issues)
	}

	// Drain both records of the first file via the
	// streaming API. Two records on one decoder.
	r1, err := dec.Record()
	if err != nil {
		t.Fatalf("first file record 1: %v", err)
	}
	if issues := v.Next(h1, r1); len(issues) != 0 {
		t.Fatalf("Next r1: %+v", issues)
	}
	r2, err := dec.Record()
	if err != nil {
		t.Fatalf("first file record 2: %v", err)
	}
	if issues := v.Next(h1, r2); len(issues) != 0 {
		t.Fatalf("Next r2 (second record of first file): %+v", issues)
	}
	if issues := v.Finish(); len(issues) != 0 {
		t.Fatalf("Finish: %+v", issues)
	}

	// Now feed the same file content to the same
	// Validator without calling Reset. The first record's
	// ID "rs1" is already in seenIDs from the first file,
	// so it MUST be flagged as a duplicate — that's the
	// footgun.
	h2, err := NewDecoder(strings.NewReader(minimalVCFWithID)).Header()
	if err != nil {
		t.Fatalf("second header (post-skip): %v", err)
	}
	if issues := v.HeaderChecks(h2); len(issues) != 0 {
		t.Fatalf("second HeaderChecks (no Reset): %+v", issues)
	}
	dec3 := NewDecoder(strings.NewReader(minimalVCFWithID))
	if _, err := dec3.Header(); err != nil {
		t.Fatalf("third header: %v", err)
	}
	rb, err := dec3.Record()
	if err != nil {
		t.Fatalf("second file record: %v", err)
	}
	issues := v.Next(h2, rb)
	sawDuplicate := false
	for _, iss := range issues {
		if strings.Contains(iss.Msg, "ID") && strings.Contains(iss.Msg, "already used") {
			sawDuplicate = true
			break
		}
	}
	if !sawDuplicate {
		t.Fatalf("expected cross-record duplicate ID without Reset, got %+v", issues)
	}

	// Calling Reset wipes the streaming state — the very
	// next pass through the same file should produce no
	// duplicate-id issues on its first record.
	v.Reset()
	if issues := v.HeaderChecks(h2); len(issues) != 0 {
		t.Fatalf("post-Reset HeaderChecks: %+v", issues)
	}
	if issues := v.Finish(); len(issues) != 0 {
		t.Fatalf("post-Reset Finish: %+v", issues)
	}
	dec4 := NewDecoder(strings.NewReader(minimalVCFWithID))
	if _, err := dec4.Header(); err != nil {
		t.Fatalf("fourth header: %v", err)
	}
	rb2, err := dec4.Record()
	if err != nil {
		t.Fatalf("third file record: %v", err)
	}
	if issues := v.Next(h2, rb2); len(issues) != 0 {
		t.Fatalf("post-Reset Next: %+v", issues)
	}
}
