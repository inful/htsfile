package main

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// specExamplePath is the VCFv4.5 spec's own worked example
// (VCFv4.5 tex section 1.6.1, lines 43-68), checked in as
// testdata. It is the only ground truth we have that does not
// come from this package's own encoder.
const specExamplePath = "../../vcf/testdata/spec_example.vcf"

// specExampleRoundTripPath is the expected result of
// vcf2bcf | bcf2vcf on the spec example.
//
// Two differences from the input are expected and spec-allowed:
//
//   - INFO keys are sorted. VCFv4.5 section 1.4.3 line 133:
//     "Implementations must not rely on the order of the fields
//     within structured lines and are not required to preserve
//     field ordering."
//   - `md5` on ##contig gains quotes. VCFv4.5 section 1.4.3
//     line 126 requires optional structured-line field values to
//     be quoted, so the quoted form is *more* spec-compliant
//     than the spec's own example.
//
// Everything else - every record, every column, every
// per-sample value including the "1/1:43:5:.,." record with
// its list-of-missing-values - must be identical.
const specExampleRoundTripPath = "../../vcf/testdata/spec_example.roundtrip.vcf"

// TestSpecExampleStrictCheck runs `check --strict` on the spec's
// own worked example. It must report zero diagnostics.
func TestSpecExampleStrictCheck(t *testing.T) {
	if _, err := os.Stat(specExamplePath); err != nil {
		t.Skipf("spec example not available: %v", err)
	}
	out, err := captureStdout(t, func() error { return runCheck([]string{"--strict", specExamplePath}) })
	if err != nil {
		t.Errorf("check --strict on the spec example: %v\n%s", err, out)
	}
	if !strings.Contains(out, "clean") {
		t.Errorf("check --strict: expected a clean report, got:\n%s", out)
	}
}

// TestSpecExampleBCFRoundTrip runs vcf2bcf | bcf2vcf on the
// spec example and compares the result field-for-field with the
// checked-in expectation.
//
// This is the regression guard for the whole BCF layer: the
// record framing, the dictionary indices, the type-byte
// encoding, the GT allele vector and the per-slot missing
// values all have to survive the trip.
func TestSpecExampleBCFRoundTrip(t *testing.T) {
	if _, err := os.Stat(specExamplePath); err != nil {
		t.Skipf("spec example not available: %v", err)
	}
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "spec.bcf")
	if err := translateToFile(specExamplePath, bcfPath, runVCFToBCF); err != nil {
		t.Fatalf("vcf2bcf: %v", err)
	}
	// The BCF must be readable by our own checker.
	out, err := captureStdout(t, func() error { return runCheck([]string{bcfPath}) })
	if err != nil {
		t.Fatalf("check on the produced BCF: %v\n%s", err, out)
	}
	if !strings.Contains(out, "clean") {
		t.Errorf("check on the produced BCF: expected a clean report, got:\n%s", out)
	}
	// Round-trip back to VCF.
	rtPath := filepath.Join(dir, "rt.vcf")
	if err := translateToFile(bcfPath, rtPath, runBCFToVCF); err != nil {
		t.Fatalf("bcf2vcf: %v", err)
	}
	got, err := os.ReadFile(rtPath)
	if err != nil {
		t.Fatal(err)
	}
	want, err := os.ReadFile(specExampleRoundTripPath)
	if err != nil {
		t.Skipf("round-trip expectation not available: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Errorf("round-trip mismatch.\n--- got ---\n%s\n--- want ---\n%s", got, want)
	}
}

// TestSpecExampleRoundTripPreservesData compares only the data
// rows, ignoring the header. This pins the property that matters
// most: no record data is lost or altered.
func TestSpecExampleRoundTripPreservesData(t *testing.T) {
	if _, err := os.Stat(specExamplePath); err != nil {
		t.Skipf("spec example not available: %v", err)
	}
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "spec.bcf")
	if err := translateToFile(specExamplePath, bcfPath, runVCFToBCF); err != nil {
		t.Fatalf("vcf2bcf: %v", err)
	}
	rtPath := filepath.Join(dir, "rt.vcf")
	if err := translateToFile(bcfPath, rtPath, runBCFToVCF); err != nil {
		t.Fatalf("bcf2vcf: %v", err)
	}
	in := dataRows(t, specExamplePath)
	out := dataRows(t, rtPath)
	if len(in) != len(out) {
		t.Fatalf("record count: got %d, want %d", len(out), len(in))
	}
	for i := range in {
		gotFields := strings.Split(out[i], "\t")
		wantFields := strings.Split(in[i], "\t")
		if len(gotFields) != len(wantFields) {
			t.Errorf("record %d: %d columns, want %d", i+1, len(gotFields), len(wantFields))
			continue
		}
		for j := range gotFields {
			if j == 7 { // INFO column
				if !sameInfoSet(gotFields[7], wantFields[7]) {
					t.Errorf("record %d INFO: got %q, want %q", i+1, gotFields[7], wantFields[7])
				}
				continue
			}
			if j >= 9 {
				// A trailing FORMAT field that is absent in
				// the input may be written as the canonical
				// "." (VCFv4.5 section 1.6.2 line 508) or
				// dropped entirely. Compare them as equal.
				if trimTrailingMissing(gotFields[j]) != trimTrailingMissing(wantFields[j]) {
					t.Errorf("record %d column %d: got %q, want %q", i+1, j+1, gotFields[j], wantFields[j])
				}
				continue
			}
			if gotFields[j] != wantFields[j] {
				t.Errorf("record %d column %d: got %q, want %q", i+1, j+1, gotFields[j], wantFields[j])
			}
		}
	}
}

// translateToFile runs one of the translation functions with the
// input file as the reader and the output file as the writer.
func translateToFile(inPath, outPath string, fn func(io.Reader, io.Writer) error) error {
	in, err := os.Open(inPath)
	if err != nil {
		return err
	}
	defer func() { _ = in.Close() }()
	out, err := os.Create(outPath)
	if err != nil {
		return err
	}
	defer func() { _ = out.Close() }()
	return fn(in, out)
}

// dataRows returns the non-header, non-comment lines of a VCF.
func dataRows(t *testing.T, path string) []string {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var out []string
	for _, line := range strings.Split(strings.TrimRight(string(b), "\n"), "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		out = append(out, line)
	}
	return out
}

// sameInfoSet compares two INFO columns as sets of key=value
// pairs, since the field order is not preserved (and is not
// required to be).
func sameInfoSet(a, b string) bool {
	if a == b {
		return true
	}
	as := strings.Split(a, ";")
	bs := strings.Split(b, ";")
	if len(as) != len(bs) {
		return false
	}
	sortStrings(as)
	sortStrings(bs)
	for i := range as {
		if as[i] != bs[i] {
			return false
		}
	}
	return true
}

// sortStrings sorts in place (insertion sort; the lists are tiny).
func sortStrings(s []string) {
	for i := 1; i < len(s); i++ {
		for j := i; j > 0 && s[j] < s[j-1]; j-- {
			s[j], s[j-1] = s[j-1], s[j]
		}
	}
}

// trimTrailingMissing removes trailing FORMAT sub-fields that
// carry no value ("." or ""). A field that is absent in the
// input may be written as the canonical "." or dropped
// entirely; both mean the same thing.
func trimTrailingMissing(field string) string {
	parts := strings.Split(field, ":")
	for len(parts) > 1 {
		last := parts[len(parts)-1]
		if last != "." && last != "" {
			break
		}
		parts = parts[:len(parts)-1]
	}
	return strings.Join(parts, ":")
}
