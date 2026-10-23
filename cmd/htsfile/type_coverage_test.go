package main

import (
	"bytes"
	"compress/flate"
	"encoding/binary"
	"encoding/hex"
	"errors"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
)

// typeCoveragePath is a VCF that exercises the INFO and FORMAT
// types the spec's own worked example does not reach: String
// (single and multi-element), Character, Float, Number=R / G / .
// variable cardinality, a String longer than the 15-byte
// descriptor nibble, and - the important one - missing values in
// every position.
const typeCoveragePath = "../../vcf/testdata/type_coverage.vcf"

// The reference bytes below were produced by `bcftools view -Ob`
// (htslib 1.24) from typeCoveragePath. They are the oracle's
// answer for this file, not ours.
//
// Record 0 covers the present-value case; record 1 covers the
// missing-value case, where every INFO and FORMAT value is
// absent. Getting record 1 right is what caught the defects
// listed in SPEC_COMPLIANCE.md T6.24 - a missing INFO value was
// written as 0, which silently turned `DBL=.,.` into `DBL=0,0`.
const (
	refCoverage0Shared = "0000000063000000010000000000f04108000200020000040717471741110011015768656c6c6f110257612c622c63110317781104f711246162636465666768696a6b6c6d6e6f707172737475767778797a303132333435363738391105250000c03f000020401106210102110731010203110831070809"
	refCoverage0Indiv  = "1109150000c03f00002040110a2778007900110b6761612c62620063632c646400110c2101020304"
	refCoverage1Shared = "00000000c7000000010000000000f0410800030002000004071747373c413e373c443e1100110107110207110307110417781105250100807f0100807f110611801107118011081180"
	refCoverage1Indiv  = "1109150100807f0100807f110a272e002e00110b272e002e00110c2180808080"
)

// TestTypeCoverageMatchesReference asserts that the records we
// write for type_coverage.vcf are byte-identical to the ones
// bcftools writes.
//
// This is the guard for the type coverage the spec's example
// misses. The v0.23.0 release passed every test in the suite
// while writing `DBL=.,.` as `DBL=0,0` and dropping the commas
// from a multi-element FORMAT String.
func TestTypeCoverageMatchesReference(t *testing.T) {
	if _, err := os.Stat(typeCoveragePath); err != nil {
		t.Skipf("type coverage fixture not available: %v", err)
	}
	bcfPath := filepath.Join(t.TempDir(), "cover.bcf")
	if err := translateToFile(typeCoveragePath, bcfPath, runVCFToBCF); err != nil {
		t.Fatalf("vcf2bcf: %v", err)
	}
	got := bcfRecords(t, bcfPath)

	want := []struct{ shared, indiv string }{
		{refCoverage0Shared, refCoverage0Indiv},
		{refCoverage1Shared, refCoverage1Indiv},
	}
	if len(got) != len(want) {
		t.Fatalf("record count: got %d, want %d", len(got), len(want))
	}
	for i, w := range want {
		if got[i].shared != w.shared {
			t.Errorf("record %d shared block differs from bcftools.\n got: %s\nwant: %s",
				i, got[i].shared, w.shared)
		}
		if got[i].indiv != w.indiv {
			t.Errorf("record %d per-sample block differs from bcftools.\n got: %s\nwant: %s",
				i, got[i].indiv, w.indiv)
		}
	}
}

// bcfRecords decompresses a BCF file and returns its records as
// hex strings: the shared block and the per-sample block of each.
func bcfRecords(t *testing.T, path string) []struct{ shared, indiv string } {
	t.Helper()
	raw, err := os.ReadFile(path) //nolint:gosec // the path comes from t.TempDir()
	if err != nil {
		t.Fatal(err)
	}
	// The file is a sequence of BGZF blocks; inflate them all.
	var plain bytes.Buffer
	br := bgzf.NewReader(bytes.NewReader(raw))
	if _, err := io.Copy(&plain, br); err != nil && !errors.Is(err, io.EOF) {
		t.Fatalf("inflate: %v", err)
	}
	_ = flate.BestSpeed
	b := plain.Bytes()
	if len(b) < 9 {
		t.Fatalf("short BCF: %d bytes", len(b))
	}
	lText := binary.LittleEndian.Uint32(b[5:9])
	p := int(9) + int(lText)
	var out []struct{ shared, indiv string }
	for p+8 <= len(b) {
		lShared := binary.LittleEndian.Uint32(b[p : p+4])
		lIndiv := binary.LittleEndian.Uint32(b[p+4 : p+8])
		end := p + 8 + int(lShared) + int(lIndiv)
		if end > len(b) {
			break
		}
		out = append(out, struct{ shared, indiv string }{
			shared: hex.EncodeToString(b[p+8 : p+8+int(lShared)]),
			indiv:  hex.EncodeToString(b[p+8+int(lShared) : end]),
		})
		p = end
	}
	return out
}
