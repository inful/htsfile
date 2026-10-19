package main

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
)

// writeVCFQueryFixture is the shared helper used by the
// query and index CLI tests. It writes a small VCF plus
// the .tbi that index-by-tabix produces, and is reused by
// index tests that exercise the round-trip.
func makeVCFIndexFixture(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "sample.vcf.gz")

	const content = "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=10000>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t.\t.\t.\n" +
		"chr1\t300\t.\tC\tT\t.\t.\t.\n" +
		"chr1\t500\t.\tG\tA\t.\t.\t.\n"

	f, err := os.Create(vcfPath) //nolint:gosec // test fixture in t.TempDir
	if err != nil {
		t.Fatalf("create vcf: %v", err)
	}
	bw := bgzf.NewWriter(f)
	if _, err := bw.Write([]byte(content)); err != nil {
		t.Fatalf("bw.Write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bw.Close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("f.Close: %v", err)
	}
	return vcfPath
}

func TestRunVCFIndex_TabixBuilds(t *testing.T) {
	vcfPath := makeVCFIndexFixture(t)
	if err := runVCFIndex([]string{vcfPath}); err != nil {
		t.Fatalf("runVCFIndex: %v", err)
	}
	if _, err := os.Stat(vcfPath + ".tbi"); err != nil {
		t.Fatalf("expected .tbi next to %s: %v", vcfPath, err)
	}
}

func TestRunVCFIndex_CSIBuilds(t *testing.T) {
	vcfPath := makeVCFIndexFixture(t)
	if err := runVCFIndex([]string{"--format=csi", vcfPath}); err != nil {
		t.Fatalf("runVCFIndex: %v", err)
	}
	if _, err := os.Stat(vcfPath + ".csi"); err != nil {
		t.Fatalf("expected .csi next to %s: %v", vcfPath, err)
	}
}

func TestRunVCFIndex_RoundTripWithQuery(t *testing.T) {
	vcfPath := makeVCFIndexFixture(t)
	if err := runVCFIndex([]string{vcfPath}); err != nil {
		t.Fatalf("runVCFIndex: %v", err)
	}
	if err := runVCFQuery([]string{vcfPath, "chr1:200-400"}); err != nil {
		t.Fatalf("runVCFQuery: %v", err)
	}
}

func TestRunVCFIndex_RejectsBadArgs(t *testing.T) {
	if err := runVCFIndex([]string{}); err == nil {
		t.Error("expected error for no args")
	}
	if err := runVCFIndex([]string{"--format=bogus", "foo.vcf.gz"}); err == nil {
		t.Error("expected error for unknown format")
	}
	if err := runVCFIndex([]string{"a", "b"}); err == nil {
		t.Error("expected error for too many args")
	}
}