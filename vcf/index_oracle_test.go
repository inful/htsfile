package vcf

import (
	"bytes"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inful/htsfile/csi"
	"github.com/inful/htsfile/tabix"
)

// The oracle test: htslib's tabix and ours must produce the
// same index and answer the same region queries.
//
// This is the check that caught the defects in SPEC_COMPLIANCE.md
// phase 7. Before it, our index had no reference names and no
// VCF preset, was written uncompressed (which htslib rejects
// outright for a .csi), used the tabix bin numbering for a CSI,
// carried a linear index that htslib's CSI does not have, and
// pointed its chunks at the wrong offsets - the last of which
// because our BGZF reader reported a stale virtual offset after
// the first block.
//
// The test is skipped when htslib's tools are not installed.

var oracleRegions = []string{
	"20",
	"20:1-100000",
	"20:14370-14370",
	"20:14000-17400",
	"20:1110696-1110696",
	"20:1230237-1230237",
	"20:1234567-1234567",
	"20:17330-17330",
	"20:1100000-1120000",
	"20:14370-1234567",
}

// oracleTool returns the path of an htslib tool, skipping the
// test when it is not installed.
func oracleTool(t *testing.T, name string) string {
	t.Helper()
	path, err := exec.LookPath(name)
	if err != nil {
		t.Skipf("%s not installed: %v", name, err)
	}
	return path
}

// oracleFixture writes the spec example, bgzipped by htslib, and
// returns the path.
func oracleFixture(t *testing.T, bgzip string) string {
	t.Helper()
	src, err := os.ReadFile("testdata/spec_example.vcf")
	if err != nil {
		t.Fatalf("read spec example: %v", err)
	}
	dir := t.TempDir()
	path := filepath.Join(dir, "spec.vcf.gz")
	cmd := exec.Command(bgzip, "-c")
	cmd.Stdin = bytes.NewReader(src)
	out, err := cmd.Output() //nolint:gosec // the arguments are the test's own paths
	if err != nil {
		t.Fatalf("bgzip: %v", err)
	}
	if err := os.WriteFile(path, out, 0o600); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	return path
}

// oracleIndex runs htslib's tabix over path and returns the
// index path.
func oracleIndex(t *testing.T, tabixBin, path, ext string, csi bool) string {
	t.Helper()
	args := []string{"-f", "-p", "vcf"}
	if csi {
		args = append(args, "-C")
	}
	args = append(args, path)
	cmd := exec.Command(tabixBin, args...)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("htslib tabix %v: %v: %s", args, err, out)
	}
	return path + ext
}

// sameChunks compares two chunk lists by content, not by order.
func sameChunks(a, b []tabix.Chunk) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// TestOracleTabixIndexMatchesHTSlib checks that our .tbi holds
// the same bins, chunks and linear index as htslib's, and that
// htslib's tabix answers the same region queries with either
// index.
func TestOracleTabixIndexMatchesHTSlib(t *testing.T) {
	tabixBin := oracleTool(t, "tabix")
	bgzipBin := oracleTool(t, "bgzip")

	path := oracleFixture(t, bgzipBin)
	oraclePath := oracleIndex(t, tabixBin, path, ".tbi", false)

	// Our index, over the same file.
	oursPath := filepath.Join(t.TempDir(), "ours.vcf.gz")
	src, err := os.ReadFile(path) //nolint:gosec // the path comes from the test's own temp dir
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	if err := os.WriteFile(oursPath, src, 0o600); err != nil { //nolint:gosec // the path comes from the test's own temp dir
		t.Fatalf("write fixture: %v", err)
	}
	if err := BuildIndex(oursPath, IndexOpts{Format: FormatTabix}); err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}

	oracle, err := os.Open(oraclePath) //nolint:gosec // the path comes from the test's own temp dir
	if err != nil {
		t.Fatalf("open oracle index: %v", err)
	}
	defer func() { _ = oracle.Close() }()
	ours, err := os.Open(oursPath + ".tbi") //nolint:gosec // the path comes from the test's own temp dir
	if err != nil {
		t.Fatalf("open our index: %v", err)
	}
	defer func() { _ = ours.Close() }()

	oi, err := tabix.NewReader(oracle).Read()
	if err != nil {
		t.Fatalf("read oracle index: %v", err)
	}
	mi, err := tabix.NewReader(ours).Read()
	if err != nil {
		t.Fatalf("read our index: %v", err)
	}

	// The header: the VCF preset and the reference names.
	if oi.Format != mi.Format || oi.ColSeq != mi.ColSeq ||
		oi.ColBeg != mi.ColBeg || oi.ColEnd != mi.ColEnd ||
		oi.MetaChar != mi.MetaChar || oi.SkipLines != mi.SkipLines {
		t.Errorf("preset: oracle %+v ours %+v", oi, mi)
	}
	if strings.Join(oi.Names, ",") != strings.Join(mi.Names, ",") {
		t.Errorf("names: oracle %v ours %v", oi.Names, mi.Names)
	}
	if len(oi.References) != len(mi.References) {
		t.Fatalf("references: oracle %d ours %d", len(oi.References), len(mi.References))
	}
	for i := range oi.References {
		or, mr := oi.References[i], mi.References[i]
		if len(or.Bins) != len(mr.Bins) {
			t.Errorf("ref %d: bins: oracle %d ours %d", i, len(or.Bins), len(mr.Bins))
			continue
		}
		for bin, ochunks := range or.Bins {
			mchunks, ok := mr.Bins[bin]
			if !ok {
				t.Errorf("ref %d: bin %d missing from our index", i, bin)
				continue
			}
			if !sameChunks(ochunks, mchunks) {
				t.Errorf("ref %d: bin %d: oracle %v ours %v", i, bin, ochunks, mchunks)
			}
		}
		if len(or.Intervals) != len(mr.Intervals) {
			t.Errorf("ref %d: intervals: oracle %d ours %d", i, len(or.Intervals), len(mr.Intervals))
			continue
		}
		for j := range or.Intervals {
			if or.Intervals[j] != mr.Intervals[j] {
				t.Errorf("ref %d: interval %d: oracle %d ours %d", i, j, or.Intervals[j], mr.Intervals[j])
			}
		}
	}

	// The same region queries, answered with either index.
	oracleDir := filepath.Dir(path)
	oursDir := filepath.Dir(oursPath)
	for _, region := range oracleRegions {
		want := oracleQuery(t, tabixBin, path, region)
		got := oracleQuery(t, tabixBin, oursPath, region)
		if want != got {
			t.Errorf("region %s: htslib's query with our index returned %q, want %q", region, got, want)
		}
	}
	_ = oracleDir
	_ = oursDir
}

// oracleQuery runs htslib's tabix over path and returns the
// records it printed.
func oracleQuery(t *testing.T, tabixBin, path, region string) string {
	t.Helper()
	out, err := exec.Command(tabixBin, path, region).Output() //nolint:gosec // the arguments are the test's own paths
	if err != nil {
		// tabix exits non-zero when a region has no records.
		var exitErr *exec.ExitError
		if errors.As(err, &exitErr) {
			return string(exitErr.Stderr)
		}
		t.Fatalf("tabix %s: %v", region, err)
	}
	return string(out)
}

// TestOracleCSIIndexMatchesHTSlib checks the same for the .csi
// format: the bin numbering, the per-bin loffset, and the
// region queries.
func TestOracleCSIIndexMatchesHTSlib(t *testing.T) {
	tabixBin := oracleTool(t, "tabix")
	bgzipBin := oracleTool(t, "bgzip")

	path := oracleFixture(t, bgzipBin)
	oraclePath := oracleIndex(t, tabixBin, path, ".csi", true)

	oursPath := filepath.Join(t.TempDir(), "ours.vcf.gz")
	src, err := os.ReadFile(path) //nolint:gosec // the path comes from the test's own temp dir
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	if err := os.WriteFile(oursPath, src, 0o600); err != nil { //nolint:gosec // the path comes from the test's own temp dir
		t.Fatalf("write fixture: %v", err)
	}
	if err := BuildIndex(oursPath, IndexOpts{Format: FormatCSI, MinShift: 14, Depth: 5}); err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}

	oracle, err := os.Open(oraclePath) //nolint:gosec // the path comes from the test's own temp dir
	if err != nil {
		t.Fatalf("open oracle index: %v", err)
	}
	defer func() { _ = oracle.Close() }()
	ours, err := os.Open(oursPath + ".csi") //nolint:gosec // the path comes from the test's own temp dir
	if err != nil {
		t.Fatalf("open our index: %v", err)
	}
	defer func() { _ = ours.Close() }()

	oi, err := csi.NewReader(oracle).Read()
	if err != nil {
		t.Fatalf("read oracle index: %v", err)
	}
	mi, err := csi.NewReader(ours).Read()
	if err != nil {
		t.Fatalf("read our index: %v", err)
	}
	if oi.MinShift != mi.MinShift || oi.Depth != mi.Depth {
		t.Errorf("min_shift/depth: oracle %d/%d ours %d/%d", oi.MinShift, oi.Depth, mi.MinShift, mi.Depth)
	}
	if !bytes.Equal(oi.Aux, mi.Aux) {
		t.Errorf("aux: oracle %x ours %x", oi.Aux, mi.Aux)
	}
	if len(oi.References) != len(mi.References) {
		t.Fatalf("references: oracle %d ours %d", len(oi.References), len(mi.References))
	}
	for i := range oi.References {
		or, mr := oi.References[i], mi.References[i]
		if len(or.Bins) != len(mr.Bins) {
			t.Errorf("ref %d: bins: oracle %d ours %d", i, len(or.Bins), len(mr.Bins))
			continue
		}
		for bin, ochunks := range or.Bins {
			mchunks, ok := mr.Bins[bin]
			if !ok {
				t.Errorf("ref %d: bin %d missing from our index", i, bin)
				continue
			}
			if len(ochunks) != len(mchunks) {
				t.Errorf("ref %d: bin %d: oracle %d chunks ours %d", i, bin, len(ochunks), len(mchunks))
				continue
			}
			for j := range ochunks {
				if ochunks[j] != mchunks[j] {
					t.Errorf("ref %d: bin %d chunk %d: oracle %+v ours %+v", i, bin, j, ochunks[j], mchunks[j])
				}
			}
			if or.Loffsets[bin] != mr.Loffsets[bin] {
				t.Errorf("ref %d: bin %d loffset: oracle %d ours %d", i, bin, or.Loffsets[bin], mr.Loffsets[bin])
			}
		}
	}

	for _, region := range oracleRegions {
		want := oracleQueryCSI(t, tabixBin, path, region)
		got := oracleQueryCSI(t, tabixBin, oursPath, region)
		if want != got {
			t.Errorf("region %s: htslib's query with our index returned %q, want %q", region, got, want)
		}
	}
}

// oracleQueryCSI runs htslib's tabix in CSI mode over path.
func oracleQueryCSI(t *testing.T, tabixBin, path, region string) string {
	t.Helper()
	out, err := exec.Command(tabixBin, "-C", path, region).Output() //nolint:gosec // the arguments are the test's own paths
	if err != nil {
		var exitErr *exec.ExitError
		if errors.As(err, &exitErr) {
			return string(exitErr.Stderr)
		}
		t.Fatalf("tabix -C %s: %v", region, err)
	}
	return string(out)
}
