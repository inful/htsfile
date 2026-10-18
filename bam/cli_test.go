package bam

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// repoRoot returns the absolute path to the htsfile repo root
// (the directory containing this package's parent's go.mod).
// We locate it by walking up from the current working directory
// until we find go.mod. The CLI build needs a stable anchor
// because go build resolves relative paths against the
// current directory, and the package's tests run with their
// own package directory as CWD.
func repoRoot(t *testing.T) string {
	t.Helper()
	wd, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	dir := wd
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("go.mod not found above %s", wd)
		}
		dir = parent
	}
}

// TestCLI_CheckAndSummary builds the cmd/htsfile binary and
// runs its `check` and `summary` subcommands against a
// synthetic BAM written to disk. The CLI's `roundtrip` for BAM
// is verified separately to return a clear "not supported"
// error.
func TestCLI_CheckAndSummary(t *testing.T) {
	// Build the CLI binary once.
	root := repoRoot(t)
	dir := t.TempDir()
	bin := filepath.Join(dir, "htsfile")
	cmd := exec.Command("go", "build", "-o", bin, "./cmd/htsfile")
	cmd.Dir = root
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build: %v\n%s", err, out)
	}

	// Build a synthetic BAM.
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n", []struct {
		name   string
		length int64
	}{{"chr1", 248956422}}); err != nil {
		t.Fatal(err)
	}
	if err := b.writeRecord(&Record{
		QName:     "smoke1",
		RefID:     0,
		Pos:       100,
		MapQ:      60,
		Bin:       4681,
		Cigar:     []uint32{cigar(CigarM, 5)}, // 5M
		Seq:       "ACGTA",
		Qual:      []byte("IIIII"),
		NextRefID: -1,
	}); err != nil {
		t.Fatal(err)
	}
	bamPath := filepath.Join(dir, "smoke.bam")
	if err := os.WriteFile(bamPath, b.finish(), 0o644); err != nil {
		t.Fatal(err)
	}

	// Run the CLI's check subcommand.
	t.Run("check", func(t *testing.T) {
		out, err := exec.Command(bin, "check", bamPath).CombinedOutput()
		if err != nil {
			t.Fatalf("check failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "clean") {
			t.Errorf("check output missing 'clean': %s", out)
		}
	})

	// Run the CLI's summary subcommand.
	t.Run("summary", func(t *testing.T) {
		out, err := exec.Command(bin, "summary", bamPath).CombinedOutput()
		if err != nil {
			t.Fatalf("summary failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "records:") {
			t.Errorf("summary output missing 'records:': %s", out)
		}
	})

	// Run the CLI's roundtrip subcommand — should fail with a
	// clear "BAM v1 is read-only" message, not silently
	// succeed.
	t.Run("roundtrip_unsupported", func(t *testing.T) {
		out, err := exec.Command(bin, "roundtrip", bamPath).CombinedOutput()
		if err == nil {
			t.Fatalf("roundtrip should have failed for BAM, but succeeded:\n%s", out)
		}
		combined := string(out)
		if !strings.Contains(combined, "read-only") {
			t.Errorf("roundtrip error doesn't mention read-only:\n%s", combined)
		}
	})
}

// TestCLI_Version runs the CLI's `version` subcommand (which
// doesn't need a file argument) to confirm the binary builds
// and reports a sensible version string.
func TestCLI_Version(t *testing.T) {
	root := repoRoot(t)
	dir := t.TempDir()
	bin := filepath.Join(dir, "htsfile")
	cmd := exec.Command("go", "build", "-o", bin, "./cmd/htsfile")
	cmd.Dir = root
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build: %v\n%s", err, out)
	}
	out, err := exec.Command(bin, "version").CombinedOutput()
	if err != nil {
		t.Fatalf("version failed: %v\n%s", err, out)
	}
	if !bytes.Contains(out, []byte("htsfile ")) {
		t.Errorf("version output unexpected: %s", out)
	}
}
