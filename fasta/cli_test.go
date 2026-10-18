package fasta

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// TestCLI_CheckAndSummary builds the cmd/htsfile binary and
// runs its `check`, `roundtrip`, and `summary` subcommands
// against a synthetic FASTA file.
func TestCLI_CheckAndSummary(t *testing.T) {
	root := repoRoot(t)
	dir := t.TempDir()
	bin := filepath.Join(dir, "htsfile")
	cmd := exec.Command("go", "build", "-o", bin, "./cmd/htsfile")
	cmd.Dir = root
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build: %v\n%s", err, out)
	}

	src := ">seq1 E. coli K-12\n" +
		"ACGTACGT\n" +
		"ACGTACGT\n" +
		">seq2 some other description\n" +
		"GATTACA\n"
	faPath := filepath.Join(dir, "smoke.fa")
	if err := os.WriteFile(faPath, []byte(src), 0o644); err != nil {
		t.Fatal(err)
	}

	t.Run("check", func(t *testing.T) {
		out, err := exec.Command(bin, "check", faPath).CombinedOutput()
		if err != nil {
			t.Fatalf("check failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "clean") {
			t.Errorf("check output missing 'clean': %s", out)
		}
	})

	t.Run("summary", func(t *testing.T) {
		out, err := exec.Command(bin, "summary", faPath).CombinedOutput()
		if err != nil {
			t.Fatalf("summary failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "records:") {
			t.Errorf("summary output missing 'records:': %s", out)
		}
		if !strings.Contains(string(out), "bases:") {
			t.Errorf("summary output missing 'bases:': %s", out)
		}
	})

	t.Run("roundtrip", func(t *testing.T) {
		out, err := exec.Command(bin, "roundtrip", faPath).CombinedOutput()
		if err != nil {
			t.Fatalf("roundtrip failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), ">seq1 E. coli K-12") {
			t.Errorf("roundtrip missing header:\n%s", out)
		}
		if !strings.Contains(string(out), "ACGTACGTACGTACGT") {
			t.Errorf("roundtrip missing sequence:\n%s", out)
		}
	})
}

// repoRoot finds the htsfile repo root by walking up from the
// current directory until it sees go.mod.
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