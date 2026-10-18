package fastq

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// TestCLI_CheckSummaryRoundTrip builds the cmd/htsfile binary
// and runs its `check`, `roundtrip`, and `summary` subcommands
// against a synthetic FASTQ file. Mirrors the FASTA package's
// cli_test.go.
func TestCLI_CheckSummaryRoundTrip(t *testing.T) {
	if os.Getenv("HT_SKIP_CLI") != "" {
		t.Skip("HT_SKIP_CLI set")
	}
	root := repoRoot(t)
	dir := t.TempDir()
	bin := filepath.Join(dir, "htsfile")
	build := exec.Command("go", "build", "-o", bin, "./cmd/htsfile")
	build.Dir = root
	if out, err := build.CombinedOutput(); err != nil {
		t.Fatalf("go build: %v\n%s", err, out)
	}

	src := "@r1\nACGT\n+\nIIII\n@r2\nGGGG\n+\nJJJJ\n"
	fqPath := filepath.Join(dir, "smoke.fq")
	if err := os.WriteFile(fqPath, []byte(src), 0o600); err != nil {
		t.Fatal(err)
	}

	t.Run("check", func(t *testing.T) {
		out, err := exec.Command(bin, "check", fqPath).CombinedOutput()
		if err != nil {
			t.Fatalf("check failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "clean (2 records, 8 bases)") {
			t.Errorf("check output: %s", out)
		}
	})

	t.Run("summary", func(t *testing.T) {
		out, err := exec.Command(bin, "summary", fqPath).CombinedOutput()
		if err != nil {
			t.Fatalf("summary failed: %v\n%s", err, out)
		}
		for _, want := range []string{"records:    2", "bases:      8", "phred+33"} {
			if !strings.Contains(string(out), want) {
				t.Errorf("summary missing %q: %s", want, out)
			}
		}
	})

	t.Run("roundtrip", func(t *testing.T) {
		out, err := exec.Command(bin, "roundtrip", fqPath).CombinedOutput()
		if err != nil {
			t.Fatalf("roundtrip failed: %v\n%s", err, out)
		}
		if string(out) != src {
			t.Errorf("round-trip mismatch:\ngot:\n%s\nwant:\n%s", out, src)
		}
	})
}

// repoRoot finds the htsfile repo root by walking up from
// the current directory until it sees go.mod.
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
