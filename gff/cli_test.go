package gff

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// TestCLI_CheckAndSummary builds the cmd/htsfile binary and
// runs its `check`, `roundtrip`, and `summary` subcommands
// against a synthetic GFF3 file. The smoke test confirms the
// extension dispatch and content-sniff paths work end-to-end.
func TestCLI_CheckAndSummary(t *testing.T) {
	root := repoRoot(t)
	dir := t.TempDir()
	bin := filepath.Join(dir, "htsfile")
	cmd := exec.Command("go", "build", "-o", bin, "./cmd/htsfile")
	cmd.Dir = root
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build: %v\n%s", err, out)
	}

	src := "##gff-version 3\n" +
		"##sequence-region ctgA 1 1497228\n" +
		"ctgA\tEMBL\tgene\t1000\t9000\t.\t+\t.\tID=gene1;Name=ABC1\n" +
		"ctgA\tEMBL\texon\t1000\t2000\t.\t+\t.\tParent=gene1\n"
	gffPath := filepath.Join(dir, "smoke.gff")
	if err := os.WriteFile(gffPath, []byte(src), 0o644); err != nil {
		t.Fatal(err)
	}

	t.Run("check", func(t *testing.T) {
		out, err := exec.Command(bin, "check", gffPath).CombinedOutput()
		if err != nil {
			t.Fatalf("check failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "clean") {
			t.Errorf("check output missing 'clean': %s", out)
		}
	})

	t.Run("summary", func(t *testing.T) {
		out, err := exec.Command(bin, "summary", gffPath).CombinedOutput()
		if err != nil {
			t.Fatalf("summary failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "records:") {
			t.Errorf("summary output missing 'records:': %s", out)
		}
	})

	t.Run("roundtrip", func(t *testing.T) {
		out, err := exec.Command(bin, "roundtrip", gffPath).CombinedOutput()
		if err != nil {
			t.Fatalf("roundtrip failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "##gff-version 3") {
			t.Errorf("roundtrip output missing version line:\n%s", out)
		}
		if !strings.Contains(string(out), "ID=gene1") {
			t.Errorf("roundtrip output missing ID=gene1:\n%s", out)
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
