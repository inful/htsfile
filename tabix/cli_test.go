package tabix

import (
	"bytes"
	"encoding/binary"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func TestCLI_CheckAndSummary(t *testing.T) {
	root := repoRoot(t)
	dir := t.TempDir()
	bin := filepath.Join(dir, "htsfile")
	cmd := exec.Command("go", "build", "-o", bin, "./cmd/htsfile")
	cmd.Dir = root
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build: %v\n%s", err, out)
	}

	// Build a synthetic Tabix index in-memory.
	ref := &Reference{
		Bins:      make(map[uint32][]Chunk),
		Intervals: []VirtualOffset{0x1000, 0x2000, 0x3000},
		NoCoor:    0,
	}
	ref.Bins[4681] = []Chunk{{Beg: 0x1000, End: 0x2000}}
	refs := []*Reference{ref}
	data := makeIndex(refs, true)
	idxPath := filepath.Join(dir, "smoke.tbi")
	if err := os.WriteFile(idxPath, data, 0o644); err != nil {
		t.Fatal(err)
	}

	t.Run("check", func(t *testing.T) {
		out, err := exec.Command(bin, "check", idxPath).CombinedOutput()
		if err != nil {
			t.Fatalf("check failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "clean") {
			t.Errorf("check output missing 'clean': %s", out)
		}
	})

	t.Run("summary", func(t *testing.T) {
		out, err := exec.Command(bin, "summary", idxPath).CombinedOutput()
		if err != nil {
			t.Fatalf("summary failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "references:") {
			t.Errorf("summary output missing 'references:': %s", out)
		}
	})

	t.Run("roundtrip_unsupported", func(t *testing.T) {
		out, err := exec.Command(bin, "roundtrip", idxPath).CombinedOutput()
		if err == nil {
			t.Fatalf("roundtrip should have failed for Tabix, but succeeded:\n%s", out)
		}
		combined := string(out)
		if !strings.Contains(combined, "read-only") {
			t.Errorf("roundtrip error doesn't mention read-only:\n%s", combined)
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

// _ keeps bytes / encoding/binary visible (used by makeIndex
// inherited from tabix_test.go).
var _ = bytes.NewReader
var _ = binary.LittleEndian
