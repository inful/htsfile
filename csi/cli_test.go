package csi

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

	// Build a synthetic CSI index in-memory.
	ref := &Reference{
		Bins:      make(map[uint32][]Chunk),
		Intervals: []VirtualOffset{0x1000, 0x2000, 0x3000},
	}
	ref.Bins[4681] = []Chunk{{Beg: 0x1000, End: 0x2000}}
	refs := []*Reference{ref}
	data := makeIndex(14, refs)
	idxPath := filepath.Join(dir, "smoke.csi")
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
		if !strings.Contains(string(out), "min_shift:") {
			t.Errorf("summary output missing 'min_shift:': %s", out)
		}
	})

	t.Run("roundtrip", func(t *testing.T) {
		// CSI now has a Writer (v2 follow-up resolved); round-trip
		// reads the synthetic index and re-emits it. The Writer
		// sorts chunks deterministically by Beg, so the output
		// re-parses into a logically-equivalent Index.
		out, err := exec.Command(bin, "roundtrip", idxPath).CombinedOutput()
		if err != nil {
			t.Fatalf("roundtrip failed: %v\n%s", err, out)
		}
		idx2, err := NewReader(bytes.NewReader(out)).Read()
		if err != nil {
			t.Fatalf("re-parse: %v", err)
		}
		if idx2.MinShift != 14 || idx2.Depth != 5 {
			t.Errorf("round-tripped min_shift/depth: got %d/%d want 14/5",
				idx2.MinShift, idx2.Depth)
		}
	})
}

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

var _ = bytes.NewReader
var _ = binary.LittleEndian
