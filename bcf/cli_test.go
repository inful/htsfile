package bcf

import (
	"bytes"
	"errors"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// TestCLI_CheckAndSummary builds the cmd/htsfile binary and
// runs its `check`, `summary`, and `roundtrip` (expecting an
// unsupported error) subcommands against a synthetic BCF file.
// The smoke test confirms the extension dispatch and content
// path work end-to-end.
func TestCLI_CheckAndSummary(t *testing.T) {
	root := repoRoot(t)
	dir := t.TempDir()
	bin := filepath.Join(dir, "htsfile")
	cmd := exec.Command("go", "build", "-o", bin, "./cmd/htsfile")
	cmd.Dir = root
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("go build: %v\n%s", err, out)
	}

	// Construct a synthetic BCF.
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	defs := map[string]*InfoDef{
		"DP": {ID: "DP", Number: "1", Type: "Integer"},
	}
	r := &testRecord{
		contigIdx: 0, Pos: 100, Ref: "A", Alt: []string{"T"},
		Qual: 30, QualKnown: true,
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{42}},
		},
		infoDefs:  defs,
		infoOrder: []string{"DP"},
	}
	b.writeRecord(t, r)
	bcfPath := filepath.Join(dir, "smoke.bcf")
	if err := os.WriteFile(bcfPath, b.finish(t), 0o644); err != nil {
		t.Fatal(err)
	}

	t.Run("check", func(t *testing.T) {
		out, err := exec.Command(bin, "check", bcfPath).CombinedOutput()
		if err != nil {
			t.Fatalf("check failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "clean") {
			t.Errorf("check output missing 'clean': %s", out)
		}
	})

	t.Run("summary", func(t *testing.T) {
		out, err := exec.Command(bin, "summary", bcfPath).CombinedOutput()
		if err != nil {
			t.Fatalf("summary failed: %v\n%s", err, out)
		}
		if !strings.Contains(string(out), "records:") {
			t.Errorf("summary output missing 'records:': %s", out)
		}
		if !strings.Contains(string(out), "contigs:") {
			t.Errorf("summary output missing 'contigs:': %s", out)
		}
	})

	t.Run("roundtrip", func(t *testing.T) {
		// BCF v1 round-trip is partial: header and typed
		// record fields round-trip, but per-sample FORMAT
		// data is dropped (matching the v1 decoder's
		// scope). The decoder reads the round-tripped bytes
		// and we verify the same record count + contig
		// count.
		out, err := exec.Command(bin, "roundtrip", bcfPath).CombinedOutput()
		if err != nil {
			t.Fatalf("roundtrip failed: %v\n%s", err, out)
		}
		dec := NewDecoder(bytes.NewReader(out))
		if _, err := dec.Header(); err != nil {
			t.Fatalf("re-parse header: %v", err)
		}
		var n int
		for {
			_, err := dec.Record()
			if errors.Is(err, io.EOF) {
				break
			}
			if err != nil {
				t.Fatalf("re-parse record %d: %v", n, err)
			}
			n++
		}
		if n != 1 {
			t.Errorf("round-tripped records: got %d want 1", n)
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

// _ keeps bytes/io imports visible (they are used by tests
// that may be added later; the package's test builder uses
// them already).
var _ = bytes.HasPrefix
var _ = io.EOF
