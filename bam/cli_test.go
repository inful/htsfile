package bam

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
// runs its `check`, `summary`, and `roundtrip` subcommands
// against a synthetic BAM written to disk. The roundtrip test
// verifies the BAM v1 writer reads back through the existing
// decoder.
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

	// Run the CLI's roundtrip subcommand. BAM v1 now has a
	// writer (added in the 'close symmetry' Phase 1), so
	// roundtrip reads and re-emits the file. The output is
	// byte-stable because records stream through verbatim.
	t.Run("roundtrip", func(t *testing.T) {
		out, err := exec.Command(bin, "roundtrip", bamPath).CombinedOutput()
		if err != nil {
			t.Fatalf("roundtrip failed: %v\n%s", err, out)
		}
		// Re-parse the round-tripped output and verify it's
		// structurally equivalent (same record count, same
		// names).
		dec := NewDecoder(bytes.NewReader(out))
		hdr2, err := dec.Header()
		if err != nil {
			t.Fatalf("re-parse header: %v", err)
		}
		if len(hdr2.SQs) != 1 || hdr2.SQs[0].Name != "chr1" {
			t.Errorf("round-tripped header: %+v", hdr2)
		}
		var n int
		for {
			rec, err := dec.Record()
			if errors.Is(err, io.EOF) {
				break
			}
			if err != nil {
				t.Fatalf("re-parse record %d: %v", n, err)
			}
			if rec.QName != "smoke1" {
				t.Errorf("round-tripped QNAME: got %q want smoke1", rec.QName)
			}
			n++
		}
		if n != 1 {
			t.Errorf("round-tripped records: got %d want 1", n)
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
