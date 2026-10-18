package fasta

import (
	"os"
	"path/filepath"
	"sort"
	"testing"
)

func TestReference_BaseAt_PlainScanned(t *testing.T) {
	// Build a small FASTA with 2 contigs of known width.
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	const content = ">chr1 ACGT\nACGTACGT\nACGTACGT\n>chr2 XYZ\nGATTACA\nGATTACA\n"
	if err := os.WriteFile(path, []byte(content), 0o600); err != nil {
		t.Fatalf("write: %v", err)
	}
	ref, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer func() { _ = ref.Close() }()

	cases := []struct {
		chrom string
		pos   int
		want  byte
	}{
		{"chr1", 1, 'A'},
		{"chr1", 2, 'C'},
		{"chr1", 4, 'T'},
		{"chr1", 5, 'A'}, // wraps to next line, first base
		{"chr1", 8, 'T'},
		{"chr1", 16, 'T'}, // last base of chr1 (8 + 8 = 16)
		{"chr2", 1, 'G'},
		{"chr2", 4, 'T'},
		{"chr2", 8, 'G'}, // start of chr2 line 2 (chr2 line 1 is 7 chars)
	}
	for _, c := range cases {
		got, err := ref.BaseAt(c.chrom, c.pos)
		if err != nil {
			t.Errorf("%s:%d: %v", c.chrom, c.pos, err)
			continue
		}
		if got != c.want {
			t.Errorf("%s:%d = %q, want %q", c.chrom, c.pos, got, c.want)
		}
	}
}

func TestReference_BaseAt_WithFAI(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	const content = ">chr1 ACGT\nACGTACGT\nACGTACGT\n>chr2 XYZ\nGATTACA\nGATTACA\n"
	if err := os.WriteFile(path, []byte(content), 0o600); err != nil {
		t.Fatalf("write: %v", err)
	}
	// Build a sibling .fai manually. The .fai offset is
	// the byte position of the first base (past the
	// header line).
	//
	//   ">chr1 ACGT\n"     11 bytes (0..10)
	//   "ACGTACGT\n"         9 bytes (11..19)
	//   "ACGTACGT\n"         9 bytes (20..28)
	//   ">chr2 XYZ\n"      10 bytes (29..38)
	//
	// So chr1 starts at byte 11 and chr2 starts at byte
	// 39. chr2 line width is 7 + 1 = 8.
	fai := "chr1\t16\t11\t8\t9\nchr2\t14\t39\t7\t8\n"
	if err := os.WriteFile(path+".fai", []byte(fai), 0o600); err != nil {
		t.Fatalf("write .fai: %v", err)
	}
	ref, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer func() { _ = ref.Close() }()
	b, err := ref.BaseAt("chr1", 9)
	if err != nil {
		t.Fatalf("BaseAt: %v", err)
	}
	if b != 'A' {
		t.Errorf("chr1:9 = %q, want A", b)
	}
	if got := ref.Contigs(); len(got) != 2 {
		t.Errorf("Contigs: %v, want 2 entries", got)
	} else {
		sort.Strings(got)
		if got[0] != "chr1" || got[1] != "chr2" {
			t.Errorf("Contigs = %v", got)
		}
	}
}

func TestReference_UnknownContig(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(path, []byte(">chr1\nACGT\n"), 0o600); err != nil {
		t.Fatalf("write: %v", err)
	}
	ref, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer func() { _ = ref.Close() }()
	if _, err := ref.BaseAt("chrX", 1); err == nil {
		t.Errorf("expected error for unknown contig")
	}
}

func TestReference_PositionOutOfRange(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(path, []byte(">chr1\nACGT\n"), 0o600); err != nil {
		t.Fatalf("write: %v", err)
	}
	ref, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer func() { _ = ref.Close() }()
	if _, err := ref.BaseAt("chr1", 100); err == nil {
		t.Errorf("expected error for out-of-range position")
	}
}

func TestReference_VariableWidthHandled(t *testing.T) {
	// 10-char lines for chr1, 8-char lines for chr2 — same
	// call must work because we use the per-contig width.
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	const content = ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nGATTACAG\nATTACAGA\n"
	if err := os.WriteFile(path, []byte(content), 0o600); err != nil {
		t.Fatalf("write: %v", err)
	}
	ref, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer func() { _ = ref.Close() }()
	b, err := ref.BaseAt("chr1", 11)
	if err != nil {
		t.Fatalf("BaseAt chr1:11: %v", err)
	}
	if b != 'G' {
		t.Errorf("chr1:11 = %q, want G", b)
	}
	b, err = ref.BaseAt("chr2", 9)
	if err != nil {
		t.Fatalf("BaseAt chr2:9: %v", err)
	}
	if b != 'A' {
		t.Errorf("chr2:9 = %q, want A", b)
	}
}