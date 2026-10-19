package main

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
	"github.com/inful/htsfile/tabix"
)

// makeVCFQueryFixture writes a BGZF-compressed VCF and
// tabix index for the vcf query CLI tests.
func makeVCFQueryFixture(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "sample.vcf.gz")

	const content = "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=10000>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t.\t.\t.\n" +
		"chr1\t300\t.\tC\tT\t.\t.\t.\n" +
		"chr1\t500\t.\tG\tA\t.\t.\t.\n"

	f, err := os.Create(vcfPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("create vcf: %v", err)
	}
	bw := bgzf.NewWriter(f)
	if _, err := bw.Write([]byte(content)); err != nil {
		t.Fatalf("bw.Write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bw.Close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("f.Close: %v", err)
	}

	// Tabix index with one chunk covering the whole file.
	idx := &tabix.Index{
		References: []*tabix.Reference{
			{
				Bins: map[uint32][]tabix.Chunk{
					0: {{Beg: 0, End: 0x10000}},
				},
				Intervals: []tabix.VirtualOffset{0, 0x10000},
			},
		},
	}
	tbiFile, err := os.Create(vcfPath + ".tbi") //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("create tbi: %v", err)
	}
	if err := tabix.NewWriter(tbiFile).Write(idx); err != nil {
		_ = tbiFile.Close()
		t.Fatalf("write tbi: %v", err)
	}
	if err := tbiFile.Close(); err != nil {
		t.Fatalf("close tbi: %v", err)
	}
	return vcfPath
}

func TestRunVCFQuery_RegionMatches(t *testing.T) {
	vcfPath := makeVCFQueryFixture(t)

	// Redirect stdout to a buffer.
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("pipe: %v", err)
	}
	oldStdout := os.Stdout
	os.Stdout = w
	defer func() { os.Stdout = oldStdout }()

	var out bytes.Buffer
	done := make(chan struct{})
	go func() {
		_, _ = out.ReadFrom(r)
		close(done)
	}()

	if err := runVCFQuery([]string{vcfPath, "chr1:200-400"}); err != nil {
		t.Fatalf("runVCFQuery: %v", err)
	}
	_ = w.Close()
	<-done

	got := out.String()
	if !strings.Contains(got, "chr1\t300\t.\tC\tT") {
		t.Errorf("output missing expected record:\n%s", got)
	}
	if strings.Contains(got, "chr1\t100\t.") {
		t.Errorf("output should not include chr1:100 (out of range):\n%s", got)
	}
	if strings.Contains(got, "chr1\t500\t.") {
		t.Errorf("output should not include chr1:500 (out of range):\n%s", got)
	}
}

func TestRunVCFQuery_WholeContig(t *testing.T) {
	vcfPath := makeVCFQueryFixture(t)

	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("pipe: %v", err)
	}
	oldStdout := os.Stdout
	os.Stdout = w
	defer func() { os.Stdout = oldStdout }()

	var out bytes.Buffer
	done := make(chan struct{})
	go func() {
		_, _ = out.ReadFrom(r)
		close(done)
	}()

	if err := runVCFQuery([]string{vcfPath, "chr1"}); err != nil {
		t.Fatalf("runVCFQuery: %v", err)
	}
	_ = w.Close()
	<-done

	got := out.String()
	if got == "" {
		t.Errorf("empty output for whole-contig query")
	}
	// Whole-contig query returns all 3 records.
	for _, want := range []string{"chr1\t100", "chr1\t300", "chr1\t500"} {
		if !strings.Contains(got, want) {
			t.Errorf("output missing %q:\n%s", want, got)
		}
	}
}

func TestRunVCFQuery_InvalidRegion(t *testing.T) {
	vcfPath := makeVCFQueryFixture(t)
	if err := runVCFQuery([]string{vcfPath, "chr1:bad-region"}); err == nil {
		t.Error("expected error for malformed region")
	}
	if err := runVCFQuery([]string{vcfPath, "chr1:200-100"}); err == nil {
		t.Error("expected error for end < start")
	}
}

func TestRunVCFQuery_NoIndex(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "no-index.vcf.gz")
	if err := os.WriteFile(path, []byte("##fileformat=VCFv4.2\n"), 0o600); err != nil {
		t.Fatalf("write: %v", err)
	}
	if err := runVCFQuery([]string{path, "chr1:1-100"}); err == nil {
		t.Error("expected error when no index is present")
	}
}

func TestParseVCFRegion(t *testing.T) {
	cases := []struct {
		in        string
		wantCh    string
		wantStart int64
		wantEnd   int64
		wantErr   bool
	}{
		{"chr1:100-200", "chr1", 100, 200, false},
		{"chr1:1-1000000", "chr1", 1, 1000000, false},
		{"chr1", "chr1", 0, 0, false}, // whole-contig
		{"chr1:bad", "", 0, 0, true},
		{"chr1:200-100", "", 0, 0, true},
		{"chr1:notnum-200", "", 0, 0, true},
	}
	for _, c := range cases {
		t.Run(c.in, func(t *testing.T) {
			ch, st, en, err := parseVCFRegion(c.in)
			if (err != nil) != c.wantErr {
				t.Fatalf("err = %v, wantErr = %v", err, c.wantErr)
			}
			if err != nil {
				return
			}
			if ch != c.wantCh || st != c.wantStart || en != c.wantEnd {
				t.Errorf("got (%q, %d, %d), want (%q, %d, %d)",
					ch, st, en, c.wantCh, c.wantStart, c.wantEnd)
			}
		})
	}
}
