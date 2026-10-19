package vcf

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
	"github.com/inful/htsfile/tabix"
)

// makeIndexedVCF writes a BGZF-compressed VCF file and a
// tabix index that covers the entire file in one chunk.
// Returns the VCF path.
func makeIndexedVCF(t *testing.T) (vcfPath string) {
	t.Helper()
	dir := t.TempDir()
	vcfPath = filepath.Join(dir, "sample.vcf.gz")
	tbiPath := vcfPath + ".tbi"
	_ = tbiPath

	// Step 1: BGZF-compress the VCF content.
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

	// Step 2: compute the BGZF virtual offset for each
	// line. Tabix chunk Beg/End are virtual offsets.
	// For a one-chunk-per-record index, we use a single
	// chunk that covers the whole file (offset 0 to EOF).
	//
	// The EOF marker of a BGZF file is at virtual offset
	// (block_count << 16) | 0. For a file with one block,
	// that's 0x10000 | 0 = 65536.
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

	// Step 3: write the index file.
	tbiFile, err := os.Create(tbiPath) //nolint:gosec // test file in t.TempDir
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

func TestTabixReader_OpenAndHeader(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	hdr, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if len(hdr.Contigs) != 1 || hdr.Contigs[0].ID != "chr1" {
		t.Errorf("Contigs: %+v, want one chr1", hdr.Contigs)
	}
}

func TestTabixReader_QueryRegion(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	// Query region that overlaps pos 100 and 300 but not
	// 500. With one chunk covering the whole file, we'll
	// fetch all three records and filter by POS.
	recs, err := r.Query("chr1", 200, 400)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 1 {
		t.Fatalf("Query: got %d records, want 1", len(recs))
	}
	if recs[0].Pos != 300 {
		t.Errorf("Query: got POS %d, want 300", recs[0].Pos)
	}
}

func TestTabixReader_QueryWholeContig(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	recs, err := r.Query("chr1", 1, 10000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 3 {
		t.Errorf("Query: got %d records, want 3", len(recs))
	}
}

func TestTabixReader_QueryUnknownContig(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	recs, err := r.Query("chrX", 1, 100)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 0 {
		t.Errorf("Query: got %d records, want 0", len(recs))
	}
}

func TestTabixReader_QueryInvalidRegion(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Query("chr1", 0, 100); err == nil {
		t.Error("Query(0, 100): expected error for 0-based start")
	}
	if _, err := r.Query("chr1", 100, 100); err == nil {
		t.Error("Query(100, 100): expected error for end <= start")
	}
}

func TestTabixReader_MissingIndex(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "no-index.vcf")
	if err := os.WriteFile(path, []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"), 0o600); err != nil {
		t.Fatalf("write: %v", err)
	}
	if _, err := OpenTabix(path); err == nil {
		t.Error("OpenTabix: expected error for missing .tbi")
	}
}

func TestTabixReader_ViaBGZFStream(t *testing.T) {
	// Confirm the BGZF stream round-trips through the
	// internal bgzf.Reader. The first two bytes must be
	// the gzip magic (1F 8B).
	vcfPath := makeIndexedVCF(t)
	f, err := os.Open(vcfPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer func() { _ = f.Close() }()
	hdr := make([]byte, 2)
	if _, err := f.Read(hdr); err != nil {
		t.Fatalf("read: %v", err)
	}
	if hdr[0] != 0x1f || hdr[1] != 0x8b {
		t.Fatalf("not BGZF: got %x %x", hdr[0], hdr[1])
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		t.Fatalf("seek: %v", err)
	}
	br := bgzf.NewReader(f)
	var buf bytes.Buffer
	if _, err := buf.ReadFrom(br); err != nil {
		t.Fatalf("bgzf read: %v", err)
	}
	if !bytes.Contains(buf.Bytes(), []byte("chr1\t100")) {
		t.Errorf("BGZF content missing chr1:100: %q", buf.String())
	}
}