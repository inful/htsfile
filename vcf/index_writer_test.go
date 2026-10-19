package vcf

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
)

// writeTestVCF builds a small BGZF-compressed VCF
// fixture with nRecords and returns the path.
func writeTestVCF(t *testing.T, dir string, nRecords int) string {
	t.Helper()
	vcfPath := filepath.Join(dir, "test.vcf.gz")

	content := make([]byte, 0, 200+nRecords*30)
	content = append(content, []byte("##fileformat=VCFv4.2\n")...)
	content = append(content, []byte("##contig=<ID=chr1,length=200000>\n")...)
	content = append(content, []byte("##contig=<ID=chr2,length=200000>\n")...)
	content = append(content, []byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")...)
	for i := 1; i <= nRecords; i++ {
		chrom := "chr1"
		if i > nRecords/2 {
			chrom = "chr2"
		}
		_, _ = fmt.Fprintf((*bytesAppender)(&content), "%s\t%d\t.\tA\tG\t.\t.\t.\n", chrom, i*100)
	}

	f, err := os.Create(vcfPath) //nolint:gosec // test fixture in t.TempDir
	if err != nil {
		t.Fatalf("create vcf: %v", err)
	}
	bw := bgzf.NewWriter(f)
	if _, err := bw.Write(content); err != nil {
		t.Fatalf("bw.Write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bw.Close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("f.Close: %v", err)
	}
	stat, err := os.Stat(vcfPath) //nolint:gosec // test fixture in t.TempDir
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("fixture %s: %d records, file size %d", vcfPath, nRecords, stat.Size())
	return vcfPath
}

// bytesAppender adapts *[]byte to io.Writer for
// fmt.Fprintf usage in tests. The Write method
// never fails because *[]byte grows in place.
type bytesAppender []byte

func (b *bytesAppender) Write(p []byte) (int, error) {
	*b = append(*b, p...)
	return len(p), nil
}

// TestBuildIndex_Tabix verifies that BuildIndex emits
// a tabix .tbi that round-trips through OpenTabix and
// yields the same record set as the file contains.
func TestBuildIndex_Tabix(t *testing.T) {
	dir := t.TempDir()
	vcfPath := writeTestVCF(t, dir, 200)

	if err := BuildIndex(vcfPath, IndexOpts{Format: FormatTabix}); err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}

	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	hdr, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if len(hdr.Contigs) != 2 {
		t.Errorf("Contigs: got %d, want 2", len(hdr.Contigs))
	}

	// Full contig query: every record returned.
	recs, err := r.Query("chr1", 1, 200000)
	if err != nil {
		t.Fatalf("Query chr1: %v", err)
	}
	if len(recs) != 100 {
		t.Errorf("chr1: got %d records, want 100", len(recs))
	}
	// First record should have POS=100 (i*100 for i=1).
	if recs[0].Pos != 100 {
		t.Errorf("chr1[0].Pos: got %d, want 100", recs[0].Pos)
	}

	recs, err = r.Query("chr2", 1, 200000)
	if err != nil {
		t.Fatalf("Query chr2: %v", err)
	}
	if len(recs) != 100 {
		t.Errorf("chr2: got %d records, want 100", len(recs))
	}
}

// TestBuildIndex_CSI verifies that BuildIndex emits
// a CSI .csi that round-trips through OpenCSI.
func TestBuildIndex_CSI(t *testing.T) {
	dir := t.TempDir()
	vcfPath := writeTestVCF(t, dir, 200)

	if err := BuildIndex(vcfPath, IndexOpts{
		Format:   FormatCSI,
		MinShift: 14,
		Depth:    5,
	}); err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}

	r, err := OpenCSI(vcfPath)
	if err != nil {
		t.Fatalf("OpenCSI: %v", err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	recs, err := r.Query("chr1", 1, 200000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 100 {
		t.Errorf("got %d records, want 100", len(recs))
	}
}

// TestBuildIndex_RejectsInvalidOpts verifies that
// nonsense min_shift and depth values are rejected.
func TestBuildIndex_RejectsInvalidOpts(t *testing.T) {
	dir := t.TempDir()
	vcfPath := writeTestVCF(t, dir, 10)

	if err := BuildIndex(vcfPath, IndexOpts{Format: FormatCSI, MinShift: -1, Depth: 5}); err == nil {
		t.Error("BuildIndex with min_shift=-1: expected error")
	}
	if err := BuildIndex(vcfPath, IndexOpts{Format: FormatCSI, MinShift: 14, Depth: 0}); err == nil {
		t.Error("BuildIndex with depth=0: expected error")
	}
}

// TestBuildIndex_SmallRegion checks that a build→query
// round-trip with a narrow region yields just the matching
// record, which exercises the chunk-boundary code path
// the chunk-End stop test covers.
func TestBuildIndex_SmallRegion(t *testing.T) {
	dir := t.TempDir()
	vcfPath := writeTestVCF(t, dir, 5000)
	if err := BuildIndex(vcfPath, IndexOpts{Format: FormatTabix}); err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	// The fixture has records at POS 100, 200, ..., 500000
	// for chr1 (records 1..2500). Query chr1:300-301 — only
	// POS=300 should match.
	recs, err := r.Query("chr1", 300, 301)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 1 {
		t.Fatalf("got %d records, want 1", len(recs))
	}
	if recs[0].Pos != 300 {
		t.Errorf("got POS=%d, want 300", recs[0].Pos)
	}
}

// TestBuildIndex_EOFBehaviour checks that the builder
// correctly closes any chunks still open at EOF. We use
// a single-record fixture so there's exactly one open
// chunk per bin.
func TestBuildIndex_EOFBehaviour(t *testing.T) {
	dir := t.TempDir()
	vcfPath := writeTestVCF(t, dir, 10)
	if err := BuildIndex(vcfPath, IndexOpts{Format: FormatTabix}); err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	recs, err := r.Query("chr1", 1, 200000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 5 {
		t.Errorf("got %d records, want 5", len(recs))
	}
}

// TestBuildIndex_BadPath checks error paths.
func TestBuildIndex_BadPath(t *testing.T) {
	if err := BuildIndex("/no/such/path.vcf.gz", IndexOpts{}); err == nil {
		t.Error("BuildIndex on missing path: expected error")
	}
}

// _ ensures io.EOF and errors.Is stay linked even if
// the error-handling code above is refactored.
var (
	_ = io.EOF
	_ = errors.Is
)