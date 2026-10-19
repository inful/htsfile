package vcf

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/inful/htsfile/csi"
	"github.com/inful/htsfile/internal/bgzf"
)

// makeIndexedCSI writes a BGZF-compressed VCF file and a
// CSI index that covers the entire file in one chunk.
func makeIndexedCSI(t *testing.T) (vcfPath string) {
	t.Helper()
	dir := t.TempDir()
	vcfPath = filepath.Join(dir, "sample.vcf.gz")
	csiPath := vcfPath + ".csi"
	_ = csiPath

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

	// CSI uses MinShift (window = 2^MinShift bp) and
	// Depth (number of bin levels). Set up for the default
	// 5-level scheme and put the chunk in bin 0.
	idx := &csi.Index{
		MinShift: 14, // 16 Kbp windows
		Depth:    5,
		References: []*csi.Reference{
			{
				Bins: map[uint32][]csi.Chunk{
					0: {{Beg: 0, End: 0x10000}},
				},
				Intervals: []csi.VirtualOffset{0, 0x10000},
			},
		},
	}

	csiFile, err := os.Create(csiPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("create csi: %v", err)
	}
	if err := csi.NewWriter(csiFile).Write(idx); err != nil {
		_ = csiFile.Close()
		t.Fatalf("write csi: %v", err)
	}
	if err := csiFile.Close(); err != nil {
		t.Fatalf("close csi: %v", err)
	}

	return vcfPath
}

func TestCSIReader_OpenAndHeader(t *testing.T) {
	vcfPath := makeIndexedCSI(t)
	r, err := OpenCSI(vcfPath)
	if err != nil {
		t.Fatalf("OpenCSI: %v", err)
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

func TestCSIReader_QueryRegion(t *testing.T) {
	vcfPath := makeIndexedCSI(t)
	r, err := OpenCSI(vcfPath)
	if err != nil {
		t.Fatalf("OpenCSI: %v", err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
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
