package vcf

import (
	"bytes"
	"strings"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
)

// plainVCF is a small, well-formed VCF used by the Open
// tests below.
const plainVCF = "##fileformat=VCFv4.2\n" +
	"##contig=<ID=chr1,length=1000>\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
	"chr1\t100\t.\tA\tC\t.\t.\t.\n" +
	"chr1\t200\t.\tG\tT\t.\t.\t.\n"

// bgzfBytes returns the BGZF-compressed bytes of the
// given VCF body.
func bgzfBytes(t *testing.T, body string) []byte {
	t.Helper()
	var buf bytes.Buffer
	w := bgzf.NewWriter(&buf)
	if _, err := w.Write([]byte(body)); err != nil {
		t.Fatalf("bgzf write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("bgzf close: %v", err)
	}
	return buf.Bytes()
}

func TestOpen_PlainVCFFromBytes(t *testing.T) {
	dec := NewDecoder(strings.NewReader(plainVCF))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr.FileFormat.String() != "VCFv4.2" {
		t.Errorf("FileFormat = %q", hdr.FileFormat.String())
	}
	r, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if r.Pos != 100 {
		t.Errorf("first record POS = %d, want 100", r.Pos)
	}
}

func TestOpenReader_DetectsBGZF(t *testing.T) {
	compressed := bgzfBytes(t, plainVCF)
	dec := openReader(bytes.NewReader(compressed))
	// Verify the BGZF wrapper passes through to the
	// underlying VCF text.
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr.FileFormat.String() != "VCFv4.2" {
		t.Errorf("FileFormat = %q", hdr.FileFormat.String())
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.Pos != 100 {
		t.Errorf("first record POS = %d, want 100", rec.Pos)
	}
}

func TestOpenReader_LeavesPlainVCFAlone(t *testing.T) {
	dec := openReader(strings.NewReader(plainVCF))
	// Plain VCF must not be wrapped in BGZF; reading
	// the header should succeed immediately and the
	// first two bytes should be "##" (the VCF meta
	// prefix), not the gzip magic.
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr.FileFormat.String() != "VCFv4.2" {
		t.Errorf("FileFormat = %q", hdr.FileFormat.String())
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.Pos != 100 {
		t.Errorf("first record POS = %d, want 100", rec.Pos)
	}
}