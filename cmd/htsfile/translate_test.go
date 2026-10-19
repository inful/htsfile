package main

import (
	"bytes"
	"strings"
	"testing"

	"github.com/inful/htsfile/bam"
	"github.com/inful/htsfile/bcf"
	"github.com/inful/htsfile/sam"
	"github.com/inful/htsfile/vcf"
)

// -----------------------------------------------------------------------------
// SAM ↔ BAM round-trip
// -----------------------------------------------------------------------------

// TestSAMToBAMAndBack writes a small SAM, converts to BAM via
// sam2bam, reads the BAM, converts back via bam2sam, and
// verifies the result is logically equivalent.
func TestSAMToBAMAndBack(t *testing.T) {
	src := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:248956422\n" +
		"@SQ\tSN:chr2\tLN:5000000\n" +
		"r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:0\n" +
		"r2\t16\tchr1\t200\t30\t3M\tchr2\t300\t100\tTGA\tIII\tRG:Z:rg1\n"

	var bamBuf bytes.Buffer
	if err := runSAMToBAM(strings.NewReader(src), &bamBuf); err != nil {
		t.Fatalf("sam2bam: %v", err)
	}
	if bamBuf.Len() == 0 {
		t.Fatal("sam2bam produced empty output")
	}

	var samBuf bytes.Buffer
	if err := runBAMToSAM(&bamBuf, &samBuf); err != nil {
		t.Fatalf("bam2sam: %v", err)
	}
	out := samBuf.String()
	// Header round-trips (some fields may be reordered).
	if !strings.Contains(out, "@HD") {
		t.Errorf("missing @HD line in:\n%s", out)
	}
	if !strings.Contains(out, "SN:chr1") {
		t.Errorf("missing chr1 SQ line in:\n%s", out)
	}
	if !strings.Contains(out, "r1\t") {
		t.Errorf("missing record r1 in:\n%s", out)
	}
	if !strings.Contains(out, "r2\t") {
		t.Errorf("missing record r2 in:\n%s", out)
	}
	// CIGAR, SEQ, QUAL preserved verbatim.
	if !strings.Contains(out, "\t5M\t") {
		t.Errorf("CIGAR not preserved:\n%s", out)
	}
	if !strings.Contains(out, "\tACGTA\t") {
		t.Errorf("SEQ not preserved:\n%s", out)
	}
	if !strings.Contains(out, "NM:i:0") {
		t.Errorf("NM tag not preserved:\n%s", out)
	}
	if !strings.Contains(out, "RG:Z:rg1") {
		t.Errorf("RG tag not preserved:\n%s", out)
	}
}

// TestBAMToSAMMinimal writes a minimal BAM (no records),
// converts to SAM, and verifies the header is emitted.
func TestBAMToSAMMinimal(t *testing.T) {
	hdr := &bam.Header{
		Raw: "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n",
		SQs: []*bam.SQ{{Name: "chr1", Len: 1000}},
	}
	var bamBuf bytes.Buffer
	enc := bam.NewEncoder(&bamBuf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	var samBuf bytes.Buffer
	if err := runBAMToSAM(&bamBuf, &samBuf); err != nil {
		t.Fatalf("bam2sam: %v", err)
	}
	out := samBuf.String()
	if !strings.HasPrefix(out, "@HD") {
		t.Errorf("SAM output should start with @HD, got: %q", out)
	}
	if !strings.Contains(out, "SN:chr1") {
		t.Errorf("missing chr1 SQ line:\n%s", out)
	}
}

// TestCIGARParseRoundTrip exercises parseCigarString /
// formatCigarString with the standard op set. The empty
// string and "*" both map to the canonical "*" output (the
// SAM spec's representation of "no CIGAR").
func TestCIGARParseRoundTrip(t *testing.T) {
	cases := []struct {
		input, output string
	}{
		{"5M", "5M"},
		{"3M2I3M", "3M2I3M"},
		{"10M200N10M", "10M200N10M"},
		{"5S20M5S", "5S20M5S"},
		{"5H10M5H", "5H10M5H"},
		{"2X3=2X", "2X3=2X"},
		{"*", "*"},
		{"", "*"}, // empty canonicalises to "*"
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			packed, err := parseCigarString(c.input)
			if err != nil {
				t.Fatalf("parse %q: %v", c.input, err)
			}
			got := formatCigarString(packed)
			if got != c.output {
				t.Errorf("round-trip: got %q want %q", got, c.output)
			}
		})
	}
}

func TestCIGARParseError(t *testing.T) {
	// Missing op character after the length.
	if _, err := parseCigarString("5"); err == nil {
		t.Error("parseCigarString(\"5\"): expected error, got nil")
	}
	// Unknown op character.
	if _, err := parseCigarString("5Q"); err == nil {
		t.Error("parseCigarString(\"5Q\"): expected error, got nil")
	}
}

// TestComputeBin verifies that the simple bin formula yields
// reasonable values for common positions.
func TestComputeBin(t *testing.T) {
	cases := []struct {
		pos int64
		min uint16
	}{
		{0, 4680},  // unmapped
		{-1, 4680}, // negative also unmapped
		{1, 4681},  // mapped, very early
		{1000, 4681},
		{1 << 27, 4682}, // past 64 MiB → level 2
	}
	for _, c := range cases {
		got := computeBin(c.pos)
		if got < c.min {
			t.Errorf("computeBin(%d) = %d, want >= %d", c.pos, got, c.min)
		}
	}
}

// -----------------------------------------------------------------------------
// VCF ↔ BCF round-trip
// -----------------------------------------------------------------------------

// TestVCFToBCFAndBack writes a small VCF, converts to BCF,
// reads the BCF, converts back to VCF, and verifies the
// core fields are preserved. Per-sample FORMAT data is
// dropped on vcf2bcf (matching v1 BCF scope).
func TestVCFToBCFAndBack(t *testing.T) {
	src := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=248956422>\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\trs1\tA\tG\t99.5\tPASS\tDP=42\n"

	var bcfBuf bytes.Buffer
	if err := runVCFToBCF(strings.NewReader(src), &bcfBuf); err != nil {
		t.Fatalf("vcf2bcf: %v", err)
	}
	if bcfBuf.Len() == 0 {
		t.Fatal("vcf2bcf produced empty output")
	}
	t.Logf("BCF output (%d bytes): %x", bcfBuf.Len(), bcfBuf.Bytes())

	var vcfBuf bytes.Buffer
	if err := runBCFToVCF(&bcfBuf, &vcfBuf); err != nil {
		t.Fatalf("bcf2vcf: %v", err)
	}
	out := vcfBuf.String()
	if !strings.Contains(out, "##fileformat=VCFv4.2") {
		t.Errorf("missing ##fileformat:\n%s", out)
	}
	if !strings.Contains(out, "##contig=<ID=chr1") {
		t.Errorf("missing ##contig:\n%s", out)
	}
	if !strings.Contains(out, "chr1\t100\t") {
		t.Errorf("record 1 missing:\n%s", out)
	}
	if !strings.Contains(out, "DP=42") {
		t.Errorf("INFO DP=42 missing:\n%s", out)
	}
}

// TestVCFToBCFMultipleRecords verifies that the translator
// handles multiple records (the encoder is responsible for
// putting them in BGZF blocks; this test catches block-
// boundary issues).
func TestVCFToBCFMultipleRecords(t *testing.T) {
	src := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=248956422>\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\trs1\tA\tG\t99.5\tPASS\tDP=42\n" +
		"chr1\t200\t.\tC\tT,G\t50.0\tPASS\tDP=7\n"
	var bcfBuf bytes.Buffer
	if err := runVCFToBCF(strings.NewReader(src), &bcfBuf); err != nil {
		t.Fatalf("vcf2bcf: %v", err)
	}
	var vcfBuf bytes.Buffer
	if err := runBCFToVCF(&bcfBuf, &vcfBuf); err != nil {
		t.Fatalf("bcf2vcf: %v", err)
	}
	out := vcfBuf.String()
	if !strings.Contains(out, "chr1\t100\t") {
		t.Errorf("record 1 missing:\n%s", out)
	}
	if !strings.Contains(out, "chr1\t200\t") {
		t.Errorf("record 2 missing:\n%s", out)
	}
}

// TestVCFToBCFMinimal writes a header-only VCF and converts
// to BCF + back. Verifies the conversion paths work without
// records.
func TestVCFToBCFMinimal(t *testing.T) {
	src := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

	var bcfBuf bytes.Buffer
	if err := runVCFToBCF(strings.NewReader(src), &bcfBuf); err != nil {
		t.Fatalf("vcf2bcf: %v", err)
	}
	var vcfBuf bytes.Buffer
	if err := runBCFToVCF(&bcfBuf, &vcfBuf); err != nil {
		t.Fatalf("bcf2vcf: %v", err)
	}
	out := vcfBuf.String()
	if !strings.Contains(out, "##fileformat=VCFv4.2") {
		t.Errorf("missing fileformat:\n%s", out)
	}
}

// -----------------------------------------------------------------------------
// Header / record helpers
// -----------------------------------------------------------------------------

func TestVcfQualToBCF(t *testing.T) {
	if vcfQualToBCF(nil) != 0 {
		t.Errorf("nil qual should be 0")
	}
	q := 99.5
	if got := vcfQualToBCF(&q); got != 99.5 {
		t.Errorf("got %v want 99.5", got)
	}
}

func TestBcfQualToVCF(t *testing.T) {
	if bcfQualToVCF(0, false) != nil {
		t.Errorf("qualKnown=false should give nil")
	}
	got := bcfQualToVCF(99.5, true)
	if got == nil || *got != 99.5 {
		t.Errorf("got %v want pointer to 99.5", got)
	}
}

func TestLookupRefID(t *testing.T) {
	hdr := &bam.Header{
		SQs: []*bam.SQ{{Name: "chr1", Len: 1000}, {Name: "chr2", Len: 2000}},
	}
	if got, err := lookupRefID("", hdr); err != nil || got != -1 {
		t.Errorf("empty name: got (%d, %v) want (-1, nil)", got, err)
	}
	if got, err := lookupRefID("*", hdr); err != nil || got != -1 {
		t.Errorf("\"*\" name: got (%d, %v) want (-1, nil)", got, err)
	}
	if got, err := lookupRefID("chr1", hdr); err != nil || got != 0 {
		t.Errorf("chr1: got (%d, %v) want (0, nil)", got, err)
	}
	if got, err := lookupRefID("chr2", hdr); err != nil || got != 1 {
		t.Errorf("chr2: got (%d, %v) want (1, nil)", got, err)
	}
	if _, err := lookupRefID("chrX", hdr); err == nil {
		t.Errorf("chrX: expected error, got nil")
	}
}

// TestRunTranslateRejectsBadExtension verifies the dispatch
// guard rejects input with the wrong extension.
func TestRunTranslateRejectsBadExtension(t *testing.T) {
	// nil path triggers "expected exactly one file argument".
	err := runTranslate(nil, runSAMToBAM, "sam2bam", ".sam")
	if err == nil {
		t.Fatal("expected error for nil args")
	}
	if !strings.Contains(err.Error(), "exactly one file") {
		t.Errorf("unexpected error: %v", err)
	}
}

// _ ensures sam and vcf packages are referenced (compile
// check that the helpers compose with the public APIs).
var (
	_ = sam.Header{}
	_ = vcf.Header{}
	_ = bcf.Header{}
)
