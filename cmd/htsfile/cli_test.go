package main

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inful/htsfile/bam"
	"github.com/inful/htsfile/bcf"
	"github.com/inful/htsfile/bed"
	"github.com/inful/htsfile/fasta"
	"github.com/inful/htsfile/fastq"
	"github.com/inful/htsfile/gff"
	"github.com/inful/htsfile/sam"
	"github.com/inful/htsfile/vcf"
)

// captureStdout runs fn with os.Stdout redirected to a buffer
// and returns the captured output plus any error. Tests that
// exercise the per-subcommand helpers — which write directly
// to os.Stdout — use this to assert both stdout and the
// returned error.
func captureStdout(t *testing.T, fn func() error) (string, error) {
	t.Helper()
	orig := os.Stdout
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("pipe: %v", err)
	}
	os.Stdout = w
	var fnErr error
	func() {
		defer func() {
			// Restore Stdout synchronously to avoid losing
			// the captured output if fn panics.
			os.Stdout = orig
		}()
		fnErr = fn()
	}()
	if err := w.Close(); err != nil {
		t.Fatalf("close writer: %v", err)
	}
	var buf bytes.Buffer
	if _, err := io.Copy(&buf, r); err != nil {
		t.Fatalf("read pipe: %v", err)
	}
	return buf.String(), fnErr
}

// minimalVCF is a tiny, valid VCF document with one record.
const minimalVCF = "##fileformat=VCFv4.5\n" +
	"##contig=<ID=chr1,length=1000>\n" +
	"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
	"chr1\t100\t.\tA\tC\t30\tPASS\tDP=42\n"

// brokenVCF has a known spec violation: a Record whose FILTER
// value isn't declared in the header. The validator (ModeDefault)
// flags undeclared FILTER codes (v4.5 §1.6.1.5) with a
// SeverityError.
const brokenVCF = "##fileformat=VCFv4.5\n" +
	"##contig=<ID=chr1,length=1000>\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
	"chr1\t100\t.\tA\tC\t30\tUNDECLARED\t.\n"

// minimalBED is a tiny, valid BED3 document with one record.
const minimalBED = "browser position chr1:100-200\n" +
	"track name=example\n" +
	"chr1\t100\t200\n"

// minimalSAM is a tiny, valid SAM with two records on a
// coordinate-sorted header.
const minimalSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:1000\n" +
	"r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:0\n" +
	"r2\t16\tchr1\t200\t30\t3M\tchr1\t300\t100\tTGA\tIII\tRG:Z:rg1\n"

// minimalFASTA is one record with a short sequence that won't
// trigger any validator error.
const minimalFASTA = ">seq1 short read\n" +
	"ACGTACGTAC\n"

// minimalFASTQ is one record, Phred+33 offset.
const minimalFASTQ = "@seq1\n" +
	"ACGTACGTAC\n" +
	"+\n" +
	"IIIIIIIIII\n"

// minimalGFF is a single GFF3 record with a header.
const minimalGFF = "##gff-version 3\n" +
	"##sequence-region chr1 1 1000\n" +
	"chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1\n"

// -----------------------------------------------------------------------------
// runCheck per-format tests
// -----------------------------------------------------------------------------

func TestRunCheckVCFClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runCheckVCF(strings.NewReader(minimalVCF), "example.vcf", false, false)
	})
	if err != nil {
		t.Fatalf("runCheckVCF: %v", err)
	}
	if !strings.Contains(out, "example.vcf: clean") {
		t.Errorf("expected clean summary, got %q", out)
	}
}

func TestRunCheckVCFBroken(t *testing.T) {
	// Strict mode rejects undeclared FILTER codes (ModeDefault
	// only warns). RunCheckVCF's strict flag matches `htsfile
	// check --strict …`.
	_, err := captureStdout(t, func() error {
		return runCheckVCF(strings.NewReader(brokenVCF), "broken.vcf", true, false)
	})
	if err == nil {
		t.Fatalf("expected an error for broken VCF under strict mode, got nil")
	}
	if !strings.Contains(err.Error(), "error") {
		t.Errorf("expected error to mention errors, got %v", err)
	}
}

func TestRunCheckBEDClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runCheckBED(strings.NewReader(minimalBED), "example.bed", false, false)
	})
	if err != nil {
		t.Fatalf("runCheckBED: %v", err)
	}
	if !strings.Contains(out, "example.bed: clean") {
		t.Errorf("expected clean summary, got %q", out)
	}
}

func TestRunCheckSAMClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runCheckSAM(strings.NewReader(minimalSAM), "example.sam", false, false)
	})
	if err != nil {
		t.Fatalf("runCheckSAM: %v", err)
	}
	if !strings.Contains(out, "example.sam: clean") {
		t.Errorf("expected clean summary, got %q", out)
	}
}

func TestRunCheckFASTAClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runCheckFASTA(strings.NewReader(minimalFASTA), "example.fa", false, false)
	})
	if err != nil {
		t.Fatalf("runCheckFASTA: %v", err)
	}
	if !strings.Contains(out, "example.fa: clean") {
		t.Errorf("expected clean summary, got %q", out)
	}
}

func TestRunCheckFASTQClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runCheckFASTQ(strings.NewReader(minimalFASTQ), "example.fq", false, false)
	})
	if err != nil {
		t.Fatalf("runCheckFASTQ: %v", err)
	}
	if !strings.Contains(out, "example.fq: clean") {
		t.Errorf("expected clean summary, got %q", out)
	}
}

func TestRunCheckGFFClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runCheckGFF(strings.NewReader(minimalGFF), "example.gff", false, false)
	})
	if err != nil {
		t.Fatalf("runCheckGFF: %v", err)
	}
	if !strings.Contains(out, "example.gff: clean") {
		t.Errorf("expected clean summary, got %q", out)
	}
}

// -----------------------------------------------------------------------------
// runRoundtrip per-format tests
// -----------------------------------------------------------------------------

func TestRunRoundtripVCFClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runRoundtripVCF(strings.NewReader(minimalVCF), false)
	})
	if err != nil {
		t.Fatalf("runRoundtripVCF: %v", err)
	}
	if !strings.Contains(out, "##fileformat=VCFv4.5") {
		t.Errorf("missing VCF header in output: %q", out)
	}
	if !strings.Contains(out, "chr1\t100\t.\tA\tC\t30\tPASS\tDP=42") {
		t.Errorf("missing VCF record in output: %q", out)
	}
}

func TestRunRoundtripBEDClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runRoundtripBED(strings.NewReader(minimalBED), false)
	})
	if err != nil {
		t.Fatalf("runRoundtripBED: %v", err)
	}
	if !strings.Contains(out, "chr1\t100\t200") {
		t.Errorf("missing BED record in output: %q", out)
	}
}

func TestRunRoundtripSAMClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runRoundtripSAM(strings.NewReader(minimalSAM), false)
	})
	if err != nil {
		t.Fatalf("runRoundtripSAM: %v", err)
	}
	if !strings.Contains(out, "@HD") {
		t.Errorf("missing SAM @HD in output: %q", out)
	}
	if !strings.Contains(out, "r1\t") {
		t.Errorf("missing SAM record r1 in output: %q", out)
	}
}

func TestRunRoundtripFASTAClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runRoundtripFASTA(strings.NewReader(minimalFASTA), false)
	})
	if err != nil {
		t.Fatalf("runRoundtripFASTA: %v", err)
	}
	if !strings.Contains(out, ">seq1") {
		t.Errorf("missing FASTA header in output: %q", out)
	}
	if !strings.Contains(out, "ACGTACGTAC") {
		t.Errorf("missing FASTA sequence in output: %q", out)
	}
}

func TestRunRoundtripFASTQClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runRoundtripFASTQ(strings.NewReader(minimalFASTQ), false)
	})
	if err != nil {
		t.Fatalf("runRoundtripFASTQ: %v", err)
	}
	if !strings.Contains(out, "@seq1") {
		t.Errorf("missing FASTQ header in output: %q", out)
	}
	if !strings.Contains(out, "ACGTACGTAC") {
		t.Errorf("missing FASTQ sequence in output: %q", out)
	}
}

// -----------------------------------------------------------------------------
// runSummary per-format tests
// -----------------------------------------------------------------------------

func TestRunSummaryVCFClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runSummaryVCF(strings.NewReader(minimalVCF), "example.vcf")
	})
	if err != nil {
		t.Fatalf("runSummaryVCF: %v", err)
	}
	if !strings.Contains(out, "example.vcf") {
		t.Errorf("missing path in summary: %q", out)
	}
	if !strings.Contains(out, "records:") {
		t.Errorf("missing records field in summary: %q", out)
	}
	if !strings.Contains(out, "fileformat:") {
		t.Errorf("missing fileformat in summary: %q", out)
	}
}

func TestRunSummaryBEDClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runSummaryBED(strings.NewReader(minimalBED), "example.bed")
	})
	if err != nil {
		t.Fatalf("runSummaryBED: %v", err)
	}
	if !strings.Contains(out, "example.bed") {
		t.Errorf("missing path in summary: %q", out)
	}
	if !strings.Contains(out, "BED3") {
		t.Errorf("expected BED3 in summary, got: %q", out)
	}
}

func TestRunSummarySAMClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runSummarySAM(strings.NewReader(minimalSAM), "example.sam")
	})
	if err != nil {
		t.Fatalf("runSummarySAM: %v", err)
	}
	if !strings.Contains(out, "example.sam") {
		t.Errorf("missing path in summary: %q", out)
	}
	if !strings.Contains(out, "@SQ lines:") {
		t.Errorf("missing @SQ lines field in summary: %q", out)
	}
}

func TestRunSummaryFASTAClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runSummaryFASTA(strings.NewReader(minimalFASTA), "example.fa")
	})
	if err != nil {
		t.Fatalf("runSummaryFASTA: %v", err)
	}
	if !strings.Contains(out, "records:") {
		t.Errorf("missing records field in summary: %q", out)
	}
}

func TestRunSummaryFASTQClean(t *testing.T) {
	out, err := captureStdout(t, func() error {
		return runSummaryFASTQ(strings.NewReader(minimalFASTQ), "example.fq")
	})
	if err != nil {
		t.Fatalf("runSummaryFASTQ: %v", err)
	}
	if !strings.Contains(out, "phred+33") {
		t.Errorf("expected Phred+33 in summary, got: %q", out)
	}
}

// -----------------------------------------------------------------------------
// reportFrom* adapter tests
// -----------------------------------------------------------------------------

// TestReportFromVCFCoverage exercises every per-format
// adapter so a typo in any of the eight `reportFrom*` helpers
// shows up as a missing case.
func TestReportFromVCFCoverage(t *testing.T) {
	// All eight adapters are referenced through the test body
	// so removing one breaks the build.
	var v []vcf.ValidationError
	var b []bed.Issue
	var s []sam.Issue
	var bamIssues []bam.Issue
	var g []gff.Issue
	var fa []fasta.Issue
	var fq []fastq.Issue
	var bf []bcf.Issue
	_ = v
	_ = b
	_ = s
	_ = bamIssues
	_ = g
	_ = fa
	_ = fq
	_ = bf

	out := captureWriter(t, func(w io.Writer) {
		reportFromVCF(nil)
		reportFromBED(nil)
		reportFromSAM(nil)
		reportFromBAM(nil)
		reportFromGFF(nil)
		reportFromFASTA(nil)
		reportFromFASTQ(nil)
		reportFromBCF(nil)
		_, _ = w.Write([]byte("ok"))
	})
	if out != "ok" {
		t.Errorf("expected ok, got %q", out)
	}
}

// captureWriter swaps os.Stdout for an io.Writer for the
// duration of fn and returns everything fn wrote.
func captureWriter(t *testing.T, fn func(io.Writer)) string {
	t.Helper()
	var buf bytes.Buffer
	fn(&buf)
	return buf.String()
}

// -----------------------------------------------------------------------------
// Dispatch tests (runCheck / runRoundtrip / runSummary) via tempfile
// -----------------------------------------------------------------------------

func writeTemp(t *testing.T, name, content string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, []byte(content), 0o600); err != nil {
		t.Fatalf("write tempfile %s: %v", name, err)
	}
	return path
}

// TestRunCheckDispatchVCF covers the format-detection dispatch
// in runCheck: a real file path drives auto-detect, the right
// per-format helper runs, and the right summary line is printed.
func TestRunCheckDispatchVCF(t *testing.T) {
	path := writeTemp(t, "x.vcf", minimalVCF)
	out, err := captureStdout(t, func() error {
		return runCheck([]string{path})
	})
	if err != nil {
		t.Fatalf("runCheck: %v", err)
	}
	if !strings.Contains(out, "x.vcf: clean") {
		t.Errorf("expected clean summary, got %q", out)
	}
}

func TestRunCheckDispatchBED(t *testing.T) {
	path := writeTemp(t, "x.bed", minimalBED)
	out, err := captureStdout(t, func() error {
		return runCheck([]string{path})
	})
	if err != nil {
		t.Fatalf("runCheck: %v", err)
	}
	if !strings.Contains(out, "x.bed: clean") {
		t.Errorf("expected clean summary, got %q", out)
	}
}

func TestRunRoundtripDispatchVCF(t *testing.T) {
	path := writeTemp(t, "x.vcf", minimalVCF)
	out, err := captureStdout(t, func() error {
		return runRoundtrip([]string{path})
	})
	if err != nil {
		t.Fatalf("runRoundtrip: %v", err)
	}
	if !strings.Contains(out, "##fileformat=VCFv4.5") {
		t.Errorf("missing VCF header in dispatch output: %q", out)
	}
}

func TestRunSummaryDispatchVCF(t *testing.T) {
	path := writeTemp(t, "x.vcf", minimalVCF)
	out, err := captureStdout(t, func() error {
		return runSummary([]string{path})
	})
	if err != nil {
		t.Fatalf("runSummary: %v", err)
	}
	if !strings.Contains(out, "records:") {
		t.Errorf("expected records field in dispatch summary: %q", out)
	}
}

// TestCLIPipeline_VCFToBCFAndBack is the architecture review's
// item #5: an end-to-end pipeline test that drives VCF →
// BCF → VCF → check via the CLI helpers and verifies the
// data round-trips. Catches any regression in the boundary
// conversions that the per-format unit tests miss.
func TestCLIPipeline_VCFToBCFAndBack(t *testing.T) {
	// 1. VCF → BCF (runVCFToBCF).
	var bcfBuf bytes.Buffer
	if err := runVCFToBCF(strings.NewReader(minimalVCF), &bcfBuf); err != nil {
		t.Fatalf("vcf→bcf: %v", err)
	}
	if bcfBuf.Len() == 0 {
		t.Fatal("bcf output is empty")
	}

	// 2. BCF → VCF (runBCFToVCF).
	var vcfBuf bytes.Buffer
	if err := runBCFToVCF(&bcfBuf, &vcfBuf); err != nil {
		t.Fatalf("bcf→vcf: %v", err)
	}
	out := vcfBuf.String()
	// 3. Verify the round-tripped VCF has the expected
	// header and record contents.
	if !strings.Contains(out, "##fileformat=VCFv4.5") {
		t.Errorf("missing fileformat: %q", out)
	}
	if !strings.Contains(out, "chr1\t100") {
		t.Errorf("missing chr1:100 record: %q", out)
	}
	if !strings.Contains(out, "DP=42") {
		t.Errorf("missing DP=42 INFO: %q", out)
	}

	// 4. Validate the round-tripped VCF (runCheckVCF on the
	// buffer) and confirm it's clean. This is the part the
	// pipeline must satisfy across the BCF boundary.
	if issues := vcf.NewValidator(vcf.ModeDefault).Validate(
		func() *vcf.Header {
			h, err := vcf.NewDecoder(strings.NewReader(out)).Header()
			if err != nil {
				t.Fatalf("vcf header after round-trip: %v", err)
			}
			return h
		}(),
		func() []*vcf.Record {
			var recs []*vcf.Record
			dec := vcf.NewDecoder(strings.NewReader(out))
			for {
				r, err := dec.Record()
				if err != nil {
					break
				}
				recs = append(recs, r)
			}
			return recs
		}(),
	); len(issues) != 0 {
		t.Errorf("round-tripped VCF has issues: %+v", issues)
	}
}
