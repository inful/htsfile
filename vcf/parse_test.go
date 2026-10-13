package vcf

import (
	"errors"
	"io"
	"os"
	"strings"
	"testing"
)

// helper: parse a VCF from a string
func parseAll(t *testing.T, src string) (*Header, []*Record) {
	t.Helper()
	dec := NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header(): %v", err)
	}
	var recs []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record(): %v", err)
		}
		recs = append(recs, rec)
	}
	return hdr, recs
}

// helper: parse from a file path
func parseFile(t *testing.T, path string) (*Header, []*Record) {
	t.Helper()
	f, err := os.Open(path) //nolint:gosec // path is test-controlled
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close() //nolint:errcheck // test helper
	dec := NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header(): %v", err)
	}
	var recs []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record(): %v", err)
		}
		recs = append(recs, rec)
	}
	return hdr, recs
}

func TestDecoder_MinimalHeader(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	hdr, recs := parseAll(t, src)

	if hdr.FileFormat.String() != "VCFv4.5" {
		t.Errorf("FileFormat = %s", hdr.FileFormat.String())
	}
	if dp, ok := hdr.Infos["DP"]; !ok {
		t.Error("DP INFO missing")
	} else if dp.Number.Number != 1 || dp.Type != TypeInteger {
		t.Errorf("DP = %+v", dp)
	}
	if len(recs) != 0 {
		t.Errorf("expected 0 records, got %d", len(recs))
	}
	if len(hdr.SampleNames) != 0 {
		t.Errorf("expected no sample names, got %v", hdr.SampleNames)
	}
}

func TestDecoder_SpecExample(t *testing.T) {
	t.Parallel()
	hdr, recs := parseFile(t, "testdata/spec_example.vcf")

	if hdr.FileFormat.String() != "VCFv4.5" {
		t.Errorf("FileFormat = %s", hdr.FileFormat.String())
	}
	if hdr.FileDate != "20090805" {
		t.Errorf("FileDate = %q", hdr.FileDate)
	}
	if hdr.Source != "myImputationProgramV3.1" {
		t.Errorf("Source = %q", hdr.Source)
	}
	if hdr.Reference != "file:///seq/references/1000GenomesPilot-NCBI36.fasta" {
		t.Errorf("Reference = %q", hdr.Reference)
	}

	// ##phasing is preserved as a generic meta line.
	var foundPhasing bool
	for _, m := range hdr.MetaLines {
		if m.Key == "phasing" {
			foundPhasing = true
			if m.Value != "partial" {
				t.Errorf("phasing value = %q", m.Value)
			}
		}
	}
	if !foundPhasing {
		t.Error("phasing meta line missing")
	}

	// Contig
	if len(hdr.Contigs) != 1 {
		t.Fatalf("Contigs = %d entries", len(hdr.Contigs))
	}
	c := hdr.Contigs[0]
	if c.ID != "20" {
		t.Errorf("Contig.ID = %q", c.ID)
	}
	if c.Length != 62435964 {
		t.Errorf("Contig.Length = %d", c.Length)
	}
	if c.MD5 != "f126cdf8a6e0c7f379d618ff66beb2da" {
		t.Errorf("Contig.MD5 = %q", c.MD5)
	}

	// Sample names
	wantSamples := []string{"NA00001", "NA00002", "NA00003"}
	if len(hdr.SampleNames) != len(wantSamples) {
		t.Fatalf("SampleNames len = %d; want %d", len(hdr.SampleNames), len(wantSamples))
	}
	for i, s := range wantSamples {
		if hdr.SampleNames[i] != s {
			t.Errorf("SampleNames[%d] = %q; want %q", i, hdr.SampleNames[i], s)
		}
	}

	// INFO defs
	if len(hdr.Infos) != 6 {
		t.Errorf("Infos count = %d; want 6", len(hdr.Infos))
	}
	dp := hdr.Infos["DP"]
	if dp == nil || dp.Type != TypeInteger || dp.Number.Number != 1 {
		t.Errorf("DP = %+v", dp)
	}
	af := hdr.Infos["AF"]
	if af == nil || af.Type != TypeFloat || af.Number.Fixed != SymbolA {
		t.Errorf("AF = %+v", af)
	}
	db := hdr.Infos["DB"]
	if db == nil || db.Type != TypeFlag || db.Number.Number != 0 {
		t.Errorf("DB = %+v", db)
	}

	// FILTER defs
	if len(hdr.Filters) != 2 {
		t.Errorf("Filters count = %d; want 2", len(hdr.Filters))
	}
	if hdr.Filters["q10"] == nil {
		t.Error("q10 filter missing")
	}

	// FORMAT defs
	if len(hdr.Formats) != 4 {
		t.Errorf("Formats count = %d; want 4", len(hdr.Formats))
	}
	if hdr.Formats["HQ"].Number.Number != 2 {
		t.Errorf("HQ.Number = %+v", hdr.Formats["HQ"].Number)
	}

	// Records
	if len(recs) != 5 {
		t.Fatalf("Records = %d; want 5", len(recs))
	}

	// First record: 20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2
	r0 := recs[0]
	if r0.Chrom != "20" || r0.Pos != 14370 {
		t.Errorf("rec0 chrom/pos = %q/%d", r0.Chrom, r0.Pos)
	}
	if len(r0.ID) != 1 || r0.ID[0] != "rs6054257" {
		t.Errorf("rec0 ID = %v", r0.ID)
	}
	if r0.Ref != "G" {
		t.Errorf("rec0 Ref = %q", r0.Ref)
	}
	if len(r0.Alt) != 1 || r0.Alt[0] != "A" {
		t.Errorf("rec0 Alt = %v", r0.Alt)
	}
	if r0.Qual == nil || *r0.Qual != 29 {
		t.Errorf("rec0 Qual = %v", r0.Qual)
	}
	if len(r0.Filter) != 1 || r0.Filter[0] != "PASS" {
		t.Errorf("rec0 Filter = %v", r0.Filter)
	}
	if len(r0.Format) != 4 {
		t.Errorf("rec0 Format = %v", r0.Format)
	}
	if len(r0.Samples) != 3 {
		t.Errorf("rec0 Samples len = %d", len(r0.Samples))
	}
	if v, ok := r0.Info["DP"]; !ok || v.Type != TypeInteger || len(v.Ints) != 1 || v.Ints[0] != 14 {
		t.Errorf("rec0 Info[DP] = %+v", v)
	}
	if v, ok := r0.Info["AF"]; !ok || v.Type != TypeFloat || len(v.Floats) != 1 || v.Floats[0] != 0.5 {
		t.Errorf("rec0 Info[AF] = %+v", v)
	}
	if v, ok := r0.Info["DB"]; !ok || v.Type != TypeFlag || !v.Flag {
		t.Errorf("rec0 Info[DB] = %+v", v)
	}
	if v, ok := r0.Info["H2"]; !ok || v.Type != TypeFlag || !v.Flag {
		t.Errorf("rec0 Info[H2] = %+v", v)
	}

	// First sample's GT
	s0 := r0.Samples[0].Fields["GT"]
	if len(s0.Strs) != 1 || s0.Strs[0] != "0|0" {
		t.Errorf("rec0 sample0 GT = %+v", s0)
	}

	// Record 2 (microsat) has GT:GQ:DP — only 3 FORMAT keys.
	r4 := recs[4]
	if len(r4.Format) != 3 {
		t.Errorf("rec4 Format = %v", r4.Format)
	}
	if r4.Format[0] != "GT" || r4.Format[1] != "GQ" || r4.Format[2] != "DP" {
		t.Errorf("rec4 Format = %v", r4.Format)
	}
	// Last sample dropped trailing field; sample still present.
	if len(r4.Samples) != 3 {
		t.Errorf("rec4 Samples len = %d", len(r4.Samples))
	}
}

func TestDecoder_CRLF(t *testing.T) {
	t.Parallel()
	// CR+LF is permitted by the spec.
	src := "##fileformat=VCFv4.5\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\n1\t2\t.\tA\tC\t.\t.\t.\r\n"
	hdr, recs := parseAll(t, src)
	if hdr.FileFormat.Major != 4 {
		t.Errorf("Major = %d", hdr.FileFormat.Major)
	}
	if len(recs) != 1 {
		t.Errorf("Records = %d", len(recs))
	}
	if recs[0].Chrom != "1" {
		t.Errorf("Chrom = %q", recs[0].Chrom)
	}
}

func TestDecoder_HeaderOnly(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	hdr, recs := parseAll(t, src)
	if hdr == nil {
		t.Fatal("Header is nil")
	}
	if len(recs) != 0 {
		t.Errorf("expected 0 records, got %d", len(recs))
	}
}

func TestDecoder_EOF(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if !errors.Is(err, io.EOF) {
		t.Fatalf("first Record: got %v; want io.EOF", err)
	}
	if rec != nil {
		t.Errorf("rec should be nil at EOF, got %+v", rec)
	}
	// Calling Record again at EOF is allowed; must keep returning io.EOF.
	rec2, err2 := dec.Record()
	if !errors.Is(err2, io.EOF) {
		t.Fatalf("second Record: got %v; want io.EOF", err2)
	}
	if rec2 != nil {
		t.Errorf("rec2 should be nil at EOF")
	}
}

func TestDecoder_CachesHeader(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(src))
	h1, err := dec.Header()
	if err != nil {
		t.Fatal(err)
	}
	h2, err := dec.Header()
	if err != nil {
		t.Fatal(err)
	}
	if h1 != h2 {
		t.Error("Header() should return the same pointer on repeated calls")
	}
}

func TestDecoder_SVExample(t *testing.T) {
	t.Parallel()
	hdr, recs := parseFile(t, "testdata/sv_example.vcf")

	if hdr.FileFormat.String() != "VCFv4.5" {
		t.Errorf("FileFormat = %s", hdr.FileFormat.String())
	}

	// ALT defs
	wantAlts := []string{"INV", "INS", "DUP", "DUP:TANDEM", "DEL", "CNV"}
	if len(hdr.Alts) != len(wantAlts) {
		t.Errorf("Alts count = %d; want %d", len(hdr.Alts), len(wantAlts))
	}
	for _, id := range wantAlts {
		if hdr.Alts[id] == nil {
			t.Errorf("ALT %s missing", id)
		}
	}

	// Symbolic ALT parsing: <DEL> on rec[1]
	r1 := recs[1]
	if len(r1.Alt) != 1 || r1.Alt[0] != "<DEL>" {
		t.Errorf("rec1 Alt = %v; want [<DEL>]", r1.Alt)
	}

	// Breakend ALT on rec[2]: T[chrA:5[
	r2 := recs[2]
	if len(r2.Alt) != 1 || r2.Alt[0] != "T[chrA:5[" {
		t.Errorf("rec2 Alt = %v", r2.Alt)
	}

	// Reverse-direction breakend: ]chrA:2]A on rec[3]
	r3 := recs[3]
	if len(r3.Alt) != 1 || r3.Alt[0] != "]chrA:2]A" {
		t.Errorf("rec3 Alt = %v", r3.Alt)
	}

	// Symbolic <INS> on rec[7] (chrA 14 T <INS> ...)
	r7 := recs[7]
	if len(r7.Alt) != 1 || r7.Alt[0] != "<INS>" {
		t.Errorf("rec7 Alt = %v", r7.Alt)
	}

	// Number=A INFO key (MATEID with Number=A) — values must align with ALT count
	if v, ok := r2.Info["MATEID"]; !ok || v.Type != TypeString || len(v.Strs) != 1 || v.Strs[0] != "delbp2" {
		t.Errorf("rec2 MATEID = %+v", v)
	}

	// CIPOS with Number=.: present on rec[6] (chrA 5 G <DUP> ...)
	r6 := recs[6]
	if v, ok := r6.Info["CIPOS"]; !ok || v.Type != TypeInteger || len(v.Ints) != 2 {
		t.Errorf("rec6 CIPOS = %+v", v)
	}

	// IMPRECISE flag
	if v, ok := r7.Info["IMPRECISE"]; !ok || v.Type != TypeFlag || !v.Flag {
		t.Errorf("rec7 IMPRECISE = %+v", v)
	}
}

func TestSplitAlleles(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in   string
		want []string
	}{
		{"A,C,G", []string{"A", "C", "G"}},
		{".", []string{"."}}, // (caller handles '.' specially; we still test splitAlleles)
		{"T[chr1:5[", []string{"T[chr1:5["}},
		{"T[chr1:5[,C", []string{"T[chr1:5[", "C"}},
		{"]chr1:2]A", []string{"]chr1:2]A"}},
		{"A,]chr1:2]A", []string{"A", "]chr1:2]A"}},
		{"<DEL>,<INS>", []string{"<DEL>", "<INS>"}},
	}
	for _, c := range cases {
		got := splitAlleles(c.in)
		if !equalStrings(got, c.want) {
			t.Errorf("splitAlleles(%q) = %v; want %v", c.in, got, c.want)
		}
	}
}

func equalStrings(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestDecoder_Errors(t *testing.T) {
	t.Parallel()
	cases := []struct {
		name string
		src  string
	}{
		{
			"missing fileformat",
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		},
		{
			"fileformat not first",
			"##source=foo\n##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		},
		{
			"missing header line",
			"##fileformat=VCFv4.5\n",
		},
		{
			"bad fileformat value",
			"##fileformat=VCFbogus\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		},
		{
			"malformed structured line (no closing >)",
			"##fileformat=VCFv4.5\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		},
		{
			"structured line missing required ID",
			"##fileformat=VCFv4.5\n##INFO=<Number=1,Type=Integer,Description=\"d\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		},
		{
			"data line too few columns",
			"##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t2\t.\tA\n",
		},
		{
			"data line non-numeric POS",
			"##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\tnot_a_number\t.\tA\tC\t.\t.\t.\n",
		},
		{
			"unstructured line with empty value",
			"##fileformat=VCFv4.5\n##source=\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		},
		{
			"unstructured line value starting with <",
			"##fileformat=VCFv4.5\n##file=<x\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			t.Parallel()
			dec := NewDecoder(strings.NewReader(c.src))
			_, err := dec.Header()
			if err == nil {
				// Header() might succeed; the error might surface in Record().
				// Try reading records.
				for {
					_, recErr := dec.Record()
					if recErr == nil {
						continue
					}
					if errors.Is(recErr, io.EOF) {
						t.Fatal("expected an error, got EOF")
					}
					err = recErr
					break
				}
			}
			if err == nil {
				t.Fatal("expected an error, got nil")
			}
		})
	}
}
