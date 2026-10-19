package vcf

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

// encodeDecode round-trips src through the Encoder and back via the
// Decoder and returns the resulting Header and Records.
func encodeDecode(t *testing.T, hdr *Header, recs []*Record) (*Header, []*Record) {
	t.Helper()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatalf("EncodeRecord: %v", err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	dec := NewDecoder(&buf)
	hdr2, err := dec.Header()
	if err != nil {
		t.Fatalf("re-decode Header: %v\nbuf:\n%s", err, buf.String())
	}
	var out []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("re-decode Record: %v\nbuf:\n%s", err, buf.String())
		}
		out = append(out, rec)
	}
	return hdr2, out
}

func TestEncoder_Minimal(t *testing.T) {
	t.Parallel()
	h := &Header{
		FileFormat: FileFormat{4, 5},
		Infos: map[string]*InfoDef{
			"DP": {ID: "DP", Number: Number{Number: 1}, Type: TypeInteger, Description: "depth"},
		},
	}
	rec := &Record{
		Chrom: "1", Pos: 100, Ref: "A",
		Info: map[string]InfoValue{"DP": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{10}}}},
	}
	hdr2, recs := encodeDecode(t, h, []*Record{rec})
	if hdr2.FileFormat.String() != "VCFv4.5" {
		t.Errorf("FileFormat = %s", hdr2.FileFormat.String())
	}
	if len(recs) != 1 {
		t.Fatalf("len(recs) = %d", len(recs))
	}
	if v, ok := recs[0].Info["DP"]; !ok || v.Ints[0] != 10 {
		t.Errorf("DP = %+v", v)
	}
}

func TestEncoder_RoundTripSpecExample(t *testing.T) {
	t.Parallel()
	hdr, recs := parseFile(t, "testdata/spec_example.vcf")
	hdr2, recs2 := encodeDecode(t, hdr, recs)

	// Header sanity
	if hdr2.FileFormat.String() != "VCFv4.5" {
		t.Errorf("FileFormat = %s", hdr2.FileFormat.String())
	}
	if len(hdr2.SampleNames) != 3 {
		t.Errorf("SampleNames = %v", hdr2.SampleNames)
	}
	if len(hdr2.Infos) != len(hdr.Infos) {
		t.Errorf("Infos: got %d, want %d", len(hdr2.Infos), len(hdr.Infos))
	}
	if hdr2.Infos["DP"].Description != hdr.Infos["DP"].Description {
		t.Errorf("DP.Description: got %q, want %q", hdr2.Infos["DP"].Description, hdr.Infos["DP"].Description)
	}

	// Records count
	if len(recs2) != len(recs) {
		t.Fatalf("Records: got %d, want %d", len(recs2), len(recs))
	}
	for i := range recs {
		if recs[i].Chrom != recs2[i].Chrom {
			t.Errorf("rec%d.Chrom: %q vs %q", i, recs[i].Chrom, recs2[i].Chrom)
		}
		if recs[i].Pos != recs2[i].Pos {
			t.Errorf("rec%d.Pos: %d vs %d", i, recs[i].Pos, recs2[i].Pos)
		}
		if recs[i].Ref != recs2[i].Ref {
			t.Errorf("rec%d.Ref: %q vs %q", i, recs[i].Ref, recs2[i].Ref)
		}
		// Check INFO[key].Type and key presence.
		for k, v := range recs[i].Info {
			v2, ok := recs2[i].Info[k]
			if !ok {
				t.Errorf("rec%d lost INFO[%s]", i, k)
				continue
			}
			if v.Type != v2.Type {
				t.Errorf("rec%d INFO[%s] type: %v vs %v", i, k, v.Type, v2.Type)
			}
		}
		// Check sample columns survive.
		if len(recs[i].Samples) != len(recs2[i].Samples) {
			t.Errorf("rec%d Samples: %d vs %d", i, len(recs[i].Samples), len(recs2[i].Samples))
			continue
		}
		for j := range recs[i].Samples {
			g1 := recs[i].Samples[j].Fields["GT"]
			g2 := recs2[i].Samples[j].Fields["GT"]
			gt := g1.String()
			gt2 := g2.String()
			if gt != gt2 {
				t.Errorf("rec%d sample%d GT: %q vs %q", i, j, gt, gt2)
			}
		}
	}
}

func TestEncoder_RoundTripSVExample(t *testing.T) {
	t.Parallel()
	hdr, recs := parseFile(t, "testdata/sv_example.vcf")
	_, recs2 := encodeDecode(t, hdr, recs)

	if len(recs2) != len(recs) {
		t.Fatalf("Records: got %d, want %d", len(recs2), len(recs))
	}
	for i := range recs {
		if !equalStrings(recs[i].Alt, recs2[i].Alt) {
			t.Errorf("rec%d Alt: %v vs %v", i, recs[i].Alt, recs2[i].Alt)
		}
	}
	// Symbolic <DEL> record should survive.
	if recs2[1].Alt[0] != "<DEL>" {
		t.Errorf("rec1 Alt[0] = %q", recs2[1].Alt[0])
	}
	// Breakend records must survive byte-for-byte.
	if recs2[2].Alt[0] != "T[chrA:2[" {
		t.Errorf("rec2 Alt[0] = %q", recs2[2].Alt[0])
	}
	if recs2[3].Alt[0] != "]chrA:2]A" {
		t.Errorf("rec3 Alt[0] = %q", recs2[3].Alt[0])
	}
}

func TestEncoder_CRLF(t *testing.T) {
	t.Parallel()
	h := &Header{FileFormat: FileFormat{4, 5}}
	var buf bytes.Buffer
	enc := NewEncoderCRLF(&buf)
	if err := enc.EncodeHeader(h); err != nil {
		t.Fatal(err)
	}
	if !bytes.Contains(buf.Bytes(), []byte("\r\n")) {
		t.Errorf("output missing CRLF: %q", buf.String())
	}
	// Round-trip must accept CRLF.
	dec := NewDecoder(&buf)
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
}

func TestFormatVCFFloat(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in   float64
		want string
	}{
		{0, "0"},
		{1, "1"},
		{-1, "-1"},
		{0.5, "0.5"},
		{1e-7, "1e-07"},
		{positiveInf, "Infinity"},
		{negativeInf, "-Infinity"},
		{nanValue, "NaN"},
	}
	for _, c := range cases {
		if got := formatVCFFloat(c.in); got != c.want {
			t.Errorf("formatVCFFloat(%v) = %q; want %q", c.in, got, c.want)
		}
	}
}

func TestEncoder_PercentEncoding(t *testing.T) {
	t.Parallel()
	// INFO value containing ";" must be encoded so it survives re-parse.
	h := &Header{
		FileFormat: FileFormat{4, 5},
		Infos: map[string]*InfoDef{
			"NOTE": {ID: "NOTE", Number: Number{Number: 1}, Type: TypeString, Description: "note"},
		},
	}
	rec := &Record{
		Chrom: "1", Pos: 1, Ref: "A",
		Info: map[string]InfoValue{
			"NOTE": {typedCore: typedCore{Type: TypeString, Strs: []string{"a;b=c,d%"}}},
		},
	}
	_, recs := encodeDecode(t, h, []*Record{rec})
	if len(recs) != 1 {
		t.Fatalf("len(recs) = %d", len(recs))
	}
	v := recs[0].Info["NOTE"]
	if len(v.Strs) != 1 || v.Strs[0] != "a;b=c,d%" {
		t.Errorf("NOTE = %v", v)
	}
}

func TestEncoder_NoDuplicateHeaders(t *testing.T) {
	t.Parallel()
	hdr, _ := parseFile(t, "testdata/spec_example.vcf")

	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatal(err)
	}
	out := buf.String()

	// Each specific ##meta line should appear exactly once.
	specific := []string{
		"##fileformat=VCFv4.5",
		"##fileDate=20090805",
		"##source=myImputationProgramV3.1",
		"##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta",
		"##INFO=<ID=NS,",
		"##INFO=<ID=DP,",
		"##INFO=<ID=AF,",
		"##INFO=<ID=AA,",
		"##INFO=<ID=DB,",
		"##INFO=<ID=H2,",
		"##FILTER=<ID=q10,",
		"##FILTER=<ID=s50,",
		"##FORMAT=<ID=GT,",
		"##FORMAT=<ID=GQ,",
		"##FORMAT=<ID=DP,",
		"##FORMAT=<ID=HQ,",
		"##contig=<ID=20,",
	}
	for _, s := range specific {
		if c := strings.Count(out, s); c != 1 {
			t.Errorf("%q appears %d times in output (expected 1)", s, c)
		}
	}
}

func TestEncoder_EncodeRecord_QualityNilAndMissing(t *testing.T) {
	t.Parallel()
	q := 5.0
	h := &Header{FileFormat: FileFormat{4, 5}}
	recs := []*Record{
		{Chrom: "1", Pos: 1, Ref: "A", Qual: nil},
		{Chrom: "1", Pos: 2, Ref: "A", Qual: &q},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(h); err != nil {
		t.Fatal(err)
	}
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatal(err)
	}
	out := buf.String()
	if !strings.Contains(out, "\t.\t") {
		t.Errorf("output should contain a missing-quality '.':\n%s", out)
	}
	if !strings.Contains(out, "\t5\t") {
		t.Errorf("output should contain QUAL=5:\n%s", out)
	}
}

func TestEncoder_EncodeRecord_NilError(t *testing.T) {
	t.Parallel()
	enc := NewEncoder(io.Discard)
	if err := enc.EncodeRecord(nil); err == nil {
		t.Error("EncodeRecord(nil) should error")
	}
}

func TestEncoder_EncodeHeader_MissingFileFormat(t *testing.T) {
	t.Parallel()
	enc := NewEncoder(io.Discard)
	if err := enc.EncodeHeader(&Header{}); err == nil {
		t.Error("EncodeHeader with empty FileFormat should error")
	}
}

func TestEncoder_PedigreeAndSample(t *testing.T) {
	t.Parallel()
	h := &Header{
		FileFormat: FileFormat{4, 5},
		Pedigrees: []*Pedigree{
			{
				ID: "TumourSample", Original: "GermlineID",
				Father: "FatherID", Mother: "MotherID",
				Names: []string{"A1", "A2"},
			},
		},
		Samples: []*Sample{
			{
				ID:          "Sample1",
				Genomes:     map[string]string{"GRCh38": "/path/to.fa"},
				Description: "European example",
				DOI:         "url",
			},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(h); err != nil {
		t.Fatal(err)
	}
	out := buf.String()
	if !strings.Contains(out, `##PEDIGREE=<ID=TumourSample,Original=GermlineID,Father=FatherID,Mother=MotherID,Name_1=A1,Name_2=A2>`) {
		t.Errorf("pedigree missing or wrong:\n%s", out)
	}
	if !strings.Contains(out, `##SAMPLE=<ID=Sample1,Genomes=GRCh38="/path/to.fa",Description="European example",DOI="url">`) {
		t.Errorf("sample missing or wrong:\n%s", out)
	}
}

func TestEncoder_EncodeSample_DroppedTrailingFields(t *testing.T) {
	t.Parallel()
	h := &Header{FileFormat: FileFormat{4, 5}, SampleNames: []string{"S"}}
	rec := &Record{
		Chrom: "1", Pos: 1, Ref: "A",
		Format: []string{"GT", "DP", "HQ"},
		Samples: []Genotype{
			{
				Fields: map[string]FieldValue{
					"GT": {typedCore: typedCore{Type: TypeString, Strs: []string{"0/1"}}},
					// DP and HQ dropped — they are NOT in Fields. Per
					// spec §1.6.1, trailing fields may be dropped, and
					// the encoder must NOT re-emit them as "." in the
					// sample column. (The FORMAT column itself stays
					// as GT:DP:HQ because other samples in the same
					// record might still have all three.)
				},
			},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(h); err != nil {
		t.Fatal(err)
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatal(err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatal(err)
	}
	out := buf.String()
	if !strings.Contains(out, "GT:DP:HQ\t0/1\n") {
		t.Errorf("FORMAT column should reflect max keys, sample column should drop trailing fields:\n%s", out)
	}
	if strings.Contains(out, "0/1:.") {
		t.Errorf("sample column must not have trailing '.':\n%s", out)
	}
}

func TestEncoder_EncodeSample_TypeFallback(t *testing.T) {
	t.Parallel()
	// Unknown FORMAT key falls back to raw.
	h := &Header{FileFormat: FileFormat{4, 5}, SampleNames: []string{"S"}}
	rec := &Record{
		Chrom: "1", Pos: 1, Ref: "A",
		Format: []string{"ZZ"},
		Samples: []Genotype{
			{
				Fields: map[string]FieldValue{
					"ZZ": {typedCore: typedCore{Type: TypeUnknown, Raw: []string{"hello"}}},
				},
			},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(h); err != nil {
		t.Fatal(err)
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatal(err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(buf.String(), "ZZ\thello") {
		t.Errorf("expected raw 'hello' in output, got:\n%s", buf.String())
	}
}
