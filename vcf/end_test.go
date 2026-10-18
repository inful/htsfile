package vcf

import (
	"bytes"
	"strings"
	"testing"
)

// parseHeaderForTest is a small helper that decodes a VCF header
// from a string. Returns the *Header so the encoder can be
// driven against it.
func parseHeaderForTest(t *testing.T, hdrText string) *Header {
	t.Helper()
	dec := NewDecoder(strings.NewReader(hdrText))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	return hdr
}

func TestComputeEND_RefOnly(t *testing.T) {
	// POS=100, REF="ACG" → end = 100 + 3 - 1 = 102.
	rec := &Record{Pos: 100, Ref: "ACG"}
	if got := computeEND(rec); got != 102 {
		t.Errorf("computeEND: got %d, want 102", got)
	}
}

func TestComputeEND_SymbolicSVWithSVLEN(t *testing.T) {
	// POS=100, REF="A", ALT=<DEL>, SVLEN=5 → end = 100 + 5 = 105.
	rec := &Record{
		Pos:  100,
		Ref:  "A",
		Alt:  []string{"<DEL>"},
		Info: map[string]InfoValue{
			"SVLEN": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{5}}},
		},
	}
	if got := computeEND(rec); got != 105 {
		t.Errorf("computeEND: got %d, want 105", got)
	}
}

func TestComputeEND_StarAlleleWithLEN(t *testing.T) {
	// POS=100, REF="A", ALT=<*>, sample LEN=10 → end = 100 + 10 - 1 = 109.
	rec := &Record{
		Pos:  100,
		Ref:  "A",
		Alt:  []string{"<*>"},
		Samples: []Genotype{{Fields: map[string]FieldValue{
			"LEN": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{10}}},
		}}},
	}
	if got := computeEND(rec); got != 109 {
		t.Errorf("computeEND: got %d, want 109", got)
	}
}

func TestComputeEND_StarAlleleMaxOverSamples(t *testing.T) {
	// Two samples with different LENs: take the max.
	rec := &Record{
		Pos: 100,
		Ref: "A",
		Alt: []string{"<*>"},
		Samples: []Genotype{
			{Fields: map[string]FieldValue{
				"LEN": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{5}}},
			}},
			{Fields: map[string]FieldValue{
				"LEN": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{12}}},
			}},
		},
	}
	if got := computeEND(rec); got != 111 {
		t.Errorf("computeEND: got %d, want 111", got)
	}
}

func TestComputeEND_MaxOfAllSources(t *testing.T) {
	// POS=100, REF="AC", ALT=<DEL> with SVLEN=50, sample LEN=20.
	// Max of:
	//   POS + len(REF) - 1   = 101
	//   POS + SVLEN           = 150
	//   POS + LEN - 1         = 119
	// → 150.
	rec := &Record{
		Pos:  100,
		Ref:  "AC",
		Alt:  []string{"<DEL>"},
		Info: map[string]InfoValue{"SVLEN": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{50}}}},
		Samples: []Genotype{{Fields: map[string]FieldValue{
			"LEN": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{20}}},
		}}},
	}
	if got := computeEND(rec); got != 150 {
		t.Errorf("computeEND: got %d, want 150", got)
	}
}

func TestComputeEND_NonSymbolicAllele(t *testing.T) {
	// Non-symbolic ALT: END is just POS + len(REF) - 1.
	rec := &Record{Pos: 50, Ref: "ATCG", Alt: []string{"A"}}
	if got := computeEND(rec); got != 53 {
		t.Errorf("computeEND: got %d, want 53", got)
	}
}

func TestEncoder_ENDIsComputed(t *testing.T) {
	// Wire test: a record with END in its Info map should emit
	// the computed value, not whatever the caller stored.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"end position\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	hdr := parseHeaderForTest(t, hdrText)
	rec := &Record{
		Chrom: "chr1",
		Pos:   100,
		Ref:   "ACGT", // end-of-REF = 103
		Info: map[string]InfoValue{
			// Caller-supplied (wrong) value: should be
			// replaced with the computed 103.
			"END": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{9999}}},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	out := buf.String()
	if !strings.Contains(out, "END=103") {
		t.Errorf("output should contain END=103 (computed), got: %s", out)
	}
	if strings.Contains(out, "END=9999") {
		t.Errorf("output should not contain caller's stale END=9999, got: %s", out)
	}
}

func TestEncoder_ENDComputedWithSVLEN(t *testing.T) {
	// POS=100, REF="A", ALT=<DEL>, SVLEN=8 → end = 108.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"end position\">\n" +
		"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"SV length\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	hdr := parseHeaderForTest(t, hdrText)
	rec := &Record{
		Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"<DEL>"},
		Info: map[string]InfoValue{
			"SVLEN": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{8}}},
			"END":   {typedCore: typedCore{Type: TypeInteger, Ints: []int64{0}}},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if !strings.Contains(buf.String(), "END=108") {
		t.Errorf("output should contain END=108, got: %s", buf.String())
	}
}

func TestEncoder_RecordUnchanged(t *testing.T) {
	// Sanity: the encoder must not mutate the caller's Record.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"end position\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	hdr := parseHeaderForTest(t, hdrText)
	rec := &Record{
		Chrom: "chr1", Pos: 100, Ref: "ACGT",
		Info: map[string]InfoValue{
			"END": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{7777}}},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if got := rec.Info["END"].Ints[0]; got != 7777 {
		t.Errorf("encoder mutated caller's Record: END now %d, want 7777", got)
	}
}