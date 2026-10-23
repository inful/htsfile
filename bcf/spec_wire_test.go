package bcf

import (
	"bytes"
	"encoding/hex"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
)

// This file is the BCF wire-format guard.
//
// Every byte sequence below was produced by the reference
// implementation (`bcftools view -Ob`, htslib 1.24) from the
// VCFv4.5 spec's own worked example, and is reproduced here as
// a literal hex string. The test asserts that OUR encoder emits
// exactly the same bytes.
//
// The point of using literal bytes rather than the package's own
// encoder is that a fixture built by the code under test cannot
// catch a change to it. That failure mode is not hypothetical:
// before 2026-10-22 the BCF record layout, the string encoding,
// the dictionary indices and the GT allele vector were all wrong
// in ways that every round-trip test in this package happily
// passed, because the writer and the reader agreed with each
// other and with nobody else. htslib rejected those files with
// "Invalid BCF2 magic string", "Invalid REF/ALT type 9",
// "Invalid ID type 0 (null)" and "Invalid FORMAT id 12".
//
// Regenerate the reference bytes with:
//
//	bcftools view -Ob -o ref.bcf vcf/testdata/spec_example.vcf
//	python3 - <<'PY'
//	import gzip, struct
//	raw = gzip.open('ref.bcf','rb').read()
//	lt = struct.unpack('<I', raw[5:9])[0]; off = 9+lt
//	ls, li = struct.unpack('<II', raw[off:off+8])
//	print(raw[off+8:off+8+ls].hex())
//	print(raw[off+8+ls:off+8+ls+li].hex())
//	PY

// refRecord1Shared is the l_shared block of the spec example's
// first record, as written by htslib.
//
// Decoded, per VCFv4.5 §BcfSiteEncoding (tex lines 2176-2197):
//
//	00000000        CHROM = 0        (contig offset)
//	21380000        POS   = 14369    (0-based)
//	01000000        rlen  = 1
//
// 0000e841        QUAL  = 29.0f
//
//	0500            n_info   = 5
//	0200            n_allele = 2  (REF + 1 ALT)
//	030000          n_sample = 3
//	04              n_fmt    = 4
//	97 727336303534323537   ID   "rs6054257"  (0x97 = size 9, type 7)
//	17 47                   REF  "G"
//	17 41                   ALT  "A"
//	11 00                   FILTER = index 0 (PASS)
//	11 01                   INFO key 1  (NS)
//	11 03                   NS = 3
//	11 02                   INFO key 2  (DP)
//	11 0e                   DP = 14
//	11 03                   INFO key 3  (AF)
//	15 0000003f             AF = 0.5f
//	11 05                   INFO key 5  (DB)
//	11 06                   DB = present (Flag)
const refRecord1Shared = "0000000021380000010000000000e841050002000300000497727336303534323537" +
	"174717411100110111031102110e1103150000003f110500110600"

// refRecord1Indiv is the l_indiv block of the spec example's
// first record, as written by htslib.
//
// Per VCFv4.5 §GenotypeEncoding (tex line 2233) each FORMAT
// field is a triplet: fmt_key (typed int, the header dictionary
// index), fmt_type (descriptor byte), fmt_values.
//
//	11 09   fmt_key 9  = GT
//	21       fmt_type  = size 2, int8
//	03 03 05 03 04 04   GT alleles for 3 samples
//	11 0a   fmt_key 10 = GQ
//	11       fmt_type  = size 1, int8
//	30 30 2b             GQ = 48,48,43
//	11 02   fmt_key 2  = DP (the same dictionary index as the
//	                      INFO DP - the index is per distinct ID)
//	11       fmt_type  = size 1, int8
//	01 08 05             DP = 1,8,5
//	11 0b   fmt_key 11 = HQ
//	21       fmt_type  = size 2, int8
//	33 33 33 33 80 80   HQ = 51,51,51,51,<missing>,<missing>
const refRecord1Indiv = "110921030305030404110a1130302b110211010805110b21333333338080"

// TestEncoderMatchesReferenceBytes asserts that our encoder
// produces the exact same record bytes as htslib for the spec
// example's first record.
//
// This is the mutation-testable guard for the whole BCF layout.
// Re-introducing any of the fixed defects makes it fail:
//
//   - swapping the n_sample/n_fmt packing
//   - writing a bare type code instead of a (size<<4)|type
//     descriptor
//   - writing a string as type 7 + uint32 length instead of a
//     char vector
//   - using the FORMAT index instead of the header dictionary
//     index
//   - writing all fmt_keys first instead of one triplet per field
//   - encoding GT as a string instead of an allele vector
func TestEncoderMatchesReferenceBytes(t *testing.T) {
	// Locate the spec example.
	exPath := filepath.Join("..", "vcf", "testdata", "spec_example.vcf")
	raw, err := os.ReadFile(exPath)
	if err != nil {
		t.Skipf("spec example not available: %v", err)
	}
	hdr := ParseHeader(string(raw))
	if hdr == nil {
		t.Fatal("ParseHeader returned nil")
	}
	// The header text the encoder writes must match the input
	// header, otherwise the dictionary indices shift.
	rec := &Record{
		CHROM: "20", Pos: 14370, ID: "rs6054257",
		Ref: "G", Alt: []string{"A"},
		Qual: 29, QualKnown: true,
		Filter: []string{"PASS"},
		Info: map[string]InfoValue{
			"NS": {Type: "Integer", Ints: []int64{3}},
			"DP": {Type: "Integer", Ints: []int64{14}},
			"AF": {Type: "Float", Floats: []float64{0.5}},
			"DB": {Type: "Flag", Flag: true},
			"H2": {Type: "Flag", Flag: true},
		},
		Format: []uint16{0, 1, 2, 3}, // GT, GQ, DP, HQ
		Samples: []Sample{
			{Fields: map[uint16]SampleValue{
				0: {Type: "String", Strs: []string{"0|0"}},
				1: {Type: "Integer", Ints: []int32{48}},
				2: {Type: "Integer", Ints: []int32{1}},
				3: {Type: "Integer", Ints: []int32{51, 51}},
			}},
			{Fields: map[uint16]SampleValue{
				0: {Type: "String", Strs: []string{"1|0"}},
				1: {Type: "Integer", Ints: []int32{48}},
				2: {Type: "Integer", Ints: []int32{8}},
				3: {Type: "Integer", Ints: []int32{51, 51}},
			}},
			{Fields: map[uint16]SampleValue{
				0: {Type: "String", Strs: []string{"1/1"}},
				1: {Type: "Integer", Ints: []int32{43}},
				2: {Type: "Integer", Ints: []int32{5}},
				3: {Type: "Integer", Ints: []int32{0, 0}, IntMissing: []bool{true, true}},
			}},
		},
	}
	shared, err := encodeRecord(rec, 0, hdr, 2)
	if err != nil {
		t.Fatalf("encodeRecord: %v", err)
	}
	indiv, err := encodeIndivid(rec, hdr, 2)
	if err != nil {
		t.Fatalf("encodeIndivid: %v", err)
	}
	if got, want := hex.EncodeToString(shared), refRecord1Shared; got != want {
		t.Errorf("shared block does not match htslib.\n got: %s\nwant: %s", got, want)
	}
	if got, want := hex.EncodeToString(indiv), refRecord1Indiv; got != want {
		t.Errorf("indiv block does not match htslib.\n got: %s\nwant: %s", got, want)
	}
}

// TestDecoderReadsReferenceBytes asserts that our decoder reads
// the reference bytes correctly - the mirror of
// TestEncoderMatchesReferenceBytes.
func TestDecoderReadsReferenceBytes(t *testing.T) {
	shared, err := hex.DecodeString(refRecord1Shared)
	if err != nil {
		t.Fatal(err)
	}
	indiv, err := hex.DecodeString(refRecord1Indiv)
	if err != nil {
		t.Fatal(err)
	}
	hdrText := "##fileformat=VCFv4.5\n" +
		"##contig=<ID=20,length=62435964>\n" +
		"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"x\">\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n" +
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">\n" +
		"##INFO=<ID=AA,Number=1,Type=String,Description=\"x\">\n" +
		"##INFO=<ID=DB,Number=0,Type=Flag,Description=\"x\">\n" +
		"##INFO=<ID=H2,Number=0,Type=Flag,Description=\"x\">\n" +
		"##FILTER=<ID=q10,Description=\"x\">\n" +
		"##FILTER=<ID=s50,Description=\"x\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">\n" +
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"x\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n" +
		"##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n"
	// Build a one-record BCF stream around the reference bytes.
	// The stream is BGZF-framed, exactly as the encoder writes
	// it: magic, l_text (counting the NUL), the header text,
	// then the record.
	var rec bytes.Buffer
	zw := bgzf.NewWriter(&rec)
	write := func(b []byte) {
		t.Helper()
		if _, err := zw.Write(b); err != nil {
			t.Fatalf("bgzf write: %v", err)
		}
	}
	write(Magic[:])
	var lt [4]byte
	putU32(lt[:], uint32(len(hdrText)+1))
	write(lt[:])
	write([]byte(hdrText))
	write([]byte{0})
	var lens [8]byte
	putU32(lens[0:4], uint32(len(shared)))
	putU32(lens[4:8], uint32(len(indiv)))
	write(lens[:])
	write(shared)
	write(indiv)
	if err := zw.Close(); err != nil {
		t.Fatalf("bgzf close: %v", err)
	}

	d := NewDecoder(bytes.NewReader(rec.Bytes()))
	if _, err := d.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	got, err := d.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if got.CHROM != "20" || got.Pos != 14370 {
		t.Errorf("CHROM/POS: got %s:%d, want 20:14370", got.CHROM, got.Pos)
	}
	if got.ID != "rs6054257" || got.Ref != "G" || len(got.Alt) != 1 || got.Alt[0] != "A" {
		t.Errorf("ID/REF/ALT: got %q %q %v", got.ID, got.Ref, got.Alt)
	}
	if got.Qual != 29 || !got.QualKnown {
		t.Errorf("QUAL: got %v (known=%v), want 29", got.Qual, got.QualKnown)
	}
	if len(got.Filter) != 1 || got.Filter[0] != "PASS" {
		t.Errorf("Filter: got %v, want [PASS]", got.Filter)
	}
	if v, ok := got.Info["NS"]; !ok || len(v.Ints) != 1 || v.Ints[0] != 3 {
		t.Errorf("INFO NS: got %+v", v)
	}
	if v, ok := got.Info["DP"]; !ok || len(v.Ints) != 1 || v.Ints[0] != 14 {
		t.Errorf("INFO DP: got %+v", v)
	}
	if v, ok := got.Info["AF"]; !ok || len(v.Floats) != 1 || v.Floats[0] != 0.5 {
		t.Errorf("INFO AF: got %+v", v)
	}
	if v, ok := got.Info["DB"]; !ok || !v.Flag {
		t.Errorf("INFO DB: got %+v", v)
	}
	if len(got.Samples) != 3 {
		t.Fatalf("Samples: got %d, want 3", len(got.Samples))
	}
	// The GT allele vector must decode back to the genotype
	// strings, and the per-slot missing HQ must survive.
	wantGT := []string{"0|0", "1|0", "1/1"}
	for i, s := range got.Samples {
		if v, ok := s.Fields[0]; !ok || len(v.Strs) != 1 || v.Strs[0] != wantGT[i] {
			t.Errorf("sample %d GT: got %+v, want %q", i, v, wantGT[i])
		}
	}
	if v, ok := got.Samples[2].Fields[3]; !ok || len(v.Ints) != 2 ||
		!v.IntMissing[0] || !v.IntMissing[1] {
		t.Errorf("sample 2 HQ: got %+v, want two missing slots", v)
	}
}

// putU32 writes a little-endian uint32.
func putU32(b []byte, v uint32) {
	b[0] = byte(v)
	b[1] = byte(v >> 8)
	b[2] = byte(v >> 16)
	b[3] = byte(v >> 24)
}

// TestReferenceBytesAreCurrent guards the reference hex strings
// against accidental edits: they must decode and must have the
// lengths the spec's example implies.
func TestReferenceBytesAreCurrent(t *testing.T) {
	shared, err := hex.DecodeString(refRecord1Shared)
	if err != nil {
		t.Fatalf("refRecord1Shared is not valid hex: %v", err)
	}
	indiv, err := hex.DecodeString(refRecord1Indiv)
	if err != nil {
		t.Fatalf("refRecord1Indiv is not valid hex: %v", err)
	}
	if len(shared) != 61 {
		t.Errorf("refRecord1Shared: %d bytes, want 61", len(shared))
	}
	if len(indiv) != 30 {
		t.Errorf("refRecord1Indiv: %d bytes, want 30", len(indiv))
	}
	if !strings.Contains(refRecord1Shared, "727336303534323537") {
		t.Error("refRecord1Shared should contain the ID \"rs6054257\"")
	}
}
