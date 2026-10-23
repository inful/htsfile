package vcf

import (
	"bytes"
	"strings"
	"testing"
)

// The spec's BCF2 header mentions a string-dictionary pragma:
//
//	"BCF2 supports encoding values in a dictionary of strings.
//	 The string map is provided by the keyword
//	 ##dictionary=S0,S1,...,SN as a comma-separate ordered list
//	 of strings." (VCFv4.5, tex line 2100)
//
// The "Dictionary of strings" section that follows (tex lines
// 2103-2143) then describes the *header* dictionary: the
// position of each FILTER / INFO / FORMAT key in the header,
// which is what a BCF2 record's typed integers refer to. That
// is the encoding we implement (see T6.11 in
// SPEC_COMPLIANCE.md); the `s0`, `s1` in the section's example
// are a visual guide, not the wire format.
//
// The pragma itself has no wire encoding in VCFv4.5, but a
// header line carrying it is valid, so it must be accepted and
// preserved rather than rejected.

// dictionaryVCF is a VCF whose header carries the pragma.
var dictionaryVCF = "##fileformat=VCFv4.5\n" +
	"##dictionary=S0,S1,S2\n" +
	"##contig=<ID=20,length=1000>\n" +
	"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"x\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
	"20\t100\t.\tG\tA\t30\tPASS\tNS=3\n"

// TestHeader_DictionaryPragmaIsAccepted checks that a header
// carrying the pragma parses and keeps the line.
func TestHeader_DictionaryPragmaIsAccepted(t *testing.T) {
	dec := NewDecoder(strings.NewReader(dictionaryVCF))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	var found bool
	for _, ml := range hdr.MetaLines {
		if ml.Key == "dictionary" {
			found = true
			if ml.Value != "S0,S1,S2" {
				t.Errorf("dictionary: got %q want %q", ml.Value, "S0,S1,S2")
			}
		}
	}
	if !found {
		t.Error("the ##dictionary line was dropped from the header")
	}
}

// TestHeader_DictionaryPragmaSurvivesEncoding checks that the
// pragma is replayed when the header is written out again.
func TestHeader_DictionaryPragmaSurvivesEncoding(t *testing.T) {
	dec := NewDecoder(strings.NewReader(dictionaryVCF))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	if !strings.Contains(buf.String(), "##dictionary=S0,S1,S2") {
		t.Errorf("the encoded header lost the pragma:\n%s", buf.String())
	}
}
