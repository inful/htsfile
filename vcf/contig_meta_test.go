package vcf

import (
	"bytes"
	"strings"
	"testing"
)

// TestParseContig_SpeciesAndTaxonomy covers the optional
// species and taxonomy sub-fields of a ##contig line
// (per VCFv4.5 §1.4.7). They must be parsed into the
// Contig struct and round-trip cleanly through the
// encoder.
func TestParseContig_SpeciesAndTaxonomy(t *testing.T) {
	in := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000,species=\"Homo sapiens\",taxonomy=\"NCBI:9606\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(in))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if len(hdr.Contigs) != 1 {
		t.Fatalf("Contigs: %d entries", len(hdr.Contigs))
	}
	c := hdr.Contigs[0]
	if c.Species != "Homo sapiens" {
		t.Errorf("Contig.Species = %q, want %q", c.Species, "Homo sapiens")
	}
	if c.Taxonomy != "NCBI:9606" {
		t.Errorf("Contig.Taxonomy = %q, want %q", c.Taxonomy, "NCBI:9606")
	}
	// Round-trip: emit and re-parse.
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	out := buf.String()
	if !strings.Contains(out, "species=\"Homo sapiens\"") {
		t.Errorf("output missing species; got %q", out)
	}
	if !strings.Contains(out, "taxonomy=\"NCBI:9606\"") {
		t.Errorf("output missing taxonomy; got %q", out)
	}
	dec2 := NewDecoder(strings.NewReader(out))
	hdr2, err := dec2.Header()
	if err != nil {
		t.Fatalf("re-parse Header: %v", err)
	}
	if len(hdr2.Contigs) != 1 {
		t.Fatalf("re-parse Contigs: %d entries", len(hdr2.Contigs))
	}
	if hdr2.Contigs[0].Species != "Homo sapiens" || hdr2.Contigs[0].Taxonomy != "NCBI:9606" {
		t.Errorf("re-parse: Species=%q Taxonomy=%q", hdr2.Contigs[0].Species, hdr2.Contigs[0].Taxonomy)
	}
}
