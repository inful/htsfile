package vcf

import (
	"testing"
)

func TestParseFileFormat(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in      string
		want    FileFormat
		wantErr bool
	}{
		{"VCFv4.5", FileFormat{4, 5}, false},
		{"VCFv4.4", FileFormat{4, 4}, false},
		{"VCFv5.0", FileFormat{5, 0}, true},     // major != 4
		{"VCFv10.20", FileFormat{10, 20}, true}, // major != 4
		{"vcfv4.5", FileFormat{4, 5}, false},    // case insensitive on prefix
		{"BCFv2.2", FileFormat{}, true},         // wrong prefix
		{"VCFv4", FileFormat{}, true},           // missing minor
		{"VCFv4.", FileFormat{}, true},          // empty minor
		{"VCFvx.y", FileFormat{}, true},         // non-numeric major
	}
	for _, c := range cases {
		t.Run(c.in, func(t *testing.T) {
			t.Parallel()
			got, err := ParseFileFormat(c.in)
			if c.wantErr {
				if err == nil {
					t.Fatalf("expected error for %q, got %+v", c.in, got)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if got != c.want {
				t.Fatalf("got %+v, want %+v", got, c.want)
			}
		})
	}
}

func TestFileFormatString(t *testing.T) {
	t.Parallel()
	if got := (FileFormat{4, 5}).String(); got != "VCFv4.5" {
		t.Errorf("got %q, want VCFv4.5", got)
	}
	if got := (FileFormat{Major: 0, Minor: 0}).String(); got != "VCFv0.0" {
		t.Errorf("got %q, want VCFv0.0", got)
	}
}

func TestHeaderLookups(t *testing.T) {
	t.Parallel()
	h := &Header{
		Infos: map[string]*InfoDef{
			"DP": {ID: "DP", Number: Number{Number: 1}, Type: TypeInteger},
		},
		Formats: map[string]*FormatDef{
			"GT": {ID: "GT", Number: Number{Number: 1}, Type: TypeString},
		},
		Filters: map[string]*FilterDef{
			"PASS": {ID: "PASS", Description: "All filters passed"},
		},
		Alts: map[string]*AltDef{
			"DEL": {ID: "DEL", Description: "Deletion"},
		},
		Contigs: []*Contig{
			{ID: "chr1", Length: 11},
		},
	}

	if d, ok := h.LookupInfo("DP"); !ok || d.ID != "DP" {
		t.Fatalf("LookupInfo(DP) = %v, %v", d, ok)
	}
	if d, ok := h.LookupInfo("ZZ"); ok {
		t.Fatalf("LookupInfo(ZZ) returned %v; want not found", d)
	}
	if d, ok := h.LookupFormat("GT"); !ok || d.ID != "GT" {
		t.Fatalf("LookupFormat(GT) = %v, %v", d, ok)
	}
	if d, ok := h.LookupFilter("PASS"); !ok || d.ID != "PASS" {
		t.Fatalf("LookupFilter(PASS) = %v, %v", d, ok)
	}
	if d, ok := h.LookupAlt("DEL"); !ok || d.ID != "DEL" {
		t.Fatalf("LookupAlt(DEL) = %v, %v", d, ok)
	}
	if c, ok := h.LookupContig("chr1"); !ok || c.Length != 11 {
		t.Fatalf("LookupContig(chr1) = %v, %v", c, ok)
	}

	// Nil header safety.
	var nh *Header
	if d, ok := nh.LookupInfo("DP"); ok {
		t.Fatalf("nil header LookupInfo returned %v", d)
	}
}
