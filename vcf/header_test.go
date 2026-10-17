package vcf

import (
	"bytes"
	"fmt"
	"strings"
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
		{"VCFv4.3", FileFormat{4, 3}, false},
		{"VCFv4.2", FileFormat{4, 2}, false},
		{"VCFv4.1", FileFormat{4, 1}, false},
		{"VCFv4.0", FileFormat{4, 0}, false},    // accepted: any v4.x minor
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

// TestRoundTripFileFormat pins the v4.1-v4.5 acceptance guarantee
// at the full document level: a header that starts with
// "##fileformat=VCFvX.Y" round-trips with the same X.Y in the
// emitted line. This is the user-facing promise of the v4.1-v4.5
// support story; per-version validation differences are not yet
// implemented (they would need a v4.1/v4.2/v4.3/v4.4 spec read to
// enumerate). The library is "liberal in what it accepts,
// conservative in what it emits" — the FileFormat is preserved
// verbatim on round-trip regardless of which v4.x it is.
func TestRoundTripFileFormat(t *testing.T) {
	t.Parallel()
	for _, minor := range []int{1, 2, 3, 4, 5} {
		t.Run(fmt.Sprintf("v4.%d", minor), func(t *testing.T) {
			t.Parallel()
			src := fmt.Sprintf("##fileformat=VCFv4.%d\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", minor) +
				"1\t2\t.\tA\tC\t.\t.\t.\n"
			dec := NewDecoder(strings.NewReader(src))
			hdr, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			if hdr.FileFormat.Major != 4 || hdr.FileFormat.Minor != minor {
				t.Errorf("FileFormat = %d.%d; want 4.%d",
					hdr.FileFormat.Major, hdr.FileFormat.Minor, minor)
			}
			var buf bytes.Buffer
			enc := NewEncoder(&buf)
			if err := enc.EncodeHeader(hdr); err != nil {
				t.Fatalf("EncodeHeader: %v", err)
			}
			if err := enc.Flush(); err != nil {
				t.Fatalf("Flush: %v", err)
			}
			wantPrefix := fmt.Sprintf("##fileformat=VCFv4.%d", minor)
			if !strings.HasPrefix(buf.String(), wantPrefix) {
				t.Errorf("output = %q; want prefix %q", buf.String(), wantPrefix)
			}
		})
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
