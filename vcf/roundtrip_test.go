package vcf

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

// TestReview_DroppedTrailingFieldsRoundTrip checks that a record with a
// trailing dropped FORMAT field (legal per spec §1.6.1) round-trips without
// re-emitting the missing field as "." — i.e. the encoded bytes match the
// input bytes for that column.
func TestReview_DroppedTrailingFieldsRoundTrip(t *testing.T) {
	t.Parallel()
	src := "##fileformat=VCFv4.5\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"h\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" +
		"1\t2\t.\tA\tC\t.\t.\t.\tGT:DP:HQ\t0/1:42,51,51\n" +
		"1\t3\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:42\n"

	// Decode, encode, decode — the second line has 3 fields, FORMAT has 2.
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatal(err)
	}
	var recs []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		recs = append(recs, rec)
	}

	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatal(err)
	}

	out := buf.String()
	t.Logf("encoded output:\n%s", out)

	// The first record has HQ dropped from the sample column. The
	// encoder must preserve the drop — emit the FORMAT keys up to the
	// last present one (HQ here) but with no trailing ".", matching the
	// input bytes.
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 data lines, got %d:\n%s", len(lines), out)
	}
	first := lines[0]
	cols := strings.Split(first, "\t")
	if len(cols) < 10 {
		t.Fatalf("expected 10 columns in first record, got %d:\n%s",
			len(cols), first)
	}
	format := cols[8]
	sample := cols[9]
	if format != "GT:DP:HQ" {
		t.Errorf("FORMAT = %q; want GT:DP:HQ", format)
	}
	// Critical assertion: the sample must NOT have a trailing '.' for
	// the dropped HQ. Original input was "0/1:42,51,51" with HQ
	// implicitly dropped; output must match.
	if sample != "0/1:42,51,51" {
		t.Errorf("sample = %q; want 0/1:42,51,51 (no trailing '.')", sample)
	}

	// Second record — GT:DP present, all present.
	second := lines[1]
	cols2 := strings.Split(second, "\t")
	if cols2[9] != "0/1:42" {
		t.Errorf("second sample = %q; want 0/1:42", cols2[9])
	}
}
