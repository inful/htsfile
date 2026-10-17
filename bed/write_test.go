package bed

import (
	"bytes"
	"os"
	"strings"
	"testing"
)

func TestEncoder_BED3_NoTrack(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	recs := []*Record{
		{Chrom: "chr1", ChromStart: 100, ChromEnd: 200},
		{Chrom: "chr2", ChromStart: 500, ChromEnd: 600},
	}
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatalf("EncodeRecord: %v", err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "chr1\t100\t200\nchr2\t500\t600\n"
	if got := buf.String(); got != want {
		t.Errorf("output = %q; want %q", got, want)
	}
}

func TestEncoder_BED6_WithTrack(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	track := &Track{
		Name:        "peaks",
		Description: "MACS2 narrow peaks",
		UseScore:    true,
	}
	if err := enc.EncodeTrack(track); err != nil {
		t.Fatalf("EncodeTrack: %v", err)
	}
	plus := StrandPlus
	minus := StrandMinus
	none := StrandNone
	recs := []*Record{
		{Chrom: "chr1", ChromStart: 100, ChromEnd: 200,
			Name: strPtr("peak_1"), Score: intPtr(500), Strand: &plus},
		{Chrom: "chr1", ChromStart: 300, ChromEnd: 400,
			Name: strPtr("peak_2"), Score: intPtr(900), Strand: &minus},
		{Chrom: "chr2", ChromStart: 500, ChromEnd: 600,
			Name: strPtr("peak_3"), Score: intPtr(750), Strand: &none},
	}
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatalf("EncodeRecord: %v", err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	wantLines := []string{
		`track name="peaks" description="MACS2 narrow peaks" useScore=1`,
		"chr1\t100\t200\tpeak_1\t500\t+",
		"chr1\t300\t400\tpeak_2\t900\t-",
		"chr2\t500\t600\tpeak_3\t750\t.",
	}
	got := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(got) != len(wantLines) {
		t.Fatalf("got %d lines, want %d:\n%s", len(got), len(wantLines), buf.String())
	}
	for i, w := range wantLines {
		if got[i] != w {
			t.Errorf("line %d = %q; want %q", i, got[i], w)
		}
	}
}

func TestEncoder_BED12_FullBlocks(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	plus := StrandPlus
	rec := &Record{
		Chrom:       "chr1",
		ChromStart:  1000,
		ChromEnd:    5000,
		Name:        strPtr("NM_001"),
		Score:       intPtr(960),
		Strand:      &plus,
		ThickStart:  int64Ptr(1200),
		ThickEnd:    int64Ptr(4800),
		ItemRgb:     &RGB{0, 0, 255},
		BlockCount:  intPtr(3),
		BlockSizes:  []int64{100, 200, 300},
		BlockStarts: []int64{0, 1500, 3700},
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "chr1\t1000\t5000\tNM_001\t960\t+\t1200\t4800\t0,0,255\t3\t100,200,300\t0,1500,3700\n"
	if got := buf.String(); got != want {
		t.Errorf("output = %q; want %q", got, want)
	}
}

func TestEncoder_RawTrackVerbatim(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	raw := `track name=foo description="an example" useScore=1 itemRgb="On" customAttr=42`
	track := &Track{Raw: raw}
	if err := enc.EncodeTrack(track); err != nil {
		t.Fatalf("EncodeTrack: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := raw + "\n"
	if got := buf.String(); got != want {
		t.Errorf("output = %q; want %q", got, want)
	}
}

func TestEncoder_CRLF(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoderCRLF(&buf)
	rec := &Record{Chrom: "chr1", ChromStart: 100, ChromEnd: 200}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "chr1\t100\t200\r\n"
	if got := buf.String(); got != want {
		t.Errorf("output = %q; want %q", got, want)
	}
}

func TestEncoder_NilTrackIsNoOp(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeTrack(nil); err != nil {
		t.Errorf("EncodeTrack(nil) = %v; want nil", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if buf.Len() != 0 {
		t.Errorf("output = %q; want empty", buf.String())
	}
}

// TestRoundTripCorpus parses each testdata file and re-encodes it,
// confirming the output is semantically equivalent. The textual
// representation may differ (e.g. trailing-field order on the track
// line) but the parsed data must match.
func TestRoundTripCorpus(t *testing.T) {
	t.Parallel()
	for _, name := range []string{"bed3.bed", "bed6.bed", "bed12.bed", "comments.bed"} {
		t.Run(name, func(t *testing.T) {
			t.Parallel()
			track, recs := readAllRecords(t, mustOpenDecoder(t, "testdata/"+name))

			var buf bytes.Buffer
			enc := NewEncoder(&buf)
			if track != nil {
				if err := enc.EncodeTrack(track); err != nil {
					t.Fatalf("EncodeTrack: %v", err)
				}
			}
			for _, r := range recs {
				if err := enc.EncodeRecord(r); err != nil {
					t.Fatalf("EncodeRecord: %v", err)
				}
			}
			if err := enc.Flush(); err != nil {
				t.Fatalf("Flush: %v", err)
			}

			// Re-parse the encoder output and compare field-by-field.
			track2, recs2 := readAllRecords(t, NewDecoder(&buf))
			if (track == nil) != (track2 == nil) {
				t.Errorf("track presence mismatch: original=%v roundtripped=%v", track, track2)
			}
			if track != nil && track2 != nil {
				// Use the raw line for comparison: round-trip is text
				// identical when the track line is preserved via Raw.
				if track.Raw != track2.Raw {
					t.Errorf("track Raw mismatch:\n  orig: %q\n  rt:   %q", track.Raw, track2.Raw)
				}
			}
			if len(recs) != len(recs2) {
				t.Fatalf("record count mismatch: orig=%d rt=%d", len(recs), len(recs2))
			}
			for i, a := range recs {
				b := recs2[i]
				if !recordEqual(a, b) {
					t.Errorf("record %d mismatch:\n  orig: %+v\n  rt:   %+v", i, a, b)
				}
			}
		})
	}
}

func mustOpenDecoder(t *testing.T, path string) *Decoder {
	t.Helper()
	//nolint:gosec // path comes from test code, not user input
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { _ = f.Close() })
	return NewDecoder(f)
}

func recordEqual(a, b *Record) bool {
	if a == nil || b == nil {
		return a == b
	}
	if a.Chrom != b.Chrom || a.ChromStart != b.ChromStart || a.ChromEnd != b.ChromEnd {
		return false
	}
	if !strPtrEq(a.Name, b.Name) || !intPtrEq(a.Score, b.Score) || !strandPtrEq(a.Strand, b.Strand) {
		return false
	}
	if !int64PtrEq(a.ThickStart, b.ThickStart) || !int64PtrEq(a.ThickEnd, b.ThickEnd) {
		return false
	}
	if !rgbPtrEq(a.ItemRgb, b.ItemRgb) || !intPtrEq(a.BlockCount, b.BlockCount) {
		return false
	}
	if !sliceEq(a.BlockSizes, b.BlockSizes) || !sliceEq(a.BlockStarts, b.BlockStarts) {
		return false
	}
	return true
}

func strPtrEq(a, b *string) bool {
	if a == nil || b == nil {
		return a == b
	}
	return *a == *b
}

func intPtrEq(a, b *int) bool {
	if a == nil || b == nil {
		return a == b
	}
	return *a == *b
}

func int64PtrEq(a, b *int64) bool {
	if a == nil || b == nil {
		return a == b
	}
	return *a == *b
}

func strandPtrEq(a, b *Strand) bool {
	if a == nil || b == nil {
		return a == b
	}
	return *a == *b
}

func rgbPtrEq(a, b *RGB) bool {
	if a == nil || b == nil {
		return a == b
	}
	return *a == *b
}
