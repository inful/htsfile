package bed

import (
	"errors"
	"io"
	"os"
	"strings"
	"testing"
)

// readAllRecords drains the decoder for testing.
func readAllRecords(t *testing.T, dec *Decoder) (*Track, []*Record) {
	t.Helper()
	track, err := dec.Track()
	if err != nil {
		t.Fatalf("Track: %v", err)
	}
	var recs []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		recs = append(recs, rec)
	}
	return track, recs
}

func openTestFile(t *testing.T, path string) (*os.File, func()) {
	t.Helper()
	//nolint:gosec // path comes from test code, not user input
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	return f, func() { _ = f.Close() }
}

func TestDecoder_BED3_NoTrack(t *testing.T) {
	t.Parallel()
	f, cleanup := openTestFile(t, "testdata/bed3.bed")
	defer cleanup()
	track, recs := readAllRecords(t, NewDecoder(f))
	if track != nil {
		t.Errorf("Track = %+v; want nil", track)
	}
	if len(recs) != 3 {
		t.Fatalf("len(recs) = %d; want 3", len(recs))
	}
	if recs[0].Chrom != "chr1" || recs[0].ChromStart != 100 || recs[0].ChromEnd != 200 {
		t.Errorf("recs[0] = %+v", recs[0])
	}
	if recs[0].ColumnCount() != 3 {
		t.Errorf("recs[0].ColumnCount = %d; want 3", recs[0].ColumnCount())
	}
	if recs[0].Name != nil {
		t.Error("BED3 record has Name; want nil")
	}
}

func TestDecoder_BED6_WithTrack(t *testing.T) {
	t.Parallel()
	f, cleanup := openTestFile(t, "testdata/bed6.bed")
	defer cleanup()
	track, recs := readAllRecords(t, NewDecoder(f))
	if track == nil {
		t.Fatal("Track = nil; want populated")
	}
	if track.Name != "peaks" {
		t.Errorf("Track.Name = %q; want peaks", track.Name)
	}
	if track.Description != "MACS2 narrow peaks" {
		t.Errorf("Track.Description = %q", track.Description)
	}
	if !track.UseScore {
		t.Error("Track.UseScore = false; want true")
	}
	if len(recs) != 3 {
		t.Fatalf("len(recs) = %d; want 3", len(recs))
	}
	if recs[0].ColumnCount() != 6 {
		t.Errorf("recs[0].ColumnCount = %d; want 6", recs[0].ColumnCount())
	}
	if recs[0].Name == nil || *recs[0].Name != "peak_1" {
		t.Errorf("recs[0].Name = %v; want peak_1", recs[0].Name)
	}
	if recs[0].Score == nil || *recs[0].Score != 500 {
		t.Errorf("recs[0].Score = %v; want 500", recs[0].Score)
	}
	if recs[0].Strand == nil || *recs[0].Strand != StrandPlus {
		t.Errorf("recs[0].Strand = %v; want StrandPlus", recs[0].Strand)
	}
}

func TestDecoder_BED12_FullBlocks(t *testing.T) {
	t.Parallel()
	f, cleanup := openTestFile(t, "testdata/bed12.bed")
	defer cleanup()
	track, recs := readAllRecords(t, NewDecoder(f))
	if track == nil {
		t.Fatal("Track = nil; want populated")
	}
	if !track.ItemRgbOn {
		t.Error("Track.ItemRgbOn = false; want true (from itemRgb=\"On\")")
	}
	if len(recs) != 2 {
		t.Fatalf("len(recs) = %d; want 2", len(recs))
	}
	r0 := recs[0]
	if r0.ColumnCount() != 12 {
		t.Errorf("r0.ColumnCount = %d; want 12", r0.ColumnCount())
	}
	if r0.BlockCount == nil || *r0.BlockCount != 3 {
		t.Errorf("r0.BlockCount = %v; want 3", r0.BlockCount)
	}
	wantSizes := []int64{100, 200, 300}
	if !sliceEq(r0.BlockSizes, wantSizes) {
		t.Errorf("r0.BlockSizes = %v; want %v", r0.BlockSizes, wantSizes)
	}
	wantStarts := []int64{0, 1500, 3700}
	if !sliceEq(r0.BlockStarts, wantStarts) {
		t.Errorf("r0.BlockStarts = %v; want %v", r0.BlockStarts, wantStarts)
	}
	if r0.ItemRgb == nil || *r0.ItemRgb != (RGB{0, 0, 255}) {
		t.Errorf("r0.ItemRgb = %v; want {0,0,255}", r0.ItemRgb)
	}
}

func TestDecoder_CommentsAndBlankLines(t *testing.T) {
	t.Parallel()
	f, cleanup := openTestFile(t, "testdata/comments.bed")
	defer cleanup()
	track, recs := readAllRecords(t, NewDecoder(f))
	if track == nil || track.Name != "with-comments" {
		t.Errorf("Track = %+v; want Name=with-comments", track)
	}
	if len(recs) != 2 {
		t.Errorf("len(recs) = %d; want 2 (comments and blanks skipped)", len(recs))
	}
}

func TestDecoder_Malformed(t *testing.T) {
	t.Parallel()
	f, cleanup := openTestFile(t, "testdata/malformed.bed")
	defer cleanup()
	dec := NewDecoder(f)
	if _, err := dec.Track(); err != nil {
		t.Fatalf("Track: %v", err)
	}
	// First record (chr1 100 200 peak_1 500 +) parses.
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("first Record: %v", err)
	}
	if rec.Chrom != "chr1" {
		t.Errorf("rec.Chrom = %q; want chr1", rec.Chrom)
	}
	// Second record has a non-numeric chromStart → error.
	_, err = dec.Record()
	if err == nil {
		t.Fatal("second Record: expected error, got nil")
	}
	if !IsParseError(err) {
		t.Errorf("error type = %T; want *ParseError", err)
	}
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatal("errors.As failed")
	}
	if pe.Field != "chromStart" {
		t.Errorf("pe.Field = %q; want chromStart", pe.Field)
	}
}

func TestDecoder_TooFewFields(t *testing.T) {
	t.Parallel()
	src := "chr1\t100\nchr1\t200\t300\n"
	dec := NewDecoder(strings.NewReader(src))
	// No track.
	_, _ = dec.Track()
	_, err := dec.Record()
	if err == nil {
		t.Fatal("expected error for line with < 3 fields")
	}
	if !IsParseError(err) {
		t.Errorf("error type = %T; want *ParseError", err)
	}
}

func TestDecoder_EOFOnEmpty(t *testing.T) {
	t.Parallel()
	dec := NewDecoder(strings.NewReader(""))
	track, err := dec.Track()
	if err != nil {
		t.Errorf("Track err = %v; want nil", err)
	}
	if track != nil {
		t.Errorf("Track = %+v; want nil on empty", track)
	}
	rec, err := dec.Record()
	if !errors.Is(err, io.EOF) {
		t.Errorf("Record err = %v; want io.EOF", err)
	}
	if rec != nil {
		t.Errorf("Record = %+v; want nil on EOF", rec)
	}
	// Subsequent calls also return io.EOF.
	rec, err = dec.Record()
	if !errors.Is(err, io.EOF) {
		t.Errorf("second Record err = %v; want io.EOF", err)
	}
	if rec != nil {
		t.Errorf("second Record = %+v; want nil", rec)
	}
}

func TestDecoder_AutoCallsTrackOnFirstRecord(t *testing.T) {
	t.Parallel()
	src := "track name=auto\nchr1\t100\t200\n"
	dec := NewDecoder(strings.NewReader(src))
	// No explicit Track() call.
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if rec.Chrom != "chr1" {
		t.Errorf("rec.Chrom = %q; want chr1", rec.Chrom)
	}
	track, err := dec.Track()
	if err != nil {
		t.Errorf("Track: %v", err)
	}
	if track == nil || track.Name != "auto" {
		t.Errorf("Track = %+v; want Name=auto", track)
	}
}

func TestParseRGB(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in      string
		want    RGB
		wantErr bool
	}{
		{"0,0,0", RGB{}, false},
		{"255,128,64", RGB{255, 128, 64}, false},
		{"  10 , 20 , 30  ", RGB{10, 20, 30}, false},
		{"0,0", RGB{}, true},     // too few components
		{"0,0,0,0", RGB{}, true}, // too many
		{"256,0,0", RGB{}, true}, // out of range
		{"abc,0,0", RGB{}, true}, // non-numeric
	}
	for _, c := range cases {
		got, err := parseRGB(c.in)
		if c.wantErr {
			if err == nil {
				t.Errorf("parseRGB(%q) = nil err; want non-nil", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("parseRGB(%q) unexpected err: %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("parseRGB(%q) = %+v; want %+v", c.in, got, c.want)
		}
	}
}

func TestParseTrack_UnterminatedQuote(t *testing.T) {
	t.Parallel()
	raw := []byte(`track name="unterminated description=foo`)
	_, err := parseTrack(raw)
	if err == nil {
		t.Fatal("expected error for unterminated quoted value")
	}
	if !strings.Contains(err.Error(), "unterminated") {
		t.Errorf("err = %v; want it to mention 'unterminated'", err)
	}
}

func TestParseIntList(t *testing.T) {
	t.Parallel()
	got, err := parseIntList("1,2,3,4")
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if !sliceEq(got, []int64{1, 2, 3, 4}) {
		t.Errorf("got = %v; want [1 2 3 4]", got)
	}
	if _, err := parseIntList(""); err == nil {
		t.Error("expected error on empty list")
	}
	if _, err := parseIntList("1,x,3"); err == nil {
		t.Error("expected error on non-numeric element")
	}
}

func sliceEq(a, b []int64) bool {
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

// TestDecoder_MaxLineLength pins the bounded-read behavior: a
// line longer than MaxLineLength must error, not silently allocate
// gigabytes.
func TestDecoder_MaxLineLength(t *testing.T) {
	t.Parallel()
	// 200 'A' characters in chromStart, with a 100-byte limit.
	longStart := strings.Repeat("9", 200)
	src := "track name=\"oops\"\n" +
		"chr1\t" + longStart + "\t300\n"

	dec := NewDecoderWithMaxLine(strings.NewReader(src), 100)
	if _, err := dec.Track(); err != nil {
		t.Fatalf("Track: %v", err)
	}
	_, err := dec.Record()
	if err == nil {
		t.Fatal("expected error for line exceeding MaxLineLength")
	}
	if !IsParseError(err) {
		t.Errorf("error type = %T; want *ParseError", err)
	}
	if !strings.Contains(err.Error(), "exceeds") {
		t.Errorf("error msg = %q; want it to mention 'exceeds'", err.Error())
	}
}
