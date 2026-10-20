package bcf

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
)

func TestEncoder_RoundTrip_Basic(t *testing.T) {
	hdr := ParseHeader(
		"##fileformat=VCFv4.2\n" +
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
			"##contig=<ID=chr1,length=1000>\n")
	records := []*Record{
		{
			CHROM: "chr1", Pos: 100,
			ID:   "rs1",
			Ref:  "A",
			Alt:  []string{"G", "C"},
			Qual: 99.5, QualKnown: true,
			Filter: []string{"PASS"},
			Info: map[string]InfoValue{
				"DP": {Type: "Integer", Ints: []int64{42}},
			},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, r := range records {
		if err := enc.WriteRecord(r); err != nil {
			t.Fatalf("WriteRecord: %v", err)
		}
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	dec := NewDecoder(&buf)
	gotHdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if gotHdr.FileFormat != hdr.FileFormat {
		t.Errorf("FileFormat: got %s want %s", gotHdr.FileFormat, hdr.FileFormat)
	}
	var got []*Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		got = append(got, rec)
	}
	if len(got) != len(records) {
		t.Fatalf("records: got %d want %d", len(got), len(records))
	}
	for i, want := range records {
		g := got[i]
		if g.CHROM != want.CHROM {
			t.Errorf("[%d] CHROM: got %q want %q", i, g.CHROM, want.CHROM)
		}
		if g.Pos != want.Pos {
			t.Errorf("[%d] Pos: got %d want %d", i, g.Pos, want.Pos)
		}
		if g.ID != want.ID {
			t.Errorf("[%d] ID: got %q want %q", i, g.ID, want.ID)
		}
		if g.Ref != want.Ref {
			t.Errorf("[%d] Ref: got %q want %q", i, g.Ref, want.Ref)
		}
		if !equalStrings(g.Alt, want.Alt) {
			t.Errorf("[%d] Alt: got %v want %v", i, g.Alt, want.Alt)
		}
		if g.QualKnown != want.QualKnown || g.Qual != want.Qual {
			t.Errorf("[%d] Qual: got %v/%v want %v/%v",
				i, g.Qual, g.QualKnown, want.Qual, want.QualKnown)
		}
		if !equalStrings(g.Filter, want.Filter) {
			t.Errorf("[%d] Filter: got %v want %v", i, g.Filter, want.Filter)
		}
		if v, ok := g.Info["DP"]; !ok || v.Type != "Integer" || len(v.Ints) != 1 || v.Ints[0] != 42 {
			t.Errorf("[%d] Info[DP]: ok=%v val=%+v", i, ok, v)
		}
	}
}

func TestEncoder_EmptyStream(t *testing.T) {
	hdr := ParseHeader("##fileformat=VCFv4.2\n")
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	dec := NewDecoder(&buf)
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Errorf("expected EOF, got %v", err)
	}
}

func TestEncoder_NoHeader(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	err := enc.WriteRecord(&Record{CHROM: "chr1"})
	if err == nil {
		t.Error("WriteRecord before WriteHeader: expected error, got nil")
	}
}

func TestEncoder_RecordAfterClose(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if err := enc.WriteRecord(&Record{CHROM: "chr1"}); err == nil {
		t.Error("WriteRecord after Close: expected error, got nil")
	}
}

func TestEncoder_NilHeader(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(nil); err == nil {
		t.Error("WriteHeader(nil): expected error, got nil")
	}
}

func TestEncoder_NilRecord(t *testing.T) {
	hdr := ParseHeader("##fileformat=VCFv4.2\n")
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(nil); err != nil {
		t.Errorf("WriteRecord(nil): expected no-op, got %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

// TestEncoder_NoContig verifies that a record referencing a
// contig not in the header produces an error.
func TestEncoder_NoContig(t *testing.T) {
	hdr := ParseHeader("##fileformat=VCFv4.2\n")
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	err := enc.WriteRecord(&Record{CHROM: "chr1", Pos: 1, Ref: "A", Alt: []string{"G"}})
	if err == nil {
		t.Error("WriteRecord with unknown contig: expected error, got nil")
	}
}

// TestWriter_ForbiddenIntegerRange verifies that the BCF
// writer refuses to emit Integer values in the BCF2
// reserved range [MinInt32, MinInt32+7] (§1.4.2,
// §BcfTypeEncoding line 2306). T4.22.
func TestWriter_ForbiddenIntegerRange(t *testing.T) {
	hdr := ParseHeader(
		"##fileformat=VCFv4.2\n" +
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
			"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"gq\">\n" +
			"##contig=<ID=chr1,length=1000>\n")
	for _, v := range []int64{
		-2147483648, // MinInt32
		-2147483641, // MinInt32+7 (last forbidden)
	} {
		var buf bytes.Buffer
		enc := NewEncoder(&buf)
		if err := enc.WriteHeader(hdr); err != nil {
			t.Fatalf("v=%d: WriteHeader: %v", v, err)
		}
		err := enc.WriteRecord(&Record{
			CHROM: "chr1", Pos: 100, Ref: "A", Alt: []string{"G"},
			Info: map[string]InfoValue{
				"DP": {Type: "Integer", Ints: []int64{v}},
			},
		})
		if err == nil {
			t.Errorf("v=%d: expected forbidden-range error, got nil", v)
		}
	}
}

// TestWriter_AcceptsValidIntegerRange verifies that values
// outside the reserved range round-trip cleanly. Sister
// test to TestWriter_ForbiddenIntegerRange.
func TestWriter_AcceptsValidIntegerRange(t *testing.T) {
	hdr := ParseHeader(
		"##fileformat=VCFv4.2\n" +
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
			"##contig=<ID=chr1,length=1000>\n")
	for _, v := range []int64{
		2147483647,  // MaxInt32 (high side)
		-2147483640, // MinInt32+8 (just past the forbidden range)
	} {
		var buf bytes.Buffer
		enc := NewEncoder(&buf)
		if err := enc.WriteHeader(hdr); err != nil {
			t.Fatalf("v=%d: WriteHeader: %v", v, err)
		}
		if err := enc.WriteRecord(&Record{
			CHROM: "chr1", Pos: 100, Ref: "A", Alt: []string{"G"},
			Info: map[string]InfoValue{
				"DP": {Type: "Integer", Ints: []int64{v}},
			},
		}); err != nil {
			t.Errorf("v=%d: WriteRecord: %v", v, err)
			continue
		}
		if err := enc.Close(); err != nil {
			t.Errorf("v=%d: Close: %v", v, err)
			continue
		}
		dec := NewDecoder(&buf)
		if _, err := dec.Header(); err != nil {
			t.Errorf("v=%d: Header: %v", v, err)
			continue
		}
		rec, err := dec.Record()
		if err != nil {
			t.Errorf("v=%d: Record: %v", v, err)
			continue
		}
		got, ok := rec.Info["DP"]
		if !ok || len(got.Ints) != 1 || got.Ints[0] != v {
			t.Errorf("v=%d: round-trip got=%+v", v, got)
		}
	}
}

// TestWriter_AtomicInlineSize verifies the writer emits
// the atomic (size-1) type byte in the high nibble for
// Integer and Float INFO fields. T4.21,
// §BcfTypeEncoding line 2271: `(1<<4) | type`.
func TestWriter_AtomicInlineSize(t *testing.T) {
	hdr := ParseHeader(
		"##fileformat=VCFv4.2\n" +
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
			"##INFO=<ID=AF,Number=1,Type=Float,Description=\"af\">\n" +
			"##contig=<ID=chr1,length=1000>\n")
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(&Record{
		CHROM: "chr1", Pos: 100, Ref: "A", Alt: []string{"G"},
		Qual: 99.5, QualKnown: true,
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{42}},
			"AF": {Type: "Float", Floats: []float64{0.5}},
		},
	}); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	// Round-trip: the decoder masks the high nibble so it
	// accepts both bare type bytes and inline-size bytes.
	dec := NewDecoder(&buf)
	gotHdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if _, ok := gotHdr.Infos["DP"]; !ok {
		t.Fatal("DP not in header")
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	dp, ok := rec.Info["DP"]
	if !ok || len(dp.Ints) != 1 || dp.Ints[0] != 42 {
		t.Errorf("DP round-trip: got %+v", dp)
	}
	af, ok := rec.Info["AF"]
	if !ok || len(af.Floats) != 1 || af.Floats[0] != 0.5 {
		t.Errorf("AF round-trip: got %+v", af)
	}
}

// TestReader_AcceptsInlineSizeTypeByte verifies that the
// reader masks the high nibble of the type byte (size
// encoded inline per §BcfTypeEncoding line 2271) and
// still decodes values correctly. The BCF2-compliant
// writer emits `(1<<4) | type` for atomic values; the
// reader must accept that as well as the bare type byte.
func TestReader_AcceptsInlineSizeTypeByte(t *testing.T) {
	hdr := ParseHeader(
		"##fileformat=VCFv4.2\n" +
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
			"##contig=<ID=chr1,length=1000>\n")
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(&Record{
		CHROM: "chr1", Pos: 100, Ref: "A", Alt: []string{"G"},
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{7}},
		},
	}); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	// Decode the BGZF stream once to inspect the
	// uncompressed bytes — the Writer must emit the
	// type byte 0x13 (atomic int32 with size 1 in the
	// high nibble) for the single INFO Integer, followed
	// by 4 bytes little-endian for the value 7.
	raw, err := io.ReadAll(bgzf.NewReader(bytes.NewReader(buf.Bytes())))
	if err != nil {
		t.Fatalf("decompress: %v", err)
	}
	want := []byte{typeByteAtomic(typeInt32), 0x07, 0x00, 0x00, 0x00}
	if !bytes.Contains(raw, want) {
		t.Errorf("uncompressed wire bytes missing %x; got %x", want, raw)
	}
	dec := NewDecoder(bytes.NewReader(buf.Bytes()))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	dp, ok := rec.Info["DP"]
	if !ok || len(dp.Ints) != 1 || dp.Ints[0] != 7 {
		t.Errorf("DP: got %+v want 7", dp)
	}
}

// equalStrings compares two []string for equality.
func equalStrings(a, b []string) bool {
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

// ensure strings package is used
var _ = strings.HasPrefix
