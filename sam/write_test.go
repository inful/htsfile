package sam

import (
	"bytes"
	"errors"
	"io"
	"os"
	"reflect"
	"strings"
	"testing"
)

func TestEncoder_RecordOnly(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	rec := &Record{
		QName: "r1", Flag: 0, RName: "chr1",
		Pos: 100, MapQ: 60, Cigar: "100M",
		RNext: "*", PNext: 0, TLen: 0,
		Seq:  "ACGTACGTAC",
		Qual: "IIIIIIIIII",
		Tags: []Tag{
			{Key: "NM", Type: Tagi, Value: "1"},
			{Key: "MD", Type: TagZ, Value: "99"},
		},
	}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "r1\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:1\tMD:Z:99\n"
	if got := buf.String(); got != want {
		t.Errorf("output = %q; want %q", got, want)
	}
}

func TestEncoder_CRLF(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoderCRLF(&buf)
	rec := &Record{QName: "r1", Seq: "ACGT", Qual: "IIII"}
	if err := enc.EncodeRecord(rec); err != nil {
		t.Fatalf("EncodeRecord: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "r1\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\r\n"
	if got := buf.String(); got != want {
		t.Errorf("output = %q; want %q", got, want)
	}
}

func TestEncoder_Header_RawVerbatim(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	hdr := &Header{
		HD: &HD{VN: "1.6", SO: "coordinate", Raw: "@HD\tVN:1.6\tSO:coordinate\tGO:none"},
	}
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "@HD\tVN:1.6\tSO:coordinate\tGO:none\n"
	if got := buf.String(); got != want {
		t.Errorf("output = %q; want %q", got, want)
	}
}

func TestEncoder_Header_Synthesised(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	hdr := &Header{
		HD: &HD{VN: "1.6", SO: "coordinate"},
		SQs: []*SQ{
			{SN: "chr1", LN: 1000},
			{SN: "chr2", LN: 2000, M5: "abc"},
		},
		RGs: []*RG{
			{ID: "RG1", SM: "sample1", LB: "lib1", PL: "ILLUMINA"},
		},
		PGs: []*PG{
			{ID: "bwa", PN: "bwa", VN: "0.7.17"},
		},
		COs: []string{"@CO\thello"},
	}
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	// Empty fields are omitted (e.g. HD.GO is empty so not emitted).
	wantLines := []string{
		"@HD\tSO:coordinate\tVN:1.6",
		"@SQ\tLN:1000\tSN:chr1",
		"@SQ\tLN:2000\tM5:abc\tSN:chr2",
		"@RG\tID:RG1\tLB:lib1\tPL:ILLUMINA\tSM:sample1",
		"@PG\tID:bwa\tPN:bwa\tVN:0.7.17",
		"@CO\thello",
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

func TestEncoder_NilHeader(t *testing.T) {
	t.Parallel()
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(nil); err != nil {
		t.Errorf("EncodeHeader(nil) = %v; want nil", err)
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if buf.Len() != 0 {
		t.Errorf("output = %q; want empty", buf.String())
	}
}

func TestRoundTripCorpus(t *testing.T) {
	t.Parallel()
	// Parse the testdata file, re-emit, re-parse, and compare
	// field-by-field. This is the strongest end-to-end test.
	original, err := os.Open("testdata/example.sam")
	if err != nil {
		t.Fatal(err)
	}
	defer original.Close() //nolint:errcheck

	dec := NewDecoder(original)
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
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

	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		t.Fatalf("EncodeHeader: %v", err)
	}
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			t.Fatalf("EncodeRecord: %v", err)
		}
	}
	if err := enc.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}

	// Re-parse and compare. Use a fresh bytes.Reader so the
	// second decode doesn't share state with the first.
	rt, err := NewDecoder(bytes.NewReader(buf.Bytes())).Header()
	if err != nil {
		t.Fatalf("re-parse Header: %v", err)
	}
	if !headerEqual(hdr, rt) {
		t.Errorf("header round-trip differs:\n orig: %+v\n rt:   %+v", hdr, rt)
	}
	var rtRecs []*Record
	dec2 := NewDecoder(bytes.NewReader(buf.Bytes()))
	if _, err := dec2.Header(); err != nil {
		t.Fatal(err)
	}
	for {
		rec, err := dec2.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("rt Record: %v", err)
		}
		rtRecs = append(rtRecs, rec)
	}
	if len(recs) != len(rtRecs) {
		t.Fatalf("record count: orig=%d rt=%d", len(recs), len(rtRecs))
	}
	for i, a := range recs {
		b := rtRecs[i]
		if !recordEqual(a, b) {
			t.Errorf("record %d differs:\n orig: %+v\n rt:   %+v", i, a, b)
		}
	}
}

func headerEqual(a, b *Header) bool {
	if a == nil || b == nil {
		return a == b
	}
	if (a.HD == nil) != (b.HD == nil) {
		return false
	}
	if a.HD != nil && !hdEqual(a.HD, b.HD) {
		return false
	}
	if len(a.SQs) != len(b.SQs) {
		return false
	}
	for i := range a.SQs {
		if !sqEqual(a.SQs[i], b.SQs[i]) {
			return false
		}
	}
	if len(a.RGs) != len(b.RGs) {
		return false
	}
	for i := range a.RGs {
		if !rgEqual(a.RGs[i], b.RGs[i]) {
			return false
		}
	}
	if len(a.PGs) != len(b.PGs) {
		return false
	}
	for i := range a.PGs {
		if !pgEqual(a.PGs[i], b.PGs[i]) {
			return false
		}
	}
	return true
}

func hdEqual(a, b *HD) bool {
	return a.VN == b.VN && a.SO == b.SO && a.GO == b.GO &&
		a.Raw == b.Raw && reflect.DeepEqual(a.Other, b.Other)
}

func sqEqual(a, b *SQ) bool {
	return a.SN == b.SN && a.LN == b.LN && a.AH == b.AH && a.AN == b.AN &&
		a.AS == b.AS && a.M5 == b.M5 && a.UR == b.UR && a.SP == b.SP &&
		a.Raw == b.Raw && reflect.DeepEqual(a.Other, b.Other)
}

func rgEqual(a, b *RG) bool {
	return a.ID == b.ID && a.SM == b.SM && a.LB == b.LB && a.DS == b.DS &&
		a.PU == b.PU && a.PI == b.PI && a.CN == b.CN && a.DT == b.DT &&
		a.PL == b.PL && a.PG == b.PG &&
		a.Raw == b.Raw && reflect.DeepEqual(a.Other, b.Other)
}

func pgEqual(a, b *PG) bool {
	return a.ID == b.ID && a.PN == b.PN && a.CL == b.CL && a.PP == b.PP &&
		a.VN == b.VN && a.Raw == b.Raw && reflect.DeepEqual(a.Other, b.Other)
}

func recordEqual(a, b *Record) bool {
	if a == nil || b == nil {
		return a == b
	}
	if a.QName != b.QName || a.Flag != b.Flag || a.RName != b.RName ||
		a.Pos != b.Pos || a.MapQ != b.MapQ || a.Cigar != b.Cigar ||
		a.RNext != b.RNext || a.PNext != b.PNext || a.TLen != b.TLen ||
		a.Seq != b.Seq || a.Qual != b.Qual {
		return false
	}
	if len(a.Tags) != len(b.Tags) {
		return false
	}
	for i := range a.Tags {
		if a.Tags[i] != b.Tags[i] {
			return false
		}
	}
	return true
}
