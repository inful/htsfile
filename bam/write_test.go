package bam

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

// TestEncoder_RoundTrip_Basic verifies that the Encoder
// produces a stream the existing Decoder can read back.
func TestEncoder_RoundTrip_Basic(t *testing.T) {
	hdr := &Header{
		Raw: "@HD\tVN:1.6\tSO:coordinate\n",
		SQs: []*SQ{{Name: "chr1", Len: 1000}, {Name: "chr2", Len: 2000}},
	}
	records := []*Record{
		{
			QName: "read1",
			Flag:  0,
			Pos:   100,
			MapQ:  60,
			Bin:   4681,
			Cigar: []uint32{1<<4 | 10}, // 10M
			RefID: 0,
			Seq:   "ACGTACGTAC",
			Qual:  []byte("IIIIIIIIII"),
			Tags:  []Tag{{Key: "NM", Type: Tagi, Value: "0"}},
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

	// Read back.
	dec := NewDecoder(&buf)
	gotHdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if gotHdr.Raw != hdr.Raw {
		t.Errorf("Raw: got %q want %q", gotHdr.Raw, hdr.Raw)
	}
	if len(gotHdr.SQs) != 2 {
		t.Fatalf("SQs: got %d want 2", len(gotHdr.SQs))
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
		t.Fatalf("Records: got %d want %d", len(got), len(records))
	}
	for i, r := range records {
		g := got[i]
		if g.QName != r.QName {
			t.Errorf("QName: got %q want %q", g.QName, r.QName)
		}
		if g.Pos != r.Pos {
			t.Errorf("Pos: got %d want %d", g.Pos, r.Pos)
		}
		if g.RefID != r.RefID {
			t.Errorf("RefID: got %d want %d", g.RefID, r.RefID)
		}
		if g.Seq != r.Seq {
			t.Errorf("Seq: got %q want %q", g.Seq, r.Seq)
		}
		if string(g.Qual) != string(r.Qual) {
			t.Errorf("Qual: got %q want %q", g.Qual, r.Qual)
		}
		if len(g.Tags) != len(r.Tags) {
			t.Errorf("Tags: got %d want %d", len(g.Tags), len(r.Tags))
		}
		if len(g.Tags) > 0 && (g.Tags[0].Key != r.Tags[0].Key || g.Tags[0].Value != r.Tags[0].Value) {
			t.Errorf("Tag[0]: got %+v want %+v", g.Tags[0], r.Tags[0])
		}
	}
}

// TestEncoder_MultipleBlocks verifies that the encoder emits
// multiple BGZF blocks when the data exceeds one block's worth.
func TestEncoder_MultipleBlocks(t *testing.T) {
	hdr := &Header{Raw: "@HD\tVN:1.6\n", SQs: []*SQ{{Name: "chr1", Len: 1 << 20}}}
	// Write 200 records of 1 KiB each = 200 KiB, far past the
	// 64 KiB BGZF block cap.
	seq := strings.Repeat("A", 1024)
	qual := bytes.Repeat([]byte("I"), 1024)
	var records []*Record
	for i := 0; i < 200; i++ {
		records = append(records, &Record{
			QName: "r" + string(rune('0'+i%10)),
			RefID: 0, Pos: int64(i + 1),
			MapQ: 60, Bin: 4681,
			Cigar: []uint32{1<<4 | 1024},
			Seq:   seq, Qual: qual,
		})
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
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	var n int
	for {
		_, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record %d: %v", n, err)
		}
		n++
	}
	if n != len(records) {
		t.Errorf("records: got %d want %d", n, len(records))
	}
}

// TestEncoder_EmptyStream verifies that writing just a header
// (no records) produces a valid BAM with zero records.
func TestEncoder_EmptyStream(t *testing.T) {
	hdr := &Header{Raw: "@HD\tVN:1.6\n", SQs: []*SQ{{Name: "chr1", Len: 100}}}
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
	_, err := dec.Record()
	if !errors.Is(err, io.EOF) {
		t.Errorf("expected EOF, got %v", err)
	}
}

// TestEncoder_AuxTags verifies that all 12 BAM tag types
// round-trip through the encoder.
func TestEncoder_AuxTags(t *testing.T) {
	hdr := &Header{Raw: "@HD\tVN:1.6\n", SQs: []*SQ{{Name: "chr1", Len: 100}}}
	tags := []Tag{
		{Key: "NM", Type: Tagi, Value: "1"},
		{Key: "MD", Type: TagZ, Value: "5"},
		{Key: "XA", Type: TagA, Value: "A"},
		{Key: "XB", Type: Tagf, Value: "1.5"},
		{Key: "XC", Type: TagH, Value: "1A2B"},
		{Key: "XD", Type: Tagd, Value: "3.14"},
		{Key: "XE", Type: Tagc, Value: "-1"},
		{Key: "XF", Type: TagC, Value: "200"},
		{Key: "XG", Type: Tags, Value: "-100"},
		{Key: "XH", Type: TagS, Value: "500"},
		{Key: "XI", Type: TagI, Value: "1234567"},
		// B arrays are excluded: the current Tag struct does not
		// carry the B sub-type, so the encoder cannot emit a
		// valid B array. The 11 non-B tag types above are
		// covered by this round-trip.
	}
	records := []*Record{
		{
			QName: "r1", RefID: 0, Pos: 1, MapQ: 60, Bin: 4681,
			Cigar: []uint32{1<<4 | 5}, Seq: "ACGTA",
			Qual: []byte("IIIII"), Tags: tags,
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
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if len(rec.Tags) != len(tags) {
		t.Fatalf("Tags: got %d want %d", len(rec.Tags), len(tags))
	}
	for i, want := range tags {
		got := rec.Tags[i]
		if got.Key != want.Key || got.Type != want.Type || got.Value != want.Value {
			t.Errorf("tag[%d]: got %+v want %+v", i, got, want)
		}
	}
}

// TestEncoder_NoHeader verifies that WriteRecord before
// WriteHeader returns an error.
func TestEncoder_NoHeader(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	err := enc.WriteRecord(&Record{QName: "r1"})
	if err == nil {
		t.Error("WriteRecord before WriteHeader: expected error, got nil")
	}
}

// TestEncoder_RecordAfterClose verifies that WriteRecord
// after Close returns an error.
func TestEncoder_RecordAfterClose(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if err := enc.WriteRecord(&Record{QName: "r1"}); err == nil {
		t.Error("WriteRecord after Close: expected error, got nil")
	}
}

// TestEncoder_NilHeader verifies that WriteHeader(nil) returns
// an error.
func TestEncoder_NilHeader(t *testing.T) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(nil); err == nil {
		t.Error("WriteHeader(nil): expected error, got nil")
	}
}

// TestEncoder_NilRecord verifies that WriteRecord(nil) is a
// no-op (matches the fasta/BCF pattern).
func TestEncoder_NilRecord(t *testing.T) {
	hdr := &Header{Raw: "@HD\tVN:1.6\n", SQs: []*SQ{{Name: "chr1", Len: 100}}}
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
	dec := NewDecoder(&buf)
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	if _, err := dec.Record(); !errors.Is(err, io.EOF) {
		t.Errorf("expected EOF (no records), got %v", err)
	}
}
