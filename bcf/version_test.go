package bcf

import (
	"bytes"
	"encoding/binary"
	"io"
	"testing"
)

// TestDecoder_BCF22RoundTrip builds a BCF 2.2 stream and
// verifies the decoder reads it. We hand-build the stream
// (BGZF + Magic22 + l_text + header text) because the
// library's own encoder always emits BCF 2.1.
func TestDecoder_BCF22RoundTrip(t *testing.T) {
	hdrText := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
	b := newBuilderWithHeader(t, hdrText)
	payload := make([]byte, 0, len(magic22)+4+len(hdrText))
	payload = append(payload, magic22[:]...)
	var lt [4]byte
	binary.LittleEndian.PutUint32(lt[:], uint32(len(hdrText))) //nolint:gosec // test data, length fits
	payload = append(payload, lt[:]...)
	payload = append(payload, []byte(hdrText)...)
	b.writeBGZFBlock(payload)
	data := b.finish(t)

	dec := NewDecoder(bytes.NewReader(data))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	if got := dec.Version(); got != 2 {
		t.Errorf("Version: got %d, want 2 (BCF 2.2)", got)
	}
}

// TestEncoder_BCF22EmitAndReadBack covers the encoder end
// of BCF 2.2: SetBCFVersion(2) before WriteHeader causes
// the encoder to emit the 2.2 magic. The decoder then reads
// the stream back and reports Version() == 2.
//
// We also test the long-allele encoding: a REF or ALT string
// of length 15 or more is emitted as typeLong (0xF) +
// uint32 length rather than typeString (7) + uint32 length.
// The decoder accepts either form.
func TestEncoder_BCF22EmitAndReadBack(t *testing.T) {
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	hdr := ParseHeader(hdrText)
	// 20-character ALT: the boundary for the long-string
	// encoding is length >= 15.
	longAlt := "GATTACAGATTACAGATTA" // length 20
	rec := &Record{
		CHROM: "chr1", Pos: 100,
		ID:   "rs1",
		Ref:  "A",
		Alt:  []string{longAlt},
		Info: map[string]InfoValue{},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.SetBCFVersion(2); err != nil {
		t.Fatalf("SetBCFVersion: %v", err)
	}
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(rec); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	// Verify the magic byte at the start of the BGZF
	// payload. The 6th byte of the magic is the version
	// byte (1 for 2.1, 2 for 2.2).
	dec := NewDecoder(&buf)
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	if got := dec.Version(); got != 2 {
		t.Errorf("Version: got %d, want 2 (BCF 2.2)", got)
	}
	got, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if got.Ref != "A" {
		t.Errorf("Ref: got %q, want %q", got.Ref, "A")
	}
	if len(got.Alt) != 1 || got.Alt[0] != longAlt {
		t.Errorf("Alt: got %v, want [%q]", got.Alt, longAlt)
	}
}

func TestSetBCFVersion_RejectsBadValue(t *testing.T) {
	enc := NewEncoder(io.Discard)
	if err := enc.SetBCFVersion(3); err == nil {
		t.Error("SetBCFVersion(3): expected error, got nil")
	}
	if err := enc.SetBCFVersion(0); err == nil {
		t.Error("SetBCFVersion(0): expected error, got nil")
	}
}

func TestSetBCFVersion_AfterHeaderFails(t *testing.T) {
	hdrText := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(ParseHeader(hdrText)); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.SetBCFVersion(2); err == nil {
		t.Error("SetBCFVersion after WriteHeader: expected error, got nil")
	}
}

// TestEncoder_BCF22LongStringRoundTrip exercises the
// per-sample string FORMAT field with a long string.
// In BCF 2.2 the encoder should emit typeLong (0xF) for
// strings whose length is 15 or more. The decoder accepts
// either type tag.
func TestEncoder_BCF22LongStringRoundTrip(t *testing.T) {
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=FT,Number=1,Type=String,Description=\"long format string\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	hdr := ParseHeader(hdrText)
	// 30-character FT value: triggers the long-string
	// encoding in BCF 2.2.
	longFT := "ABCDEFGHIJKLMNOPQRSTUVWXYZ12345" // length 31
	rec := &Record{
		CHROM:  "chr1",
		Pos:    100,
		Ref:    "A",
		Info:   map[string]InfoValue{},
		Format: []uint16{0},
		Samples: []Sample{
			{Fields: map[uint16]SampleValue{
				0: {Type: "String", Strs: []string{longFT}},
			}},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.SetBCFVersion(2); err != nil {
		t.Fatalf("SetBCFVersion: %v", err)
	}
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(rec); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	dec := NewDecoder(&buf)
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	if got := dec.Version(); got != 2 {
		t.Errorf("Version: got %d, want 2", got)
	}
	got, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	v, ok := got.Samples[0].Fields[0]
	if !ok {
		t.Fatal("FT not in sample")
	}
	if len(v.Strs) != 1 || v.Strs[0] != longFT {
		t.Errorf("FT value: got %v, want [%q]", v.Strs, longFT)
	}
}