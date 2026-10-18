package bcf

import (
	"bytes"
	"testing"
)

func TestFormatRoundTrip(t *testing.T) {
	// Header declares one FORMAT key (DP) plus one sample.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	rec := &Record{
		CHROM: "chr1",
		Pos:   100,
		Ref:   "A",
		Alt:   []string{"C"},
		Info:  map[string]InfoValue{},
		Format: []uint16{0}, // DP is the only FORMAT key
		Samples: []Sample{
			{Fields: map[uint16]SampleValue{
				0: {Type: "Integer", Ints: []int32{42}},
			}},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(ParseHeader(hdrText)); err != nil {
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
	got, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if len(got.Format) != 1 {
		t.Fatalf("Format: got %d, want 1", len(got.Format))
	}
	if len(got.Samples) != 1 {
		t.Fatalf("Samples: got %d, want 1", len(got.Samples))
	}
	v, ok := got.Samples[0].Fields[0]
	if !ok {
		t.Fatal("DP not in sample")
	}
	if v.Type != "Integer" || len(v.Ints) != 1 || v.Ints[0] != 42 {
		t.Errorf("DP value: got %+v, want Integer [42]", v)
	}
}

func TestFormatMissingRoundTrip(t *testing.T) {
	// A per-sample missing (".") must round-trip as Missing=true.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	rec := &Record{
		CHROM: "chr1", Pos: 100, Ref: "A",
		Info:   map[string]InfoValue{},
		Format: []uint16{0},
		Samples: []Sample{
			{Fields: map[uint16]SampleValue{
				0: {Missing: true},
			}},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(ParseHeader(hdrText)); err != nil {
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
	got, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	v := got.Samples[0].Fields[0]
	if !v.Missing {
		t.Errorf("DP: got %+v, want Missing=true", v)
	}
}

func TestFormatVariableCardinality(t *testing.T) {
	// Number=. means the value count is on the wire.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"phred-scaled likelihoods\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	rec := &Record{
		CHROM: "chr1", Pos: 100, Ref: "A", Alt: []string{"C"},
		Info:   map[string]InfoValue{},
		Format: []uint16{0},
		Samples: []Sample{
			{Fields: map[uint16]SampleValue{
				0: {Type: "Integer", Ints: []int32{0, 100, 200}},
			}},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(ParseHeader(hdrText)); err != nil {
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
	got, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	v := got.Samples[0].Fields[0]
	if len(v.Ints) != 3 || v.Ints[0] != 0 || v.Ints[1] != 100 || v.Ints[2] != 200 {
		t.Errorf("PL: got %+v, want [0 100 200]", v.Ints)
	}
}

func TestFormatNoSamples(t *testing.T) {
	// A record with FORMAT declared but no samples is the
	// "site-only" VCF case: the per-sample block is empty
	// but n_fmt > 0. The decoder should still return a
	// record with Format set and Samples empty.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
	rec := &Record{
		CHROM: "chr1", Pos: 100, Ref: "A",
		Info: map[string]InfoValue{},
		// Format set but Samples empty: site-only.
		Format:  []uint16{0},
		Samples: nil,
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(ParseHeader(hdrText)); err != nil {
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
	got, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if len(got.Format) != 1 {
		t.Errorf("Format: got %d, want 1", len(got.Format))
	}
	if len(got.Samples) != 0 {
		t.Errorf("Samples: got %d, want 0", len(got.Samples))
	}
}