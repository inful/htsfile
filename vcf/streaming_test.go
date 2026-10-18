package vcf

import (
	"fmt"
	"strings"
	"testing"
)

// streamingTestHeader is a small header suitable for streaming
// tests. It declares GT and M5mC so tests can exercise both
// the GT structural check and the M-family value range check.
const streamingTestHeader = "##fileformat=VCFv4.2\n" +
	"##contig=<ID=chr1,length=1000>\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
	"##FORMAT=<ID=M5mC,Number=.,Type=Float,Description=\"5mC fraction\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"

// streamRecord builds a record by feeding a one-line VCF
// through the parser. This gives us the same Raw/Strs/Floats
// layout that real records have, without hand-rolling it.
// pos is the 1-based record position; gt is the GT string;
// m5mc is the M5mC value (empty = unset). The REF is always
// "C"; tests do not need to vary it.
func streamRecord(t *testing.T, pos int64, gt, m5mc string) *Record {
	t.Helper()
	sample := gt
	if m5mc != "" {
		sample = gt + ":" + m5mc
	}
	line := "chr1\t" + itoa(pos) + "\t.\tC\tA\t.\t.\t.\tGT:M5mC\t" + sample + "\n"
	dec := NewDecoder(strings.NewReader(streamingTestHeader + line))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	return rec
}

func itoa(n int64) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}

func TestValidator_Streaming_HeaderOnce(t *testing.T) {
	dec := NewDecoder(strings.NewReader(streamingTestHeader))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	v := NewValidator(ModeDefault)
	// First call runs the header checks (none of the
	// well-formed-header issues should fire here).
	if issues := v.HeaderChecks(hdr); len(issues) != 0 {
		t.Errorf("first HeaderChecks: got %d issues, want 0: %+v", len(issues), issues)
	}
	// Second call must be a no-op (headerDone flag set).
	if issues := v.HeaderChecks(hdr); issues != nil {
		t.Errorf("second HeaderChecks should be a no-op, got: %+v", issues)
	}
}

func TestValidator_Streaming_NextReturnsPerRecordIssues(t *testing.T) {
	dec := NewDecoder(strings.NewReader(streamingTestHeader))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	v := NewValidator(ModeDefault)
	_ = v.HeaderChecks(hdr)

	// Record 1: clean.
	clean := streamRecord(t, 100, "0/0", "0.5")
	if issues := v.Next(hdr, clean); len(issues) != 0 {
		t.Errorf("clean record: got %d issues, want 0: %+v", len(issues), issues)
	}
	// Record 2: M5mC out of range (1.5).
	bad := streamRecord(t, 200, "0/0", "1.5")
	issues := v.Next(hdr, bad)
	found := false
	for _, iss := range issues {
		if iss.Field == "FORMAT/M5mC/sample[0]" {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("bad M5mC record: expected M5mC issue, got: %+v", issues)
	}
	// Record 3: GT mixes separators.
	mixed := streamRecord(t, 300, "0/1|2", "0.5")
	issues = v.Next(hdr, mixed)
	found = false
	for _, iss := range issues {
		if iss.Field == "FORMAT/GT/sample[0]" {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("mixed GT record: expected GT issue, got: %+v", issues)
	}
}

func TestValidator_Streaming_POSMonotonicityIncremental(t *testing.T) {
	dec := NewDecoder(strings.NewReader(streamingTestHeader))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	v := NewValidator(ModeStrict)
	_ = v.HeaderChecks(hdr)

	// First record: POS 100, no previous.
	first := streamRecord(t, 100, "0/0", "0.5")
	if issues := v.Next(hdr, first); len(issues) != 0 {
		t.Errorf("first record: got %d issues, want 0: %+v", len(issues), issues)
	}
	// Second record: POS 50 (out of order). The streaming
	// check should fire here, not deferred to Finish.
	bad := streamRecord(t, 50, "0/0", "0.5")
	issues := v.Next(hdr, bad)
	found := false
	for _, iss := range issues {
		if iss.Field == "POS" {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("out-of-order POS: expected POS issue at Next, got: %+v", issues)
	}
}

func TestValidator_Streaming_Reset(t *testing.T) {
	dec := NewDecoder(strings.NewReader(streamingTestHeader))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	v := NewValidator(ModeStrict)
	_ = v.HeaderChecks(hdr)
	// Feed a record that bumps recordIdx and prevPos.
	_ = v.Next(hdr, streamRecord(t, 200, "0/0", "0.5"))
	// Reset and feed a new first record. The monotonicity
	// check should not fire on the very first record of the
	// new run.
	v.Reset()
	_ = v.HeaderChecks(hdr)
	first := streamRecord(t, 100, "0/0", "0.5")
	if issues := v.Next(hdr, first); len(issues) != 0 {
		t.Errorf("after Reset, first record: got %d issues, want 0: %+v", len(issues), issues)
	}
}

func TestValidator_Streaming_EquivalenceWithValidate(t *testing.T) {
	// The streaming API must produce the same issues as the
	// batch Validate API for the same inputs.
	dec := NewDecoder(strings.NewReader(streamingTestHeader))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	recs := []*Record{
		streamRecord(t, 100, "0/0", "0.5"),
		streamRecord(t, 200, "0/1|2", "0.5"), // mixed separator
		streamRecord(t, 300, "0/0", "1.5"),  // out of range
		streamRecord(t, 50, "0/0", "0.5"),   // out of order
	}
	// Batch.
	batch := NewValidator(ModeStrict).Validate(hdr, recs)
	// Streaming.
	s := NewValidator(ModeStrict)
	stream := make([]ValidationError, 0, 8)
	stream = append(stream, s.HeaderChecks(hdr)...)
	for _, r := range recs {
		stream = append(stream, s.Next(hdr, r)...)
	}
	stream = append(stream, s.Finish()...)

	if len(batch) != len(stream) {
		t.Fatalf("issue count mismatch: batch=%d, stream=%d\nbatch: %+v\nstream: %+v",
			len(batch), len(stream), batch, stream)
	}
	for i := range batch {
		if batch[i].Field != stream[i].Field || batch[i].Msg != stream[i].Msg {
			t.Errorf("issue %d differs: batch=%+v, stream=%+v", i, batch[i], stream[i])
		}
	}
}

// streamingTestHeaderPS declares GT and PS in the FORMAT
// column so tests can exercise the cross-record PS check.
const streamingTestHeaderPS = "##fileformat=VCFv4.2\n" +
	"##contig=<ID=chr1,length=1000>\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
	"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"phase set\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"

// streamRecordPS builds a record that includes PS for the
// streaming phase-style test.
func streamRecordPS(t *testing.T, pos int64, gt string, ps int64) *Record {
	t.Helper()
	sample := gt
	if ps > 0 {
		sample = gt + ":" + fmt.Sprintf("%d", ps)
	}
	line := "chr1\t" + itoa(pos) + "\t.\tC\tA\t.\t.\t.\tGT:PS\t" + sample + "\n"
	dec := NewDecoder(strings.NewReader(streamingTestHeaderPS + line))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	return rec
}

func TestValidator_Streaming_PhaseStyleMismatch(t *testing.T) {
	// Two records on the same chrom with the same PS value
	// but different phasing styles (one phased, one
	// unphased) must produce a cross-record issue in
	// ModeStrict.
	dec := NewDecoder(strings.NewReader(streamingTestHeaderPS))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	v := NewValidator(ModeStrict)
	_ = v.HeaderChecks(hdr)
	issues := v.Next(hdr, streamRecordPS(t, 100, "0|1", 42))
	if len(issues) != 0 {
		t.Errorf("first record: got %d issues, want 0: %+v", len(issues), issues)
	}
	issues = v.Next(hdr, streamRecordPS(t, 200, "0/1", 42))
	found := false
	for _, iss := range issues {
		if iss.Field == "FORMAT/GT/sample[0]" && strings.Contains(iss.Msg, "separator style") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("second record: expected phase-style issue, got: %+v", issues)
	}
	// A third record with a different PS (or no PS) on the
	// same chrom must not fire.
	issues = v.Next(hdr, streamRecordPS(t, 300, "0|1", 99))
	for _, iss := range issues {
		if iss.Field == "FORMAT/GT/sample[0]" && strings.Contains(iss.Msg, "separator style") {
			t.Errorf("third record (different PS): unexpected phase-style issue: %+v", iss)
		}
	}
}

func TestValidator_Streaming_PhaseStyleConsistent(t *testing.T) {
	// Two records with the same PS and consistent phasing
	// style must not produce any phase-style issue.
	dec := NewDecoder(strings.NewReader(streamingTestHeaderPS))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	v := NewValidator(ModeStrict)
	_ = v.HeaderChecks(hdr)
	for _, gt := range []string{"0|1", "1|0", "0|1|2"} {
		issues := v.Next(hdr, streamRecordPS(t, 100, gt, 42))
		for _, iss := range issues {
			if iss.Field == "FORMAT/GT/sample[0]" && strings.Contains(iss.Msg, "separator style") {
				t.Errorf("GT=%q: unexpected phase-style issue: %+v", gt, iss)
			}
		}
	}
}