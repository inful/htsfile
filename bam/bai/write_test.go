package bai

import (
	"bytes"
	"testing"
)

func TestWriter_RoundTrip_Basic(t *testing.T) {
	original := &Index{
		References: []*Reference{
			{
				Bins: []*Bin{
					{Number: 4681, Lo: 0x1000, Hi: 0x2000},
				},
				Intervals: []uint64{0x1000, 0x2000, 0x3000},
				NoCoor:    0,
			},
		},
	}
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	idx, err := NewReader(&buf).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References: got %d want 1", len(idx.References))
	}
	ref := idx.References[0]
	if len(ref.Bins) != 1 || ref.Bins[0].Number != 4681 ||
		ref.Bins[0].Lo != 0x1000 || ref.Bins[0].Hi != 0x2000 {
		t.Errorf("Bins: got %+v", ref.Bins)
	}
	if len(ref.Intervals) != 3 {
		t.Errorf("Intervals: got %v", ref.Intervals)
	}
	if ref.NoCoor != 0 {
		t.Errorf("NoCoor: got %d want 0", ref.NoCoor)
	}
}

func TestWriter_RoundTrip_MultipleRefs(t *testing.T) {
	original := &Index{
		References: []*Reference{
			{
				Bins: []*Bin{
					{Number: 1, Lo: 0x1000, Hi: 0x2000},
				},
				Intervals: []uint64{0x1000},
				NoCoor:    5,
			},
			{
				Bins: []*Bin{
					{Number: 2, Lo: 0x3000, Hi: 0x4000},
					{Number: 100, Lo: 0x4000, Hi: 0x5000},
				},
				Intervals: []uint64{0x3000, 0x4000, 0x5000},
				NoCoor:    7,
			},
		},
	}
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	idx, err := NewReader(&buf).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 2 {
		t.Fatalf("References: got %d want 2", len(idx.References))
	}
	if len(idx.References[1].Bins) != 2 {
		t.Errorf("ref[1].Bins: got %d want 2", len(idx.References[1].Bins))
	}
	if idx.References[0].NoCoor != 5 || idx.References[1].NoCoor != 7 {
		t.Errorf("NoCoor values wrong: %d %d", idx.References[0].NoCoor, idx.References[1].NoCoor)
	}
}

func TestWriter_EmptyIndex(t *testing.T) {
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(&Index{}); err != nil {
		t.Fatalf("Write: %v", err)
	}
	idx, err := NewReader(&buf).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 0 {
		t.Errorf("References: got %d want 0", len(idx.References))
	}
}

func TestWriter_MagicAndHeaderBytes(t *testing.T) {
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(&Index{}); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if got := buf.Bytes()[:4]; !bytes.Equal(got, Magic[:]) {
		t.Errorf("magic: got %q want %q", got, Magic[:])
	}
}

func TestWriter_BinsSorted(t *testing.T) {
	// Bins within a reference must be in ascending order in the
	// emitted file (the Reader assumes it). The Writer
	// sorts them by Number.
	original := &Index{
		References: []*Reference{
			{
				Bins: []*Bin{
					{Number: 100, Lo: 0x4000, Hi: 0x5000},
					{Number: 1, Lo: 0x1000, Hi: 0x2000},
					{Number: 50, Lo: 0x3000, Hi: 0x4000},
				},
			},
		},
	}
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	idx, err := NewReader(&buf).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	bins := idx.References[0].Bins
	if bins[0].Number != 1 || bins[1].Number != 50 || bins[2].Number != 100 {
		t.Errorf("bins not sorted by Number: %+v", bins)
	}
}

func TestWriter_NilIndex(t *testing.T) {
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(nil); err == nil {
		t.Error("Write(nil): expected error, got nil")
	}
}

func TestWriter_FailingWriter(t *testing.T) {
	if err := NewWriter(errWriter{}).Write(&Index{}); err == nil {
		t.Error("Write to failing writer: expected error, got nil")
	}
}

// errWriter is an io.Writer that always returns an error.
type errWriter struct{}

func (errWriter) Write(p []byte) (int, error) { return 0, &writeError{} }

// writeError is a sentinel error for failingWriter.
type writeError struct{}

func (*writeError) Error() string { return "synthetic write failure" }
