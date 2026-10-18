package tabix

import (
	"bytes"
	"errors"
	"testing"
)

func TestWriter_RoundTrip_Basic(t *testing.T) {
	original := &Index{
		References: []*Reference{
			{
				Bins: map[uint32][]Chunk{
					4681: {{Beg: 0x1000, End: 0x2000}},
				},
				Intervals: []VirtualOffset{0x1000, 0x2000, 0x3000},
				NoCoor:    0,
			},
		},
	}

	// Write to a buffer.
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}

	// Read it back.
	idx, err := NewReader(&buf).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References: got %d want 1", len(idx.References))
	}
	ref := idx.References[0]
	if got := ref.Bins[4681]; len(got) != 1 || got[0].Beg != 0x1000 || got[0].End != 0x2000 {
		t.Errorf("Bins[4681]: got %+v", got)
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
				Bins: map[uint32][]Chunk{
					1: {{Beg: VirtualOffset(0x1000), End: VirtualOffset(0x2000)}},
				},
				Intervals: asVO([]uint64{0x1000}),
				NoCoor:    3,
			},
			{
				Bins: map[uint32][]Chunk{
					2: {
						{Beg: VirtualOffset(0x3000), End: VirtualOffset(0x4000)},
						{Beg: VirtualOffset(0x4000), End: VirtualOffset(0x5000)},
					},
				},
				Intervals: asVO([]uint64{0x3000, 0x4000, 0x5000}),
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
	if idx.References[0].NoCoor != 3 {
		t.Errorf("ref[0].NoCoor: got %d want 3", idx.References[0].NoCoor)
	}
	if len(idx.References[1].Bins[2]) != 2 {
		t.Errorf("ref[1].Bins[2]: got %d chunks want 2", len(idx.References[1].Bins[2]))
	}
}

func TestWriter_EmptyIndex(t *testing.T) {
	// No references — just the magic + n_ref=0 header.
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
	// The first 4 bytes of any Tabix file are the magic.
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(&Index{}); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if got := buf.Bytes()[:4]; !bytes.Equal(got, Magic[:]) {
		t.Errorf("magic: got %q want %q", got, Magic[:])
	}
}

func TestWriter_NilIndex(t *testing.T) {
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(nil); err == nil {
		t.Error("Write(nil): expected error, got nil")
	}
}

func TestWriter_FailingWriter(t *testing.T) {
	// A writer that always errors surfaces the error from Write.
	if err := NewWriter(errWriter{}).Write(&Index{}); err == nil {
		t.Error("Write to failing writer: expected error, got nil")
	}
}

func TestWriter_EmptyBinsAndIntervals(t *testing.T) {
	// Reference with no bins and no intervals — degenerate but
	// legal (writer emits n_bin=0 and n_intv=0).
	original := &Index{
		References: []*Reference{
			{
				Bins:      map[uint32][]Chunk{},
				Intervals: nil,
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
		t.Fatalf("References: got %d", len(idx.References))
	}
	if len(idx.References[0].Bins) != 0 || len(idx.References[0].Intervals) != 0 {
		t.Errorf("expected empty bins + intervals")
	}
}

func TestWriter_ChunkOrderingPreserved(t *testing.T) {
	// Chunk order within a bin matters — round-trip must
	// preserve it. The map iteration is non-deterministic,
	// but within one bin we sort by Beg for deterministic
	// output (matching how the Reader sorts).
	original := &Index{
		References: []*Reference{
			{
				Bins: map[uint32][]Chunk{
					100: {
						{Beg: 0x1000, End: 0x2000},
						{Beg: 0x0500, End: 0x0800},
						{Beg: 0x2000, End: 0x3000},
					},
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
	got := idx.References[0].Bins[100]
	if len(got) != 3 {
		t.Fatalf("got %d chunks want 3", len(got))
	}
	// Sorted by Beg ascending.
	if got[0].Beg != 0x0500 || got[1].Beg != 0x1000 || got[2].Beg != 0x2000 {
		t.Errorf("chunks not sorted by Beg: %+v", got)
	}
}

// errWriter is an io.Writer that always returns an error.
type errWriter struct{}

func (errWriter) Write(p []byte) (int, error) { return 0, errors.New("synthetic") }
