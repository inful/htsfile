package csi

import (
	"io"

	"bytes"
	"errors"
	"github.com/inful/htsfile/internal/bgzf"
	"testing"
)

func TestWriter_RoundTrip_Basic(t *testing.T) {
	original := &Index{
		MinShift: 14,
		Depth:    5,
		References: []*Reference{
			{
				Bins: map[uint32][]Chunk{
					4681: {{Beg: 0x1000, End: 0x2000}},
				},
				Loffsets: map[uint32]VirtualOffset{4681: 0x1000},
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
	if idx.MinShift != 14 || idx.Depth != 5 {
		t.Errorf("min_shift/depth: got %d/%d want 14/5", idx.MinShift, idx.Depth)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References: got %d want 1", len(idx.References))
	}
	ref := idx.References[0]
	if got := ref.Bins[4681]; len(got) != 1 || got[0].Beg != 0x1000 || got[0].End != 0x2000 {
		t.Errorf("Bins[4681]: got %+v", got)
	}
	// A CSI has no linear index: the per-bin loffset is what
	// says where the bin's records start.
	if got := ref.Loffsets[4681]; got != 0x1000 {
		t.Errorf("Loffsets[4681]: got %#x want 0x1000", got)
	}
}

func TestWriter_RoundTrip_TunableParams(t *testing.T) {
	// Non-default min_shift / depth — verifies the header
	// bytes are written as little-endian int32s.
	original := &Index{
		MinShift: 18,
		Depth:    6,
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
	if idx.MinShift != 18 || idx.Depth != 6 {
		t.Errorf("min_shift/depth: got %d/%d want 18/6", idx.MinShift, idx.Depth)
	}
}

func TestWriter_RoundTrip_MultipleRefs(t *testing.T) {
	original := &Index{
		MinShift: 14,
		Depth:    5,
		References: []*Reference{
			{
				Bins: map[uint32][]Chunk{
					1: {{Beg: VirtualOffset(0x1000), End: VirtualOffset(0x2000)}},
				},
				Intervals: []VirtualOffset{0x1000},
			},
			{
				Bins: map[uint32][]Chunk{
					2: {
						{Beg: VirtualOffset(0x3000), End: VirtualOffset(0x4000)},
						{Beg: VirtualOffset(0x4000), End: VirtualOffset(0x5000)},
					},
				},
				Intervals: []VirtualOffset{0x3000, 0x4000, 0x5000},
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
	if len(idx.References[1].Bins[2]) != 2 {
		t.Errorf("ref[1].Bins[2]: got %d chunks want 2", len(idx.References[1].Bins[2]))
	}
}

func TestWriter_EmptyIndex(t *testing.T) {
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(&Index{MinShift: 14, Depth: 5}); err != nil {
		t.Fatalf("Write: %v", err)
	}
	idx, err := NewReader(&buf).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 0 {
		t.Errorf("References: got %d want 0", len(idx.References))
	}
	if idx.MinShift != 14 || idx.Depth != 5 {
		t.Errorf("min_shift/depth: got %d/%d want 14/5", idx.MinShift, idx.Depth)
	}
}

func TestWriter_MagicAndHeaderBytes(t *testing.T) {
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(&Index{MinShift: 14, Depth: 5}); err != nil {
		t.Fatalf("Write: %v", err)
	}
	// The index is BGZF-compressed, the way htslib writes
	// it, so the magic sits inside the compressed stream.
	zr := bgzf.NewReader(bytes.NewReader(buf.Bytes()))
	got := make([]byte, 4)
	if _, err := io.ReadFull(zr, got); err != nil {
		t.Fatalf("read magic: %v", err)
	}
	if !bytes.Equal(got, Magic[:]) {
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
	if err := NewWriter(errWriter{}).Write(&Index{MinShift: 14, Depth: 5}); err == nil {
		t.Error("Write to failing writer: expected error, got nil")
	}
}

func TestWriter_EmptyBinsAndIntervals(t *testing.T) {
	original := &Index{
		MinShift: 14,
		Depth:    5,
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
	if len(idx.References) != 1 || len(idx.References[0].Bins) != 0 ||
		len(idx.References[0].Intervals) != 0 {
		t.Errorf("expected one ref with empty bins + intervals, got %+v", idx.References[0])
	}
}

func TestWriter_ChunkOrderingPreserved(t *testing.T) {
	original := &Index{
		MinShift: 14,
		Depth:    5,
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
	if got[0].Beg != 0x0500 || got[1].Beg != 0x1000 || got[2].Beg != 0x2000 {
		t.Errorf("chunks not sorted by Beg: %+v", got)
	}
}

// errWriter is an io.Writer that always returns an error.
type errWriter struct{}

func (errWriter) Write(p []byte) (int, error) { return 0, errors.New("synthetic") }
