package bai

import "testing"

func TestReference_Query_Empty(t *testing.T) {
	ref := &Reference{}
	q := ref.Query(0, 100)
	if len(q.Chunks) != 0 {
		t.Errorf("empty ref: got %d chunks, want 0", len(q.Chunks))
	}
	if q.MinOffset != 0 {
		t.Errorf("empty ref: MinOffset=%d, want 0", q.MinOffset)
	}
}

func TestReference_Query_ZeroRange(t *testing.T) {
	ref := &Reference{Bins: []*Bin{{Number: 0, Lo: 100, Hi: 200}}}
	q := ref.Query(100, 100)
	if len(q.Chunks) != 0 {
		t.Errorf("zero range: got %d chunks, want 0", len(q.Chunks))
	}
}

func TestReference_Query_InvertedRange(t *testing.T) {
	ref := &Reference{Bins: []*Bin{{Number: 0, Lo: 100, Hi: 200}}}
	q := ref.Query(200, 100)
	if len(q.Chunks) != 0 {
		t.Errorf("inverted range: got %d chunks, want 0", len(q.Chunks))
	}
}

func TestReference_Query_SingleBin(t *testing.T) {
	ref := &Reference{
		Bins: []*Bin{
			{Number: 0, Lo: 100, Hi: 200},
		},
		Intervals: []uint64{150, 0, 0, 0},
	}
	// Region 0..100 hits bin 0 (covers [0, 2^29)).
	q := ref.Query(0, 100)
	if len(q.Chunks) != 1 {
		t.Fatalf("got %d chunks, want 1", len(q.Chunks))
	}
	if q.Chunks[0].Lo != 100 || q.Chunks[0].Hi != 200 {
		t.Errorf("chunk: got (%d,%d), want (100,200)", q.Chunks[0].Lo, q.Chunks[0].Hi)
	}
	if q.MinOffset != 150 {
		t.Errorf("MinOffset: got %d, want 150", q.MinOffset)
	}
}

func TestReference_Query_MultipleBinsDeduped(t *testing.T) {
	// Two regions that hit overlapping bin sets should
	// produce deduplicated chunks.
	ref := &Reference{
		Bins: []*Bin{
			{Number: 0, Lo: 100, Hi: 200},
			{Number: 1, Lo: 50, Hi: 150},
		},
	}
	q := ref.Query(0, 100)
	// A small region hits the leaf bin (level 0) and a
	// handful of parent bins; we should see at most 2
	// distinct chunks here.
	seen := make(map[Chunk]bool)
	for _, c := range q.Chunks {
		seen[c] = true
	}
	if len(seen) != len(q.Chunks) {
		t.Errorf("duplicate chunks in result: %+v", q.Chunks)
	}
}

func TestReference_Query_NoMatchingBins(t *testing.T) {
	// Bins for a different region of the reference.
	ref := &Reference{
		Bins: []*Bin{
			{Number: 4681 + 100, Lo: 100, Hi: 200}, // far bin
		},
	}
	q := ref.Query(0, 100)
	// Region 0..100 hits bins at the bottom of the tree,
	// not 4681+100, so no chunks should be returned.
	if len(q.Chunks) != 0 {
		t.Errorf("no matching bins: got %d chunks, want 0", len(q.Chunks))
	}
}

func TestReference_FindOffset(t *testing.T) {
	ref := &Reference{
		Intervals: []uint64{100, 0, 200, 0, 300},
	}
	// pos=0 -> window 0 -> 100
	if off := ref.FindOffset(0); off != 100 {
		t.Errorf("pos=0: got %d, want 100", off)
	}
	// pos=16384 -> walk 0..1, first non-zero is 100
	if off := ref.FindOffset(16384); off != 100 {
		t.Errorf("pos=16384: got %d, want 100", off)
	}
	// pos=32768 -> walk 0..2, first non-zero is 100
	if off := ref.FindOffset(32768); off != 100 {
		t.Errorf("pos=32768: got %d, want 100", off)
	}
	// pos=80000 (window 4) -> walk 0..4, first non-zero is 100
	if off := ref.FindOffset(80000); off != 100 {
		t.Errorf("pos=80000: got %d, want 100", off)
	}
}

func TestReference_FindOffset_PastEnd(t *testing.T) {
	ref := &Reference{Intervals: []uint64{100, 200}}
	// 2 intervals cover windows 0..1 (positions up to 32768).
	// pos=100000 is well past the end -> 0.
	if off := ref.FindOffset(100000); off != 0 {
		t.Errorf("past end: got %d, want 0", off)
	}
}

func TestReference_FindOffset_Empty(t *testing.T) {
	ref := &Reference{}
	if off := ref.FindOffset(100); off != 0 {
		t.Errorf("empty ref: got %d, want 0", off)
	}
}

func TestReference_Query_LinearIndexSkipZeros(t *testing.T) {
	// Window 0 has no alignment; the first one starts at
	// window 3. A query that reaches into window 3 should
	// report the offset of the first non-zero window.
	ref := &Reference{
		Intervals: []uint64{0, 0, 0, 50, 60},
	}
	// Query [0, 200000) reaches into windows 0..12, so
	// the first non-zero interval is 50.
	q := ref.Query(0, 200000)
	if q.MinOffset != 50 {
		t.Errorf("MinOffset: got %d, want 50", q.MinOffset)
	}
	// Query [0, 100) covers only window 0 — no alignment
	// starts there, MinOffset stays 0.
	q = ref.Query(0, 100)
	if q.MinOffset != 0 {
		t.Errorf("narrow query: MinOffset=%d, want 0", q.MinOffset)
	}
}