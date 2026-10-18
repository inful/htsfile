package tabix

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"io"
	"testing"
)

// makeIndex builds a synthetic Tabix index file in memory.
func makeIndex(refs []*Reference, withNoCoor bool) []byte {
	var buf bytes.Buffer
	buf.Write(Magic[:])
	var nRef [4]byte
	binary.LittleEndian.PutUint32(nRef[:], uint32(len(refs)))
	buf.Write(nRef[:])
	for _, ref := range refs {
		var nBin [4]byte
		binary.LittleEndian.PutUint32(nBin[:], uint32(len(ref.Bins)))
		buf.Write(nBin[:])
		for bin, chunks := range ref.Bins {
			var bBytes [4]byte
			binary.LittleEndian.PutUint32(bBytes[:], bin)
			buf.Write(bBytes[:])
			var nChunk [4]byte
			binary.LittleEndian.PutUint32(nChunk[:], uint32(len(chunks)))
			buf.Write(nChunk[:])
			for _, c := range chunks {
				var begBytes [8]byte
				binary.LittleEndian.PutUint64(begBytes[:], uint64(c.Beg))
				buf.Write(begBytes[:])
				var endBytes [8]byte
				binary.LittleEndian.PutUint64(endBytes[:], uint64(c.End))
				buf.Write(endBytes[:])
			}
		}
		var nIntv [4]byte
		binary.LittleEndian.PutUint32(nIntv[:], uint32(len(ref.Intervals)))
		buf.Write(nIntv[:])
		for _, iv := range ref.Intervals {
			var ivBytes [8]byte
			binary.LittleEndian.PutUint64(ivBytes[:], uint64(iv))
			buf.Write(ivBytes[:])
		}
		if withNoCoor {
			var ncBytes [8]byte
			binary.LittleEndian.PutUint64(ncBytes[:], ref.NoCoor)
			buf.Write(ncBytes[:])
		}
	}
	return buf.Bytes()
}

func TestRead_Basic(t *testing.T) {
	refs := []*Reference{
		{
			Bins: map[uint32][]Chunk{
				4681: {{Beg: 0x1000, End: 0x2000}},
			},
			Intervals: []VirtualOffset{0x1000, 0x2000, 0x3000},
		},
	}
	data := makeIndex(refs, true)
	idx, err := NewReader(bufio.NewReader(bytes.NewReader(data))).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References: %d", len(idx.References))
	}
	ref := idx.References[0]
	if got := ref.Bins[4681]; len(got) != 1 || got[0].Beg != 0x1000 || got[0].End != 0x2000 {
		t.Errorf("Bins[4681]: %+v", got)
	}
	if len(ref.Intervals) != 3 {
		t.Errorf("Intervals: %v", ref.Intervals)
	}
	if ref.NoCoor != 0 {
		t.Errorf("NoCoor: %d (expected 0 when builder skips)", ref.NoCoor)
	}
}

func TestRead_MultipleRefs(t *testing.T) {
	mkRef := func(bin uint32, beg, end uint64, ivs []uint64) *Reference {
		return &Reference{
			Bins:      map[uint32][]Chunk{bin: {{Beg: VirtualOffset(beg), End: VirtualOffset(end)}}},
			Intervals: asVO(ivs),
		}
	}
	refs := []*Reference{
		mkRef(1, 0x1000, 0x2000, []uint64{0x1000}),
		mkRef(2, 0x3000, 0x4000, []uint64{0x3000, 0x4000, 0x5000}),
	}
	// Use withNoCoor=true (the modern Tabix format includes
	// n_no_coor after every reference's intervals).
	data := makeIndex(refs, true)
	idx, err := NewReader(bufio.NewReader(bytes.NewReader(data))).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 2 {
		t.Fatalf("References: %d", len(idx.References))
	}
	if idx.References[1].Bins[2][0].Beg != 0x3000 {
		t.Errorf("ref[1] bin 2: %+v", idx.References[1].Bins[2])
	}
	// Verify NoCoor values.
	if idx.References[0].NoCoor != 0 {
		t.Errorf("ref[0].NoCoor: %d", idx.References[0].NoCoor)
	}
	if idx.References[1].NoCoor != 0 {
		t.Errorf("ref[1].NoCoor: %d", idx.References[1].NoCoor)
	}
}

func TestRead_BadMagic(t *testing.T) {
	data := []byte{'X', 'B', 'I', 0x01, 0, 0, 0, 0}
	_, err := NewReader(bufio.NewReader(bytes.NewReader(data))).Read()
	if err == nil {
		t.Fatal("expected bad-magic error")
	}
	var re *ReadError
	if !errors.As(err, &re) {
		t.Errorf("error type: %v", err)
	}
}

func TestRead_Truncated(t *testing.T) {
	_, err := NewReader(bytes.NewReader([]byte{'T', 'B', 'I'})).Read()
	if err == nil {
		t.Fatal("expected truncated-magic error")
	}
}

func TestRead_EmptyStream(t *testing.T) {
	_, err := NewReader(bytes.NewReader(nil)).Read()
	if err == nil {
		t.Fatal("expected empty-stream error")
	}
}

func TestReg2Bins(t *testing.T) {
	// Single base: reg2bins returns one bin per level
	// (5 bins for a 1-bp region in samtools 5-level scheme).
	bins := reg2bins(100, 101)
	if len(bins) != 5 {
		t.Errorf("reg2bins(100, 101): got %d bins, want 5", len(bins))
	}
	// Always includes bin 0 (covers the entire reference).
	if bins[0] != 0 {
		t.Errorf("reg2bins(100, 101): first bin %d, want 0", bins[0])
	}
	// Empty range: no bins.
	if bins := reg2bins(100, 100); len(bins) != 0 {
		t.Errorf("reg2bins(100, 100): got %v, want empty", bins)
	}
	// Wide range: many bins across multiple levels.
	bins = reg2bins(0, 100_000)
	if len(bins) == 0 {
		t.Error("reg2bins(0, 100000): got empty")
	}
	// All bins should be in 0..4680 (samtools 5-level scheme).
	for _, b := range bins {
		if b > 4680 {
			t.Errorf("reg2bins(0, 100000): bin %d > 4680", b)
		}
	}
}

func TestLinearIndexOffset(t *testing.T) {
	tests := []struct {
		pos  int64
		want int64
	}{
		{0, 0},
		{16383, 0},
		{16384, 1},
		{16385, 1},
		{100000, 6},
	}
	for _, tc := range tests {
		if got := LinearIndexOffset(tc.pos); got != tc.want {
			t.Errorf("LinearIndexOffset(%d): got %d want %d", tc.pos, got, tc.want)
		}
	}
}

func TestBinFromRange(t *testing.T) {
	// Wrapper around reg2bins; just smoke-test the public path.
	bins := BinFromRange(1000, 5000)
	if len(bins) == 0 {
		t.Error("BinFromRange returned empty")
	}
}

func TestReference_LookupChunks(t *testing.T) {
	ref := &Reference{
		Bins: map[uint32][]Chunk{
			1: {{Beg: 0x1000, End: 0x2000}},
			2: {{Beg: 0x2000, End: 0x3000}},
		},
	}
	// Empty bins → return all.
	all := ref.LookupChunks(nil)
	if len(all) != 2 {
		t.Errorf("all chunks: got %d", len(all))
	}
	// Specific bins → return those.
	some := ref.LookupChunks([]uint32{1})
	if len(some) != 1 || some[0].Beg != 0x1000 {
		t.Errorf("bin 1 chunks: %+v", some)
	}
}

func TestReadError_Error(t *testing.T) {
	e := &ReadError{Err: "bad"}
	if msg := e.Error(); !contains(msg, "tabix") || !contains(msg, "bad") {
		t.Errorf("Error: %q", msg)
	}
	if e.Unwrap() != nil {
		t.Errorf("Unwrap: %v", e.Unwrap())
	}
}

// contains is a tiny strings.Contains shim.
func contains(s, sub string) bool {
	for i := 0; i+len(sub) <= len(s); i++ {
		if s[i:i+len(sub)] == sub {
			return true
		}
	}
	return false
}

// asVO converts a []uint64 to []VirtualOffset.
func asVO(in []uint64) []VirtualOffset {
	out := make([]VirtualOffset, len(in))
	for i, v := range in {
		out[i] = VirtualOffset(v)
	}
	return out
}

// _ keeps io import visible (used by io.EOF in other tests).
var _ = io.EOF
