package csi

import (
	"bytes"
	"encoding/binary"
	"errors"
	"testing"
)

// makeIndex builds a synthetic CSI index.
func makeIndex(minShift int32, refs []*Reference) []byte {
	var buf bytes.Buffer
	buf.Write(Magic[:])
	var m [4]byte
	binary.LittleEndian.PutUint32(m[:], uint32(minShift))
	buf.Write(m[:])
	// depth defaults to 5 in the samtools/htslib CSI writer.
	binary.LittleEndian.PutUint32(m[:], 5)
	buf.Write(m[:])
	// l_aux = 0 (no auxiliary data).
	binary.LittleEndian.PutUint32(m[:], 0)
	buf.Write(m[:])
	var nr [4]byte
	binary.LittleEndian.PutUint32(nr[:], uint32(len(refs)))
	buf.Write(nr[:])
	for _, ref := range refs {
		var nBin [4]byte
		binary.LittleEndian.PutUint32(nBin[:], uint32(len(ref.Bins)))
		buf.Write(nBin[:])
		for bin, chunks := range ref.Bins {
			var bBytes [4]byte
			binary.LittleEndian.PutUint32(bBytes[:], bin)
			buf.Write(bBytes[:])
			// Per-bin loffset (uint64), written as 0 — the
			// reader discards the value but must consume it.
			var zero [8]byte
			buf.Write(zero[:])
			var nChunk [4]byte
			binary.LittleEndian.PutUint32(nChunk[:], uint32(len(chunks)))
			buf.Write(nChunk[:])
			for _, c := range chunks {
				var begBytes [8]byte
				binary.LittleEndian.PutUint64(begBytes[:], c.Beg)
				buf.Write(begBytes[:])
				var endBytes [8]byte
				binary.LittleEndian.PutUint64(endBytes[:], c.End)
				buf.Write(endBytes[:])
			}
		}
		var nIntv [4]byte
		binary.LittleEndian.PutUint32(nIntv[:], uint32(len(ref.Intervals)))
		buf.Write(nIntv[:])
		for _, iv := range ref.Intervals {
			var ivBytes [8]byte
			binary.LittleEndian.PutUint64(ivBytes[:], iv)
			buf.Write(ivBytes[:])
		}
		// n_no_coor (uint64) per CSI spec.
		var ncBytes [8]byte
		binary.LittleEndian.PutUint64(ncBytes[:], ref.NoCoor)
		buf.Write(ncBytes[:])
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
	data := makeIndex(14, refs)
	idx, err := NewReader(bytes.NewReader(data)).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if idx.MinShift != 14 {
		t.Errorf("MinShift: %d", idx.MinShift)
	}
	if idx.Depth != 5 {
		t.Errorf("Depth: %d", idx.Depth)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References: %d", len(idx.References))
	}
	if got := idx.WindowSize(); got != 16384 {
		t.Errorf("WindowSize: %d want 16384", got)
	}
}

func TestRead_CustomMinShift(t *testing.T) {
	ref := &Reference{
		Bins:      make(map[uint32][]Chunk),
		Intervals: []VirtualOffset{0x1000},
	}
	ref.Bins[1] = []Chunk{{Beg: 0x1000, End: 0x2000}}
	// min_shift=20 → 1 MiB windows.
	data := makeIndex(20, []*Reference{ref})
	idx, err := NewReader(bytes.NewReader(data)).Read()
	if err != nil {
		t.Fatal(err)
	}
	if got := idx.WindowSize(); got != 1<<20 {
		t.Errorf("WindowSize: %d", got)
	}
}

func TestRead_BadMagic(t *testing.T) {
	data := []byte{'X', 'S', 'I', 0x01, 0, 0, 0, 0}
	_, err := NewReader(bytes.NewReader(data)).Read()
	if err == nil {
		t.Fatal("expected bad-magic error")
	}
	var re *ReadError
	if !errors.As(err, &re) {
		t.Errorf("error type: %v", err)
	}
}

func TestRead_Truncated(t *testing.T) {
	_, err := NewReader(bytes.NewReader([]byte{'C', 'S', 'I'})).Read()
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

func TestRead_MultipleRefs(t *testing.T) {
	mkRef := func(bin uint32, beg, end uint64, ivs []uint64) *Reference {
		ref := &Reference{
			Bins:      make(map[uint32][]Chunk),
			Intervals: asVO(ivs),
		}
		ref.Bins[bin] = []Chunk{{Beg: beg, End: end}}
		return ref
	}
	refs := []*Reference{
		mkRef(1, 0x1000, 0x2000, []uint64{0x1000}),
		mkRef(2, 0x3000, 0x4000, []uint64{0x3000, 0x4000, 0x5000}),
	}
	data := makeIndex(14, refs)
	idx, err := NewReader(bytes.NewReader(data)).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 2 {
		t.Fatalf("References: %d", len(idx.References))
	}
}

func TestReg2Bins(t *testing.T) {
	// Single base: 5 bins (one per level).
	bins := reg2bins(100, 101, 14, 5)
	// 6-level scheme: one bin per level, 6 bins for a 1-bp region.
	if len(bins) != 6 {
		t.Errorf("reg2bins(100, 101): got %d bins, want 6", len(bins))
	}
	if bins[0] != 0 {
		t.Errorf("reg2bins(100, 101): first bin %d, want 0", bins[0])
	}
	if bins := reg2bins(100, 100, 14, 5); len(bins) != 0 {
		t.Errorf("reg2bins(100, 100): got %v", bins)
	}
	bins = reg2bins(0, 100_000, 14, 5)
	if len(bins) == 0 {
		t.Error("reg2bins(0, 100000): got empty")
	}
	// samtools 6-level scheme includes 16 Kb bins (>= 4681) for a
	// 100 Kb region. Require that we emit at least one such bin
	// and never exceed BinLimit.
	var saw16K bool
	for _, b := range bins {
		if b >= 4681 {
			saw16K = true
		}
		if b > 37448 {
			t.Errorf("reg2bins(0, 100000): bin %d > 37448 (bin_limit)", b)
		}
	}
	if !saw16K {
		t.Error("reg2bins(0, 100000): no 16 Kb-level bin emitted (>= 4681)")
	}
}

func TestBinFromRange(t *testing.T) {
	// Smoke-test the public API.
	if len(BinFromRange(1000, 5000)) == 0 {
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
	all := ref.LookupChunks(nil)
	if len(all) != 2 {
		t.Errorf("all chunks: got %d", len(all))
	}
	some := ref.LookupChunks([]uint32{1})
	if len(some) != 1 || some[0].Beg != 0x1000 {
		t.Errorf("bin 1 chunks: %+v", some)
	}
}

func TestReadError_Error(t *testing.T) {
	cause := errors.New("bad")
	e := &ReadError{Err: cause}
	if msg := e.Error(); !contains(msg, "csi") || !contains(msg, "bad") {
		t.Errorf("Error: %q", msg)
	}
	if e.Unwrap() != cause { //nolint:errorlint
		t.Errorf("Unwrap: %v, want %v", e.Unwrap(), cause)
	}
}

func contains(s, sub string) bool {
	for i := 0; i+len(sub) <= len(s); i++ {
		if s[i:i+len(sub)] == sub {
			return true
		}
	}
	return false
}

func asVO(in []uint64) []VirtualOffset {
	out := make([]VirtualOffset, len(in))
	copy(out, in)
	return out
}
