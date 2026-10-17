package bai

import (
	"bytes"
	"encoding/binary"
	"errors"
	"io"
	"testing"
)

// makeBAI builds a synthetic BAI file with the given bins and
// intervals and returns the byte slice.
func makeBAI(refs []*Reference, withNoCoor bool) []byte {
	var buf bytes.Buffer
	buf.Write(Magic[:])
	var nr [4]byte
	binary.LittleEndian.PutUint32(nr[:], uint32(len(refs)))
	buf.Write(nr[:])
	for _, ref := range refs {
		var nb [4]byte
		binary.LittleEndian.PutUint32(nb[:], uint32(len(ref.Bins)))
		buf.Write(nb[:])
		for _, b := range ref.Bins {
			var bn [4]byte
			binary.LittleEndian.PutUint32(bn[:], b.Number)
			buf.Write(bn[:])
			var lo [8]byte
			binary.LittleEndian.PutUint64(lo[:], b.Lo)
			buf.Write(lo[:])
			var hi [8]byte
			binary.LittleEndian.PutUint64(hi[:], b.Hi)
			buf.Write(hi[:])
		}
		var ni [4]byte
		binary.LittleEndian.PutUint32(ni[:], uint32(len(ref.Intervals)))
		buf.Write(ni[:])
		for _, iv := range ref.Intervals {
			var ib [8]byte
			binary.LittleEndian.PutUint64(ib[:], iv)
			buf.Write(ib[:])
		}
		if withNoCoor {
			var nc [8]byte
			binary.LittleEndian.PutUint64(nc[:], ref.NoCoor)
			buf.Write(nc[:])
		}
	}
	return buf.Bytes()
}

func TestRead_Basic(t *testing.T) {
	refs := []*Reference{
		{
			Bins: []*Bin{
				{Number: 4681, Lo: 100, Hi: 500},
				{Number: 4682, Lo: 500, Hi: 1000},
			},
			Intervals: []uint64{10, 20, 30, 40},
			NoCoor:    5,
		},
	}
	data := makeBAI(refs, true)
	idx, err := NewReader(bytes.NewReader(data)).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References: %d", len(idx.References))
	}
	ref := idx.References[0]
	if len(ref.Bins) != 2 {
		t.Errorf("Bins: %d", len(ref.Bins))
	}
	if ref.Bins[0].Number != 4681 || ref.Bins[0].Lo != 100 || ref.Bins[0].Hi != 500 {
		t.Errorf("Bin[0]: %+v", ref.Bins[0])
	}
	if len(ref.Intervals) != 4 {
		t.Errorf("Intervals: %d", len(ref.Intervals))
	}
	if ref.NoCoor != 5 {
		t.Errorf("NoCoor: %d", ref.NoCoor)
	}
}

func TestRead_NoNoCoor(t *testing.T) {
	refs := []*Reference{
		{
			Bins:      []*Bin{{Number: 4681, Lo: 100, Hi: 500}},
			Intervals: []uint64{10},
		},
	}
	data := makeBAI(refs, false) // no n_no_coor field
	idx, err := NewReader(bytes.NewReader(data)).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if idx.References[0].NoCoor != 0 {
		t.Errorf("NoCoor: %d (want 0 when field absent)", idx.References[0].NoCoor)
	}
}

func TestRead_BadMagic(t *testing.T) {
	data := []byte{'B', 'A', 'I', 0x02, 0, 0, 0, 0}
	_, err := NewReader(bytes.NewReader(data)).Read()
	if err == nil {
		t.Fatal("expected error on bad magic")
	}
}

func TestRead_TruncatedMagic(t *testing.T) {
	_, err := NewReader(bytes.NewReader([]byte{'B', 'A'})).Read()
	if !errors.Is(err, io.ErrUnexpectedEOF) {
		t.Errorf("expected ErrUnexpectedEOF, got %v", err)
	}
}

func TestRead_Empty(t *testing.T) {
	// n_ref = 0
	refs := []*Reference{}
	data := makeBAI(refs, false)
	idx, err := NewReader(bytes.NewReader(data)).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 0 {
		t.Errorf("References: %d", len(idx.References))
	}
}

func TestReg2Bins(t *testing.T) {
	tests := []struct {
		beg, end int64
		want     []uint32
	}{
		// Single-position region hits bin 4681 (level 5).
		{100, 101, []uint32{0, 1, 9, 73, 585, 4681}},
		// Empty region yields nil.
		{100, 100, nil},
		// 256 KiB range: spans 2 level-4 bins (2^17 each)
		// and 16 level-5 bins (2^14 each).
		{0, 1 << 18, append(
			[]uint32{0, 1, 9, 73, 585, 586},
			range16(4681, 4697)...,
		)},
	}
	for _, tc := range tests {
		got := Reg2Bins(tc.beg, tc.end)
		if !equalBins(got, tc.want) {
			t.Errorf("Reg2Bins(%d, %d):\n  got %v\n want %v", tc.beg, tc.end, got, tc.want)
		}
	}
}

// range16 returns the half-open integer range [lo, hi).
func range16(lo, hi uint32) []uint32 {
	out := make([]uint32, 0, hi-lo)
	for i := lo; i < hi; i++ {
		out = append(out, i)
	}
	return out
}

func equalBins(a, b []uint32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}