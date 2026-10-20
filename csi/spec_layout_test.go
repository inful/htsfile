package csi

import (
	"bytes"
	"encoding/binary"
	"testing"
)

// A real CSI v1 header (per hts-specs CSIv1.tex):
//
//	magic      char[4]  "CSI\1"
//	min_shift  int32
//	depth      int32
//	l_aux      int32
//	aux        uint8[l_aux]
//	n_ref      int32
//	ref_i:
//	  n_bin      int32
//	  bin_i:
//	    bin      uint32
//	    loffset  uint64
//	    n_chunk  int32
//	      chunk_beg uint64, chunk_end uint64
//	  n_intv     int32
//	    ioffset  uint64 * n_intv
//	  n_no_coor  uint64
func buildRealShapeCSI(t *testing.T) []byte {
	t.Helper()
	var buf bytes.Buffer
	// bytes.Buffer.Write never errors; the helpers discard
	// the errcheck warning.
	w32 := func(v uint32) { _ = binary.Write(&buf, binary.LittleEndian, v) }
	w32s := func(v int32) { _ = binary.Write(&buf, binary.LittleEndian, v) }
	w64 := func(v uint64) { _ = binary.Write(&buf, binary.LittleEndian, v) }
	buf.Write(Magic[:])
	w32s(14) // min_shift
	w32s(5)  // depth
	w32s(0)  // l_aux = 0
	w32(1)   // n_ref
	// reference 0: 1 bin (bin 0), 1 chunk, 0 intervals, n_no_coor=0
	w32(1)   // n_bin
	w32(0)   // bin
	w64(0)   // loffset
	w32(1)   // n_chunk
	w64(100) // chunk.beg
	w64(200) // chunk.end
	w32(0)   // n_intv
	w64(0)   // n_no_coor
	return buf.Bytes()
}

func TestRead_RealSpecShapedCSIFile(t *testing.T) {
	data := buildRealShapeCSI(t)
	idx, err := NewReader(bytes.NewReader(data)).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if got := idx.MinShift; got != 14 {
		t.Errorf("MinShift = %d, want 14", got)
	}
	if got := idx.Depth; got != 5 {
		t.Errorf("Depth = %d, want 5", got)
	}
	if got := idx.BinLimit(); got != 37449 {
		t.Errorf("BinLimit = %d, want 37449", got)
	}
	if len(idx.Aux) != 0 {
		t.Errorf("Aux = %v, want empty", idx.Aux)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References = %d, want 1", len(idx.References))
	}
	if got := idx.References[0].Bins[0]; len(got) != 1 || got[0].Beg != 100 || got[0].End != 200 {
		t.Errorf("Bins[0] = %+v", got)
	}
}

func TestReg2Bins_EmitsAllSixLevels(t *testing.T) {
	// 1-bp region, 6 levels at min_shift=14, depth=5.
	want := []uint32{0, 1, 9, 73, 585, 4681}
	got := BinFromRange(100, 101)
	if len(got) != len(want) {
		t.Fatalf("BinFromRange(100, 101) = %v, want %v", got, want)
	}
	for i, w := range want {
		if got[i] != w {
			t.Errorf("BinFromRange[%d] = %d, want %d", i, got[i], w)
		}
	}
}
