package tabix

import (
	"bytes"
	"encoding/binary"
	"testing"
)

// A real Tabix header (per hts-specs tabix.tex): magic +
// n_ref + format + col_seq + col_beg + col_end + meta + skip
// + l_nm + names, then per-reference data. For a VCF index
// the header is typically format=2, col_seq=1, col_beg=2,
// col_end=0, meta=35 ('#'), skip=0.
func buildRealShapeTabix(t *testing.T) []byte {
	t.Helper()
	var buf bytes.Buffer
	// bytes.Buffer.Write never errors, but errcheck sees the
	// pattern. The local helpers discard the errcheck warning.
	w32 := func(v uint32) { _ = binary.Write(&buf, binary.LittleEndian, v) }
	w32s := func(v int32) { _ = binary.Write(&buf, binary.LittleEndian, v) }
	w64 := func(v uint64) { _ = binary.Write(&buf, binary.LittleEndian, v) }
	buf.Write(Magic[:])
	w32(1)   // n_ref
	w32(2)   // format = VCF
	w32(1)   // col_seq
	w32(2)   // col_beg
	w32(0)   // col_end
	w32s(35) // meta = '#'
	w32s(0)  // skip
	// l_nm + names (chr1\0)
	name := []byte{'c', 'h', 'r', '1', 0}
	w32s(int32(len(name)))
	buf.Write(name)
	// reference 0: bin 0, 1 chunk, 0 intervals, n_no_coor=0
	w32(1)   // n_bin
	w32(0)   // bin 0
	w32(1)   // n_chunk
	w64(100) // chunk.beg
	w64(200) // chunk.end
	w32(0)   // n_intv
	w64(0)   // n_no_coor
	return buf.Bytes()
}

func TestRead_RealSpecShapedTabixFile(t *testing.T) {
	data := buildRealShapeTabix(t)
	idx, err := NewReader(bytes.NewReader(data)).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if got := idx.Format; got != FormatVCF {
		t.Errorf("Format = %d, want FormatVCF (%d)", got, FormatVCF)
	}
	if got := idx.ColSeq; got != 1 {
		t.Errorf("ColSeq = %d, want 1", got)
	}
	if got := idx.MetaChar; got != '#' {
		t.Errorf("MetaChar = %d, want '#' (%d)", got, '#')
	}
	if got := idx.Names; len(got) != 1 || got[0] != "chr1" {
		t.Errorf("Names = %v, want [chr1]", got)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References = %d, want 1", len(idx.References))
	}
	if got := idx.References[0].Name; got != "chr1" {
		t.Errorf("Reference[0].Name = %q, want chr1", got)
	}
	if got := idx.References[0].Bins[0]; len(got) != 1 || got[0].Beg != 100 || got[0].End != 200 {
		t.Errorf("Reference[0].Bins[0] = %+v", got)
	}
}

func TestReg2Bins_EmitsAllSixLevels(t *testing.T) {
	// A 1-bp region must produce exactly 6 bins (one per
	// level) for the canonical 14/5 config: 0, 1, 9, 73, 585,
	// 4681.
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
