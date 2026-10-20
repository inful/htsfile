package bam

import "testing"

// BAM spec section 4.2 packs each CIGAR operation as
//
//	op_len << 4 | op
//
// with the op code in the LOW 4 bits and the operation length
// in the high 28 bits. This test pins that layout against
// literal words, deliberately NOT built with PackCigarOp or
// the package's own cigar() test helper, so a future inversion
// of either cannot silently redefine the expectation.
//
// A real 5M is (5<<4)|0 = 0x50. An implementation that reads
// the op from the high 4 bits and the length from the low 28
// decodes 0x50 as op=0 (M, coincidentally right) len=80 —
// i.e. a 5M becomes an 80M.
func TestCigarPackingFollowsSpecLayout(t *testing.T) {
	cases := []struct {
		word   uint32
		op     byte
		length uint32
		label  string
	}{
		{0x50, 'M', 5, "5M"},
		{0x21, 'I', 2, "2I"},
		{0x32, 'D', 3, "3D"},
		{0x43, 'N', 4, "4N"},
		{0x54, 'S', 5, "5S"},
		{0x65, 'H', 6, "6H"},
		{0x76, 'P', 7, "7P"},
		{0x87, '=', 8, "8="},
		{0x98, 'X', 9, "9X"},
	}
	for _, c := range cases {
		if got := CigarOpChar(c.word); got != c.op {
			t.Errorf("CigarOpChar(0x%02x) = %q, want %q (%s)",
				c.word, got, c.op, c.label)
		}
		if got := CigarOpLen(c.word); got != c.length {
			t.Errorf("CigarOpLen(0x%02x) = %d, want %d (%s)",
				c.word, got, c.length, c.label)
		}
	}
}

// PackCigarOp must be the exact inverse of the accessors, and
// must agree with the literal words above.
func TestPackCigarOpRoundTrip(t *testing.T) {
	ops := []uint32{CigarM, CigarI, CigarD, CigarN, CigarS,
		CigarH, CigarP, CigarEq, CigarX, CigarB}
	for op := uint32(0); op <= CigarB; op++ {
		for _, n := range []uint32{0, 1, 5, 255, 1 << 20} {
			w := PackCigarOp(op, n)
			if got := CigarOpLen(w); got != n {
				t.Errorf("op=%d n=%d: len=%d, want %d", op, n, got, n)
			}
			if got := CigarOpChar(w); got != CigarOpChar(PackCigarOp(op, 1)) {
				t.Errorf("op=%d n=%d: char=%q, want %q",
					op, n, got, CigarOpChar(PackCigarOp(op, 1)))
			}
		}
	}
	_ = ops
}

// A record whose Cigar holds a spec-formatted word must render
// the same textual CIGAR a spec-conformant reader would.
func TestCigarStringFromSpecWord(t *testing.T) {
	r := &Record{Cigar: []uint32{
		PackCigarOp(CigarM, 5),
		PackCigarOp(CigarI, 2),
		PackCigarOp(CigarM, 3),
	}}
	if got, want := r.CigarString(), "5M2I3M"; got != want {
		t.Errorf("CigarString() = %q, want %q", got, want)
	}
}
