package main

import (
	"bytes"
	"strings"
	"testing"

	"github.com/inful/htsfile/bam"
)

// sam2bam must emit CIGAR words in the BAM spec's layout,
// op_len<<4 | op, with the op code in the LOW 4 bits and the
// operation length in the high 28.
//
// This is the end-to-end interop guard. The unit-level round
// trip (SAM -> BAM -> SAM) cannot catch a layout inversion on
// its own, because both the packer and the unpacker were
// wrong in the same direction and agreed with each other. The
// literal words below are what samtools, htsjdk and htslib
// actually put on disk, so they pin the layout to the outside
// world rather than to us.
func TestSamToBamCigarUsesSpecLayout(t *testing.T) {
	src := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:1000\n" +
		"r1\t0\tchr1\t100\t60\t5M2I3M\t*\t0\t0\tACGTACGTA\tIIIIIIIII\tNM:i:0\n" +
		"r2\t0\tchr1\t200\t60\t2D\t*\t0\t0\tAC\tII\tNM:i:0\n"

	var buf bytes.Buffer
	if err := runSAMToBAM(strings.NewReader(src), &buf); err != nil {
		t.Fatalf("sam2bam: %v", err)
	}

	dec := bam.NewDecoder(bytes.NewReader(buf.Bytes()))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("bam header: %v", err)
	}
	if hdr == nil {
		t.Fatal("nil bam header")
	}

	// r1: 5M2I3M -> 0x50, 0x21, 0x33
	want1 := []uint32{5 << 4, 2<<4 | 1, 3 << 4}
	rec1, err := dec.Record()
	if err != nil {
		t.Fatalf("bam record 1: %v", err)
	}
	assertCigar(t, rec1, want1, "5M2I3M")

	// r2: 2D -> 0x22
	rec2, err := dec.Record()
	if err != nil {
		t.Fatalf("bam record 2: %v", err)
	}
	assertCigar(t, rec2, []uint32{2<<4 | 2}, "2D")
}

func assertCigar(t *testing.T, rec *bam.Record, want []uint32, label string) {
	t.Helper()
	if len(rec.Cigar) != len(want) {
		t.Fatalf("%s: got %d cigar words %v, want %d %v",
			label, len(rec.Cigar), rec.Cigar, len(want), want)
	}
	for i := range want {
		if rec.Cigar[i] != want[i] {
			t.Errorf("%s: cigar[%d] = %#08x, want %#08x (CigarString %q)",
				label, i, rec.Cigar[i], want[i], rec.CigarString())
		}
	}
	// And the library's own renderer must agree with the
	// literal words we just asserted.
	if got := rec.CigarString(); got != label {
		t.Errorf("CigarString() = %q, want %q", got, label)
	}
}
