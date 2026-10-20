package bcf

import (
	"reflect"
	"testing"
)

// TestEncodeGT covers the spec's per-allele int encoding
// for GT (§5 line 2463). Each example is from the spec
// table at line 2474.
func TestEncodeGT(t *testing.T) {
	cases := []struct {
		gt   string
		want []uint8
	}{
		{"0/1", []uint8{0x02, 0x04}},
		{"0|1", []uint8{0x03, 0x05}},
		{"./.", []uint8{0x00, 0x00}},
		{"0", []uint8{0x03}}, // single-value haploid → phased
		{"1", []uint8{0x05}}, // single-value haploid → phased
		{"0/1/2", []uint8{0x02, 0x04, 0x06}},
		{"0/1|2", []uint8{0x02, 0x04, 0x07}},
	}
	for _, c := range cases {
		t.Run(c.gt, func(t *testing.T) {
			got, err := encodeGT(c.gt)
			if err != nil {
				t.Fatalf("encodeGT: %v", err)
			}
			if !reflect.DeepEqual(got, c.want) {
				t.Errorf("got %v, want %v", got, c.want)
			}
		})
	}
}

// TestDecodeGT verifies the inverse mapping for a few
// canonical examples. The decoder is lossless on the
// round-tripped forms above.
func TestDecodeGT(t *testing.T) {
	cases := []struct {
		ints []uint8
		want string
	}{
		{[]uint8{0x02, 0x04}, "0/1"},
		{[]uint8{0x03, 0x05}, "0|1"},
		{[]uint8{0x00, 0x00}, "./."},
		{[]uint8{0x03}, "0"},
		{[]uint8{0x05}, "1"},
		{[]uint8{0x02, 0x04, 0x06}, "0/1/2"},
		{[]uint8{0x02, 0x04, 0x07}, "0/1|2"},
		{nil, ""},
	}
	for _, c := range cases {
		got, err := decodeGT(c.ints)
		if err != nil {
			t.Errorf("decodeGT(%v): %v", c.ints, err)
			continue
		}
		if got != c.want {
			t.Errorf("decodeGT(%v): got %q, want %q", c.ints, got, c.want)
		}
	}
}

// TestGTRoundTrip ensures encodeGT ∘ decodeGT is the
// identity for a variety of forms.
func TestGTRoundTrip(t *testing.T) {
	cases := []string{
		"0/0", "0/1", "1/1", "./.", ".",
		"0|1", "0|0", "1|0",
		"0/1/2", "0|1|2", "0/1|2", "0|1/2",
		"0", "1", ".",
	}
	for _, gt := range cases {
		t.Run(gt, func(t *testing.T) {
			enc, err := encodeGT(gt)
			if err != nil {
				t.Fatalf("encodeGT(%q): %v", gt, err)
			}
			// Strip trailing END_OF_VECTOR padding if
			// any, then decode.
			stripped := stripEndOfVectorGT(enc)
			got, err := decodeGT(stripped)
			if err != nil {
				t.Fatalf("decodeGT: %v", err)
			}
			if got != gt {
				t.Errorf("round-trip: got %q, want %q", got, gt)
			}
		})
	}
}

// stripEndOfVectorGT removes trailing 0x81 bytes from a
// GT int vector. Used by TestGTRoundTrip to handle
// mixed-ploidy GTs that the encoder pads out to the
// record's max ploidy with END_OF_VECTOR.
func stripEndOfVectorGT(v []uint8) []uint8 {
	for len(v) > 0 && v[len(v)-1] == bcfEndOfVectorInt8 {
		v = v[:len(v)-1]
	}
	return v
}
