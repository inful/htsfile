package bcf

import (
	"bytes"
	"encoding/binary"
	"testing"
)

// TestDecoder_BCF22RoundTrip builds a BCF 2.2 stream and
// verifies the decoder reads it. We hand-build the stream
// (BGZF + Magic22 + l_text + header text) because the
// library's own encoder always emits BCF 2.1.
func TestDecoder_BCF22RoundTrip(t *testing.T) {
	hdrText := "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
	b := newBuilderWithHeader(t, hdrText)
	payload := make([]byte, 0, len(magic22)+4+len(hdrText))
	payload = append(payload, magic22[:]...)
	var lt [4]byte
	binary.LittleEndian.PutUint32(lt[:], uint32(len(hdrText))) //nolint:gosec // test data, length fits
	payload = append(payload, lt[:]...)
	payload = append(payload, []byte(hdrText)...)
	b.writeBGZFBlock(payload)
	data := b.finish(t)

	dec := NewDecoder(bytes.NewReader(data))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	if got := dec.Version(); got != 2 {
		t.Errorf("Version: got %d, want 2 (BCF 2.2)", got)
	}
}