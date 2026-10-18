package bgzf

import (
	"bytes"
	"compress/flate"
	"encoding/binary"
	"errors"
	"hash/crc32"
	"io"
	"testing"
)

// bgzfMagic is the gzip magic + deflate + FEXTRA flag header prefix.
var bgzfMagic = []byte{0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00}

// makeBlock returns a single BGZF block containing payload, with
// the correct BSIZE, CRC32, and ISIZE trailer.
func makeBlock(payload []byte) []byte {
	// CDATA = raw deflate of payload.
	var cdata bytes.Buffer
	w, _ := flate.NewWriter(&cdata, flate.DefaultCompression)
	_, _ = w.Write(payload)
	_ = w.Close()
	cd := cdata.Bytes()
	// Block size = 18 (header) + len(cd) + 4 (CRC32) + 4 (ISIZE) = 26 + len(cd).
	bsize := uint16(26 + len(cd) - 1) // BSIZE = total - 1
	buf := bytes.Buffer{}
	buf.Write(bgzfMagic[:16])
	_ = binary.Write(&buf, binary.LittleEndian, bsize)
	buf.Write(cd)
	_ = binary.Write(&buf, binary.LittleEndian, crc32.ChecksumIEEE(payload))
	_ = binary.Write(&buf, binary.LittleEndian, uint32(len(payload)))
	return buf.Bytes()
}

// makeEOFFenceMarker returns the BGZF EOF marker block (BSIZE = 27,
// CDATA = empty-deflate-stream, CRC32 = 0, ISIZE = 0).
func makeEOFFenceMarker() []byte {
	// The standard samtools EOF marker.
	return []byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00,
	}
}

func TestReader_SingleBlock(t *testing.T) {
	want := []byte("hello, BGZF world")
	stream := makeBlock(want)
	r := NewReader(bytes.NewReader(stream))

	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("payload mismatch:\n got=%q\nwant=%q", got, want)
	}
}

func TestReader_MultipleBlocks(t *testing.T) {
	chunks := [][]byte{
		[]byte("first "),
		[]byte("middle "),
		[]byte("third "),
		[]byte("fourth"),
	}
	var stream bytes.Buffer
	for _, c := range chunks {
		stream.Write(makeBlock(c))
	}
	r := NewReader(&stream)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	want := []byte("first middle third fourth")
	if !bytes.Equal(got, want) {
		t.Fatalf("payload mismatch:\n got=%q\nwant=%q", got, want)
	}
}

func TestReader_EOFMarker(t *testing.T) {
	want := []byte("only one block")
	stream := bytes.Buffer{}
	stream.Write(makeBlock(want))
	stream.Write(makeEOFFenceMarker())
	r := NewReader(&stream)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("payload mismatch:\n got=%q\nwant=%q", got, want)
	}
}

func TestReader_BadMagic(t *testing.T) {
	bad := []byte{
		0x00, 0x00, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00,
	}
	r := NewReader(bytes.NewReader(bad))
	_, err := io.ReadAll(r)
	if err == nil {
		t.Fatal("expected error on bad magic, got nil")
	}
	if !errors.Is(err, ErrBadMagic) {
		t.Fatalf("expected ErrBadMagic, got %v", err)
	}
}

func TestReader_BadBSIZE(t *testing.T) {
	// A block whose BSIZE claims a CDATA length different from the
	// actual stream length should fail CRC or ISIZE validation.
	bad := []byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00,
		// BSIZE = 100 (claims 75 bytes of CDATA), but we only have 1.
		0x64, 0x00,
		// CDATA: 1 byte (but BSIZE says 75)
		0x00,
		// pad with zeros up to 75+8 bytes total after the BSIZE so
		// the reader can attempt to read.
	}
	for range 82 {
		bad = append(bad, 0x00)
	}
	r := NewReader(bytes.NewReader(bad))
	_, err := io.ReadAll(r)
	if err == nil {
		t.Fatal("expected error on bad BSIZE, got nil")
	}
}

func TestReader_Empty(t *testing.T) {
	r := NewReader(bytes.NewReader(nil))
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(got) != 0 {
		t.Fatalf("expected empty, got %q", got)
	}
}

// TestReader_BinaryPattern ensures the reader handles bytes that
// include nulls and the full 0x00..0xff range (no string-style
// line handling).
func TestReader_BinaryPattern(t *testing.T) {
	want := make([]byte, 1024)
	for i := range want {
		want[i] = byte(i)
	}
	stream := makeBlock(want)
	r := NewReader(bytes.NewReader(stream))
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatal("binary payload mismatch")
	}
}

// TestReader_LargePayload verifies a payload spanning multiple
// 64KB BGZF blocks round-trips.
func TestReader_LargePayload(t *testing.T) {
	want := bytes.Repeat([]byte("ABCD"), 40000) // 160 KB → multiple blocks
	var stream bytes.Buffer
	// Chunk into pieces small enough for a single BGZF block.
	for i := 0; i < len(want); i += 60000 {
		end := i + 60000
		if end > len(want) {
			end = len(want)
		}
		stream.Write(makeBlock(want[i:end]))
	}
	r := NewReader(&stream)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatal("large payload mismatch")
	}
}
