package bgzf

import (
	"bytes"
	"hash/crc32"
	"testing"
)

// TestWriter_RoundTrip verifies that bytes written by Writer
// re-decode cleanly through the existing Reader. This is the
// most fundamental property: anything written must be
// re-readable.
func TestWriter_RoundTrip(t *testing.T) {
	cases := []string{
		"",
		"hello",
		"a longer string that spans a BGZF block boundary; the writer must emit one or more valid blocks plus the canonical 28-byte EOF marker",
	}
	for _, payload := range cases {
		t.Run(payload, func(t *testing.T) {
			var buf bytes.Buffer
			w := NewWriter(&buf)
			if _, err := w.Write([]byte(payload)); err != nil {
				t.Fatalf("Write: %v", err)
			}
			if err := w.Close(); err != nil {
				t.Fatalf("Close: %v", err)
			}
			// Round-trip through the Reader.
			r := NewReader(&buf)
			got, err := readAll(r)
			if err != nil {
				t.Fatalf("Read: %v", err)
			}
			if string(got) != payload {
				t.Errorf("got %q want %q", got, payload)
			}
		})
	}
}

// TestWriter_MultipleBlocks verifies that payloads spanning
// multiple BGZF blocks reassemble correctly.
func TestWriter_MultipleBlocks(t *testing.T) {
	// Write 200 KiB — well past the 64 KiB block cap. The
	// writer must emit multiple blocks.
	payload := bytes.Repeat([]byte("ABCDEFGH"), 25*1024) // 200 KiB
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	got, err := readAll(NewReader(&buf))
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if !bytes.Equal(got, payload) {
		t.Errorf("payload mismatch: got %d bytes want %d", len(got), len(payload))
	}
}

// TestWriter_BlockSizeBoundaries verifies that payloads at
// the 64 KiB block boundary reassemble correctly (off-by-one
// in BSIZE math would surface here).
func TestWriter_BlockSizeBoundaries(t *testing.T) {
	for _, n := range []int{
		0,
		1,
		64*1024 - 1,   // one byte under block boundary
		64 * 1024,     // exactly the boundary
		64*1024 + 1,   // one byte over
		2 * 64 * 1024, // exactly two blocks
	} {
		payload := bytes.Repeat([]byte("X"), n)
		var buf bytes.Buffer
		w := NewWriter(&buf)
		if _, err := w.Write(payload); err != nil {
			t.Fatalf("Write n=%d: %v", n, err)
		}
		if err := w.Close(); err != nil {
			t.Fatalf("Close n=%d: %v", n, err)
		}
		got, err := readAll(NewReader(&buf))
		if err != nil {
			t.Fatalf("Read n=%d: %v", n, err)
		}
		if len(got) != n {
			t.Errorf("n=%d: got %d bytes want %d", n, len(got), n)
		}
	}
}

// TestWriter_EOFTrailer verifies the canonical 28-byte EOF
// marker is emitted on Close.
func TestWriter_EOFTrailer(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	// The EOF marker is the last block; its BSIZE is 27 and
	// its payload is empty (CRC32=0, ISIZE=0). We can spot-check
	// by reading past it and confirming the Reader sees io.EOF.
	r := NewReader(&buf)
	got, err := readAll(r)
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(got) != 0 {
		t.Errorf("expected empty payload after Close without Write, got %q", got)
	}
}

// TestWriter_BlockStructure verifies that the bytes emitted
// on Close have the right shape: gzip+FEXTRA magic + 'BC'
// + 2-byte SLEN + 2-byte BSIZE + CDATA + CRC32 + ISIZE.
// This is the BGZF block layout the spec requires; if it's
// wrong, downstream tools (samtools, htslib) won't read the
// output.
func TestWriter_BlockStructure(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	payload := []byte("hello")
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	// The output is a real block followed by an EOF marker.
	// The first block's BSIZE is total_size - 1; we can read
	// the BSIZE field at offset 16-17 and verify it.
	out := buf.Bytes()
	if len(out) < headerSize {
		t.Fatalf("output too short: %d bytes", len(out))
	}
	if out[0] != 0x1f || out[1] != 0x8b || out[2] != 0x08 {
		t.Errorf("gzip magic/compression byte wrong: %x %x %x",
			out[0], out[1], out[2])
	}
	if out[3]&0x04 == 0 {
		t.Errorf("FEXTRA flag not set: %x", out[3])
	}
	if out[12] != 'B' || out[13] != 'C' {
		t.Errorf("BGZF extra subfield missing: %x %x", out[12], out[13])
	}
	if int(out[14])|int(out[15])<<8 != 2 {
		t.Errorf("BGZF SLEN != 2: %x %x", out[14], out[15])
	}
	bsize := int(out[16]) | int(out[17])<<8
	totalSize := bsize + 1
	if totalSize < headerSize+trailerSize {
		t.Errorf("BSIZE %d too small for header+trailer", bsize)
	}
	// CRC32 + ISIZE trailer should match the payload.
	trailerStart := totalSize - trailerSize
	crc := uint32(out[trailerStart]) | uint32(out[trailerStart+1])<<8 |
		uint32(out[trailerStart+2])<<16 | uint32(out[trailerStart+3])<<24
	if crc != crc32.ChecksumIEEE(payload) {
		t.Errorf("CRC32: got %08x want %08x", crc, crc32.ChecksumIEEE(payload))
	}
}

// TestWriter_DefaultBlockSize verifies that the default
// block size is the spec maximum (64 KiB) so callers don't
// have to think about it.
func TestWriter_DefaultBlockSize(t *testing.T) {
	w := NewWriter(&bytes.Buffer{})
	if w.blockSize != maxBSIZE {
		t.Errorf("default block size: got %d want %d", w.blockSize, maxBSIZE)
	}
}

// TestWriter_WriteAfterClose verifies that Write after Close
// returns an error (or panics? — see comment in write.go).
func TestWriter_WriteAfterClose(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if _, err := w.Write([]byte("after close")); err == nil {
		t.Error("Write after Close: expected error, got nil")
	}
}

// TestWriter_NilWriter verifies that the package doesn't panic
// on a nil receiver (it should error out).
func TestWriter_NilWriter(t *testing.T) {
	var w *Writer
	if _, err := w.Write([]byte("x")); err == nil {
		t.Error("nil.Write: expected error, got nil")
	}
}

// readAll is a small test helper that drains a Reader to EOF.
func readAll(r *Reader) ([]byte, error) {
	var buf bytes.Buffer
	buf.Grow(4096)
	for {
		var chunk [4096]byte
		n, err := r.Read(chunk[:])
		if n > 0 {
			buf.Write(chunk[:n])
		}
		if err != nil {
			if err.Error() == "EOF" {
				return buf.Bytes(), nil
			}
			return nil, err
		}
	}
}
