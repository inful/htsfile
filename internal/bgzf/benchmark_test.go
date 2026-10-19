package bgzf

import (
	"bytes"
	"io"
	"testing"
)

// BenchmarkReader_Decompress measures raw BGZF
// decompression throughput on an in-memory stream of
// one block of incompressible data.
func BenchmarkReader_Decompress(b *testing.B) {
	const payloadSize = 64 * 1024
	payload := make([]byte, payloadSize)
	for i := range payload {
		payload[i] = byte(i)
	}
	stream := append(makeBlock(payload), makeEOFFenceMarker()...)
	b.SetBytes(int64(payloadSize))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r := NewReader(bytes.NewReader(stream))
		buf := make([]byte, 4096)
		for {
			_, err := r.Read(buf)
			if err != nil {
				break
			}
		}
	}
	_ = io.EOF
}

// BenchmarkSeekReader_RandomSeek measures the cost of
// repeatedly seeking to a fixed virtual offset (the start
// of the second BGZF block in a two-block stream).
func BenchmarkSeekReader_RandomSeek(b *testing.B) {
	const payloadSize = 64 * 1024
	payload := make([]byte, payloadSize)
	for i := range payload {
		payload[i] = byte(i)
	}
	block1 := makeBlock(payload)
	block2 := makeBlock(payload)
	stream := append(append(block1, block2...), makeEOFFenceMarker()...)

	// Virtual offset of the second block's first byte.
	secondBlockVirt := uint64(len(block1)) << 16

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		rd := &readSeekerBytes{buf: stream}
		sr, err := NewSeekReader(rd)
		if err != nil {
			b.Fatal(err)
		}
		if err := sr.Seek(secondBlockVirt); err != nil {
			b.Fatal(err)
		}
		buf := make([]byte, 1024)
		if _, err := sr.Read(buf); err != nil {
			b.Fatal(err)
		}
		_ = sr.Close()
	}
}

// readSeekerBytes wraps a byte slice as io.ReadSeekCloser
// so SeekReader can mount it.
type readSeekerBytes struct {
	buf []byte
	off int
}

func (r *readSeekerBytes) Read(p []byte) (int, error) {
	if r.off >= len(r.buf) {
		return 0, io.EOF
	}
	n := copy(p, r.buf[r.off:])
	r.off += n
	return n, nil
}

func (r *readSeekerBytes) Seek(offset int64, whence int) (int64, error) {
	var abs int64
	switch whence {
	case io.SeekStart:
		abs = offset
	case io.SeekCurrent:
		abs = int64(r.off) + offset
	case io.SeekEnd:
		abs = int64(len(r.buf)) + offset
	default:
		return 0, io.ErrNoProgress //nolint:goerr113
	}
	if abs < 0 {
		return 0, io.ErrNoProgress //nolint:goerr113
	}
	r.off = int(abs)
	return abs, nil
}

func (r *readSeekerBytes) Close() error { return nil }