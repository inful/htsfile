package bgzf

import (
	"bufio"
	"io"
)

// gzipMagic is the two-byte magic every gzip and BGZF stream
// starts with.
var gzipMagic = [2]byte{0x1f, 0x8b}

// InflateIfBGZF returns a reader over the decompressed contents
// of r when r is a gzip/BGZF stream, and r itself otherwise.
//
// It is used by the index readers: htslib writes BGZF-compressed
// .tbi and .csi files, and a raw index file is also accepted.
func InflateIfBGZF(r io.Reader) io.Reader {
	br, ok := r.(*bufio.Reader)
	if !ok {
		br = bufio.NewReader(r)
	}
	var magic [2]byte
	if _, err := io.ReadFull(br, magic[:]); err != nil {
		return r
	}
	// Hand the magic back to the reader in both cases: the
	// caller needs the whole stream, starting with the gzip
	// header.
	p := &prefixReader{r: br}
	copy(p.prefix[:], magic[:])
	if magic != gzipMagic {
		// Not compressed.
		return p
	}
	return NewReader(p)
}

// prefixReader re-reads bytes that were consumed while
// sniffing the stream.
type prefixReader struct {
	prefix [2]byte
	n      int
	r      io.Reader
}

func (p *prefixReader) Read(b []byte) (int, error) {
	if p.n < len(p.prefix) {
		n := copy(b, p.prefix[p.n:])
		p.n += n
		return n, nil
	}
	return p.r.Read(b)
}
