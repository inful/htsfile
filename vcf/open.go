package vcf

import (
	"fmt"
	"io"
	"os"

	"github.com/inful/htsfile/internal/bgzf"
)

// Open returns a streaming [*Decoder] reading from path.
// The path may be a plain .vcf file or a BGZF-compressed
// .vcf.gz file. The format is auto-detected: the BGZF
// reader peeks at the first two bytes and matches the
// gzip magic (1F 8B) before wrapping the reader; plain
// files are read directly.
//
// Calling code should always close the returned decoder
// when finished (it owns no resources but tests and
// callers expect symmetry).
//
// Example:
//
//	dec, err := vcf.Open("sample.vcf.gz")
//	if err != nil { ... }
//	defer /* nothing */ // Decoder has no Close
//	hdr, err := dec.Header()
//	for { rec, err := dec.Record(); ... }
func Open(path string) (*Decoder, error) {
	// path is supplied by the caller; gosec's G304
	// warning about file-inclusion does not apply to
	// a deliberately-passed filename.
	f, err := os.Open(path) //nolint:gosec // path is caller-supplied
	if err != nil {
		return nil, fmt.Errorf("vcf.Open: %w", err)
	}
	return openReader(f), nil
}

// openReader wraps r in a BGZF decoder iff the first two
// bytes are the gzip magic (1F 8B). Plain VCF files pass
// through untouched.
func openReader(r io.Reader) *Decoder {
	return NewDecoder(autoDetectBGZF(r))
}

// autoDetectBGZF peeks at the first two bytes and returns
// an io.Reader that produces the BGZF-decoded stream if
// the magic matches, or the original stream otherwise.
func autoDetectBGZF(r io.Reader) io.Reader {
	// We need to peek without consuming. bufio.Reader.Peek
	// only works for n bytes already buffered; for the
	// common case the first 2 bytes may be at the very
	// start of the file, so we read into a small buffer
	// and tee.
	const magic = "\x1f\x8b"
	pr := newPeekReader(r, 2)
	head, _ := pr.Peek(2)
	if len(head) >= 2 && string(head[:2]) == magic {
		return bgzf.NewReader(pr)
	}
	return pr
}

// peekReader wraps an io.Reader and exposes a Peek
// method backed by a small buffer. After Peek returns
// bytes, subsequent Read calls return those bytes before
// pulling from the underlying reader.
type peekReader struct {
	r   io.Reader
	buf []byte
	off int
}

func newPeekReader(r io.Reader, n int) *peekReader {
	return &peekReader{r: r, buf: make([]byte, 0, n)}
}

// Peek returns the next n bytes without advancing the
// reader. Returns fewer than n bytes if EOF is hit.
func (p *peekReader) Peek(n int) ([]byte, error) {
	for len(p.buf)-p.off < n {
		var b [1]byte
		nr, err := p.r.Read(b[:])
		if nr > 0 {
			p.buf = append(p.buf, b[0])
		}
		if err != nil {
			if err == io.EOF {
				break
			}
			return p.buf[p.off:], err
		}
	}
	return p.buf[p.off : p.off+n], nil
}

func (p *peekReader) Read(dst []byte) (int, error) {
	if p.off < len(p.buf) {
		n := copy(dst, p.buf[p.off:])
		p.off += n
		return n, nil
	}
	return p.r.Read(dst)
}
