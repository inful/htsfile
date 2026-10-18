package fai

import (
	"errors"
	"fmt"
	"io"
)

// Build reads a FASTA document from r and returns an [*Index]
// describing it, suitable for serialising with [Writer.Write].
// The returned index is in source order.
//
// The input must be a valid FASTA file: records begin with a
// '>' header line; sequence data follows on subsequent lines
// until the next '>' or EOF. Blank lines between records are
// tolerated.
//
// The index captures, per record:
//
//   - Name      — the first whitespace-delimited token of the
//     header.
//   - Length    — total number of sequence bases (whitespace
//     stripped).
//   - Offset    — byte offset of the first base in the input.
//   - LineBases — number of bases on the first sequence line
//     of the record.
//   - LineWidth — byte width of that first sequence line
//     including its trailing newline (lineBases + 1 for LF,
//     lineBases + 2 for CR-LF).
//
// Build assumes a uniform line width throughout each record
// (the standard FASTA convention). Mixed-width records
// produce an index whose seek math is undefined downstream.
//
// Build never buffers the whole FASTA file: it reads one line
// at a time, so multi-gigabyte inputs work in O(1) memory
// beyond the index itself.
//
// Implementation note: Build reads through a countingReader
// byte-by-byte rather than via bufio.Reader. A bufio.Reader
// prefetches in chunks, which makes the byte count of the
// underlying reader diverge from the line offset the parser
// is interested in. Byte-at-a-time is O(N) syscalls but
// correct; index generation is a one-shot operation.
func Build(r io.Reader) (*Index, error) {
	idx := &Index{}
	cr := &countingReader{r: r}
	lr := &lineReader{r: cr}

	var cur *Reference
	var bases int64
	var sawHeader bool
	var firstBaseSeen bool

	flush := func() {
		if cur == nil {
			return
		}
		cur.Length = bases
		idx.References = append(idx.References, cur)
	}

	for {
		line, err := lr.readLine()
		if len(line) > 0 {
			// Determine whether the line ended in LF or CR-LF.
			n := len(line)
			hasLF := n > 0 && line[n-1] == '\n'
			hasCR := n > 1 && line[n-2] == '\r'
			contentEnd := n
			if hasLF {
				contentEnd--
			}
			if hasCR {
				contentEnd--
			}
			content := line[:contentEnd]

			if len(content) == 0 {
				goto next
			}

			if content[0] == '>' {
				flush()
				cur = &Reference{Name: parseHeaderName(string(content[1:]))}
				bases = 0
				firstBaseSeen = false
				sawHeader = true
				goto next
			}

			if cur == nil {
				if !sawHeader {
					return nil, errMissingGT
				}
				goto next
			}
			if !firstBaseSeen {
				// cr.n now sits at the byte just past
				// the line's trailing newline. The
				// first byte of this line is therefore
				// cr.n - n.
				cur.Offset = cr.n - int64(n)
				cur.LineBases = countBases(content)
				cur.LineWidth = n
				firstBaseSeen = true
			}
			bases += int64(countBases(content))
		}
	next:
		if err != nil {
			if errors.Is(err, io.EOF) {
				flush()
				if !sawHeader {
					return nil, errMissingGT
				}
				return idx, nil
			}
			return nil, fmt.Errorf("fai: read: %w", err)
		}
	}
}

// parseHeaderName returns the first whitespace-delimited token
// of a FASTA header (the "ID" portion).
func parseHeaderName(header string) string {
	for i := 0; i < len(header); i++ {
		if header[i] == ' ' || header[i] == '\t' {
			return header[:i]
		}
	}
	return header
}

// countBases returns the number of non-whitespace bytes in b
// (the FASTA spec allows spaces and tabs inside wrapped
// sequence lines; the .fai Length field counts bases only).
func countBases(b []byte) int {
	n := 0
	for _, c := range b {
		if c == ' ' || c == '\t' || c == '\r' || c == '\n' {
			continue
		}
		n++
	}
	return n
}

// countingReader is an io.Reader wrapper that tallies the
// total number of bytes consumed. We need the offset of the
// first base of each record, and the FASTA parser does not
// expose one.
type countingReader struct {
	r io.Reader
	n int64
}

func (c *countingReader) Read(p []byte) (int, error) {
	n, err := c.r.Read(p)
	c.n += int64(n)
	return n, err
}

// lineReader returns one line at a time (including the
// trailing newline) by reading byte-by-byte from the wrapped
// reader. It exists to keep the line offset (relative to the
// start of the stream) trivially computable from the wrapped
// reader's byte count.
type lineReader struct {
	r io.Reader
	buf []byte
}

func (l *lineReader) readLine() ([]byte, error) {
	l.buf = l.buf[:0]
	var b [1]byte
	for {
		n, err := l.r.Read(b[:])
		if n > 0 {
			l.buf = append(l.buf, b[0])
			if b[0] == '\n' {
				return l.buf, nil
			}
			continue
		}
		if err != nil {
			if len(l.buf) > 0 {
				// Last line without trailing newline;
				// surface the buffered bytes alongside
				// the EOF so the caller still sees the
				// final line.
				return l.buf, err
			}
			return nil, err
		}
	}
}
