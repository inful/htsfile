package fasta

import (
	"bufio"
	"errors"
	"fmt"
	"io"
)

// Decoder reads a FASTA document incrementally from an
// underlying [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("sequences.fa")
//	defer f.Close()
//	dec := fasta.NewDecoder(f)
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	}
//
// A Decoder is not safe for concurrent use; create one per
// goroutine.
type Decoder struct {
	scanner *bufio.Scanner

	maxLineLength int

	line int
	eof  bool

	// If we've read a header line and parked it before the
	// caller asked for Record, return that line first.
	pendingHeader []byte
	pendingLine   int
	pendingSet    bool
}

// DefaultMaxLineLength is the per-line size cap applied by
// [NewDecoder]. Lines longer than this cause a parse error
// rather than an unbounded allocation; the value is large
// enough for any real-world FASTA sequence line (typical
// wrapping is 60–80 chars) and small enough to prevent OOM
// on hostile input.
const DefaultMaxLineLength = 10 << 20 // 10 MB

// NewDecoder wraps r in a [Decoder] with the default
// [DefaultMaxLineLength] per-line size cap.
func NewDecoder(r io.Reader) *Decoder {
	return NewDecoderWithMaxLine(r, DefaultMaxLineLength)
}

// NewDecoderWithMaxLine is like [NewDecoder] but lets the
// caller set a custom per-line size cap. Pass maxLineLength
// <= 0 to keep the default cap.
func NewDecoderWithMaxLine(r io.Reader, maxLineLength int) *Decoder {
	if maxLineLength <= 0 {
		maxLineLength = DefaultMaxLineLength
	}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, maxLineLength), maxLineLength)
	return &Decoder{scanner: sc, maxLineLength: maxLineLength}
}

// Record returns the next record, or (nil, [io.EOF]) when no
// more records are available. Calling Record after EOF
// continues to return (nil, [io.EOF]).
//
// FASTA's structure has no header section that needs to be
// read separately — records are self-delimiting. The
// streaming pattern (Record before Header / no Header at all)
// matches the other packages' Decoder shapes; there is no
// [Decoder.Header] method.
func (d *Decoder) Record() (*Record, error) {
	var header []byte
	var headerLine int
	if d.pendingSet {
		header = d.pendingHeader
		headerLine = d.pendingLine
		d.pendingHeader = nil
		d.pendingSet = false
	} else {
		hdr, l, err := d.readLine()
		if err != nil {
			return nil, err
		}
		if len(hdr) == 0 || hdr[0] != '>' {
			return nil, &ParseError{Line: l, Field: "header",
				Msg: "FASTA record must start with '>'"}
		}
		header = hdr[1:]
		headerLine = l
	}
	id, desc := parseHeader(string(header))
	rec := &Record{
		Line:        headerLine,
		ID:          id,
		Description: desc,
		Raw:         string(header),
	}
	// Read sequence lines until next '>' or EOF.
	var seq []byte
	for {
		raw, l, err := d.readLine()
		if err != nil {
			if errors.Is(err, io.EOF) {
				rec.Sequence = stripWhitespace(seq)
				return rec, nil
			}
			return nil, err
		}
		if len(raw) == 0 {
			continue
		}
		if raw[0] == '>' {
			// Park for the next Record call (already
			// stripped of the leading ">" so the next
			// Record entry can take the same path as a
			// fresh header).
			d.pendingHeader = raw[1:]
			d.pendingLine = l
			d.pendingSet = true
			rec.Sequence = stripWhitespace(seq)
			return rec, nil
		}
		seq = append(seq, raw...)
	}
}

// readLine reads one line. Trailing newline / CR is stripped.
// The Scanner internally enforces [DefaultMaxLineLength] (or
// whatever the caller set via [NewDecoderWithMaxLine]); lines
// longer than that return a [ParseError] without unbounded
// allocation.
func (d *Decoder) readLine() (raw []byte, lineNo int, err error) {
	if d.eof {
		return nil, d.line, io.EOF
	}
	lineNo = d.line + 1
	if !d.scanner.Scan() {
		if scanErr := d.scanner.Err(); scanErr != nil {
			if errors.Is(scanErr, bufio.ErrTooLong) {
				return nil, lineNo, &ParseError{Line: lineNo,
					Msg: fmt.Sprintf("line exceeds max length %d", d.maxLineLength)}
			}
			return nil, lineNo, &ParseError{Line: lineNo, Msg: scanErr.Error()}
		}
		d.eof = true
		d.line = lineNo
		return nil, d.line, io.EOF
	}
	src := d.scanner.Bytes()
	raw = make([]byte, len(src))
	copy(raw, src)
	if n := len(raw); n > 0 && raw[n-1] == '\r' {
		raw = raw[:n-1]
	}
	d.line = lineNo
	return raw, lineNo, nil
}
