package gff

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"strconv"
)

// Decoder reads a GFF3 document incrementally from an
// underlying [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("annotations.gff")
//	defer f.Close()
//	dec := gff.NewDecoder(f)
//	hdr, _ := dec.Header()
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

	// Header state.
	hdrRead bool
	hdr     *Header
	hdrErr  error

	// If the scanner is past the header and parked a data
	// line before the caller asked for Record, return that
	// line first.
	pending     []byte
	pendingLine int
	pendingSet  bool

	// If we've seen `###`, the feature section is over;
	// subsequent Record calls return io.EOF.
	terminated bool
}

// DefaultMaxLineLength is the per-line size cap applied by
// [NewDecoder]. Lines longer than this cause a parse error
// rather than an unbounded allocation; the value is large
// enough for any real-world GFF3 line (attributes can be
// long but rarely exceed a few KiB) and small enough to
// prevent OOM on hostile input.
const DefaultMaxLineLength = 10 << 20 // 10 MB

// NewDecoder wraps r in a [Decoder] with the default
// [DefaultMaxLineLength] per-line size cap. Use
// [NewDecoderWithMaxLine] to override the cap.
func NewDecoder(r io.Reader) *Decoder {
	return NewDecoderWithMaxLine(r, DefaultMaxLineLength)
}

// NewDecoderWithMaxLine is like [NewDecoder] but lets the
// caller set a custom per-line size cap. Pass maxLineLength <= 0
// to keep the default cap.
func NewDecoderWithMaxLine(r io.Reader, maxLineLength int) *Decoder {
	if maxLineLength <= 0 {
		maxLineLength = DefaultMaxLineLength
	}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, maxLineLength), maxLineLength)
	return &Decoder{scanner: sc, maxLineLength: maxLineLength}
}

// Header returns the file's header, or (nil, nil) if the file
// has no header pragmas. The first call reads through blank
// lines and bare `#` comments until it finds the first `##`
// pragma or the first feature line (signalling "no header").
// The result is cached; subsequent calls return the same value.
//
// If Record is called before Header, Header is invoked
// implicitly so the header/data boundary is correctly handled.
func (d *Decoder) Header() (*Header, error) {
	if d.hdrRead {
		return d.hdr, d.hdrErr
	}
	d.hdrRead = true
	hdr := &Header{}
	for {
		raw, lineNo, err := d.readLine()
		if err != nil {
			if errors.Is(err, io.EOF) {
				// Empty file with no header.
				d.hdr = hdr
				return hdr, nil
			}
			d.hdr = hdr
			d.hdrErr = err
			return hdr, err
		}
		if isBlank(raw) {
			continue
		}
		if raw[0] != '#' {
			// Data line; park for Record and report empty header.
			d.pending = raw
			d.pendingLine = lineNo
			d.pendingSet = true
			d.hdr = hdr
			return hdr, nil
		}
		if len(raw) == 1 || raw[1] != '#' {
			// Bare `#` comment line — preserve in Comments.
			hdr.Comments = append(hdr.Comments, string(raw))
			continue
		}
		// `##` pragma.
		if err := parsePragma(hdr, string(raw[2:]), lineNo); err != nil {
			d.hdr = hdr
			d.hdrErr = err
			return hdr, err
		}
		// We've seen the first `##` pragma; we now switch to
		// reading features. The next non-pragma, non-comment
		// line is the first record.
		break
	}
	// Continue reading: blanks/comments/pragmas + first data.
	for {
		raw, lineNo, err := d.readLine()
		if err != nil {
			if errors.Is(err, io.EOF) {
				d.hdr = hdr
				return hdr, nil
			}
			d.hdr = hdr
			d.hdrErr = err
			return hdr, err
		}
		if isBlank(raw) {
			continue
		}
		if raw[0] == '#' {
			if len(raw) == 1 || raw[1] != '#' {
				hdr.Comments = append(hdr.Comments, string(raw))
				continue
			}
			// `##` after the first pragma is unusual but
			// allowed (some producers emit the spec line
			// last).
			if err := parsePragma(hdr, string(raw[2:]), lineNo); err != nil {
				d.hdr = hdr
				d.hdrErr = err
				return hdr, err
			}
			continue
		}
		// First feature line; park it for Record.
		d.pending = raw
		d.pendingLine = lineNo
		d.pendingSet = true
		d.hdr = hdr
		return hdr, nil
	}
}

// Record returns the next feature record. When no more records
// are available it returns (nil, [io.EOF]). Calling Record
// after EOF continues to return (nil, [io.EOF]).
//
// If Header has not yet been called, Record calls it implicitly.
// A `###` line terminates the feature section; subsequent Record
// calls return io.EOF.
func (d *Decoder) Record() (*Record, error) {
	if !d.hdrRead {
		if _, err := d.Header(); err != nil {
			return nil, err
		}
	}
	if d.terminated {
		// The `###` line ended the feature section; any
		// subsequent lines (in v1, just whitespace; FASTA
		// after `###` is out of scope) are ignored.
		return nil, io.EOF
	}
	for {
		var raw []byte
		var lineNo int
		if d.pendingSet {
			raw = d.pending
			lineNo = d.pendingLine
			d.pending = nil
			d.pendingSet = false
		} else {
			r, l, err := d.readLine()
			if err != nil {
				return nil, err
			}
			raw = r
			lineNo = l
		}
		if isBlank(raw) {
			continue
		}
		// `###` record terminator: end the feature section.
		if bytes.Equal(raw, []byte("###")) {
			d.terminated = true
			return nil, io.EOF
		}
		// Bare `#` comments and `##` pragmas are silently
		// skipped in the data section (the spec only
		// requires pragmas in the header; some tools emit
		// tool-specific pragmas mid-file).
		if raw[0] == '#' {
			continue
		}
		return parseRecord(raw, lineNo)
	}
}

// readLine reads one line. Trailing newline / CR is stripped.
// The Scanner internally enforces the per-line cap; lines longer
// than that return a [ParseError] without unbounded allocation.
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

func isBlank(b []byte) bool {
	for _, c := range b {
		if c != ' ' && c != '\t' {
			return false
		}
	}
	return true
}

// parseRecord parses one tab-separated feature line.
func parseRecord(raw []byte, lineNo int) (*Record, error) {
	fields := bytes.Split(raw, []byte{'\t'})
	if len(fields) != 9 {
		return nil, &ParseError{Line: lineNo, Field: "record",
			Msg: "expected 9 tab-separated columns, got " + strconv.Itoa(len(fields))}
	}
	rec := &Record{
		Line:  lineNo,
		SeqID: string(fields[0]),
		// Source: "." means "not available"; preserve as-is.
		Source: string(fields[1]),
		Type:   string(fields[2]),
		Raw:    string(raw),
	}
	// start / end
	start, err := strconv.ParseInt(string(fields[3]), 10, 64)
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "start",
			Msg: "non-integer: " + string(fields[3])}
	}
	end, err := strconv.ParseInt(string(fields[4]), 10, 64)
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "end",
			Msg: "non-integer: " + string(fields[4])}
	}
	if start < 1 {
		return nil, &ParseError{Line: lineNo, Field: "start",
			Msg: "GFF3 coordinates are 1-based; got " + strconv.FormatInt(start, 10)}
	}
	if end < start {
		return nil, &ParseError{Line: lineNo, Field: "end",
			Msg: fmt.Sprintf("end %d < start %d", end, start)}
	}
	rec.Start = start
	rec.End = end
	// score
	scoreStr := string(fields[5])
	if scoreStr == "." {
		rec.ScoreKnown = false
	} else {
		s, err := strconv.ParseFloat(scoreStr, 64)
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "score",
				Msg: "non-numeric: " + scoreStr}
		}
		rec.Score = s
		rec.ScoreKnown = true
	}
	// strand
	switch s := string(fields[6]); s {
	case "+", "-", ".", "?":
		rec.Strand = s[0]
	default:
		return nil, &ParseError{Line: lineNo, Field: "strand",
			Msg: "must be one of +, -, ., ?; got " + s}
	}
	// phase
	phaseStr := string(fields[7])
	if phaseStr == "." {
		rec.Phase = PhaseUnknown
	} else {
		p, err := strconv.Atoi(phaseStr)
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "phase",
				Msg: "non-integer: " + phaseStr}
		}
		if p < 0 || p > 2 {
			return nil, &ParseError{Line: lineNo, Field: "phase",
				Msg: "must be 0, 1, 2, or .; got " + phaseStr}
		}
		rec.Phase = p
	}
	// attributes
	rec.Attributes = parseAttributes(string(fields[8]))
	return rec, nil
}