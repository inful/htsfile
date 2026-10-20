package bed

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Decoder reads a BED document incrementally.
//
// Typical usage:
//
//	dec := bed.NewDecoder(r)
//	track, _ := dec.Track()    // optional, may be nil
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    ...
//	}
//
// A Decoder is not safe for concurrent use; create one per goroutine.
type Decoder struct {
	// scanner reads lines from the input with a bounded buffer
	// (DefaultMaxLineLength or whatever was set via
	// NewDecoderWithMaxLine).
	scanner *bufio.Scanner

	// maxLineLength is the per-line byte cap used to configure the
	// scanner. The scanner returns bufio.ErrTooLong when a line
	// exceeds this; readLine translates that into a ParseError.
	maxLineLength int

	line int
	eof  bool

	// Track state.
	trackRead bool
	track     *Track
	trackErr  error

	// If we read past the headers and parked a data line before the
	// caller asked for a Record, we return that line first.
	pending     []byte
	pendingLine int
	pendingSet  bool
}

// DefaultMaxLineLength is the per-line size cap applied by
// [NewDecoder]. Lines longer than this cause a parse error rather
// than an unbounded allocation; the value is large enough for any
// real-world BED line and small enough to prevent OOM on hostile
// input.
const DefaultMaxLineLength = 10 << 20 // 10 MB

// NewDecoder wraps r in a [Decoder] with the default
// [DefaultMaxLineLength] per-line size cap. Use
// [NewDecoderWithMaxLine] to override the cap.
func NewDecoder(r io.Reader) *Decoder {
	return NewDecoderWithMaxLine(r, DefaultMaxLineLength)
}

// NewDecoderWithMaxLine is like [NewDecoder] but lets the caller
// set a custom per-line size cap. Pass maxLineLength <= 0 to keep
// the default cap.
func NewDecoderWithMaxLine(r io.Reader, maxLineLength int) *Decoder {
	if maxLineLength <= 0 {
		maxLineLength = DefaultMaxLineLength
	}
	sc := bufio.NewScanner(r)
	// Both the initial buffer cap and the explicit max must equal
	// the limit; otherwise the Scanner silently grows the buffer
	// up to cap(buf) regardless of max.
	sc.Buffer(make([]byte, 0, maxLineLength), maxLineLength)
	return &Decoder{scanner: sc, maxLineLength: maxLineLength}
}

// Track returns the file's "track" header line, or (nil, nil) if
// there is no track line. The first call reads through browser /
// comment / blank lines until it finds the track line or the first
// data line. The result is cached; subsequent calls return the same
// value.
//
// Browser lines (e.g. "browser position chr1:1-1000") are display
// hints and are skipped here. They are not preserved for round-trip
// in this version.
func (d *Decoder) Track() (*Track, error) {
	if d.trackRead {
		return d.track, d.trackErr
	}
	d.trackRead = true
	for {
		raw, lineNo, err := d.readLine()
		if err != nil {
			if errors.Is(err, io.EOF) {
				return nil, nil
			}
			return nil, err
		}
		if isBlank(raw) || isComment(raw) {
			continue
		}
		if bytes.HasPrefix(raw, []byte("browser")) {
			continue
		}
		if bytes.HasPrefix(raw, []byte("track")) {
			t, err := parseTrack(raw)
			d.track = t
			d.trackErr = err
			return t, err
		}
		// Data line: park it for Record() and report no track.
		d.pending = raw
		d.pendingLine = lineNo
		d.pendingSet = true
		return nil, nil
	}
}

// Record returns the next data record. When no more records are
// available it returns (nil, [io.EOF]). Calling Record after EOF
// continues to return (nil, [io.EOF]).
//
// If Track has not yet been called, Record calls it implicitly so the
// header / first-data-line boundary is correctly handled.
func (d *Decoder) Record() (*Record, error) {
	if !d.trackRead {
		if _, err := d.Track(); err != nil {
			return nil, err
		}
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
			if isBlank(r) || isComment(r) {
				continue
			}
			if bytes.HasPrefix(r, []byte("track")) || bytes.HasPrefix(r, []byte("browser")) {
				return nil, &ParseError{Line: l, Msg: "unexpected header line in data section: " + string(r)}
			}
			raw = r
			lineNo = l
		}
		return parseRecord(raw, lineNo)
	}
}

// readLine reads one line. Returns (raw, lineNo, nil) on success,
// (nil, lineNo, io.EOF) at end of file. Trailing newline / CR is
// stripped. The Scanner internally enforces [DefaultMaxLineLength]
// (or whatever the caller set via [NewDecoderWithMaxLine]); lines
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
	// Strip a trailing CR (Scanner doesn't include it; be safe
	// anyway in case of mixed line endings).
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

func isComment(b []byte) bool {
	if len(b) == 0 || b[0] != '#' {
		return false
	}
	// The "#" must not be followed by an alpha/digit, to avoid
	// treating an actual data field as a comment. BED has no formal
	// comment syntax, but the de-facto convention is a leading "#".
	// We treat any line starting with "#" as a comment.
	return true
}

// parseTrack parses a "track name=... description=..." header line.
// The Raw field is always set to the original line for round-trip
// fidelity; known fields are extracted into the struct.
func parseTrack(raw []byte) (*Track, error) {
	t := &Track{Raw: string(raw)}
	// Strip the leading "track" keyword and any whitespace after it.
	rest := bytes.TrimSpace(raw[5:])
	// Split into name=value pairs. Use a simple state machine that
	// handles double-quoted values.
	type kv struct{ k, v string }
	var pairs []kv
	i := 0
	for i < len(rest) {
		// skip whitespace
		for i < len(rest) && (rest[i] == ' ' || rest[i] == '\t') {
			i++
		}
		if i >= len(rest) {
			break
		}
		// read key
		kstart := i
		for i < len(rest) && rest[i] != '=' && rest[i] != ' ' && rest[i] != '\t' {
			i++
		}
		key := string(rest[kstart:i])
		if i >= len(rest) || rest[i] != '=' {
			// Malformed: no '=' for this key. Skip to next whitespace.
			for i < len(rest) && rest[i] != ' ' && rest[i] != '\t' {
				i++
			}
			continue
		}
		i++ // consume '='
		// read value (quoted or unquoted)
		if i < len(rest) && rest[i] == '"' {
			i++ // consume opening quote
			vstart := i
			terminated := false
			for i < len(rest) {
				if rest[i] == '\\' && i+1 < len(rest) {
					i += 2
					continue
				}
				if rest[i] == '"' {
					terminated = true
					break
				}
				i++
			}
			if !terminated {
				return nil, fmt.Errorf("unterminated quoted value in attribute %q", key)
			}
			val := string(rest[vstart:i])
			i++ // consume closing quote
			pairs = append(pairs, kv{key, val})
		} else {
			vstart := i
			for i < len(rest) && rest[i] != ' ' && rest[i] != '\t' {
				i++
			}
			pairs = append(pairs, kv{key, string(rest[vstart:i])})
		}
	}
	for _, p := range pairs {
		switch p.k {
		case "name":
			t.Name = p.v
		case "description":
			t.Description = p.v
		case "useScore":
			t.UseScore = p.v == "1"
		case "itemRgb":
			t.ItemRgbOn = strings.EqualFold(p.v, "On")
		case "type":
			t.Type = p.v
		case "visibility":
			t.Visibility = p.v
		}
	}
	return t, nil
}

// parseRecord parses one tab-delimited data line into a Record. The
// column count is auto-detected from the number of fields (3..=12).
// Extra trailing fields (>12) are ignored with no error.
func parseRecord(raw []byte, lineNo int) (*Record, error) {
	fields := strings.Split(string(raw), "\t")
	if len(fields) < 3 {
		return nil, &ParseError{Line: lineNo, Msg: fmt.Sprintf(
			"data line has %d fields; need at least 3", len(fields))}
	}
	rec := &Record{Line: lineNo}
	if fields[0] == "" {
		return nil, &ParseError{Line: lineNo, Field: "chrom", Msg: "empty chromosome"}
	}
	rec.Chrom = fields[0]

	start, err := strconv.ParseInt(fields[1], 10, 64)
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "chromStart",
			Msg: "non-numeric chromStart: " + fields[1]}
	}
	rec.ChromStart = start

	end, err := strconv.ParseInt(fields[2], 10, 64)
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "chromEnd",
			Msg: "non-numeric chromEnd: " + fields[2]}
	}
	rec.ChromEnd = end

	// Optional fields. The spec requires a contiguous prefix; we
	// detect and surface that via ColumnCount. Here we populate
	// whatever is present; the validator (later) reports non-contiguous
	// prefixes.
	if len(fields) > 3 {
		s := fields[3]
		rec.Name = &s
	}
	if len(fields) > 4 {
		n, err := strconv.Atoi(fields[4])
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "score",
				Msg: "non-numeric score: " + fields[4]}
		}
		rec.Score = &n
	}
	if len(fields) > 5 {
		if fields[5] != "" {
			st, err := ParseStrand(fields[5][0])
			if err != nil {
				return nil, &ParseError{Line: lineNo, Field: "strand", Msg: err.Error()}
			}
			rec.Strand = &st
		}
	}
	if len(fields) > 6 {
		n, err := strconv.ParseInt(fields[6], 10, 64)
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "thickStart",
				Msg: "non-numeric thickStart: " + fields[6]}
		}
		rec.ThickStart = &n
	}
	if len(fields) > 7 {
		n, err := strconv.ParseInt(fields[7], 10, 64)
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "thickEnd",
				Msg: "non-numeric thickEnd: " + fields[7]}
		}
		rec.ThickEnd = &n
	}
	if len(fields) > 8 {
		rgb, err := parseRGB(fields[8])
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "itemRgb", Msg: err.Error()}
		}
		rec.ItemRgb = &rgb
	}
	if len(fields) > 9 {
		n, err := strconv.Atoi(fields[9])
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "blockCount",
				Msg: "non-numeric blockCount: " + fields[9]}
		}
		rec.BlockCount = &n
	}
	if len(fields) > 10 {
		sizes, err := parseIntList(fields[10])
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "blockSizes", Msg: err.Error()}
		}
		rec.BlockSizes = sizes
	}
	if len(fields) > 11 {
		starts, err := parseIntList(fields[11])
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "blockStarts", Msg: err.Error()}
		}
		rec.BlockStarts = starts
	}
	return rec, nil
}

// parseRGB parses "R,G,B" into an RGB. Each component must be 0..=255.
func parseRGB(s string) (RGB, error) {
	parts := strings.Split(s, ",")
	if len(parts) != 3 {
		return RGB{}, fmt.Errorf("itemRgb must be R,G,B; got %q", s)
	}
	var rgb RGB
	for i, p := range parts {
		n, err := strconv.Atoi(strings.TrimSpace(p))
		if err != nil {
			return RGB{}, fmt.Errorf("non-numeric RGB component %q: %w", p, err)
		}
		if n < 0 || n > 255 {
			return RGB{}, fmt.Errorf("RGB component out of range [0,255]: %d", n)
		}
		switch i {
		case 0:
			rgb.R = uint8(n)
		case 1:
			rgb.G = uint8(n)
		case 2:
			rgb.B = uint8(n)
		}
	}
	return rgb, nil
}

// parseIntList parses a comma-separated list of integers.
func parseIntList(s string) ([]int64, error) {
	if s == "" {
		return nil, fmt.Errorf("empty comma-separated list")
	}
	parts := strings.Split(s, ",")
	out := make([]int64, 0, len(parts))
	for _, p := range parts {
		n, err := strconv.ParseInt(p, 10, 64)
		if err != nil {
			return nil, fmt.Errorf("non-numeric element %q: %w", p, err)
		}
		out = append(out, n)
	}
	return out, nil
}
