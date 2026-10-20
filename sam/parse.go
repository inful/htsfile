package sam

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Decoder reads a SAM document incrementally.
//
// Typical usage:
//
//	dec := sam.NewDecoder(r)
//	hdr, _ := dec.Header()    // optional; nil if absent
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    ...
//	}
//
// A Decoder is not safe for concurrent use; create one per
// goroutine.
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

	// Header state.
	hdrRead bool
	hdr     *Header
	hdrErr  error

	// If we read past the header and parked a data line before the
	// caller asked for Record, we return that line first.
	pending     []byte
	pendingLine int
	pendingSet  bool
}

// DefaultMaxLineLength is the per-line size cap applied by
// [NewDecoder]. Lines longer than this cause a parse error rather
// than an unbounded allocation; the value is large enough for any
// real-world SAM line (SEQ is the largest, bounded by sequencing
// technology) and small enough to prevent OOM on hostile input.
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

// Header returns the file's header, or (nil, nil) if the file has
// no "@" lines. The first call reads through blank lines and
// non-header lines until it finds the first "@" line (or the first
// data line, signalling "no header"). The result is cached;
// subsequent calls return the same value.
//
// If Record is called before Header, Header is invoked implicitly
// so the header / first-data-line boundary is correctly handled.
func (d *Decoder) Header() (*Header, error) {
	if d.hdrRead {
		return d.hdr, d.hdrErr
	}
	d.hdrRead = true
	for {
		raw, lineNo, err := d.readLine()
		if err != nil {
			if errors.Is(err, io.EOF) {
				return nil, nil
			}
			return nil, err
		}
		if isBlank(raw) {
			continue
		}
		if raw[0] != '@' {
			// Data line; park it for Record and report no header.
			d.pending = raw
			d.pendingLine = lineNo
			d.pendingSet = true
			return nil, nil
		}
		err = d.parseHeader(raw, lineNo)
		if err != nil {
			d.hdrErr = err
			// Return d.hdr (the partial header) alongside the error
			// so callers can inspect what parsed before the failure.
			return d.hdr, err
		}
		// Continue reading until a non-header line.
		for {
			raw, lineNo, err := d.readLine()
			if err != nil {
				if errors.Is(err, io.EOF) {
					return d.hdr, nil
				}
				return nil, err
			}
			if isBlank(raw) {
				continue
			}
			if raw[0] != '@' {
				// Data line; park it for Record.
				d.pending = raw
				d.pendingLine = lineNo
				d.pendingSet = true
				return d.hdr, nil
			}
			if err := d.parseHeader(raw, lineNo); err != nil {
				d.hdrErr = err
				return d.hdr, err
			}
		}
	}
}

// Record returns the next alignment record. When no more records
// are available it returns (nil, [io.EOF]). Calling Record after
// EOF continues to return (nil, [io.EOF]).
//
// If Header has not yet been called, Record calls it implicitly.
func (d *Decoder) Record() (*Record, error) {
	if !d.hdrRead {
		if _, err := d.Header(); err != nil {
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
			if isBlank(r) {
				continue
			}
			if r[0] == '@' {
				return nil, &ParseError{Line: l, Msg: "unexpected header line in data section: " + string(r)}
			}
			raw = r
			lineNo = l
		}
		return parseRecord(raw, lineNo)
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

// parseHeader parses one "@" header line and appends it to the
// appropriate section of d.hdr. The caller passes raw with the
// leading "@" and WITHOUT the trailing newline.
func (d *Decoder) parseHeader(raw []byte, lineNo int) error {
	// raw starts with '@' (e.g. "@HD\tVN:1.6\tSO:coordinate").
	// Reject lines that are too short to contain the two-character
	// type tag, and lines where raw[3] is non-tab/non-EOF (which
	// means the type tag is more than two characters).
	if len(raw) < 3 || raw[0] != '@' {
		return &ParseError{Line: lineNo, Field: "header",
			Msg: "header line must start with '@' followed by a two-character type tag"}
	}
	if len(raw) > 3 && raw[3] != '\t' {
		return &ParseError{Line: lineNo, Field: "header",
			Msg: "header type tag must be exactly two characters"}
	}
	// The header type is the two characters after '@'.
	tag := string(raw[1:3])
	rest := strings.TrimPrefix(string(raw[3:]), "\t")
	pairs := parseTagPairs(rest)
	hdr := d.hdr
	switch tag {
	case "HD":
		if hdr == nil {
			hdr = &Header{}
		}
		if hdr.HD != nil {
			return &ParseError{Line: lineNo, Field: "header",
				Msg: "duplicate @HD line (only one @HD is allowed per file)"}
		}
		hdr.HD = &HD{
			VN:    pairs["VN"],
			SO:    pairs["SO"],
			GO:    pairs["GO"],
			SS:    pairs["SS"],
			Other: pairs.otherExcept("VN", "SO", "GO", "SS"),
			Raw:   string(raw),
		}
		d.hdr = hdr
		return nil
	case "SQ":
		if hdr == nil {
			hdr = &Header{}
		}
		ln, err := parseInt64Field(pairs["LN"], lineNo, "SQ/LN")
		if err != nil {
			return err
		}
		hdr.SQs = append(hdr.SQs, &SQ{
			SN:    pairs["SN"],
			LN:    ln,
			AH:    pairs["AH"],
			AN:    pairs["AN"],
			AS:    pairs["AS"],
			DS:    pairs["DS"],
			M5:    pairs["M5"],
			TP:    pairs["TP"],
			UR:    pairs["UR"],
			SP:    pairs["SP"],
			Other: pairs.otherExcept("SN", "LN", "AH", "AN", "AS", "DS", "M5", "TP", "UR", "SP"),
			Raw:   string(raw),
		})
		d.hdr = hdr
		return nil
	case "RG":
		if hdr == nil {
			hdr = &Header{}
		}
		hdr.RGs = append(hdr.RGs, &RG{
			ID:    pairs["ID"],
			SM:    pairs["SM"],
			LB:    pairs["LB"],
			DS:    pairs["DS"],
			PU:    pairs["PU"],
			PI:    pairs["PI"],
			CN:    pairs["CN"],
			DT:    pairs["DT"],
			PL:    pairs["PL"],
			PG:    pairs["PG"],
			Other: pairs.otherExcept("ID", "SM", "LB", "DS", "PU", "PI", "CN", "DT", "PL", "PG"),
			Raw:   string(raw),
		})
		d.hdr = hdr
		return nil
	case "PG":
		if hdr == nil {
			hdr = &Header{}
		}
		hdr.PGs = append(hdr.PGs, &PG{
			ID:    pairs["ID"],
			PN:    pairs["PN"],
			CL:    pairs["CL"],
			PP:    pairs["PP"],
			VN:    pairs["VN"],
			Other: pairs.otherExcept("ID", "PN", "CL", "PP", "VN"),
			Raw:   string(raw),
		})
		d.hdr = hdr
		return nil
	case "CO":
		if hdr == nil {
			hdr = &Header{}
		}
		hdr.COs = append(hdr.COs, string(raw))
		d.hdr = hdr
		return nil
	default:
		// Unknown @XX line — preserve verbatim in COs so the
		// encoder can round-trip it.
		if hdr == nil {
			hdr = &Header{}
		}
		hdr.COs = append(hdr.COs, string(raw))
		d.hdr = hdr
		return nil
	}
}

// parseRecord parses one tab-delimited alignment line.
func parseRecord(raw []byte, lineNo int) (*Record, error) {
	fields := bytes.Split(raw, []byte{'\t'})
	if len(fields) < 11 {
		return nil, &ParseError{Line: lineNo, Msg: "data line has " +
			strconv.Itoa(len(fields)) + " fields; need at least 11"}
	}
	rec := &Record{Line: lineNo}
	rec.QName = string(fields[0])
	flag, err := strconv.Atoi(string(fields[1]))
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "FLAG", Msg: "non-numeric: " + string(fields[1])}
	}
	rec.Flag = flag
	rec.RName = string(fields[2])
	pos, err := strconv.ParseInt(string(fields[3]), 10, 64)
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "POS", Msg: "non-numeric: " + string(fields[3])}
	}
	rec.Pos = pos
	mapq, err := strconv.Atoi(string(fields[4]))
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "MAPQ", Msg: "non-numeric: " + string(fields[4])}
	}
	rec.MapQ = mapq
	rec.Cigar = string(fields[5])
	rec.RNext = string(fields[6])
	pnext, err := strconv.ParseInt(string(fields[7]), 10, 64)
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "PNEXT", Msg: "non-numeric: " + string(fields[7])}
	}
	rec.PNext = pnext
	tlen, err := strconv.ParseInt(string(fields[8]), 10, 64)
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "TLEN", Msg: "non-numeric: " + string(fields[8])}
	}
	rec.TLen = tlen
	rec.Seq = string(fields[9])
	rec.Qual = string(fields[10])

	// Optional tag fields: TAG:TYPE:VALUE
	for i := 11; i < len(fields); i++ {
		tag, err := parseTag(string(fields[i]), lineNo)
		if err != nil {
			return nil, err
		}
		rec.Tags = append(rec.Tags, tag)
	}
	return rec, nil
}

// parseTag parses one "TAG:TYPE:VALUE" field. The value may itself
// contain colons (for Z-type strings) or commas (for B-type
// arrays); we only split on the first two colons.
func parseTag(s string, lineNo int) (Tag, error) {
	// Find the first ':' (key separator) and the second ':' (type
	// separator). Anything after the second ':' is the value.
	first := strings.IndexByte(s, ':')
	if first < 0 {
		return Tag{}, &ParseError{Line: lineNo, Field: "tag", Msg: "missing type separator: " + s}
	}
	key := s[:first]
	if len(key) != 2 {
		return Tag{}, &ParseError{Line: lineNo, Field: "tag", Msg: "tag key must be 2 characters: " + key}
	}
	rest := s[first+1:]
	second := strings.IndexByte(rest, ':')
	if second < 0 {
		return Tag{}, &ParseError{Line: lineNo, Field: "tag", Msg: "missing value separator: " + s}
	}
	typeChar := rest[0]
	tt, err := ParseTagType(typeChar)
	if err != nil {
		return Tag{}, &ParseError{Line: lineNo, Field: "tag", Msg: err.Error()}
	}
	return Tag{
		Key:   key,
		Type:  tt,
		Value: rest[second+1:],
	}, nil
}

// tagPairs is a key=value map with an "other" view that excludes
// well-known fields. The SAM header format is a flat tag list; we
// keep all keys in a single map and let callers extract "Other"
// via otherExcept when they want the well-known fields stripped.
type tagPairs map[string]string

func (m tagPairs) otherExcept(keys ...string) map[string]string {
	excluded := make(map[string]struct{}, len(keys))
	for _, k := range keys {
		excluded[k] = struct{}{}
	}
	out := make(map[string]string)
	for k, v := range m {
		if _, ok := excluded[k]; ok {
			continue
		}
		out[k] = v
	}
	return out
}

// parseTagPairs parses a tab-separated (or space-separated) list
// of key:value pairs. The value may contain colons (we split on
// the first colon after the key).
func parseTagPairs(s string) tagPairs {
	out := make(tagPairs)
	i := 0
	for i < len(s) {
		// skip whitespace
		for i < len(s) && (s[i] == ' ' || s[i] == '\t') {
			i++
		}
		if i >= len(s) {
			break
		}
		// read key
		kstart := i
		for i < len(s) && s[i] != ':' && s[i] != ' ' && s[i] != '\t' {
			i++
		}
		key := s[kstart:i]
		if i >= len(s) || s[i] != ':' {
			// malformed: no ':' for this key; skip to next ws
			for i < len(s) && s[i] != ' ' && s[i] != '\t' {
				i++
			}
			continue
		}
		i++ // consume ':'
		vstart := i
		// value: until next whitespace, OR for quoted-like values
		// (Z type may contain colons) — but the SAM header format
		// doesn't use quoting. So read until next whitespace.
		for i < len(s) && s[i] != ' ' && s[i] != '\t' {
			i++
		}
		out[key] = s[vstart:i]
	}
	return out
}

func parseInt64Field(s string, lineNo int, field string) (int64, error) {
	if s == "" {
		return 0, &ParseError{Line: lineNo, Field: field, Msg: "missing"}
	}
	n, err := strconv.ParseInt(s, 10, 64)
	if err != nil {
		return 0, &ParseError{Line: lineNo, Field: field, Msg: "non-numeric: " + s}
	}
	return n, nil
}
