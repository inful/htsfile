package fastq

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
)

// DefaultPhredOffset is the offset used by [NewDecoder] for
// Phred quality decoding. 33 corresponds to Sanger / Illumina
// 1.8+ (the dominant modern convention); 64 corresponds to
// older Illumina 1.3–1.7.
const DefaultPhredOffset uint8 = 33

// Illumina13PhredOffset is the offset used by older Illumina
// 1.3–1.7 encoders. Pass this to [Decoder.SetPhredOffset] when
// parsing data whose quality bytes are in the [64, 104] range.
const Illumina13PhredOffset uint8 = 64

// DefaultMaxLineLength is the per-line size cap applied by
// [NewDecoder]. Lines longer than this cause a parse error
// rather than an unbounded allocation; the value is large
// enough for any real-world FASTQ record (typical sequence
// lines are <200 chars and quality lines match them) and
// small enough to prevent OOM on hostile input.
const DefaultMaxLineLength = 10 << 20 // 10 MB

// Decoder reads a FASTQ document incrementally from an
// underlying [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("reads.fq")
//	defer f.Close()
//	dec := fastq.NewDecoder(f)
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	}
//
// The first record's quality bytes drive Phred-offset
// auto-detection: if every byte is < 64 the data is treated as
// Phred+33 (the default); if any byte is ≥ 64 the data is
// treated as Phred+64. Subsequent records use the same offset.
// Callers can disable auto-detect and force an offset via
// [Decoder.SetPhredOffset].
//
// Solexa Solexa+64 data uses the same byte range as Phred+64
// but a different formula. Callers with known Solexa+64 data
// should call [Decoder.SetPhredEncoding] (Solexa) before
// reading. The library does not auto-detect Solexa.
//
// A Decoder is not safe for concurrent use; create one per
// goroutine.
type Decoder struct {
	scanner *bufio.Scanner

	maxLineLength int

	// Phred offset config.
	phredOffset uint8
	autoDetect  bool

	// Phred formula config.
	phredEncoding PhredEncoding

	// After the first record, this holds the detected (or
	// configured) offset used for all records.
	offsetDecided bool

	line int
	eof  bool
}

// NewDecoder wraps r in a [Decoder] with the default
// [DefaultMaxLineLength] per-line size cap and Phred+33
// default. Auto-detect is enabled: the first record's quality
// bytes decide whether the offset is actually 33 or 64.
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
	return &Decoder{
		scanner:       sc,
		maxLineLength: maxLineLength,
		phredOffset:   DefaultPhredOffset,
		phredEncoding: PhredEncodingPhred,
		autoDetect:    true,
	}
}

// SetPhredOffset overrides the Phred offset (and disables
// auto-detect for the lifetime of this Decoder). The standard
// offsets are [DefaultPhredOffset] (33, Sanger / Illumina
// 1.8+) and [Illumina13PhredOffset] (64, older Illumina).
func (d *Decoder) SetPhredOffset(offset uint8) {
	d.phredOffset = offset
	d.autoDetect = false
}

// SetPhredEncoding sets the formula the Qual bytes are
// intended to be decoded with. The default is
// PhredEncodingPhred; Solexa data should set
// PhredEncodingSolexa. The byte range of the Qual line is
// unaffected by this setting — see [PhredEncoding] for the
// formula difference.
func (d *Decoder) SetPhredEncoding(e PhredEncoding) {
	d.phredEncoding = e
}

// PhredOffset returns the offset in effect for the most
// recently decoded record (or the configured default if no
// record has been decoded yet). After the first record, this
// is stable for the lifetime of the Decoder.
func (d *Decoder) PhredOffset() uint8 {
	return d.phredOffset
}

// PhredEncoding returns the encoding in effect for the
// Decoder (default [PhredEncodingPhred]).
func (d *Decoder) PhredEncoding() PhredEncoding {
	return d.phredEncoding
}

// Record returns the next FASTQ record, or (nil, [io.EOF])
// when no more records are available. Calling Record after
// EOF continues to return (nil, [io.EOF]).
//
// FASTQ has no header section that needs to be read
// separately — records are self-delimiting. The streaming
// pattern matches the other packages' Decoder shapes; there
// is no [Decoder.Header] method.
func (d *Decoder) Record() (*Record, error) {
	hdr, hdrLine, err := d.readNonBlankLine()
	if err != nil {
		return nil, err
	}
	if len(hdr) == 0 || hdr[0] != '@' {
		return nil, &ParseError{Line: hdrLine, Field: "header",
			Msg: "FASTQ record must start with '@'"}
	}
	rec := &Record{
		Line:   hdrLine,
		Header: string(hdr[1:]),
	}
	// Sequence line.
	seqLine, _, err := d.readLine()
	if err != nil {
		if errors.Is(err, io.EOF) {
			return nil, &ParseError{Line: hdrLine, Field: "sequence",
				Msg: "unexpected EOF after header"}
		}
		return nil, &ParseError{Line: hdrLine, Field: "sequence",
			Msg: err.Error()}
	}
	rec.Seq = append(rec.Seq, seqLine...)
	// Plus line.
	plusLine, plusLineNo, err := d.readLine()
	if err != nil {
		if errors.Is(err, io.EOF) {
			return nil, &ParseError{Line: hdrLine, Field: "plus",
				Msg: "unexpected EOF after sequence"}
		}
		return nil, &ParseError{Line: plusLineNo, Field: "plus",
			Msg: err.Error()}
	}
	if len(plusLine) == 0 || plusLine[0] != '+' {
		return nil, &ParseError{Line: plusLineNo, Field: "plus",
			Msg: "FASTQ record separator must start with '+'"}
	}
	rec.PlusLine = string(bytes.TrimLeft(plusLine[1:], " \t"))
	// Qualities line.
	qual, qualLineNo, err := d.readLine()
	if err != nil {
		if errors.Is(err, io.EOF) {
			return nil, &ParseError{Line: plusLineNo, Field: "qualities",
				Msg: "unexpected EOF after '+' line"}
		}
		return nil, &ParseError{Line: qualLineNo, Field: "qualities",
			Msg: err.Error()}
	}
	// Phred offset detection (first record only).
	if d.autoDetect && !d.offsetDecided {
		// Heuristic: Phred+33 quality bytes for Q0-Q40 are in
		// the [33, 73] ASCII range; Phred+64 spans [64, 104].
		// A max byte > 73 means the data cannot be Phred+33
		// (Q41+ is unrepresented) and must be Phred+64. A
		// max byte <= 73 and min byte >= 33 is Phred+33.
		// Min byte < 33 is an error either way (Q values < 0
		// don't exist in either encoding).
		if len(qual) > 0 {
			minQ, maxQ := qualRange(qual)
			// Doc-promised rules: a byte < 33 is invalid for
			// Phred+33 (Q-1 doesn't exist) and a byte < 64
			// in the 64..73 range can still be valid Phred+64
			// with Q0. The lowest valid byte across the two
			// encodings is 33 (Phred+33 Q0).
			if minQ < 33 {
				return nil, &ParseError{Line: qualLineNo, Field: "qualities",
					Msg: fmt.Sprintf("quality byte %d below Phred minimum 33", minQ)}
			}
			if maxQ > 73 {
				d.phredOffset = Illumina13PhredOffset
			} else if minQ >= 64 {
				// Documented rule: min >= 64 unambiguously
				// indicates Phred+64.
				d.phredOffset = Illumina13PhredOffset
			}
		}
		d.offsetDecided = true
	}
	rec.Qual = append(rec.Qual, qual...)
	rec.PhredOffset = d.phredOffset
	rec.PhredEncoding = d.phredEncoding
	return rec, nil
}

// readNonBlankLine reads lines until one is non-empty or EOF.
// Returns the line bytes (without trailing CR/LF), the 1-based
// line number, and any error. Used for the @header line so
// that records separated by blank lines parse cleanly.
//
// mirrors readLine's so callers can pass through either.
//
//nolint:gocritic // unnamed results are intentional; the signature
func (d *Decoder) readNonBlankLine() ([]byte, int, error) {
	for {
		raw, lineNo, err := d.readLine()
		if err != nil {
			return nil, lineNo, err
		}
		if len(raw) > 0 {
			return raw, lineNo, nil
		}
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

// qualRange returns the (min, max) byte values in b. If b is
// empty, both are 0.
//
// single-use helper used inline during auto-detect.
//
//nolint:gocritic // unnamed results are intentional; this is a
func qualRange(b []byte) (byte, byte) {
	if len(b) == 0 {
		return 0, 0
	}
	mn, mx := b[0], b[0]
	for _, c := range b[1:] {
		if c < mn {
			mn = c
		}
		if c > mx {
			mx = c
		}
	}
	return mn, mx
}
