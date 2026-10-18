package vcf

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"
)

// Decoder reads a VCF document incrementally.
//
// Typical usage:
//
//	dec := vcf.NewDecoder(r)
//	hdr, err := dec.Header()
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    // ...
//	}
//
// A [Decoder] is not safe for concurrent use; create one per goroutine.
type Decoder struct {
	// scanner reads lines from the input with a bounded buffer
	// (DefaultMaxLineLength or whatever was set via
	// NewDecoderWithMaxLine).
	scanner *bufio.Scanner

	// maxLineLength is the per-line byte cap used to configure the
	// scanner. The scanner returns bufio.ErrTooLong when a line
	// exceeds this; readLine translates that into a *ParseError.
	maxLineLength int

	// line is the 1-based line number of the most recently read line.
	line int

	// hdr is cached after the first successful [Decoder.Header] call.
	hdr *Header
	// hdrErr, when non-nil, is returned by every subsequent [Decoder.Header]
	// call. We surface header parse failures persistently so callers don't
	// re-read malformed input and end up with different results.
	hdrErr error
	// hdrParsed is true once Header() has run to completion (with or without
	// success). It prevents re-parsing.
	hdrParsed bool

	// eof is true once the underlying reader has been exhausted.
	eof bool
}

// DefaultMaxLineLength is the per-line size cap applied by
// [NewDecoder]. Lines longer than this cause a parse error rather
// than an unbounded allocation; the value is large enough for any
// real-world VCF line (even a long INFO string rarely exceeds a few
// KB) and small enough to prevent OOM on hostile input.
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

// Header parses the VCF header. The first call reads ##meta lines until the
// #CHROM header line is seen; subsequent calls return the same [*Header]
// without re-reading the stream.
func (d *Decoder) Header() (*Header, error) {
	if d.hdrParsed {
		return d.hdr, d.hdrErr
	}
	d.hdrParsed = true
	h, err := d.parseHeader()
	d.hdr = h
	d.hdrErr = err
	return h, err
}

// Record returns the next data record. When no more records are available,
// it returns (nil, [io.EOF]). Calling Record after EOF is allowed and
// continues to return (nil, [io.EOF]).
func (d *Decoder) Record() (*Record, error) {
	if !d.hdrParsed {
		if _, err := d.Header(); err != nil {
			return nil, err
		}
	}
	for {
		raw, lineNo, err := d.readLine()
		if errors.Is(err, io.EOF) {
			return nil, io.EOF
		}
		if err != nil {
			return nil, err
		}
		// Allow empty lines in the data section.
		if len(raw) == 0 {
			continue
		}
		rec, err := parseRecord(raw, lineNo, d.hdr)
		if err != nil {
			return nil, err
		}
		return rec, nil
	}
}

// readLine reads one line from the underlying reader. The returned
// slice is the line without its trailing newline. The Scanner
// internally enforces [DefaultMaxLineLength] (or whatever the
// caller set via [NewDecoderWithMaxLine]); lines longer than that
// return a [ParseError] without unbounded allocation.
//
// On the first call, a leading UTF-8 BOM (EF BB BF) is detected
// and rejected with a [ParseError] — VCF files must not contain
// a BOM per the spec.
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
	// Copy the scanner's bytes because the slice is reused on the
	// next Scan call.
	src := d.scanner.Bytes()
	raw = make([]byte, len(src))
	copy(raw, src)
	// Strip a trailing CR from the copy (Scanner preserves the
	// newline but not a leading CR; be safe anyway).
	if n := len(raw); n > 0 && raw[n-1] == '\r' {
		raw = raw[:n-1]
	}
	// BOM detection on the first line.
	if d.line == 0 && len(raw) >= 3 && bytes.Equal(raw[0:3], []byte{0xEF, 0xBB, 0xBF}) {
		return nil, lineNo, &ParseError{Line: lineNo, Field: "<BOM>",
			Msg: "VCF files must not contain a UTF-8 byte order mark"}
	}
	// Non-printable character rejection. The spec forbids
	// the ASCII control ranges U+0000..U+0008, U+000B..U+000C,
	// and U+000E..U+001F. For UTF-8 input the single-byte
	// (ASCII) form is the only relevant one; multi-byte
	// sequences start at 0x80+ and are out of these ranges.
	for i, b := range raw {
		if isDisallowedControl(b) {
			return nil, lineNo, &ParseError{Line: lineNo,
				Field: "<non-printable>",
				Msg: fmt.Sprintf("non-printable byte 0x%02X at column %d is not allowed in VCF", b, i+1)}
		}
	}
	d.line = lineNo
	return raw, lineNo, nil
}

// isDisallowedControl reports whether b is in one of the
// non-printable byte ranges the VCF v4.5 spec forbids:
// U+0000..U+0008, U+000B..U+000C, U+000E..U+001F.
func isDisallowedControl(b byte) bool {
	switch {
	case b <= 0x08:
		return true
	case b == 0x0B, b == 0x0C:
		return true
	case b >= 0x0E && b <= 0x1F:
		return true
	}
	return false
}

// parseHeader reads ## meta lines and the #CHROM header line.
func (d *Decoder) parseHeader() (*Header, error) {
	h := &Header{
		Infos:   map[string]*InfoDef{},
		Filters: map[string]*FilterDef{},
		Formats: map[string]*FormatDef{},
		Alts:    map[string]*AltDef{},
	}

	// First line must be ##fileformat=...
	firstLine, lineNo, err := d.readLine()
	if err != nil {
		if errors.Is(err, io.EOF) {
			return nil, &ParseError{Line: lineNo, Msg: "empty document"}
		}
		return nil, err
	}
	if err := h.parseFirstLine(firstLine, lineNo); err != nil {
		return nil, err
	}

	// Subsequent ## / # lines.
	for {
		raw, ln, err := d.readLine()
		if errors.Is(err, io.EOF) {
			return nil, &ParseError{Line: ln, Field: "header",
				Msg: "expected #CHROM header line before end of file"}
		}
		if err != nil {
			return nil, err
		}
		if len(raw) == 0 {
			// Allow blank lines in the meta section.
			continue
		}
		if bytes.HasPrefix(raw, []byte("##")) {
			if err := h.parseMetaLine(raw, ln); err != nil {
				return nil, err
			}
			continue
		}
		if bytes.HasPrefix(raw, []byte("#")) && !bytes.HasPrefix(raw, []byte("##")) {
			// Header line.
			if err := h.parseHeaderLine(raw, ln); err != nil {
				return nil, err
			}
			break
		}
		return nil, &ParseError{Line: ln,
			Msg: "expected '##' or '#CHROM' line, got: " + string(raw)}
	}

	return h, nil
}

// parseFirstLine parses the mandatory first ##fileformat line.
func (h *Header) parseFirstLine(raw []byte, lineNo int) error {
	const prefix = "##fileformat="
	if !bytes.HasPrefix(raw, []byte(prefix)) {
		return &ParseError{Line: lineNo, Field: "fileformat",
			Msg: "first line must be '##fileformat=...', got: " + string(raw)}
	}
	val := string(raw[len(prefix):])
	ff, err := ParseFileFormat(val)
	if err != nil {
		var pe *ParseError
		if errors.As(err, &pe) {
			pe.Line = lineNo
			pe.Field = "fileformat"
			return pe
		}
		return &ParseError{Line: lineNo, Field: "fileformat", Msg: err.Error()}
	}
	h.FileFormat = ff
	h.MetaLines = append(h.MetaLines, MetaLine{Key: "fileformat", Value: val, Raw: string(raw)})
	return nil
}

// parseMetaLine dispatches a single ##key=value or ##key=<...> line.
func (h *Header) parseMetaLine(raw []byte, lineNo int) error {
	rest := string(raw[2:]) // strip "##"

	// Split on first '='.
	eq := strings.IndexByte(rest, '=')
	if eq < 0 {
		return &ParseError{Line: lineNo, Msg: "##meta line missing '=': " + rest}
	}
	key := rest[:eq]
	value := rest[eq+1:]

	if key == "" {
		return &ParseError{Line: lineNo, Msg: "##meta line has empty key"}
	}

	// Common case: "fileformat" is only legal as the first line; reject
	// duplicates here.
	if key == "fileformat" {
		return &ParseError{Line: lineNo, Field: "fileformat",
			Msg: "duplicate ##fileformat line"}
	}

	if !strings.HasPrefix(value, "<") {
		// Unstructured line. Value must be non-empty and not start with '<'.
		if value == "" {
			return &ParseError{Line: lineNo, Field: key, Msg: "empty value"}
		}
		if strings.HasPrefix(value, "<") {
			return &ParseError{Line: lineNo, Field: key,
				Msg: "unstructured ##meta value may not start with '<'"}
		}
		switch key {
		case "fileDate":
			h.FileDate = value
		case "source":
			h.Source = value
		case "reference":
			h.Reference = value
		case "assembly":
			h.Assembly = value
		case "pedigreeDB":
			h.PedigreeDB = value
		case "PEDIGREE":
			// Disambiguator: PEDIGREE is usually structured but is
			// reserved as a key in both forms by the spec. We accept only
			// the structured form below; reject this case.
			return &ParseError{Line: lineNo, Field: "PEDIGREE",
				Msg: "##PEDIGREE must be structured (<...>)"}
		default:
			// Generic meta line: just record it.
		}
		h.MetaLines = append(h.MetaLines, MetaLine{Key: key, Value: value, Raw: string(raw)})
		return nil
	}

	// Structured line: <...>
	if !strings.HasSuffix(value, ">") {
		return &ParseError{Line: lineNo, Field: key,
			Msg: "structured ##meta line missing closing '>'"}
	}
	inner := value[1 : len(value)-1]

	pairs, err := parseStructuredBody(inner)
	if err != nil {
		return &ParseError{Line: lineNo, Field: key, Msg: err.Error()}
	}

	h.MetaLines = append(h.MetaLines, MetaLine{Key: key, Value: value, Raw: string(raw)})
	return h.applyStructured(key, pairs, lineNo)
}

// parseHeaderLine parses "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO[\tFORMAT[\t...samples]]".
func (h *Header) parseHeaderLine(raw []byte, lineNo int) error {
	// Trailing-tab rejection: the spec says "there must be
	// no tab characters at the end of the line". The
	// readLine helper strips a trailing CR but not a
	// trailing tab, so we check here. The split below
	// would silently drop an empty trailing field.
	if n := len(raw); n > 0 && raw[n-1] == '\t' {
		return &ParseError{Line: lineNo, Field: "header",
			Msg: "header line ends with a tab character"}
	}
	fields := bytes.Split(raw, []byte{'\t'})
	if len(fields) < 8 {
		return &ParseError{Line: lineNo, Field: "header",
			Msg: fmt.Sprintf("header line has %d fields; need at least 8", len(fields))}
	}
	expected := []string{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}
	for i, want := range expected {
		if string(fields[i]) != want {
			return &ParseError{Line: lineNo, Field: want,
				Msg: fmt.Sprintf("expected %q at column %d, got %q", want, i+1, fields[i])}
		}
	}
	if len(fields) == 8 {
		return nil
	}
	if string(fields[8]) != "FORMAT" {
		return &ParseError{Line: lineNo, Field: "FORMAT",
			Msg: fmt.Sprintf("expected FORMAT at column 9, got %q", fields[8])}
	}
	samples := make([]string, 0, len(fields)-9)
	seen := map[string]bool{}
	for i := 9; i < len(fields); i++ {
		name := string(fields[i])
		if name == "" {
			return &ParseError{Line: lineNo, Field: "FORMAT",
				Msg: fmt.Sprintf("empty sample name at column %d", i+1)}
		}
		// Sample names must not contain whitespace or tabs;
		// the spec calls this out for several fields and
		// the parser is the natural place to enforce it.
		if strings.ContainsAny(name, " \t") {
			return &ParseError{Line: lineNo, Field: "FORMAT",
				Msg: fmt.Sprintf("sample name %q contains whitespace", name)}
		}
		if seen[name] {
			return &ParseError{Line: lineNo, Field: "FORMAT",
				Msg: fmt.Sprintf("duplicate sample name %q", name)}
		}
		seen[name] = true
		samples = append(samples, name)
	}
	h.SampleNames = samples
	return nil
}

// applyStructured decodes the key=value pairs of a structured meta line into
// the corresponding header field.
func (h *Header) applyStructured(key string, pairs []fieldPair, lineNo int) error {
	switch key {
	case "contig":
		c, err := parseContigAttrs(pairs, lineNo)
		if err != nil {
			return err
		}
		h.Contigs = append(h.Contigs, c)

	case "INFO":
		d := &InfoDef{Other: map[string]string{}}
		for _, p := range pairs {
			switch p.key {
			case "ID":
				d.ID = unquote(p.value)
			case "Number":
				n, err := ParseNumber(p.value)
				if err != nil {
					return &ParseError{Line: lineNo, Field: "INFO/Number", Msg: err.Error()}
				}
				d.Number = n
			case "Type":
				t, err := ParseFieldType(p.value)
				if err != nil {
					return &ParseError{Line: lineNo, Field: "INFO/Type", Msg: err.Error()}
				}
				d.Type = t
			case "Description":
				d.Description = unquote(p.value)
			case "Source":
				d.Source = unquote(p.value)
			case "Version":
				d.Version = unquote(p.value)
			default:
				d.Other[p.key] = unquote(p.value)
			}
		}
		if d.ID == "" {
			return &ParseError{Line: lineNo, Field: "INFO", Msg: "missing required ID"}
		}
		if _, exists := h.Infos[d.ID]; exists {
			return &ParseError{Line: lineNo, Field: "INFO",
				Msg: "duplicate ID " + d.ID}
		}
		if d.Type == TypeFlag && d.Number.Number != 0 {
			return &ParseError{Line: lineNo, Field: "INFO/Number",
				Msg: "Flag fields must have Number=0"}
		}
		h.Infos[d.ID] = d

	case "FILTER":
		d := &FilterDef{Other: map[string]string{}}
		for _, p := range pairs {
			switch p.key {
			case "ID":
				d.ID = unquote(p.value)
			case "Description":
				d.Description = unquote(p.value)
			default:
				d.Other[p.key] = unquote(p.value)
			}
		}
		if d.ID == "" {
			return &ParseError{Line: lineNo, Field: "FILTER",
				Msg: "missing required ID"}
		}
		if _, exists := h.Filters[d.ID]; exists {
			return &ParseError{Line: lineNo, Field: "FILTER",
				Msg: "duplicate ID " + d.ID}
		}
		h.Filters[d.ID] = d

	case "FORMAT":
		d := &FormatDef{Other: map[string]string{}}
		for _, p := range pairs {
			switch p.key {
			case "ID":
				d.ID = unquote(p.value)
			case "Number":
				n, err := ParseNumber(p.value)
				if err != nil {
					return &ParseError{Line: lineNo, Field: "FORMAT/Number", Msg: err.Error()}
				}
				d.Number = n
			case "Type":
				t, err := ParseFieldType(p.value)
				if err != nil {
					return &ParseError{Line: lineNo, Field: "FORMAT/Type", Msg: err.Error()}
				}
				d.Type = t
			case "Description":
				d.Description = unquote(p.value)
			default:
				d.Other[p.key] = unquote(p.value)
			}
		}
		if d.ID == "" {
			return &ParseError{Line: lineNo, Field: "FORMAT",
				Msg: "missing required ID"}
		}
		if _, exists := h.Formats[d.ID]; exists {
			return &ParseError{Line: lineNo, Field: "FORMAT",
				Msg: "duplicate ID " + d.ID}
		}
		h.Formats[d.ID] = d

	case "ALT":
		d := &AltDef{Other: map[string]string{}}
		for _, p := range pairs {
			switch p.key {
			case "ID":
				d.ID = unquote(p.value)
			case "Description":
				d.Description = unquote(p.value)
			default:
				d.Other[p.key] = unquote(p.value)
			}
		}
		if d.ID == "" {
			return &ParseError{Line: lineNo, Field: "ALT",
				Msg: "missing required ID"}
		}
		// VCFv4.5 §1.4.4: the first level type of a
		// symbolic structural variant ALT ID must be
		// DEL, INS, DUP, INV, or CNV. Implementations
		// are free to add subtypes via ":" but the
		// first segment must be one of these five.
		// IUPAC ambiguity codes (R, Y, S, W, K, M,
		// B, D, H, V) are also accepted as bare
		// symbolic-allele IDs.
		first := d.ID
		if i := strings.IndexByte(d.ID, ':'); i >= 0 {
			first = d.ID[:i]
		}
		if iupacIDs[first] {
			// IUPAC code as a bare symbolic-allele
			// ID is accepted (no colon-separated
			// subtypes for IUPAC codes).
		} else {
			switch first {
			case "DEL", "INS", "DUP", "INV", "CNV":
				// ok
			default:
				return &ParseError{Line: lineNo, Field: "ALT",
					Msg: fmt.Sprintf("ALT ID %q has first-level type %q; must be one of DEL, INS, DUP, INV, CNV, or an IUPAC code", d.ID, first)}
			}
		}
		if _, exists := h.Alts[d.ID]; exists {
			return &ParseError{Line: lineNo, Field: "ALT",
				Msg: "duplicate ID " + d.ID}
		}
		h.Alts[d.ID] = d

	case "PEDIGREE":
		d := &Pedigree{Other: map[string]string{}}
		for _, p := range pairs {
			switch p.key {
			case "ID":
				d.ID = unquote(p.value)
			case "Original":
				d.Original = unquote(p.value)
			case "Father":
				d.Father = unquote(p.value)
			case "Mother":
				d.Mother = unquote(p.value)
			default:
				if strings.HasPrefix(p.key, "Name_") {
					d.Names = append(d.Names, unquote(p.value))
					continue
				}
				d.Other[p.key] = unquote(p.value)
			}
		}
		if d.ID == "" {
			return &ParseError{Line: lineNo, Field: "PEDIGREE",
				Msg: "missing required ID"}
		}
		h.Pedigrees = append(h.Pedigrees, d)

	case "SAMPLE":
		d := &Sample{Other: map[string]string{}}
		for _, p := range pairs {
			switch p.key {
			case "ID":
				d.ID = unquote(p.value)
			case "Genomes":
				d.Genomes = map[string]string{}
				inner, _ := parseStructuredBody(p.value)
				for _, ip := range inner {
					d.Genomes[ip.key] = unquote(ip.value)
				}
			case "MolConvs":
				d.MolConvs = map[string]string{}
				inner, _ := parseStructuredBody(p.value)
				for _, ip := range inner {
					d.MolConvs[ip.key] = unquote(ip.value)
				}
			case "Description":
				d.Description = unquote(p.value)
			case "DOI":
				d.DOI = unquote(p.value)
			default:
				d.Other[p.key] = unquote(p.value)
			}
		}
		if d.ID == "" {
			return &ParseError{Line: lineNo, Field: "SAMPLE",
				Msg: "missing required ID"}
		}
		h.Samples = append(h.Samples, d)

	default:
		// Unknown structured line — already recorded in MetaLines.
	}
	return nil
}

// fieldPair is one key=value entry extracted from a structured meta body.
type fieldPair struct {
	key, value string
}

// parseStructuredBody parses the inner text of a structured meta line — a
// comma-separated list of key=value pairs, where values may be enclosed in
// double quotes and contain escaped characters (\" → ", \\ → \).
//
// Quoting is preserved in the returned value strings (including the
// surrounding double quotes, if any); [unquote] strips them.
func parseStructuredBody(body string) ([]fieldPair, error) {
	var out []fieldPair
	var current strings.Builder
	var key string
	var seenEquals bool
	inQuotes := false
	escape := false
	fieldStart := true

	flush := func() {
		out = append(out, fieldPair{key: key, value: current.String()})
		key = ""
		current.Reset()
		seenEquals = false
		fieldStart = true
	}

	for i := 0; i < len(body); i++ {
		c := body[i]
		if escape {
			current.WriteByte(c)
			escape = false
			continue
		}
		if inQuotes {
			if c == '\\' {
				escape = true
				continue
			}
			if c == '"' {
				inQuotes = false
				current.WriteByte(c)
				continue
			}
			current.WriteByte(c)
			continue
		}
		switch c {
		case '"':
			inQuotes = true
			current.WriteByte(c)
			fieldStart = false
		case ',':
			flush()
		case '=':
			if !seenEquals {
				key = current.String()
				current.Reset()
				seenEquals = true
				fieldStart = false
			} else {
				current.WriteByte(c)
			}
		default:
			if fieldStart && (c == ' ' || c == '\t') {
				continue
			}
			fieldStart = false
			current.WriteByte(c)
		}
	}
	if inQuotes {
		return nil, fmt.Errorf("unterminated quoted value")
	}
	if current.Len() > 0 || key != "" {
		flush()
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("empty structured body")
	}
	return out, nil
}

// unquote strips surrounding double quotes from a structured-body value and
// reverses the backslash escapes the spec mandates.
//
// Inputs that are not quoted are returned verbatim; that handles values
// like "1" (unquoted Number) or "Integer" (unquoted Type).
func unquote(s string) string {
	if len(s) < 2 || s[0] != '"' || s[len(s)-1] != '"' {
		return s
	}
	inner := s[1 : len(s)-1]
	var b strings.Builder
	b.Grow(len(inner))
	escape := false
	for i := 0; i < len(inner); i++ {
		c := inner[i]
		if escape {
			switch c {
			case '"', '\\':
				b.WriteByte(c)
			default:
				b.WriteByte('\\')
				b.WriteByte(c)
			}
			escape = false
			continue
		}
		if c == '\\' {
			escape = true
			continue
		}
		b.WriteByte(c)
	}
	return b.String()
}

// parseRecord parses one tab-delimited data line using the header context
// to look up INFO/FORMAT key declarations and typed value decoding.
func parseRecord(raw []byte, lineNo int, h *Header) (*Record, error) {
	fields := bytes.Split(raw, []byte{'\t'})
	if len(fields) < 8 {
		return nil, &ParseError{Line: lineNo, Msg: fmt.Sprintf(
			"data line has %d fields; need at least 8", len(fields))}
	}
	rec := &Record{Line: lineNo, Info: map[string]InfoValue{}}

	// CHROM (col 1)
	rec.Chrom = string(fields[0])
	if rec.Chrom == "" {
		return nil, &ParseError{Line: lineNo, Field: "CHROM",
			Msg: "empty chromosome"}
	}

	// POS (col 2)
	pos, err := strconv.ParseInt(string(fields[1]), 10, 64)
	if err != nil {
		return nil, &ParseError{Line: lineNo, Field: "POS",
			Msg: "non-numeric position: " + string(fields[1])}
	}
	rec.Pos = pos

	// ID (col 3) — semicolon-separated. Per VCFv4.5 / BCF
// v2.2 the ID field must not contain whitespace, comma,
	// semicolon, or equals signs (any of these would break
	// downstream parsing).
	idStr := string(fields[2])
	if idStr != "." {
		for _, id := range strings.Split(idStr, ";") {
			if strings.ContainsAny(id, " \t,;=") {
				return nil, &ParseError{Line: lineNo, Field: "ID",
					Msg: "identifier contains whitespace, ',', ';', or '=': " + id}
			}
			rec.ID = append(rec.ID, id)
		}
	}

	// REF (col 4) — never percent-decode; the spec allows only A/C/G/T/N
	// bases here, plus the leading base for symbolic SV padding. Percent
	// encoding only applies inside INFO/FORMAT values.
	rec.Ref = string(fields[3])
	// Reduce IUPAC ambiguity codes to one of their
	// constituent bases (per VCFv4.5 §1.6.1.1 REF must be
	// A/C/G/T/N). This mutates rec.Ref in place; callers
	// that need the original sequence can re-derive it
	// from REF+ALT.
	rec.Ref = reduceIUPACBases(rec.Ref)
	if rec.Ref == "" {
		return nil, &ParseError{Line: lineNo, Field: "REF",
			Msg: "empty reference allele"}
	}

	// ALT (col 5) — comma-separated. We do not percent-decode alleles:
	// the spec lists : as a reserved character, but the spec example for
	// breakends uses literal ':' in position specifiers (e.g. "T[chrA:5[")
	// without encoding, and REF/ALT allele strings are not value strings.
	altStr := string(fields[4])
	if altStr != "." {
		rec.Alt = splitAlleles(altStr)
		// <NON_REF> is an alias for <*> (per VCFv4.5
		// §1.7, "Representing unspecified alleles") for
		// gVCF backwards compatibility. Normalise so
		// downstream code only sees the canonical form.
		for i, a := range rec.Alt {
			if a == "<NON_REF>" {
				rec.Alt[i] = "<*>"
			}
		}
		// IUPAC ambiguity codes in ALT are reduced to
		// one of their constituent bases when the ID is
		// declared via ##ALT (per VCFv4.5 §1.4.4). An
		// undeclared IUPAC code in ALT is rejected by the
		// existing symbolic-ALT declaration check
		// (Validator.checkRecord).
		for i, a := range rec.Alt {
			if id := iupacIDFromAllele(a); id != "" {
				if _, declared := h.Alts[id]; declared {
					rec.Alt[i] = reduceIUPACBases(a)
				}
			}
		}
	}

	// QUAL (col 6)
	qualStr := string(fields[5])
	if qualStr != "." {
		q, err := parseVCFFloat(qualStr)
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "QUAL",
				Msg: err.Error()}
		}
		rec.Qual = &q
	}

	// FILTER (col 7)
	filterStr := string(fields[6])
	if filterStr != "." {
		seen := make(map[string]bool, 2)
		for _, f := range strings.Split(filterStr, ";") {
			if strings.ContainsAny(f, " \t;") {
				return nil, &ParseError{Line: lineNo, Field: "FILTER",
					Msg: "filter contains whitespace or ';': " + f}
			}
			if f == "0" {
				return nil, &ParseError{Line: lineNo, Field: "FILTER",
					Msg: "'0' is reserved and must not be used as a filter String"}
			}
			if seen[f] {
				return nil, &ParseError{Line: lineNo, Field: "FILTER",
					Msg: "duplicate filter value: " + f}
			}
			seen[f] = true
			rec.Filter = append(rec.Filter, f)
		}
	}

	// INFO (col 8)
	infoStr := string(fields[7])
	if infoStr != "." {
		rec.Info, err = parseInfoColumn(infoStr, h, lineNo)
		if err != nil {
			return nil, err
		}
	}

	// Optional FORMAT (col 9) + sample columns.
	if len(fields) > 8 {
		formatStr := string(fields[8])
		if formatStr == "" {
			return nil, &ParseError{Line: lineNo, Field: "FORMAT",
				Msg: "empty FORMAT column"}
		}
		rec.Format = strings.Split(formatStr, ":")
		// GT-first is enforced by the validator, not here, so the parser
		// remains useful for inspecting malformed files.
		// We require exactly len(SampleNames) sample columns.
		wantSamples := len(h.SampleNames)
		gotSamples := len(fields) - 9
		if wantSamples != gotSamples {
			return nil, &ParseError{Line: lineNo, Field: "FORMAT",
				Msg: fmt.Sprintf("expected %d sample columns, got %d",
					wantSamples, gotSamples)}
		}
		rec.Samples = make([]Genotype, gotSamples)
		for i := 0; i < gotSamples; i++ {
			rec.Samples[i] = parseSample(fields[9+i], rec.Format, h)
		}
	} else if len(h.SampleNames) > 0 {
		return nil, &ParseError{Line: lineNo, Field: "FORMAT",
			Msg: "header declares sample columns but record has no FORMAT column"}
	}

	return rec, nil
}

// splitAlleles splits a comma-separated ALT string while respecting
// breakend brackets.
//
// Breakend replacement strings (VCFv4.5 §1.6.5) are encoded as one of:
//
//	t[p[   piece extending right of p is joined before t
//	]p]t   piece extending left of p is joined before t
//	t]p]   piece extending right of p is joined after t
//	[p[t   piece extending left of p is joined after t
//
// Each form opens and closes with the same bracket character (`[` or `]`).
// The opener enters the sequence; a comma outside any breakend is an allele
// delimiter. Commas inside a breakend (e.g. inside a contig reference) are
// preserved.
func splitAlleles(s string) []string {
	var out []string
	var current strings.Builder
	depth := 0               // 0 = outside any breakend, 1 = inside one
	expectedClose := byte(0) // the bracket that will close the current breakend

	flush := func() {
		out = append(out, current.String())
		current.Reset()
	}

	for i := 0; i < len(s); i++ {
		c := s[i]
		switch c {
		case '[', ']':
			if depth == 0 {
				depth = 1
				expectedClose = c
			} else if c == expectedClose {
				// Closing bracket — exit the breakend.
				depth = 0
				expectedClose = 0
			}
			// Any other bracket inside an open breakend is unusual; we
			// keep it as data and let downstream parsing complain if
			// needed.
			current.WriteByte(c)
		case ',':
			if depth == 0 {
				flush()
				continue
			}
			current.WriteByte(c)
		default:
			current.WriteByte(c)
		}
	}
	flush()
	return out
}

// parseVCFFloat parses a Float value, accepting "INF", "INFINITY", "NAN"
// case-insensitively in addition to the standard regular expression.
func parseVCFFloat(s string) (float64, error) {
	switch strings.ToUpper(s) {
	case "INF", "INFINITY":
		return positiveInf, nil
	case "-INF", "-INFINITY":
		return negativeInf, nil
	case "NAN":
		return nanValue, nil
	}
	return strconv.ParseFloat(s, 32)
}

// positiveInf, negativeInf, and nanValue use math.Inf / math.NaN so they
// work at compile time.
var (
	positiveInf = math.Inf(+1)
	negativeInf = math.Inf(-1)
	nanValue    = math.NaN()
)

// parseInfoColumn parses the contents of the INFO column into a typed map
// keyed by INFO ID. The declaration in h.Infos is consulted for type and
// cardinality.
func parseInfoColumn(s string, h *Header, lineNo int) (map[string]InfoValue, error) {
	out := map[string]InfoValue{}
	for _, entry := range strings.Split(s, ";") {
		if entry == "" {
			continue
		}
		key := entry
		var rawValues []string
		if eq := strings.IndexByte(entry, '='); eq >= 0 {
			key = entry[:eq]
			rawValues = strings.Split(entry[eq+1:], ",")
		}
		if !validInfoKey(key) {
			return nil, &ParseError{Line: lineNo, Field: "INFO",
				Msg: "invalid key: " + key}
		}
		def, declared := h.LookupInfo(key)
		if !declared {
			// Unknown INFO keys are permitted by the spec but with
			// TypeUnknown. We still store the raw values.
			v := InfoValue{
				typedCore: typedCore{Type: TypeUnknown, Raw: rawValues},
			}
			if len(rawValues) == 0 {
				v.Flag = true
			}
			out[key] = v
			continue
		}
		v, err := decodeInfoValue(def, rawValues)
		if err != nil {
			return nil, &ParseError{Line: lineNo, Field: "INFO/" + key,
				Msg: err.Error()}
		}
		out[key] = v
	}
	return out, nil
}

// contigNamePattern is the SAM-compatible contig-name
// regex from VCFv4.5 §1.4.7:
//
//	[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*
//
// The first character excludes '*' and '='. Subsequent
// characters may include them. ',' is excluded entirely.
var contigNamePattern = `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`

// validContigName reports whether s matches contigNamePattern.
// iupacReduction maps each uppercase IUPAC ambiguity code
// to one of its constituent A/C/G/T bases. Lower case input
// is reduced by folding case and then copying the result
// back to the input case. Per VCFv4.5 §1.6.1.1 REF must be
// only A/C/G/T/N; we reduce IUPAC codes to satisfy the rule
// without rejecting real-world records that contain them.
var iupacReduction = map[byte]byte{
	'R': 'A', // A or G
	'Y': 'C', // C or T
	'S': 'C', // G or C
	'W': 'A', // A or T
	'K': 'G', // G or T
	'M': 'A', // A or C
	'B': 'C', // C, G or T
	'D': 'A', // A, G or T
	'H': 'A', // A, C or T
	'V': 'A', // A, C or G
}

// reduceIUPACBases rewrites each IUPAC ambiguity code in s
// to one of its constituent bases, preserving the case of
// the input. A, C, G, T, N (any case) are returned
// unchanged.
func reduceIUPACBases(s string) string {
	if !strings.ContainsAny(s, "RYSWKMBDHVryswkmbdhv") {
		return s
	}
	out := []byte(s)
	for i, c := range out {
		// Fold ASCII case without going through rune to
		// keep gosec happy.
		var upper byte
		switch {
		case c >= 'a' && c <= 'z':
			upper = c - ('a' - 'A')
		default:
			upper = c
		}
		base, ok := iupacReduction[upper]
		if !ok {
			continue
		}
		if c >= 'a' && c <= 'z' {
			base += 'a' - 'A'
		}
		out[i] = base
	}
	return string(out)
}

// iupacIDs lists the IUPAC ambiguity codes the spec
// treats as eligible symbolic-allele IDs (§1.4.4).
var iupacIDs = map[string]bool{
	"R": true, "Y": true, "S": true, "W": true, "K": true,
	"M": true, "B": true, "D": true, "H": true, "V": true,
}

// iupacIDFromAllele returns the IUPAC code ID when a is
// a bare IUPAC allele (no surrounding angle brackets).
// Returns "" otherwise.
func iupacIDFromAllele(a string) string {
	if !iupacIDs[a] {
		return ""
	}
	return a
}

// validMD5Hex reports whether s is 32 lowercase hex characters,
// the format VCFv4.5 §1.4.7 prescribes for ##contig MD5. We
// do not recompute the digest because the sequence is not in
// the VCF header.
func validMD5Hex(s string) bool {
	if len(s) != 32 {
		return false
	}
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c < '0' || (c > '9' && c < 'a') || c > 'f' {
			return false
		}
	}
	return true
}

// reservedSymbolicAlleles lists the symbolic ALT IDs the
// spec treats as reserved. A ##contig ID must not match
// one of these names (VCFv4.5 §1.4.7).
var reservedSymbolicAlleles = map[string]bool{
	"DEL": true, "INS": true, "DUP": true, "INV": true,
	"CNV": true, "NON_REF": true,
	"DEL:ME": true, "INS:ME": true,
	"DUP:TANDEM": true, "DUP:DISPERSED": true, "CNV:TR": true,
}

// parseContigAttrs parses a ##contig structured line's
// key/value pairs and validates them per VCFv4.5 §1.4.7.
// Returns the Contig or an error.
func parseContigAttrs(pairs []fieldPair, lineNo int) (*Contig, error) {
	c := &Contig{Other: map[string]string{}}
	for _, p := range pairs {
		switch p.key {
		case "ID":
			c.ID = unquote(p.value)
		case "length":
			n, err := strconv.ParseInt(p.value, 10, 64)
			if err != nil {
				return nil, &ParseError{Line: lineNo, Field: "contig/length",
					Msg: "invalid length: " + p.value}
			}
			c.Length = n
		case "md5":
			if err := setContigMD5(c, p.value, lineNo); err != nil {
				return nil, err
			}
		case "URL":
			c.URL = unquote(p.value)
		case "species":
			c.Species = unquote(p.value)
		case "taxonomy":
			c.Taxonomy = unquote(p.value)
		default:
			c.Other[p.key] = unquote(p.value)
		}
	}
	if c.ID == "" {
		return nil, &ParseError{Line: lineNo, Field: "contig",
			Msg: "missing required ID"}
	}
	// Contig name regex (SAM-compatible per VCFv4.5
	// §1.4.7 "Contig field format"). The first character
	// excludes '*' and '='; subsequent characters may
	// include them.
	if !validContigName(c.ID) {
		return nil, &ParseError{Line: lineNo, Field: "contig/ID",
			Msg: fmt.Sprintf("contig name %q does not match the required regex %s",
				c.ID, contigNamePattern)}
	}
	// Contig ID must not be a reserved symbolic allele
	// name (VCFv4.5 §1.4.7).
	if reservedSymbolicAlleles[c.ID] {
		return nil, &ParseError{Line: lineNo, Field: "contig/ID",
			Msg: fmt.Sprintf("contig name %q is a reserved symbolic allele name and cannot be used as a contig ID", c.ID)}
	}
	return c, nil
}

// setContigMD5 parses and validates a ##contig MD5 attribute
// value. The actual digest is not recomputed (the sequence is
// not in the header) but the 32 lowercase hex format is
// checked.
func setContigMD5(c *Contig, raw string, lineNo int) error {
	c.MD5 = unquote(raw)
	if !validMD5Hex(c.MD5) {
		return &ParseError{Line: lineNo, Field: "contig/md5",
			Msg: fmt.Sprintf("contig MD5 %q must be 32 lowercase hex characters", c.MD5)}
	}
	return nil
}

func validContigName(s string) bool {
	if s == "" {
		return false
	}
	if !isContigNameStart(s[0]) {
		return false
	}
	for i := 1; i < len(s); i++ {
		if !isContigNameCont(s[i]) {
			return false
		}
	}
	return true
}

func isContigNameStart(c byte) bool {
	return (c >= '0' && c <= '9') ||
		(c >= 'A' && c <= 'Z') ||
		(c >= 'a' && c <= 'z') ||
		c == '!' || c == '#' || c == '$' || c == '%' || c == '&' ||
		c == '+' || c == '.' || c == '/' || c == ':' || c == ';' ||
		c == '?' || c == '@' || c == '^' || c == '_' || c == '|' ||
		c == '~' || c == '-'
}

func isContigNameCont(c byte) bool {
	return isContigNameStart(c) || c == '*' || c == '='
}

// validInfoKey matches the regex ^([A-Za-z_][0-9A-Za-z_.]*|1000G)$ from
// the spec.
func validInfoKey(s string) bool {
	if s == "1000G" {
		return true
	}
	if s == "" {
		return false
	}
	if !isKeyStart(s[0]) {
		return false
	}
	for i := 1; i < len(s); i++ {
		if !isKeyCont(s[i]) {
			return false
		}
	}
	return true
}

func isKeyStart(c byte) bool {
	return c == '_' || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')
}

func isKeyCont(c byte) bool {
	return c == '_' || c == '.' ||
		(c >= '0' && c <= '9') ||
		(c >= 'A' && c <= 'Z') ||
		(c >= 'a' && c <= 'z')
}

// decodeTyped populates the typed slices on c from raw according to c.Type.
// It handles the "." missing sentinel and the Integer / Float / String
// cases. Flag, Character, and Unknown have type-specific semantics that
// callers ([decodeInfoValue], [decodeFieldValue]) handle themselves.
//
// Decoding is lenient: a token that fails to parse (e.g. "abc" for an
// Integer field) is recorded as a zero in the typed slice while the raw
// token remains in c.Raw for the validator to surface. The parser never
// rejects on type errors so semantically-broken VCFs can still be
// inspected.
func decodeTyped(c *typedCore, raw []string) {
	if len(raw) == 1 && raw[0] == "." {
		c.Missing = true
		return
	}
	switch c.Type {
	case TypeInteger:
		for _, r := range raw {
			if r == "." {
				c.Ints = append(c.Ints, 0)
				continue
			}
			n, err := strconv.ParseInt(r, 10, 32)
			if err != nil {
				// Lenient: keep raw but mark type mismatch for the
				// validator to surface.
				c.Ints = append(c.Ints, 0)
				continue
			}
			c.Ints = append(c.Ints, n)
		}
	case TypeFloat:
		for _, r := range raw {
			if r == "." {
				c.Floats = append(c.Floats, 0)
				continue
			}
			f, err := parseVCFFloat(r)
			if err != nil {
				c.Floats = append(c.Floats, 0)
				continue
			}
			c.Floats = append(c.Floats, f)
		}
	case TypeString:
		for _, r := range raw {
			c.Strs = append(c.Strs, percentDecode(r))
		}
	}
}

// decodeInfoValue converts the raw string values for one INFO entry to
// typed slices based on the declared FieldType. Flag entries must have
// no values (Number == 0); other entries may have any number.
func decodeInfoValue(def *InfoDef, raw []string) (InfoValue, error) {
	v := InfoValue{typedCore: typedCore{Type: def.Type, Raw: raw}}

	if def.Type == TypeFlag {
		if len(raw) > 0 {
			return v, fmt.Errorf("flag field must not have values; got %d",
				len(raw))
		}
		v.Flag = true
		return v, nil
	}

	decodeTyped(&v.typedCore, raw)

	// INFO Character values are percent-decoded per the spec; FORMAT
	// Character values are not (see decodeFieldValue). The "." case is
	// already handled by decodeTyped setting Missing.
	if def.Type == TypeCharacter && len(raw) == 1 {
		s := percentDecode(raw[0])
		if len(s) == 1 {
			v.Chars = []byte{s[0]}
		}
	}
	return v, nil
}

// parseSample splits a sample column into per-FORMAT-key FieldValues.
func parseSample(raw []byte, format []string, h *Header) Genotype {
	fields := bytes.Split(raw, []byte{':'})
	out := Genotype{Fields: map[string]FieldValue{}}

	// Trailing fields may be dropped (except GT, when present). We do not add
	// dropped fields to the Fields map — the validator distinguishes "dropped"
	// from "explicitly missing" by checking whether the FORMAT key is past the
	// last present field.
	n := min(len(fields), len(format))
	for i := 0; i < n; i++ {
		key := format[i]
		rawVals := strings.Split(string(fields[i]), ",")
		def, declared := h.LookupFormat(key)
		var v FieldValue
		if !declared {
			v = FieldValue{typedCore: typedCore{Type: TypeUnknown, Raw: rawVals}}
		} else {
			v = decodeFieldValue(def, rawVals)
		}
		out.Fields[key] = v
	}
	return out
}

func decodeFieldValue(def *FormatDef, raw []string) FieldValue {
	v := FieldValue{typedCore: typedCore{Type: def.Type, Raw: raw}}

	decodeTyped(&v.typedCore, raw)

	// FORMAT Character values are NOT percent-decoded (see decodeInfoValue
	// for the INFO counterpart). The "." case is already handled by
	// decodeTyped setting Missing, so only the single-char branch remains.
	if def.Type == TypeCharacter && len(raw) == 1 && len(raw[0]) == 1 {
		v.Chars = []byte{raw[0][0]}
	}
	return v
}
