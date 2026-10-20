package bcf

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strings"

	"github.com/inful/htsfile/internal/bgzf"
)

// Magic is the six-byte BCF file signature. The first
// three bytes are the ASCII "BCF". Byte 3 is the major
// version (2). Byte 4 is the minor version (1 for BCF 2.1,
// 2 for BCF 2.2). Byte 5 is reserved and must be zero.
var Magic = [6]byte{'B', 'C', 'F', 0x02, 0x01, 0x00}

// magic22 is the BCF 2.2 signature. The decoder accepts
// either 2.1 or 2.2. The encoder emits 2.1.
var magic22 = [6]byte{'B', 'C', 'F', 0x02, 0x02, 0x00}

// Decoder reads a BCF v2.1 document incrementally from a
// BGZF-compressed [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("sample.bcf")
//	defer f.Close()
//	dec := bcf.NewDecoder(f)
//	hdr, _ := dec.Header()
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	}
//
// A Decoder is not safe for concurrent use.
type Decoder struct {
	bgzf *bgzf.Reader
	br   *bufio.Reader

	hdrRead bool
	hdr     *Header
	hdrErr  error

	// version is the minor-version byte from the BCF magic
	// (1 for BCF 2.1, 2 for BCF 2.2). Zero before Header()
	// is called.
	version uint8

	recordLine int
}

// NewDecoder wraps r in a [Decoder]. r is expected to be a
// BGZF-compressed BCF stream.
func NewDecoder(r io.Reader) *Decoder {
	gz := bgzf.NewReader(r)
	return &Decoder{
		bgzf: gz,
		br:   bufio.NewReader(gz),
	}
}

// Version returns the BCF minor-version byte from the magic
// (1 for BCF 2.1, 2 for BCF 2.2). Returns 0 before Header()
// has been called.
func (d *Decoder) Version() uint8 {
	return d.version
}

// Header reads the BCF magic + length-prefixed VCF-style
// header text + parses the header into a [Header] value. The
// first call performs the read; subsequent calls return the
// cached value.
//
// If Record is called before Header, Header is invoked
// implicitly.
func (d *Decoder) Header() (*Header, error) {
	if d.hdrRead {
		return d.hdr, d.hdrErr
	}
	d.hdrRead = true
	// Magic check.
	var magic [6]byte
	if _, err := io.ReadFull(d.br, magic[:]); err != nil {
		if errors.Is(err, io.EOF) || errors.Is(err, io.ErrUnexpectedEOF) {
			d.hdrErr = &ParseError{Field: "magic", Msg: "BCF stream too short for magic"}
			return nil, d.hdrErr
		}
		d.hdrErr = err
		return nil, err
	}
	if magic != Magic && magic != magic22 {
		d.hdrErr = &ParseError{Field: "magic",
			Msg: fmt.Sprintf("not a BCF file: got %q", magic[:])}
		return nil, d.hdrErr
	}
	d.version = magic[4]
	// Length-prefixed header text.
	var lText uint32
	if err := binary.Read(d.br, binary.LittleEndian, &lText); err != nil {
		d.hdrErr = &ParseError{Field: "header",
			Msg: "failed to read l_text: " + err.Error()}
		return nil, d.hdrErr
	}
	textBytes := make([]byte, lText)
	if _, err := io.ReadFull(d.br, textBytes); err != nil {
		d.hdrErr = &ParseError{Field: "header",
			Msg: "failed to read header text: " + err.Error()}
		return nil, d.hdrErr
	}
	d.hdr = ParseHeader(string(textBytes))
	return d.hdr, nil
}

// Record returns the next BCF record, or (nil, [io.EOF]) when
// no more records are available. Calling Record after EOF
// continues to return (nil, [io.EOF]).
//
// v2 supports the shared per-record fields (CHROM, POS, ID,
// REF, ALT, FILTER, QUAL, INFO) and the per-sample FORMAT
// block. Records with n_fmt == 0 (no per-sample data) round
// trip without a FORMAT section.
func (d *Decoder) Record() (*Record, error) {
	if !d.hdrRead {
		if _, err := d.Header(); err != nil {
			return nil, err
		}
	}
	d.recordLine++
	// CHROM (int32) + POS (int32) + rlen (int32) + n_info (uint32)
	var prefix [16]byte
	if _, err := io.ReadFull(d.br, prefix[:]); err != nil {
		if errors.Is(err, io.EOF) || errors.Is(err, io.ErrUnexpectedEOF) {
			return nil, io.EOF
		}
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: "failed to read record prefix: " + err.Error()}
	}
	chromIdx := int32(binary.LittleEndian.Uint32(prefix[0:4]))
	pos0 := int32(binary.LittleEndian.Uint32(prefix[4:8]))
	_ = int32(binary.LittleEndian.Uint32(prefix[8:12])) // rlen (unused)
	nInfo := binary.LittleEndian.Uint32(prefix[12:16])
	// n_fmt (uint16) — we read it but skip the per-sample
	// data entirely (v1 limitation).
	var nFmtBytes [2]byte
	if _, err := io.ReadFull(d.br, nFmtBytes[:]); err != nil {
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: "failed to read n_fmt: " + err.Error()}
	}
	nFmt := binary.LittleEndian.Uint16(nFmtBytes[:])
	_ = nFmt // consumed later in the per-sample block
	// ID (1 byte type + optional 4-byte length + bytes).
	id, err := d.readTypedString()
	if err != nil {
		return nil, err
	}
	// REF (same shape as ID).
	ref, err := d.readTypedString()
	if err != nil {
		return nil, err
	}
	// ALT: a sequence of typed-string entries (each with
	// type byte 7 or, in BCF 2.2, 0xF for long strings).
	// The list is terminated by the first non-string type
	// byte, which begins the FILTER field. An ALT count of
	// zero is implicit (no typed-string bytes appear between
	// REF and FILTER).
	var alt []string
	for {
		peek, err := d.peekByte()
		if err != nil {
			return nil, err
		}
		if peek != typeString && peek != typeLong {
			break
		}
		a, err := d.readTypedString()
		if err != nil {
			return nil, err
		}
		alt = append(alt, a)
	}
	// FILTER (typed: missing | single int32 | array of int32).
	filters, err := d.readTypedFilter()
	if err != nil {
		return nil, err
	}
	// QUAL (typed: missing | float64).
	qual, qualKnown, err := d.readTypedQual()
	if err != nil {
		return nil, err
	}
	// INFO: n_info typed entries. Each is keyed by a uint16
	// header index plus a typed value.
	info, err := d.readTypedInfo(nInfo)
	if err != nil {
		return nil, err
	}
	rec := &Record{
		Line:      d.recordLine,
		CHROM:     d.hdr.LookupContig(int(chromIdx)),
		Pos:       int64(pos0) + 1,
		ID:        id,
		Ref:       ref,
		Alt:       alt,
		Qual:      qual,
		QualKnown: qualKnown,
		Filter:    filters,
		Info:      info,
	}
	// Per-sample FORMAT block.
	if nFmt > 0 {
		formatKeys, err := d.readFormatKeys(nFmt)
		if err != nil {
			return nil, err
		}
		rec.Format = formatKeys
		sampleNames := d.hdr.SampleNames
		rec.Samples = make([]Sample, len(sampleNames))
		for i := range sampleNames {
			s, err := d.readSample(formatKeys)
			if err != nil {
				return nil, err
			}
			rec.Samples[i] = *s
		}
	}
	return rec, nil
}

// readFormatKeys reads the n_fmt uint16 header indices that
// form the per-record FORMAT list.
func (d *Decoder) readFormatKeys(nFmt uint16) ([]uint16, error) {
	keys := make([]uint16, nFmt)
	for i := uint16(0); i < nFmt; i++ {
		var b [2]byte
		if _, err := io.ReadFull(d.br, b[:]); err != nil {
			return nil, fmt.Errorf("bcf: read FORMAT[%d]: %w", i, err)
		}
		keys[i] = binary.LittleEndian.Uint16(b[:])
	}
	return keys, nil
}

// readTypedString reads a type byte + optional length-prefixed
// string body. Returns "" when the field is missing.
//
// BCF 2.1 strings use typeString (7) + uint32 length + bytes.
// BCF 2.2 adds typeLong (0xF) as an equivalent form, used by
// the encoder for strings of length 15 or more. The decoder
// accepts either.
func (d *Decoder) readTypedString() (string, error) {
	t, err := d.readByte()
	if err != nil {
		return "", err
	}
	switch t {
	case typeMissing:
		return "", nil
	case typeString, typeLong:
		var l uint32
		if err := binary.Read(d.br, binary.LittleEndian, &l); err != nil {
			return "", err
		}
		if l == bcfInt32Missing {
			return "", nil
		}
		buf := make([]byte, l)
		if _, err := io.ReadFull(d.br, buf); err != nil {
			return "", err
		}
		return string(buf), nil
	default:
		return "", &ParseError{Field: "string",
			Msg: fmt.Sprintf("unexpected type byte %d for string field", t)}
	}
}

// readTypedFilter reads the FILTER typed array. BCF FILTER
// stores filter indices (int32) keyed by header.Filters.
// Index 0 traditionally means "PASS" — the resolver in
// the header returns "PASS" for index 0 when no filters
// were declared.
func (d *Decoder) readTypedFilter() ([]string, error) {
	t, err := d.readByte()
	if err != nil {
		return nil, err
	}
	switch t {
	case typeMissing:
		return nil, nil
	case typeInt8:
		// Spec says "single int8" — but the filter list is
		// indexed by int32; in practice the type is 1 (int8
		// per the type-tag spec) for the single-filter case
		// and the int32 follows. We treat type 1 as "single
		// int32" for compatibility with real-world files.
		var idx uint32
		if err := binary.Read(d.br, binary.LittleEndian, &idx); err != nil {
			return nil, err
		}
		if idx == bcfFilterMissing {
			return nil, nil
		}
		if idx == 0 && d.hdr != nil && len(d.hdr.Filters) == 0 {
			// Traditional convention: no ##FILTER lines
			// means index 0 is PASS.
			return []string{"PASS"}, nil
		}
		name := d.hdr.LookupFilter(int(idx))
		if name == "" {
			return []string{fmt.Sprintf("FILTER#%d", idx)}, nil
		}
		return []string{name}, nil
	case typeInt32:
		// Array form: int32 length + that many int32 indices.
		var n uint32
		if err := binary.Read(d.br, binary.LittleEndian, &n); err != nil {
			return nil, err
		}
		if n == 0 {
			return nil, nil
		}
		out := make([]string, 0, n)
		for i := uint32(0); i < n; i++ {
			var idx uint32
			if err := binary.Read(d.br, binary.LittleEndian, &idx); err != nil {
				return nil, err
			}
			if idx == bcfFilterMissing {
				continue
			}
			if idx == 0 && d.hdr != nil && len(d.hdr.Filters) == 0 {
				out = append(out, "PASS")
				continue
			}
			name := d.hdr.LookupFilter(int(idx))
			if name == "" {
				out = append(out, fmt.Sprintf("FILTER#%d", idx))
				continue
			}
			out = append(out, name)
		}
		return out, nil
	default:
		return nil, &ParseError{Field: "FILTER",
			Msg: fmt.Sprintf("unexpected FILTER type byte %d", t)}
	}
}

// readTypedQual reads the QUAL typed value (float64).
func (d *Decoder) readTypedQual() (qual float64, known bool, err error) {
	t, err := d.readByte()
	if err != nil {
		return 0, false, err
	}
	switch t {
	case typeMissing:
		return 0, false, nil
	case typeFloat:
		var bits uint64
		if err := binary.Read(d.br, binary.LittleEndian, &bits); err != nil {
			return 0, false, err
		}
		return math.Float64frombits(bits), true, nil
	default:
		return 0, false, &ParseError{Field: "QUAL",
			Msg: fmt.Sprintf("unexpected QUAL type byte %d", t)}
	}
}

// readTypedInfo reads the n_info typed INFO entries and
// returns the resolved Info map.
func (d *Decoder) readTypedInfo(n uint32) (map[string]InfoValue, error) {
	out := make(map[string]InfoValue)
	for i := uint32(0); i < n; i++ {
		// Key: 2 bytes uint16 = index into the header's
		// ordered INFO declarations. We need an ordered
		// list — Header.Infos is a map; we look up by
		// iterating to find the i-th key.
		var keyBytes [2]byte
		if _, err := io.ReadFull(d.br, keyBytes[:]); err != nil {
			return nil, err
		}
		keyIdx := int(binary.LittleEndian.Uint16(keyBytes[:]))
		def, ok := d.hdr.InfoAt(keyIdx)
		if !ok {
			// Unknown INFO key — read past the typed value
			// so we stay aligned, but skip the assignment.
			if err := d.skipTypedInfoValue(); err != nil {
				return nil, err
			}
			continue
		}
		val, err := d.readTypedInfoValue(def)
		if err != nil {
			return nil, err
		}
		out[def.ID] = val
	}
	return out, nil
}

// readTypedInfoValue reads one typed INFO body. The format
// varies by declared type:
//
//   - Flag: one type byte; typeMissing = absent, anything
//     else = present. No length prefix, no body.
//   - All other types: optional length prefix (for A/R/G/.
//     Number modes) followed by a value type byte and count
//     typed values.
//
// Per the BCF v2.1 spec, the length-prefix's own type byte
// is *separate* from the value's type byte — the decoder
// must consume the length-prefix type byte before reading
// the value type byte.
func (d *Decoder) readTypedInfoValue(def *InfoDef) (InfoValue, error) {
	//nolint:dupl // not actually duplicate; S1006 placeholder
	v := InfoValue{Type: def.Type}

	// Special case: Flag has no length prefix; just one type byte.
	if def.Type == "Flag" {
		t, err := d.readByte()
		if err != nil {
			return v, err
		}
		if t != typeMissing {
			v.Flag = true
		}
		return v, nil
	}

	// All other types: optional length prefix, then a value
	// type byte, then the values.
	count, isVar := ParseNumber(def.Number)
	if isVar {
		// Length prefix: 1 type byte + 4-byte uint32 count.
		// The type byte is conventionally typeInt32 (1).
		if _, err := d.readByte(); err != nil {
			return v, err
		}
		var l uint32
		if err := binary.Read(d.br, binary.LittleEndian, &l); err != nil {
			return v, err
		}
		count = int(l)
	}

	// Value type byte.
	t, err := d.readByte()
	if err != nil {
		return v, err
	}
	if t == typeMissing {
		return v, nil
	}
	if count == 0 {
		return v, nil
	}

	switch def.Type {
	case "Integer":
		v.Ints = make([]int64, count)
		switch t & 0x0F { // mask the high nibble (inline size)
		case typeInt8:
			for i := 0; i < count; i++ {
				b, err := d.readByte()
				if err != nil {
					return InfoValue{}, err
				}
				v.Ints[i] = int64(int8(b))
			}
		case typeInt16:
			for i := 0; i < count; i++ {
				var x uint16
				if err := binary.Read(d.br, binary.LittleEndian, &x); err != nil {
					return InfoValue{}, err
				}
				v.Ints[i] = int64(int16(x))
			}
		case typeInt32:
			for i := 0; i < count; i++ {
				var x uint32
				if err := binary.Read(d.br, binary.LittleEndian, &x); err != nil {
					return InfoValue{}, err
				}
				v.Ints[i] = int64(int32(x))
			}
		default:
			return InfoValue{}, &ParseError{Field: "INFO",
				Msg: fmt.Sprintf("Integer field has type byte %d", t)}
		}
	case "Float":
		v.Floats = make([]float64, count)
		switch t & 0x0F { // mask the high nibble (inline size)
		case typeFloat:
			for i := 0; i < count; i++ {
				var bits uint32
				if err := binary.Read(d.br, binary.LittleEndian, &bits); err != nil {
					return InfoValue{}, err
				}
				v.Floats[i] = float64(math.Float32frombits(bits))
			}
		default:
			return InfoValue{}, &ParseError{Field: "INFO",
				Msg: fmt.Sprintf("Float field has type byte %d", t)}
		}
	case "Character":
		if count != 1 {
			return InfoValue{}, &ParseError{Field: "INFO",
				Msg: "Character field must have Number=1"}
		}
		b, err := d.readByte()
		if err != nil {
			return InfoValue{}, err
		}
		v.Str = string(rune(b))
	case "String":
		v.Strs = make([]string, count)
		for i := 0; i < count; i++ {
			var l uint32
			if err := binary.Read(d.br, binary.LittleEndian, &l); err != nil {
				return InfoValue{}, err
			}
			buf := make([]byte, l)
			if _, err := io.ReadFull(d.br, buf); err != nil {
				return InfoValue{}, err
			}
			v.Strs[i] = string(buf)
		}
	default:
		return InfoValue{}, &ParseError{Field: "INFO",
			Msg: "unknown declared type " + def.Type}
	}
	return v, nil
}

// skipTypedInfoValue reads past one typed INFO body (used
// when the key index points to an unknown declaration).
func (d *Decoder) skipTypedInfoValue() error {
	t, err := d.readByte()
	if err != nil {
		return err
	}
	if t == typeMissing {
		return nil
	}
	// Variable-length prefix.
	var l uint32
	if err := binary.Read(d.br, binary.LittleEndian, &l); err != nil {
		return err
	}
	// Skip l typed elements. We don't know the declared
	// type, so we make a best-effort guess from the type
	// byte: int8=1, int16=2, int32/float=4, string=variable.
	var elemSize int
	switch t {
	case typeInt8:
		elemSize = 1
	case typeInt16:
		elemSize = 2
	case typeInt32, typeFloat:
		elemSize = 4
	default:
		// Unknown — best effort: read 1 byte per element.
		elemSize = 1
	}
	total := int(l) * elemSize
	buf := make([]byte, total)
	if _, err := io.ReadFull(d.br, buf); err != nil {
		return err
	}
	return nil
}

// skipFmtBlock is no longer used; the v2 FORMAT reader
// consumes the per-sample block above. Kept for backwards
// compatibility with any external callers.
func (d *Decoder) skipFmtBlock() error {
	return &ParseError{Field: "FORMAT",
		Msg: "skipFmtBlock is no longer used; v2 reads the per-sample block directly"}
}

// peekByte reads one byte without consuming it from the
// stream.
func (d *Decoder) peekByte() (byte, error) {
	b, err := d.br.Peek(1)
	if err != nil {
		return 0, err
	}
	if len(b) == 0 {
		return 0, io.EOF
	}
	return b[0], nil
}

// readByte reads one byte (helper for the Decoder).
func (d *Decoder) readByte() (byte, error) {
	var b [1]byte
	if _, err := io.ReadFull(d.br, b[:]); err != nil {
		return 0, err
	}
	return b[0], nil
}

// InfoAt returns the n-th INFO declaration in source order.
// Used by the Decoder to map BCF key indices (which are
// positionally ordered, not named) back to declarations.
// We expose this on Header because the BCF record's key
// index is a position, not a name.
func (h *Header) InfoAt(idx int) (*InfoDef, bool) {
	if h == nil || idx < 0 || idx >= len(h.InfosOrder) {
		return nil, false
	}
	id := h.InfosOrder[idx]
	def, ok := h.Infos[id]
	return def, ok
}

// LookupInfoIndex returns the positional index of id in the
// header's ##INFO declaration order. BCF records store
// INFO keys as uint16 indices into that order, not as names;
// a writer (future) and the test builder (today) must use
// this index to avoid the off-by-N bug where a local
// declaration order differs from the header's.
//
// Returns (-1, false) if id is not declared.
func (h *Header) LookupInfoIndex(id string) (int, bool) {
	if h == nil {
		return -1, false
	}
	for i, declaredID := range h.InfosOrder {
		if declaredID == id {
			return i, true
		}
	}
	return -1, false
}

// FormatAt returns the n-th FORMAT declaration in source order.
func (h *Header) FormatAt(idx int) (*FormatDef, bool) {
	if h == nil || idx < 0 || idx >= len(h.FormatsOrder) {
		return nil, false
	}
	id := h.FormatsOrder[idx]
	def, ok := h.Formats[id]
	return def, ok
}

// LookupFormatIndex returns the positional index of id in the
// header's ##FORMAT declaration order. See [Header.LookupInfoIndex]
// for the rationale.
func (h *Header) LookupFormatIndex(id string) (int, bool) {
	if h == nil {
		return -1, false
	}
	for i, declaredID := range h.FormatsOrder {
		if declaredID == id {
			return i, true
		}
	}
	return -1, false
}

// _ = strings.HasPrefix keeps the strings import visible
// (it is used by parseHeader).
var _ = strings.HasPrefix
