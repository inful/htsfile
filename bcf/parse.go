package bcf

import (
	"bufio"
	"bytes"
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

	// pending holds bytes that were read ahead of the
	// structured parse (l_shared / l_indiv blocks are read
	// in one io.ReadFull, then handed back to the typed
	// helpers). ReadByte and peekByte drain it first.
	pending []byte
}

// pushBack returns b to the front of the read stream. Used
// after a bulk io.ReadFull of a length-delimited block, so
// the typed-value helpers can consume it field by field.
func (d *Decoder) pushBack(b []byte) {
	if len(b) == 0 {
		return
	}
	d.pending = append(append([]byte(nil), b...), d.pending...)
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
	if err := d.readLe(&lText); err != nil {
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
	// l_shared + l_indiv (uint32 each), spec BcfSiteEncoding.
	var lens [8]byte
	if _, err := io.ReadFull(d.br, lens[:]); err != nil {
		if errors.Is(err, io.EOF) || errors.Is(err, io.ErrUnexpectedEOF) {
			return nil, io.EOF
		}
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: "failed to read record lengths: " + err.Error()}
	}
	lShared := binary.LittleEndian.Uint32(lens[0:4])
	lIndiv := binary.LittleEndian.Uint32(lens[4:8])

	// 24-byte site prefix: CHROM, POS, rlen, QUAL, n_info,
	// n_allele, n_sample (24 bits), n_fmt (8 bits).
	shared := make([]byte, lShared)
	if _, err := io.ReadFull(d.br, shared); err != nil {
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: "failed to read record shared block: " + err.Error()}
	}
	if len(shared) < 24 {
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: fmt.Sprintf("l_shared=%d too short for the 24-byte site prefix", lShared)}
	}
	chromIdx := int32(binary.LittleEndian.Uint32(shared[0:4]))
	pos0 := int32(binary.LittleEndian.Uint32(shared[4:8]))
	// rlen at [8:12] is informational; the decoder does not need it.
	qualBits := binary.LittleEndian.Uint32(shared[12:16])
	nInfo := binary.LittleEndian.Uint16(shared[16:18])
	nAllele := binary.LittleEndian.Uint16(shared[18:20])
	nFmt := uint8(binary.LittleEndian.Uint32(shared[20:24]))

	// QUAL is a plain float32; the per-type MISSING sentinel
	// 0x7F800001 means "unknown".
	var qual float64
	qualKnown := qualBits != bcfFloatMissing
	if qualKnown {
		qual = float64(math.Float32frombits(qualBits))
	}

	// Rest of the shared block, read through the buffered
	// reader so the typed helpers stay in one place.
	// Hand the remainder of the shared block back so the
	// typed helpers below consume it field by field.
	d.pushBack(shared[24:])

	// ID (typed string).
	id, err := d.readTypedString()
	if err != nil {
		return nil, err
	}
	// REF is the first of n_allele typed strings; the rest
	// are ALT.
	if nAllele == 0 {
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: "n_allele=0 but REF is mandatory"}
	}
	ref, err := d.readTypedString()
	if err != nil {
		return nil, err
	}
	alt := make([]string, 0, int(nAllele)-1)
	for i := 1; i < int(nAllele); i++ {
		a, err := d.readTypedString()
		if err != nil {
			return nil, err
		}
		alt = append(alt, a)
	}
	// FILTER.
	filters, err := d.readTypedFilter()
	if err != nil {
		return nil, err
	}
	// INFO.
	info, err := d.readTypedInfo(uint32(nInfo))
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

	// Per-sample block (l_indiv bytes).
	if lIndiv > 0 && nFmt > 0 {
		indiv := make([]byte, lIndiv)
		if _, err := io.ReadFull(d.br, indiv); err != nil {
			return nil, &ParseError{Line: d.recordLine, Field: "record",
				Msg: "failed to read per-sample block: " + err.Error()}
		}
		d.pushBack(indiv)
		formatKeys, err := d.readFormatKeys(nFmt)
		if err != nil {
			return nil, err
		}
		rec.Format = formatKeys
		samples, err := d.readFormatBlock(formatKeys, len(d.hdr.SampleNames))
		if err != nil {
			return nil, err
		}
		rec.Samples = samples
	}
	return rec, nil
}

// readFormatKeys reads the n_fmt typed atomic integers that
// form the per-record FORMAT key list (spec section 5,
// GenotypeEncoding: "fmt_key typed int").
func (d *Decoder) readFormatKeys(nFmt uint8) ([]uint16, error) {
	keys := make([]uint16, nFmt)
	for i := uint8(0); i < nFmt; i++ {
		v, err := d.readTypedInt()
		if err != nil {
			return nil, fmt.Errorf("bcf: read FORMAT[%d]: %w", i, err)
		}
		keys[i] = uint16(v)
	}
	return keys, nil
}

// readLe reads a little-endian value of any fixed width,
// draining the pushback buffer first.
func (d *Decoder) readLe(v any) error {
	sz := binary.Size(v)
	if sz <= 0 {
		return nil
	}
	buf := make([]byte, sz)
	if err := d.readPending(buf); err != nil {
		return err
	}
	return binary.Read(bytes.NewReader(buf), binary.LittleEndian, v)
}

// readTypedInt reads one typed atomic integer: a type byte
// (any of int8/int16/int32) followed by the little-endian
// value. Used for the INFO dictionary key and the FORMAT key
// list.
func (d *Decoder) readTypedInt() (int32, error) {
	t, err := d.readByte()
	if err != nil {
		return 0, err
	}
	if t == typeMissing {
		return 0, nil
	}
	switch t & 0x0F {
	case typeInt8:
		b, err := d.readByte()
		if err != nil {
			return 0, err
		}
		return int32(int8(b)), nil
	case typeInt16:
		var b [2]byte
		if err := d.readPending(b[:]); err != nil {
			return 0, err
		}
		return int32(int16(binary.LittleEndian.Uint16(b[:]))), nil
	case typeInt32:
		var b [4]byte
		if err := d.readPending(b[:]); err != nil {
			return 0, err
		}
		return int32(binary.LittleEndian.Uint32(b[:])), nil
	}
	return 0, &ParseError{Field: "int",
		Msg: fmt.Sprintf("unexpected type byte %d for integer", t)}
}

// readPending fills b, draining the pushback buffer first.
func (d *Decoder) readPending(b []byte) error {
	for len(b) > 0 && len(d.pending) > 0 {
		b[0] = d.pending[0]
		d.pending = d.pending[1:]
		b = b[1:]
	}
	if len(b) == 0 {
		return nil
	}
	_, err := io.ReadFull(d.br, b)
	return err
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
		if err := d.readLe(&l); err != nil {
			return "", err
		}
		if l == bcfInt32Missing {
			return "", nil
		}
		buf := make([]byte, l)
		if err := d.readPending(buf); err != nil {
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
		if err := d.readLe(&idx); err != nil {
			return nil, err
		}
		if idx == bcfFilterMissing {
			return nil, nil
		}
		// LookupFilter encodes the implicit PASS at 0,
		// so we don't need a special case here.
		name := d.hdr.LookupFilter(int(idx))
		if name == "" {
			return []string{fmt.Sprintf("FILTER#%d", idx)}, nil
		}
		return []string{name}, nil
	case typeInt32:
		// Array form: int32 length + that many int32 indices.
		var n uint32
		if err := d.readLe(&n); err != nil {
			return nil, err
		}
		if n == 0 {
			return nil, nil
		}
		out := make([]string, 0, n)
		for i := uint32(0); i < n; i++ {
			var idx uint32
			if err := d.readLe(&idx); err != nil {
				return nil, err
			}
			if idx == bcfFilterMissing {
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

// readTypedInfo reads the n_info typed INFO entries and
// returns the resolved Info map.
func (d *Decoder) readTypedInfo(n uint32) (map[string]InfoValue, error) {
	out := make(map[string]InfoValue)
	for i := uint32(0); i < n; i++ {
		// Key: 2 bytes uint16 = index into the header's
		// ordered INFO declarations. We need an ordered
		// list — Header.Infos is a map; we look up by
		// iterating to find the i-th key.
		// Key: a typed atomic integer giving the offset of
		// the INFO field key into the dictionary (spec
		// section 5).
		keyVal, err := d.readTypedInt()
		if err != nil {
			return nil, err
		}
		keyIdx := int(keyVal)
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
	//nolint:dupl // INFO and FORMAT readers parallel by design
	v := InfoValue{Type: def.Type}

	// Flag has no vector body; one type byte.
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

	// One type descriptor byte per spec §BcfTypeEncoding:
	// high nibble = size (1..14 inline, 15 = overflow follows,
	// 0 = missing), low nibble = BCF type code.
	// If size=15, read a typed Integer from the byte stream
	// giving the true vector length.
	t, err := d.readByte()
	if err != nil {
		return v, err
	}
	if t == typeMissing {
		return v, nil
	}
	size := int(t >> 4)
	typeCode := t & 0x0F
	count := size
	if size == 15 {
		// Overflow: typed Integer giving the true count.
		// Read the int's type byte (typically 0x11=int8) and
		// its value.
		ct, err := d.readByte()
		if err != nil {
			return v, err
		}
		switch ct & 0x0F {
		case typeInt8:
			b, err := d.readByte()
			if err != nil {
				return v, err
			}
			count = int(int8(b))
		case typeInt16:
			var x uint16
			if err := d.readLe(&x); err != nil {
				return v, err
			}
			count = int(int16(x))
		case typeInt32:
			var x uint32
			if err := d.readLe(&x); err != nil {
				return v, err
			}
			count = int(int32(x))
		default:
			return v, &ParseError{Field: "INFO",
				Msg: fmt.Sprintf("overflow Integer type byte %d", ct)}
		}
	}
	if count == 0 {
		return v, nil
	}

	switch def.Type {
	case "Integer":
		v.Ints = make([]int64, count)
		switch typeCode {
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
				if err := d.readLe(&x); err != nil {
					return InfoValue{}, err
				}
				v.Ints[i] = int64(int16(x))
			}
		case typeInt32:
			for i := 0; i < count; i++ {
				var x uint32
				if err := d.readLe(&x); err != nil {
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
		if typeCode != typeFloat {
			return InfoValue{}, &ParseError{Field: "INFO",
				Msg: fmt.Sprintf("Float field has type byte %d", t)}
		}
		for i := 0; i < count; i++ {
			var bits uint32
			if err := d.readLe(&bits); err != nil {
				return InfoValue{}, err
			}
			v.Floats[i] = float64(math.Float32frombits(bits))
		}
	case "String":
		v.Strs = make([]string, count)
		if typeCode != typeString {
			return InfoValue{}, &ParseError{Field: "INFO",
				Msg: fmt.Sprintf("String field has type byte %d", t)}
		}
		for i := 0; i < count; i++ {
			var l uint32
			if err := d.readLe(&l); err != nil {
				return InfoValue{}, err
			}
			if l == bcfIntMissingInt32 {
				continue
			}
			buf := make([]byte, l)
			if err := d.readPending(buf); err != nil {
				return InfoValue{}, err
			}
			v.Strs[i] = string(buf)
		}
	case "Character":
		if typeCode != typeString || count != 1 {
			return InfoValue{}, &ParseError{Field: "INFO",
				Msg: fmt.Sprintf("Character field has type byte %d", t)}
		}
		var l uint32
		if err := d.readLe(&l); err != nil {
			return InfoValue{}, err
		}
		buf := make([]byte, l)
		if err := d.readPending(buf); err != nil {
			return InfoValue{}, err
		}
		v.Str = string(buf)
	default:
		return InfoValue{}, &ParseError{Field: "INFO",
			Msg: fmt.Sprintf("unsupported field type %q", def.Type)}
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
	if err := d.readLe(&l); err != nil {
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
	return d.readPending(buf)
}

// skipFmtBlock is no longer used; the v2 FORMAT reader
// consumes the per-sample block above. Kept for backwards
// compatibility with any external callers.
func (d *Decoder) skipFmtBlock() error {
	return &ParseError{Field: "FORMAT",
		Msg: "skipFmtBlock is no longer used; v2 reads the per-sample block directly"}
}

// readByte reads one byte (helper for the Decoder).
func (d *Decoder) readByte() (byte, error) {
	if len(d.pending) > 0 {
		b := d.pending[0]
		d.pending = d.pending[1:]
		return b, nil
	}
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
	if def, ok := h.Formats[id]; ok && def != nil {
		return def, true
	}
	// Synthesise a default for FORMAT keys that exist
	// in declaration order but lack a typed ##FORMAT
	// line (rare, but allowed). Integer is a safe
	// default for the common case (GT-as-int, DP, GQ).
	return &FormatDef{ID: id, Number: "1", Type: "Integer"}, true
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
