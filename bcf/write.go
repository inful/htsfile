package bcf

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"

	"github.com/inful/htsfile/internal/bgzf"
)

// Encoder writes a BCF v2.1 document to an underlying
// [io.Writer]. The output is BGZF-compressed per the BCF
// spec.
//
// Typical usage:
//
//	f, _ := os.Create("out.bcf")
//	defer f.Close()
//	enc := bcf.NewEncoder(f)
//	if err := enc.WriteHeader(hdr); err != nil { ... }
//	for _, r := range recs { enc.WriteRecord(r) }
//	enc.Close()
//
// v1 scope: records are written without per-sample FORMAT
// data (matching what the v1 decoder accepts). n_fmt is
// always emitted as 0; the decoder skips the FORMAT block
// when n_fmt == 0, which means the produced stream round-
// trips cleanly through the existing decoder.
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w *bgzf.Writer

	hdr        *Header
	hdrWritten bool
	closed     bool

	// Cached contig index for fast CHROM → int32 mapping.
	contigIdx map[string]int32
}

// NewEncoder wraps w in an [Encoder].
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{w: bgzf.NewWriter(w)}
}

// WriteHeader emits the BCF file header: magic + l_text +
// header text. The header is round-tripped through ParseHeader
// on read; emitting the original text bytes (rather than
// re-serialising from the parsed structure) preserves any
// unknown pragmas verbatim.
//
// A nil header returns an error. WriteHeader must be called
// before any WriteRecord.
func (e *Encoder) WriteHeader(hdr *Header) error {
	if hdr == nil {
		return errors.New("bcf: WriteHeader: nil header")
	}
	if e.hdrWritten {
		return errors.New("bcf: WriteHeader called twice")
	}
	e.hdrWritten = true

	text := hdr.Raw
	if text == "" {
		// Fall back to re-emitting the parsed header text.
		// This path is hit only if a caller constructs a
		// Header by hand (no header text round-tripped from
		// a parsed file).
		text = serialiseHeader(hdr)
	}
	var payload bytes.Buffer
	payload.Write(Magic[:])
	var lt [4]byte
	binary.LittleEndian.PutUint32(lt[:], uint32(len(text)))
	payload.Write(lt[:])
	payload.WriteString(text)
	if _, err := e.w.Write(payload.Bytes()); err != nil {
		return fmt.Errorf("bcf: writing header: %w", err)
	}

	// Build contig lookup. (The binary reference dictionary
	// — n_ref + per-ref entries — is part of the BCF v2.1
	// spec but the v1 decoder only parses the VCF text. We
	// keep the contig map in the encoder for the record-time
	// chromIdx resolution but don't emit the dictionary
	// bytes; this matches what the v1 test fixtures do and
	// what the v1 decoder expects.)
	e.contigIdx = make(map[string]int32, len(hdr.Contigs))
	for i, c := range hdr.Contigs {
		e.contigIdx[c] = int32(i)
	}
	e.hdr = hdr
	return nil
}

// WriteRecord emits one BCF record: CHROM (int32) + POS
// (int32) + rlen (int32, unused) + n_info (uint32) +
// n_fmt (uint16, always 0 in v1) + ID + REF + ALT + FILTER
// + QUAL + INFO. The per-sample FORMAT block is skipped
// (n_fmt = 0).
//
// A nil record is a no-op. CHROM must appear in the header's
// Contig dictionary (LookupContig returns it); an unknown
// CHROM produces an error so callers get a clear failure
// rather than silent data corruption.
func (e *Encoder) WriteRecord(r *Record) error {
	if r == nil {
		return nil
	}
	if !e.hdrWritten {
		return errors.New("bcf: WriteRecord before WriteHeader")
	}
	if e.closed {
		return errors.New("bcf: WriteRecord after Close")
	}
	chromIdx, ok := e.contigIdx[r.CHROM]
	if !ok {
		return fmt.Errorf("bcf: record CHROM %q not in header", r.CHROM)
	}
	payload, err := encodeRecord(r, chromIdx, e.hdr)
	if err != nil {
		return fmt.Errorf("bcf: encoding record: %w", err)
	}
	if _, err := e.w.Write(payload); err != nil {
		return fmt.Errorf("bcf: writing record: %w", err)
	}
	return nil
}

// Close flushes any buffered BGZF data and emits the
// canonical 28-byte EOF marker.
func (e *Encoder) Close() error {
	if e.closed {
		return nil
	}
	e.closed = true
	return e.w.Close()
}

// encodeRecord serializes one Record into the BGZF block
// payload: prefix (CHROM + POS + rlen + n_info + n_fmt) +
// ID + REF + ALT + FILTER + QUAL + INFO.
func encodeRecord(r *Record, chromIdx int32, hdr *Header) ([]byte, error) {
	var payload bytes.Buffer
	// Prefix: CHROM (int32) + POS (int32) + rlen (int32) +
	// n_info (uint32) + n_fmt (uint16) = 18 bytes.
	prefix := make([]byte, 18)
	binary.LittleEndian.PutUint32(prefix[0:4], uint32(chromIdx))
	pos0 := uint32(int32(r.Pos) - 1)
	if r.Pos <= 0 {
		pos0 = 0xFFFFFFFF // -1 in int32 (BCF's missing sentinel)
	}
	binary.LittleEndian.PutUint32(prefix[4:8], pos0)
	// rlen (uint32) — unused in v1; emit 0.
	binary.LittleEndian.PutUint32(prefix[8:12], 0)
	nInfo := uint32(len(r.Info))
	binary.LittleEndian.PutUint32(prefix[12:16], nInfo)
	// n_fmt (uint16 LE) — always 0 in v1 (skip FORMAT block).
	payload.Write(prefix[:18])
	// ID (typed string).
	writeTypedString(&payload, r.ID)
	// REF (typed string).
	writeTypedString(&payload, r.Ref)
	// ALT (zero or more typed strings, terminated by FILTER).
	for _, alt := range r.Alt {
		writeTypedString(&payload, alt)
	}
	// FILTER (typed: missing | single int32 | array of int32).
	writeTypedFilter(&payload, r.Filter, hdr)
	// QUAL (typed: missing | float64).
	writeTypedQual(&payload, r.Qual, r.QualKnown)
	// INFO (n_info typed entries; each is keyed by uint16
	// header index plus a typed value).
	if err := writeTypedInfo(&payload, r.Info, hdr); err != nil {
		return nil, err
	}
	return payload.Bytes(), nil
}

// writeTypedString emits a type byte + optional length prefix
// + string body. Missing (empty) is typeMissing.
func writeTypedString(w *bytes.Buffer, s string) {
	if s == "" {
		w.WriteByte(typeMissing)
		return
	}
	w.WriteByte(typeString)
	var lb [4]byte
	binary.LittleEndian.PutUint32(lb[:], uint32(len(s)))
	w.Write(lb[:])
	w.WriteString(s)
}

// writeTypedFilter emits the FILTER typed array. BCF FILTER
// stores int32 indices. An empty filter list means typeMissing.
// A single filter is written as typeInt8 + single int32
// (the v1 decoder's preferred single-filter form); multiple
// filters are written as typeInt32 + count + elements.
//
// Filter names not declared in hdr.Filters are stored with
// the bcfFilterMissing sentinel (0xFFFFFFFF). The decoder's
// LookupFilter will return empty for that sentinel — a clear
// signal that the filter name didn't round-trip, not silent
// data corruption.
func writeTypedFilter(w *bytes.Buffer, filters []string, hdr *Header) {
	if len(filters) == 0 {
		w.WriteByte(typeMissing)
		return
	}
	if len(filters) == 1 {
		w.WriteByte(typeInt8)
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], filterIndex(filters[0], hdr))
		w.Write(b[:])
		return
	}
	w.WriteByte(typeInt32)
	var count [4]byte
	binary.LittleEndian.PutUint32(count[:], uint32(len(filters)))
	w.Write(count[:])
	for _, f := range filters {
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], filterIndex(f, hdr))
		w.Write(b[:])
	}
}

// filterIndex returns the int32 index of filter name in the
// header's ##FILTER declaration list. If no filters are
// declared, "PASS" maps to index 0 (the BCF v2.1 convention;
// the decoder treats 0 as PASS when len(hdr.Filters) == 0).
// Undeclared filters get the bcfFilterMissing sentinel so
// the decoder surfaces a clear "unknown filter" instead of
// silently re-mapping.
func filterIndex(name string, hdr *Header) uint32 {
	if hdr == nil {
		return bcfFilterMissing
	}
	for i, f := range hdr.Filters {
		if f == name {
			return uint32(i)
		}
	}
	if len(hdr.Filters) == 0 && name == "PASS" {
		return 0
	}
	return bcfFilterMissing
}

// writeTypedQual emits the QUAL typed float64. Missing is
// typeMissing (qualKnown=false → emit sentinel).
func writeTypedQual(w *bytes.Buffer, qual float64, qualKnown bool) {
	if !qualKnown {
		w.WriteByte(typeMissing)
		return
	}
	w.WriteByte(typeFloat)
	var b [8]byte
	binary.LittleEndian.PutUint64(b[:], math.Float64bits(qual))
	w.Write(b[:])
}

// writeTypedInfo emits one typed entry per key in info. Each
// entry is the header's uint16 index of the key followed by
// the typed value. The encoder emits keys in the order they
// appear in hdr.InfosOrder (the header's ##INFO declaration
// order) so the decoder's lookup-by-index works cleanly.
//
// INFO keys present in the record but missing from the
// header produce an error — the decoder would reject them
// anyway (header doesn't match key index 0xFFFF).
func writeTypedInfo(w *bytes.Buffer, info map[string]InfoValue, hdr *Header) error {
	if len(info) == 0 || hdr == nil {
		return nil
	}
	for _, key := range hdr.InfosOrder {
		val, ok := info[key]
		if !ok {
			continue
		}
		idx, ok := hdr.LookupInfoIndex(key)
		if !ok {
			continue // unreachable, but defensive
		}
		v := val
		var ib [2]byte
		binary.LittleEndian.PutUint16(ib[:], uint16(idx))
		w.Write(ib[:])
		if err := writeTypedInfoValue(w, &v, hdr); err != nil {
			return err
		}
	}
	// Any keys in info that aren't in the header's
	// InfosOrder are reported as errors so the caller
	// doesn't silently lose fields.
	for k := range info {
		if _, ok := hdr.Infos[k]; !ok {
			return fmt.Errorf("bcf: INFO key %q not declared in header", k)
		}
	}
	return nil
}

// writeTypedInfoValue serializes a single INFO value per
// the header-declared type. The format matches the v1 decoder:
//
//   - Flag: typeMissing = absent, anything else = present. No
//     length prefix, no body.
//   - Fixed-length (Number=1/2/...): type byte + value bytes.
//     No count uint32 (the decoder reads count values from the
//     header declaration, not the stream).
//   - Variable-length (Number=A/R/G/.): length prefix (type byte
//   - count uint32) + value type byte + value bytes.
func writeTypedInfoValue(w *bytes.Buffer, val *InfoValue, hdr *Header) error {
	def, ok := hdr.Infos[val.Type]
	_ = def // type info is in val.Type; def is for Number
	isVar := false
	if ok {
		_, isVar = ParseNumber(def.Number)
	}
	switch val.Type {
	case "Flag":
		// Flag: typeMissing = absent, anything else = present.
		if val.Flag {
			w.WriteByte(typeInt8)
		} else {
			w.WriteByte(typeMissing)
		}
	case "Integer":
		if isVar {
			// Variable-length: length prefix + value type +
			// count + values.
			w.WriteByte(typeInt32) // length-prefix type byte
			var lb [4]byte
			binary.LittleEndian.PutUint32(lb[:], uint32(len(val.Ints)))
			w.Write(lb[:])
		}
		w.WriteByte(typeInt32) // value type byte
		for _, n := range val.Ints {
			var b [4]byte
			binary.LittleEndian.PutUint32(b[:], uint32(int32(n)))
			w.Write(b[:])
		}
	case "Float":
		if isVar {
			w.WriteByte(typeInt32)
			var lb [4]byte
			binary.LittleEndian.PutUint32(lb[:], uint32(len(val.Floats)))
			w.Write(lb[:])
		}
		w.WriteByte(typeFloat)
		for _, f := range val.Floats {
			var b [4]byte
			binary.LittleEndian.PutUint32(b[:], math.Float32bits(float32(f)))
			w.Write(b[:])
		}
	case "String":
		if isVar {
			w.WriteByte(typeInt32)
			var lb [4]byte
			binary.LittleEndian.PutUint32(lb[:], uint32(len(val.Strs)))
			w.Write(lb[:])
		}
		w.WriteByte(typeString)
		for _, s := range val.Strs {
			var lb [4]byte
			binary.LittleEndian.PutUint32(lb[:], uint32(len(s)))
			w.Write(lb[:])
			w.WriteString(s)
		}
	case "Character":
		// Character is always fixed-length 1.
		w.WriteByte(typeString)
		var lb [4]byte
		binary.LittleEndian.PutUint32(lb[:], uint32(len(val.Str)))
		w.Write(lb[:])
		w.WriteString(val.Str)
	default:
		return fmt.Errorf("bcf: unknown INFO type %q", val.Type)
	}
	return nil
}

// serialiseHeader is a no-op fallback for the rare case of
// a Header constructed by hand without Raw text. v1 callers
// always go through ParseHeader or build a Header with Raw
// set, so this path is unreachable in practice. Kept as a
// forward-compatible stub in case future callers construct
// Header from structured fields.
func serialiseHeader(hdr *Header) string {
	if hdr.Raw != "" {
		return hdr.Raw
	}
	// Without Raw text the best we can do is emit a stub
	// fileformat pragma. Callers relying on this should
	// ensure hdr.Raw is set.
	return "##fileformat=VCFv4.2\n"
}
