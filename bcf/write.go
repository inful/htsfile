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

// typeByteAtomic returns the BCF2 type-descriptor byte for
// an atomic (size-1) value: the size in the high nibble
// (must be 1 for atomic values) and the type in the low
// nibble. Per §BcfTypeEncoding (line 2271):
//
//	"5,6,7,8 bits & The number of elements of the
//	upcoming type. For atomic values, the size must be
//	1."
//
// yEncodes as `(1<<4) | type`, e.g. `0x13` for a single
// int32, `0x15` for a single float32. The decoder ignores
// the high nibble (it reads the count from the header's
// declared Number), so this is wire-compatible with the
// existing reader.
func typeByteAtomic(t uint8) uint8 {
	return (1 << 4) | t
}

// bcfIntReservedMin is the lowest of the eight Integer
// values that BCF2 (§BcfTypeEncoding, line 2306 / §1.4.2)
// reserves for its type-prefix encoding. The eight values
// [MinInt32, MinInt32+7] are forbidden in both VCF and
// BCF; the writer refuses to emit them so a BCF file
// round-trips through bcftools / htslib cleanly.
const (
	bcfIntReservedMin = int64(math.MinInt32)
	bcfIntReservedMax = bcfIntReservedMin + 7 // -2147483641
)

// bcfIntegerReserved reports whether v is one of the
// eight BCF2-reserved Integer values. See [T4.22].
func bcfIntegerReserved(v int64) bool {
	return v >= bcfIntReservedMin && v <= bcfIntReservedMax
}

// checkBcfIntegerRange returns an error if any element of
// ints is in the BCF2-reserved Integer range
// [MinInt32, MinInt32+7]. Used by the writer to refuse
// forbidden values before they hit the wire. Format
// errors are kept under one message so the writer
// returns one error per field rather than one per
// value.
func checkBcfIntegerRange(field string, ints []int64) error {
	for _, v := range ints {
		if bcfIntegerReserved(v) {
			return fmt.Errorf("bcf: %s value %d is in the BCF2 forbidden range [%d, %d] (§1.4.2)",
				field, v, bcfIntReservedMin, bcfIntReservedMax)
		}
	}
	return nil
}

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
// v2 scope: per-sample FORMAT data is encoded when the
// record's Format field is non-empty. n_fmt is set from
// len(r.Format); one uint16 index per FORMAT key, then one
// typed-value block per sample, in the same order as the
// header's sample list.
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w *bgzf.Writer

	hdr        *Header
	hdrWritten bool
	closed     bool

	// version is the BCF minor version to emit (1 for 2.1,
	// 2 for 2.2). Defaults to 1.
	version uint8

	// Cached contig index for fast CHROM → int32 mapping.
	contigIdx map[string]int32
}

// NewEncoder wraps w in an [Encoder]. The default BCF
// version is 2.1; call [Encoder.SetBCFVersion] before
// [Encoder.WriteHeader] to emit 2.2.
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{w: bgzf.NewWriter(w), version: 1}
}

// SetBCFVersion selects the BCF minor version to emit.
// Pass 1 for BCF 2.1, 2 for BCF 2.2. Must be called before
// [Encoder.WriteHeader]; calling it after the header is
// written returns an error.
func (e *Encoder) SetBCFVersion(v uint8) error {
	if e.hdrWritten {
		return errors.New("bcf: SetBCFVersion after WriteHeader")
	}
	if v != 1 && v != 2 {
		return fmt.Errorf("bcf: unsupported BCF version %d (want 1 or 2)", v)
	}
	e.version = v
	return nil
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
	// Select the magic for the configured version.
	magic := Magic
	if e.version == 2 {
		magic = magic22
	}
	payload.Write(magic[:])
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
// n_fmt (uint16) + ID + REF + ALT + FILTER + QUAL + INFO +
// (per-sample FORMAT block when n_fmt > 0).
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
	shared, err := encodeRecord(r, chromIdx, e.hdr, e.version)
	if err != nil {
		return fmt.Errorf("bcf: encoding record: %w", err)
	}
	individ, err := encodeIndivid(r, e.hdr, e.version)
	if err != nil {
		return fmt.Errorf("bcf: encoding per-sample data: %w", err)
	}
	// Spec section 5, BcfSiteEncoding: every record is
	// preceded by l_shared (CHROM..INFO) and l_indiv
	// (FORMAT + per-sample genotype fields).
	var lens [8]byte
	binary.LittleEndian.PutUint32(lens[0:4], uint32(len(shared)))
	binary.LittleEndian.PutUint32(lens[4:8], uint32(len(individ)))
	var rec bytes.Buffer
	rec.Write(lens[:])
	rec.Write(shared)
	rec.Write(individ)
	if _, err := e.w.Write(rec.Bytes()); err != nil {
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
// ID + REF + ALT + FILTER + QUAL + INFO + per-sample FORMAT.
// encodeRecord serializes one Record into the BCF2 "site"
// representation defined in VCFv4.5 section 5, "BCF2 site
// information encoding":
//
//	l_shared  uint32   (written by the caller as a length prefix)
//	l_indiv   uint32   (ditto)
//	CHROM     int32
//	POS       int32    0-based
//	rlen      int32
//	QUAL      float32  0x7F800001 when missing
//	n_info    uint16
//	n_allele  uint16   REF + ALTs
//	n_sample  uint24
//	n_fmt     uint8
//	ID        typed string
//	REF+ALT   n_allele typed strings
//	FILTER    typed vector of int
//	INFO      n_info x (typed atomic int key, typed vector value)
//
// The per-sample FORMAT block is appended by the caller path
// in writeRecord; this function returns the shared part only.
func encodeRecord(r *Record, chromIdx int32, hdr *Header, version uint8) ([]byte, error) {
	var payload bytes.Buffer
	nInfo := uint16(len(r.Info))
	nAllele := uint16(1 + len(r.Alt)) // REF plus the ALTs
	nFmt := uint8(len(r.Format))

	// Fixed 24-byte site prefix.
	var prefix [24]byte
	binary.LittleEndian.PutUint32(prefix[0:4], uint32(chromIdx))
	pos0 := uint32(int32(r.Pos) - 1)
	if r.Pos <= 0 {
		pos0 = 0xFFFFFFFF // -1 in int32 (BCF's missing sentinel)
	}
	binary.LittleEndian.PutUint32(prefix[4:8], pos0)
	// rlen: length of the record projected onto the reference.
	// Spec (BcfSiteEncoding): the maximum of len(REF) and the
	// lengths inferred from the SVLEN/LEN of any symbolic
	// alleles. We have no reference or symbolic-length metadata
	// on bcf.Record, so len(REF) is the best available answer
	// and is what htslib writes for simple records.
	rlen := int32(len(r.Ref))
	if rlen < 1 {
		rlen = 1
	}
	binary.LittleEndian.PutUint32(prefix[8:12], uint32(rlen))
	// QUAL: a plain float32 in the prefix. The per-type
	// MISSING sentinel 0x7F800001 marks "unknown".
	qualBits := bcfFloatMissing
	if r.QualKnown {
		qualBits = math.Float32bits(float32(r.Qual))
	}
	binary.LittleEndian.PutUint32(prefix[12:16], qualBits)
	binary.LittleEndian.PutUint16(prefix[16:18], nInfo)
	binary.LittleEndian.PutUint16(prefix[18:20], nAllele)
	// n_sample occupies 24 bits at offset 20, n_fmt the final
	// byte.
	if hdr != nil {
		binary.LittleEndian.PutUint32(prefix[20:24],
			uint32(len(hdr.SampleNames))<<8|uint32(nFmt))
	} else {
		binary.LittleEndian.PutUint32(prefix[20:24], uint32(nFmt))
	}
	payload.Write(prefix[:])

	// ID (typed string; 0x07 marks missing).
	writeTypedString(&payload, r.ID, version)
	// REF then ALT, all typed strings.
	writeTypedString(&payload, r.Ref, version)
	for _, alt := range r.Alt {
		writeTypedString(&payload, alt, version)
	}
	// FILTER: typed vector of dictionary offsets.
	writeTypedFilter(&payload, r.Filter, hdr)
	// INFO: n_info pairs.
	if err := writeTypedInfo(&payload, r.Info, hdr); err != nil {
		return nil, err
	}
	return payload.Bytes(), nil
}

// encodeIndivid writes the per-sample (l_indiv) part: n_fmt
// typed atomic integer FORMAT-key offsets followed by one
// typed vector per key (spec section 5, GenotypeEncoding).
func encodeIndivid(r *Record, hdr *Header, version uint8) ([]byte, error) {
	var buf bytes.Buffer
	if len(r.Format) == 0 {
		return nil, nil
	}
	for _, keyIdx := range r.Format {
		writeTypedInt(&buf, int32(keyIdx))
	}
	if err := writeFormatBlock(&buf, r, hdr, version); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// writeTypedString emits a type byte + optional length prefix
// + string body. Missing (empty) is typeMissing. For BCF 2.2
// (version == 2) strings of length 15 or more use the typeLong
// form (type byte 0xF + uint32 length) instead of the standard
// typeString (7) form. The reader accepts both.
func writeTypedString(w *bytes.Buffer, s string, version uint8) {
	if s == "" {
		w.WriteByte(typeMissing)
		return
	}
	if version == 2 && len(s) >= 15 {
		w.WriteByte(typeLong)
	} else {
		w.WriteByte(typeString)
	}
	var lb [4]byte
	binary.LittleEndian.PutUint32(lb[:], uint32(len(s)))
	w.Write(lb[:])
	w.WriteString(s)
}

// writeTypedInt emits a typed atomic 32-bit integer: the type
// byte (size 1 in the high nibble, int32 in the low nibble)
// followed by the little-endian value. Spec section 5 uses
// this form for the INFO dictionary key and for the FORMAT
// key list.
func writeTypedInt(w *bytes.Buffer, v int32) {
	w.WriteByte(typeByteAtomic(typeInt32))
	var b [4]byte
	binary.LittleEndian.PutUint32(b[:], uint32(v))
	w.Write(b[:])
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
			// Declared filters sit at +1 (index 0 is the
			// implicit PASS).
			return uint32(i + 1)
		}
	}
	if name == "PASS" {
		return 0
	}
	return bcfFilterMissing
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
		def, _ := hdr.InfoAt(idx)
		v := val
		// Spec section 5: "The first value must be a typed
		// atomic integer giving the offset of the INFO field
		// key into the dictionary."
		writeTypedInt(w, int32(idx))
		if err := writeTypedInfoValue(w, &v, def); err != nil {
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

// writeTypedInfoValue serializes a single INFO value per
// the BCF 2 type-byte encoding (spec §BcfTypeEncoding,
// lines 2253-2309). The four high bits of the type byte
// hold the vector size: 1..14 are inline, 15 means the true
// count follows as a typed Integer in the byte stream, 0
// means MISSING. The four low bits hold the BCF type code
// (1=int8, 2=int16, 3=int32, 5=float32, 7=char/string).
//
// Flag has no vector body: typeMissing marks absence, any
// other byte marks presence.
func writeTypedInfoValue(w *bytes.Buffer, val *InfoValue, def *InfoDef) error {
	switch val.Type {
	case "Flag":
		if val.Flag {
			w.WriteByte(typeInt8)
		} else {
			w.WriteByte(typeMissing)
		}
	case "Integer":
		count, isVar := len(val.Ints), false
		if def != nil {
			_, isVar = ParseNumber(def.Number)
		}
		size, overflow := pickBcfSize(count, isVar)
		w.WriteByte(typeByteSize(uint8(size), typeInt32))
		if overflow {
			// Typed Integer overflow: type byte + value.
			w.WriteByte(typeByteAtomic(typeInt8))
			w.WriteByte(byte(count))
		}
		if err := checkBcfIntegerRange("INFO", val.Ints); err != nil {
			return err
		}
		for _, n := range val.Ints {
			var b [4]byte
			binary.LittleEndian.PutUint32(b[:], uint32(int32(n)))
			w.Write(b[:])
		}
	case "Float":
		count, isVar := len(val.Floats), false
		if def != nil {
			_, isVar = ParseNumber(def.Number)
		}
		size, overflow := pickBcfSize(count, isVar)
		w.WriteByte(typeByteSize(uint8(size), typeFloat))
		if overflow {
			w.WriteByte(typeByteAtomic(typeInt8))
			w.WriteByte(byte(count))
		}
		for _, f := range val.Floats {
			var b [4]byte
			binary.LittleEndian.PutUint32(b[:], math.Float32bits(float32(f)))
			w.Write(b[:])
		}
	case "String":
		count, isVar := len(val.Strs), false
		if def != nil {
			_, isVar = ParseNumber(def.Number)
		}
		size, overflow := pickBcfSize(count, isVar)
		w.WriteByte(typeByteSize(uint8(size), typeString))
		if overflow {
			w.WriteByte(typeByteAtomic(typeInt8))
			w.WriteByte(byte(count))
		}
		for _, s := range val.Strs {
			var lb [4]byte
			binary.LittleEndian.PutUint32(lb[:], uint32(len(s)))
			w.Write(lb[:])
			w.WriteString(s)
		}
	case "Character":
		w.WriteByte(typeByteAtomic(typeString))
		var lb [4]byte
		binary.LittleEndian.PutUint32(lb[:], uint32(len(val.Str)))
		w.Write(lb[:])
		w.WriteString(val.Str)
	default:
		return fmt.Errorf("bcf: unknown INFO type %q", val.Type)
	}
	return nil
}

// pickBcfSize picks the size nibble for a BCF type byte.
// For variable-length and counts >= 15 we use the overflow
// marker and emit the true count as a typed int8 right after.
func pickBcfSize(count int, isVar bool) (size int, overflow bool) {
	if isVar {
		return 15, true
	}
	if count >= 15 {
		return 15, true
	}
	return count, false
}

// typeByteAtomic is a convenience for a size-1 typed value:
// (1<<4)|lo.
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
