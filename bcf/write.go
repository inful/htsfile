package bcf

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strings"

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

	// version is the BCF minor version to emit. Defaults to 2
	// (BCF 2.2) because that is what the current spec fixes
	// (tex line 2088) and the only version htslib 1.24+
	// will read.
	version uint8

	// Cached contig index for fast CHROM → int32 mapping.
	contigIdx map[string]int32
}

// NewEncoder wraps w in an [Encoder]. The default BCF
// version is 2.1; call [Encoder.SetBCFVersion] before
// [Encoder.WriteHeader] to emit 2.2.
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{w: bgzf.NewWriter(w), version: 2}
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
	// Select the magic for the configured version. BCF 2.1 is
	// emitted only on explicit SetBCFVersion(1); the default is
	// 2.2, which is the version the spec fixes and the only one
	// current htslib reads.
	magic := Magic
	if e.version == 1 {
		magic = magic21
	}
	payload.Write(magic[:])
	// The header text is NUL-terminated and l_text counts the
	// terminator (spec §BCF2 header, tex lines 2086-2091:
	// "l_text - Length of the 'text' field, including the
	// terminating NUL character"). htslib writes it the same
	// way; our previous output omitted the NUL and set l_text
	// to the bare text length.
	var lt [4]byte
	binary.LittleEndian.PutUint32(lt[:], uint32(len(text)+1))
	payload.Write(lt[:])
	payload.WriteString(text)
	payload.WriteByte(0)
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
	// n_sample is a uint24 occupying bytes [20:22] and n_fmt a
	// uint8 at byte [23] (spec §BcfSiteEncoding, tex line 2192).
	// So the packed uint32 is n_sample | n_fmt<<24 - NOT
	// n_sample<<8 | n_fmt. The old expression put n_fmt in the
	// low byte and n_sample in bytes 21-23, which a conformant
	// reader decodes as n_sample=772, n_fmt=0.
	var nSample uint32
	if hdr != nil {
		nSample = uint32(len(hdr.SampleNames))
	}
	binary.LittleEndian.PutUint32(prefix[20:24], nSample|uint32(nFmt)<<24)
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
	// Spec §GenotypeEncoding (tex line 2233): "Suppose there
	// are i genotype fields in a specific record. Each i is
	// encoded by a triplet" - fmt_key, fmt_type, fmt_values.
	//
	// An earlier layout wrote every fmt_key first and then
	// every fmt_type/fmt_values pair. htslib rejected that
	// with "individuals section malformed or too short".
	for pos, keyIdx := range r.Format {
		name := ""
		if hdr != nil && int(keyIdx) < len(hdr.FormatsOrder) {
			name = hdr.FormatsOrder[keyIdx]
		}
		writeTypedInt(&buf, int32(hdr.IDXFor(name)))
		if err := writeFormatField(&buf, r, hdr, version, keyIdx, pos); err != nil {
			return nil, err
		}
	}
	return buf.Bytes(), nil
}

// writeTypedString emits a type byte + optional length prefix
// + string body. Missing (empty) is typeMissing. For BCF 2.2
// (version == 2) strings of length 15 or more use the typeLong
// form (type byte 0xF + uint32 length) instead of the standard
// typeString (7) form. The reader accepts both.
func writeTypedString(w *bytes.Buffer, s string, _ uint8) {
	if s == "" {
		// Spec §BcfSiteEncoding line 2196: "ID - typed string -
		// Variant identifier; 0x07 for a missing value". 0x07 is
		// also the descriptor for a zero-length char vector, so
		// htslib reads it as an empty ID. Writing 0x00 here
		// made htslib reject the record with "Invalid ID type 0
		// (null)".
		w.WriteByte(typeString)
		return
	}
	writeBCFString(w, s)
}

// writeTypedInt emits a typed atomic integer: the type byte
// (size 1 in the high nibble, the integer type in the low
// nibble) followed by the little-endian value. Spec section 5
// uses this form for the INFO dictionary key and for the FORMAT
// key list.
//
// The width is the narrowest that fits, matching htslib: a
// dictionary index is almost always small enough for int8.
func writeTypedInt(w *bytes.Buffer, v int32) {
	width := typeInt32
	if v >= -128 && v <= 127 {
		width = typeInt8
	}
	w.WriteByte(typeByteAtomic(width))
	writeIntValue(w, v, width)
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
	// FILTER is a typed *vector* of dictionary offsets, so the
	// descriptor carries the element count in its high nibble
	// (spec §BcfTypeEncoding). Writing a bare type code such as
	// 0x01 declares size 0 (MISSING) and desynchronises the
	// reader.
	vals := make([]int32, 0, len(filters))
	for _, f := range filters {
		vals = append(vals, int32(hdr.IDXFor(f)))
	}
	width := intWidthFor32(vals)
	size, overflow := pickBcfSize(len(vals))
	w.WriteByte(typeByteSize(uint8(size), width))
	if overflow {
		w.WriteByte(typeByteAtomic(typeInt8))
		w.WriteByte(byte(len(vals)))
	}
	writeTypedIntVector(w, vals)
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
		v := val
		// Spec section 5: "The first value must be a typed
		// atomic integer giving the offset of the INFO field
		// key into the dictionary."
		writeTypedInt(w, int32(hdr.IDXFor(key)))
		if err := writeTypedInfoValue(w, &v); err != nil {
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
func writeTypedInfoValue(w *bytes.Buffer, val *InfoValue) error {
	switch val.Type {
	case "Flag":
		// htslib writes a single 0x00 byte for a present
		// Flag; the INFO entry's presence is what carries
		// the information. Verified against `bcftools
		// view -Ob`.
		if val.Flag {
			w.WriteByte(0)
		} else {
			w.WriteByte(typeMissing)
		}
	case "Integer":
		ints, missing := val.Ints, val.IntMissing
		if len(ints) == 0 {
			// An INFO key with an empty value list is a single
			// missing element, which is what htslib writes
			// (`bcftools view -Ob` on "VR=" gives "11 06 11 80").
			ints = []int64{0}
			missing = []bool{true}
		}
		count := len(ints)
		size, overflow := pickBcfSize(count)
		if err := checkBcfIntegerRange("INFO", ints); err != nil {
			return err
		}
		width := intWidthFor(ints)
		w.WriteByte(typeByteSize(uint8(size), width))
		if overflow {
			// Typed Integer overflow: type byte + value.
			w.WriteByte(typeByteAtomic(typeInt8))
			w.WriteByte(byte(count))
		}
		for i, n := range ints {
			if i < len(missing) && missing[i] {
				writeIntValue(w, int32(intMissingFor(width)), width)
				continue
			}
			writeIntValue(w, int32(n), width)
		}
	case "Float":
		floats, missing := val.Floats, val.FloatMissing
		if len(floats) == 0 {
			floats = []float64{0}
			missing = []bool{true}
		}
		count := len(floats)
		size, overflow := pickBcfSize(count)
		w.WriteByte(typeByteSize(uint8(size), typeFloat))
		if overflow {
			w.WriteByte(typeByteAtomic(typeInt8))
			w.WriteByte(byte(count))
		}
		for i, f := range floats {
			var b [4]byte
			if i < len(missing) && missing[i] {
				// Per-type MISSING: the float32 signaling NaN
				// 0x7F800001 (T4.24 / line 2333).
				binary.LittleEndian.PutUint32(b[:], bcfFloatMissing)
			} else {
				binary.LittleEndian.PutUint32(b[:], math.Float32bits(float32(f)))
			}
			w.Write(b[:])
		}
	case "String":
		// A String-typed INFO field is a char vector (spec
		// §BcfTypeEncoding, tex line 2276). The list is written
		// as ONE char run with the elements joined by commas,
		// and the descriptor nibble is the run's total character
		// count. Verified byte-for-byte against `bcftools view
		// -Ob`: CSQ=alpha,beta,gamma becomes
		// f7 11 10 "alpha,beta,gamma".
		run := strings.Join(val.Strs, ",")
		if run == "" {
			// An absent or empty String is a zero-length char
			// vector (0x07), the same form as a missing ID
			// (spec line 2196) - not typeMissing.
			w.WriteByte(typeString)
			return nil
		}
		size, overflow := pickBcfSize(len(run))
		w.WriteByte(typeByteSize(uint8(size), typeString))
		if overflow {
			w.WriteByte(typeByteAtomic(typeInt8))
			w.WriteByte(byte(len(run)))
		}
		w.WriteString(run)
	case "Character":
		writeBCFString(w, val.Str)
	default:
		return fmt.Errorf("bcf: unknown INFO type %q", val.Type)
	}
	return nil
}

// intWidthFor picks the narrowest BCF integer type that can hold
// every value in vals: int8, then int16, then int32.
//
// The spec (§BcfTypeEncoding, tex line 2304) allows any of
// "8, 16, or 32 bit values" and says "It is up to the encoder
// to determine the appropriate ranged value to use". htslib
// makes the same choice, so matching it keeps our output
// byte-identical to the reference implementation.
func intWidthFor(vals []int64) uint8 {
	fits8 := true
	fits16 := true
	for _, v := range vals {
		if v < -128 || v > 127 {
			fits8 = false
		}
		if v < -32768 || v > 32767 {
			fits16 = false
		}
	}
	switch {
	case fits8:
		return typeInt8
	case fits16:
		return typeInt16
	default:
		return typeInt32
	}
}

// intWidthFor32 is intWidthFor for a []int32.
func intWidthFor32(vals []int32) uint8 {
	fits8 := true
	fits16 := true
	for _, v := range vals {
		if v < -128 || v > 127 {
			fits8 = false
		}
		if v < -32768 || v > 32767 {
			fits16 = false
		}
	}
	switch {
	case fits8:
		return typeInt8
	case fits16:
		return typeInt16
	default:
		return typeInt32
	}
}

// writeIntValue writes one integer in the given width.
func writeIntValue(w *bytes.Buffer, v int32, width uint8) {
	switch width {
	case typeInt8:
		w.WriteByte(byte(int8(v)))
	case typeInt16:
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(int16(v)))
		w.Write(b[:])
	default:
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(v))
		w.Write(b[:])
	}
}

// writeTypedIntVector writes a vector of integers using the
// narrowest width that fits all of them.
func writeTypedIntVector(w *bytes.Buffer, vals []int32) {
	width := intWidthFor32(vals)
	for _, v := range vals {
		writeIntValue(w, v, width)
	}
}

// pickBcfSize picks the size nibble for a BCF type byte.
// For variable-length and counts >= 15 we use the overflow
// marker and emit the true count as a typed int8 right after.
func pickBcfSize(count int) (size int, overflow bool) {
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
