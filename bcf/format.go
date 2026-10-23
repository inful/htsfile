package bcf

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"math"
	"strings"
)

// Sample is one row of the per-sample FORMAT table. Fields
// maps the format-key header index (the same int16 that
// appears in Record.Format) to its typed value. The
// indirection via header index matches the BCF wire format
// and keeps the in-memory representation compact when many
// records share the same FORMAT layout.
type Sample struct {
	Fields map[uint16]SampleValue
}

// SampleValue is the typed value of one FORMAT field for one
// sample. It mirrors the vcf package's FieldValue shape so
// that downstream code can apply the same per-type accessors.
type SampleValue struct {
	// Type is the header-declared type (Integer / Float /
	// Character / String). Flag is not permitted in FORMAT
	// per the VCF spec; the encoder refuses to emit it.
	Type string

	// Ints holds the typed value(s) for Integer-typed FORMAT
	// keys. GT is also stored here as the spec's int vector.
	Ints []int32

	// Floats holds the typed value(s) for Float-typed FORMAT
	// keys. Stored as float32 because BCF FORMAT floats
	// are 32-bit on the wire.
	Floats []float32

	// Str is the single Character value for Character-typed
	// FORMAT keys. Character fields always have length 1.
	Str string

	// Strs holds the typed value(s) for String-typed FORMAT
	// keys. BCF encodes these as length-prefixed strings.
	Strs []string

	// Missing is true when the field was emitted as the BCF
	// typeMissing sentinel. Per the VCF spec, "." in a
	// sample's column means "not applicable" — distinct
	// from an absent field, which the v1 BCF encoder does
	// not emit.
	Missing bool

	// IntMissing parallels Ints slot-by-slot: IntMissing[i]
	// marks Ints[i] as the per-slot VCF missing value ".".
	// Output round-trip must emit "." for those slots
	// (VCF §1.6.2 line 508).
	IntMissing []bool

	// FloatMissing parallels Floats slot-by-slot. A "." in
	// a Float FORMAT list stays "." on output.
	FloatMissing []bool
}

// readFormatBlockTriplets consumes the per-field FORMAT block
// (spec §GenotypeEncoding, tex line 2233). Each field is a
// triplet: fmt_key (typed int, the header dictionary index),
// fmt_type (a descriptor byte, optionally followed by a typed
// Integer vector length), then fmt_values (one value per
// sample).
func (d *Decoder) readFormatBlockTriplets(nFmt uint8, sampleCount int) ([]uint16, []Sample, error) {
	samples := make([]Sample, sampleCount)
	for i := range samples {
		samples[i].Fields = make(map[uint16]SampleValue, int(nFmt))
	}
	keys := make([]uint16, 0, nFmt)
	for i := uint8(0); i < nFmt; i++ {
		keyIdx, err := d.readFormatKey()
		if err != nil {
			return nil, nil, err
		}
		keys = append(keys, keyIdx)
		def, ok := d.hdr.FormatAt(int(keyIdx))
		if !ok {
			return nil, nil, fmt.Errorf("bcf: FORMAT key index %d not in header", keyIdx)
		}
		vals, err := d.readFieldForAllSamples(def, sampleCount)
		if err != nil {
			return nil, nil, fmt.Errorf("bcf: FORMAT key %d (%s): %w", keyIdx, def.ID, err)
		}
		for sIdx := range vals {
			samples[sIdx].Fields[keyIdx] = vals[sIdx]
		}
	}
	return keys, samples, nil
}

// intWidthMask returns the bit mask for an integer width, used
// to compare raw wire bits against the sentinel constants
// without sign extension getting in the way.
func intWidthMask(width uint8) uint32 {
	switch width {
	case typeInt8:
		return 0xFF
	case typeInt16:
		return 0xFFFF
	default:
		return 0xFFFFFFFF
	}
}

// given integer width. The spec (tex line 2304) defines them
// per width: 0x80, 0x8000, 0x80000000.
func intMissingFor(width uint8) uint32 {
	switch width {
	case typeInt8:
		return 0x80
	case typeInt16:
		return 0x8000
	default:
		return 0x80000000
	}
}

// intEndOfVectorFor returns the per-type END_OF_VECTOR sentinel
// for the given integer width: 0x81, 0x8001, 0x80000001.
func intEndOfVectorFor(width uint8) uint32 {
	switch width {
	case typeInt8:
		return 0x81
	case typeInt16:
		return 0x8001
	default:
		return 0x80000001
	}
}

// readIntValue reads one integer in the width named by the
// descriptor's type nibble (int8, int16 or int32).
func (d *Decoder) readIntValue(width uint8) (int32, error) {
	switch width {
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
	default:
		var b [4]byte
		if err := d.readPending(b[:]); err != nil {
			return 0, err
		}
		return int32(binary.LittleEndian.Uint32(b[:])), nil
	}
}

// readFormatKey reads one fmt_key: a typed atomic integer
// holding the header dictionary index of the FORMAT key.
func (d *Decoder) readFormatKey() (uint16, error) {
	idx, err := d.readTypedInt()
	if err != nil {
		return 0, err
	}
	name, ok := d.hdr.IDAt(int(idx))
	if !ok {
		return 0, fmt.Errorf("bcf: FORMAT dictionary index %d out of range", idx)
	}
	for i, n := range d.hdr.FormatsOrder {
		if n == name {
			return uint16(i), nil
		}
	}
	return 0, fmt.Errorf("bcf: FORMAT key %q not declared in header", name)
}

// readFieldForAllSamples decodes one FORMAT field across N
// samples. The wire shape is the type byte (size in high
// nibble, type in low) followed by N typed values. For
// variable-length vectors, END_OF_VECTOR-padded values are
// dropped so each sample sees only its own cardinality.
//
// Character fields are written as a single byte per sample
// (with a typeMissing escape); the high-nibble size is
// ignored on read.
func (d *Decoder) readFieldForAllSamples(def *FormatDef, nSamples int) ([]SampleValue, error) {
	t, err := d.readByte()
	if err != nil {
		return nil, err
	}
	if t == typeMissing {
		// Whole field is missing — every sample sees
		// Missing=true with the declared type.
		out := make([]SampleValue, nSamples)
		for i := range out {
			out[i] = SampleValue{Missing: true, Type: def.Type}
		}
		return out, nil
	}

	// Per the spec, size goes in the high nibble of the
	// type byte. Size 15 = overflow; the real size follows
	// as a typed INT8.
	size := int((t >> 4) & 0x0F)
	if size == 15 {
		n, err := d.readTypedLen()
		if err != nil {
			return nil, err
		}
		size = n
	}
	typeLo := t & 0x0F

	out := make([]SampleValue, nSamples)

	if def.ID == "GT" {
		// GT is an integer allele vector; rebuild the VCF
		// genotype string from it.
		for s := 0; s < nSamples; s++ {
			alleles := make([]int32, 0, size)
			for i := 0; i < size; i++ {
				v, err := d.readIntValue(typeLo)
				if err != nil {
					return nil, fmt.Errorf("read GT allele: %w", err)
				}
				if uint32(v)&intWidthMask(typeLo) == intEndOfVectorFor(typeLo) {
					break
				}
				alleles = append(alleles, v)
			}
			out[s] = SampleValue{Type: "String", Strs: []string{decodeGenotype(alleles)}}
		}
		return out, nil
	}
	switch def.Type {
	case "Integer":
		// Per-sample int32 values, with END_OF_VECTOR
		// marker dropped. For GT-as-int the size is
		// the field's max ploidy across samples.
		//
		// The spec's per-type MISSING sentinel for int32
		// is `0x80000000` (line 2304) — the first of the
		// eight reserved values just above MinInt32
		// (T4.22). When all values for a sample slot read
		// the sentinel, the slot is treated as Missing.
		// The integer width comes from the descriptor's type
		// nibble: int8, int16 or int32. htslib picks the
		// narrowest that fits, and so do we.
		missingU := intMissingFor(typeLo)
		endU := intEndOfVectorFor(typeLo)
		for s := 0; s < nSamples; s++ {
			ints := make([]int32, 0, size)
			missing := true
			for i := 0; i < size; i++ {
				v, err := d.readIntValue(typeLo)
				if err != nil {
					return nil, fmt.Errorf("read int: %w", err)
				}
				raw := uint32(v) & intWidthMask(typeLo)
				if raw == endU && size > 0 {
					continue
				}
				if raw == missingU {
					ints = append(ints, v)
					out[s].IntMissing = append(out[s].IntMissing, true)
					continue
				}
				missing = false
				ints = append(ints, v)
				out[s].IntMissing = append(out[s].IntMissing, false)
			}
			out[s].Type = "Integer"
			out[s].Ints = ints
			out[s].Missing = missing && size > 0
		}
	case "Float":
		// Per-type MISSING is signaling NaN `0x7F800001`
		// (T4.24, line 2333). END_OF_VECTOR is the
		// related `0x7F800002` (line 2325). When all size
		// values for a sample read the MISSING sentinel,
		// the slot is treated as Missing.
		for s := 0; s < nSamples; s++ {
			floats := make([]float32, 0, size)
			missing := true
			for i := 0; i < size; i++ {
				var b [4]byte
				if err := d.readPending(b[:]); err != nil {
					return nil, fmt.Errorf("read float32: %w", err)
				}
				bits := binary.LittleEndian.Uint32(b[:])
				if bits == bcfEndOfVectorFloat32 && size > 0 {
					continue
				}
				if bits == bcfFloatMissing && size > 0 {
					floats = append(floats, math.Float32frombits(bits))
					out[s].FloatMissing = append(out[s].FloatMissing, true)
					continue
				}
				missing = false
				floats = append(floats, math.Float32frombits(bits))
				out[s].FloatMissing = append(out[s].FloatMissing, false)
			}
			out[s].Type = "Float"
			out[s].Floats = floats
			out[s].Missing = missing && size > 0
		}
	case "Character":
		// A Character FORMAT value is the raw character, NUL
		// padded to the field size.
		for s := 0; s < nSamples; s++ {
			buf := make([]byte, size)
			if err := d.readPending(buf); err != nil {
				return nil, fmt.Errorf("read char: %w", err)
			}
			text := strings.TrimRight(string(buf), "\x00")
			out[s] = SampleValue{Type: "Character", Str: text}
		}
	case "String":
		// A FORMAT string is the raw field text - the
		// comma-joined elements - NUL padded to the field
		// size. The commas are part of the text, so the
		// element boundaries survive the round trip.
		for s := 0; s < nSamples; s++ {
			buf := make([]byte, size)
			if err := d.readPending(buf); err != nil {
				return nil, fmt.Errorf("read string body: %w", err)
			}
			text := strings.TrimRight(string(buf), "\x00")
			var strs []string
			if text != "" {
				strs = strings.Split(text, ",")
			}
			out[s] = SampleValue{Type: "String", Strs: strs}
		}
	default:
		return nil, fmt.Errorf("unsupported FORMAT type %q", def.Type)
	}
	_ = typeLo
	return out, nil
}

// int64Slice converts a []int32 to a []int64 so the
// shared checkBcfIntegerRange helper accepts the FORMAT
// Ints slice. A copy avoids aliasing of the caller's
// storage.
func int64Slice(xs []int32) []int64 {
	out := make([]int64, len(xs))
	for i, v := range xs {
		out[i] = int64(v)
	}
	return out
}

// writeFormatBlock emits the per-field FORMAT block (§5
// line 2208–2255) for a record. The encoder's caller is
// expected to have already written the n_fmt uint16
// header indices; writeFormatBlock writes only the per-
// field type byte + N values for each.
//
//	type byte         (vector size in high nibble; type in low)
//	typed values      (one value per sample, padded with
//	                   END_OF_VECTOR for shorter entries)
//
// When all samples are missing for a field, the whole
// field collapses to a single typeMissing byte.
// writeFormatField emits one FORMAT field's fmt_type descriptor
// followed by its fmt_values (one value per sample). The
// fmt_key that precedes it is written by encodeIndivid.
func writeFormatField(w *bytes.Buffer, rec *Record, hdr *Header, version uint8, keyIdx uint16, pos int) error {
	nSamples := len(rec.Samples)
	def, ok := hdr.FormatAt(int(keyIdx))
	if !ok {
		return fmt.Errorf("bcf: FORMAT key index %d not in header", keyIdx)
	}
	// Vector size = max cardinality across samples
	// for this field. Atomic fields (Number=1) →
	// size 1. Variable-length fields (Number=.,
	// A, R, G) → max across samples.
	fieldSize, hasValues := fieldMaxSize(def, rec, keyIdx)
	if def.Type == "String" || def.Type == "Character" {
		// htslib leaves room for a NUL terminator in every
		// FORMAT field except the first one in the record
		// (vcf.c vcf_parse_format_max3 counts the field
		// length from the preceding delimiter, so every
		// field but the first is one byte longer). This is
		// an htslib implementation detail rather than a
		// spec requirement, but matching it keeps our
		// output byte-identical to the reference.
		if pos > 0 {
			fieldSize++
		}
	}
	// Emit the type byte (size in high nibble). We
	// always emit at least size 0 even when nothing
	// is present so the decoder can read past the
	// field cleanly.
	tb := formatTypeByte(def, version, hasValues)
	if def.Type == "Integer" && def.ID != "GT" {
		// The descriptor's type nibble must match the width
		// the values are actually written in.
		var vals []int32
		for s := range rec.Samples {
			v, present := rec.Samples[s].Fields[keyIdx]
			if present && !v.Missing {
				vals = v.Ints
				break
			}
		}
		tb = intWidthFor32(vals)
	}
	if fieldSize >= 15 {
		// Overflow: descriptor nibble 15 followed by a
		// typed Integer holding the true vector length
		// (spec §BcfTypeEncoding). Without this the reader
		// cannot tell a 15+ element vector from a short one.
		w.WriteByte(typeByteSize(15, tb))
		w.WriteByte(typeByteSize(1, typeInt8))
		w.WriteByte(byte(fieldSize)) //nolint:gosec // fieldSize < 256 here
	} else {
		w.WriteByte(typeByteSize(uint8(fieldSize), tb)) //nolint:gosec // guarded by fieldSize < 15
	}

	for s := 0; s < nSamples; s++ {
		v, present := rec.Samples[s].Fields[keyIdx]
		if err := writeFormatSampleValue(w, &v, def, version, fieldSize, present, tb); err != nil {
			return fmt.Errorf("bcf: FORMAT key %d (%s) sample %d: %w", keyIdx, def.ID, s, err)
		}
	}
	return nil
}

// fieldMaxSize returns the maximum cardinality observed
// for a FORMAT field across the record's samples, plus a
// has-values flag. Variable-length fields use each
// sample's slice length; atomic fields use the header's
// declared Number; Character fields use 1; String fields
// use each sample's Strs length.
func fieldMaxSize(def *FormatDef, rec *Record, keyIdx uint16) (int, bool) {
	count, isVar := ParseNumber(def.Number)
	// GT is an allele vector on the wire; its size is the
	// ploidy, so it must not take the String/Character path
	// below (which measures characters).
	if def.ID != "GT" && !isVar {
		// Atomic (Number=N) or Character/String with a
		// fixed declared count.
		switch def.Type {
		case "Character":
			return 1, true
		case "String":
			// A String FORMAT field is a *char vector*, so the
			// size nibble carries the character count, not the
			// element count. A Number=1 String holding 31
			// characters needs size 31.
			hasStr := false
			maxChars := 0
			for s := range rec.Samples {
				v, present := rec.Samples[s].Fields[keyIdx]
				if !present || v.Missing {
					continue
				}
				n := len(strings.Join(v.Strs, ","))
				if n > maxChars {
					maxChars = n
				}
				hasStr = true
			}
			if hasStr {
				return maxChars, true
			}
			if count == 1 {
				return 0, false
			}
		}
	}
	fieldMax := 0
	has := false
	for s := range rec.Samples {
		v, present := rec.Samples[s].Fields[keyIdx]
		if !present || v.Missing {
			// htslib still writes a slot for a sample with no
			// value: the per-type MISSING sentinel for
			// numbers, the literal "." text for strings and
			// characters.
			if !has {
				fieldMax = 1
				has = true
			}
			continue
		}
		has = true
		var n int
		switch {
		case def.ID == "GT":
			// GT is an allele vector on the wire; its size is
			// the ploidy, not the character count.
			if len(v.Ints) > 0 {
				n = len(v.Ints)
			} else if len(v.Strs) > 0 {
				n = len(encodeGenotype(v.Strs[0]))
			}
		case def.Type == "Integer":
			n = len(v.Ints)
		case def.Type == "Float":
			n = len(v.Floats)
		case def.Type == "Character":
			n = len(v.Str)
			if n == 0 {
				n = 1
			}
		case def.Type == "String":
			n = len(strings.Join(v.Strs, ","))
			if n == 0 {
				n = 1
			}
		}
		if n > fieldMax {
			fieldMax = n
		}
	}
	return fieldMax, has
}

// formatTypeByte returns the low-nibble type byte for the
// per-field type-byte header. Strings are returned as
// either typeString (7) or typeLong (0xF) depending on
// the longest string in this field, in BCF 2.2 mode.
// `version` is reserved for future expansion; the long
// string promotion is handled inside writeFormatBlock
// rather than here.
func formatTypeByte(def *FormatDef, _ uint8, _ bool) uint8 {
	// GT is stored in BCF as an integer allele vector even
	// though the VCF header declares it Type=String (spec
	// §GenotypeEncoding; htslib's bcf_format_gt_v2 rejects
	// any other type).
	if def.ID == "GT" {
		return typeInt8
	}
	switch def.Type {
	case "Integer":
		return typeInt32
	case "Float":
		return typeFloat
	case "Character":
		return typeString
	case "String":
		return typeString
	}
	return 0
}

// typeByteSize composes the per-field type byte: high
// nibble = vector size (1..14 inline; 15 = overflow),
// low nibble = type. We use inline-size 1..14 for
// atomic and known-length vectors; sizes 15+ are
// overflow-encoded.
func typeByteSize(size, lo uint8) uint8 {
	if size >= 15 {
		return (15 << 4) | lo
	}
	return (size << 4) | lo
}

// writeFormatSampleValue writes one sample's worth of
// typed values, padded with END_OF_VECTOR for shorter
// entries. A missing-or-absent sample pads the whole
// field with END_OF_VECTOR (so the stream stays
// uniform-length across samples).
func writeFormatSampleValue(w *bytes.Buffer, v *SampleValue, def *FormatDef, _ uint8, fieldSize int, present bool, width uint8) error {
	// Whole-field-missing collapses are only emitted
	// when the caller didn't provide per-slot missing
	// tracking. If IntMissing/FloatMissing are set,
	// the missing flags take precedence per slot
	// (VCF §1.6.2 line 508): a value like "." in a
	// Number=2 Integer list emits ".,.", not ".".
	if (!present || v.Missing) &&
		len(v.IntMissing) == 0 && len(v.FloatMissing) == 0 {
		switch def.Type {
		case "Integer":
			// Per-type MISSING: emit the sentinel for the
			// field's integer width (T4.24 / line 2304:
			// 0x80, 0x8000, 0x80000000). Matches the width
			// declared by the field type byte so the
			// reader's value loop stays aligned.
			//
			// An absent value is ONE missing element
			// followed by END_OF_VECTOR padding, which is
			// what htslib writes. An explicit ".,." is a
			// missing element per slot and is handled by
			// the per-slot path below.
			for i := 0; i < fieldSize; i++ {
				if i == 0 {
					writeIntValue(w, int32(intMissingFor(width)), width)
					continue
				}
				writeIntValue(w, int32(intEndOfVectorFor(width)), width)
			}
			return nil
		case "Float":
			// Per-type MISSING: emit the float32 signaling
			// NaN 0x7F800001 (T4.24 / line 2333), then
			// END_OF_VECTOR padding.
			for i := 0; i < fieldSize; i++ {
				if i == 0 {
					var b [4]byte
					binary.LittleEndian.PutUint32(b[:], bcfFloatMissing)
					w.Write(b[:])
					continue
				}
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], bcfEndOfVectorFloat32)
				w.Write(b[:])
			}
			return nil
		}
		// Character and String fall through: htslib writes the
		// literal "." text for a sample with no value, not a
		// MISSING sentinel.
	}
	if def.ID == "GT" {
		// GT is an integer allele vector on the wire, encoded
		// as int8 (htslib's bcf_int8).
		alleles := v.Ints
		if len(alleles) == 0 && len(v.Strs) > 0 {
			alleles = encodeGenotype(v.Strs[0])
		}
		for i := 0; i < fieldSize; i++ {
			if i < len(alleles) {
				w.WriteByte(byte(alleles[i]))
			} else {
				w.WriteByte(bcfEndOfVectorInt8)
			}
		}
		return nil
	}
	switch def.Type {
	case "Integer":
		if err := checkBcfIntegerRange(fmt.Sprintf("FORMAT/%s", def.ID), int64Slice(v.Ints)); err != nil {
			return err
		}
		// Pad with END_OF_VECTOR for shorter entries;
		// emit bcfIntMissingInt32 for per-slot missing.
		width := intWidthFor32(v.Ints)
		// Route through variables: these are untyped constants
		// whose high bit would overflow a constant int32
		// conversion at compile time.
		missingU := intMissingFor(width)
		endU := intEndOfVectorFor(width)
		missingVal := int32(missingU)
		endVal := int32(endU)
		for i := 0; i < fieldSize; i++ {
			if i < len(v.IntMissing) && v.IntMissing[i] {
				writeIntValue(w, missingVal, width)
				continue
			}
			if i < len(v.Ints) {
				writeIntValue(w, v.Ints[i], width)
			} else {
				writeIntValue(w, endVal, width)
			}
		}
		return nil
	case "Float":
		for i := 0; i < fieldSize; i++ {
			if i < len(v.FloatMissing) && v.FloatMissing[i] {
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], bcfFloatMissing)
				w.Write(b[:])
				continue
			}
			if i < len(v.Floats) {
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], math.Float32bits(v.Floats[i]))
				w.Write(b[:])
			} else {
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], bcfEndOfVectorFloat32)
				w.Write(b[:])
			}
		}
		return nil
	case "Character":
		// A Character FORMAT value is the raw character, NUL
		// padded to the field size. A sample with no value
		// gets the literal "." text, as htslib writes.
		text := "."
		if present && !v.Missing && v.Str != "" {
			text = v.Str[:1]
		}
		w.WriteString(text)
		writeNulPad(w, fieldSize-len(text))
		return nil
	case "String":
		// A FORMAT string is the raw field text - the
		// comma-joined elements - NUL padded to the field
		// size. htslib copies the text verbatim, commas
		// included, so the element boundaries survive. A
		// sample with no value gets the literal ".".
		text := "."
		if present && !v.Missing {
			text = strings.Join(v.Strs, ",")
			if text == "" {
				text = "."
			}
		}
		w.WriteString(text)
		writeNulPad(w, fieldSize-len(text))
		return nil
	}
	return fmt.Errorf("bcf: unknown FORMAT type %q", def.Type)
}

// writeNulPad writes n NUL bytes.
func writeNulPad(w *bytes.Buffer, n int) {
	for i := 0; i < n; i++ {
		w.WriteByte(0)
	}
}
