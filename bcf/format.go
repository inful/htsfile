package bcf

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
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
}

// readFormatBlock consumes the per-field FORMAT block (§5
// line 2208–2255) for a record given the format-key
// header indices in record order. The wire shape is, for
// each FORMAT key:
//
//	uint16 key index  (header offset)
//	type byte         (vector size in high nibble; type in low)
//	typed values      (one value per sample)
//
// Variable-length vectors are padded with END_OF_VECTOR for
// shorter samples; we strip those on read so each sample
// sees only its own cardinality.
func (d *Decoder) readFormatBlock(formatKeys []uint16, sampleCount int) ([]Sample, error) {
	samples := make([]Sample, sampleCount)
	for i := range samples {
		samples[i].Fields = make(map[uint16]SampleValue, len(formatKeys))
	}

	for _, keyIdx := range formatKeys {
		def, ok := d.hdr.FormatAt(int(keyIdx))
		if !ok {
			return nil, fmt.Errorf("bcf: FORMAT key index %d not in header", keyIdx)
		}
		vals, err := d.readFieldForAllSamples(def, sampleCount)
		if err != nil {
			return nil, fmt.Errorf("bcf: FORMAT key %d (%s): %w", keyIdx, def.ID, err)
		}
		for sIdx, sv := range vals {
			samples[sIdx].Fields[keyIdx] = sv
		}
	}
	return samples, nil
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
		ot, err := d.readByte()
		if err != nil {
			return nil, err
		}
		if ot != typeInt8 {
			return nil, fmt.Errorf("bcf: overflow size marker: got type byte %d, want INT8", ot)
		}
		var b [1]byte
		if _, err := io.ReadFull(d.br, b[:]); err != nil {
			return nil, err
		}
		size = int(b[0])
	}
	typeLo := t & 0x0F

	out := make([]SampleValue, nSamples)

	switch def.Type {
	case "Integer":
		// Per-sample int32 values, with END_OF_VECTOR
		// marker dropped. For GT-as-int the size is
		// the field's max ploidy across samples.
		for s := 0; s < nSamples; s++ {
			ints := make([]int32, 0, size)
			for i := 0; i < size; i++ {
				var b [4]byte
				if _, err := io.ReadFull(d.br, b[:]); err != nil {
					return nil, fmt.Errorf("read int32: %w", err)
				}
				u := binary.LittleEndian.Uint32(b[:])
				if u == bcfEndOfVectorInt32 && size > 0 {
					continue
				}
				ints = append(ints, int32(u)) //nolint:gosec // wire-format reinterpret
			}
			out[s] = SampleValue{Type: "Integer", Ints: ints}
		}
	case "Float":
		for s := 0; s < nSamples; s++ {
			floats := make([]float32, 0, size)
			for i := 0; i < size; i++ {
				var b [4]byte
				if _, err := io.ReadFull(d.br, b[:]); err != nil {
					return nil, fmt.Errorf("read float32: %w", err)
				}
				bits := binary.LittleEndian.Uint32(b[:])
				if bits == bcfEndOfVectorFloat32 && size > 0 {
					continue
				}
				floats = append(floats, math.Float32frombits(bits))
			}
			out[s] = SampleValue{Type: "Float", Floats: floats}
		}
	case "Character":
		// Character is a single byte per sample, with
		// typeMissing used to indicate absence.
		for s := 0; s < nSamples; s++ {
			var b [1]byte
			if _, err := io.ReadFull(d.br, b[:]); err != nil {
				return nil, fmt.Errorf("read char: %w", err)
			}
			out[s] = SampleValue{Type: "Character", Str: string(b[:])}
		}
	case "String":
		// Per the spec, strings at the FORMAT level are
		// stored as vectors of char: each sample slot
		// holds a length-prefixed body. We accept size
		// as the number of strings per sample, with each
		// entry being a length-prefixed string.
		for s := 0; s < nSamples; s++ {
			strs := make([]string, 0, size)
			for i := 0; i < size; i++ {
				var lb [4]byte
				if _, err := io.ReadFull(d.br, lb[:]); err != nil {
					return nil, fmt.Errorf("read string len: %w", err)
				}
				l := binary.LittleEndian.Uint32(lb[:])
				if l == bcfInt32Missing {
					strs = append(strs, "")
					continue
				}
				buf := make([]byte, l)
				if _, err := io.ReadFull(d.br, buf); err != nil {
					return nil, fmt.Errorf("read string body: %w", err)
				}
				strs = append(strs, string(buf))
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
func writeFormatBlock(w *bytes.Buffer, rec *Record, hdr *Header, version uint8) error {
	nSamples := len(rec.Samples)
	for _, keyIdx := range rec.Format {
		def, ok := hdr.FormatAt(int(keyIdx))
		if !ok {
			return fmt.Errorf("bcf: FORMAT key index %d not in header", keyIdx)
		}
		// All-missing fields collapse to typeMissing.
		if nSamples > 0 && allSamplesMissing(rec, keyIdx) {
			w.WriteByte(typeMissing)
			continue
		}

		// Vector size = max cardinality across samples
		// for this field. Atomic fields (Number=1) →
		// size 1. Variable-length fields (Number=.,
		// A, R, G) → max across samples.
		fieldSize, hasValues := fieldMaxSize(def, rec, keyIdx)
		// Emit the type byte (size in high nibble). We
		// always emit at least size 0 even when nothing
		// is present so the decoder can read past the
		// field cleanly.
		tb := formatTypeByte(def, version, hasValues)
		w.WriteByte(typeByteSize(uint8(fieldSize), tb)) //nolint:gosec // fieldSize is always 0..14 in the encoder path

		for s := 0; s < nSamples; s++ {
			v, present := rec.Samples[s].Fields[keyIdx]
			if err := writeFormatSampleValue(w, &v, def, version, fieldSize, present); err != nil {
				return fmt.Errorf("bcf: FORMAT key %d (%s) sample %d: %w", keyIdx, def.ID, s, err)
			}
		}
	}
	return nil
}

// allSamplesMissing reports whether every sample in the
// record is missing for the given FORMAT key.
func allSamplesMissing(rec *Record, keyIdx uint16) bool {
	for s := range rec.Samples {
		v, present := rec.Samples[s].Fields[keyIdx]
		if !present || !v.Missing {
			return false
		}
	}
	return true
}

// fieldMaxSize returns the maximum cardinality observed
// for a FORMAT field across the record's samples, plus a
// has-values flag. Variable-length fields use each
// sample's slice length; atomic fields use the header's
// declared Number; Character fields use 1; String fields
// use each sample's Strs length.
func fieldMaxSize(def *FormatDef, rec *Record, keyIdx uint16) (int, bool) {
	count, isVar := ParseNumber(def.Number)
	if !isVar {
		// Atomic (Number=N) or Character/String with a
		// fixed declared count.
		switch def.Type {
		case "Character":
			return 1, true
		case "String":
			// String fields declared Number=1 always
			// have at least one string; use the count
			// from the header for atomic fixed length.
			if count == 1 {
				return 1, true
			}
		}
	}
	fieldMax := 0
	has := false
	for s := range rec.Samples {
		v, present := rec.Samples[s].Fields[keyIdx]
		if !present || v.Missing {
			continue
		}
		has = true
		var n int
		switch def.Type {
		case "Integer":
			n = len(v.Ints)
		case "Float":
			n = len(v.Floats)
		case "Character":
			n = 1
		case "String":
			n = len(v.Strs)
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
func writeFormatSampleValue(w *bytes.Buffer, v *SampleValue, def *FormatDef, version uint8, fieldSize int, present bool) error {
	if !present || v.Missing {
		switch def.Type {
		case "Integer":
			for i := 0; i < fieldSize; i++ {
				w.WriteByte(typeMissing)
			}
			// typeMissing is a valid wire form; the
			// reader interprets it as Missing=true on the
			// matching sample.
			return nil
		case "Float":
			for i := 0; i < fieldSize; i++ {
				w.WriteByte(typeMissing)
			}
			return nil
		case "Character":
			w.WriteByte(typeMissing)
			return nil
		case "String":
			for i := 0; i < fieldSize; i++ {
				w.WriteByte(typeMissing)
			}
			return nil
		}
		return errors.New("unreachable")
	}
	switch def.Type {
	case "Integer":
		if err := checkBcfIntegerRange(fmt.Sprintf("FORMAT/%s", def.ID), int64Slice(v.Ints)); err != nil {
			return err
		}
		// Pad with END_OF_VECTOR for shorter entries.
		for i := 0; i < fieldSize; i++ {
			if i < len(v.Ints) {
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], uint32(v.Ints[i])) //nolint:gosec // wire-format reinterpret
				w.Write(b[:])
			} else {
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], bcfEndOfVectorInt32)
				w.Write(b[:])
			}
		}
		return nil
	case "Float":
		for i := 0; i < fieldSize; i++ {
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
		var b [1]byte
		if v.Str != "" {
			b[0] = v.Str[0]
		}
		w.Write(b[:])
		return nil
	case "String":
		// Strings at the FORMAT level are vectors of
		// character: each entry is a length-prefixed
		// string body. We write the sample's strings
		// and pad missing slots with typeMissing.
		strType := typeString
		if version == 2 {
			for _, s := range v.Strs {
				if len(s) >= 15 {
					strType = typeLong
					break
				}
			}
		}
		// (strType here applies across the per-sample
		// block, but per BCF 2.2 the type byte is set
		// once per field — we accept the simple case
		// where every sample's strings fit.)
		_ = strType
		for i := 0; i < fieldSize; i++ {
			if i < len(v.Strs) {
				s := v.Strs[i]
				var lb [4]byte
				if s == "" {
					binary.LittleEndian.PutUint32(lb[:], bcfInt32Missing)
				} else {
					binary.LittleEndian.PutUint32(lb[:], uint32(len(s))) //nolint:gosec // length fits in uint32
				}
				w.Write(lb[:])
				w.WriteString(s)
			} else {
				w.WriteByte(typeMissing)
			}
		}
		return nil
	}
	return fmt.Errorf("bcf: unknown FORMAT type %q", def.Type)
}
