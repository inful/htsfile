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
	// keys.
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

// readSample reads one per-sample block from the BCF stream
// given the FORMAT layout (the list of format-key header
// indices for this record). The sample populates s.Fields;
// existing fields are not cleared.
func (d *Decoder) readSample(formatKeys []uint16) (*Sample, error) {
	s := &Sample{Fields: make(map[uint16]SampleValue, len(formatKeys))}
	for _, keyIdx := range formatKeys {
		def, ok := d.hdr.FormatAt(int(keyIdx))
		if !ok {
			return nil, fmt.Errorf("bcf: FORMAT key index %d not in header", keyIdx)
		}
		val, err := d.readSampleTypedValue(def)
		if err != nil {
			return nil, fmt.Errorf("bcf: FORMAT key %d (%s): %w", keyIdx, def.ID, err)
		}
		s.Fields[keyIdx] = val
	}
	return s, nil
}

// readSampleTypedValue reads one typed value for a per-sample
// FORMAT field, honouring the FORMAT key's declared type and
// Number. The cardinality rules are:
//
//   - Fixed (Number=1, 2, 3, …): read exactly that many
//     values of the declared type.
//   - Variable (Number=., A, R, G): read a uint32 count then
//     that many values of the declared type.
//
// GT (Number=G) and other ploidy-dependent fields are
// encoded with the variable-length form using the count
// implied by the sample's ploidy.
func (d *Decoder) readSampleTypedValue(def *FormatDef) (SampleValue, error) {
	_, isVar := ParseNumber(def.Number)
	var count int
	if isVar {
		var c [4]byte
		if _, err := io.ReadFull(d.br, c[:]); err != nil {
			return SampleValue{}, fmt.Errorf("read format count: %w", err)
		}
		count = int(binary.LittleEndian.Uint32(c[:]))
	} else {
		count, _ = ParseNumber(def.Number)
	}
	var val SampleValue
	t, err := d.readByte()
	if err != nil {
		return SampleValue{}, err
	}
	if t == typeMissing {
		val.Missing = true
		return val, nil
	}
	val.Type = def.Type
	switch def.Type {
	case "Integer":
		// Per-sample integers are 32-bit (the spec uses
		// int32, not int8, for FORMAT).
		val.Ints = make([]int32, count)
		for i := 0; i < count; i++ {
			var b [4]byte
			if _, err := io.ReadFull(d.br, b[:]); err != nil {
				return SampleValue{}, fmt.Errorf("read int32: %w", err)
			}
			val.Ints[i] = int32(binary.LittleEndian.Uint32(b[:])) //nolint:gosec // wire-format reinterpret
		}
	case "Float":
		val.Floats = make([]float32, count)
		for i := 0; i < count; i++ {
			var b [4]byte
			if _, err := io.ReadFull(d.br, b[:]); err != nil {
				return SampleValue{}, fmt.Errorf("read float32: %w", err)
			}
			val.Floats[i] = math.Float32frombits(binary.LittleEndian.Uint32(b[:]))
		}
	case "Character":
		// Character is always length 1; the type tag is
		// typeString with a 1-byte body.
		if t != typeString || count != 1 {
			return SampleValue{}, fmt.Errorf("character field: unexpected type=%d count=%d", t, count)
		}
		var b [1]byte
		if _, err := io.ReadFull(d.br, b[:]); err != nil {
			return SampleValue{}, fmt.Errorf("read char: %w", err)
		}
		val.Str = string(b[:])
	case "String":
		if t != typeString && t != typeLong {
			return SampleValue{}, fmt.Errorf("string field: unexpected type=%d", t)
		}
		val.Strs = make([]string, count)
		for i := 0; i < count; i++ {
			var lb [4]byte
			if _, err := io.ReadFull(d.br, lb[:]); err != nil {
				return SampleValue{}, fmt.Errorf("read string len: %w", err)
			}
			l := binary.LittleEndian.Uint32(lb[:])
			if l == bcfInt32Missing {
				val.Strs[i] = ""
				continue
			}
			buf := make([]byte, l)
			if _, err := io.ReadFull(d.br, buf); err != nil {
				return SampleValue{}, fmt.Errorf("read string body: %w", err)
			}
			val.Strs[i] = string(buf)
		}
	default:
		return SampleValue{}, fmt.Errorf("unsupported FORMAT type %q", def.Type)
	}
	return val, nil
}

// writeSample emits one per-sample block given the FORMAT
// layout. The version selects the long-string type tag
// (0xF) for BCF 2.2; pass 1 for BCF 2.1.
func writeSample(w *bytes.Buffer, s *Sample, formatKeys []uint16, hdr *Header, version uint8) error {
	for _, keyIdx := range formatKeys {
		def, ok := hdr.FormatAt(int(keyIdx))
		if !ok {
			return fmt.Errorf("bcf: FORMAT key index %d not in header", keyIdx)
		}
		val, present := s.Fields[keyIdx]
		if !present || val.Missing {
			w.WriteByte(typeMissing)
			continue
		}
		if err := writeSampleTypedValue(w, &val, def, version); err != nil {
			return fmt.Errorf("bcf: FORMAT key %d (%s): %w", keyIdx, def.ID, err)
		}
	}
	return nil
}

// writeSampleTypedValue is the encoder counterpart to
// readSampleTypedValue. The version selects the long-string
// type tag (0xF) for BCF 2.2.
func writeSampleTypedValue(w *bytes.Buffer, val *SampleValue, def *FormatDef, version uint8) error {
	_, isVar := ParseNumber(def.Number)
	if isVar {
		// Variable-length: emit a uint32 count followed by
		// that many values.
		var count int
		switch def.Type {
		case "Integer":
			count = len(val.Ints)
		case "Float":
			count = len(val.Floats)
		case "Character":
			count = 1
		case "String":
			count = len(val.Strs)
		default:
			return fmt.Errorf("unsupported FORMAT type %q", def.Type)
		}
		var cb [4]byte
		binary.LittleEndian.PutUint32(cb[:], uint32(count))
		w.Write(cb[:])
	}
	switch def.Type {
	case "Integer":
		w.WriteByte(typeInt32)
		for _, n := range val.Ints {
			var b [4]byte
			binary.LittleEndian.PutUint32(b[:], uint32(n)) //nolint:gosec // wire-format reinterpret
			w.Write(b[:])
		}
	case "Float":
		w.WriteByte(typeFloat)
		for _, f := range val.Floats {
			var b [4]byte
			binary.LittleEndian.PutUint32(b[:], math.Float32bits(f))
			w.Write(b[:])
		}
	case "Character":
		if val.Str == "" {
			w.WriteByte(typeMissing)
			return nil
		}
		w.WriteByte(typeString)
		w.WriteByte(val.Str[0])
	case "String":
		// For BCF 2.2, strings of length 15 or more use
		// the typeLong (0xF) tag instead of typeString (7).
		strType := typeString
		if version == 2 {
			for _, s := range val.Strs {
				if len(s) >= 15 {
					strType = typeLong
					break
				}
			}
		}
		w.WriteByte(strType)
		for _, s := range val.Strs {
			var lb [4]byte
			if s == "" {
				binary.LittleEndian.PutUint32(lb[:], bcfInt32Missing)
			} else {
				binary.LittleEndian.PutUint32(lb[:], uint32(len(s))) //nolint:gosec // length fits in uint32
			}
			w.Write(lb[:])
			w.WriteString(s)
		}
	default:
		return errors.New("bcf: unknown FORMAT type " + def.Type)
	}
	return nil
}
