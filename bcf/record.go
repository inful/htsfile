package bcf

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"strconv"
)

// Type tags for the typed arrays that appear in BCF records.
// The same byte is reused with different meanings across
// fields (e.g. ID type 1 means "string present"; QUAL type 1
// means "missing") — the parser is field-aware, not just
// byte-aware.
const (
	// typeMissing marks a typed field that has no value.
	typeMissing uint8 = 0
	// typeInt8 marks an int8 (or, in some contexts, a filter
	// index) payload.
	typeInt8 uint8 = 1
	// typeInt16 marks an int16 payload.
	typeInt16 uint8 = 2
	// typeInt32 marks an int32 payload.
	typeInt32 uint8 = 3
	// typeFloat marks a float32 / float64 payload (the
	// width is field-dependent: QUAL uses float64, INFO
	// fields use float32).
	typeFloat uint8 = 5
	// typeString marks a BCF2 string, which is encoded as a
	// *vector of chars*: the type descriptor byte carries the
	// length in its high nibble and there is no separate length
	// prefix. A 9-character ID is therefore 0x97 followed by the
	// 9 bytes - NOT 0x07 followed by a uint32 length, which a
	// conformant reader decodes as size 0 (MISSING).
	//
	// Spec §BcfTypeEncoding, tex line 2276: "In BCF2 strings are
	// type=7 with an appropriate length. Single characters are
	// represented as strings with length=1."
	typeString uint8 = 7

	// bcfFilterMissing is the sentinel value stored in the
	// FILTER typed array when no filter applies (PASS or
	// unspecified).
	bcfFilterMissing uint32 = 0xFFFFFFFF
)

// Record is one BCF alignment record. The shared per-record
// fields (CHROM, POS, ID, REF, ALT, FILTER, QUAL, INFO,
// FORMAT list, and per-sample data) are populated by the
// Decoder.
//
// Zero value is not meaningful; obtain Records from
// [Decoder.Record].
type Record struct {
	// Line is the 1-based record number within the alignment
	// section (the first BCF record is line 1).
	Line int

	// CHROM is the resolved contig name (header.Contigs[i]
	// for the int32 CHROM index in the binary record).
	CHROM string

	// Pos is the 1-based leftmost position. BCF stores 0-based
	// positions internally; the decoder adds 1.
	Pos int64

	// ID is the record identifier. Empty when missing.
	ID string

	// Ref is the reference allele. Empty when missing.
	Ref string

	// Alt is the list of alternate alleles. The length is the
	// number of ALTs; empty slice when missing.
	Alt []string

	// Qual is the phred-scaled quality score. NaN (specifically
	// the value stored as the "missing" sentinel) is surfaced
	// as 0 with QualKnown=false; use QualKnown to distinguish
	// "0.0" from "missing".
	Qual float64

	// QualKnown reports whether Qual was read from a real
	// value (true) or from the "missing" sentinel (false).
	QualKnown bool

	// Filter is the list of filter IDs (resolved names) for
	// this record. Empty when no filter applies (PASS or
	// unspecified).
	Filter []string

	// Info maps header INFO ID → typed value. The map is
	// built field-by-field as the BCF INFO typed array is
	// parsed; field types are coerced according to the
	// header's ##INFO declarations.
	Info map[string]InfoValue

	// Format is the list of FORMAT-key header indices
	// (uint16) in the order the encoder emitted / the
	// decoder read them. Each entry also appears as a key
	// in every Sample's Fields map. Empty when the record
	// has no per-sample data.
	Format []uint16

	// Samples is the per-sample data, one entry per sample
	// declared in the header (in header order). Empty when
	// the record has no per-sample data.
	Samples []Sample
}

// InfoValue holds the typed value(s) for a single BCF INFO
// field. BCF INFO values are arrays (Number may be 1, A, R,
// G, or an integer); the parser stores all values in a Go
// slice per the most natural Go type.
//
// The Type field records the VCF/BCF declared type
// (Integer / Float / Flag / Character / String).
type InfoValue struct {
	Type   string // "Integer", "Float", "Flag", "Character", "String"
	Ints   []int64
	Floats []float64
	Strs   []string // for String type (one per Number element)
	Str    string   // for Character type (single value)
	Flag   bool     // true when Type == "Flag" and the field is present
}

// String returns a SAM/VCF-style textual representation of
// the value (used by tests; not part of the public API in
// the same way as vcf.InfoValue.String).
func (v InfoValue) String() string {
	switch v.Type {
	case "Flag":
		if v.Flag {
			return "."
		}
		return ""
	case "Integer":
		out := ""
		for i, n := range v.Ints {
			if i > 0 {
				out += ","
			}
			out += strconv.FormatInt(n, 10)
		}
		return out
	case "Float":
		out := ""
		for i, f := range v.Floats {
			if i > 0 {
				out += ","
			}
			out += strconv.FormatFloat(f, 'g', -1, 32)
		}
		return out
	default:
		return v.Str
	}
}

// _ keeps the fmt import used (the parser calls
// fmt.Errorf in several places).
var _ = fmt.Errorf

// writeBCFString emits a BCF2 string as a char vector: the type
// descriptor byte carries the length in its high nibble
// (spec §BcfTypeEncoding, tex line 2276), followed by the raw
// bytes with no separate length prefix. Lengths of 15 or more
// use the overflow form - descriptor nibble 15 followed by a
// typed Integer holding the true length.
//
// This is the single correct place to encode a BCF string. The
// earlier implementation wrote `0x07` + a uint32 length, which a
// conformant reader decodes as size 0 (MISSING) followed by
// garbage; bcftools rejected such records with "Invalid REF/ALT
// type 9".
func writeBCFString(w *bytes.Buffer, s string) {
	n := len(s)
	if n == 0 {
		w.WriteByte(typeMissing)
		return
	}
	if n <= 14 {
		w.WriteByte(byte(n)<<4 | typeString)
		w.WriteString(s)
		return
	}
	// Overflow: size nibble 15, then a typed Integer length.
	w.WriteByte(15<<4 | typeString)
	switch {
	case n <= 0xFF:
		w.WriteByte(1<<4 | typeInt8)
		w.WriteByte(byte(n))
	case n <= 0xFFFF:
		w.WriteByte(1<<4 | typeInt16)
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(n))
		w.Write(b[:])
	default:
		w.WriteByte(1<<4 | typeInt32)
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(n)) //nolint:gosec // n <= 0xFFFF
		w.Write(b[:])
	}
	w.WriteString(s)
}
