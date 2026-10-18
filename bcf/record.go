package bcf

import (
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
	// typeString marks a string payload (length-prefixed,
	// uint32 length). For BCF 2.1 the length is the only
	// prefix. For BCF 2.2 strings whose length exceeds 14
	// characters use typeLong (0xF) instead, with a uint32
	// length prefix.
	typeString uint8 = 7
	// typeLong marks a long string (BCF 2.2). The payload
	// is a uint32 length followed by that many bytes; this
	// form is required for strings whose length is 15 or
	// more, where the type 7 + length encoding runs out of
	// space in the type-byte's high nibble.
	typeLong uint8 = 0xF

	// bcfFilterMissing is the sentinel value stored in the
	// FILTER typed array when no filter applies (PASS or
	// unspecified).
	bcfFilterMissing uint32 = 0xFFFFFFFF

	// bcfInt32Missing is the sentinel value stored in
	// length-prefixed string fields when the string is
	// missing.
	bcfInt32Missing uint32 = 0xFFFFFFFF
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
