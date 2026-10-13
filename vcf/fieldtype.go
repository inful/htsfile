package vcf

import "fmt"

// FieldType is the data type of an INFO or FORMAT field, as declared in a
// ##INFO=<ID=...,Type=...,...> or ##FORMAT=<ID=...,Type=...,...> meta line.
//
// The spec lists Integer (32-bit signed), Float (32-bit IEEE-754), Flag,
// Character, and String.
type FieldType int8

const (
	// TypeUnknown is the zero value, used when a field type is missing or
	// could not be parsed. [Validator] treats TypeUnknown fields as opaque.
	TypeUnknown FieldType = iota
	// TypeInteger — 32-bit signed integer.
	TypeInteger
	// TypeFloat — IEEE-754 32-bit float. VCF permits the textual values
	// "INF", "INFINITY", "NAN" (case insensitive).
	TypeFloat
	// TypeFlag — present/absent toggle; the field's Number must be 0 and
	// the record column contains the key alone (e.g. "DB") or the missing
	// value ".".
	TypeFlag
	// TypeCharacter — single ASCII character.
	TypeCharacter
	// TypeString — free-form text. Reserved characters must be percent-
	// encoded on the wire.
	TypeString
)

// String returns the spec text representation of t.
func (t FieldType) String() string {
	switch t {
	case TypeInteger:
		return "Integer"
	case TypeFloat:
		return "Float"
	case TypeFlag:
		return "Flag"
	case TypeCharacter:
		return "Character"
	case TypeString:
		return "String"
	default:
		return ""
	}
}

// ParseFieldType parses a VCF Type= field.
func ParseFieldType(s string) (FieldType, error) {
	switch s {
	case "Integer":
		return TypeInteger, nil
	case "Float":
		return TypeFloat, nil
	case "Flag":
		return TypeFlag, nil
	case "Character":
		return TypeCharacter, nil
	case "String":
		return TypeString, nil
	}
	return TypeUnknown, fmt.Errorf("vcf: invalid Type %q (must be one of Integer, Float, Flag, Character, String)", s)
}
