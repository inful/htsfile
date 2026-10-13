package vcf

// InfoValue is the typed representation of one INFO field's value on a
// record. The exact field that should be read depends on the declared
// [FieldType] for the INFO key (see [Header.LookupInfo]):
//
//   - [TypeInteger] → use Ints
//   - [TypeFloat]   → use Floats
//   - [TypeString]  → use Strs (a list of one entry for Number=1)
//   - [TypeCharacter] → use Chars (always exactly one entry)
//   - [TypeFlag]    → the field is either present ([Flag] = true) or absent
//     from the Info map.
//
// The MISSING value "." is represented by [InfoValue.Missing] = true.
//
// Raw preserves the original token(s) from the VCF file before type
// decoding. For Number=1 fields the slice has length 1; for Number=A and
// similar it has length matching the ALT allele count (or fewer if some
// values are missing).
type InfoValue struct {
	Type    FieldType
	Ints    []int64
	Floats  []float64
	Strs    []string
	Chars   []byte
	Flag    bool
	Missing bool
	Raw     []string
}

// HasValue reports whether v represents a present, non-missing value.
func (v *InfoValue) HasValue() bool {
	if v.Missing {
		return false
	}
	if v.Type == TypeFlag {
		return v.Flag
	}
	return true
}

// String returns the first string value, or "" if no string values are
// present. Intended for Callers that know the field is Number=1 and
// Type=String.
func (v *InfoValue) String() string {
	if len(v.Strs) == 0 {
		return ""
	}
	return v.Strs[0]
}

// Int returns the first integer value, or 0 if no integer values are present.
func (v *InfoValue) Int() int64 {
	if len(v.Ints) == 0 {
		return 0
	}
	return v.Ints[0]
}

// Float returns the first float value, or 0 if no float values are present.
func (v *InfoValue) Float() float64 {
	if len(v.Floats) == 0 {
		return 0
	}
	return v.Floats[0]
}

// FieldValue is the typed representation of one FORMAT field's value for a
// single sample. Its fields parallel those of [InfoValue] except that no
// Flag type is permitted for FORMAT columns (a FORMAT field is always
// present per-sample; absence is encoded by [Missing] = true or ".").
type FieldValue struct {
	Type    FieldType
	Ints    []int64
	Floats  []float64
	Strs    []string
	Chars   []byte
	Missing bool
	Raw     []string
}

// HasValue reports whether v represents a present, non-missing value.
func (v *FieldValue) HasValue() bool { return !v.Missing }

// String returns the first string value, or "" if no string values are
// present.
func (v *FieldValue) String() string {
	if len(v.Strs) == 0 {
		return ""
	}
	return v.Strs[0]
}

// Int returns the first integer value, or 0 if no integer values are present.
func (v *FieldValue) Int() int64 {
	if len(v.Ints) == 0 {
		return 0
	}
	return v.Ints[0]
}

// Float returns the first float value, or 0 if no float values are present.
func (v *FieldValue) Float() float64 {
	if len(v.Floats) == 0 {
		return 0
	}
	return v.Floats[0]
}
