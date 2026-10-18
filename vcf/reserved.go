package vcf

// reservedInfoKeys enumerates the VCF v4.5 reserved INFO keys
// whose expected Number/Type the Validator can compare against
// user-declared ##INFO lines. The full v4.5 spec lists many
// more; this table is the v4.5-new subset (base modifications,
// tandem repeats, local-allele-related). The other reserved
// INFO keys from earlier versions are already covered by the
// generic Number/Type validation path.
//
// Each entry: {key, number string, type code}. Number=.
// denotes "any" per spec (no cardinality check).
//
// Source: VCFv4.5 §1.6.1.4 (reserved keys), §1.6.1.5 (base
// modification), §1.6.1.6 (tandem repeats).
var reservedInfoKeys = map[string]reservedKey{
	// Base modification: M followed by a non-negative integer
	// and an optional single nucleotide (e.g. M0A, M5mC,
	// M5hmC). Number=1, Type=String. The key is treated
	// literally by the validator — callers can match the
	// pattern with reservedBaseModKeys below.
	// Tandem-repeat fields (VCFv4.5 §1.6.1.7 lines 1070
	// and 1084-1110): RN is per-allele; RUS/RUL/RUC/RB are
	// flat lists whose length is determined by RN; CIRUC /
	// CIRB are flat lists whose length is 2×RN.
	"RN":    {Number: "A", Type: TypeInteger},
	"RUS":   {Number: ".", Type: TypeString},
	"RUL":   {Number: ".", Type: TypeInteger},
	"RUC":   {Number: ".", Type: TypeFloat},
	"RB":    {Number: ".", Type: TypeInteger},
	"CIRUC": {Number: ".", Type: TypeFloat},
	"CIRB":  {Number: ".", Type: TypeInteger},
	"RUB":   {Number: ".", Type: TypeInteger},
	// Local-allele related (also appear as FORMAT keys below).
	"LAA": {Number: ".", Type: TypeString},
	"LAD": {Number: ".", Type: TypeString},
}

// reservedFormatKeys enumerates the VCF v4.5 reserved FORMAT
// keys. LAA is the only one validated structurally today
// (it must come right after GT; see Validator.checkSamples).
// The others are documented here so callers can detect them
// and the validator can check declared Number/Type against
// the spec.
var reservedFormatKeys = map[string]reservedKey{
	// Local-allele family (LAA must follow GT; rest may
	// follow any of the local-allele keys).
	"LAA":  {Number: ".", Type: TypeString},
	"LAD":  {Number: ".", Type: TypeString},
	"LADF": {Number: ".", Type: TypeString},
	"LADR": {Number: ".", Type: TypeString},
	"LEC":  {Number: ".", Type: TypeInteger},
	"LGL":  {Number: "G", Type: TypeString},
	"LGP":  {Number: "G", Type: TypeString},
	"LPL":  {Number: "P", Type: TypeString},
	"LPP":  {Number: "G", Type: TypeInteger},
	// Phase-set list trio (v4.5). All three use Number=P
	// (one value per allele value in GT).
	"PSL": {Number: "P", Type: TypeString},
	"PSO": {Number: "P", Type: TypeInteger},
	"PSQ": {Number: "P", Type: TypeInteger},
	// Genotype likelihood / posterior trio (per VCFv4.5
	// §1.6.2). All three use Number=G; PL is Integer,
	// GL and GP are Float.
	"PL": {Number: "G", Type: TypeInteger},
	"GL": {Number: "G", Type: TypeFloat},
	"GP": {Number: "G", Type: TypeFloat},
}

// reservedKey is the expected (Number, Type) pair for a
// reserved key. Number is a string ("1", ".", "G", "P", etc.)
// because the v4.5 spec expresses reserved-key cardinality
// symbolically and we don't want to convert at import time.
type reservedKey struct {
	Number string
	Type   FieldType
}

// IsReservedInfoKey reports whether key is a v4.5-reserved
// INFO key whose metadata the validator compares against.
// Matches the exact key names in reservedInfoKeys; the
// base-modification family (M0A, M5mC, M5hmC, …) is NOT
// covered by the literal lookup — use IsBaseModKey for
// that.
func IsReservedInfoKey(key string) bool {
	_, ok := reservedInfoKeys[key]
	return ok
}

// IsReservedFormatKey reports whether key is a v4.5-reserved
// FORMAT key.
func IsReservedFormatKey(key string) bool {
	_, ok := reservedFormatKeys[key]
	return ok
}

// MFamily returns the family prefix for a v4.5 base-modification
// FORMAT key: "M" for the fraction field, "DPM" for the
// total-read-depth field, "ADM" for the modified-read-depth
// field. Returns "" if the key does not match any of these
// patterns. The patterns are:
//
//	M<digits>[modifier]
//	DPM<digits>[modifier]
//	ADM<digits>[modifier]
//
// where <digits> is one or more decimal digits and [modifier]
// is the optional ChEBI-derived code (e.g. "5mC", "0A") or
// the trailing nucleotide (A/C/G/T/U/N).
func MFamily(key string) string {
	for _, prefix := range []string{"DPM", "ADM", "M"} {
		if len(key) <= len(prefix) {
			continue
		}
		if key[:len(prefix)] != prefix {
			continue
		}
		// After the prefix, require at least one digit and
		// no whitespace. The optional modifier code can be
		// anything non-whitespace (the v4.5 spec lists a
		// finite set; we are permissive).
		rest := key[len(prefix):]
		hasDigit := false
		for i := 0; i < len(rest); i++ {
			c := rest[i]
			if c == ' ' || c == '\t' || c == '\n' || c == '\r' {
				return ""
			}
			if !hasDigit && c >= '0' && c <= '9' {
				hasDigit = true
			}
		}
		if hasDigit {
			return prefix
		}
		return ""
	}
	return ""
}

// IsBaseModKey reports whether key matches the v4.5 base-
// modification pattern: "M" followed by a non-negative
// decimal integer and an optional single nucleotide code
// from the IUPAC alphabet (A, C, G, T, U, N) or one of the
// documented modification codes (5mC, 5hmC, 5fC, 5caC, 5hmU,
// 6mA, …). The validator does not enforce a specific
// modification code, only the shape.
//
// Examples that match: M0A, M5mC, M5hmC, M5fC, M6mA.
// Examples that do not match: MN, MR, M, MC (no digit).
func IsBaseModKey(key string) bool {
	if len(key) < 2 || key[0] != 'M' {
		return false
	}
	// Walk digits.
	i := 1
	for i < len(key) && key[i] >= '0' && key[i] <= '9' {
		i++
	}
	if i == 1 {
		return false // no digit after M
	}
	// Remainder is the modification code. It may be empty
	// (just the digit) or a non-empty code. We accept any
	// non-empty, non-whitespace tail here — the v4.5 spec
	// lists a finite set of codes, but the parser is
	// permissive.
	if i == len(key) {
		return true
	}
	for j := i; j < len(key); j++ {
		c := key[j]
		if c == ' ' || c == '\t' || c == '\n' || c == '\r' {
			return false
		}
	}
	return true
}