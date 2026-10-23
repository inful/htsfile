package vcf

// reservedInfoKeys enumerates the VCF v4.5 reserved INFO keys
// whose expected Number/Type the Validator compares against
// user-declared ##INFO lines.
//
// Two groups:
//
//  1. The Reserved INFO keys table (spec lines 460-480):
//     AA, AC, AD, ADF, ADR, AF, AN, BQ, CIGAR, DB, DP, END,
//     H2, H3, MQ, MQ0, NS, SB, SOMATIC, VALIDATED, 1000G.
//  2. The v4.5-new structural-variant keys (spec lines
//     1070-1110): RN, RUS, RUL, RUC, RB, CIRUC, CIRB, RUB,
//     plus the INFO-side local-allele fields LAA and LAD.
//
// Each entry: {key, number string, type code}. Number=.
// denotes "any" per spec (no cardinality check).
//
// Source: VCFv4.5 section 1.6.1 Reserved INFO keys table
// (lines 460-480) and section 1.6.1.6 (tandem repeats).
var reservedInfoKeys = map[string]reservedKey{
	// -- Reserved INFO keys table, spec lines 460-480. --
	"AA":        {Number: "1", Type: TypeString},
	"AC":        {Number: "A", Type: TypeInteger},
	"AD":        {Number: "R", Type: TypeInteger},
	"ADF":       {Number: "R", Type: TypeInteger},
	"ADR":       {Number: "R", Type: TypeInteger},
	"AF":        {Number: "A", Type: TypeFloat},
	"AN":        {Number: "1", Type: TypeInteger},
	"BQ":        {Number: "1", Type: TypeFloat},
	"CIGAR":     {Number: "A", Type: TypeString},
	"DB":        {Number: "0", Type: TypeFlag},
	"DP":        {Number: "1", Type: TypeInteger},
	"END":       {Number: "1", Type: TypeInteger},
	"H2":        {Number: "0", Type: TypeFlag},
	"H3":        {Number: "0", Type: TypeFlag},
	"MQ":        {Number: "1", Type: TypeFloat},
	"MQ0":       {Number: "1", Type: TypeInteger},
	"NS":        {Number: "1", Type: TypeInteger},
	"SB":        {Number: "4", Type: TypeInteger},
	"SOMATIC":   {Number: "0", Type: TypeFlag},
	"VALIDATED": {Number: "0", Type: TypeFlag},
	// 1000G is both a reserved Flag and the one INFO key the
	// spec's key regex (line 433) carves out as a legacy
	// exception: ^([A-Za-z_][0-9A-Za-z_.]*|1000G)$.
	"1000G": {Number: "0", Type: TypeFlag},

	// -- v4.5-new tandem-repeat fields (spec lines 1070 and
	//    1084-1110). RN is per-allele; RUS/RUL/RUC/RB are
	//    flat lists whose length is determined by RN;
	//    CIRUC/CIRB are flat lists whose length is twice RN. --
	"RN":    {Number: "A", Type: TypeInteger},
	"RUS":   {Number: ".", Type: TypeString},
	"RUL":   {Number: ".", Type: TypeInteger},
	"RUC":   {Number: ".", Type: TypeFloat},
	"RB":    {Number: ".", Type: TypeInteger},
	"CIRUC": {Number: ".", Type: TypeFloat},
	"CIRB":  {Number: ".", Type: TypeInteger},
	"RUB":   {Number: ".", Type: TypeInteger},
	// Local-allele related (also appear as FORMAT keys below).
	"LAA": {Number: ".", Type: TypeInteger},
	"LAD": {Number: "LR", Type: TypeInteger},
}

// reservedFormatKeys enumerates the VCF v4.5 reserved FORMAT
// keys. Two groups:
//
//  1. The Reserved genotype keys table (spec lines 534-596):
//     AD, ADF, ADR, DP, EC, LEN, FT, GL, GP, GQ, GT, HQ, LA,
//     LAA, LAD, LADF, LADR, LEC, LGL, LGP, LPL, LPP, MQ, PL,
//     PP, PQ, PS, PSL, PSO, PSQ.
//  2. The local-allele family (LAA, LAD, LADF, LADR, LEC,
//     LGL, LGP, LPL, LPP) and the phase-set list trio (PSL,
//     PSO, PSQ), which the v4.5 local-allele and phase-set
//     work introduced.
//
// The base-modification families (M, DPM, ADM - spec lines
// 556-588, all Number=M) are deliberately NOT literal entries:
// they are resolved by [MFamily] so the table stays finite.
var reservedFormatKeys = map[string]reservedKey{
	// -- Reserved genotype keys table, spec lines 534-596. --
	"AD":  {Number: "R", Type: TypeInteger},
	"ADF": {Number: "R", Type: TypeInteger},
	"ADR": {Number: "R", Type: TypeInteger},
	"DP":  {Number: "1", Type: TypeInteger},
	"EC":  {Number: "A", Type: TypeInteger},
	"LEN": {Number: "1", Type: TypeInteger},
	"FT":  {Number: "1", Type: TypeString},
	"GL":  {Number: "G", Type: TypeFloat},
	"GP":  {Number: "G", Type: TypeFloat},
	"GQ":  {Number: "1", Type: TypeInteger},
	"GT":  {Number: "1", Type: TypeString},
	"HQ":  {Number: "2", Type: TypeInteger},
	"LA":  {Number: ".", Type: TypeInteger},
	"MQ":  {Number: "1", Type: TypeInteger},
	"PL":  {Number: "G", Type: TypeInteger},
	"PP":  {Number: "G", Type: TypeInteger},
	"PQ":  {Number: "1", Type: TypeInteger},
	"PS":  {Number: "1", Type: TypeInteger},

	// -- Local-allele family. LAA must follow GT; the rest may
	//    follow any of the local-allele keys.
	//
	//    LAA is Integer, not String: spec line 547 reads
	//    `LAA & . & Integer & 1-based indices into ALT`.
	//    laaLenFor / checkLAA already accept either wire form
	//    (the values may parse as Ints or as Strs), so only
	//    the header-metadata cross-check needed correcting. --
	"LAA":  {Number: ".", Type: TypeInteger},
	"LAD":  {Number: "LR", Type: TypeInteger},
	"LADF": {Number: "LR", Type: TypeInteger},
	"LADR": {Number: "LR", Type: TypeInteger},
	"LEC":  {Number: "LA", Type: TypeInteger},
	"LGL":  {Number: "LG", Type: TypeFloat},
	"LGP":  {Number: "LG", Type: TypeFloat},
	// LPL and LPP are the local-allele forms of the
	// phred-scaled likelihood/posterior, so Number=LG with
	// Type=Integer - not LG Float (that is LGL/LGP) and not
	// Number=G. Spec lines 554-555.
	"LPL": {Number: "LG", Type: TypeInteger},
	"LPP": {Number: "LG", Type: TypeInteger},

	// -- Phase-set list trio (spec lines 594-596). All three
	//    use Number=P: one value per allele value in GT. --
	"PSL": {Number: "P", Type: TypeString},
	"PSO": {Number: "P", Type: TypeInteger},
	"PSQ": {Number: "P", Type: TypeInteger},
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

// ReservedFormatDef returns the reserved Number and Type of a
// v4.5-reserved FORMAT key, and whether the key is reserved.
// It lets a consumer that has to synthesise a header
// declaration (the BCF writer, which must declare every key a
// record uses) use the spec's own definitions instead of
// guessing.
func ReservedFormatDef(key string) (number string, fieldType FieldType, ok bool) {
	r, ok := reservedFormatKeys[key]
	if !ok {
		return "", TypeUnknown, false
	}
	return r.Number, r.Type, true
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
