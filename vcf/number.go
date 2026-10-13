package vcf

import (
	"fmt"
	"strconv"
)

// Number describes the cardinality of an INFO or FORMAT field.
//
// A [Number] is one of:
//
//   - a non-negative integer literal count (0, 1, 2, …), where 0 is reserved
//     for [TypeFlag] fields that do not carry a value;
//   - one of the symbolic cardinalities ([SymbolA], [SymbolR], [SymbolG],
//     [SymbolDot], [SymbolLA], [SymbolLR], [SymbolLG]), which expand at
//     validation time based on record context.
//
// See section 1.4.2 (Information field format) and section 1.4.4 (Individual
// format field field format) of VCFv4.5 for the original definitions.
type Number struct {
	// Fixed holds the symbolic cardinality when non-zero.
	Fixed Symbol
	// Number holds the literal count when Fixed is zero.
	Number int
}

// Symbol enumerates the symbolic cardinalities defined by the spec.
type Symbol int8

const (
	// SymbolNone is the zero value of [Symbol]; it indicates the [Number]
	// is a literal count rather than a symbolic cardinality.
	SymbolNone Symbol = iota
	// SymbolA — one value per alternate allele. The values must appear in
	// the same order as the ALT column.
	SymbolA
	// SymbolR — one value per possible allele (REF followed by each ALT).
	SymbolR
	// SymbolG — one value per possible genotype, ordered as described in
	// section 1.6.1.2 (Genotype Ordering).
	SymbolG
	// SymbolDot — the number of possible values varies, is unknown, or is
	// unbounded.
	SymbolDot
	// SymbolLA — identical to A but only the alternate alleles listed in
	// the FORMAT/LAA field are considered present (local-allele variant).
	SymbolLA
	// SymbolLR — identical to R for local alleles.
	SymbolLR
	// SymbolLG — identical to G for local alleles.
	SymbolLG
)

// String returns the spec text representation of n.
func (n Number) String() string {
	if n.Fixed == SymbolNone {
		return strconv.Itoa(n.Number)
	}
	switch n.Fixed {
	case SymbolDot:
		return "."
	case SymbolA:
		return "A"
	case SymbolR:
		return "R"
	case SymbolG:
		return "G"
	case SymbolLA:
		return "LA"
	case SymbolLR:
		return "LR"
	case SymbolLG:
		return "LG"
	}
	return ""
}

// ParseNumber parses a VCF Number= field.
func ParseNumber(s string) (Number, error) {
	switch s {
	case "A":
		return Number{Fixed: SymbolA}, nil
	case "R":
		return Number{Fixed: SymbolR}, nil
	case "G":
		return Number{Fixed: SymbolG}, nil
	case ".":
		return Number{Fixed: SymbolDot}, nil
	case "LA":
		return Number{Fixed: SymbolLA}, nil
	case "LR":
		return Number{Fixed: SymbolLR}, nil
	case "LG":
		return Number{Fixed: SymbolLG}, nil
	case "P", "M":
		// P (FORMAT only): one value per allele in GT.
		// M (FORMAT only): variable, depends on genotype + allele bases.
		// Defer detailed interpretation; treat as variable for cardinality
		// checks.
		return Number{Fixed: SymbolDot}, nil
	case "0":
		return Number{Number: 0}, nil
	}
	n, err := strconv.Atoi(s)
	if err != nil {
		return Number{}, fmt.Errorf("vcf: invalid Number %q: %w", s, err)
	}
	if n < 0 {
		return Number{}, fmt.Errorf("vcf: invalid Number %q (must be non-negative)", s)
	}
	return Number{Number: n}, nil
}

// IsFixed reports whether n is a symbolic cardinality
// (A, R, G, ., LA, LR, LG).
func (n Number) IsFixed() bool { return n.Fixed != SymbolNone }

// IsLocal reports whether n is a local-allele cardinality.
func (n Number) IsLocal() bool {
	return n.Fixed == SymbolLA || n.Fixed == SymbolLR || n.Fixed == SymbolLG
}

// Cardinality returns the number of values that a field with this Number is
// expected to carry on a record with the given alternate-allele count and
// (for G) ploidy. For symbolic cardinalities whose cardinality depends on
// context that the caller does not supply, -1 is returned to mean "variable".
//
// altCount is the number of entries in the ALT column (0 means no ALT
// alleles; "." is recorded as one entry per spec section 1.4.5).
//
// ploidy is the per-sample ploidy derived from the GT field for that sample;
// pass 0 when no GT is available or when the cardinality does not depend on
// ploidy.
func (n Number) Cardinality(altCount, ploidy int) int {
	if !n.IsFixed() {
		return n.Number
	}
	switch n.Fixed {
	case SymbolA:
		return altCount
	case SymbolR:
		return altCount + 1 // REF + ALTs
	case SymbolG:
		if ploidy <= 0 {
			return -1
		}
		return genotypeCount(altCount+1, ploidy)
	case SymbolDot, SymbolLA, SymbolLR, SymbolLG:
		return -1
	}
	return -1
}

// genotypeCount returns the number of ordered genotypes with ploidy p and
// alleleCount alleles (0..N), per VCFv4.5 §"Genotype Ordering".
//
// The closed-form is C(alleleCount + ploidy − 1, ploidy).
func genotypeCount(alleleCount, ploidy int) int {
	if ploidy <= 0 || alleleCount <= 0 {
		return -1
	}
	n := alleleCount + ploidy - 1
	k := ploidy
	if k > n-k {
		k = n - k
	}
	num, den := 1, 1
	for i := 0; i < k; i++ {
		num *= n - i
		den *= i + 1
	}
	return num / den
}
