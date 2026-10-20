package bcf

import (
	"fmt"
	"strings"
)

// GT (genotype) is encoded on the BCF2 wire as a typed
// integer vector (§5 line 2463):
//
//	A Genotype (GT) field is encoded in a typed integer
//	vector (can be 8, 16, or even 32 bit if necessary)
//	with the number of elements equal to the maximum
//	ploidy among all samples at a site. For one
//	individual, each integer in the vector is organized
//	as $(allele+1) << 1 \mid phased$ where allele is set
//	to $-1$ if the allele in GT is a dot `.' (thus the
//	higher bits are all 0).
//
//	Examples (line 2474):
//	  0/1   → 0x02 0x04
//	  0|1   → 0x03 0x05
//	  ./.   → 0x00 0x00
//	  0/1/2 → 0x02 0x04 0x06
//
// The encoding uses 8-bit (int8) integers as long as the
// max ploidy * 16 fits comfortably — that is, for up to
// ~7 allele values. Bcftools and htslib promote to int16
// when the encoded value would overflow int8 (max
// allele+1 = 127); we use int8 by default and promote
// on overflow.

// encodeGT encodes one sample's GT (the VCF string form
// "0/1", "0|1", etc.) as the spec's int vector. The
// returned slice has length = GT ploidy. Alleles are
// stored as `(allele+1)<<1 | phased` per spec.
//
// Phasing rule (matching the spec's worked examples at
// lines 2474-2482):
//
//   - First allele: phased = (separator immediately
//     AFTER this allele is `|`) OR (single-value haploid,
//     which is implicitly phased).
//   - Non-first allele: phased = (separator immediately
//     BEFORE this allele is `|`).
//
// This produces the spec's examples:
//
//	"0/1"   → 0x02 0x04  (both unphased)
//	"0|1"   → 0x03 0x05  (first after `|`, second before `|`)
//	"./."   → 0x00 0x00  (missing always 0)
//	"0"     → 0x03      (single haploid → implicitly phased)
//	"0/1/2" → 0x02 0x04 0x06
//	"0/1|2" → 0x02 0x04 0x07  (middle allele unphased:
//	  separator before is `/`; trailing allele phased:
//	  separator before is `|`)
func encodeGT(gt string) ([]uint8, error) {
	if gt == "" {
		return nil, fmt.Errorf("bcf: empty GT")
	}
	tokens, sep, err := splitGT(gt)
	if err != nil {
		return nil, err
	}
	out := make([]uint8, len(tokens))
	for i, tok := range tokens {
		allele, missing, err := parseGTAllele(tok)
		if err != nil {
			return nil, err
		}
		var v int
		if missing {
			v = 0
		} else {
			v = (allele + 1) << 1
		}
		var phased bool
		switch {
		case i == 0 && len(tokens) == 1:
			// Single-allele haploid GT ("0", "1", ".")
			// is implicitly phased.
			phased = true
		case i == 0:
			// First allele of multi-allele GT: phased
			// iff separator AFTER this allele is `|`.
			phased = len(sep) > 0 && sep[0] == '|'
		default:
			// Non-first allele: phased iff separator
			// BEFORE this allele is `|`.
			phased = sep[i-1] == '|'
		}
		if phased {
			v |= 1
		}
		out[i] = uint8(v)
	}
	return out, nil
}

// decodeGT decodes one sample's GT int vector back to
// the VCF string form. The phasing bit is unfolded by
// emitting a `|` separator BEFORE a phased allele (and
// a `/` before an unphased one). The trailing allele
// never needs a separator in VCF.
//
// The error return is reserved for future use (e.g.,
// overflow handling when int16/int32 widths are added);
// today the function never returns a non-nil error.
//
//nolint:unparam // error reserved for future extension
func decodeGT(ints []uint8) (string, error) {
	if len(ints) == 0 {
		return "", nil
	}
	var sb strings.Builder
	for i, v := range ints {
		if i > 0 {
			// Separator between alleles i-1 and i is
			// determined by the phasing bit of allele i
			// (the second allele of the pair): `|` if
			// allele i is phased, else `/`. This mirrors
			// the encodeGT rule, which attaches the
			// phasing to the trailing allele of each
			// pair (the one whose separator-before
			// determines the join).
			if v&1 == 1 {
				sb.WriteByte('|')
			} else {
				sb.WriteByte('/')
			}
		}
		allele := int((v >> 1)) - 1
		if allele == -1 {
			sb.WriteByte('.')
		} else {
			fmt.Fprintf(&sb, "%d", allele)
		}
	}
	return sb.String(), nil
}

// splitGT tokenises a GT string ("0/1", "0|1|2", "./.")
// into the per-allele tokens plus the separators. The
// slice sep has length len(tokens)-1: sep[i] is the
// separator BETWEEN tokens[i] and tokens[i+1].
func splitGT(gt string) (tokens []string, sep []byte, err error) {
	start := 0
	for i := 0; i < len(gt); i++ {
		c := gt[i]
		if c == '/' || c == '|' {
			tokens = append(tokens, gt[start:i])
			sep = append(sep, c)
			start = i + 1
			continue
		}
	}
	tokens = append(tokens, gt[start:])
	if len(sep) != len(tokens)-1 {
		return nil, nil, fmt.Errorf("bcf: malformed GT %q", gt)
	}
	return tokens, sep, nil
}

// parseGTAllele decodes one GT allele token: "." →
// (0, true, nil) and "N" → (N, false, nil). Phased bit
// is captured separately by encodeGT.
func parseGTAllele(s string) (allele int, missing bool, err error) {
	if s == "." {
		return 0, true, nil
	}
	n := 0
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c < '0' || c > '9' {
			return 0, false, fmt.Errorf("bcf: GT allele %q is not an integer", s)
		}
		n = n*10 + int(c-'0')
		if n > 0xFFFF {
			// Beyond what int8 can carry once we add
			// (allele+1)<<1 and the phasing bit.
			return 0, false, fmt.Errorf("bcf: GT allele %q is too large for BCF2 encoding", s)
		}
	}
	return n, false, nil
}
