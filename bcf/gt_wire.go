package bcf

import (
	"strconv"
	"strings"
)

// BCF2 stores the GT field as an integer allele vector, not as
// the VCF genotype string. Each allele is encoded as
//
//	(allele+1) << 1 | phased
//
// where a missing allele is 0 and `phased` is the phasing of
// the separator *before* this allele. The first allele's bit
// is the "prefixed phasing" flag.
//
// The rules below are a direct transcription of htslib's
// vcf_parse_format() GT branch (htslib 1.24, vcf.c ~line
// 3437), verified byte-for-byte against `bcftools view -Ob`
// output for "0|1", "0/1", "1|0", "1/0", "./.", "0", "0|1/2"
// and "0/1|2".
//
// encodeGenotype converts a VCF genotype string to the BCF
// allele vector.
func encodeGenotype(gt string) []int32 {
	if gt == "" {
		return nil
	}
	// A leading '|' or '/' is VCF 4.4 "prefixed phasing".
	prefixed := false
	phased := false
	t := gt
	if t[0] == '|' || t[0] == '/' {
		phased = t[0] == '|'
		prefixed = true
		t = t[1:]
	}
	out := []int32{}
	ploidy := 0
	anyUnphased := false
	unknown1 := false
	for {
		ploidy++
		if t == "" {
			break
		}
		if t[0] == '.' {
			t = t[1:]
			out = append(out, int32(boolToInt(phased)))
			if len(out) == 1 {
				unknown1 = true
			}
		} else {
			// Read the allele index.
			i := 0
			for i < len(t) && t[i] >= '0' && t[i] <= '9' {
				i++
			}
			if i == 0 {
				// Not a number: treat as missing.
				out = append(out, int32(boolToInt(phased)))
				if len(out) == 1 {
					unknown1 = true
				}
			} else {
				v, _ := strconv.Atoi(t[:i])
				t = t[i:]
				out = append(out, int32((v+1)<<1|boolToInt(phased)))
			}
		}
		if ploidy != 1 && !phased {
			anyUnphased = true
		}
		phased = strings.HasPrefix(t, "|")
		if !strings.HasPrefix(t, "|") && !strings.HasPrefix(t, "/") {
			break
		}
		t = t[1:]
	}
	if !prefixed {
		if ploidy == 1 {
			// Haploid: implicitly phased unless the allele is
			// missing.
			if !unknown1 {
				out[0] |= 1
			}
		} else {
			// The first allele is phased only if every other
			// allele is phased.
			if !anyUnphased {
				out[0] |= 1
			}
		}
	}
	return out
}

// decodeGenotype converts a BCF allele vector back to a VCF
// genotype string. It mirrors htslib's bcf_format_gt_v2.
func decodeGenotype(alleles []int32) string {
	if len(alleles) == 0 {
		return "."
	}
	var sb strings.Builder
	anyUnphased := false
	for i, v := range alleles {
		if i > 0 {
			if v&1 == 1 {
				sb.WriteByte('|')
			} else {
				sb.WriteByte('/')
				anyUnphased = true
			}
		}
		if v>>1 == 0 {
			sb.WriteByte('.')
		} else {
			sb.WriteString(strconv.Itoa(int(v>>1) - 1))
		}
	}
	val0 := alleles[0]
	ploidy := len(alleles)
	switch val0 & 1 {
	case 1:
		// The first allele is phased. VCF 4.4+ needs an
		// explicit leading '|' when some other allele is
		// unphased, or when the genotype is a haploid
		// missing call.
		if (ploidy > 1 && anyUnphased) || (ploidy <= 1 && val0>>1 == 0) {
			sb.WriteByte('|')
			s := sb.String()
			return s[1:]
		}
	case 0:
		// The first allele is unphased. VCF 4.4+ needs an
		// explicit leading '/' for a haploid call or when
		// every allele is phased.
		if (ploidy <= 1 && val0 != 0) || (ploidy > 1 && !anyUnphased) {
			sb.WriteByte('/')
			s := sb.String()
			return s[1:]
		}
	}
	return sb.String()
}

// boolToInt converts a bool to 0 or 1.
func boolToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}
