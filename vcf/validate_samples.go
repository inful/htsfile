package vcf

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"
)

// architecture cleanup (item #1 of the gortex review). The
// per-sample/FORMAT block and the small helper layer
// (hasGT, alleleIndex, wantedCardinality, etc.) live here.
func (v *Validator) checkSamples(h *Header, rec *Record, idx int) []ValidationError {
	var issues []ValidationError
	line := rec.Line

	if len(rec.Samples) == 0 {
		return issues
	}

	// GT must be the first FORMAT key when present.
	for i, k := range rec.Format {
		if k == "GT" && i != 0 {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "FORMAT",
				Msg:   "GT must be the first FORMAT key when present",
			})
			break
		}
	}

	// LAA must precede all fields other than GT when present.
	for i, k := range rec.Format {
		if k == "LAA" && i > 1 {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "FORMAT",
				Msg:   "LAA must precede all FORMAT fields other than GT",
			})
			break
		}
	}

	// Per-sample checks.
	for sIdx, g := range rec.Samples {
		// PS and PSL are mutually exclusive per sample
		// (VCFv4.5 §1.6.2: "A given sample-genotype must
		// not have values for both PS and PSL").
		if hasPSandPSL(g) {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/PS,sample[%d]", sIdx),
				Msg:   "PS and PSL cannot both be present on the same sample",
			})
		}

		// FT per-sample filter must not contain whitespace
		// or semicolons (VCFv4.5 §1.6.2 FT field: "No
		// whitespace or semicolons permitted").
		for _, tok := range ftBadTokens(g) {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/FT/sample[%d]", sIdx),
				Msg:   fmt.Sprintf("FT token %q must not contain whitespace or ';'", tok),
			})
		}

		// GT structural checks (separator consistency, allele
		// index range) — emit one issue per problem so the
		// caller can see all errors at once.
		if v, ok := g.Fields["GT"]; ok && !v.Missing && len(v.Strs) > 0 {
			if err := checkGT(v.Strs[0], len(rec.Alt)); err != nil {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: fmt.Sprintf("FORMAT/GT/sample[%d]", sIdx),
					Msg:   err.Error(),
				})
			}
			// Haploid chromosomes (mtDNA) must have GT
			// with exactly one allele (§1.6.2:
			// "Haploid calls, e.g. on Y, male non-
			// pseudoautosomal X, or mitochondria, should
			// be indicated by having only one allele
			// value."). We only check mtDNA — Y is
			// ambiguous because PAR regions can carry
			// diploid calls.
			if isMitochondrial(rec.Chrom) && ploidyFromGT(g) != 1 {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: fmt.Sprintf("FORMAT/GT/sample[%d]", sIdx),
					Msg:   fmt.Sprintf("mitochondrial chromosome %s should have a haploid GT (one allele) but found ploidy %d", rec.Chrom, ploidyFromGT(g)),
				})
			}
		}

		// PSL per-allele check (§1.6.2: "Unphased alleles
		// (without a | separator before them) must have
		// the value '$.$' in their corresponding position
		// in the list"). For each unphased allele in GT,
		// the PSL entry must be exactly ".".
		issues = append(issues, checkPSLUnphased(g, sIdx, idx, line)...)

		// GP values must be in [0, 1] (§1.6.2 GP field).
		issues = append(issues, checkGPRange(g, sIdx, idx, line)...)

		// M-key present requires GT defined for the sample
		// (§1.6.1.5: "If any base modification key is present
		// for a sample, GT must be defined for that sample").
		issues = append(issues, checkMRequiresGT(g, sIdx, idx, line)...)

		// PSL / PSO / PSQ trio consistency per allele
		// (§1.6.2 PSQ: "The missing value '$.$' should be
		// used when the corresponding PSL value is missing").
		// When PSL[i] is ".", PSO[i] and PSQ[i] must also
		// be "."; the same holds in reverse.
		issues = append(issues, checkPSLTrio(g, sIdx, idx, line)...)

		// LAA must precede all FORMAT keys other than GT
		// (per §1.6.2: "If any local-allele field is present,
		// LAA must also be present and precede all fields
		// other than GT").
		if _, hasLAA := g.Fields["LAA"]; hasLAA {
			issues = append(issues, checkLAAPosition(rec, sIdx, idx, line)...)
		}
		// LAA value range and distinctness (per VCFv4.5
		// §1.6.1.4). LAA must be distinct 1-based indices
		// into ALT.
		if v, ok := g.Fields["LAA"]; ok && !v.Missing {
			issues = append(issues, checkLAA(&v, rec, idx, sIdx, line)...)
		}

		// If any local-allele field is present, LAA must
		// also be present (VCFv4.5 §1.6.2: "If any
		// local-allele field is present, LAA must also be
		// present"). The local-allele FORMAT keys are the
		// L-prefixed ones listed in the spec: LAD, LADF,
		// LADR, LEC, LGL, LGP, LPL, LPP. Per the spec, an
		// empty LAA is equivalent to missing for the
		// purpose of this rule, so the check fires only
		// when LAA is genuinely non-empty.
		for _, key := range missingLAAFields(g) {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/%s/sample[%d]", key, sIdx),
				Msg:   fmt.Sprintf("local-allele field %q requires LAA to be present on the same sample", key),
			})
		}

		// Determine the last present FORMAT index for this sample — keys
		// past that index are considered dropped (per spec).
		lastPresent := -1
		for i, key := range rec.Format {
			if _, ok := g.Fields[key]; ok {
				lastPresent = i
			}
		}

		for i, key := range rec.Format {
			fv, present := g.Fields[key]
			if !present {
				// Trailing keys past the last present one are dropped;
				// otherwise this is a truly missing field.
				if i > lastPresent {
					continue
				}
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: fmt.Sprintf("FORMAT/%s/sample[%d]", key, sIdx),
					Msg:   "FORMAT key missing from sample data",
				})
				continue
			}
			def, declared := h.Formats[key]
			if !declared {
				iss := ValidationError{
					Severity: SeverityWarning,
					Record:   idx + 1,
					Line:     line,
					Field:    fmt.Sprintf("FORMAT/%s/sample[%d]", key, sIdx),
					Msg:      "FORMAT key is not declared in the header",
				}
				if v.Mode == ModeStrict {
					iss.Severity = SeverityError
				}
				issues = append(issues, iss)
				continue
			}
			// Number checking: depends on ploidy for G / LG. Use GT to
			// derive ploidy when possible. Local-allele
			// FORMAT fields (LAD, LEC, LGL, …) derive their
			// cardinality from LAA when present, not from
			// the full ALT list (§1.6.2); without LAA they
			// fall back to ALT count.
			want := wantedCardinality(def.Number, key, g, rec)
			got := len(fv.Raw)
			if fv.Missing {
				got = 0
			}
			if want >= 0 && got != want {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: fmt.Sprintf("FORMAT/%s/sample[%d]", key, sIdx),
					Msg: fmt.Sprintf("expected %d value(s) for Number=%s, got %d",
						want, def.Number, got),
				})
			}
			// Type checking.
			for _, r := range fv.Raw {
				if !validTypedValue(def.Type, r) {
					issues = append(issues, ValidationError{
						Severity: SeverityError, Record: idx + 1, Line: line,
						Field: fmt.Sprintf("FORMAT/%s/sample[%d]", key, sIdx),
						Msg:   fmt.Sprintf("value %q does not match declared Type=%q", r, def.Type),
					})
				}
			}

			// Base-modification family value range (VCF v4.5).
			//   M   : float in [0, 1]
			//   DPM  : non-negative integer (read depth)
			//   ADM  : non-negative integer (read depth)
			// We rely on fv.Ints / fv.Floats, which decodeTyped
			// populates from the raw values.
			if fam := MFamily(key); fam != "" {
				if fam == "M" {
					for _, f := range fv.Floats {
						if f < 0 || f > 1 {
							issues = append(issues, ValidationError{
								Severity: SeverityError, Record: idx + 1, Line: line,
								Field: fmt.Sprintf("FORMAT/%s/sample[%d]", key, sIdx),
								Msg:   fmt.Sprintf("%s value %g is out of range [0, 1]", key, f),
							})
						}
					}
				} else {
					// DPM / ADM: non-negative integer.
					for _, n := range fv.Ints {
						if n < 0 {
							issues = append(issues, ValidationError{
								Severity: SeverityError, Record: idx + 1, Line: line,
								Field: fmt.Sprintf("FORMAT/%s/sample[%d]", key, sIdx),
								Msg:   fmt.Sprintf("%s value %d is negative", key, n),
							})
						}
					}
				}
			}
		}
	}
	return issues
}

// validRefBases reports whether s contains only A/C/G/T/N bases (case
// insensitive). An empty string is not valid.
func validRefBases(s string) bool {
	if s == "" {
		return false
	}
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch c {
		case 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n':
		default:
			return false
		}
	}
	return true
}

// validAltAllele reports whether a is a syntactically valid ALT allele.
// Acceptable forms: a non-empty ACGTN string, "*", "<ID>" where ID is any
// string without whitespace/angle brackets, "<*>", or a breakend
// replacement string.
func validAltAllele(a string) bool {
	if a == "" {
		return false
	}
	if a == "*" || a == "<*>" {
		return true
	}
	if strings.HasPrefix(a, "<") && strings.HasSuffix(a, ">") {
		// <ID> symbolic allele — ID must not contain whitespace or
		// angle brackets.
		inner := a[1 : len(a)-1]
		if inner == "" {
			return false
		}
		for i := 0; i < len(inner); i++ {
			c := inner[i]
			if c == '<' || c == '>' || c == ' ' || c == '\t' || c == '\n' || c == '\r' {
				return false
			}
		}
		return true
	}
	// Breakend: bracket-balanced form starting or ending with [ or ].
	if isBreakend(a) {
		return true
	}
	// Bare IUPAC ambiguity code (R, Y, S, W, K, M, B, D,
	// H, V) is a valid symbolic-allele ALT form per
	// VCFv4.5 §1.4.4.
	if iupacIDs[a] {
		return true
	}
	// Plain nucleotide string.
	return validRefBases(a)
}

// isBreakend reports whether a is a breakend replacement string of one of
// the four forms: t]p], t[p[, ]p]t, [p[t. Per VCFv4.5 §1.6.5.1, the
// position specifier is bracketed by two of the same kind of bracket character
// (`[` or `]`). Leading/trailing base strings ("t") are optional.
//
// We do a structural check: the string must contain two matching brackets
// of the same kind and the substring between them must include a colon
// (the position specifier uses "chr:pos" form).
func isBreakend(a string) bool {
	if len(a) < 4 {
		return false
	}
	// Find the first bracket character.
	start := -1
	for i := 0; i < len(a); i++ {
		if a[i] == '[' || a[i] == ']' {
			start = i
			break
		}
	}
	if start < 0 || start+2 >= len(a) {
		return false
	}
	open := a[start]
	// Find the matching closer (same bracket character) at a later position.
	for i := start + 1; i < len(a); i++ {
		if a[i] == open {
			// Position specifier sits between start+1 and i (exclusive).
			inner := a[start+1 : i]
			return strings.Contains(inner, ":")
		}
	}
	return false
}

// vcfFloatPattern is the VCF v4.5 §1.4.3 "Data types" regex for
// Float values:
//
//	^[-+]?\d*\.?\d+([eE][-+]?\d+)?$|^[-+]?(INF|INFINITY|NAN)$
//	case-insensitively
//
// The regex rejects "1." (trailing dot) and "1e" (trailing
// exponent with no digits) which strconv.ParseFloat accepts.
var vcfFloatPattern = regexp.MustCompile(`^[-+]?\d*\.?\d+([eE][-+]?\d+)?$`)

// validFloatString reports whether s matches the VCF Float
// format. The spec allows "INF", "INFINITY", "NAN" (and
// negative variants) case-insensitively; everything else
// must match vcfFloatPattern.
func validFloatString(s string) bool {
	upper := strings.ToUpper(s)
	switch upper {
	case "INF", "INFINITY", "-INF", "-INFINITY", "NAN":
		return true
	}
	return vcfFloatPattern.MatchString(s)
}

// validTypedValue reports whether s is a valid representation of a value
// of the declared type.
func validTypedValue(t FieldType, s string) bool {
	if s == "." {
		return true // MISSING is always permitted
	}
	switch t {
	case TypeInteger:
		_, err := strconv.ParseInt(s, 10, 32)
		return err == nil
	case TypeFloat:
		return validFloatString(s)
	case TypeFlag:
		return s == ""
	case TypeCharacter:
		return len(s) == 1
	case TypeString:
		return true
	}
	return true
}

// empty reports whether the field carries no values
// (missing value, zero entries, or a single empty
// string entry — all of which the spec treats as "no
// local alleles").
func (fv *FieldValue) empty() bool {
	if fv.Missing {
		return true
	}
	if len(fv.Raw) == 0 {
		return true
	}
	if len(fv.Raw) == 1 && fv.Raw[0] == "" {
		return true
	}
	return false
}

// hasPSandPSL reports whether the sample carries both
// PS and PSL with non-empty values.
func hasPSandPSL(g Genotype) bool {
	psFV, hasPS := g.Fields["PS"]
	if !hasPS || psFV.empty() {
		return false
	}
	pslFV, hasPSL := g.Fields["PSL"]
	return hasPSL && !pslFV.empty()
}

// ftBadTokens returns the FT tokens on this sample that
// contain whitespace or semicolons (VCFv4.5 §1.6.2 FT).
func ftBadTokens(g Genotype) []string {
	fv, ok := g.Fields["FT"]
	if !ok || fv.Missing {
		return nil
	}
	var bad []string
	for _, tok := range fv.Strs {
		if strings.ContainsAny(tok, " \t\n;") {
			bad = append(bad, tok)
		}
	}
	return bad
}

// localAlleleKeys lists the FORMAT keys that depend on
// LAA being present (VCFv4.5 §1.6.2). If any of these is
// declared on a sample without LAA, the sample is invalid.
var localAlleleKeys = [...]string{
	"LAD", "LADF", "LADR", "LEC",
	"LGL", "LGP", "LPL", "LPP",
}

// isLocalAlleleKey reports whether k is a local-allele
// FORMAT key. Matches the keys listed in
// localAlleleKeys; "LAA" itself is not considered a
// local-allele key here (it is the container).
func isLocalAlleleKey(k string) bool {
	for _, lk := range localAlleleKeys {
		if k == lk {
			return true
		}
	}
	return false
}

// missingLAAFields returns the local-allele FORMAT keys
// declared on the sample whose LAA field is missing or
// empty.
func missingLAAFields(g Genotype) []string {
	laaFV, hasLAA := g.Fields["LAA"]
	laaPresent := hasLAA && !laaFV.empty()
	if laaPresent {
		return nil
	}
	var missing []string
	for _, key := range localAlleleKeys {
		fv, ok := g.Fields[key]
		if !ok || fv.empty() {
			continue
		}
		missing = append(missing, key)
	}
	return missing
}

// laaLenFor returns the effective allele count for
// computing the cardinality of a local-allele FORMAT
// field. When LAA is present and non-empty, that's the
// count; otherwise fall back to the total ALT count
// (matching the spec's behaviour when LAA is missing).
//
// Per the spec LAA's declared Type is String (§1.6.2),
// so the values are stored in Strs (or fall back to Raw
// for non-numeric lengths). Use whatever length is
// non-zero so LAA works regardless of declared Type.
func laaLenFor(g Genotype, altCount int) int {
	laaFV, has := g.Fields["LAA"]
	if !has || laaFV.empty() {
		return altCount
	}
	if n := len(laaFV.Ints); n > 0 {
		return n
	}
	if n := len(laaFV.Strs); n > 0 {
		return n
	}
	return len(laaFV.Raw)
}

// wantedCardinality returns the expected element count
// for a FORMAT field. Local-allele FORMAT fields
// (LAD, LEC, LGL, …) derive their cardinality from LAA
// when present (§1.6.2); without LAA they fall back to
// the ALT count. The cross-check ensures LA-field
// values are indexed into LAA, not the full ALT list —
// the only meaningful interpretation of a per-allele
// (LAD/LADR) or per-genotype (LGL/LGP/LPL/LPP) value
// when LAA is present.
func wantedCardinality(n Number, key string, g Genotype, rec *Record) int {
	cardinalAlt := len(rec.Alt)
	if isLocalAlleleKey(key) {
		// For local-allele keys the cardinality is
		// computed against LAA, not REF+ALTs.
		// Number=G (= SymbolG) is the only case where
		// the conventional formula adds REF back in
		// (altCount+1); for LA semantics we want
		// genotypeCount(len(LAA), ploidy) instead, so
		// subtract 1 before calling Number.Cardinality.
		cardinalAlt = laaLenFor(g, cardinalAlt)
		if n.Fixed == SymbolG {
			cardinalAlt--
		}
		return n.Cardinality(cardinalAlt, ploidyFromGT(g))
	}
	return n.Cardinality(cardinalAlt, ploidyFromGT(g))
}

// samePSL reports whether two PSL slices match element
// by element. Two empty/nil slices match.
func samePSL(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// samePSO reports whether two PSO slices match element
// by element. Two empty/nil slices match.
func samePSO(a, b []int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// checkPSLUnphased enforces VCFv4.5 §1.6.2: "Unphased
// alleles (without a | separator before them) must have
// the value '$.$' in their corresponding position in the
// list". For each allele in GT separated by '/' (or
// appearing as the first allele), the corresponding PSL
// entry must be ".".
//
// Extracted from checkSamples to keep the per-sample
// loop under the cognitive-complexity lint cap.
func checkPSLUnphased(g Genotype, sIdx, idx, line int) []ValidationError {
	gtFV, hasGT := g.Fields["GT"]
	pslFV, hasPSL := g.Fields["PSL"]
	if !hasGT || !hasPSL {
		return nil
	}
	if gtFV.Missing || pslFV.Missing || len(gtFV.Strs) == 0 || len(pslFV.Strs) == 0 {
		return nil
	}
	if gtFV.Strs[0] == "" {
		return nil
	}
	gt := gtFV.Strs[0]
	var issues []ValidationError
	for i, c := range gt {
		if c == '/' || c == '|' || c == ' ' || c == '\t' {
			continue
		}
		// `c` is an allele digit. The separator that
		// preceded it determines phasing: '|' = phased
		// (no constraint on PSL); anything else (or
		// start of string) = unphased (PSL must be ".").
		phased := i > 0 && gt[i-1] == '|'
		if phased {
			continue
		}
		alleleIdx := alleleIndex(gt, i)
		if alleleIdx < len(pslFV.Strs) && pslFV.Strs[alleleIdx] != "." {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/PSL/sample[%d]", sIdx),
				Msg:   fmt.Sprintf("PSL[%d]=%q must be %q because the corresponding GT allele is unphased (§1.6.2)", alleleIdx, pslFV.Strs[alleleIdx], "."),
			})
		}
	}
	return issues
}

// alleleIndex returns the index of the allele at
// position `pos` in the GT string, ignoring
// separators. Position 0 (start of string) corresponds
// to allele index 0.
func alleleIndex(gt string, pos int) int {
	n := 0
	for i := 0; i < pos; i++ {
		c := gt[i]
		if c != '/' && c != '|' && c != ' ' && c != '\t' {
			n++
		}
	}
	return n
}

// checkGPRange enforces VCFv4.5 §1.6.2 GP field:
// "Genotype posterior probabilities in the range 0 to 1".
// Emits one issue per out-of-range value.
func checkGPRange(g Genotype, sIdx, idx, line int) []ValidationError {
	gp, ok := g.Fields["GP"]
	if !ok || gp.Missing {
		return nil
	}
	var issues []ValidationError
	for j, f := range gp.Floats {
		if j < len(gp.Raw) && gp.Raw[j] == "." {
			continue
		}
		if f < 0 || f > 1 {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/GP/sample[%d]", sIdx),
				Msg:   fmt.Sprintf("GP[%d]=%g is out of range; GP values must be between 0 and 1 (§1.6.2)", j, f),
			})
		}
	}
	return issues
}

// checkMRequiresGT enforces VCFv4.5 §1.6.1.5: "If any
// base modification key is present for a sample, GT must
// be defined for that sample." A present (non-missing)
// M-key for a sample whose GT is missing fires one
// SeverityError.
func checkMRequiresGT(g Genotype, sIdx, idx, line int) []ValidationError {
	gt, hasGT := g.Fields["GT"]
	if !hasGT || !gt.Missing {
		return nil
	}
	for k := range g.Fields {
		fv := g.Fields[k]
		if !IsBaseModKey(k) {
			continue
		}
		if fv.Missing {
			continue
		}
		return []ValidationError{{
			Severity: SeverityError, Record: idx + 1, Line: line,
			Field: fmt.Sprintf("FORMAT/%s/sample[%d]", k, sIdx),
			Msg:   fmt.Sprintf("base-modification field %q is present but GT is missing; §1.6.1.5 requires GT when any M-key is present", k),
		}}
	}
	return nil
}

// checkPSLTrio enforces VCFv4.5 §1.6.2 PSQ field: "The
// missing value '$.$' should be used when the corresponding
// PSL value is missing, or when the phasing is of unknown
// quality." We treat the rule as bidirectional: PSL[i]
// is missing iff PSO[i] is missing iff PSQ[i] is missing.
// Each mismatch fires one SeverityError. Spec uses "should";
// we still fire in ModeStrict (matching the existing
// PSL-unphased check).
func checkPSLTrio(g Genotype, sIdx, idx, line int) []ValidationError {
	pslFV, hasPSL := g.Fields["PSL"]
	if !hasPSL {
		return nil
	}
	if pslFV.Missing {
		return nil
	}
	pslCount := len(pslFV.Raw)
	var issues []ValidationError
	for _, k := range []string{"PSO", "PSQ"} {
		fv, ok := g.Fields[k]
		if !ok || fv.Missing {
			continue
		}
		for i := 0; i < pslCount && i < len(fv.Raw); i++ {
			pslMiss := pslFV.Raw[i] == "."
			otherMiss := fv.Raw[i] == "."
			if pslMiss == otherMiss {
				continue
			}
			var msg string
			if pslMiss {
				msg = fmt.Sprintf("%s[%d]=%q is non-missing but the corresponding PSL[%d] is missing (§1.6.2 PSQ)", k, i, fv.Raw[i], i)
			} else {
				msg = fmt.Sprintf("PSL[%d]=%q is non-missing but the corresponding %s[%d] is missing (§1.6.2 PSQ)", i, pslFV.Raw[i], k, i)
			}
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/%s/sample[%d]", k, sIdx),
				Msg:   msg,
			})
		}
	}
	return issues
}

// checkMKeyRange enforces VCFv4.5 §1.6.1.5: "Values
// must be between 0 and 1". Only applies to keys that
// match the M-key shape AND are declared as TypeFloat.
// Extracted from checkInfo to keep the per-key loop under
// the cognitive-complexity lint cap.
func checkMKeyRange(key string, val *InfoValue, idx, line int) []ValidationError {
	if !IsBaseModKey(key) {
		return nil
	}
	defType := val.Type
	if defType != TypeFloat {
		return nil
	}
	var issues []ValidationError
	for j, f := range val.Floats {
		if j < len(val.Raw) && val.Raw[j] == "." {
			continue
		}
		if f < 0 || f > 1 {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/" + key,
				Msg:   fmt.Sprintf("%s[%d]=%g is out of range; M-key values must be between 0 and 1 (§1.6.1.5)", key, j, f),
			})
		}
	}
	return issues
}

// parseBreakendPosition extracts the remote contig and
// 1-based position from a breakend replacement string of
// the forms t]p], t[p[, ]p]t, or [p[t. Per §1.6.5.1, the
// position specifier is bracketed by two matching brackets
// (`[` or `]`). The inner content is "chr[:pos[:strand]]"
// where the trailing bracket character (if present)
// encodes the strand of the remote end.
//
// Returns (contig, pos, true) on success. Returns
// ("", 0, false) when the string is not a breakend,
// when the position spec omits the position (a telomere
// breakpoint), or when the position cannot be parsed.
func parseBreakendPosition(a string) (contig string, pos int64, ok bool) {
	if len(a) < 4 {
		return "", 0, false
	}
	start := -1
	for i := 0; i < len(a); i++ {
		if a[i] == '[' || a[i] == ']' {
			start = i
			break
		}
	}
	if start < 0 {
		return "", 0, false
	}
	open := a[start]
	for i := start + 1; i < len(a); i++ {
		if a[i] != open {
			continue
		}
		inner := a[start+1 : i]
		colonIdx := strings.LastIndex(inner, ":")
		if colonIdx < 0 {
			return "", 0, false
		}
		contig := inner[:colonIdx]
		posStr := inner[colonIdx+1:]
		if posStr != "" && (posStr[len(posStr)-1] == '[' || posStr[len(posStr)-1] == ']') {
			posStr = posStr[:len(posStr)-1]
		}
		pos, err := strconv.ParseInt(posStr, 10, 64)
		if err != nil || pos < 1 {
			return "", 0, false
		}
		return contig, pos, true
	}
	return "", 0, false
}

// checkBreakendPositions fires when a FASTA reference is
// attached. For every breakend ALT, the remote position
// must be within the remote contig's length. ALTs that
// omit the position (a telomere breakpoint) are skipped
// because there is no position to bound-check. Unknown
// remote contigs degrade to SeverityWarning in ModeStrict
// and are skipped in ModeDefault.
func (v *Validator) checkBreakendPositions(rec *Record, idx, line int) []ValidationError {
	if v.Reference == nil {
		return nil
	}
	var issues []ValidationError
	for i, alt := range rec.Alt {
		if !isBreakend(alt) {
			continue
		}
		contig, pos, ok := parseBreakendPosition(alt)
		if !ok {
			continue
		}
		contigLen, err := v.Reference.Length(contig)
		if err != nil {
			iss := ValidationError{
				Severity: SeverityWarning,
				Record:   idx + 1, Line: line,
				Field: fmt.Sprintf("ALT[%d]", i),
				Msg:   fmt.Sprintf("breakend references unknown contig %q: %v", contig, err),
			}
			if v.Mode == ModeStrict {
				iss.Severity = SeverityError
			}
			issues = append(issues, iss)
			continue
		}
		if pos > contigLen {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("ALT[%d]", i),
				Msg:   fmt.Sprintf("breakend position %d on %q exceeds contig length %d", pos, contig, contigLen),
			})
		}
	}
	return issues
}

// checkMateMatch compares two breakend mate snapshots
// and appends SeverityErrors on mismatch. Two checks
// run:
//
//  1. PSL and PSO (per-allele) must be equal (§1.6.2).
//  2. Each mate's breakend remote positions must
//     point at the partner's local (Chrom, POS) —
//     §1.6.5.1 breakend round-trip.
func checkMateMatch(existing, snap *mateSnapshot, mid string, idx, line int, issues *[]ValidationError) {
	if !samePSL(existing.psl, snap.psl) {
		*issues = append(*issues, ValidationError{
			Severity: SeverityError, Record: idx + 1, Line: line,
			Field: "INFO/PSL",
			Msg:   fmt.Sprintf("PSL for MATEID %q does not match the mate record", mid),
		})
	}
	if !samePSO(existing.pso, snap.pso) {
		*issues = append(*issues, ValidationError{
			Severity: SeverityError, Record: idx + 1, Line: line,
			Field: "INFO/PSO",
			Msg:   fmt.Sprintf("PSO for MATEID %q does not match the mate record", mid),
		})
	}
	// §1.6.5.1 breakend round-trip: this mate's
	// breakend remote positions must equal the
	// partner's local (Chrom, POS). Only fires when
	// both mates carry breakend ALTs; otherwise one
	// side might be a paired SNV (MATEID linking two
	// distinct records), which has no breakend
	// position to round-trip.
	if len(existing.breakends) > 0 && len(snap.breakends) > 0 {
		for _, bp := range snap.breakends {
			if bp.contig != existing.chrom || bp.pos != existing.pos {
				*issues = append(*issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "ALT",
					Msg: fmt.Sprintf("breakend mate %q points to %s:%d but partner is at %s:%d (§1.6.5.1)",
						mid, bp.contig, bp.pos, existing.chrom, existing.pos),
				})
			}
		}
		for _, bp := range existing.breakends {
			if bp.contig != snap.chrom || bp.pos != snap.pos {
				*issues = append(*issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "ALT",
					Msg: fmt.Sprintf("breakend mate %q (already registered) points to %s:%d but partner is at %s:%d (§1.6.5.1)",
						mid, bp.contig, bp.pos, snap.chrom, snap.pos),
				})
			}
		}
	}
}

// checkCNVAllelesShareSVLEN enforces the VCFv4.5
// §"FORMAT keys used for structural variants" rule:
//
//	When FORMAT CN is present, all <CNV>, <DEL> and
//	<DUP> alleles must have the same SVLEN.
//
// The check fires in ModeStrict.
func checkCNVAllelesShareSVLEN(rec *Record, svlens []int64, idx, line int, issues *[]ValidationError) {
	// Collect SVLENs from <CNV>, <DEL>, <DUP> alleles,
	// skipping missing values (.). The first non-missing
	// SVLEN is the reference; subsequent ones must match.
	var first int64
	firstSeen := false
	for i, alt := range rec.Alt {
		if !isCNVAllele(alt) {
			continue
		}
		if i >= len(svlens) || svlens[i] == -1 || isMissing(svlens[i]) {
			continue
		}
		if !firstSeen {
			first = svlens[i]
			firstSeen = true
			continue
		}
		if svlens[i] != first {
			*issues = append(*issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/SVLEN",
				Msg:   fmt.Sprintf("all <CNV>, <DEL>, and <DUP> alleles must share the same SVLEN when FORMAT CN is declared; %q has SVLEN=%d, expected %d", alt, svlens[i], first),
			})
		}
	}
}

// isCNVAllele reports whether a is a CNV-family symbolic
// allele (<CNV>, <CNV:TR>, <DEL>, <DUP>, including their
// subtypes).
func isCNVAllele(a string) bool {
	if !strings.HasPrefix(a, "<") || !strings.HasSuffix(a, ">") {
		return false
	}
	inner := a[1 : len(a)-1]
	for _, prefix := range []string{"CNV", "DEL", "DUP"} {
		if inner == prefix || strings.HasPrefix(inner, prefix+":") {
			return true
		}
	}
	return false
}

// isCNVTRAllelesPresent reports whether the record
// carries any <CNV:TR> allele.
func isCNVTRAllelesPresent(rec *Record) bool {
	for _, a := range rec.Alt {
		if strings.HasPrefix(a, "<CNV:TR") && strings.HasSuffix(a, ">") {
			return true
		}
	}
	return false
}

// hasAnyInfoField reports whether the INFO map carries
// any of the named keys with a present value (Flag or
// non-missing).
func hasAnyInfoField(info map[string]InfoValue, keys ...string) bool {
	for _, k := range keys {
		if fv, ok := info[k]; ok && fv.HasValue() {
			return true
		}
	}
	return false
}

// isMissing is a small helper that treats a sentinel
// value (the parser uses -1 to mean "missing" for the
// Ints slice when the field has fewer entries than ALT
// alleles) as missing. The actual VALUE MISSING (".")
// round-trips through as 0 in the typedCore, but for
// SVLEN a 0 entry is already caught by the
// "must be positive" check; we skip the SVLEN-
// equality check on zero values to avoid false
// positives.
func isMissing(n int64) bool { return n <= 0 }

// impreciseFlagged reports whether the record carries
// the ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,...>
// declaration and the IMPRECISE flag is set.
func impreciseFlagged(rec *Record) bool {
	iv, ok := rec.Info["IMPRECISE"]
	if !ok {
		return false
	}
	return iv.Flag
}

// checkIntegerBCFRange enforces VCFv4.5 §1.4.2: Integer
// values from -2147483648 to -2147483641 (the eight
// values just above math.MinInt32) cannot be stored in
// BCF2 and are therefore disallowed in both VCF and
// BCF (§"BcfTypeEncoding"). The check walks every
// Integer-typed INFO and FORMAT field across the
// record and emits a SeverityError per offending
// value.
func checkIntegerBCFRange(rec *Record, idx, line int) []ValidationError {
	const bcfMin = int64(math.MinInt32)
	const bcfMaxForbidden = bcfMin + 7 // -2147483641
	var issues []ValidationError
	checkField := func(field string, ints []int64) {
		for i, v := range ints {
			if v < bcfMin || v > bcfMaxForbidden {
				continue
			}
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("%s[%d]", field, i),
				Msg: fmt.Sprintf("Integer value %d is in the BCF2 forbidden range [%d, %d] (§1.4.2 / §BcfTypeEncoding)",
					v, bcfMin, bcfMaxForbidden),
			})
		}
	}
	for key := range rec.Info {
		val := rec.Info[key]
		if len(val.Ints) == 0 {
			continue
		}
		checkField("INFO/"+key, val.Ints)
	}
	for sIdx, sample := range rec.Samples {
		for key := range sample.Fields {
			val := sample.Fields[key]
			if len(val.Ints) == 0 {
				continue
			}
			checkField(fmt.Sprintf("FORMAT/%s/sample[%d]", key, sIdx), val.Ints)
		}
	}
	return issues
}

// ciposCovers reports whether the CIPOS INFO field
// carries an entry for the i-th ALT allele. The
// IMPRECISE rule cares only that the field has the
// expected slot — [0,0] is a valid explicit value per
// the spec ("CIPOS must be explicitly set for that
// allele, even if it is 0,0"). What is missing is
// the absence of the pair (i.e. the field has fewer
// than i*2+2 entries, or the field is missing).
func ciposCovers(rec *Record, i int) bool {
	iv, ok := rec.Info["CIPOS"]
	if !ok || iv.Missing {
		return false
	}
	return i*2+1 < len(iv.Ints)
}

// isSymbolicSV reports whether a is a symbolic
// structural-variant ALT allele (not <*>, not a
// breakend, not a bare IUPAC code).
func isSymbolicSV(a string) bool {
	if !strings.HasPrefix(a, "<") || !strings.HasSuffix(a, ">") {
		return false
	}
	if a == "<*>" || strings.ContainsAny(a, "[]") {
		return false
	}
	return isCNVAllele(a)
}

// svlenPresent reports whether the SVLEN INFO field
// carries a present (non-missing) entry for the i-th
// ALT allele.
func svlenPresent(rec *Record, i int) bool {
	iv, ok := rec.Info["SVLEN"]
	if !ok || iv.Missing {
		return false
	}
	if i >= len(iv.Ints) {
		return false
	}
	v := iv.Ints[i]
	return !isMissing(v)
}

// checkCIGAR enforces VCFv4.5 §1.6.1.1: the sum of REF-
// consuming op lengths in a CIGAR string must equal
// len(REF), and the sum of ALT-consuming op lengths
// must equal len(ALT) for the corresponding ALT.
// Extracted from checkInfo to keep the per-record
// check under the cognitive-complexity lint cap.
func checkCIGAR(rec *Record, idx, line int) []ValidationError {
	cigar, ok := rec.Info["CIGAR"]
	if !ok || cigar.Missing || cigar.Type != TypeString {
		return nil
	}
	var issues []ValidationError
	for i, c := range cigar.Strs {
		refLen, altLen, perr := parseCIGAR(c)
		if perr != nil {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/CIGAR",
				Msg:   fmt.Sprintf("CIGAR[%d]=%q is malformed: %v", i, c, perr),
			})
			continue
		}
		if refLen != len(rec.Ref) {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/CIGAR",
				Msg:   fmt.Sprintf("CIGAR[%d] consumes %d reference bases but REF has length %d", i, refLen, len(rec.Ref)),
			})
		}
		if i < len(rec.Alt) {
			want := len(rec.Alt[i])
			if altLen != want {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/CIGAR",
					Msg:   fmt.Sprintf("CIGAR[%d] consumes %d query bases but ALT[%d] has length %d", i, altLen, i, want),
				})
			}
		}
	}
	return issues
}

// parseCIGAR parses a SAM/BAM-style CIGAR string and
// returns the number of reference bases and query
// (ALT) bases it consumes, plus an error if the string
// is malformed. SAM CIGAR ops:
//
//	M, =, X: consume 1 ref + 1 query
//	I, S, H: consume 0 ref + N query
//	D, N    : consume N ref + 0 query
//	P       : consume 0 ref + 0 query (padding)
func parseCIGAR(cigar string) (refBases, altBases int, err error) {
	i := 0
	for i < len(cigar) {
		// Read the run length (decimal digits).
		n := 0
		hadDigit := false
		for i < len(cigar) && cigar[i] >= '0' && cigar[i] <= '9' {
			hadDigit = true
			n = n*10 + int(cigar[i]-'0')
			if n > 1<<24 {
				return 0, 0, fmt.Errorf("run length overflow")
			}
			i++
		}
		if !hadDigit {
			return 0, 0, fmt.Errorf("expected digit at offset %d", i)
		}
		if i >= len(cigar) {
			return 0, 0, fmt.Errorf("missing op at end of CIGAR")
		}
		op := cigar[i]
		i++
		switch op {
		case 'M', '=', 'X':
			refBases += n
			altBases += n
		case 'I', 'S', 'H':
			altBases += n
		case 'D', 'N':
			refBases += n
		case 'P':
			// Padding consumes nothing.
		default:
			return 0, 0, fmt.Errorf("unknown CIGAR op %q at offset %d", op, i-1)
		}
	}
	return refBases, altBases, nil
}

// checkRefMatchesReference compares every byte of REF
// against the corresponding reference base at
// POS+i. Returns one SeverityError per mismatch, or one
// SeverityWarning if the reference cannot answer
// (unknown contig, position out of range). The warning
// is dropped in ModeDefault.
//
// This is the v0.11.0 successor to the older single-base
// padding check (checkPaddingBase, removed). It strictly
// subsumes that check while also catching wrong contig,
// wrong POS, transcription errors, and undeclared IUPAC
// reductions.
func (v *Validator) checkRefMatchesReference(rec *Record, idx, line int) []ValidationError {
	pos := int(rec.Pos)
	if pos < 1 {
		return nil
	}
	end := pos + len(rec.Ref)
	bases, err := v.Reference.BaseRange(rec.Chrom, pos, end)
	if err != nil {
		iss := ValidationError{
			Severity: SeverityWarning,
			Record:   idx + 1, Line: line,
			Field: "REF",
			Msg:   fmt.Sprintf("reference lookup for %s:%d-%d failed: %v; skipping full-REF check", rec.Chrom, pos, end-1, err),
		}
		if v.Mode == ModeStrict {
			iss.Severity = SeverityError
		}
		return []ValidationError{iss}
	}
	var issues []ValidationError
	for i, refBase := range bases {
		want := strings.ToUpper(string(rec.Ref[i]))
		got := strings.ToUpper(string(refBase))
		if want != got {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "REF",
				Msg:   fmt.Sprintf("REF[%d]=%q does not match reference base %q at %s:%d (§1.6.1.1)", i, string(rec.Ref[i]), string(refBase), rec.Chrom, pos+i),
			})
		}
	}
	return issues
}

// checkLAAPosition enforces VCFv4.5 §1.6.2: "LAA must
// also be present and precede all fields other than GT".
// Extracted from checkSamples to keep the per-sample
// loop under the cognitive-complexity lint cap.
func checkLAAPosition(rec *Record, sIdx, idx, line int) []ValidationError {
	laaIdx := -1
	for i, k := range rec.Format {
		if k == "LAA" {
			laaIdx = i
			break
		}
	}
	if laaIdx < 0 {
		return nil
	}
	var issues []ValidationError
	for i, k := range rec.Format {
		if i == laaIdx || k == "GT" {
			continue
		}
		if !isLocalAlleleKey(k) {
			continue
		}
		if i < laaIdx {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/%s/sample[%d]", k, sIdx),
				Msg:   fmt.Sprintf("local-allele field %q must appear after LAA in the FORMAT column", k),
			})
		}
	}
	return issues
}

// checkLAA validates a sample's LAA field. Per VCFv4.5
// §1.6.1.4, LAA is a list of 1-based indices into ALT and
// the values must be distinct.
func checkLAA(vv *FieldValue, rec *Record, idx, sIdx, line int) []ValidationError {
	if len(vv.Ints) == 0 {
		return nil
	}
	nAlts := int64(len(rec.Alt))
	var issues []ValidationError
	seen := make(map[int64]int, len(vv.Ints))
	for i, n := range vv.Ints {
		if n < 1 || n > nAlts {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/LAA/sample[%d]", sIdx),
				Msg: fmt.Sprintf("LAA[%d]=%d is out of range (must be a 1-based ALT index in [1, %d])",
					i, n, nAlts),
			})
		}
		if prev, dup := seen[n]; dup {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: fmt.Sprintf("FORMAT/LAA/sample[%d]", sIdx),
				Msg: fmt.Sprintf("LAA[%d]=%d duplicates LAA[%d] (values must be distinct)",
					i, n, prev),
			})
		} else {
			seen[n] = i
		}
	}
	return issues
}

// ploidyFromGT derives the per-sample ploidy from a GT string by counting
// allele indices separated by / or |. Returns 0 when GT is missing.
func ploidyFromGT(g Genotype) int {
	v, ok := g.Fields["GT"]
	if !ok || v.Missing || len(v.Strs) == 0 {
		return 0
	}
	gt := v.Strs[0]
	if gt == "" {
		return 0
	}
	count := 0
	for _, c := range gt {
		if c == '/' || c == '|' {
			count++
		}
	}
	return count + 1
}

// checkGT performs the two structural checks the spec
// requires for a GT field:
//  1. A single record may not mix / and | separators (the
//     phasing must be consistent within a record).
//  2. Each allele index must be a valid 0-based index into
//     REF+ALT. The missing value "." is always allowed.
//
// Returns nil for a clean GT. nAlts is the number of ALT
// alleles in the record (0 means REF-only).
func checkGT(gt string, nAlts int) error {
	if gt == "" {
		return nil
	}
	var sawPipe, sawSlash bool
	for i := 0; i < len(gt); {
		c := gt[i]
		switch c {
		case '|':
			sawPipe = true
			i++
		case '/':
			sawSlash = true
			i++
		default:
			// Read the allele token up to the next separator.
			j := i
			for j < len(gt) && gt[j] != '|' && gt[j] != '/' {
				j++
			}
			tok := gt[i:j]
			if tok != "." {
				n, err := strconv.Atoi(tok)
				if err != nil {
					return fmt.Errorf("GT allele %q is not an integer", tok)
				}
				if n < 0 || n > nAlts {
					return fmt.Errorf("GT allele %d is out of range (REF + %d ALT alleles)", n, nAlts)
				}
			}
			i = j
		}
	}
	if sawPipe && sawSlash {
		return fmt.Errorf("GT mixes unphased (/) and phased (|) separators within a single record")
	}
	return nil
}
