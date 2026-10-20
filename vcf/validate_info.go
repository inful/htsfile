package vcf

import (
	"fmt"
	"strings"
)

// (Validator lifecycle) is in validate.go.
func (v *Validator) checkInfo(h *Header, rec *Record, idx int) []ValidationError {
	var issues []ValidationError
	altCount := len(rec.Alt)
	line := rec.Line

	// BCF Integer range (§1.4.2): the values from
	// $-2^{31}$ to $-2^{31}+7$ cannot be stored in
	// BCF2 and are therefore disallowed in both VCF
	// and BCF (§"BcfTypeEncoding"). The eight
	// forbidden values are MIN_INT32..MIN_INT32+7.
	issues = append(issues, checkIntegerBCFRange(rec, idx, line)...)

	// CIGAR consistency with REF and ALT (§1.6.1.1).
	issues = append(issues, checkCIGAR(rec, idx, line)...)

	// IMPRECISE + precise ALT requires CIPOS (§1.6.1.7):
	// "If a precise ALT allele is present in a record
	// with the IMPRECISE flag, CIPOS must be explicitly
	// set for that allele, even if it is 0,0."
	if impreciseFlagged(rec) {
		for i, alt := range rec.Alt {
			if isBreakend(alt) {
				continue
			}
			if !ciposCovers(rec, i) {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/CIPOS",
					Msg:   fmt.Sprintf("IMPRECISE flag is set but CIPOS is missing for precise ALT %q; §1.6.1.7 requires CIPOS for every precise allele when IMPRECISE is set", alt),
				})
			}
		}
	}

	// SVLEN required for symbolic structural variant
	// alleles (§1.6.1.7): "SVLEN must be specified for
	// symbolic structural variant alleles." A missing
	// SVLEN entry is -1 in the typed-core Ints slice.
	for i, alt := range rec.Alt {
		if !isSymbolicSV(alt) {
			continue
		}
		if !svlenPresent(rec, i) {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/SVLEN",
				Msg:   fmt.Sprintf("symbolic ALT %q requires SVLEN to be specified (§1.6.1.7)", alt),
			})
		}
	}

	for key := range rec.Info {
		val := rec.Info[key]
		def, declared := h.Infos[key]
		if !declared {
			iss := ValidationError{
				Severity: SeverityWarning,
				Record:   idx + 1,
				Line:     line,
				Field:    "INFO/" + key,
				Msg:      "INFO key is not declared in the header",
			}
			if v.Mode == ModeStrict {
				iss.Severity = SeverityError
			}
			issues = append(issues, iss)
			continue
		}

		// Flag fields must not carry values; conversely, non-flag fields
		// must carry at least one value.
		if def.Type == TypeFlag {
			if !val.Flag {
				// Defensive: the parser already sets Flag=true on Flag
				// entries, but downstream tools might construct a Record
				// programmatically.
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/" + key,
					Msg:   "Flag field not marked present",
				})
			}
			continue
		}

		// Number checking.
		want := def.Number.Cardinality(altCount, 0)
		got := len(val.Raw)
		if val.Missing {
			got = 0
		}
		if want >= 0 && got != want {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/" + key,
				Msg: fmt.Sprintf("expected %d value(s) for Number=%s, got %d",
					want, def.Number, got),
			})
		}

		// Type checking.
		for _, r := range val.Raw {
			if !validTypedValue(def.Type, r) {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/" + key,
					Msg:   fmt.Sprintf("value %q does not match declared Type=%q", r, def.Type),
				})
			}
		}

		// Integer range: the BCF binary format cannot store
		// the 8 values in the range [-2^31, -2^31+7]. The
		// spec forbids these in both VCF and BCF.
		if def.Type == TypeInteger {
			for _, n := range val.Ints {
				if n >= -2147483648 && n <= -2147483641 {
					issues = append(issues, ValidationError{
						Severity: SeverityError, Record: idx + 1, Line: line,
						Field: "INFO/" + key,
						Msg:   fmt.Sprintf("Integer value %d is in the BCF-unsafe range [-2^31, -2^31+7]", n),
					})
				}
			}
		}

		// M-key value range (§1.6.1.5).
		issues = append(issues, checkMKeyRange(key, &val, idx, line)...)

		// SVLEN sign rule: for symbolic SV alleles the length
		// must be positive. (The spec says SVLEN is always
		// positive, even for INV.)
		if key == "SVLEN" {
			// When FORMAT CN is present, all <CNV>,
			// <DEL>, and <DUP> alleles must share the
			// same SVLEN (VCFv4.5 §1.6.1.7 SV FORMAT
			// keys).
			if _, hasCN := h.Formats["CN"]; hasCN {
				checkCNVAllelesShareSVLEN(rec, val.Ints, idx, line, &issues)
			}
			for i, n := range val.Ints {
				if n <= 0 {
					alt := "."
					if i < len(rec.Alt) {
						alt = rec.Alt[i]
					}
					issues = append(issues, ValidationError{
						Severity: SeverityError, Record: idx + 1, Line: line,
						Field: "INFO/SVLEN",
						Msg:   fmt.Sprintf("SVLEN[%d] for ALT %q must be positive (got %d)", i, alt, n),
					})
				}
			}
		}

		// CIPOS / CIEND / CILEN confidence intervals: the
		// interval must span 0. The field is a list of pairs
		// (low, high) per ALT allele. The lower bound cannot
		// be greater than 0, and the upper bound cannot be
		// less than 0. The field must also contain exactly
		// 2 * len(ALT) entries (VCFv4.5 §1.6.1.7) — each
		// ALT allele gets one pair.
		if key == "CIPOS" || key == "CIEND" || key == "CILEN" {
			// Cardinality: 2 * len(ALT). A trailing odd
			// entry or any count other than 2*N is
			// rejected. A field of length 0 means the
			// field is absent (the spec says "If
			// present, the number of entries must be
			// twice the number of ALT alleles"); we skip
			// the cardinality check when the field is
			// effectively empty.
			if len(val.Ints) > 0 && len(val.Ints) != 2*altCount {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/" + key,
					Msg: fmt.Sprintf("%s must contain exactly %d entries (2 * %d ALT alleles); got %d",
						key, 2*altCount, altCount, len(val.Ints)),
				})
			}
			for i := 0; i+1 < len(val.Ints); i += 2 {
				lo, hi := val.Ints[i], val.Ints[i+1]
				if lo > 0 || hi < 0 {
					issues = append(issues, ValidationError{
						Severity: SeverityError, Record: idx + 1, Line: line,
						Field: "INFO/" + key,
						Msg: fmt.Sprintf("%s pair [%d, %d] must span 0 (lower ≤ 0, upper ≥ 0)",
							key, lo, hi),
					})
				}
			}
		}

		// MEINFO / METRANS: 4 entries per ALT allele
		// (VCFv4.5 §1.6.1.7).
		if key == "MEINFO" || key == "METRANS" {
			if len(val.Strs) > 0 && len(val.Strs) != 4*altCount {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/" + key,
					Msg: fmt.Sprintf("%s must contain exactly %d entries (4 * %d ALT alleles); got %d",
						key, 4*altCount, altCount, len(val.Strs)),
				})
			}
		}

		// CICN: 2 entries per ALT allele (VCFv4.5
		// §1.6.1.7).
		if key == "CICN" {
			if len(val.Floats) > 0 && len(val.Floats) != 2*altCount {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/" + key,
					Msg: fmt.Sprintf("CICN must contain exactly %d entries (2 * %d ALT alleles); got %d",
						2*altCount, altCount, len(val.Floats)),
				})
			}
		}
	}

	// <CNV:TR> alleles must declare either RUS or RUL
	// (VCFv4.5 §"Representing tandem repeats": "RUS or
	// RUL must be specified for each <CNV:TR>"). The
	// check fires per-record after the per-key loop so
	// INFO="." records (which never enter the loop)
	// still trigger the rule.
	if isCNVTRAllelesPresent(rec) && !hasAnyInfoField(rec.Info, "RUS", "RUL") {
		issues = append(issues, ValidationError{
			Severity: SeverityError, Record: idx + 1, Line: line,
			Field: "INFO",
			Msg:   "<CNV:TR> allele is present but neither RUS nor RUL is declared",
		})
	}

	// RN-determined fields (RUS / RUL / RUC / RB / CIRUC /
	// CIRB) must have a total length equal to the sum of
	// RN values (VCFv4.5 §1.6.1.7 RN). Missing RN entries
	// are treated as 0.
	if rnSum := computeRNSum(rec); rnSum >= 0 {
		checkRNDeterminedFields(rec, rnSum, idx, line, &issues)
	}
	return issues
}

// checkRNDeterminedFields enforces the cardinality
// invariants on RUS / RUL / RUC / RB / CIRUC / CIRB
// relative to the sum of RN values (VCFv4.5 §1.6.1.7)
// and the RUS length vs RUL value rule. Extracted to
// keep checkInfo under the cognitive-complexity lint cap.
func checkRNDeterminedFields(rec *Record, rnSum, idx, line int, issues *[]ValidationError) {
	for _, key := range rnDeterminedKeys {
		fv, has := rec.Info[key]
		if !has || !fv.HasValue() {
			continue
		}
		got := valueLen(&fv)
		if got != rnSum {
			*issues = append(*issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/" + key,
				Msg:   fmt.Sprintf("%s length %d does not match sum of RN values %d", key, got, rnSum),
			})
		}
	}
	// CIRUC/CIRB are paired confidence intervals around
	// RUC/RB; their total length is 2 × sum(RN).
	for _, key := range rnDeterminedKeysDoubled {
		fv, has := rec.Info[key]
		if !has || !fv.HasValue() {
			continue
		}
		got := valueLen(&fv)
		want := 2 * rnSum
		if got != want {
			*issues = append(*issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/" + key,
				Msg:   fmt.Sprintf("%s length %d does not match 2 × sum of RN values (%d)", key, got, want),
			})
		}
	}
	// CIRUC/CIRB sub-allele missingness must match RUC/RB
	// (VCFv4.5 §1.6.1.7: "CIRUC/CIRB must not be non-missing
	// for any alleles with no corresponding RUC/RB value").
	// RUC and RB are integer vectors; their raw "." entries
	// mark the missing values that the parser collapses to
	// 0 in Ints. CIRUC/CIRB are length 2*RN, paired against
	// RUC/RB by position: CIRUC[2i] and CIRUC[2i+1] are the
	// confidence bounds around RUC[i].
	if ruc, ok := rec.Info["RUC"]; ok && ruc.HasValue() {
		checkCIRUCMissing(&ruc, rec, "CIRUC", idx, line, issues)
	}
	if rb, ok := rec.Info["RB"]; ok && rb.HasValue() {
		checkCIRUCMissing(&rb, rec, "CIRB", idx, line, issues)
	}
	// CIRUC/CIRB bounds span zero
	// (VCFv4.5 §1.6.1.7: "The lower bound must be less
	// than or equal to zero and the upper bound must be
	// greater than or equal to zero."). Missing values
	// for the lower bound are treated as 0; for the
	// upper bound they're treated as unbounded. The
	// check runs per (CIRUC, CIRB) key.
	for _, key := range rnDeterminedKeysDoubled {
		fv, has := rec.Info[key]
		if !has || !fv.HasValue() {
			continue
		}
		// CIRUC is Float, CIRB is Integer. Walk both
		// representations and dispatch by the value
		// slice that's actually populated.
		var loVals, hiVals []int64
		switch key {
		case "CIRUC":
			if len(fv.Floats) >= 2 {
				loVals = make([]int64, 0, len(fv.Floats)/2)
				hiVals = make([]int64, 0, len(fv.Floats)/2)
				for i := 0; i+1 < len(fv.Floats); i += 2 {
					loVals = append(loVals, int64(fv.Floats[i]))
					hiVals = append(hiVals, int64(fv.Floats[i+1]))
				}
			}
		case "CIRB":
			loVals = fv.Ints
			hiVals = fv.Ints
		}
		for i := 0; i < len(loVals); i++ {
			lo, hi := loVals[i], hiVals[i]
			if lo > 0 {
				*issues = append(*issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/" + key,
					Msg:   fmt.Sprintf("%s lower bound must be ≤ 0 (§1.6.1.7); got %d at pair %d", key, lo, i),
				})
			}
			if hi < 0 {
				*issues = append(*issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/" + key,
					Msg:   fmt.Sprintf("%s upper bound must be ≥ 0 (§1.6.1.7); got %d at pair %d", key, hi, i),
				})
			}
		}
	}

	// RUS length must match RUL per element (VCFv4.5
	// §1.6.1.7 RUS).
	rus, hasRUS := rec.Info["RUS"]
	rul, hasRUL := rec.Info["RUL"]
	if !hasRUS || !hasRUL || !rus.HasValue() || !rul.HasValue() {
		return
	}
	if len(rus.Strs) != len(rul.Ints) {
		return
	}
	for i, s := range rus.Strs {
		rlen := len(s)
		if rlen != int(rul.Ints[i]) {
			*issues = append(*issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/RUS",
				Msg:   fmt.Sprintf("RUS[%d] length %d does not match RUL[%d]=%d", i, rlen, i, rul.Ints[i]),
			})
		}
	}
}

// pairKey returns the base-key name for a confidence key:
// "CIRUC" -> "RUC", "CIRB" -> "RB".
func pairKey(confKey string) string {
	switch confKey {
	case "CIRUC":
		return "RUC"
	case "CIRB":
		return "RB"
	}
	return ""
}

// checkCIRUCMissing emits one SeverityError per non-missing
// CIRUC/CIRB pair whose corresponding RUC/RB entry is the
// missing sentinel. The pairing is: confidence key entry j
// corresponds to confidence key entries [2j, 2j+1] and base
// key entry j (so a missing base entry forces the pair to
// be missing too).
func checkCIRUCMissing(base *InfoValue, rec *Record, confKey string, idx, line int, issues *[]ValidationError) {
	cf, ok := rec.Info[confKey]
	if !ok || !cf.HasValue() {
		return
	}
	for j := 0; j < len(base.Raw); j++ {
		if base.Raw[j] != "." {
			continue
		}
		for k := 0; k < 2; k++ {
			cj := 2*j + k
			if cj >= len(cf.Raw) {
				break
			}
			if cf.Raw[cj] != "." {
				*issues = append(*issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "INFO/" + confKey,
					Msg:   fmt.Sprintf("%s[%d] is non-missing but the corresponding %s[%d] is missing (§1.6.1.7)", confKey, cj, pairKey(confKey), j),
				})
			}
		}
	}
}

// computeRNSum returns the sum of all RN values in the
// record (treating missing as 0 for non-<CNV:TR>
// alleles and 1 for <CNV:TR> alleles when RN is absent).
// Returns -1 when no <CNV:TR> allele is present so the
// caller can skip the RN-determined checks.
func computeRNSum(rec *Record) int {
	hasCNVTR := false
	for _, a := range rec.Alt {
		if strings.HasPrefix(a, "<CNV:TR") && strings.HasSuffix(a, ">") {
			hasCNVTR = true
			break
		}
	}
	if !hasCNVTR {
		return -1
	}
	rnFV, hasRN := rec.Info["RN"]
	if !hasRN {
		// Default RN per allele (1 for CNV:TR, 0
		// otherwise) — but we only care about CNV:TR
		// alleles in this branch.
		var sum int
		for _, a := range rec.Alt {
			if strings.HasPrefix(a, "<CNV:TR") && strings.HasSuffix(a, ">") {
				sum++
			}
		}
		return sum
	}
	var sum int
	for i, n := range rnFV.Ints {
		if i < len(rec.Alt) && strings.HasPrefix(rec.Alt[i], "<CNV:TR") {
			if n > 0 {
				sum += int(n)
			}
		}
	}
	return sum
}

// valueLen returns the cardinality of a typedCore for
// the purpose of RN-summed checks: number of entries
// for strings/floats, number of integer values for ints.
func valueLen(v *InfoValue) int {
	switch v.Type {
	case TypeInteger:
		return len(v.Ints)
	case TypeFloat:
		return len(v.Floats)
	case TypeString, TypeCharacter:
		return len(v.Strs)
	case TypeFlag:
		if v.Flag {
			return 1
		}
	}
	return len(v.Raw)
}

// rnDeterminedKeys lists the INFO keys whose total length
// must equal the sum of RN values for the record
// (VCFv4.5 §1.6.1.7 RN).
// rnDeterminedKeys lists the INFO keys whose total length
// must equal the sum of RN values for the record
// (VCFv4.5 §1.6.1.7 RN).
var rnDeterminedKeys = [...]string{
	"RUS", "RUL", "RUC", "RB",
}

// rnDeterminedKeysDoubled lists the INFO keys whose total
// length must equal twice the sum of RN values
// (VCFv4.5 §1.6.1.7: "The length of these fields must be
// twice that of their corresponding fields" — applies to
// the confidence-interval fields CIRUC and CIRB).
var rnDeterminedKeysDoubled = [...]string{
	"CIRUC", "CIRB",
}

// reservedEventTypes is the set of EVENTTYPE values the spec
// reserves (VCFv4.5 §1.6.1.7 "The following EVENTTYPE values
// are reserved and should be used when appropriate"). The
// spec also says "The semantics of other EVENTTYPE values
// is implementation-defined" — the validator warns in
// ModeStrict for non-reserved values, but they are not
// strictly invalid.
var reservedEventTypes = map[string]bool{
	"DEL":            true,
	"DEL:ME":         true,
	"INS":            true,
	"INS:ME":         true,
	"DUP":            true,
	"DUP:TANDEM":     true,
	"DUP:DISPERSED":  true,
	"INV":            true,
	"TRA":            true,
	"TRA:BALANCED":   true,
	"TRA:UNBALANCED": true,
	"CHROMOTHRIPSIS": true,
	"CHROMOPLEXY":    true,
	"BFB":            true,
	"DOUBLEMINUTE":   true,
}

// IsReservedEventType reports whether value is one of the
// spec-reserved EVENTTYPE values.
func IsReservedEventType(value string) bool {
	return reservedEventTypes[value]
}

// checkEventType warns when a record's EVENTTYPE value is
// not in the reserved list. Implementation-defined values
// are allowed (not an error) but flagged in ModeStrict so
// callers can spot typos. The check is a no-op when the
// record has no EVENTTYPE key.
