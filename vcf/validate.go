package vcf

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"
)

// Mode controls the strictness of the [Validator].
//
//	ModeDefault — the parsed record values must agree with the header-declared
//	              Number/Type, declared keys must exist, declared FILTER
//	              codes must be valid. This is the recommended starting
//	              point and is what most users want.
//	ModeStrict   — additionally enforces spec-required rules like
//	              sorted-by-POS-within-CHROM, REF alphabet, ALT not equal
//	              to REF, and ##fileformat-first. May reject some real-world
//	              VCFs that the spec marks "recommended" but does not strictly
//	              require.
type Mode int

const (
	// ModeDefault is the well-formed + semantically-consistent mode.
	ModeDefault Mode = iota
	// ModeStrict adds spec-compliance rules.
	ModeStrict
)

// Validator runs semantic checks against a parsed VCF.
//
// Typical usage:
//
//	v := &vcf.Validator{Mode: vcf.ModeDefault}
//	issues := v.Validate(hdr, records)
//	for _, iss := range issues {
//	    log.Printf("%v\n", iss)
//	}
//
// A Validator is safe to re-use across multiple Validate calls. The
// streaming API ([Validator.HeaderChecks], [Validator.Next],
// [Validator.Finish]) holds cross-record state in the
// Validator itself, so do not share a Validator across
// concurrent streams.
type Validator struct {
	Mode Mode

	// Streaming state (used by Next / Finish). Zero value is
	// fine; the validator initialises lazily.
	recordIdx   int     // 0-based index of the most recently fed record
	prevChrom   string  // CHROM of the most recently fed record
	prevPos     int64   // POS of the most recently fed record; -1 means "no previous"
	headerDone  bool    // true once HeaderChecks has been called
	finished    bool    // true once Finish has been called

	// phaseStyle tracks the phasing separator ("|" or "/")
	// observed in each (sample, PS) pair so far. The key is
	// sample-index and PS-string; the value is one of "phased",
	// "unphased", or "" (no GT). The map resets across
	// chromosomes because phase sets are chromosome-bound per
	// the spec.
	phaseStyle map[phaseKey]string
}

// phaseKey identifies a (sample, PS) pair. PS of "" means
// the implicit default phase set.
type phaseKey struct {
	sampleIdx int
	ps        string
}

// NewValidator constructs a Validator with the given mode.
func NewValidator(m Mode) *Validator { return &Validator{Mode: m} }

// Validate returns a list of issues found by applying the validator to a
// batch of records under the given header. Issues are returned in
// (Line, Record) order; nil is returned for a clean file.
//
// Validate is a thin wrapper around the streaming API: it
// runs the header checks once, feeds every record through
// Next, and then runs Finish. Callers that want to validate
// files larger than memory should use the streaming API
// directly instead.
//
// Validate never modifies the input.
func (v *Validator) Validate(h *Header, records []*Record) []ValidationError {
	v.Reset()
	issues := make([]ValidationError, 0, 16)
	issues = append(issues, v.HeaderChecks(h)...)
	for _, rec := range records {
		issues = append(issues, v.Next(h, rec)...)
	}
	issues = append(issues, v.Finish()...)
	return issues
}

// Reset clears any state from a previous streaming run. It
// does not change the Mode. Call this before reusing a
// Validator for a new file when you used the streaming API.
func (v *Validator) Reset() {
	v.recordIdx = 0
	v.prevChrom = ""
	v.prevPos = -1
	v.headerDone = false
	v.finished = false
	v.phaseStyle = nil
}

// HeaderChecks runs the header-only checks and the
// reserved-key metadata check (ModeStrict). Call this once
// at the start of a streaming run. Subsequent calls are
// no-ops.
func (v *Validator) HeaderChecks(h *Header) []ValidationError {
	if v.headerDone {
		return nil
	}
	v.headerDone = true
	var issues []ValidationError
	issues = append(issues, v.checkHeader(h)...)
	if v.Mode == ModeStrict {
		issues = append(issues, v.checkReservedKeys(h)...)
	}
	return issues
}

// Next runs the per-record checks for a single record and
// updates the cross-record state used by Finish. Returns
// the issues for this record only. Callers that want issues
// in (line, record) order should append the returned slice
// to their own output.
//
// Records must be supplied in the order they appear in the
// file. Records within a chromosome must be sorted by POS
// for the strict-mode POS-monotonicity check to be accurate.
func (v *Validator) Next(h *Header, rec *Record) []ValidationError {
	idx := v.recordIdx
	v.recordIdx++
	issues := v.checkRecord(h, rec, idx)
	// Update cross-record state. We do this even in
	// non-strict mode so a future mode switch is consistent.
	if rec.Chrom != v.prevChrom {
		v.prevChrom = rec.Chrom
		v.prevPos = -1
	}
	if v.Mode == ModeStrict && v.prevPos >= 0 && rec.Pos < v.prevPos {
		issues = append(issues, ValidationError{
			Severity: SeverityError, Record: idx + 1, Line: rec.Line,
			Field: "POS",
			Msg: fmt.Sprintf("POS %d < previous POS %d within CHROM %s",
				rec.Pos, v.prevPos, rec.Chrom),
		})
	}
	if rec.Pos > v.prevPos {
		v.prevPos = rec.Pos
	}

	// Per-sample PS phasing consistency. The spec says all
	// phased genotypes in the same phase set use the same
	// separator style. We track the first observed style per
	// (sample, PS) and flag any later inconsistency.
	issues = append(issues, v.checkPhaseStyle(rec, idx)...)

	return issues
}

// Finish runs any checks that need to see the full record
// stream (currently nothing cross-record — the POS
// monotonicity check runs incrementally in Next). Call this
// once at the end of a streaming run. Subsequent calls are
// no-ops.
func (v *Validator) Finish() []ValidationError {
	if v.finished {
		return nil
	}
	v.finished = true
	return nil
}

// checkHeader runs header-only checks. Records are inspected separately.
func (v *Validator) checkHeader(h *Header) []ValidationError {
	var issues []ValidationError

	if h.FileFormat.Major == 0 && h.FileFormat.Minor == 0 {
		issues = append(issues, ValidationError{
			Severity: SeverityError,
			Field:    "fileformat",
			Msg:      "missing ##fileformat line",
		})
	}

	if v.Mode == ModeStrict && h.FileFormat.Major != 4 {
		issues = append(issues, ValidationError{
			Severity: SeverityError,
			Field:    "fileformat",
			Msg: fmt.Sprintf("unsupported VCF major version %d (only v4.x is supported)",
				h.FileFormat.Major),
		})
	}

	// Duplicate sample IDs (defensive — parser should already catch this).
	seen := map[string]bool{}
	for _, name := range h.SampleNames {
		if seen[name] {
			issues = append(issues, ValidationError{
				Severity: SeverityError,
				Field:    "header/sample",
				Msg:      fmt.Sprintf("duplicate sample name %q", name),
			})
		}
		seen[name] = true
	}

	return issues
}

// checkReservedKeys warns when a reserved INFO/FORMAT key
// (per VCF v4.5 §1.6.1.4–1.6.1.7) is declared in the header
// with a Number or Type that disagrees with the spec. The
// check is non-fatal: the spec lists these as reservations,
// not strict requirements, so we surface a warning the
// caller can promote to an error by switching to ModeStrict.
func (v *Validator) checkReservedKeys(h *Header) []ValidationError {
	var issues []ValidationError
	for key, def := range h.Infos {
		want, ok := reservedInfoKeys[key]
		if !ok {
			continue
		}
		if def.Number.String() != want.Number || def.Type != want.Type {
			issues = append(issues, ValidationError{
				Severity: SeverityWarning,
				Field:    "header/INFO/" + key,
				Msg: fmt.Sprintf("reserved key %q should be Number=%s,Type=%q; got Number=%s,Type=%q",
					key, want.Number, want.Type, def.Number, def.Type),
			})
		}
	}
	for key, def := range h.Formats {
		want, ok := reservedFormatKeys[key]
		if !ok {
			continue
		}
		if def.Number.String() != want.Number || def.Type != want.Type {
			issues = append(issues, ValidationError{
				Severity: SeverityWarning,
				Field:    "header/FORMAT/" + key,
				Msg: fmt.Sprintf("reserved key %q should be Number=%s,Type=%q; got Number=%s,Type=%q",
					key, want.Number, want.Type, def.Number, def.Type),
			})
		}
	}
	return issues
}

// checkRecord runs per-record checks.
func (v *Validator) checkRecord(h *Header, rec *Record, idx int) []ValidationError {
	var issues []ValidationError
	line := rec.Line

	// CHROM may be required to match a declared contig. Spec says it's
	// "recommended" for VCF and "required" for BCF; we emit a warning in
	// default mode, an error in strict.
	if len(h.Contigs) > 0 {
		if _, ok := h.LookupContig(rec.Chrom); !ok {
			iss := ValidationError{
				Severity: SeverityWarning,
				Record:   idx + 1, Line: line,
				Field: "CHROM",
				Msg:   fmt.Sprintf("chromosome %q is not declared in any ##contig line", rec.Chrom),
			}
			if v.Mode == ModeStrict {
				iss.Severity = SeverityError
			}
			issues = append(issues, iss)
		}
	}

	// POS validation (only strict — POS=0 or POS=N+1 are valid telomere
	// markers, and may be intentional).
	if v.Mode == ModeStrict && rec.Pos < 0 {
		issues = append(issues, ValidationError{
			Severity: SeverityError, Record: idx + 1, Line: line,
			Field: "POS",
			Msg:   "negative position",
		})
	}
	// POS / POS+len(REF)-1 must not extend past a declared
	// contig's length (per VCFv4.5 §1.4.7 / §1.6.1.2). We
	// only check when the contig is declared AND has a
	// length; otherwise we cannot constrain it.
	if len(h.Contigs) > 0 {
		if c, ok := h.LookupContig(rec.Chrom); ok && c.Length > 0 {
			end := rec.Pos + int64(len(rec.Ref)) - 1
			if rec.Pos > c.Length {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "POS",
					Msg:   fmt.Sprintf("POS %d exceeds declared contig length %d for %s", rec.Pos, c.Length, rec.Chrom),
				})
			} else if end > c.Length {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: "REF",
					Msg:   fmt.Sprintf("REF extends past declared contig length %d for %s (POS+len(REF)-1 = %d)", c.Length, rec.Chrom, end),
				})
			}
		}
	}

	// REF must contain only A/C/G/T/N (case insensitive) per spec.
	if v.Mode == ModeStrict && !validRefBases(rec.Ref) {
		issues = append(issues, ValidationError{
			Severity: SeverityError, Record: idx + 1, Line: line,
			Field: "REF",
			Msg:   fmt.Sprintf("invalid REF allele %q (must contain only A/C/G/T/N, case insensitive)", rec.Ref),
		})
	}

	// ALT alleles: spec limits them to ACGTN, "*", "<ID>", "<* >", ".",
	// or breakend forms. We only check for forbidden characters in default
	// mode; stricter checks are left to callers that understand SV.
	if v.Mode == ModeStrict {
		for i, a := range rec.Alt {
			if !validAltAllele(a) {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: fmt.Sprintf("ALT[%d]", i),
					Msg:   fmt.Sprintf("invalid alternate allele %q", a),
				})
			}
			if a == rec.Ref {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: line,
					Field: fmt.Sprintf("ALT[%d]", i),
					Msg:   "ALT allele equals REF",
				})
			}
			// Symbolic ALT <ID> must be declared in
			// ##ALT. Breakend forms (containing '[' or
			// ']') and the special <*> allele are
			// exempt — they are defined by the spec
			// rather than by a ##ALT entry.
			if strings.HasPrefix(a, "<") && strings.HasSuffix(a, ">") &&
				a != "<*>" && !strings.ContainsAny(a, "[]") {
				id := a[1 : len(a)-1]
				if _, ok := h.Alts[id]; !ok {
					issues = append(issues, ValidationError{
						Severity: SeverityError, Record: idx + 1, Line: line,
						Field: fmt.Sprintf("ALT[%d]", i),
						Msg:   fmt.Sprintf("symbolic allele %q is not declared in any ##ALT line", a),
					})
				}
			}
		}
	}

	// INFO field checks.
	issues = append(issues, v.checkInfo(h, rec, idx)...)

	// SVCLAIM rules per symbolic allele type. Per
	// VCFv4.5 §1.6.1.7:
	//   - DEL/DUP: SVCLAIM must be specified, allowed D/J/DJ.
	//   - CNV: optional, but if specified must be D or ".".
	//   - INV/INS/breakend: optional, but if specified must
	//     be J, DJ, or ".".
	//   - non-SV: must be "." if present.
	issues = append(issues, v.checkSVClaim(rec, idx)...)

	// EVENTTYPE reserved-value check (§1.6.1.7). Non-
	// reserved values are implementation-defined, so this
	// warns by default and errors in ModeStrict.
	issues = append(issues, v.checkEventType(rec, idx)...)

	// END is a computed field per VCFv4.5 §1.6.1.6. When
	// declared it must match the computed value.
	issues = append(issues, v.checkEND(rec, idx)...)

	// FILTER values must be declared (or "PASS"/".").
	issues = append(issues, v.checkFilter(h, rec, idx)...)

	// FORMAT / sample checks.
	issues = append(issues, v.checkSamples(h, rec, idx)...)

	// GT must be the first FORMAT field if it is declared
	// at all (VCFv4.5 §1.6.2). This is a hard requirement;
	// the parser intentionally does not enforce it so
	// malformed files can still be inspected.
	if v.Mode == ModeStrict && len(rec.Format) > 0 && hasGT(rec.Format) {
		if rec.Format[0] != "GT" {
			iss := ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "FORMAT",
				Msg:   fmt.Sprintf("GT must be the first FORMAT field when present (got %q)", rec.Format[0]),
			}
			issues = append(issues, iss)
		}
	}

	return issues
}

// checkInfo validates each INFO key against the header declaration.
func (v *Validator) checkInfo(h *Header, rec *Record, idx int) []ValidationError {
	var issues []ValidationError
	altCount := len(rec.Alt)
	line := rec.Line

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
						Msg: fmt.Sprintf("Integer value %d is in the BCF-unsafe range [-2^31, -2^31+7]", n),
					})
				}
			}
		}

		// SVLEN sign rule: for symbolic SV alleles the length
		// must be positive. (The spec says SVLEN is always
		// positive, even for INV.)
		if key == "SVLEN" {
			for i, n := range val.Ints {
				if n <= 0 {
					alt := "."
					if i < len(rec.Alt) {
						alt = rec.Alt[i]
					}
					issues = append(issues, ValidationError{
						Severity: SeverityError, Record: idx + 1, Line: line,
						Field: "INFO/SVLEN",
						Msg: fmt.Sprintf("SVLEN[%d] for ALT %q must be positive (got %d)", i, alt, n),
					})
				}
			}
		}

		// CIPOS / CIEND / CILEN confidence intervals: the
		// interval must span 0. The field is a list of pairs
		// (low, high) per ALT allele. The lower bound cannot
		// be greater than 0, and the upper bound cannot be
		// less than 0.
		if key == "CIPOS" || key == "CIEND" || key == "CILEN" {
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
	}
	return issues
}

// reservedEventTypes is the set of EVENTTYPE values the spec
// reserves (VCFv4.5 §1.6.1.7 "The following EVENTTYPE values
// are reserved and should be used when appropriate"). The
// spec also says "The semantics of other EVENTTYPE values
// is implementation-defined" — the validator warns in
// ModeStrict for non-reserved values, but they are not
// strictly invalid.
var reservedEventTypes = map[string]bool{
	"DEL":             true,
	"DEL:ME":          true,
	"INS":             true,
	"INS:ME":          true,
	"DUP":             true,
	"DUP:TANDEM":      true,
	"DUP:DISPERSED":   true,
	"INV":             true,
	"TRA":             true,
	"TRA:BALANCED":    true,
	"TRA:UNBALANCED":  true,
	"CHROMOTHRIPSIS":  true,
	"CHROMOPLEXY":     true,
	"BFB":             true,
	"DOUBLEMINUTE":    true,
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
func (v *Validator) checkEventType(rec *Record, idx int) []ValidationError {
	raw, ok := rec.Info["EVENTTYPE"]
	if !ok {
		return nil
	}
	var issues []ValidationError
	for i, et := range raw.Strs {
		if et == "" || et == "." {
			continue
		}
		if !IsReservedEventType(et) {
			iss := ValidationError{
				Severity: SeverityWarning,
				Record:    idx + 1, Line: rec.Line,
				Field:    "INFO/EVENTTYPE",
				Msg:      fmt.Sprintf("EVENTTYPE[%d]=%q is not in the spec's reserved set (semantics are implementation-defined)", i, et),
			}
			if v.Mode == ModeStrict {
				iss.Severity = SeverityError
			}
			issues = append(issues, iss)
		}
	}
	return issues
}

// checkPhaseStyle enforces the PS phasing-consistency rule
// from VCFv4.5 §"PS": records that share a PS value (and
// the same sample + chromosome) must use the same GT
// separator style. We only check records that declare an
// explicit PS — the spec is silent on records without PS
// (the implicit default phase set), and the spec example
// itself mixes separators on records without PS.
func (v *Validator) checkPhaseStyle(rec *Record, idx int) []ValidationError {
	if v.Mode != ModeStrict {
		return nil
	}
	if rec.Chrom != v.prevChrom {
		v.phaseStyle = nil
	}
	if v.phaseStyle == nil {
		v.phaseStyle = map[phaseKey]string{}
	}
	var issues []ValidationError
	for sIdx, g := range rec.Samples {
		psv, ok := g.Fields["PS"]
		if !ok || psv.Missing || len(psv.Ints) == 0 {
			continue
		}
		ps := fmt.Sprintf("%d", psv.Ints[0])
		vv, ok := g.Fields["GT"]
		if !ok || vv.Missing || len(vv.Strs) == 0 {
			continue
		}
		gt := vv.Strs[0]
		if !strings.Contains(gt, "|") && !strings.Contains(gt, "/") {
			continue
		}
		style := "unphased"
		if strings.Contains(gt, "|") {
			style = "phased"
		}
		key := phaseKey{sampleIdx: sIdx, ps: ps}
		if prev, seen := v.phaseStyle[key]; seen && prev != style {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: rec.Line,
				Field: fmt.Sprintf("FORMAT/GT/sample[%d]", sIdx),
				Msg: fmt.Sprintf("GT separator style for sample %d, PS %q changed from %s to %s within CHROM %s",
					sIdx, ps, prev, style, rec.Chrom),
			})
		}
		v.phaseStyle[key] = style
	}
	return issues
}

// hasGT reports whether s contains a "GT" entry.
func hasGT(s []string) bool {
	for _, k := range s {
		if k == "GT" {
			return true
		}
	}
	return false
}

// checkEND enforces VCFv4.5 §1.6.1.6: END is a computed
// field. When the header declares ##INFO=<ID=END,...> and
// the record carries an END value, the value must match
// the computed END (POS+len(REF)-1, plus SVLEN / LEN
// adjustments for symbolic ALT alleles).
func (v *Validator) checkEND(rec *Record, idx int) []ValidationError {
	if v.Mode != ModeStrict {
		return nil
	}
	raw, ok := rec.Info["END"]
	if !ok || raw.Missing || len(raw.Ints) == 0 {
		return nil
	}
	want := computeEND(rec)
	if want <= 0 {
		return nil
	}
	got := raw.Ints[0]
	if got != want {
		return []ValidationError{{
			Severity: SeverityError, Record: idx + 1, Line: rec.Line,
			Field:    "INFO/END",
			Msg:      fmt.Sprintf("declared END %d does not match the computed END %d (POS+len(REF)-1, plus SVLEN/LEN adjustments)", got, want),
		}}
	}
	return nil
}

// checkSVClaim enforces the per-symbolic-allele SVCLAIM
// rules from VCFv4.5 §1.6.1.7.
//
//   - DEL/DUP: SVCLAIM must be specified and one of D, J, DJ.
//   - CNV: optional, but if specified must be D or ".".
//   - INV/INS/breakend alleles: optional, but if specified
//     must be J, DJ, or ".".
//   - non-SV (e.g. SNV) alleles: must be "." if present.
//
// The check is a no-op when the record has no SVCLAIM key
// in its INFO column.
func (v *Validator) checkSVClaim(rec *Record, idx int) []ValidationError {
	raw, ok := rec.Info["SVCLAIM"]
	// SVCLAIM is stored on Strs (the spec declares it as
	// Type=String). Fall back to Ints for hand-built records
	// (BCF-encoded claims).
	claims := raw.Strs
	if len(claims) == 0 {
		for _, n := range raw.Ints {
			claims = append(claims, intToSVClaim(n))
		}
	}
	if !ok || len(claims) == 0 {
		// SVCLAIM is missing: spec says it must be present
		// for DEL/DUP records. Flag here only when those
		// alleles exist.
		for _, alt := range rec.Alt {
			if alt == "<DEL>" || alt == "<DUP>" {
				return []ValidationError{{
					Severity: SeverityError, Record: idx + 1, Line: rec.Line,
					Field: "INFO/SVCLAIM",
					Msg:   "SVCLAIM must be specified for DEL and DUP symbolic alleles",
				}}
			}
		}
		return nil
	}
	var issues []ValidationError
	for i, alt := range rec.Alt {
		if i >= len(claims) {
			break
		}
		claim := claims[i]
		want := svClaimAllowed(alt)
		if want == "" {
			// Non-SV allele: claim must be "." if present.
			if claim != "." {
				issues = append(issues, ValidationError{
					Severity: SeverityError, Record: idx + 1, Line: rec.Line,
					Field: "INFO/SVCLAIM",
					Msg: fmt.Sprintf("non-SV allele %q must have SVCLAIM \".\" (got %q)", alt, claim),
				})
			}
			continue
		}
		// Symbolic SV: claim must be one of want.
		if !svClaimInSet(claim, want) {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: rec.Line,
				Field: "INFO/SVCLAIM",
				Msg: fmt.Sprintf("SVCLAIM[%d] for ALT %q must be one of %s (got %q)", i, alt, want, claim),
			})
		}
	}
	return issues
}

// intToSVClaim decodes a BCF-encoded SVCLAIM int back into
// the spec's letter form. The BCF v1 spec uses the values
// 0 (PASS / no claim), 1 (D), 2 (J), 3 (DJ). Anything else
// is reported as a "?(<n>)" string the validator surfaces
// verbatim.
func intToSVClaim(n int64) string {
	switch n {
	case 0:
		return "."
	case 1:
		return "D"
	case 2:
		return "J"
	case 3:
		return "DJ"
	}
	return fmt.Sprintf("?(%d)", n)
}

// svClaimAllowed returns the set of spec-allowed SVCLAIM
// letters for the given ALT allele. Empty string means the
// allele is non-SV and only "." is allowed.
func svClaimAllowed(alt string) string {
	// Symbolic structural variants: <ID> form.
	if !strings.HasPrefix(alt, "<") || !strings.HasSuffix(alt, ">") {
		return ""
	}
	id := alt[1 : len(alt)-1]
	switch id {
	case "DEL", "DUP":
		return "{D, J, DJ}"
	case "CNV":
		return "{D, .}"
	case "INV", "INS", "*":
		return "{J, DJ, .}"
	}
	// Breakend alleles: <DEL> / <DUP> with qualifiers, or
	// the bracket-form breakends. Per spec they are
	// implicit J claims.
	if isBreakend(alt) {
		return "{J, .}"
	}
	// User-defined symbolic alleles: spec is silent; allow
	// all four codes plus "." to avoid false positives.
	return "{D, J, DJ, .}"
}

func svClaimInSet(claim, allowed string) bool {
	if claim == "." {
		return strings.Contains(allowed, ".")
	}
	if claim == "D" {
		return strings.Contains(allowed, "D")
	}
	if claim == "J" {
		return strings.Contains(allowed, "J")
	}
	if claim == "DJ" {
		return strings.Contains(allowed, "DJ")
	}
	return false
}

// checkFilter validates FILTER column values against declared filters.
//
// Per spec, "0" is reserved and must not be used as a filter string; we
// flag it as an error in both modes.
func (v *Validator) checkFilter(h *Header, rec *Record, idx int) []ValidationError {
	var issues []ValidationError
	line := rec.Line
	for _, f := range rec.Filter {
		if f == "PASS" || f == "." {
			continue
		}
		if f == "0" {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "FILTER",
				Msg:   `"0" is not a valid filter value (reserved)`,
			})
			continue
		}
		if _, ok := h.Filters[f]; !ok {
			iss := ValidationError{
				Severity: SeverityWarning, Record: idx + 1, Line: line,
				Field: "FILTER",
				Msg:   fmt.Sprintf("filter %q is not declared in any ##FILTER line", f),
			}
			if v.Mode == ModeStrict {
				iss.Severity = SeverityError
			}
			issues = append(issues, iss)
		}
	}
	return issues
}

// checkSamples validates the FORMAT column and per-sample values.
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
		}
		// LAA value range and distinctness (per VCFv4.5
		// §1.6.1.4). LAA must be distinct 1-based indices
		// into ALT.
		if v, ok := g.Fields["LAA"]; ok && !v.Missing {
			issues = append(issues, checkLAA(&v, rec, idx, sIdx, line)...)
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
			// derive ploidy when possible.
			ploidy := ploidyFromGT(g)
			want := def.Number.Cardinality(len(rec.Alt), ploidy)
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
//   1. A single record may not mix / and | separators (the
//      phasing must be consistent within a record).
//   2. Each allele index must be a valid 0-based index into
//      REF+ALT. The missing value "." is always allowed.
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
