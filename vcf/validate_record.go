package vcf

import (
	"fmt"
	"strings"
)

// Record-level rules (header, EVENTTYPE, ploidy, END bounds,
// GVCF block ref, SVLEN bounds, SVCLAIM, FILTER, breakend
// positions, reference matching) and their helpers live here.
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
				Record:   idx + 1, Line: rec.Line,
				Field: "INFO/EVENTTYPE",
				Msg:   fmt.Sprintf("EVENTTYPE[%d]=%q is not in the spec's reserved set (semantics are implementation-defined)", i, et),
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

// isMitochondrial reports whether chrom names a
// mitochondrial contig. VCFv4.5 §1.6.1.6 calls out
// mitochondrial DNA (chrM / MT / M / etc.) as a special
// case: the reference is circular, POS can exceed the
// contig length, and END can be smaller than POS. The
// list covers the common conventions; callers using
// non-standard names can set a ##contig attribute in
// the future.
func isMitochondrial(chrom string) bool {
	switch chrom {
	case "MT", "Mt", "mt", "chrM", "chrm", "M", "chrMT", "chrMt", "MTDNA", "mtdna":
		return true
	}
	return false
}

// checkPloidy enforces VCFv4.5 §1.6.2: all samples in a
// record must have the same ploidy. Mixed ploidies in a
// single record make Number=G and Number=P cardinality
// ambiguous, so the spec forbids the combination. The
// check is a no-op when fewer than two samples have a
// present GT, or when every present GT has the same
// ploidy.
func (v *Validator) checkPloidy(rec *Record, idx int) []ValidationError {
	if v.Mode != ModeStrict {
		return nil
	}
	if len(rec.Samples) < 2 {
		return nil
	}
	var firstPloidy int
	firstSeen := false
	for sIdx, g := range rec.Samples {
		p := ploidyFromGT(g)
		if p <= 0 {
			continue
		}
		if !firstSeen {
			firstPloidy = p
			firstSeen = true
			continue
		}
		if p != firstPloidy {
			return []ValidationError{{
				Severity: SeverityError, Record: idx + 1, Line: rec.Line,
				Field: "FORMAT/GT",
				Msg:   fmt.Sprintf("samples have different ploidies in a single record (sample 0 ploidy %d, sample %d ploidy %d)", firstPloidy, sIdx, p),
			}}
		}
	}
	return nil
}

// checkEND enforces VCFv4.5 §1.6.1.6: END is a computed
// field. When the header declares ##INFO=<ID=END,...> and
// the record carries an END value, the value must match
// the computed END (POS+len(REF)-1, plus SVLEN / LEN
// adjustments for symbolic ALT alleles). Mitochondrial
// chromosomes are exempt: the reference is circular and
// END may legitimately wrap past the origin.
func (v *Validator) checkEND(rec *Record, idx int) []ValidationError {
	if v.Mode != ModeStrict {
		return nil
	}
	if isMitochondrial(rec.Chrom) {
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
	// A gVCF reference block's END is POS+LEN-1, but the
	// spec's own §1.7 example uses POS+LEN in three of its
	// seven rows and POS+LEN-1 in the rest, so both forms are
	// accepted for a <*> block. (The example also has a
	// semicolon where the FORMAT separator should be a colon
	// in three rows, and a missing PL in one FORMAT column -
	// typos in the spec, not defects on our side.)
	if got != want && (got != want+1 || !hasReferenceBlock(rec)) {
		return []ValidationError{{
			Severity: SeverityError, Record: idx + 1, Line: rec.Line,
			Field: "INFO/END",
			Msg:   fmt.Sprintf("declared END %d does not match the computed END %d (POS+len(REF)-1, plus SVLEN/LEN adjustments)", got, want),
		}}
	}
	return nil
}

// checkENDBounds fires when a FASTA reference is attached
// and the record declares an INFO END value. END must be
// within the contig (END <= ref.Length). Without a
// reference, this check is skipped (the existing
// checkEND only enforces value consistency, not bounds).
//
// For mitochondrial chromosomes END may legitimately
// wrap past the origin, so the check is skipped there
// (matching checkEND's existing exception).
func (v *Validator) checkENDBounds(rec *Record, idx, line int) []ValidationError {
	if isMitochondrial(rec.Chrom) {
		return nil
	}
	raw, ok := rec.Info["END"]
	if !ok || raw.Missing || len(raw.Ints) == 0 {
		return nil
	}
	end := raw.Ints[0]
	contigLen, err := v.Reference.Length(rec.Chrom)
	if err != nil {
		return nil
	}
	if end > contigLen {
		return []ValidationError{{
			Severity: SeverityError, Record: idx + 1, Line: line,
			Field: "INFO/END",
			Msg:   fmt.Sprintf("END %d exceeds contig length %d for %s", end, contigLen, rec.Chrom),
		}}
	}
	if end < rec.Pos {
		return []ValidationError{{
			Severity: SeverityError, Record: idx + 1, Line: line,
			Field: "INFO/END",
			Msg:   fmt.Sprintf("END %d is before POS %d for %s", end, rec.Pos, rec.Chrom),
		}}
	}
	return nil
}

// isNonVariantBlock reports whether rec is a gVCF
// non-variant block per §1.7: ALT is exactly the
// "<*>" symbolic allele (which is what the parser
// normalises "<NON_REF>" to) AND the record carries
// an END value.
//
// ALT containing "<*>" alongside real ALTs is NOT a
// non-variant block — that's a normal variant record
// that lists the unspecified allele as one of its
// alternatives.
func isNonVariantBlock(rec *Record) bool {
	if len(rec.Alt) != 1 || rec.Alt[0] != "<*>" {
		return false
	}
	raw, ok := rec.Info["END"]
	if !ok || raw.Missing || len(raw.Ints) == 0 {
		return false
	}
	return true
}

// infoEnd returns the END value of rec, or 0 if
// END is missing or non-integer. Callers should gate
// on isNonVariantBlock first.
func infoEnd(rec *Record) int64 {
	raw, ok := rec.Info["END"]
	if !ok || raw.Missing || len(raw.Ints) == 0 {
		return 0
	}
	return raw.Ints[0]
}

// checkGVCFBlock enforces the per-record gVCF rules:
//
//   - ALT is "<*>" but no END declared.
//   - END < POS.
//
// Cross-record overlap (next record's POS inside the
// previous block) is checked in the Next() state
// machine using prevEnd.
func (v *Validator) checkGVCFBlock(rec *Record, idx int) []ValidationError {
	if len(rec.Alt) != 1 || rec.Alt[0] != "<*>" {
		return nil
	}
	// Single-element "<*>" ALT — check that END is
	// declared. Without END the record is just a
	// poorly-formed gVCF block.
	raw, hasEND := rec.Info["END"]
	if !hasEND || raw.Missing || len(raw.Ints) == 0 {
		if v.Mode == ModeStrict {
			return []ValidationError{{
				Severity: SeverityError, Record: idx + 1, Line: rec.Line,
				Field: "ALT",
				Msg:   `ALT="<*>" (gVCF non-variant block) requires an INFO END field`,
			}}
		}
		// ModeDefault: surface as a warning so
		// callers can see the issue without
		// changing the default-mode behaviour.
		return []ValidationError{{
			Severity: SeverityWarning, Record: idx + 1, Line: rec.Line,
			Field: "ALT",
			Msg:   `ALT="<*>" (gVCF non-variant block) requires an INFO END field`,
		}}
	}
	if raw.Ints[0] < rec.Pos {
		return []ValidationError{{
			Severity: SeverityError, Record: idx + 1, Line: rec.Line,
			Field: "INFO/END",
			Msg:   fmt.Sprintf("END %d is before POS %d (non-variant block cannot go backwards)", raw.Ints[0], rec.Pos),
		}}
	}
	return nil
}

// checkGVCFBlockRef fires when a FASTA reference is
// attached and rec is a gVCF non-variant block. The
// REF column must match the reference sequence over
// the entire [POS, END] block — every position in the
// block is reference-matching by definition.
func (v *Validator) checkGVCFBlockRef(rec *Record, idx, line int) []ValidationError {
	if !isNonVariantBlock(rec) {
		return nil
	}
	end := infoEnd(rec)
	if end <= rec.Pos {
		// Already flagged by checkGVCFBlock.
		return nil
	}
	pos := int(rec.Pos)
	bases, err := v.Reference.BaseRange(rec.Chrom, pos, int(end)+1)
	if err != nil {
		// Defer to checkRefMatchesReference-style
		// behaviour: warn (or error in strict)
		// and skip.
		sev := SeverityWarning
		if v.Mode == ModeStrict {
			sev = SeverityError
		}
		return []ValidationError{{
			Severity: sev, Record: idx + 1, Line: line,
			Field: "REF",
			Msg:   fmt.Sprintf("reference lookup for %s:%d-%d failed: %v; skipping gVCF block-REF check", rec.Chrom, pos, end, err),
		}}
	}
	wantLen := int(end) - pos + 1
	if len(bases) != wantLen {
		sev := SeverityWarning
		if v.Mode == ModeStrict {
			sev = SeverityError
		}
		return []ValidationError{{
			Severity: sev, Record: idx + 1, Line: line,
			Field: "REF",
			Msg:   fmt.Sprintf("reference lookup for %s:%d-%d returned %d bases, want %d; skipping gVCF block-REF check", rec.Chrom, pos, end, len(bases), wantLen),
		}}
	}
	var issues []ValidationError
	for i, refBase := range bases {
		if i >= len(rec.Ref) {
			// REF is shorter than the block — that
			// would imply a deletion within the
			// block, which a gVCF block cannot
			// carry.
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "REF",
				Msg:   fmt.Sprintf("gVCF block REF has length %d, but block range [%d,%d] requires %d bases", len(rec.Ref), pos, end, wantLen),
			})
			break
		}
		want := strings.ToUpper(string(rec.Ref[i]))
		got := strings.ToUpper(string(refBase))
		if want != got {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "REF",
				Msg:   fmt.Sprintf("REF[%d]=%q does not match reference base %q at %s:%d (gVCF block must be reference-matching)", i, string(rec.Ref[i]), string(refBase), rec.Chrom, pos+i),
			})
		}
	}
	return issues
}

// checkSVLENBounds fires when a FASTA reference is attached
// and the record declares INFO SVLEN. For each non-missing
// entry, the implied end position (POS + |SVLEN| - 1)
// must be within the contig. The absolute-value rule
// follows §1.6.1.7 ("the absolute value of SVLEN should
// be taken and a negative SVLEN should be treated as
// positive").
func (v *Validator) checkSVLENBounds(rec *Record, idx, line int) []ValidationError {
	if isMitochondrial(rec.Chrom) {
		return nil
	}
	iv, ok := rec.Info["SVLEN"]
	if !ok || iv.Missing || len(iv.Ints) == 0 {
		return nil
	}
	contigLen, err := v.Reference.Length(rec.Chrom)
	if err != nil {
		return nil
	}
	var issues []ValidationError
	for i, svlen := range iv.Ints {
		if isMissing(svlen) {
			continue
		}
		absVal := svlen
		if absVal < 0 {
			absVal = -absVal
		}
		impliedEnd := rec.Pos + absVal - 1
		if impliedEnd > contigLen {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: line,
				Field: "INFO/SVLEN",
				Msg: fmt.Sprintf("SVLEN[%d]=%d implies end position %d which exceeds contig length %d for %s",
					i, svlen, impliedEnd, contigLen, rec.Chrom),
			})
		}
	}
	return issues
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
					Msg:   fmt.Sprintf("non-SV allele %q must have SVCLAIM \".\" (got %q)", alt, claim),
				})
			}
			continue
		}
		// Symbolic SV: claim must be one of want.
		if !svClaimInSet(claim, want) {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: rec.Line,
				Field: "INFO/SVCLAIM",
				Msg:   fmt.Sprintf("SVCLAIM[%d] for ALT %q must be one of %s (got %q)", i, alt, want, claim),
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

// hasReferenceBlock reports whether rec carries a <*> (or
// <NON_REF>) ALT allele, i.e. a gVCF reference block.
func hasReferenceBlock(rec *Record) bool {
	for _, alt := range rec.Alt {
		if !strings.HasPrefix(alt, "<") || !strings.HasSuffix(alt, ">") {
			continue
		}
		switch alt[1 : len(alt)-1] {
		case "*", "NON_REF":
			return true
		}
	}
	return false
}
