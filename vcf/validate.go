package vcf

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"

	"github.com/inful/htsfile/fasta"
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

	// Reference is an optional FASTA handle. When set, the
	// validator enforces VCFv4.5 §1.6.1.1 "For simple
	// insertions and deletions in which either the REF or
	// one of the ALT alleles would otherwise be
	// null/empty, the REF and ALT Strings must include the
	// base before the variant" by checking that REF[0]
	// matches the reference base at POS-1 (or POS when
	// POS=1).
	Reference *fasta.Reference

	// Streaming state (used by Next / Finish). Zero value is
	// fine; the validator initialises lazily.
	recordIdx  int    // 0-based index of the most recently fed record
	prevChrom  string // CHROM of the most recently fed record
	prevPos    int64  // POS of the most recently fed record; -1 means "no previous"
	headerDone bool   // true once HeaderChecks has been called
	finished   bool   // true once Finish has been called

	// phaseStyle tracks the phasing separator ("|" or "/")
	// observed in each (sample, PS) pair so far. The key is
	// sample-index and PS-string; the value is one of "phased",
	// "unphased", or "" (no GT). The map resets across
	// chromosomes because phase sets are chromosome-bound per
	// the spec.
	phaseStyle map[phaseKey]string

	// seenChroms tracks chromosomes observed so far so the
	// validator can enforce the §1.6.1 "contiguous
	// block" rule. The map's keys are CHROM strings; the
	// values are unused.
	seenChroms map[string]struct{}

	// seenIDs tracks the data-record IDs observed so far
	// so the validator can enforce the §1.6.1 "No
	// identifier should be present in more than one data
	// record" rule. The map keys are ID strings; the
	// values are unused.
	seenIDs map[string]struct{}

	// mateRecords maps a MATEID to the first record
	// that declared it, so we can later compare the
	// PSL/PSO of breakend mate pairs (§1.6.2 PSO).
	mateRecords map[string]*mateSnapshot
}

// mateSnapshot captures the PSL and PSO of the first
// record encountered with a given MATEID. The second
// record with the same MATEID is checked against this
// snapshot for the §1.6.2 rule:
//
//	"Variants in breakpoint notation must have the
//	 same PSL and PSO on both records."
type mateSnapshot struct {
	psl []string
	pso []int64
}

// phaseKey identifies a (sample, PS) pair. PS of "" means
// the implicit default phase set.
type phaseKey struct {
	sampleIdx int
	ps        string
}

// NewValidator constructs a Validator with the given mode.
func NewValidator(m Mode) *Validator { return &Validator{Mode: m} }

// WithReference attaches a FASTA reference to the validator
// for the §1.6.1.1 REF/ALT padding-base check. Returns the
// receiver for chaining. The reference must remain open for
// the lifetime of validation.
func (v *Validator) WithReference(ref *fasta.Reference) *Validator {
	v.Reference = ref
	return v
}

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
	v.seenChroms = nil
	v.seenIDs = nil
	v.mateRecords = nil
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

	// CHROM contiguous block (VCFv4.5 §1.6.1: "All
	// entries for a specific CHROM must form a contiguous
	// block within the VCF file"). Checked before
	// prevChrom is updated so we know the actual previous
	// CHROM. A CHROM seen earlier and not the immediate
	// previous chromosome is a reorder.
	if v.seenChroms == nil {
		v.seenChroms = make(map[string]struct{})
	}
	if _, seen := v.seenChroms[rec.Chrom]; seen && rec.Chrom != v.prevChrom {
		issues = append(issues, ValidationError{
			Severity: SeverityError, Record: idx + 1, Line: rec.Line,
			Field: "CHROM",
			Msg:   fmt.Sprintf("chromosome %q was already seen earlier in the file; CHROM blocks must be contiguous", rec.Chrom),
		})
	}
	v.seenChroms[rec.Chrom] = struct{}{}

	// ID uniqueness across records (§1.6.1: "No
	// identifier should be present in more than one data
	// record"). The ID column is a list of IDs separated
	// by ';' within a single record; '.' in any position
	// is the missing value. We track every non-missing ID
	// seen so far; the second is an error.
	if v.seenIDs == nil {
		v.seenIDs = make(map[string]struct{})
	}
	for _, id := range rec.ID {
		if id == "." {
			continue
		}
		if _, seen := v.seenIDs[id]; seen {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: rec.Line,
				Field: "ID",
				Msg:   fmt.Sprintf("ID %q was already used in an earlier data record", id),
			})
		}
		v.seenIDs[id] = struct{}{}
	}

	// Breakend mate PSL / PSO consistency (§1.6.2 PSO:
	// "Variants in breakpoint notation must have the
	// same PSL and PSO on both records"). Records must
	// carry a MATEID; only records whose ALT is a
	// breakend form participate. We key by BOTH the
	// record's own IDs and its MATEID so the first mate
	// to arrive registers under both keys, and the
	// second mate to arrive finds the first regardless
	// of arrival order.
	for _, a := range rec.Alt {
		if !isBreakend(a) {
			continue
		}
		mateID, ok := rec.Info["MATEID"]
		if !ok || mateID.Missing || len(mateID.Strs) == 0 {
			continue
		}
		mid := mateID.Strs[0]
		if mid == "" {
			continue
		}
		if v.mateRecords == nil {
			v.mateRecords = make(map[string]*mateSnapshot)
		}
		pslFV := rec.Info["PSL"]
		psoFV := rec.Info["PSO"]
		var psl []string
		var pso []int64
		if pslFV.HasValue() {
			psl = pslFV.Strs
		}
		if psoFV.HasValue() {
			pso = psoFV.Ints
		}
		snap := &mateSnapshot{psl: psl, pso: pso}
		// Check against any of the record's own IDs.
		matched := false
		for _, id := range rec.ID {
			if id == "." {
				continue
			}
			if existing, seen := v.mateRecords[id]; seen {
				checkMateMatch(existing, snap, mid, idx, rec.Line, &issues)
				matched = true
				break
			}
		}
		if !matched {
			if existing, seen := v.mateRecords[mid]; seen {
				checkMateMatch(existing, snap, mid, idx, rec.Line, &issues)
			}
		}
		// Register the record under its own IDs and its
		// MATEID so the mate finds it whichever arrives
		// first.
		for _, id := range rec.ID {
			if id == "." {
				continue
			}
			if _, already := v.mateRecords[id]; !already {
				v.mateRecords[id] = snap
			}
		}
		if _, already := v.mateRecords[mid]; !already {
			v.mateRecords[mid] = snap
		}
	}

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

	// ##META Values enforcement (§1.4.7 "Sample field
	// format"): every ##SAMPLE field whose key matches
	// a ##META ID that declares a Values list must have
	// a value drawn from that list. Fires in ModeStrict;
	// skipped in ModeDefault to avoid surprising callers
	// who don't want strict header validation.
	if v.Mode == ModeStrict && len(h.Metas) > 0 {
		for _, sample := range h.Samples {
			for key, rawValue := range sample.Other {
				meta, ok := h.Metas[key]
				if !ok || meta.Values == nil {
					continue
				}
				if !valuesContain(meta.Values, rawValue) {
					issues = append(issues, ValidationError{
						Severity: SeverityWarning,
						Field:    fmt.Sprintf("header/META/%s", key),
						Msg: fmt.Sprintf("##SAMPLE %s uses value %q which is not in the ##META %s.Values list (%v)",
							sample.ID, rawValue, key, meta.Values),
					})
				}
			}
		}
	}

	return issues
}

// valuesContain reports whether v contains s (string
// equality; empty-string and missing-value are skipped
// at the call site).
func valuesContain(values []string, s string) bool {
	for _, v := range values {
		if v == s {
			return true
		}
	}
	return false
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
	// Mitochondrial chromosomes are exempt: the reference is
	// circular (§1.6.1.6), so POS can wrap past the declared
	// length.
	if len(h.Contigs) > 0 && !isMitochondrial(rec.Chrom) {
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

	// Full REF match against the reference (§1.6.1.1).
	// For every byte of REF, the validator compares it
	// against ref[POS+i]. This catches wrong contig,
	// wrong POS, transcription errors, and undeclared
	// IUPAC codes — all in one check. It strictly
	// subsumes the older single-base padding check.
	if v.Reference != nil && rec.Ref != "" {
		issues = append(issues, v.checkRefMatchesReference(rec, idx, line)...)
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
			// rather than by a ##ALT entry. Bare IUPAC
			// ambiguity codes (R, Y, S, W, K, M, B, D,
			// H, V) are checked the same way; they are
			// also only valid when declared via ##ALT.
			switch {
			case a == "<*>" || strings.ContainsAny(a, "[]"):
				// skip
			case strings.HasPrefix(a, "<") && strings.HasSuffix(a, ">"):
				id := a[1 : len(a)-1]
				if _, ok := h.Alts[id]; !ok {
					issues = append(issues, ValidationError{
						Severity: SeverityError, Record: idx + 1, Line: line,
						Field: fmt.Sprintf("ALT[%d]", i),
						Msg:   fmt.Sprintf("symbolic allele %q is not declared in any ##ALT line", a),
					})
				}
			case iupacIDs[a]:
				if _, ok := h.Alts[a]; !ok {
					issues = append(issues, ValidationError{
						Severity: SeverityError, Record: idx + 1, Line: line,
						Field: fmt.Sprintf("ALT[%d]", i),
						Msg:   fmt.Sprintf("IUPAC code %q is not declared in any ##ALT line", a),
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

	// Reference-aware bounds: with a FASTA reference
	// attached, verify END and SVLEN stay inside the
	// contig. Without a reference, these checks are
	// skipped (the existing checks rely on declared
	// ##contig lengths).
	if v.Reference != nil {
		issues = append(issues, v.checkENDBounds(rec, idx, line)...)
		issues = append(issues, v.checkSVLENBounds(rec, idx, line)...)
		// Breakend remote positions must lie within the
		// remote contig (§1.6.5.1).
		issues = append(issues, v.checkBreakendPositions(rec, idx, line)...)
	}

	// FILTER values must be declared (or "PASS"/".").
	issues = append(issues, v.checkFilter(h, rec, idx)...)

	// FORMAT / sample checks.
	issues = append(issues, v.checkPloidy(rec, idx)...)
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
	if got != want {
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
// and appends SeverityErrors on mismatch.
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
