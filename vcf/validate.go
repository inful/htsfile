package vcf

import (
	"fmt"
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
	recordIdx int    // 0-based index of the most recently fed record
	prevChrom string // CHROM of the most recently fed record
	prevPos   int64  // POS of the most recently fed record; -1 means "no previous"
	// prevEnd is the END of the previous record if it
	// was a gVCF non-variant block, else 0. Combined
	// with prevPos it lets the gVCF overlap check fire
	// when the next record's POS would land inside a
	// previous block.
	prevEnd    int64
	headerDone bool // true once HeaderChecks has been called
	finished   bool // true once Finish has been called

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

// mateSnapshot captures PSL/PSO plus the local
// position and parsed breakend remote positions of
// the first record encountered with a given MATEID.
// The second record with the same MATEID is checked
// against this snapshot for the §1.6.2 PSL/PSO rule
// and the §1.6.5.1 breakend position round-trip
// rule.
type mateSnapshot struct {
	psl       []string
	pso       []int64
	chrom     string        // local CHROM of the first mate seen
	pos       int64         // local POS of the first mate seen
	breakends []breakendPos // parsed remote positions of the first mate's breakend ALTs
}

// breakendPos is one parsed "chr:pos" from a
// breakend replacement string (§1.6.5.1).
type breakendPos struct {
	contig string
	pos    int64
}

// recordBreakends returns the parsed (contig, pos)
// for each breakend ALT in rec. ALTs that aren't
// breakends, or whose position is omitted (a telomere
// breakpoint), are skipped.
func recordBreakends(rec *Record) []breakendPos {
	var out []breakendPos
	for _, a := range rec.Alt {
		if !isBreakend(a) {
			continue
		}
		contig, pos, ok := parseBreakendPosition(a)
		if !ok {
			continue
		}
		out = append(out, breakendPos{contig: contig, pos: pos})
	}
	return out
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
	v.prevEnd = 0
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
		snap := &mateSnapshot{
			psl:       psl,
			pso:       pso,
			chrom:     rec.Chrom,
			pos:       rec.Pos,
			breakends: recordBreakends(rec),
		}
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
		v.prevEnd = 0
	}
	// gVCF non-variant block overlap (§1.7): if the
	// previous record was a block, the new record's
	// POS must be > previous block's END. This is a
	// tighter check than the per-contig POS monotonic
	// rule, which only ensures POS > prevPos.
	if v.prevEnd > 0 && rec.Pos <= v.prevEnd {
		if v.Mode == ModeStrict {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: idx + 1, Line: rec.Line,
				Field: "POS",
				Msg: fmt.Sprintf("POS %d overlaps previous non-variant block ending at %d within CHROM %s",
					rec.Pos, v.prevEnd, rec.Chrom),
			})
		}
	}
	if v.Mode == ModeStrict && v.prevPos >= 0 && rec.Pos < v.prevPos {
		issues = append(issues, ValidationError{
			Severity: SeverityError, Record: idx + 1, Line: rec.Line,
			Field: "POS",
			Msg: fmt.Sprintf("POS %d < previous POS %d within CHROM %s",
				rec.Pos, v.prevPos, rec.Chrom),
		})
	}
	// For non-variant blocks, advance prevPos to the
	// END (the effective "next free position" on this
	// chromosome) and set prevEnd so subsequent
	// records can be checked against the block. For
	// other records, just advance prevPos to POS.
	if isNonVariantBlock(rec) {
		endVal := infoEnd(rec)
		if endVal > v.prevPos {
			v.prevPos = endVal
		}
		v.prevEnd = endVal
	} else {
		if rec.Pos > v.prevPos {
			v.prevPos = rec.Pos
		}
		v.prevEnd = 0
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

	// PEDIGREE Father / Mother / Original references must
	// resolve to declared samples (or to other PEDIGREE
	// IDs). §"PedigreeInDetail" leaves the exact
	// semantics open; we enforce the simpler "names
	// must be known" rule in ModeStrict.
	if v.Mode == ModeStrict && len(h.SampleNames) > 0 {
		known := make(map[string]struct{}, len(h.SampleNames)+len(h.Pedigrees))
		for _, n := range h.SampleNames {
			known[n] = struct{}{}
		}
		for _, p := range h.Pedigrees {
			known[p.ID] = struct{}{}
		}
		for _, p := range h.Pedigrees {
			for _, ref := range []string{p.Original, p.Father, p.Mother} {
				if ref == "" {
					continue
				}
				if _, ok := known[ref]; !ok {
					issues = append(issues, ValidationError{
						Severity: SeverityError,
						Field:    "header/PEDIGREE",
						Msg:      fmt.Sprintf("PEDIGREE ID %q references unknown sample %q", p.ID, ref),
					})
				}
			}
		}
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

	// gVCF non-variant block rules (§1.7): <*>
	// requires END ≥ POS, and the block must not
	// overlap the previous record. Cross-record
	// overlap is checked here (in the per-record
	// path) using prevPos / prevEnd.
	issues = append(issues, v.checkGVCFBlock(rec, idx)...)
	if v.Reference != nil {
		issues = append(issues, v.checkGVCFBlockRef(rec, idx, line)...)
	}

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
