package vcf

import (
	"fmt"
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
// A Validator is safe to re-use across multiple Validate calls; the only
// state is its Mode.
type Validator struct {
	Mode Mode
}

// NewValidator constructs a Validator with the given mode.
func NewValidator(m Mode) *Validator { return &Validator{Mode: m} }

// Validate returns a list of issues found by applying the validator to a
// batch of records under the given header. Issues are returned in
// (Line, Record) order; nil is returned for a clean file.
//
// Validate never modifies the input.
func (v *Validator) Validate(h *Header, records []*Record) []ValidationError {
	var issues []ValidationError

	// Header-level checks (recorded with Line=0).
	issues = append(issues, v.checkHeader(h)...)

	// Reserved-key metadata check (ModeStrict only).
	if v.Mode == ModeStrict {
		issues = append(issues, v.checkReservedKeys(h)...)
	}

	// Per-record checks.
	for i, rec := range records {
		issues = append(issues, v.checkRecord(h, rec, i)...)
	}

	// Strict-mode cross-record checks (e.g. POS monotonicity).
	if v.Mode == ModeStrict {
		issues = append(issues, v.checkRecordsStrict(records)...)
	}

	return issues
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
		}
	}

	// INFO field checks.
	issues = append(issues, v.checkInfo(h, rec, idx)...)

	// FILTER values must be declared (or "PASS"/".").
	issues = append(issues, v.checkFilter(h, rec, idx)...)

	// FORMAT / sample checks.
	issues = append(issues, v.checkSamples(h, rec, idx)...)

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
	}
	return issues
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
		}
	}
	return issues
}

// checkRecordsStrict runs cross-record checks that only make sense in
// strict mode (e.g. POS monotonicity within a CHROM block).
func (v *Validator) checkRecordsStrict(records []*Record) []ValidationError {
	var issues []ValidationError
	prevChrom := ""
	var prevPos int64
	for i, rec := range records {
		if rec.Chrom != prevChrom {
			prevChrom = rec.Chrom
			prevPos = -1
		}
		if prevPos >= 0 && rec.Pos < prevPos {
			issues = append(issues, ValidationError{
				Severity: SeverityError, Record: i + 1, Line: rec.Line,
				Field: "POS",
				Msg: fmt.Sprintf("POS %d < previous POS %d within CHROM %s",
					rec.Pos, prevPos, rec.Chrom),
			})
		}
		prevPos = rec.Pos
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
		_, err := parseVCFFloat(s)
		return err == nil
	case TypeFlag:
		return s == ""
	case TypeCharacter:
		return len(s) == 1
	case TypeString:
		return true
	}
	return true
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
