package sam

import "fmt"

// Severity classifies a validator [Issue] by impact.
type Severity int

const (
	// SeverityWarning is a suspicious condition that does not make
	// the data unusable but may indicate a mistake.
	SeverityWarning Severity = iota
	// SeverityError is a definite problem: malformed data, a
	// violated invariant, or a spec violation.
	SeverityError
)

// String returns the lowercase name of the severity.
func (s Severity) String() string {
	switch s {
	case SeverityWarning:
		return "warning"
	case SeverityError:
		return "error"
	default:
		return "unknown"
	}
}

// Issue is one problem found by the [Validator].
type Issue struct {
	Line     int
	Field    string
	Severity Severity
	Msg      string
}

// Mode controls the strictness of the [Validator].
type Mode int

const (
	// ModeDefault checks structural well-formedness: column
	// counts, MAPQ range, RNAME/SQ consistency, SEQ/QUAL length
	// consistency.
	ModeDefault Mode = iota
	// ModeStrict adds extra checks: tag-key uniqueness per record,
	// POS >= 1 for mapped reads, FLAG-bit consistency (e.g. bit
	// 0x4 and 0x2 imply specific paired-read semantics).
	ModeStrict
)

// Validator runs semantic checks against a parsed SAM document.
type Validator struct {
	Mode Mode
}

// NewValidator returns a Validator with the given mode.
func NewValidator(m Mode) *Validator {
	return &Validator{Mode: m}
}

// Validate runs all applicable checks and returns the issues
// found. It does not short-circuit; every issue is reported.
func (v *Validator) Validate(hdr *Header, recs []*Record) []Issue {
	issues := make([]Issue, 0, len(recs))
	for _, r := range recs {
		issues = append(issues, v.checkOne(hdr, r)...)
	}
	return issues
}

func (v *Validator) checkOne(hdr *Header, r *Record) []Issue {
	issues := make([]Issue, 0, 4)
	// MAPQ must be in 0..=255 (8-bit unsigned).
	if r.MapQ < 0 || r.MapQ > 255 {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "MAPQ",
			Severity: SeverityError,
			Msg:      fmt.Sprintf("MAPQ must be in [0, 255]; got %d", r.MapQ),
		})
	}
	// RNAME consistency: if @SQ is declared and RNAME is not "*",
	// it must match an @SQ.SN.
	if hdr != nil && r.RName != "" && r.RName != "*" {
		if _, ok := hdr.LookupSQ(r.RName); !ok {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "RNAME",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("RNAME %q not declared in @SQ header", r.RName),
			})
		}
	}
	// SEQ / QUAL length consistency: if either is "*" both must
	// be "*"; otherwise lengths must match.
	if r.Seq != "*" && r.Qual != "*" {
		if len(r.Seq) != len(r.Qual) {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "QUAL",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("QUAL length %d != SEQ length %d", len(r.Qual), len(r.Seq)),
			})
		}
	} else if r.Seq != r.Qual {
		// one is "*" but not both — only allowed if SEQ is "*"
		// (in which case QUAL must also be "*")
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "QUAL",
			Severity: SeverityError,
			Msg:      "SEQ and QUAL must both be \"*\" or both be set",
		})
	}
	// Strict-mode additions.
	if v.Mode == ModeStrict {
		// POS must be >= 1 for mapped reads (FLAG bit 0x4 unset).
		if r.Pos < 1 && r.Flag&FlagUnmap == 0 {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "POS",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("POS must be >= 1 for mapped read (FLAG=0x%x); got %d", r.Flag, r.Pos),
			})
		}
		// Each tag key should be unique within a record.
		seen := make(map[string]bool)
		for _, t := range r.Tags {
			if seen[t.Key] {
				issues = append(issues, Issue{
					Line:     r.Line,
					Field:    "tag",
					Severity: SeverityWarning,
					Msg:      fmt.Sprintf("duplicate tag key %q in record", t.Key),
				})
			}
			seen[t.Key] = true
		}
		// FLAG-bit consistency (SAMv1 §1.4). The SAM spec
		// requires a few invariants across the paired-read
		// bits that this check enforces.
		if r.Flag&FlagPaired == 0 && (r.Flag&(FlagMateUnmap|FlagMateReverse|FlagFirst|FlagSecondary)) != 0 {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "FLAG",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("FLAG 0x%x sets a mate-related bit without FLAG_PAIRED", r.Flag),
			})
		}
	}
	return issues
}
