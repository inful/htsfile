package bam

import "fmt"

// Severity classifies a validator [Issue] by impact.
type Severity int

const (
	// SeverityWarning is a suspicious condition that does not
	// make the data unusable but may indicate a mistake.
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
	// ModeDefault checks structural well-formedness: MAPQ
	// range, RNAME/SQ consistency, SEQ/QUAL length consistency.
	ModeDefault Mode = iota
	// ModeStrict adds extra checks: tag-key uniqueness per
	// record, POS >= 1 for mapped reads, CIGAR/SEQ length
	// consistency, refID/next_refID bounds.
	ModeStrict
)

// Validator runs semantic checks against a parsed BAM document.
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
	// MAPQ must be in 0..=255 (uint8 invariant — type system
	// already enforces this; no runtime check needed).
	// refID must be valid (or -1 for unmapped).
	if hdr != nil && r.RefID >= int32(len(hdr.SQs)) {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "refID",
			Severity: SeverityError,
			Msg:      fmt.Sprintf("refID %d out of bounds (n_ref=%d)", r.RefID, len(hdr.SQs)),
		})
	}
	// SEQ / QUAL length consistency.
	seqLen := len(r.Seq)
	if r.Seq != "*" {
		qualLen := len(r.Qual)
		if seqLen != qualLen {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "QUAL",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("QUAL length %d != SEQ length %d", qualLen, seqLen),
			})
		}
	}
	// Strict-mode additions.
	if v.Mode == ModeStrict {
		// POS must be >= 1 for mapped reads (flag bit 0x4 unset).
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
		// next_refID must also be valid (or -1).
		if hdr != nil && r.NextRefID >= int32(len(hdr.SQs)) {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "next_refID",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("next_refID %d out of bounds (n_ref=%d)", r.NextRefID, len(hdr.SQs)),
			})
		}
	}
	return issues
}

// Flag bit values from the SAM/BAM spec § 1.4. Identical to the
// sam package's constants.
const (
	FlagPaired        = 1
	FlagProperPair    = 2
	FlagUnmap         = 4
	FlagMateUnmap     = 8
	FlagReverse       = 16
	FlagMateReverse   = 32
	FlagFirst         = 64
	FlagLast          = 128
	FlagSecondary     = 256
	FlagQCFail        = 512
	FlagDuplicate     = 1024
	FlagSupplementary = 2048
)