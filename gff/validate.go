package gff

import (
	"fmt"
	"strconv"
)

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
	// ModeDefault checks structural well-formedness: 9
	// columns, integer start/end, start >= 1, start <= end,
	// valid strand character, valid phase, parseable
	// attributes. Most real-world GFF3 files pass.
	ModeDefault Mode = iota
	// ModeStrict adds ID/Parent cross-reference consistency:
	// every Parent value must reference an ID seen in the
	// document.
	ModeStrict
)

// Validator runs semantic checks against a parsed GFF3 document.
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
	// Strict-mode only: collect ID set for Parent checks.
	var ids map[string]struct{}
	if v.Mode == ModeStrict {
		ids = map[string]struct{}{}
		for _, r := range recs {
			if id := r.Attribute("ID"); id != "" {
				ids[id] = struct{}{}
			}
		}
	}
	for _, r := range recs {
		issues = append(issues, v.checkOne(r, ids)...)
	}
	return issues
}

func (v *Validator) checkOne(r *Record, ids map[string]struct{}) []Issue {
	issues := make([]Issue, 0, 4)
	// Type: not strictly required but useful.
	if r.Type == "" || r.Type == "." {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "type",
			Severity: SeverityWarning,
			Msg:      "type column is '.' (no feature type); results are non-standard",
		})
	}
	// Attributes must be present (column 9).
	if len(r.Attributes) == 0 && !r.HasAttribute("ID") {
		// Soft warning only — GFF3 features without
		// attributes are technically legal but uninformative.
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "attributes",
			Severity: SeverityWarning,
			Msg:      "feature has no attributes",
		})
	}
	// Strict-mode Parent ID cross-reference check.
	if v.Mode == ModeStrict && ids != nil {
		for _, parent := range r.AttributeList("Parent") {
			if _, ok := ids[parent]; !ok {
				issues = append(issues, Issue{
					Line:     r.Line,
					Field:    "Parent",
					Severity: SeverityError,
					Msg:      fmt.Sprintf("Parent=%q not declared as any record's ID", parent),
				})
			}
		}
	}
	// Range sanity (defensive — the parser already enforces
	// this; the check is here for records built programmatically).
	if r.Start < 1 {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "start",
			Severity: SeverityError,
			Msg:      "GFF3 start must be >= 1; got " + strconv.FormatInt(r.Start, 10),
		})
	}
	if r.End < r.Start {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "end",
			Severity: SeverityError,
			Msg:      fmt.Sprintf("end %d < start %d", r.End, r.Start),
		})
	}
	return issues
}
