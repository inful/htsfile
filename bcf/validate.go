package bcf

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
	// ModeDefault checks structural well-formedness: REF not
	// empty, ALT not equal to REF.
	ModeDefault Mode = iota
	// ModeStrict additionally enforces header-declared
	// constraints (e.g. an INFO key seen on a record must be
	// declared in ##INFO; declared Type matches the parsed
	// value's coerced type).
	ModeStrict
)

// Validator runs semantic checks against a parsed BCF
// document.
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
	if r.Ref == "" {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "REF",
			Severity: SeverityError,
			Msg:      "REF is empty",
		})
	}
	for i, alt := range r.Alt {
		if alt == r.Ref {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "ALT",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("ALT[%d] %q equals REF", i, alt),
			})
		}
	}
	if v.Mode == ModeStrict && hdr != nil {
		for k := range r.Info {
			if _, ok := hdr.Infos[k]; !ok {
				issues = append(issues, Issue{
					Line:     r.Line,
					Field:    "INFO",
					Severity: SeverityError,
					Msg:      fmt.Sprintf("INFO key %q is not declared in the header", k),
				})
			}
		}
	}
	return issues
}