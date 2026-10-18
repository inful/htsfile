package fastq

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
	// ModeDefault checks structural well-formedness: non-empty
	// header, non-empty sequence, sequence and quality same
	// length, every quality byte >= the record's PhredOffset.
	// Most real-world FASTQ files pass.
	ModeDefault Mode = iota
	// ModeStrict adds alphabet validation: every sequence
	// byte must be one of the documented IUPAC characters
	// (or the configured custom alphabet). The check is
	// case-insensitive for nucleic-acid sequences; the
	// default alphabet covers both upper- and lower-case.
	ModeStrict
)

// DefaultAlphabet is the alphabet accepted by strict-mode
// validation: IUPAC nucleotide codes (uppercase and lowercase)
// plus the '-' gap character, which appears in some
// protein-style FASTQ files. Any byte not in this set
// surfaces as a strict-mode validation error.
const DefaultAlphabet = "ACGTURYSWKMBDHVNacgturyswkmbdhvn-"

// Validator runs semantic checks against a parsed FASTQ
// document.
type Validator struct {
	// Mode controls strictness (see [Mode] constants).
	Mode Mode

	// Alphabet overrides the default allowed byte set in
	// strict mode. Empty means use [DefaultAlphabet].
	Alphabet string
}

// NewValidator returns a Validator with the given mode.
func NewValidator(m Mode) *Validator {
	return &Validator{Mode: m, Alphabet: DefaultAlphabet}
}

// Validate runs all applicable checks and returns the issues
// found. It does not short-circuit; every issue is reported.
func (v *Validator) Validate(recs []*Record) []Issue {
	issues := make([]Issue, 0, len(recs))
	for _, r := range recs {
		issues = append(issues, v.checkOne(r)...)
	}
	return issues
}

func (v *Validator) checkOne(r *Record) []Issue {
	issues := make([]Issue, 0, 4)
	if r.Header == "" {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "header",
			Severity: SeverityError,
			Msg:      "FASTQ header is empty",
		})
	}
	if len(r.Seq) == 0 {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "sequence",
			Severity: SeverityError,
			Msg:      "record has no sequence data",
		})
	}
	if len(r.Qual) != len(r.Seq) {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "quality",
			Severity: SeverityError,
			Msg: fmt.Sprintf("quality length %d != sequence length %d",
				len(r.Qual), len(r.Seq)),
		})
	}
	for i, q := range r.Qual {
		if q < r.PhredOffset {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "quality",
				Severity: SeverityError,
				Msg: fmt.Sprintf("quality byte %d at index %d is below PhredOffset %d",
					q, i, r.PhredOffset),
			})
			break // one qual-below-offset issue per record is enough
		}
	}
	if v.Mode == ModeStrict {
		alpha := v.Alphabet
		if alpha == "" {
			alpha = DefaultAlphabet
		}
		for i, c := range r.Seq {
			if !isInAlphabet(c, alpha) {
				issues = append(issues, Issue{
					Line:     r.Line,
					Field:    "sequence",
					Severity: SeverityError,
					Msg:      fmt.Sprintf("byte %q at position %d is not in the configured alphabet", c, i),
				})
				break // one alphabet error per record is enough
			}
		}
	}
	return issues
}

func isInAlphabet(c byte, alphabet string) bool {
	// Small alphabets (<=128 entries) — linear scan is fine
	// and avoids the allocations of a map.
	for i := 0; i < len(alphabet); i++ {
		if alphabet[i] == c {
			return true
		}
	}
	return false
}
