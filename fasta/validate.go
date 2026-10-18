package fasta

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
	// ModeDefault checks structural well-formedness:
	// non-empty header ID, non-empty sequence. Most real-world
	// FASTA files pass.
	ModeDefault Mode = iota
	// ModeStrict adds alphabet validation: every sequence
	// byte must be one of the documented IUPAC characters
	// for nucleic acid or protein alphabets. The check is
	// permissive by default — extended IUPAC codes
	// (ambiguity codes like R, Y, N) are accepted — but
	// whitespace and control characters inside Sequence are
	// already stripped by the parser.
	ModeStrict
)

// DefaultAlphabet is the alphabet accepted by strict-mode
// validation: IUPAC nucleotide codes (uppercase and lowercase)
// plus the '-' gap character and '*' translation stop, which
// appear in protein sequences. Any byte not in this set
// surfaces as a strict-mode validation error.
const DefaultAlphabet = "ACGTURYSWKMBDHVNACGTURYSWKMBDHVN-*"

// Validator runs semantic checks against a parsed FASTA
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
	seenIDs := make(map[string]int)
	for _, r := range recs {
		issues = append(issues, v.checkOne(r)...)
		if r.ID != "" {
			seenIDs[r.ID]++
		}
	}
	// Strict-mode only: warn on duplicate IDs.
	if v.Mode == ModeStrict {
		for _, r := range recs {
			if r.ID != "" && seenIDs[r.ID] > 1 {
				issues = append(issues, Issue{
					Line:     r.Line,
					Field:    "ID",
					Severity: SeverityWarning,
					Msg:      fmt.Sprintf("duplicate ID %q (appears %d times)", r.ID, seenIDs[r.ID]),
				})
			}
		}
	}
	return issues
}

func (v *Validator) checkOne(r *Record) []Issue {
	issues := make([]Issue, 0, 4)
	if r.ID == "" {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "ID",
			Severity: SeverityError,
			Msg:      "header has no ID (empty token before first whitespace)",
		})
	}
	if len(r.Sequence) == 0 {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "sequence",
			Severity: SeverityError,
			Msg:      "record has no sequence data",
		})
	}
	if v.Mode == ModeStrict {
		alpha := v.Alphabet
		if alpha == "" {
			alpha = DefaultAlphabet
		}
		for i, c := range r.Sequence {
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
