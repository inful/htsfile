package bed

import (
	"errors"
	"fmt"
)

// ParseError describes a failure to parse a BED line. The Line field
// is 1-based; Field identifies which column (or "track" / "browser"
// for header lines) caused the error; Msg is a human-readable
// explanation.
type ParseError struct {
	Line  int
	Field string
	Msg   string
}

// Error implements the error interface.
func (e *ParseError) Error() string {
	if e.Field == "" {
		return fmt.Sprintf("line %d: %s", e.Line, e.Msg)
	}
	return fmt.Sprintf("line %d %s: %s", e.Line, e.Field, e.Msg)
}

// Unwrap returns nil: ParseError is a leaf error.
func (e *ParseError) Unwrap() error { return nil }

// IsParseError reports whether err is or wraps a *ParseError.
func IsParseError(err error) bool {
	if err == nil {
		return false
	}
	var pe *ParseError
	return errors.As(err, &pe)
}
