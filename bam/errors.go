// Package bam provides a streaming parser and validator for the
// BAM (Binary Alignment/Map) format used to wrap SAM records
// inside a BGZF-compressed binary container.
//
// The package-level documentation is in [doc.go].
package bam

import (
	"errors"
	"fmt"
)

// ParseError describes a failure to parse a BAM stream.
//
// The Line field is set to the 1-based offset of the offending
// record within the alignment section (the first record is line 1).
// For header-level errors Line is 0. Field identifies which
// structural element failed ("magic", "header", "ref", "record",
// "tag", "cigar", "bin", etc.).
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
