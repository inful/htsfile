package vcf

import (
	"errors"
	"fmt"
	"strings"
)

// ParseError describes a problem encountered while reading a VCF document.
//
// Line is the 1-based line number in the input stream; Field, when set,
// identifies the column or meta-line token the error concerns.
type ParseError struct {
	Line  int
	Field string
	Msg   string
}

// Error implements [error].
func (e *ParseError) Error() string {
	var b strings.Builder
	fmt.Fprintf(&b, "vcf: parse error")
	if e.Line > 0 {
		fmt.Fprintf(&b, " at line %d", e.Line)
	}
	if e.Field != "" {
		fmt.Fprintf(&b, " (%s)", e.Field)
	}
	b.WriteString(": ")
	b.WriteString(e.Msg)
	return b.String()
}

// Unwrap allows [errors.Is]/[errors.As] on the underlying message.
func (e *ParseError) Unwrap() error { return nil }

// ValidationError describes a semantic problem detected by [Validator].
//
// Severity indicates whether the issue is an Error (the file is invalid) or
// a Warning (the file is technically valid but unusual).
//
// The same fields as errors already exist: Line, Record, Field, Msg.
// Line and Record are 1-based; Line is the absolute file line and Record is
// the zero-based record index.
type ValidationError struct {
	Severity Severity
	Line     int
	Record   int
	Field    string
	Msg      string
}

// Error implements [error].
func (e *ValidationError) Error() string {
	var b strings.Builder
	fmt.Fprintf(&b, "vcf: %s", e.Severity)
	if e.Record > 0 {
		fmt.Fprintf(&b, " record %d", e.Record)
	}
	if e.Line > 0 {
		fmt.Fprintf(&b, " (line %d)", e.Line)
	}
	if e.Field != "" {
		fmt.Fprintf(&b, " %s", e.Field)
	}
	b.WriteString(": ")
	b.WriteString(e.Msg)
	return b.String()
}

// Severity classifies a [ValidationError].
type Severity int8

const (
	// SeverityError means the file violates a hard requirement of the spec.
	// The file is invalid.
	SeverityError Severity = iota
	// SeverityWarning means the file is technically valid but exhibits an
	// unusual or discouraged pattern (e.g. a sample missing from a column
	// when others are present).
	SeverityWarning
)

// String returns the spec text representation of the severity.
func (s Severity) String() string {
	switch s {
	case SeverityError:
		return "error"
	case SeverityWarning:
		return "warning"
	}
	return "unknown"
}

// IsParseError reports whether err is or wraps a [*ParseError].
func IsParseError(err error) bool {
	var pe *ParseError
	return errors.As(err, &pe)
}

// IsValidationError reports whether err is or wraps a [*ValidationError].
func IsValidationError(err error) bool {
	var ve *ValidationError
	return errors.As(err, &ve)
}

// percentEncode escapes the bytes that must be percent-encoded when
// appearing in the textual representation of a String / Character value.
// Per VCFv4.5 §1.2 the encoded forms are:
//
//	%3A  :  (colon)
//	%3B  ;  (semicolon)
//	%3D  =  (equal sign)
//	%25  %  (percent sign)
//	%2C  ,  (comma)
//	%0D  CR
//	%0A  LF
//	%09  TAB
func percentEncode(s string) string {
	if !needsEscape(s) {
		return s
	}
	var b strings.Builder
	b.Grow(len(s) + 4)
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch c {
		case ':':
			b.WriteString("%3A")
		case ';':
			b.WriteString("%3B")
		case '=':
			b.WriteString("%3D")
		case '%':
			b.WriteString("%25")
		case ',':
			b.WriteString("%2C")
		case '\r':
			b.WriteString("%0D")
		case '\n':
			b.WriteString("%0A")
		case '\t':
			b.WriteString("%09")
		default:
			b.WriteByte(c)
		}
	}
	return b.String()
}

func needsEscape(s string) bool {
	return strings.ContainsAny(s, ":;=%,\r\n\t")
}

// percentDecode reverses [percentEncode]. Unknown % sequences are passed
// through unchanged so that downstream code can detect them.
func percentDecode(s string) string {
	if !strings.ContainsRune(s, '%') {
		return s
	}
	var b strings.Builder
	b.Grow(len(s))
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c != '%' || i+2 >= len(s) {
			b.WriteByte(c)
			continue
		}
		hh := s[i+1 : i+3]
		var v byte
		_, scanErr := fmt.Sscanf(hh, "%02X", &v)
		if scanErr != nil {
			b.WriteByte(c)
			continue
		}
		switch hh {
		case "3A", "3a":
			b.WriteByte(':')
		case "3B", "3b":
			b.WriteByte(';')
		case "3D", "3d":
			b.WriteByte('=')
		case "25":
			b.WriteByte('%')
		case "2C", "2c":
			b.WriteByte(',')
		case "0D", "0d":
			b.WriteByte('\r')
		case "0A", "0a":
			b.WriteByte('\n')
		case "09":
			b.WriteByte('\t')
		default:
			b.WriteByte(c)
		}
		i += 2
	}
	return b.String()
}
