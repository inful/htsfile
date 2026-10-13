package vcf

import (
	"errors"
	"testing"
)

func TestPercentEncode(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in, want string
	}{
		{"", ""},
		{"hello", "hello"},                   // no specials
		{"a:b", "a%3Ab"},                     // colon
		{"a;b", "a%3Bb"},                     // semicolon
		{"a=b", "a%3Db"},                     // equals
		{"a,b", "a%2Cb"},                     // comma
		{"a%b", "a%25b"},                     // percent
		{"a\rb", "a%0Db"},                    // CR
		{"a\nb", "a%0Ab"},                    // LF
		{"a\tb", "a%09b"},                    // TAB
		{"a;b=c,d\n", "a%3Bb%3Dc%2Cd%0A"},    // multiple
		{"no_change", "no_change"},           // safe
		{":::::::", "%3A%3A%3A%3A%3A%3A%3A"}, // all colons
		{"already %3A encoded", "already %253A encoded"}, // % becomes %25
	}
	for _, c := range cases {
		if got := percentEncode(c.in); got != c.want {
			t.Errorf("percentEncode(%q) = %q; want %q", c.in, got, c.want)
		}
	}
}

func TestPercentDecode(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in, want string
	}{
		{"", ""},
		{"hello", "hello"},
		{"a%3Ab", "a:b"},
		{"a%3Bb", "a;b"},
		{"a%3Db", "a=b"},
		{"a%2Cb", "a,b"},
		{"a%25b", "a%b"},
		{"a%0Db", "a\rb"},
		{"a%0Ab", "a\nb"},
		{"a%09b", "a\tb"},
		{"a%3Bb%3Dc%2Cd%0A", "a;b=c,d\n"},
		{"no%XXchange", "no%XXchange"}, // unknown sequence passed through
	}
	for _, c := range cases {
		if got := percentDecode(c.in); got != c.want {
			t.Errorf("percentDecode(%q) = %q; want %q", c.in, got, c.want)
		}
	}
}

func TestPercentRoundTrip(t *testing.T) {
	t.Parallel()
	// Anything that needs escaping round-trips losslessly.
	inputs := []string{
		"hello",
		"key=value;foo:bar,baz",
		"a\tb\nc\rd",
		"100%",
		";=,",
	}
	for _, s := range inputs {
		if got := percentDecode(percentEncode(s)); got != s {
			t.Errorf("round-trip(%q) = %q", s, got)
		}
	}
}

func TestParseError_Error(t *testing.T) {
	t.Parallel()
	e := &ParseError{Line: 12, Field: "INFO", Msg: "missing key"}
	// Just exercise the message; specifics are not asserted.
	if e.Error() == "" {
		t.Fatal("Error() returned empty string")
	}
	if e.Unwrap() != nil {
		t.Fatal("Unwrap should return nil")
	}
}

func TestValidationError_Error(t *testing.T) {
	t.Parallel()
	e := &ValidationError{Severity: SeverityError, Record: 3, Line: 42, Field: "INFO/DP", Msg: "wrong number of values"}
	if e.Error() == "" {
		t.Fatal("Error() returned empty string")
	}
}

func TestErrorWrappingHelpers(t *testing.T) {
	t.Parallel()
	pe := &ParseError{Line: 1, Msg: "x"}
	ve := &ValidationError{Line: 1, Msg: "y"}
	other := errors.New("other")

	if !IsParseError(pe) {
		t.Error("IsParseError(pe) = false")
	}
	if IsParseError(other) {
		t.Error("IsParseError(other) = true")
	}
	if !IsValidationError(ve) {
		t.Error("IsValidationError(ve) = false")
	}
	if IsValidationError(other) {
		t.Error("IsValidationError(other) = true")
	}
}

func TestSeverityString(t *testing.T) {
	t.Parallel()
	if got := SeverityError.String(); got != "error" {
		t.Errorf("SeverityError.String() = %q", got)
	}
	if got := SeverityWarning.String(); got != "warning" {
		t.Errorf("SeverityWarning.String() = %q", got)
	}
}
