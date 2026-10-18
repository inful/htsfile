package fasta

import (
	"errors"
	"testing"
)

// ---------------------------------------------------------------------------
// Coverage pin — ParseError helpers and Severity
// ---------------------------------------------------------------------------

func TestParseError_ErrorAndHelpers(t *testing.T) {
	e := &ParseError{Line: 5, Msg: "missing field"}
	if msg := e.Error(); msg != "line 5: missing field" {
		t.Errorf("Error: got %q", msg)
	}
	eWithField := &ParseError{Line: 7, Field: "header", Msg: "must start with >"}
	if msg := eWithField.Error(); msg != "line 7 header: must start with >" {
		t.Errorf("Error (with field): got %q", msg)
	}
	if e.Unwrap() != nil {
		t.Errorf("Unwrap: got %v want nil", e.Unwrap())
	}
	if IsParseError(nil) {
		t.Error("IsParseError(nil) = true")
	}
	if IsParseError(errors.New("plain")) {
		t.Error("IsParseError(plain) = true")
	}
	if !IsParseError(&ParseError{}) {
		t.Error("IsParseError(*ParseError) = false")
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Severity.String
// ---------------------------------------------------------------------------

func TestSeverity_String(t *testing.T) {
	cases := []struct {
		sev  Severity
		want string
	}{
		{SeverityWarning, "warning"},
		{SeverityError, "error"},
		{Severity(99), "unknown"},
		{Severity(0), "warning"}, // zero value is SeverityWarning
	}
	for _, tc := range cases {
		if got := tc.sev.String(); got != tc.want {
			t.Errorf("String(%v): got %q want %q", tc.sev, got, tc.want)
		}
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Record.Length (nil and populated receivers)
// ---------------------------------------------------------------------------

func TestRecord_Length(t *testing.T) {
	if got := (*Record)(nil).Length(); got != 0 {
		t.Errorf("nil Length: got %d want 0", got)
	}
	rec := &Record{Sequence: []byte("ACGT")}
	if got := rec.Length(); got != 4 {
		t.Errorf("Length: got %d want 4", got)
	}
	// Empty sequence.
	rec = &Record{Sequence: []byte{}}
	if got := rec.Length(); got != 0 {
		t.Errorf("empty Length: got %d want 0", got)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — NewEncoderCRLF
// ---------------------------------------------------------------------------

func TestNewEncoderCRLF(t *testing.T) {
	// Just verify the constructor returns a usable encoder with
	// the expected CR+LF line ending.
	e := NewEncoderCRLF(&errWriter{})
	if e == nil {
		t.Fatal("NewEncoderCRLF returned nil")
	}
	if string(e.lineEnding) != "\r\n" {
		t.Errorf("lineEnding: got %q want %q", e.lineEnding, "\r\n")
	}
	if e.lineWidth != DefaultLineWidth {
		t.Errorf("lineWidth: got %d want %d", e.lineWidth, DefaultLineWidth)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — stripWhitespace edge cases (embedded tab / CR / multiple spaces)
// ---------------------------------------------------------------------------

func TestStripWhitespace_Embedded(t *testing.T) {
	cases := []struct {
		in, want string
	}{
		{"AC GT", "ACGT"},
		{"AC\tGT", "ACGT"},
		{"A\rB\nC", "ABC"},
		{"A   B   C", "ABC"},
		{"  ABC  ", "ABC"},
		{"ABC", "ABC"}, // no whitespace
	}
	for _, tc := range cases {
		got := stripWhitespace([]byte(tc.in))
		if string(got) != tc.want {
			t.Errorf("stripWhitespace(%q): got %q want %q", tc.in, got, tc.want)
		}
	}
}

// errWriter is an io.Writer that always returns an error. Used to
// exercise Encoder's write-error paths if needed in future tests.
type errWriter struct{}

func (errWriter) Write(p []byte) (int, error) { return 0, errors.New("synthetic") }
