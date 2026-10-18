package gff

import (
	"bytes"
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
	eWithField := &ParseError{Line: 7, Field: "feature", Msg: "9 columns required"}
	if msg := eWithField.Error(); msg != "line 7 feature: 9 columns required" {
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
// Coverage pin — Record.HasAttribute and AttributeList
// ---------------------------------------------------------------------------

func TestRecord_HasAttribute(t *testing.T) {
	rec := &Record{
		Attributes: map[string][]string{
			"ID":    {"gene1"},
			"Alias": {"ABC1", "ABC2"},
		},
	}
	if !rec.HasAttribute("ID") {
		t.Error("HasAttribute(ID) = false, want true")
	}
	if !rec.HasAttribute("Alias") {
		t.Error("HasAttribute(Alias) = false, want true")
	}
	if rec.HasAttribute("Missing") {
		t.Error("HasAttribute(Missing) = true, want false")
	}
	// Nil receiver.
	if (*Record)(nil).HasAttribute("ID") {
		t.Error("nil.HasAttribute(ID) = true, want false")
	}
}

func TestRecord_AttributeList(t *testing.T) {
	rec := &Record{
		Attributes: map[string][]string{
			"ID":    {"gene1"},
			"Alias": {"ABC1", "ABC2"},
		},
	}
	if got := rec.AttributeList("ID"); len(got) != 1 || got[0] != "gene1" {
		t.Errorf("AttributeList(ID): got %v", got)
	}
	if got := rec.AttributeList("Alias"); len(got) != 2 {
		t.Errorf("AttributeList(Alias): got %v want 2 entries", got)
	}
	if got := rec.AttributeList("Missing"); got != nil {
		t.Errorf("AttributeList(Missing): got %v want nil", got)
	}
	// Nil receiver.
	if got := (*Record)(nil).AttributeList("ID"); got != nil {
		t.Errorf("nil.AttributeList: got %v want nil", got)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — NewEncoderCRLF
// ---------------------------------------------------------------------------

func TestNewEncoderCRLF(t *testing.T) {
	var buf bytes.Buffer
	e := NewEncoderCRLF(&buf)
	if e == nil {
		t.Fatal("NewEncoderCRLF returned nil")
	}
	if string(e.lineEnding) != "\r\n" {
		t.Errorf("lineEnding: got %q want %q", e.lineEnding, "\r\n")
	}
}
