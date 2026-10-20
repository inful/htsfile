package fastq

import (
	"errors"
	"reflect"
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
	eWithField := &ParseError{Line: 7, Field: "header", Msg: "must start with @"}
	if msg := eWithField.Error(); msg != "line 7 header: must start with @" {
		t.Errorf("Error (with field): got %q", msg)
	}
	if e.Unwrap() != nil {
		t.Errorf("Unwrap: got %v want nil", e.Unwrap())
	}
	// IsParseError: nil, plain error, *ParseError.
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
// Coverage pin — Decoder.PhredOffset accessor and
// qualRange empty-input branch
// ---------------------------------------------------------------------------

func TestDecoder_PhredOffsetAccessor(t *testing.T) {
	// Default before any record is decoded.
	dec := NewDecoder(nil) //nolint:staticcheck // nil reader; we only check the config
	if got := dec.PhredOffset(); got != 33 {
		t.Errorf("default PhredOffset: got %d want 33", got)
	}
	dec.SetPhredOffset(64)
	if got := dec.PhredOffset(); got != 64 {
		t.Errorf("after SetPhredOffset: got %d want 64", got)
	}
}

func TestQualRange_EmptyInput(t *testing.T) {
	mn, mx := qualRange(nil)
	if mn != 0 || mx != 0 {
		t.Errorf("empty input: got (%d, %d) want (0, 0)", mn, mx)
	}
	mn, mx = qualRange([]byte{})
	if mn != 0 || mx != 0 {
		t.Errorf("empty slice: got (%d, %d) want (0, 0)", mn, mx)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Record.Length (nil and populated receivers)
// ---------------------------------------------------------------------------

func TestRecord_Length(t *testing.T) {
	if got := (*Record)(nil).Length(); got != 0 {
		t.Errorf("nil Length: got %d want 0", got)
	}
	rec := &Record{Seq: []byte("ACGT")}
	if got := rec.Length(); got != 4 {
		t.Errorf("Length: got %d want 4", got)
	}
	// Make sure Length uses Seq, not Qual.
	rec.Seq = []byte("AC")
	rec.Qual = []byte("!!!!!")
	if got := rec.Length(); got != 2 {
		t.Errorf("Length uses Seq: got %d want 2", got)
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — Validate ModeStrict alphabet-on-empty-seq path
// ---------------------------------------------------------------------------

func TestValidator_Strict_EmptySeqStillChecked(t *testing.T) {
	// Strict mode on a record with empty Seq: no alphabet
	// errors are reported (nothing to check), but the
	// empty-seq error from ModeDefault should still surface.
	v := NewValidator(ModeStrict)
	rec := &Record{Header: "r", Seq: []byte{}, Qual: []byte{}, PhredOffset: 33}
	issues := v.Validate([]*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("got %d issues, want 1: %v", len(issues), issues)
	}
	if issues[0].Field != "sequence" {
		t.Errorf("issue: %+v", issues[0])
	}
}

// ---------------------------------------------------------------------------
// Coverage pin — NewDecoderWithMaxLine max<=0 fallback
// ---------------------------------------------------------------------------

func TestNewDecoderWithMaxLine_ZeroFallsBackToDefault(t *testing.T) {
	dec := NewDecoderWithMaxLine(nil, 0) //nolint:staticcheck
	if dec.maxLineLength != DefaultMaxLineLength {
		t.Errorf("max=0: got %d want %d", dec.maxLineLength, DefaultMaxLineLength)
	}
	dec = NewDecoderWithMaxLine(nil, -5) //nolint:staticcheck
	if dec.maxLineLength != DefaultMaxLineLength {
		t.Errorf("max=-5: got %d want %d", dec.maxLineLength, DefaultMaxLineLength)
	}
}

// _ keeps reflect imports reflectively hidden if the only call site goes away.
var _ = reflect.DeepEqual
