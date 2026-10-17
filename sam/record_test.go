package sam

import "testing"

func TestParseTagType(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in      byte
		want    TagType
		wantErr bool
	}{
		{'A', TagA, false},
		{'C', TagC, false},
		{'f', Tagf, false},
		{'H', TagH, false},
		{'i', Tagi, false},
		{'Z', TagZ, false},
		{'B', TagB, false},
		{'X', 0, true},
		{'a', 0, true}, // case sensitive
		{0, 0, true},
	}
	for _, c := range cases {
		got, err := ParseTagType(c.in)
		if c.wantErr {
			if err == nil {
				t.Errorf("ParseTagType(%q) = nil err; want non-nil", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseTagType(%q) unexpected err: %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseTagType(%q) = %q; want %q", c.in, got, c.want)
		}
		if got.String() != string(c.in) {
			t.Errorf("TagType(%q).String() = %q; want %q", c.in, got.String(), string(c.in))
		}
	}
}

func TestTagType_IsValid(t *testing.T) {
	t.Parallel()
	for _, tt := range []TagType{TagA, TagC, Tagf, TagH, Tagi, TagZ, TagB} {
		if !tt.IsValid() {
			t.Errorf("%q.IsValid = false; want true", tt)
		}
	}
	for _, tt := range []TagType{0, 'x', TagType(255)} {
		if tt.IsValid() {
			t.Errorf("TagType(%d).IsValid = true; want false", tt)
		}
	}
}

func TestParseError_Message(t *testing.T) {
	t.Parallel()

	t.Run("with field", func(t *testing.T) {
		t.Parallel()
		e := &ParseError{Line: 42, Field: "QNAME", Msg: "missing"}
		want := "line 42 QNAME: missing"
		if got := e.Error(); got != want {
			t.Errorf("Error() = %q; want %q", got, want)
		}
	})

	t.Run("without field", func(t *testing.T) {
		t.Parallel()
		e := &ParseError{Line: 7, Msg: "empty document"}
		want := "line 7: empty document"
		if got := e.Error(); got != want {
			t.Errorf("Error() = %q; want %q", got, want)
		}
	})
}

func TestIsParseError(t *testing.T) {
	t.Parallel()
	if !IsParseError(&ParseError{Line: 1, Msg: "x"}) {
		t.Error("IsParseError(*ParseError) = false; want true")
	}
	if !IsParseError(error(&ParseError{Line: 1, Msg: "x"})) {
		t.Error("IsParseError(error(*ParseError)) = false; want true")
	}
	if IsParseError(nil) {
		t.Error("IsParseError(nil) = true; want false")
	}
	plain := errorString("not a parse error")
	if IsParseError(plain) {
		t.Error("IsParseError(string error) = true; want false")
	}
}

type errorString string

func (e errorString) Error() string { return string(e) }
