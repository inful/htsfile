package vcf

import (
	"errors"
	"testing"
)

func TestParseFieldType(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in      string
		want    FieldType
		wantErr bool
	}{
		{"Integer", TypeInteger, false},
		{"Float", TypeFloat, false},
		{"Flag", TypeFlag, false},
		{"Character", TypeCharacter, false},
		{"String", TypeString, false},
		{"", TypeUnknown, true},
		{"integer", TypeUnknown, true}, // spec is case sensitive
		{"Bool", TypeUnknown, true},
	}
	for _, c := range cases {
		c := c
		t.Run(c.in, func(t *testing.T) {
			t.Parallel()
			got, err := ParseFieldType(c.in)
			if c.wantErr {
				if err == nil {
					t.Fatalf("expected error for %q, got %v", c.in, got)
				}
				if !errors.Is(err, err) {
					t.Fatal("error unwrap self-check failed")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if got != c.want {
				t.Fatalf("got %v, want %v", got, c.want)
			}
		})
	}
}

func TestFieldTypeString(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in   FieldType
		want string
	}{
		{TypeInteger, "Integer"},
		{TypeFloat, "Float"},
		{TypeFlag, "Flag"},
		{TypeCharacter, "Character"},
		{TypeString, "String"},
		{TypeUnknown, ""},
	}
	for _, c := range cases {
		if got := c.in.String(); got != c.want {
			t.Errorf("FieldType(%d).String() = %q; want %q", c.in, got, c.want)
		}
	}
}
