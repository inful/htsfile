package vcf

import (
	"bytes"
	"math"
	"slices"
	"testing"
)

func TestInfoValue_Accessors(t *testing.T) {
	t.Parallel()
	t.Run("String", func(t *testing.T) {
		t.Parallel()
		v := InfoValue{typedCore: typedCore{Type: TypeString, Strs: []string{"hello", "world"}}}
		if got := v.String(); got != "hello" {
			t.Errorf("String() = %q; want %q", got, "hello")
		}
		zero := InfoValue{}
		if zero.String() != "" {
			t.Error("zero-value should yield empty string")
		}
	})
	t.Run("Int", func(t *testing.T) {
		t.Parallel()
		v := InfoValue{typedCore: typedCore{Type: TypeInteger, Ints: []int64{42, 7}}}
		if got := v.Int(); got != 42 {
			t.Errorf("Int() = %d; want 42", got)
		}
		zero := InfoValue{}
		if zero.Int() != 0 {
			t.Error("zero-value should yield 0")
		}
	})
	t.Run("Float", func(t *testing.T) {
		t.Parallel()
		v := InfoValue{typedCore: typedCore{Type: TypeFloat, Floats: []float64{3.14, 2.72}}}
		if got := v.Float(); got != 3.14 {
			t.Errorf("Float() = %v; want 3.14", got)
		}
		zero := InfoValue{}
		if zero.Float() != 0 {
			t.Error("zero-value should yield 0")
		}
	})
	t.Run("HasValue", func(t *testing.T) {
		t.Parallel()
		v := InfoValue{typedCore: typedCore{Missing: true}}
		if v.HasValue() {
			t.Error("missing should not have value")
		}
		v = InfoValue{typedCore: typedCore{Type: TypeFlag}}
		if v.HasValue() {
			t.Error("flag not set should not have value")
		}
		v = InfoValue{typedCore: typedCore{Type: TypeFlag}, Flag: true}
		if !v.HasValue() {
			t.Error("flag set should have value")
		}
	})
}

func TestFieldValue_Accessors(t *testing.T) {
	t.Parallel()
	t.Run("String", func(t *testing.T) {
		t.Parallel()
		v := FieldValue{typedCore: typedCore{Type: TypeString, Strs: []string{"0/1"}}}
		if got := v.String(); got != "0/1" {
			t.Errorf("String() = %q", got)
		}
	})
	t.Run("Int", func(t *testing.T) {
		t.Parallel()
		v := FieldValue{typedCore: typedCore{Type: TypeInteger, Ints: []int64{99}}}
		if got := v.Int(); got != 99 {
			t.Errorf("Int() = %d", got)
		}
	})
	t.Run("Float", func(t *testing.T) {
		t.Parallel()
		v := FieldValue{typedCore: typedCore{Type: TypeFloat, Floats: []float64{0.95}}}
		if got := v.Float(); got != 0.95 {
			t.Errorf("Float() = %v", got)
		}
	})
	t.Run("HasValue", func(t *testing.T) {
		t.Parallel()
		v := FieldValue{typedCore: typedCore{Missing: true}}
		if v.HasValue() {
			t.Error("missing should not have value")
		}
		v = FieldValue{typedCore: typedCore{Type: TypeString, Strs: []string{"x"}}}
		if !v.HasValue() {
			t.Error("present should have value")
		}
	})
}

func TestNumberIsLocal(t *testing.T) {
	t.Parallel()
	cases := []struct {
		n    Number
		want bool
	}{
		{Number{Number: 1}, false},
		{Number{Fixed: SymbolA}, false},
		{Number{Fixed: SymbolLA}, true},
		{Number{Fixed: SymbolLR}, true},
		{Number{Fixed: SymbolLG}, true},
		{Number{Fixed: SymbolG}, false},
		{Number{Fixed: SymbolDot}, false},
	}
	for _, c := range cases {
		if got := c.n.IsLocal(); got != c.want {
			t.Errorf("Number%+v.IsLocal() = %v; want %v", c.n, got, c.want)
		}
	}
}

// TestDecodeTyped pins the behaviour of decodeTyped for every case the
// per-type entry points (decodeInfoValue, decodeFieldValue) depend on.
// This is the single source of truth for the Integer/Float/String/Missing
// coercion rules that previously lived as four near-identical type
// switches in parse.go.
func TestDecodeTyped(t *testing.T) {
	t.Parallel()

	t.Run("MissingDot", func(t *testing.T) {
		t.Parallel()
		// A single "." raw sets Missing and populates no typed slices,
		// regardless of Type.
		for _, ft := range []FieldType{TypeInteger, TypeFloat, TypeString, TypeCharacter, TypeFlag} {
			c := typedCore{Type: ft, Raw: []string{"."}}
			decodeTyped(&c, []string{"."})
			if !c.Missing {
				t.Errorf("%v: Missing = false; want true", ft)
			}
			if len(c.Ints) != 0 || len(c.Floats) != 0 || len(c.Strs) != 0 || len(c.Chars) != 0 {
				t.Errorf("%v: typed slices populated for \".\": %+v", ft, c)
			}
		}
	})

	t.Run("Integer", func(t *testing.T) {
		t.Parallel()
		c := typedCore{Type: TypeInteger, Raw: []string{"1", "2", "3"}}
		decodeTyped(&c, []string{"1", "2", "3"})
		if !slices.Equal(c.Ints, []int64{1, 2, 3}) {
			t.Errorf("Ints = %v; want [1 2 3]", c.Ints)
		}
	})

	t.Run("IntegerLenient", func(t *testing.T) {
		t.Parallel()
		// Non-numeric tokens become 0 and Raw is preserved for the
		// validator to surface.
		c := typedCore{Type: TypeInteger, Raw: []string{"7", "abc", "9"}}
		decodeTyped(&c, []string{"7", "abc", "9"})
		if !slices.Equal(c.Ints, []int64{7, 0, 9}) {
			t.Errorf("Ints = %v; want [7 0 9]", c.Ints)
		}
		if c.Missing {
			t.Error("Missing set unexpectedly on non-dot raw")
		}
	})

	t.Run("IntegerDotInList", func(t *testing.T) {
		t.Parallel()
		// A "." mid-list becomes 0 (per spec §1.4.5).
		c := typedCore{Type: TypeInteger, Raw: []string{"1", ".", "3"}}
		decodeTyped(&c, []string{"1", ".", "3"})
		if !slices.Equal(c.Ints, []int64{1, 0, 3}) {
			t.Errorf("Ints = %v; want [1 0 3]", c.Ints)
		}
	})

	t.Run("Float", func(t *testing.T) {
		t.Parallel()
		c := typedCore{Type: TypeFloat, Raw: []string{"1.5", "2.5"}}
		decodeTyped(&c, []string{"1.5", "2.5"})
		if !slices.Equal(c.Floats, []float64{1.5, 2.5}) {
			t.Errorf("Floats = %v; want [1.5 2.5]", c.Floats)
		}
	})

	t.Run("FloatLenient", func(t *testing.T) {
		t.Parallel()
		c := typedCore{Type: TypeFloat, Raw: []string{"3.14", "xyz", "NaN"}}
		decodeTyped(&c, []string{"3.14", "xyz", "NaN"})
		if len(c.Floats) != 3 {
			t.Fatalf("len(Floats) = %d; want 3", len(c.Floats))
		}
		if math.Abs(c.Floats[0]-3.14) > 1e-6 {
			t.Errorf("Floats[0] = %v; want ~3.14", c.Floats[0])
		}
		if c.Floats[1] != 0 {
			t.Errorf("Floats[1] = %v; want 0", c.Floats[1])
		}
		if !math.IsNaN(c.Floats[2]) {
			t.Errorf("Floats[2] = %v; want NaN", c.Floats[2])
		}
	})

	t.Run("StringPercentDecoded", func(t *testing.T) {
		t.Parallel()
		c := typedCore{Type: TypeString, Raw: []string{"a%2Cb", "x"}}
		decodeTyped(&c, []string{"a%2Cb", "x"})
		if !slices.Equal(c.Strs, []string{"a,b", "x"}) {
			t.Errorf("Strs = %v; want [a,b x]", c.Strs)
		}
	})

	t.Run("FlagAndCharacterAndUnknownAreNoOps", func(t *testing.T) {
		t.Parallel()
		// Flag, Character, and Unknown are handled by the entry points
		// (decodeInfoValue / decodeFieldValue); decodeTyped must not
		// touch their typed slices.
		for _, ft := range []FieldType{TypeFlag, TypeCharacter, TypeUnknown} {
			c := typedCore{Type: ft, Raw: []string{"x", "y"}}
			decodeTyped(&c, []string{"x", "y"})
			if len(c.Ints) != 0 || len(c.Floats) != 0 || len(c.Strs) != 0 || len(c.Chars) != 0 {
				t.Errorf("%v: typed slices populated: %+v", ft, c)
			}
			if c.Missing {
				t.Errorf("%v: Missing set unexpectedly", ft)
			}
		}
	})
}

// TestWriteTyped pins the per-type emission rules that the encoder
// (writeInfoColumn, encodeSample) both rely on. Missing is rendered as
// "."; Integer/Float/String/Character follow the spec formatting; an
// Unknown type falls back to the Raw tokens.
func TestWriteTyped(t *testing.T) {
	t.Parallel()

	cases := []struct {
		name string
		c    typedCore
		want string
	}{
		{
			name: "Missing",
			c:    typedCore{Type: TypeInteger, Missing: true},
			want: ".",
		},
		{
			name: "Integer",
			c:    typedCore{Type: TypeInteger, Ints: []int64{1, 2, 3}},
			want: "1,2,3",
		},
		{
			name: "Float",
			c:    typedCore{Type: TypeFloat, Floats: []float64{1.5, 2.5}},
			want: "1.5,2.5",
		},
		{
			name: "FloatSpecial",
			c:    typedCore{Type: TypeFloat, Floats: []float64{math.NaN()}},
			want: "NaN",
		},
		{
			name: "StringPercentEncoded",
			c:    typedCore{Type: TypeString, Strs: []string{"a,b", "x"}},
			want: "a%2Cb,x",
		},
		{
			name: "Character",
			c:    typedCore{Type: TypeCharacter, Chars: []byte{'Q'}},
			want: "Q",
		},
		{
			name: "CharacterEmpty",
			c:    typedCore{Type: TypeCharacter},
			want: "",
		},
		{
			name: "UnknownFallsBackToRaw",
			c:    typedCore{Type: TypeUnknown, Raw: []string{"hello", "world"}},
			want: "hello,world",
		},
		{
			name: "UnknownEmptyRawEmitsNothing",
			c:    typedCore{Type: TypeUnknown},
			want: "",
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			t.Parallel()
			var buf bytes.Buffer
			writeTyped(&buf, &tc.c)
			if got := buf.String(); got != tc.want {
				t.Errorf("writeTyped = %q; want %q", got, tc.want)
			}
		})
	}
}
