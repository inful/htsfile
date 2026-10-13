package vcf

import "testing"

func TestInfoValue_Accessors(t *testing.T) {
	t.Parallel()
	t.Run("String", func(t *testing.T) {
		t.Parallel()
		v := InfoValue{Type: TypeString, Strs: []string{"hello", "world"}}
		if got := v.String(); got != "hello" {
			t.Errorf("String() = %q; want %q", got, "hello")
		}
		if (InfoValue{}).String() != "" {
			t.Error("zero-value should yield empty string")
		}
	})
	t.Run("Int", func(t *testing.T) {
		t.Parallel()
		v := InfoValue{Type: TypeInteger, Ints: []int64{42, 7}}
		if got := v.Int(); got != 42 {
			t.Errorf("Int() = %d; want 42", got)
		}
		if (InfoValue{}).Int() != 0 {
			t.Error("zero-value should yield 0")
		}
	})
	t.Run("Float", func(t *testing.T) {
		t.Parallel()
		v := InfoValue{Type: TypeFloat, Floats: []float64{3.14, 2.72}}
		if got := v.Float(); got != 3.14 {
			t.Errorf("Float() = %v; want 3.14", got)
		}
		if (InfoValue{}).Float() != 0 {
			t.Error("zero-value should yield 0")
		}
	})
	t.Run("HasValue", func(t *testing.T) {
		t.Parallel()
		v := InfoValue{Missing: true}
		if v.HasValue() {
			t.Error("missing should not have value")
		}
		v = InfoValue{Type: TypeFlag}
		if v.HasValue() {
			t.Error("flag not set should not have value")
		}
		v = InfoValue{Type: TypeFlag, Flag: true}
		if !v.HasValue() {
			t.Error("flag set should have value")
		}
	})
}

func TestFieldValue_Accessors(t *testing.T) {
	t.Parallel()
	t.Run("String", func(t *testing.T) {
		t.Parallel()
		v := FieldValue{Type: TypeString, Strs: []string{"0/1"}}
		if got := v.String(); got != "0/1" {
			t.Errorf("String() = %q", got)
		}
	})
	t.Run("Int", func(t *testing.T) {
		t.Parallel()
		v := FieldValue{Type: TypeInteger, Ints: []int64{99}}
		if got := v.Int(); got != 99 {
			t.Errorf("Int() = %d", got)
		}
	})
	t.Run("Float", func(t *testing.T) {
		t.Parallel()
		v := FieldValue{Type: TypeFloat, Floats: []float64{0.95}}
		if got := v.Float(); got != 0.95 {
			t.Errorf("Float() = %v", got)
		}
	})
	t.Run("HasValue", func(t *testing.T) {
		t.Parallel()
		v := FieldValue{Missing: true}
		if v.HasValue() {
			t.Error("missing should not have value")
		}
		v = FieldValue{Type: TypeString, Strs: []string{"x"}}
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