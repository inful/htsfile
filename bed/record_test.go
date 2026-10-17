package bed

import "testing"

func TestRecord_ColumnCount(t *testing.T) {
	t.Parallel()
	cases := []struct {
		name string
		r    Record
		want int
	}{
		{"BED3 (zero value)", Record{}, 3},
		{"BED4 with name", Record{Name: strPtr("peak1")}, 4},
		{"BED5 with score", Record{Score: intPtr(500)}, 5},
		{"BED6 with strand", Record{Strand: strandPtr(StrandPlus)}, 6},
		{"BED7 with thickStart", Record{ThickStart: int64Ptr(100)}, 7},
		{"BED8 with thickEnd", Record{ThickEnd: int64Ptr(200)}, 8},
		{"BED9 with itemRgb", Record{ItemRgb: &RGB{255, 0, 0}}, 9},
		{
			"BED10 with blockCount",
			Record{BlockCount: intPtr(2)},
			10,
		},
		{
			"BED11 with blockSizes",
			Record{BlockSizes: []int64{100, 200}},
			11,
		},
		{
			"BED12 full",
			Record{
				Name:        strPtr("gene"),
				Score:       intPtr(900),
				Strand:      strandPtr(StrandMinus),
				ThickStart:  int64Ptr(50),
				ThickEnd:    int64Ptr(450),
				ItemRgb:     &RGB{0, 255, 0},
				BlockCount:  intPtr(2),
				BlockSizes:  []int64{100, 200},
				BlockStarts: []int64{0, 150},
			},
			12,
		},
		{
			// Non-contiguous prefix (BlockStarts set without BlockCount)
			// — ColumnCount reflects the highest set field (12); the
			// Validator will report this as an error.
			"BED12 with non-contiguous prefix → 12",
			Record{
				Name:        strPtr("g"),
				BlockStarts: []int64{0},
			},
			12,
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			t.Parallel()
			if got := c.r.ColumnCount(); got != c.want {
				t.Errorf("ColumnCount = %d; want %d", got, c.want)
			}
		})
	}
}

func TestRecord_HasBlocks(t *testing.T) {
	t.Parallel()
	if (&Record{BlockCount: intPtr(2), BlockSizes: []int64{1, 2}, BlockStarts: []int64{0, 1}}).HasBlocks() != true {
		t.Error("HasBlocks = false; want true when all three block fields are set")
	}
	if (&Record{BlockCount: intPtr(2)}).HasBlocks() != false {
		t.Error("HasBlocks = true; want false when only BlockCount is set")
	}
	if (&Record{}).HasBlocks() != false {
		t.Error("HasBlocks = true; want false on zero value")
	}
}

func TestStrand_Parse(t *testing.T) {
	t.Parallel()
	cases := []struct {
		in      byte
		want    Strand
		wantErr bool
	}{
		{'+', StrandPlus, false},
		{'-', StrandMinus, false},
		{'.', StrandNone, false},
		{'*', 0, true},
		{'x', 0, true},
		{0, 0, true},
	}
	for _, c := range cases {
		got, err := ParseStrand(c.in)
		if c.wantErr {
			if err == nil {
				t.Errorf("ParseStrand(%q) = nil err; want non-nil", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseStrand(%q) unexpected err: %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseStrand(%q) = %q; want %q", c.in, got, c.want)
		}
		if got.String() != string(c.in) {
			t.Errorf("Strand(%q).String() = %q; want %q", c.in, got.String(), string(c.in))
		}
	}
}

func TestParseError_Message(t *testing.T) {
	t.Parallel()

	t.Run("with field", func(t *testing.T) {
		t.Parallel()
		e := &ParseError{Line: 42, Field: "chromEnd", Msg: "non-numeric: foo"}
		want := "line 42 chromEnd: non-numeric: foo"
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

// strPtr, intPtr, int64Ptr, strandPtr are tiny helpers that return
// the address of a literal so tests can construct optional-field
// records without ceremony.
func strPtr(s string) *string    { return &s }
func intPtr(i int) *int          { return &i }
func int64Ptr(i int64) *int64    { return &i }
func strandPtr(s Strand) *Strand { return &s }
