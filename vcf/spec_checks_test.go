package vcf

import (
	"strings"
	"testing"
)

func TestParser_RejectsDisallowedControlByte(t *testing.T) {
	// Build a VCF where the second line contains a NUL byte
	// (U+0000), which is in the forbidden range.
	src := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\x00\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err == nil {
		t.Fatal("expected error for NUL byte in header, got nil")
	} else if !strings.Contains(err.Error(), "non-printable") {
		t.Errorf("error should mention non-printable, got: %v", err)
	}
}

func TestParser_RejectsAllDisallowedRanges(t *testing.T) {
	// Exercise each forbidden range on a header line so the
	//   Header() call surfaces the error:
	//   U+0000..U+0008
	//   U+000B..U+000C
	//   U+000E..U+001F
	for _, b := range []byte{0x00, 0x05, 0x08, 0x0B, 0x0C, 0x0E, 0x1F} {
		src := "##fileformat=VCFv4.2\n" +
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=" + string([]byte{b}) + "d\">\n" +
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
		dec := NewDecoder(strings.NewReader(src))
		_, err := dec.Header()
		if err == nil {
			t.Errorf("byte 0x%02X: expected error, got nil", b)
			continue
		}
		if !strings.Contains(err.Error(), "non-printable") {
			t.Errorf("byte 0x%02X: error %q should mention non-printable", b, err)
		}
	}
}

func TestParser_AcceptsTabAndOtherAllowedBytes(t *testing.T) {
	// Bytes NOT in the forbidden range must pass. The
	// common case is the tab character (0x09) which is the
	// VCF field separator.
	src := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t1\t.\tA\t.\t.\t.\tDP=1\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	if _, err := dec.Record(); err != nil {
		t.Errorf("tab-containing record: unexpected err: %v", err)
	}
}

func TestValidator_IntegerBCFUnsafeRange(t *testing.T) {
	// The range [-2^31, -2^31+7] is forbidden per spec
	// because BCF cannot store those 8 values.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdrText))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	for _, n := range []int64{-2147483648, -2147483645, -2147483641} {
		rec := &Record{
			Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"."},
			Info: map[string]InfoValue{
				"DP": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{n}}},
			},
		}
		v := NewValidator(ModeStrict)
		issues := v.Validate(hdr, []*Record{rec})
		found := false
		for _, iss := range issues {
			if iss.Field == "INFO/DP" && strings.Contains(iss.Msg, "BCF-unsafe") {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("DP=%d: expected BCF-unsafe issue, got: %+v", n, issues)
		}
	}
	// In-range value: no issue.
	rec := &Record{
		Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"."},
		Info: map[string]InfoValue{
			"DP": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{-2147483640}}},
		},
	}
	issues := NewValidator(ModeStrict).Validate(hdr, []*Record{rec})
	for _, iss := range issues {
		if iss.Field == "INFO/DP" && strings.Contains(iss.Msg, "BCF-unsafe") {
			t.Errorf("DP=-2147483640: unexpected BCF-unsafe issue: %+v", iss)
		}
	}
}

func TestValidator_SVLENMustBePositive(t *testing.T) {
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"SV length\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

	// Negative SVLEN: spec says SVLEN must be positive even
	// for deletions.
	dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t<DEL>\t.\t.\tSVLEN=-5\n"))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "INFO/SVLEN" && strings.Contains(iss.Msg, "must be positive") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("negative SVLEN: expected positive issue, got: %+v", issues)
	}

	// Positive SVLEN: no SVLEN-positive issue.
	dec2 := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t<DEL>\t.\t.\tSVLEN=5\n"))
	if _, err := dec2.Header(); err != nil {
		t.Fatalf("Header2: %v", err)
	}
	rec2, err := dec2.Record()
	if err != nil {
		t.Fatalf("Record2: %v", err)
	}
	issues = v.Validate(hdr, []*Record{rec2})
	for _, iss := range issues {
		if iss.Field == "INFO/SVLEN" && strings.Contains(iss.Msg, "must be positive") {
			t.Errorf("positive SVLEN: unexpected issue: %+v", iss)
		}
	}
}

func TestValidator_CIPOSMustSpanZero(t *testing.T) {
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=CIPOS,Number=.,Type=Integer,Description=\"CI around POS\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdrText))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	// Lower bound > 0: must fail.
	rec := &Record{
		Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"<DEL>"},
		Info: map[string]InfoValue{
			"CIPOS": {typedCore: typedCore{Type: TypeInteger, Ints: []int64{1, 5}}},
		},
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "INFO/CIPOS" && strings.Contains(iss.Msg, "span 0") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("CIPOS lower>0: expected span-0 issue, got: %+v", issues)
	}
	// Upper bound < 0: must fail.
	rec.Info["CIPOS"] = InfoValue{typedCore: typedCore{Type: TypeInteger, Ints: []int64{-5, -1}}}
	issues = v.Validate(hdr, []*Record{rec})
	found = false
	for _, iss := range issues {
		if iss.Field == "INFO/CIPOS" && strings.Contains(iss.Msg, "span 0") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("CIPOS upper<0: expected span-0 issue, got: %+v", issues)
	}
	// Valid pair [-5, 5]: no issue.
	rec.Info["CIPOS"] = InfoValue{typedCore: typedCore{Type: TypeInteger, Ints: []int64{-5, 5}}}
	issues = v.Validate(hdr, []*Record{rec})
	for _, iss := range issues {
		if iss.Field == "INFO/CIPOS" {
			t.Errorf("valid CIPOS pair: unexpected issue: %+v", iss)
		}
	}
}

func TestParser_RejectsInvalidContigName(t *testing.T) {
	// The SAM-compatible regex forbids leading '*' and '='.
	// (A comma in a contig ID would also be invalid, but the
	// current parser splits on comma before validating the
	// ID, so that case is exercised separately if at all.)
	cases := []struct {
		name string
		id   string
	}{
		{"leading-asterisk", "*chr1"},
		{"leading-equals", "=chr1"},
	}
	for _, c := range cases {
		src := "##fileformat=VCFv4.2\n" +
			"##contig=<ID=" + c.id + ",length=1000>\n" +
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
		dec := NewDecoder(strings.NewReader(src))
		_, err := dec.Header()
		if err == nil {
			t.Errorf("%s (%q): expected error, got nil", c.name, c.id)
			continue
		}
		if !strings.Contains(err.Error(), "contig") {
			t.Errorf("%s (%q): error %q should mention contig", c.name, c.id, err)
		}
	}
}

func TestParser_AcceptsValidContigName(t *testing.T) {
	// Common contig name forms: "1", "chr1", "NC_000001.1",
	// names with dots, dashes, pluses.
	for _, id := range []string{"chr1", "1", "NC_000001.1", "chrX", "scaffold_42", "MT", "HLA-DRB1*01:01"} {
		src := "##fileformat=VCFv4.2\n" +
			"##contig=<ID=" + id + ",length=1000>\n" +
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
		dec := NewDecoder(strings.NewReader(src))
		if _, err := dec.Header(); err != nil {
			t.Errorf("ID %q: unexpected err: %v", id, err)
		}
	}
}

func TestParser_RejectsTrailingTabInHeaderLine(t *testing.T) {
	// The spec says the header line must not end with a tab.
	// We append one here.
	src := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t\n"
	dec := NewDecoder(strings.NewReader(src))
	_, err := dec.Header()
	if err == nil {
		t.Fatal("expected error for trailing tab, got nil")
	}
	if !strings.Contains(err.Error(), "tab") {
		t.Errorf("error should mention tab, got: %v", err)
	}
}

func TestParser_RejectsSampleNameWithWhitespace(t *testing.T) {
	src := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA 00001\n"
	dec := NewDecoder(strings.NewReader(src))
	_, err := dec.Header()
	if err == nil {
		t.Fatal("expected error for whitespace in sample name, got nil")
	}
	if !strings.Contains(err.Error(), "whitespace") {
		t.Errorf("error should mention whitespace, got: %v", err)
	}
}

func TestValidator_FloatFormatStrict(t *testing.T) {
	// The spec regex rejects "1." (trailing dot) and "1e"
	// (trailing exponent with no digits), which
	// strconv.ParseFloat accepts.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"freq\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	for _, bad := range []string{"1.", "1e", "1e+", ".e5"} {
		// The parser itself is lenient (it doesn't reject bad
		// float formats at parse time); the validator catches
		// them.
		rec := &Record{
			Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"C"},
			Info: map[string]InfoValue{
				"AF": {typedCore: typedCore{
					Type: TypeFloat,
					Strs: []string{bad},
					Raw:  []string{bad},
				}},
			},
		}
		v := NewValidator(ModeStrict)
		issues := v.Validate(parseHeaderForTest(t, hdrText), []*Record{rec})
		found := false
		for _, iss := range issues {
			if iss.Field == "INFO/AF" && strings.Contains(iss.Msg, "Type=\"Float\"") {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("AF=%q: expected Float format issue, got: %+v", bad, issues)
		}
	}
	// Valid formats: no issue.
	for _, good := range []string{"0.5", "1", "0.5e-3", "1.5E+10", ".5", "INF", "NAN"} {
		rec := &Record{
			Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"C"},
			Info: map[string]InfoValue{
				"AF": {typedCore: typedCore{
					Type: TypeFloat,
					Strs: []string{good},
					Raw:  []string{good},
				}},
			},
		}
		issues := NewValidator(ModeStrict).Validate(parseHeaderForTest(t, hdrText), []*Record{rec})
		for _, iss := range issues {
			if iss.Field == "INFO/AF" && strings.Contains(iss.Msg, "Type=\"Float\"") {
				t.Errorf("AF=%q: unexpected Float format issue: %+v", good, iss)
			}
		}
	}
}

func TestValidator_SVCLAIMMissingForDEL(t *testing.T) {
	// A DEL record without SVCLAIM is a spec violation.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=SVCLAIM,Number=A,Type=String,Description=\"claim\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t<DEL>\t.\t.\t.\n"))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	issues := NewValidator(ModeStrict).Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "INFO/SVCLAIM" && strings.Contains(iss.Msg, "must be specified") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected SVCLAIM-required issue for DEL, got: %+v", issues)
	}
}

func TestValidator_SVCLAIMValidForDEL(t *testing.T) {
	// A DEL record with SVCLAIM=DJ passes.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=SVCLAIM,Number=A,Type=String,Description=\"claim\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t<DEL>\t.\t.\tSVCLAIM=DJ\n"))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	issues := NewValidator(ModeStrict).Validate(hdr, []*Record{rec})
	for _, iss := range issues {
		if iss.Field == "INFO/SVCLAIM" {
			t.Errorf("unexpected SVCLAIM issue for valid DEL+DJ: %+v", iss)
		}
	}
}

func TestValidator_SVCLAIMInvalidForDEL(t *testing.T) {
	// A DEL with SVCLAIM=X (not in {D,J,DJ}) is invalid.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=SVCLAIM,Number=A,Type=String,Description=\"claim\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t<DEL>\t.\t.\tSVCLAIM=X\n"))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	issues := NewValidator(ModeStrict).Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "INFO/SVCLAIM" && strings.Contains(iss.Msg, "must be one of") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected SVCLAIM invalid-claim issue, got: %+v", issues)
	}
}

func TestValidator_SVCLAIMNonSVAlleleMustBeDot(t *testing.T) {
	// A non-SV (SNV) record with SVCLAIM=A is invalid.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=SVCLAIM,Number=A,Type=String,Description=\"claim\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC\t.\t.\tSVCLAIM=A\n"))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	issues := NewValidator(ModeStrict).Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "INFO/SVCLAIM" && strings.Contains(iss.Msg, "non-SV") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected non-SV SVCLAIM issue, got: %+v", issues)
	}
}