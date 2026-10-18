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

func TestValidator_CIPOSCIENDCILENCardinality(t *testing.T) {
	// VCFv4.5 §1.6.1.7: "If present, the number of
	// entries must be twice the number of ALT alleles."
	hdrText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=CIPOS,Number=.,Type=Integer,Description=\"x\">\n" +
		"##INFO=<ID=CIEND,Number=.,Type=Integer,Description=\"x\">\n" +
		"##INFO=<ID=CILEN,Number=.,Type=Integer,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		{"biallelic 4 entries ok", "CIPOS=-5,5,0,0", false},
		{"biallelic 2 entries wrong", "CIPOS=-5,5", true},
		{"biallelic CIEND ok", "CIEND=-5,5,0,0", false},
		{"biallelic CILEN ok", "CILEN=-5,5,0,0", false},
		{"biallelic CILEN wrong", "CILEN=-5,5,3", true},
		{"single alt 2 entries ok", "CIPOS=-5,5", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			alt := "C,G"
			if c.name == "single alt 2 entries ok" {
				alt = "C"
			}
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + alt + "\t.\t.\t" + c.info + "\n"))
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
				if strings.Contains(iss.Msg, "must contain exactly") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("info=%q alt=%s: found=%v, want %v; issues=%+v", c.info, alt, found, c.want, issues)
			}
		})
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

func TestIsReservedEventType(t *testing.T) {
	// All 15 reserved values must return true.
	for _, v := range []string{
		"DEL", "DEL:ME", "INS", "INS:ME",
		"DUP", "DUP:TANDEM", "DUP:DISPERSED",
		"INV", "TRA", "TRA:BALANCED", "TRA:UNBALANCED",
		"CHROMOTHRIPSIS", "CHROMOPLEXY", "BFB", "DOUBLEMINUTE",
	} {
		if !IsReservedEventType(v) {
			t.Errorf("IsReservedEventType(%q): should be true", v)
		}
	}
	// Implementation-defined (non-reserved) values return false.
	for _, v := range []string{"CUSTOM", "MY_SV_TYPE", "DEL_CIRCULAR", ""} {
		if IsReservedEventType(v) {
			t.Errorf("IsReservedEventType(%q): should be false", v)
		}
	}
}

func TestValidator_EventTypeReserved(t *testing.T) {
	// A reserved EVENTTYPE value must not produce an issue.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=EVENTTYPE,Number=A,Type=String,Description=\"event type\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	for _, et := range []string{"DEL", "INS", "DUP", "TRA", "DEL:ME", "DUP:TANDEM", "CHROMOTHRIPSIS"} {
		dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC\t.\t.\tEVENTTYPE=" + et + "\n"))
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
			if iss.Field == "INFO/EVENTTYPE" {
				t.Errorf("EVENTTYPE=%q: unexpected issue: %+v", et, iss)
			}
		}
	}
}

func TestValidator_EventTypeUnreserved(t *testing.T) {
	// An unreserved EVENTTYPE value: warning in default
	// mode, error in ModeStrict.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=EVENTTYPE,Number=A,Type=String,Description=\"event type\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

	dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC\t.\t.\tEVENTTYPE=MY_SV\n"))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	// Default mode: warning.
	issues := NewValidator(ModeDefault).Validate(hdr, []*Record{rec})
	found := false
	for _, iss := range issues {
		if iss.Field == "INFO/EVENTTYPE" {
			if iss.Severity != SeverityWarning {
				t.Errorf("default mode: severity %v, want Warning", iss.Severity)
			}
			found = true
		}
	}
	if !found {
		t.Errorf("default mode: expected EVENTTYPE warning, got: %+v", issues)
	}
	// Strict mode: error.
	issues = NewValidator(ModeStrict).Validate(hdr, []*Record{rec})
	found = false
	for _, iss := range issues {
		if iss.Field == "INFO/EVENTTYPE" {
			if iss.Severity != SeverityError {
				t.Errorf("strict mode: severity %v, want Error", iss.Severity)
			}
			found = true
		}
	}
	if !found {
		t.Errorf("strict mode: expected EVENTTYPE error, got: %+v", issues)
	}
}

func TestValidator_LAAValueRange(t *testing.T) {
	// LAA is a 1-based index into ALT (per VCFv4.5 §1.6.1.4).
	// Out-of-range values (0, or > len(ALT)) are an error.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"local alleles\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	for _, bad := range []string{"0", "4", "-1"} {
		// 2 ALTs, so valid LAA is {1, 2}. 0, 4, -1 all fail.
		dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC,G\t.\t.\t.\tGT:LAA\t0/1:" + bad + "\n"))
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
			if iss.Field == "FORMAT/LAA/sample[0]" && strings.Contains(iss.Msg, "out of range") {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("LAA=%q: expected out-of-range issue, got: %+v", bad, issues)
		}
	}
	// Valid LAA values: no issue.
	for _, good := range []string{"1", "2", "1,2", "2,1"} {
		dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC,G\t.\t.\t.\tGT:LAA\t0/1:" + good + "\n"))
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
			if iss.Field == "FORMAT/LAA/sample[0]" {
				t.Errorf("LAA=%q: unexpected issue: %+v", good, iss)
			}
		}
	}
}

func TestValidator_LAADistinctness(t *testing.T) {
	// LAA values must be distinct. A duplicate fires a
	// SeverityError.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"local alleles\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC,G\t.\t.\t.\tGT:LAA\t0/1:1,1\n"))
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
		if iss.Field == "FORMAT/LAA/sample[0]" && strings.Contains(iss.Msg, "duplicates") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected LAA duplicate issue, got: %+v", issues)
	}
}

func TestValidator_POSWithinContigLength(t *testing.T) {
	// VCFv4.5 §1.6.1.2: a record's POS must not exceed the
	// declared ##contig length, and POS+len(REF)-1 must
	// not extend past it either.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name  string
		line  string
		field string // expected field on the issue ("POS" or "REF")
		want  bool
	}{
		{"pos within", "chr1\t100\t.\tA\tC\t.\t.\t.", "POS", false},
		{"ref end within", "chr1\t999\t.\tA\tC\t.\t.\t.", "POS", false},
		{"pos exceeds", "chr1\t2000\t.\tA\tC\t.\t.\t.", "POS", true},
		{"ref extends past", "chr1\t1000\t.\tAA\tC\t.\t.\t.", "REF", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + c.line + "\n"))
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
				if iss.Field == c.field && (strings.Contains(iss.Msg, "contig length") || strings.Contains(iss.Msg, "extends past")) {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("issue found=%v, want %v; issues=%+v", found, c.want, issues)
			}
		})
	}
}

func TestValidator_POSWithinContigLength_Undeclared(t *testing.T) {
	// When the CHROM is not declared, the position
	// constraint must not fire (no length to compare
	// against).
	hdrText := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC\t.\t.\t.\n"))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	for _, iss := range NewValidator(ModeStrict).Validate(hdr, []*Record{rec}) {
		if strings.Contains(iss.Msg, "contig length") {
			t.Errorf("unexpected contig-length issue for undeclared contig: %+v", iss)
		}
	}
}

func TestValidator_SymbolicALTDeclaredInHeader(t *testing.T) {
	// A symbolic ALT allele like <DEL> must be declared in
	// ##ALT before it is used in a record. <*>, breakend
	// forms, and a declared ##ALT pass without issue.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=DEL,Description=\"deletion\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, alt string
		want      bool
	}{
		{"declared", "<DEL>", false},
		{"undeclared", "<UNKNOWN>", true},
		{"star allele is allowed", "<*>", false},
		{"breakend is allowed", "A[chr1:200[", false},
		{"non-symbolic", "C", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + c.alt + "\t.\t.\t.\n"))
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
				if strings.HasPrefix(iss.Field, "ALT[") && strings.Contains(iss.Msg, "not declared in any ##ALT") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("ALT=%q: found=%v, want %v; issues=%+v", c.alt, found, c.want, issues)
			}
		})
	}
}

func TestValidator_METAValuesEnforced(t *testing.T) {
	// Per VCFv4.5 §1.4.7 "Sample field format": when a
	// ##META entry declares a Values list, every
	// ##SAMPLE field with that key must take a value
	// from the list. Fires as a SeverityWarning in
	// ModeStrict; skipped in ModeDefault.
	cases := []struct {
		name  string
		samp  string
		want  bool
	}{
		{"allowed", "Assay=WholeGenome", false},
		{"not in list", "Assay=Bogus", true},
		{"no values declared", "OtherField=anything", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			hdr := "##fileformat=VCFv4.2\n" +
				"##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>\n" +
				"##META=<ID=OtherField,Type=String,Number=.>\n" +
				"##SAMPLE=<ID=Sample1," + c.samp + ">\n" +
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
			dec := NewDecoder(strings.NewReader(hdr))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			issues := NewValidator(ModeStrict).HeaderChecks(h)
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "##SAMPLE") && strings.Contains(iss.Msg, "not in the ##META") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("%s: found=%v, want %v; issues=%+v", c.name, found, c.want, issues)
			}
		})
	}
}

func TestValidator_METAValuesDefaultModeSkipsCheck(t *testing.T) {
	// The check fires only in ModeStrict. In ModeDefault
	// the same input passes without issues.
	hdr := "##fileformat=VCFv4.2\n" +
		"##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>\n" +
		"##SAMPLE=<ID=Sample1,Assay=Bogus>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdr))
	h, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if issues := NewValidator(ModeDefault).HeaderChecks(h); len(issues) != 0 {
		t.Errorf("ModeDefault should skip META Values check; got: %+v", issues)
	}
}

func TestValidator_LocalAlleleRequiresLAA(t *testing.T) {
	// Per VCFv4.5 §1.6.2: "If any local-allele field is
	// present, LAA must also be present". A sample that
	// declares LAD / LGL / LPL / etc. without LAA fires
	// a SeverityError in ModeStrict.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"l\">\n" +
		"##FORMAT=<ID=LAD,Number=.,Type=Integer,Description=\"l\">\n" +
		"##FORMAT=<ID=LGL,Number=LG,Type=Float,Description=\"l\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, sample string
		want         bool
	}{
		{"LAA and LAD ok", "0/1:2:5,10", false},
		{"LAD without LAA", "0/1::5,10", true},
		{"LGL without LAA", "0/1::0.0,1.0,2.0", true},
		{"LAA alone ok", "0/1:2", false},
		{"no local allele ok", "0/1", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC,G\t.\t.\t.\tGT:LAA:LAD:LGL\t" + c.sample + "\n"))
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
				if strings.Contains(iss.Msg, "requires LAA") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("sample=%q: found=%v, want %v; issues=%+v", c.sample, found, c.want, issues)
			}
		})
	}
}

func TestParser_NonRefAlias(t *testing.T) {
	// VCFv4.5 §1.7: <NON_REF> is an alias for <*> for
	// gVCF backwards compatibility. The parser must
	// normalise the alias to <*> so downstream code
	// only sees the canonical form.
	src := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tG\t<NON_REF>\t.\t.\t.\n"
	dec := NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if len(rec.Alt) != 1 || rec.Alt[0] != "<*>" {
		t.Errorf("rec.Alt = %q, want [\"<*>\"]", rec.Alt)
	}
}

func TestValidator_PSAndPSLMutuallyExclusive(t *testing.T) {
	// Per VCFv4.5 §1.6.2: "A given sample-genotype must
	// not have values for both PS and PSL."
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"p\">\n" +
		"##FORMAT=<ID=PSL,Number=P,Type=String,Description=\"l\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, sample string
		want         bool
	}{
		{"only PS", "0|1:100", false},
		{"only PSL", "0|1:chr1*100*1,chr1*200*1", false},
		{"both PS and PSL", "0|1:100:chr1*100*1,chr1*200*1", true},
		{"neither", "0|1", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC,G\t.\t.\t.\tGT:PS:PSL\t" + c.sample + "\n"))
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
				if strings.Contains(iss.Msg, "PS and PSL cannot both be present") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("sample=%q: found=%v, want %v; issues=%+v", c.sample, found, c.want, issues)
			}
		})
	}
}

func TestValidator_ENDMatchesComputed(t *testing.T) {
	// END is a computed field per VCFv4.5 §1.6.1.6. When
	// declared, the value must match the computed END
	// (POS + len(REF) - 1 for a non-symbolic ALT).
	hdrText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"end position\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		{"matches computed", "END=103", false},        // POS=100, REF=ACGT => END=103
		{"does not match", "END=200", true},
		{"missing END", ".", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tACGT\tC\t.\t.\t" + c.info + "\n"))
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
				if iss.Field == "INFO/END" {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("info=%q: found=%v, want %v; issues=%+v", c.info, found, c.want, issues)
			}
		})
	}
}

func TestValidator_CNVAllelesShareSVLENWhenFormatCN(t *testing.T) {
	// Per VCFv4.5 §"FORMAT keys used for structural
	// variants": "When FORMAT CN is present, all <CNV>,
	// <DEL> and <DUP> alleles must have the same SVLEN."
	hdrText := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=CNV,Description=\"x\">\n" +
		"##ALT=<ID=DEL,Description=\"x\">\n" +
		"##ALT=<ID=DUP,Description=\"x\">\n" +
		"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"x\">\n" +
		"##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"x\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, alt, svlen string
		want             bool
	}{
		{"CNV all share", "<CNV>,<CNV>", "100,100", false},
		{"CNV different", "<CNV>,<CNV>", "100,200", true},
		{"DEL DUP same", "<DEL>,<DUP>", "50,50", false},
		{"DEL DUP different", "<DEL>,<DUP>", "50,100", true},
		{"mixed INV same ok", "<DEL>,<INV>,<DUP>", "30,30,30", false},
		{"mixed INV different", "<DEL>,<INV>,<DUP>", "30,30,40", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + c.alt + "\t.\t.\tSVLEN=" + c.svlen + "\tGT:CN\t0/0:1.0\n"))
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
				if strings.Contains(iss.Msg, "must share the same SVLEN") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("alt=%s svlen=%s: found=%v, want %v; issues=%+v", c.alt, c.svlen, found, c.want, issues)
			}
		})
	}
}

func TestValidator_FTNoWhitespaceOrSemicolons(t *testing.T) {
	// Per VCFv4.5 §1.6.2 (FT field): "No whitespace or
	// semicolons permitted." A per-sample check fires
	// when an FT token contains a space, tab, newline,
	// or semicolon.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=FT,Number=1,Type=String,Description=\"f\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, sample string
		want         bool
	}{
		{"clean", "0/1:PASS", false},
		{"space inside token", "0/1:q 10", true},
		{"semicolon inside token", "0/1:q10;s50", true},
		{"multiple tokens one bad", "0/1:PASS,q 10", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC\t.\t.\t.\tGT:FT\t" + c.sample + "\n"))
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
				if iss.Field == "FORMAT/FT/sample[0]" {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("sample=%q: found=%v, want %v; issues=%+v", c.sample, found, c.want, issues)
			}
		})
	}
}

func TestValidator_RNDeterminedFieldsCardinality(t *testing.T) {
	// Per VCFv4.5 §1.6.1.7: "The length of these fields
	// [RUS, RUL, RB, CIRUC, RUC, CIRB] must match the
	// sum of all RN values for the record."
	hdrText := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=CNV:TR,Description=\"x\">\n" +
		"##INFO=<ID=RN,Number=A,Type=Integer,Description=\"x\">\n" +
		"##INFO=<ID=RUS,Number=.,Type=String,Description=\"x\">\n" +
		"##INFO=<ID=RUL,Number=.,Type=Integer,Description=\"x\">\n" +
		"##INFO=<ID=RB,Number=.,Type=Integer,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		// Two CNV:TR alleles, RN per allele = 2,2 -> sum 4
		{"RN=2,2 RUS=4 ok", "RN=2,2;RUS=CAG,TG,CAGG,CAA", false},
		{"RN=2,2 RUS=3 wrong", "RN=2,2;RUS=CAG,TG,CAGG", true},
		{"RN=2,2 RUL=4 ok", "RN=2,2;RUL=3,2,4,1", false},
		{"RN=2,2 RUL=3 wrong", "RN=2,2;RUL=3,2,4", true},
		{"RN=2,2 RB=4 ok", "RN=2,2;RB=30,15,12,20", false},
		{"RN=2,2 RB=3 wrong", "RN=2,2;RB=30,15,12", true},
		// RUS length matches RUL per element.
		{"RUS RUL aligned ok", "RN=2,2;RUS=CAG,TG,CAA,TT;RUL=3,2,3,2", false},
		{"RUS RUL mismatch", "RN=2,2;RUS=CAG,TG,CAA,TT;RUL=3,2,3,5", true},
		// No <CNV:TR> → skip the RN checks.
		{"no CNV:TR RUS mismatch", "RUS=CAG", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			alt := "<CNV:TR>,<CNV:TR>"
			if c.name == "no CNV:TR RUS mismatch" {
				alt = "C"
			}
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + alt + "\t.\t.\t" + c.info + "\n"))
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
				if strings.Contains(iss.Msg, "does not match") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("info=%q alt=%s: found=%v, want %v; issues=%+v", c.info, alt, found, c.want, issues)
			}
		})
	}
}

func TestValidator_CNVTRRequiresRUSorRUL(t *testing.T) {
	// Per VCFv4.5 §"Representing tandem repeats":
	// "RUS or RUL must be specified for each <CNV:TR>."
	hdrText := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=CNV:TR,Description=\"x\">\n" +
		"##INFO=<ID=RUS,Number=.,Type=String,Description=\"x\">\n" +
		"##INFO=<ID=RUL,Number=.,Type=Integer,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		{"RUS declared", "RUS=CAG", false},
		{"RUL declared", "RUL=3", false},
		{"neither", ".", true},
		{"empty RUS only", "RUS=.", true},
		{"no CNV:TR in record", "RUS=CAG", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			alt := "<CNV:TR>"
			if c.name == "no CNV:TR in record" {
				alt = "C"
			}
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + alt + "\t.\t.\t" + c.info + "\n"))
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
				if strings.Contains(iss.Msg, "RUS nor RUL") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("info=%q alt=%s: found=%v, want %v; issues=%+v", c.info, alt, found, c.want, issues)
			}
		})
	}
}

func TestValidator_CICNCardinality(t *testing.T) {
	// VCFv4.5 §1.6.1.7: "CICN ... If present, the
	// number of entries must be twice the number of
	// ALT alleles." The field is a pair of
	// (lower, upper) copy-number bounds per allele.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=CICN,Number=.,Type=Float,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		{"biallelic 4 entries ok", "CICN=1.5,2.0,1.0,0.5", false},
		{"biallelic 2 entries wrong", "CICN=1.5,2.0", true},
		{"single alt 2 entries ok", "CICN=1.5,2.0", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			alt := "C,G"
			if c.name == "single alt 2 entries ok" {
				alt = "C"
			}
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + alt + "\t.\t.\t" + c.info + "\n"))
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
				if strings.Contains(iss.Msg, "CICN must contain exactly") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("info=%q alt=%s: found=%v, want %v; issues=%+v", c.info, alt, found, c.want, issues)
			}
		})
	}
}

func TestValidator_MEINFOMETRANSCardinality(t *testing.T) {
	// VCFv4.5 §1.6.1.7: MEINFO and METRANS each contain
	// 4 entries per ALT allele (NAME, START, END,
	// POLARITY).
	hdrText := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=MEINFO,Number=.,Type=String,Description=\"x\">\n" +
		"##INFO=<ID=METRANS,Number=.,Type=String,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		{"biallelic MEINFO 8 ok", "MEINFO=AluY,100,200,+,SVA,300,400,-", false},
		{"biallelic MEINFO wrong", "MEINFO=AluY,100,200,+", true},
		{"biallelic METRANS 8 ok", "METRANS=chr1,100,200,+,chr2,300,400,-", false},
		{"single alt MEINFO 4 ok", "MEINFO=AluY,100,200,+", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			alt := "C,G"
			if c.name == "single alt MEINFO 4 ok" {
				alt = "C"
			}
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + alt + "\t.\t.\t" + c.info + "\n"))
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
				if strings.Contains(iss.Msg, "must contain exactly") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("info=%q alt=%s: found=%v, want %v; issues=%+v", c.info, alt, found, c.want, issues)
			}
		})
	}
}

func TestValidator_NumberGKeysMustHaveFloatOrInteger(t *testing.T) {
	// PL, GL, GP all use Number=G. The Type must be
	// Integer for PL and Float for GL/GP. Declaring
	// them with the wrong Type fires a header-reserved-
	// key warning in ModeDefault, an error in
	// ModeStrict.
	cases := []struct {
		name, hdrKey, hdrType string
		wantIssue             bool
	}{
		{"PL integer", "PL", "Integer", false},
		{"PL float wrong", "PL", "Float", true},
		{"GL float", "GL", "Float", false},
		{"GL integer wrong", "GL", "Integer", true},
		{"GP float", "GP", "Float", false},
		{"GP integer wrong", "GP", "Integer", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			hdrText := "##fileformat=VCFv4.2\n" +
				"##FORMAT=<ID=" + c.hdrKey + ",Number=G,Type=" + c.hdrType + ",Description=\"x\">\n" +
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
				"chr1\t100\t.\tA\tC\t.\t.\t.\n"
			dec := NewDecoder(strings.NewReader(hdrText))
			hdr, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			issues := NewValidator(ModeStrict).HeaderChecks(hdr)
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "should be Number=G,Type=") {
					found = true
					break
				}
			}
			if found != c.wantIssue {
				t.Errorf("hdrKey=%s,Type=%s: issue found=%v, want %v; issues=%+v",
					c.hdrKey, c.hdrType, found, c.wantIssue, issues)
			}
		})
	}
}

func TestValidator_PLCardinalityMatchesPloidy(t *testing.T) {
	// PL is Number=G; for a diploid with biallelic ALT
	// (REF + 1 ALT = 2 alleles) the expected count is
	// (2+2-1 choose 2) = 3. For triploid with biallelic
	// it is (2+3-1 choose 3) = 4.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"pl\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, gt, alt, pl string
		want              bool
	}{
		{"diploid biallelic ok", "0/1", "C", "0,10,100", false},
		{"diploid biallelic wrong", "0/1", "C", "0,10", true},
		{"haploid biallelic ok", "1", "C", "0,100", false},
		{"haploid biallelic wrong", "1", "C", "0,10,100", true},
		{"diploid with 2 alts ok", "1/2", "C,G", "0,10,100,1000,10000,100000", false},
		{"diploid with 2 alts wrong", "1/2", "C,G", "0,10,100", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + c.alt + "\t.\t.\t.\tGT:PL\t" + c.gt + ":" + c.pl + "\n"))
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
				if iss.Field == "FORMAT/PL/sample[0]" && strings.Contains(iss.Msg, "Number=G") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("gt=%s alt=%s pl=%s: found=%v, want %v; issues=%+v", c.gt, c.alt, c.pl, found, c.want, issues)
			}
		})
	}
}

func TestValidator_ADCardinalityMatchesAlleleCount(t *testing.T) {
	// AD is Number=R; the expected count per sample is
	// len(REF) + len(ALT). Diploid biallelic gives 2,
	// diploid with 2 ALTs gives 3.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"allele depths\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, alt, ad string
		want          bool
	}{
		{"biallelic ok", "C", "10,20", false},
		{"biallelic wrong", "C", "10", true},
		{"two alts ok", "C,G", "10,20,30", false},
		{"two alts wrong", "C,G", "10,20", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\t" + c.alt + "\t.\t.\t.\tGT:AD\t0/1:" + c.ad + "\n"))
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
				if iss.Field == "FORMAT/AD/sample[0]" && strings.Contains(iss.Msg, "Number=R") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("alt=%s ad=%s: found=%v, want %v; issues=%+v", c.alt, c.ad, found, c.want, issues)
			}
		})
	}
}

func TestValidator_PloidyConsistencyAcrossRecord(t *testing.T) {
	// Per VCFv4.5 §1.6.2 all samples in a record must
	// share the same ploidy. Mixed ploidies in a single
	// record fire a SeverityError in ModeStrict.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
	cases := []struct {
		name, samples string
		want          bool
	}{
		{"all diploid", "0/1\t0/0\t1/1", false},
		{"all haploid", "1\t0\t1", false},
		{"haploid + diploid mixed", "1\t0/1\t0", true},
		{"diploid + triploid", "0/1\t0/0/1\t0", true},
	}
	// Single-sample case: skip the header check (the
	// shared header declares 3 sample columns).
	hdrOne := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
	cases = append(cases, struct {
		name, samples string
		want          bool
	}{"single sample", "0/1", false})
	one := cases[len(cases)-1]
	cases = cases[:len(cases)-1]
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t" + c.samples + "\n"))
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
				if iss.Field == "FORMAT/GT" && strings.Contains(iss.Msg, "different ploidies") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("samples=%q: found=%v, want %v; issues=%+v", c.samples, found, c.want, issues)
			}
		})
	}

	t.Run(one.name, func(t *testing.T) {
		dec := NewDecoder(strings.NewReader(hdrOne + "chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t" + one.samples + "\n"))
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
			if iss.Field == "FORMAT/GT" && strings.Contains(iss.Msg, "different ploidies") {
				found = true
				break
			}
		}
		if found != one.want {
			t.Errorf("samples=%q: found=%v, want %v; issues=%+v", one.samples, found, one.want, issues)
		}
	})
}

func TestParser_IUPACALTReduction(t *testing.T) {
	// Per VCFv4.5 §1.4.4: symbolic alleles can be used
	// to represent genuinely ambiguous IUPAC codes
	// (e.g. ##ALT=<ID=R,Description="IUPAC code R =
	// A/G">). When such an ALT allele appears in a record
	// and the corresponding ##ALT is declared, the parser
	// reduces the allele to one of its constituent bases.
	in := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=R,Description=\"IUPAC code R = A/G\">\n" +
		"##ALT=<ID=M,Description=\"IUPAC code M = A/C\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tR,M\t.\t.\t.\n"
	dec := NewDecoder(strings.NewReader(in))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if len(rec.Alt) != 2 || rec.Alt[0] != "A" || rec.Alt[1] != "A" {
		t.Errorf("rec.Alt = %q, want [A A]", rec.Alt)
	}
}

func TestParser_IUPACALTUndeclaredRejected(t *testing.T) {
	// An IUPAC code used as an ALT without a matching
	// ##ALT declaration is rejected by the
	// symbolic-ALT-declaration check in the validator
	// (ModeStrict).
	in := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tR\t.\t.\t.\n"
	dec := NewDecoder(strings.NewReader(in))
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
		if iss.Field == "ALT[0]" && strings.Contains(iss.Msg, "not declared in any ##ALT") {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("expected validation issue for undeclared IUPAC ALT, got: %+v", issues)
	}
}

func TestValidator_MitochondrialChromExemptions(t *testing.T) {
	// VCFv4.5 §1.6.1.6: the mitochondrial reference is
	// circular. POS may exceed the contig length and END
	// may not match the simple POS+len(REF)-1 formula
	// because of origin-spanning variants. Both checks
	// must skip mitochondrial chromosomes.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chrM,length=16569>\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"end position\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	// POS > length: a normal chromosome would be
	// rejected; chrM must be accepted.
	t.Run("pos exceeds length", func(t *testing.T) {
		dec := NewDecoder(strings.NewReader(hdrText + "chrM\t20000\t.\tA\tC\t.\t.\t.\tGT\t0/1\n"))
		hdr, err := dec.Header()
		if err != nil {
			t.Fatalf("Header: %v", err)
		}
		rec, err := dec.Record()
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		for _, iss := range NewValidator(ModeStrict).Validate(hdr, []*Record{rec}) {
			if strings.Contains(iss.Msg, "contig length") {
				t.Errorf("mtDNA: unexpected length issue: %+v", iss)
			}
		}
	})
	// END that doesn't match the linear formula: a
	// normal chromosome would be rejected; chrM must
	// be accepted.
	t.Run("end wraps past origin", func(t *testing.T) {
		dec := NewDecoder(strings.NewReader(hdrText + "chrM\t16000\t.\tACGT\tA\t.\t.\tEND=50\tGT\t0/1\n"))
		hdr, err := dec.Header()
		if err != nil {
			t.Fatalf("Header: %v", err)
		}
		rec, err := dec.Record()
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		for _, iss := range NewValidator(ModeStrict).Validate(hdr, []*Record{rec}) {
			if iss.Field == "INFO/END" {
				t.Errorf("mtDNA: unexpected END issue: %+v", iss)
			}
		}
	})
	// Other mitochondrial name conventions must also
	// be exempt.
	for _, name := range []string{"MT", "M", "chrMT"} {
		t.Run("alt name "+name, func(t *testing.T) {
			hdr := "##fileformat=VCFv4.2\n" +
				"##contig=<ID=" + name + ",length=100>\n" +
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
			dec := NewDecoder(strings.NewReader(hdr + name + "\t200\t.\tA\tC\t.\t.\t.\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			for _, iss := range NewValidator(ModeStrict).Validate(h, []*Record{rec}) {
				if strings.Contains(iss.Msg, "contig length") {
					t.Errorf("mtDNA %q: unexpected length issue: %+v", name, iss)
				}
			}
		})
	}
	// Non-mitochondrial chromosomes are still checked.
	t.Run("non-mt still checked", func(t *testing.T) {
		hdr := "##fileformat=VCFv4.2\n" +
			"##contig=<ID=chr1,length=1000>\n" +
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
		dec := NewDecoder(strings.NewReader(hdr + "chr1\t2000\t.\tA\tC\t.\t.\t.\n"))
		h, err := dec.Header()
		if err != nil {
			t.Fatalf("Header: %v", err)
		}
		rec, err := dec.Record()
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		found := false
		for _, iss := range NewValidator(ModeStrict).Validate(h, []*Record{rec}) {
			if strings.Contains(iss.Msg, "contig length") {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("chr1 POS=2000 should still fire the length check")
		}
	})
}

func TestValidator_GTMustBeFirstFormatField(t *testing.T) {
	// Per VCFv4.5 §1.6.2, GT must be the first FORMAT
	// key when it is declared. The check fires in
	// ModeStrict.
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, format string
		wantIssue    bool
	}{
		{"GT first", "GT:DP", false},
		{"GT second", "DP:GT", true},
		{"GT alone", "GT", false},
		{"DP without GT", "DP", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			sample := "0/1:30"
			if c.name == "DP without GT" {
				sample = "30"
			}
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC\t.\t.\t.\t" + c.format + "\t" + sample + "\n"))
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
				if iss.Field == "FORMAT" && strings.Contains(iss.Msg, "GT must be the first") {
					found = true
					break
				}
			}
			if found != c.wantIssue {
				t.Errorf("format=%q: found=%v, want %v; issues=%+v", c.format, found, c.wantIssue, issues)
			}
		})
	}
}
