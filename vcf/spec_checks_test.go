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
