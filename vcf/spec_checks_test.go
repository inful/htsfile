package vcf

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/inful/htsfile/fasta"
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

func TestValidator_FullRefMatchesReference(t *testing.T) {
	// Per VCFv4.5 §1.6.1.1: REF must match the reference
	// base by base. The full-REF check (v0.11.0) replaces
	// the older single-base padding check and strictly
	// subsumes it.
	refText := ">chr1\nACGTACGTACGTACGT\n"
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(refText), 0o600); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	ref, err := fasta.Open(refPath)
	if err != nil {
		t.Fatalf("fasta.Open: %v", err)
	}
	defer func() { _ = ref.Close() }()

	hdr := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=DEL,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, line string
		want       bool
	}{
		// SNV at pos 3: REF="G" matches ref[3]='G'.
		{"SNV ok", "chr1\t3\t.\tG\tA\t.\t.\t.", false},
		// SNV with wrong REF: ref[3]='G' but REF="T".
		// Fires — the full check catches this where the
		// old padding check would not.
		{"SNV REF mismatch", "chr1\t3\t.\tT\tA\t.\t.\t.", true},
		// Insertion at pos 3: REF="G" matches ref[3]='G'.
		{"Insertion right", "chr1\t3\t.\tG\tGX\t.\t.\t.", false},
		// Insertion with wrong padding: REF="A" but
		// ref[3]='G'.
		{"Insertion wrong", "chr1\t3\t.\tA\tAX\t.\t.\t.", true},
		// Deletion of pos 4..5: REF="TA", REF[0]='T'
		// matches ref[4]='T' and REF[1]='A' matches
		// ref[5]='A'.
		{"Deletion ok", "chr1\t4\t.\tTA\tT\t.\t.\t.", false},
		// Deletion with one base wrong: REF="TA" but
		// REF[1] should be 'A'. ref[5]='A' so this
		// passes. Use REF="TG" to make REF[1]='G'
		// mismatch ref[5]='A'.
		{"Deletion REF mismatch", "chr1\t4\t.\tTG\tT\t.\t.\t.", true},
		// Multi-base substitution at pos 3..4: REF="GT"
		// matches ref[3..5]='GT'.
		{"Substitution ok", "chr1\t3\t.\tGT\tAC\t.\t.\t.", false},
		// Wrong contig: warning, not error.
		{"Unknown contig", "chrX\t3\t.\tG\tA\t.\t.\t.", false},
		// Symbolic <DEL> at pos 2 with wrong padding.
		{"Symbolic REF mismatch", "chr1\t2\t.\tAG\t<DEL>\t.\t.\t.", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			v := NewValidator(ModeStrict).WithReference(ref)
			issues := v.Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "does not match reference base") {
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

func TestValidator_NoReference_NoCheck(t *testing.T) {
	// Without a reference attached, no full-REF check fires.
	hdr := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	dec := NewDecoder(strings.NewReader(hdr + "chr1\t2\t.\tC\tT\t.\t.\t.\n"))
	h, _ := dec.Header()
	rec, _ := dec.Record()
	issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
	for _, iss := range issues {
		if strings.Contains(iss.Msg, "does not match reference base") {
			t.Errorf("unexpected full-REF check without reference: %+v", iss)
		}
	}
}

func TestValidator_ENDBoundsWithReference(t *testing.T) {
	// With a FASTA reference, END must be within the
	// contig length even when no ##contig=<length>
	// declaration is present.
	refText := ">chr1\nACGTACGTACGTACGT\n" // 16 bases
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(refText), 0o600); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	ref, err := fasta.Open(refPath)
	if err != nil {
		t.Fatalf("fasta.Open: %v", err)
	}
	defer func() { _ = ref.Close() }()

	hdr := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=DEL,Description=\"x\">\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n" +
		"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"s\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, line string
		want       bool
	}{
		// END = 6 (within contig). No bounds error from
		// the reference check. (Other checks may fire.)
		{"END within contig", "chr1\t3\t.\tAC\t<DEL>\t.\t.\tEND=6;SVLEN=4", false},
		// END = 100 (past contig). Bounds error fires.
		{"END past contig", "chr1\t3\t.\tAC\t<DEL>\t.\t.\tEND=100;SVLEN=98", true},
		// END = 2, POS = 3. END < POS. Bounds error fires.
		{"END before POS", "chr1\t3\t.\tAC\t<DEL>\t.\t.\tEND=2;SVLEN=-1", true},
		// END = 16 (last base). Within contig. No error.
		{"END at last base", "chr1\t3\t.\tAC\t<DEL>\t.\t.\tEND=16;SVLEN=14", false},
		// No END declared: skip.
		{"no END", "chr1\t3\t.\tAC\t<DEL>\t.\t.\tSVLEN=4", false},
		// Mitochondrial: skipped even when END is past.
		{"mtDNA skipped", "chrM\t3\t.\tAC\t<DEL>\t.\t.\tEND=100;SVLEN=98", false},
		// Unknown contig: skipped (Length lookup fails).
		{"unknown contig", "chrX\t3\t.\tAC\t<DEL>\t.\t.\tEND=100;SVLEN=98", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			v := NewValidator(ModeStrict).WithReference(ref)
			issues := v.Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "exceeds contig length") || strings.Contains(iss.Msg, "before POS") {
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

func TestValidator_BreakendPositionsWithReference(t *testing.T) {
	// Per VCFv4.5 §1.6.5.1, a breakend replacement string
	// has the form t]p], t[p[, ]p]t, or [p[t, where p
	// names a remote contig and position. With a FASTA
	// reference attached, the remote position must be
	// within the remote contig's length.
	refText := ">chr1\nACGTACGTACGTACGT\n>chr2\nGATTACAGATTACA\n" // 16 + 14 bases
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(refText), 0o600); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	ref, err := fasta.Open(refPath)
	if err != nil {
		t.Fatalf("fasta.Open: %v", err)
	}
	defer func() { _ = ref.Close() }()

	hdr := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, line string
		want       bool
	}{
		// Breakend with remote position within chr1
		// (length 16). Forms: t]p] and t[p[.
		{"t]p] within contig", "chr1\t3\t.\tA\tA]chr1:10]\t.\t.\t.", false},
		{"t[p[ within contig", "chr1\t3\t.\tA\tA[chr1:10[\t.\t.\t.", false},
		// Cross-contig breakend within chr2 (length 14).
		{"]p]t within remote", "chr1\t3\t.\tA\t]chr2:5]A\t.\t.\t.", false},
		{"[p[t within remote", "chr1\t3\t.\tA\t[chr2:5[A\t.\t.\t.", false},
		// Remote position past contig length.
		{"remote past contig", "chr1\t3\t.\tA\tA]chr1:100]\t.\t.\t.", true},
		// Remote position past chr2 length (chr2 = 14).
		{"remote past chr2", "chr1\t3\t.\tA\tA]chr2:100]\t.\t.\t.", true},
		// Unknown remote contig: warning in ModeDefault,
		// promoted to error in ModeStrict.
		{"unknown remote contig", "chr1\t3\t.\tA\tA]chrX:100]\t.\t.\t.", true},
		// Telomere breakpoint (no position): skip.
		{"telomere (no position)", "chr1\t3\t.\tA\t]chr1]A\t.\t.\t.", false},
		// Non-breakend ALT: skip.
		{"non-breakend ALT", "chr1\t3\t.\tA\tT\t.\t.\t.", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			v := NewValidator(ModeStrict).WithReference(ref)
			issues := v.Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "exceeds contig length") ||
					strings.Contains(iss.Msg, "unknown contig") {
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

func TestValidator_SVLENBoundsWithReference(t *testing.T) {
	// With a FASTA reference, POS + |SVLEN| - 1 must be
	// within the contig.
	refText := ">chr1\nACGTACGTACGTACGT\n" // 16 bases
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(refText), 0o600); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	ref, err := fasta.Open(refPath)
	if err != nil {
		t.Fatalf("fasta.Open: %v", err)
	}
	defer func() { _ = ref.Close() }()

	hdr := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=DEL,Description=\"x\">\n" +
		"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"s\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, line string
		want       bool
	}{
		// POS=3, SVLEN=4 → end=6 (within contig).
		{"SVLEN within contig", "chr1\t3\t.\tAC\t<DEL>\t.\t.\tSVLEN=4", false},
		// POS=3, SVLEN=100 → end=102 (past contig).
		{"SVLEN past contig", "chr1\t3\t.\tAC\t<DEL>\t.\t.\tSVLEN=100", true},
		// Negative SVLEN is rejected by the existing
		// positive-value check (T2.6) before the bounds
		// check fires — both checks agree it's invalid.
		// (No separate sub-test for negative here.)
		// POS=3, SVLEN=14 → end=16 (last base). No error.
		{"SVLEN at last base", "chr1\t3\t.\tAC\t<DEL>\t.\t.\tSVLEN=14", false},
		// POS=1, SVLEN=16 → end=16 (last base). No error.
		{"POS=1 SVLEN=length", "chr1\t1\t.\tAC\t<DEL>\t.\t.\tSVLEN=16", false},
		// No SVLEN: skip.
		{"no SVLEN", "chr1\t3\t.\tAC\t<DEL>\t.\t.\t.", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			v := NewValidator(ModeStrict).WithReference(ref)
			issues := v.Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "implies end position") {
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

func TestValidator_CIGARConsistency(t *testing.T) {
	// Per VCFv4.5 §1.6.1.1: CIGAR describes the
	// alignment between REF and ALT; the sum of REF-
	// consuming op lengths must equal len(REF) and the
	// sum of ALT-consuming op lengths must equal len(ALT).
	hdr := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"c\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, line string
		want       bool
	}{
		// SNV: REF=A, ALT=T, CIGAR=1M. refLen=1, altLen=1.
		{"SNV ok", "chr1\t3\t.\tA\tT\t.\t.\tCIGAR=1M", false},
		// SNV with wrong REF length.
		{"SNV wrong REF length", "chr1\t3\t.\tA\tT\t.\t.\tCIGAR=2M", true},
		// Deletion: REF=ACGT, ALT=A. CIGAR=1M3D.
		{"Deletion ok", "chr1\t4\t.\tACGT\tA\t.\t.\tCIGAR=1M3D", false},
		// Insertion: REF=A, ALT=ACGT. CIGAR=1M3I.
		{"Insertion ok", "chr1\t3\t.\tA\tACGT\t.\t.\tCIGAR=1M3I", false},
		// Insertion wrong ALT length.
		{"Insertion wrong ALT", "chr1\t3\t.\tA\tACGT\t.\t.\tCIGAR=1M2I", true},
		// Malformed CIGAR (missing op).
		{"Malformed CIGAR", "chr1\t3\t.\tA\tT\t.\t.\tCIGAR=1", true},
		// Malformed CIGAR (unknown op).
		{"Unknown CIGAR op", "chr1\t3\t.\tA\tT\t.\t.\tCIGAR=1Z", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "CIGAR") {
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

func TestValidator_CIRUCMissingMatchesRUC(t *testing.T) {
	// Per VCFv4.5 §1.6.1.7: "CIRUC/CIRB must not be
	// non-missing for any alleles with no corresponding
	// RUC/RB value." When RUC[i] is missing (raw "."),
	// CIRUC[2i] and CIRUC[2i+1] must also be ".".
	hdr := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=CNV:TR,Description=\"x\">\n" +
		"##INFO=<ID=RN,Number=A,Type=Integer,Description=\"x\">\n" +
		"##INFO=<ID=RUC,Number=.,Type=Float,Description=\"x\">\n" +
		"##INFO=<ID=RB,Number=.,Type=Integer,Description=\"x\">\n" +
		"##INFO=<ID=CIRUC,Number=.,Type=Float,Description=\"x\">\n" +
		"##INFO=<ID=CIRB,Number=.,Type=Integer,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		{"RUC=., CIRUC=.,. ok", "RN=2;RUC=.,5;CIRUC=.,.,0,2", false},
		{"RUC=., CIRUC[1]=0 (non-missing)", "RN=2;RUC=.,5;CIRUC=.,0,0,2", true},
		{"RUC=., CIRUC[0]=0 (non-missing)", "RN=2;RUC=.,5;CIRUC=0,.,0,2", true},
		{"RUC all present, CIRUC all missing", "RN=2;RUC=5,5;CIRUC=.,.,.,.", false},
		{"RB missing, CIRB[1] non-missing", "RN=2;RB=.,5;CIRB=.,0,0,5", true},
		{"all RUC present, all CIRUC present", "RN=2;RUC=5,5;CIRUC=0,1,0,2", false},
		{"no CIRUC", "RN=2;RUC=.,5", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + "chr1\t100\t.\tA\t<CNV:TR>\t.\t.\t" + c.info + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "non-missing but the corresponding") {
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

func TestValidator_PSLTrioConsistency(t *testing.T) {
	// Per VCFv4.5 §1.6.2 PSQ: "The missing value '$.$'
	// should be used when the corresponding PSL value is
	// missing." We enforce PSL[i] missing iff PSO[i]
	// missing iff PSQ[i] missing per allele.
	hdr := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=PSL,Number=P,Type=String,Description=\"l\">\n" +
		"##FORMAT=<ID=PSO,Number=P,Type=Integer,Description=\"o\">\n" +
		"##FORMAT=<ID=PSQ,Number=P,Type=Integer,Description=\"q\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, gt, psl, pso, psq string
		want                    bool
	}{
		{"all consistent", "|0|1", "chr1*5*1,chr1*5*2", "1,2", "10,20", false},
		{"PSL miss, PSO has value", "|0|1", ".,chr1*5*2", "1,2", "10,20", true},
		{"PSL has value, PSO miss", "|0|1", "chr1*5*1,chr1*5*2", ".,2", "10,20", true},
		{"PSL has value, PSQ miss", "|0|1", "chr1*5*1,chr1*5*2", "1,2", ".,20", true},
		{"all missing", "|0|1", ".,.", ".,.", ".,.", false},
		{"no PSO/PSQ", "|0|1", "chr1*5*1,chr1*5*2", "", "", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			line := "chr1\t100\t.\tA\tC,G\t.\t.\t.\tGT:PSL"
			sample := c.gt + ":" + c.psl
			if c.pso != "" {
				line += ":PSO"
				sample += ":" + c.pso
			}
			if c.psq != "" {
				line += ":PSQ"
				sample += ":" + c.psq
			}
			dec := NewDecoder(strings.NewReader(hdr + line + "\t" + sample + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "missing") && (strings.Contains(iss.Msg, "PSL[") || strings.Contains(iss.Msg, "PSO[") || strings.Contains(iss.Msg, "PSQ[")) {
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

func TestValidator_MKeyValueRange(t *testing.T) {
	// Per VCFv4.5 §1.6.1.5: "Values must be between 0 and 1".
	hdr := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=M5mC,Number=.,Type=Float,Description=\"x\">\n" +
		"##INFO=<ID=M5hmC,Number=.,Type=Float,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		{"in range", "M5mC=0.5", false},
		{"at boundary 0", "M5mC=0", false},
		{"at boundary 1", "M5mC=1", false},
		{"above 1", "M5mC=1.5", true},
		{"below 0", "M5mC=-0.1", true},
		{"missing ok", "M5mC=.", false},
		{"mixed (one bad)", "M5mC=0.5,2", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + "chr1\t100\t.\tA\tC\t.\t.\t" + c.info + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "M-key values must be between 0 and 1") {
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

func TestValidator_GPRange(t *testing.T) {
	// Per VCFv4.5 §1.6.2 GP field: "Genotype posterior
	// probabilities in the range 0 to 1".
	hdr := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"g\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, gt, gp string
		want         bool
	}{
		{"in range", "0/1", "0.9,0.05,0.05", false},
		{"above 1", "0/1", "1.5,0.05,0.05", true},
		{"below 0", "0/1", "-0.1,0.5,0.5", true},
		{"missing ok", "0/1", ".,.,.", false},
		{"no GP ok", "0/1", ".,.,.", false}, // no GP, no error
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + "chr1\t100\t.\tA\tC\t.\t.\t.\tGT:GP\t" + c.gt + ":" + c.gp + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "GP values must be between 0 and 1") {
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

func TestValidator_MRequiresGT(t *testing.T) {
	// Per VCFv4.5 §1.6.1.5: "If any base modification key
	// is present for a sample, GT must be defined for that
	// sample."
	hdr := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=M5mC,Number=1,Type=Float,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, sample string
		want         bool
	}{
		{"GT present, M present", "0/1:0.5", false},
		{"GT present, M missing", "0/1:.", false},
		{"GT missing, M present", ".:0.5", true},
		{"GT missing, M missing", ".:.", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + "chr1\t100\t.\tA\tC\t.\t.\t.\tGT:M5mC\t" + c.sample + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "requires GT") {
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

func TestValidator_INFOGeneralCardinality(t *testing.T) {
	// Per VCFv4.5 §1.6.1: each INFO key with Number=A or
	// Number=R must carry one value per declared allele
	// (A = each ALT, R = REF + each ALT).
	hdr := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">\n" +
		"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, info string
		want       bool
	}{
		{"Number=A correct", "AF=0.5", false},
		{"Number=A wrong count", "AF=0.5,0.3", true},
		{"Number=R correct", "AD=10,20", false},
		{"Number=R wrong count", "AD=10,20,30", true},
		{"Number=R missing", ".", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + "chr1\t100\t.\tA\tC\t.\t.\t" + c.info + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "expected") && strings.Contains(iss.Msg, "value(s) for Number=") {
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

func TestValidator_METAValuesEnforced(t *testing.T) {
	// Per VCFv4.5 §1.4.7 "Sample field format": when a
	// ##META entry declares a Values list, every
	// ##SAMPLE field with that key must take a value
	// from the list. Fires as a SeverityWarning in
	// ModeStrict; skipped in ModeDefault.
	cases := []struct {
		name string
		samp string
		want bool
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

func TestValidator_ImpreciseRequiresCIPOSForPreciseAlt(t *testing.T) {
	// Per VCFv4.5 §1.6.1.7: "If a precise ALT allele
	// is present in a record with the IMPRECISE flag,
	// CIPOS must be explicitly set for that allele, even
	// if it is 0,0."
	hdr := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"x\">\n" +
		"##INFO=<ID=CIPOS,Number=.,Type=Integer,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, line string
		want       bool
	}{
		{"imp imprecise + cipos present", "chr1\t100\t.\tA\tC\t.\t.\tIMPRECISE;CIPOS=0,0", false},
		{"imp imprecise + cipos missing", "chr1\t100\t.\tA\tC\t.\t.\tIMPRECISE", true},
		{"no imprecise flag + no cipos", "chr1\t100\t.\tA\tC\t.\t.\t.", false},
		{"breakend allele + imprecise", "chr1\t100\t.\tA\tG]chr1:200]\t.\t.\tIMPRECISE", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "IMPRECISE flag is set but CIPOS is missing") {
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

func TestValidator_SVLENRequiredForSymbolicAlt(t *testing.T) {
	// Per VCFv4.5 §1.6.1.7: "SVLEN must be specified
	// for symbolic structural variant alleles."
	hdr := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=DEL,Description=\"x\">\n" +
		"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	cases := []struct {
		name, alt, info string
		want            bool
	}{
		{"svlen present", "<DEL>", "SVLEN=100", false},
		{"svlen missing (only 1 entry)", "<DEL>,<DUP>", "SVLEN=100", true},
		{"svlen missing entirely", "<DEL>", ".", true},
		{"non-symbolic alt", "C", ".", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + "chr1\t100\t.\tA\t" + c.alt + "\t.\t.\t" + c.info + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "requires SVLEN to be specified") {
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

func TestValidator_LAAPrecedesLocalAlleleKeys(t *testing.T) {
	// Per VCFv4.5 §1.6.2: "LAA must also be present and
	// precede all fields other than GT."
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"l\">\n" +
		"##FORMAT=<ID=LAD,Number=.,Type=Integer,Description=\"l\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, format, sample string
		want                 bool
	}{
		{"GT LAA LAD", "GT:LAA:LAD", "0/1:2:5,10", false},
		{"GT LAD before LAD", "GT:LAA:LAD", "0/1:2:5,10", false},
		{"LAD before LAA", "GT:LAD:LAA", "0/1:5,10:2", true},
		{"GT DP LAA LAD", "GT:DP:LAA:LAD", "0/1:30:2:5,10", false},
		{"GT LAA DP LAD (non-local allele after LAA is fine)", "GT:LAA:DP:LAD", "0/1:2:30:5,10", false},
		{"LAA missing", "GT:LAD", "0/1:5,10", false}, // caught by T2.1 already
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC,G\t.\t.\t.\t" + c.format + "\t" + c.sample + "\n"))
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
				if strings.Contains(iss.Msg, "local-allele field") && strings.Contains(iss.Msg, "must appear after LAA") {
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

func TestValidator_PSLUnphasedMustBeDot(t *testing.T) {
	// Per VCFv4.5 §1.6.2: "Unphased alleles (without a |
	// separator before them) must have the value '$.$' in
	// their corresponding position in the list."
	hdrText := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"##FORMAT=<ID=PSL,Number=P,Type=String,Description=\"l\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n"
	cases := []struct {
		name, gt, psl string
		want          bool
	}{
		{"all unphased ok", "0/1", ".,.", false},
		{"unphased wrong", "0/1", "chr1*100*1,chr1*200*1", true},
		{"all phased ok", "|0|1", "chr1*100*1,chr1*200*1", false},
		{"first allele (no separator) ok", "0|1", ".,chr1*200*1", false},
		{"first allele (no separator) wrong", "0|1", "chr1*100*1,chr1*200*1", true},
		{"one of two unphased wrong", "0|1", "chr1*100*1,chr1*200*1", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdrText + "chr1\t100\t.\tA\tC,G\t.\t.\t.\tGT:PSL\t" + c.gt + ":" + c.psl + "\n"))
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
				if strings.Contains(iss.Msg, "unphased") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("%s (%s/%s): found=%v, want %v; issues=%+v", c.name, c.gt, c.psl, found, c.want, issues)
			}
		})
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
		{"matches computed", "END=103", false}, // POS=100, REF=ACGT => END=103
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

// TestValidator_CIRUCBIRBBoundsSpanZero covers
// VCFv4.5 §1.6.1.7: "The lower bound must be less
// than or equal to zero and the upper bound must be
// greater than or equal to zero." Fires for both
// CIRUC (paired with RUC) and CIRB (paired with RB).
func TestValidator_CIRUCBIRBBoundsSpanZero(t *testing.T) {
	hdr := "##fileformat=VCFv4.2\n" +
		"##ALT=<ID=CNV:TR,Description=\"tr\">\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n" +
		"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"s\">\n" +
		"##INFO=<ID=CN,Number=A,Type=Float,Description=\"c\">\n" +
		"##INFO=<ID=RN,Number=A,Type=Integer,Description=\"r\">\n" +
		"##INFO=<ID=RUS,Number=.,Type=String,Description=\"us\">\n" +
		"##INFO=<ID=RUL,Number=.,Type=Integer,Description=\"ul\">\n" +
		"##INFO=<ID=RUC,Number=.,Type=Float,Description=\"uc\">\n" +
		"##INFO=<ID=RB,Number=.,Type=Integer,Description=\"rb\">\n" +
		"##INFO=<ID=CIRUC,Number=.,Type=Float,Description=\"ciruc\">\n" +
		"##INFO=<ID=CIRB,Number=.,Type=Integer,Description=\"cirb\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

	cases := []struct {
		name, line string
		want       bool // true = expect an issue
	}{
		{
			// RUC=10, CIRUC lower=-5, upper=5: spans 0,
			// valid.
			name: "valid CIRUC",
			line: "chr1\t100\t.\tA\t<CNV:TR>\t.\t.\tSVLEN=1;RN=1;RUS=CAG;RUL=3;RUC=10;CIRUC=-5,5\n",
			want: false,
		},
		{
			// CIRUC lower=+1 (must be ≤ 0): error.
			name: "CIRUC lower bound positive",
			line: "chr1\t100\t.\tA\t<CNV:TR>\t.\t.\tSVLEN=1;RN=1;RUS=CAG;RUL=3;RUC=10;CIRUC=1,5\n",
			want: true,
		},
		{
			// CIRUC upper=-1 (must be ≥ 0): error.
			name: "CIRUC upper bound negative",
			line: "chr1\t100\t.\tA\t<CNV:TR>\t.\t.\tSVLEN=1;RN=1;RUS=CAG;RUL=3;RUC=10;CIRUC=-5,-1\n",
			want: true,
		},
		{
			// CIRB lower=+2: error.
			name: "CIRB lower bound positive",
			line: "chr1\t100\t.\tA\t<CNV:TR>\t.\t.\tSVLEN=1;RN=1;RUS=CAG;RUL=3;RB=30;CIRB=2,5\n",
			want: true,
		},
		{
			// CIRB upper=-2: error.
			name: "CIRB upper bound negative",
			line: "chr1\t100\t.\tA\t<CNV:TR>\t.\t.\tSVLEN=1;RN=1;RUS=CAG;RUL=3;RB=30;CIRB=-5,-2\n",
			want: true,
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			// Other checks (RN cardinality, CIRUC
			// missingness) may also fire; we want the
			// bounds-span issue specifically.
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "must be ≤ 0") ||
					strings.Contains(iss.Msg, "must be ≥ 0") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("found=%v, want %v; issues=%+v", found, c.want, issues)
			}
		})
	}
}

// TestValidator_GVCFBlockRules covers §1.7 gVCF
// non-variant block validation. The four rules are:
//
//  1. <*>-ALT must have END in INFO.
//  2. END ≥ POS.
//  3. No overlap with the previous record (POS >
//     prevEnd when the previous record was a block).
//  4. (With a FASTA reference) REF must match the
//     reference over the [POS, END] block range.
//
// Cases run both with and without a reference, in
// both ModeDefault and ModeStrict where the rule
// differs.
func TestValidator_GVCFBlockRules(t *testing.T) {
	hdr := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

	type wantCheck struct {
		fieldSubstr string // substring of Field+Msg we expect to find
		mode        Mode
	}
	cases := []struct {
		name, line string
		wants      []wantCheck
	}{
		{
			// Valid block: ALT="<*>" with END=200.
			// No errors expected.
			name:  "valid block passes",
			line:  "chr1\t100\t.\tA\t<*>\t.\t.\tEND=200",
			wants: nil,
		},
		{
			// <NON_REF> is the parser-canonicalised form
			// (T1.6). It is rewritten to <*>.
			name:  "NON_REF normalised",
			line:  "chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200",
			wants: nil,
		},
		{
			// <*>-ALT with no END: error in ModeStrict,
			// warning in ModeDefault.
			name: "missing END strict",
			line: "chr1\t100\t.\tA\t<*>\t.\t.\t.",
			wants: []wantCheck{
				{mode: ModeStrict, fieldSubstr: `ALT="<*>"`},
				{mode: ModeDefault, fieldSubstr: `ALT="<*>"`},
			},
		},
		{
			// END < POS: error in both modes (the END
			// < POS rule is unambiguous).
			name: "END before POS",
			line: "chr1\t100\t.\tA\t<*>\t.\t.\tEND=50",
			wants: []wantCheck{
				{mode: ModeStrict, fieldSubstr: "before POS"},
				{mode: ModeDefault, fieldSubstr: "before POS"},
			},
		},
		{
			// ALT=<*> alongside a real ALT: NOT a
			// non-variant block (the spec is strict
			// about "ALT contains exactly <*>"). The
			// checkGVCFBlock short-circuits — no
			// gVCF-related issues. The existing
			// checkALT / checkHaplotag logic still
			// runs and may flag unrelated things,
			// but for this case it's a normal
			// variant record.
			name:  "non-star ALT ignored by gVCF",
			line:  "chr1\t100\t.\tA\tG,<*>\t.\t.\tEND=200",
			wants: nil,
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			for _, w := range c.wants {
				t.Run(fmt.Sprintf("mode=%d", w.mode), func(t *testing.T) {
					v := NewValidator(w.mode)
					issues := v.Validate(h, []*Record{rec})
					if !hasIssueWithField(issues, w.fieldSubstr) {
						t.Errorf("expected issue containing %q; got %+v", w.fieldSubstr, issues)
					}
				})
			}
		})
	}
}

// TestValidator_GVCFBlockOverlap covers the cross-
// record overlap rule: a non-variant block followed
// by another record whose POS lands inside the block
// is an error in ModeStrict.
func TestValidator_GVCFBlockOverlap(t *testing.T) {
	hdr := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	lines := []string{
		"chr1\t100\t.\tA\t<*>\t.\t.\tEND=200", // block 100..200
		"chr1\t150\t.\tA\tG\t.\t.\t.",         // overlap
		"chr1\t201\t.\tA\tG\t.\t.\t.",         // OK
	}
	dec := NewDecoder(strings.NewReader(hdr + strings.Join(lines, "\n") + "\n"))
	h, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	var recs []*Record
	for {
		rec, err := dec.Record()
		if err != nil {
			break
		}
		recs = append(recs, rec)
	}
	v := NewValidator(ModeStrict)
	issues := v.Validate(h, recs)
	if !hasIssueWithField(issues, "overlaps previous non-variant block") {
		t.Errorf("expected overlap issue; got %+v", issues)
	}
}

// TestValidator_GVCFBlockRefMatchesReference covers the
// reference-aware rule 4: when a FASTA reference is
// attached, REF must match the reference sequence
// over [POS, END] for non-variant blocks.
func TestValidator_GVCFBlockRefMatchesReference(t *testing.T) {
	refText := ">chr1\nACGTACGTACGTACGTACGT\n" // 20 A/C/G/T bases
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(refText), 0o600); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	ref, err := fasta.Open(refPath)
	if err != nil {
		t.Fatalf("fasta.Open: %v", err)
	}
	defer func() { _ = ref.Close() }()

	hdr := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	// Reference bases 1..20 are A C G T A C G T A C G T A C G T A C G T.
	// Block 5..9 covers reference bases 5..9 = A C G T A (5 bases),
	// so REF must be ACGTA.
	goodLine := "chr1\t5\t.\tACGTA\t<*>\t.\t.\tEND=9"
	// Mismatched REF (G instead of A at position 5).
	badLine := "chr1\t5\t.\tGCGTA\t<*>\t.\t.\tEND=9"

	dec := NewDecoder(strings.NewReader(hdr + goodLine + "\n"))
	h, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	v := NewValidator(ModeStrict).WithReference(ref)
	issues := v.Validate(h, []*Record{rec})
	if hasIssueWithField(issues, "does not match reference base") {
		t.Errorf("good block fired a REF-mismatch issue; got %+v", issues)
	}

	dec2 := NewDecoder(strings.NewReader(hdr + badLine + "\n"))
	if _, err := dec2.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	rec2, err := dec2.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	v.Reset()
	issues = v.Validate(h, []*Record{rec2})
	if !hasIssueWithField(issues, "does not match reference base") {
		t.Errorf("bad block did not fire a REF-mismatch issue; got %+v", issues)
	}
}

// hasIssueWithField reports whether any issue in the
// slice contains substr in its Field+Msg concatenated
// form. The substring match is case-sensitive and
// matches across Field and Msg concatenation.
func hasIssueWithField(issues []ValidationError, substr string) bool {
	for _, iss := range issues {
		if strings.Contains(iss.Field+": "+iss.Msg, substr) {
			return true
		}
	}
	return false
}

// TestValidator_BreakendMatePositionRoundTrip
// covers §1.6.5.1: the two records in a breakend pair
// must round-trip — each mate's breakend remote
// position must equal the partner's local (Chrom, POS).
//
// Cases:
//   - Two records forming a balanced mate pair (round-trip).
//   - One mate's POS wrong (round-trip broken).
//   - One mate's breakend remote position wrong.
//   - One mate's breakend pointing to a third contig.
//   - Records without breakends (e.g. SNV with MATEID) skipped.
//   - Mate pair missing MATEID: each record evaluated alone.
func TestValidator_BreakendMatePositionRoundTrip(t *testing.T) {
	hdr := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=MATEID,Number=A,Type=String,Description=\"ID of mate breakend\">\n" +
		"##contig=<ID=chrA,length=10000>\n" +
		"##contig=<ID=chrB,length=10000>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

	// Build a balanced pair: A at chrA:100 with breakend
	// pointing to chrA:200, B at chrA:200 with breakend
	// pointing to chrA:100.
	lines := []struct {
		name      string
		aLine     string
		bLine     string
		wantIssue string // empty = no issue expected
	}{
		{
			name:      "balanced round-trip",
			aLine:     "chrA\t100\tbndA\tA\tA]chrA:200]\t.\t.\tMATEID=bndB\n",
			bLine:     "chrA\t200\tbndB\tA\t]chrA:100]A\t.\t.\tMATEID=bndA\n",
			wantIssue: "",
		},
		{
			name:      "second mate POS wrong",
			aLine:     "chrA\t100\tbndA\tA\tA]chrA:200]\t.\t.\tMATEID=bndB\n",
			bLine:     "chrA\t201\tbndB\tA\t]chrA:100]A\t.\t.\tMATEID=bndA\n",
			wantIssue: "points to chrA:200 but partner is at chrA:201",
		},
		{
			name:      "second mate breakend remote wrong",
			aLine:     "chrA\t100\tbndA\tA\tA]chrA:200]\t.\t.\tMATEID=bndB\n",
			bLine:     "chrA\t200\tbndB\tA\t]chrA:999]A\t.\t.\tMATEID=bndA\n",
			wantIssue: "points to chrA:999 but partner is at chrA:100",
		},
		{
			name:  "second mate breakend wrong contig",
			aLine: "chrA\t100\tbndA\tA\tA]chrB:200]\t.\t.\tMATEID=bndB\n",
			bLine: "chrA\t200\tbndB\tA\t]chrB:100]A\t.\t.\tMATEID=bndA\n",
			// First mate's remote (chrB:200) doesn't match partner's
			// local (chrA:200): the partner should round-trip
			// back to chrA:100. Same for chrB.
			wantIssue: "partner is at",
		},
		{
			name:  "missing MATEID on first mate is skipped",
			aLine: "chrA\t100\tbndA\tA\tA]chrA:200]\t.\t.\t.\n",
			bLine: "chrA\t200\tbndB\tA\t]chrA:100]A\t.\t.\tMATEID=bndA\n",
			// Without MATEID, no round-trip check fires; both
			// records validated independently.
			wantIssue: "",
		},
		{
			name:  "non-breakend ALT with MATEID skipped",
			aLine: "chrA\t100\tbndA\tA\tG\t.\t.\tMATEID=bndB\n",
			bLine: "chrA\t200\tbndB\tA\tC\t.\t.\tMATEID=bndA\n",
			// Records without breakend ALTs are paired SNVs;
			// no position round-trip applies.
			wantIssue: "",
		},
	}
	for _, c := range lines {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.aLine + c.bLine))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			var recs []*Record
			for {
				rec, err := dec.Record()
				if err != nil {
					break
				}
				recs = append(recs, rec)
			}
			v := NewValidator(ModeStrict)
			issues := v.Validate(h, recs)
			got := hasIssueWithField(issues, c.wantIssue)
			if c.wantIssue == "" {
				if got {
					t.Errorf("unexpected issue found; got %+v", issues)
				}
			} else if !got {
				t.Errorf("expected issue containing %q; got %+v", c.wantIssue, issues)
			}
		})
	}
}

// TestValidator_IntegerBCFRange covers VCFv4.5 §1.4.2:
// the eight values from -2147483648 to -2147483641 are
// reserved by BCF2's integer encoding and therefore
// disallowed in both VCF and BCF. Rejecting them at
// parse-time prevents a round-trip VCF→BCF→VCF from
// silently rewriting the values.
func TestValidator_IntegerBCFRange(t *testing.T) {
	hdr := "##fileformat=VCFv4.5\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n" +
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

	cases := []struct {
		name, line string
		want       bool
	}{
		{
			// DP=2147483647 (math.MaxInt32) is fine.
			name: "MaxInt32 ok",
			line: "chr1\t100\t.\tA\tC\t.\t.\tDP=2147483647\tGT\t0/1\n",
			want: false,
		},
		{
			// -2147483648 (MinInt32) is forbidden.
			name: "MinInt32 in INFO forbidden",
			line: "chr1\t100\t.\tA\tC\t.\t.\tDP=-2147483648\tGT\t0/1\n",
			want: true,
		},
		{
			// -2147483641 (MinInt32+7) is the last
			// forbidden value.
			name: "MinInt32+7 in INFO forbidden",
			line: "chr1\t100\t.\tA\tC\t.\t.\tDP=-2147483641\tGT\t0/1\n",
			want: true,
		},
		{
			// -2147483640 (MinInt32+8) is fine.
			name: "MinInt32+8 ok",
			line: "chr1\t100\t.\tA\tC\t.\t.\tDP=-2147483640\tGT\t0/1\n",
			want: false,
		},
		{
			// FORMAT/GQ with MinInt32 also forbidden.
			name: "MinInt32 in FORMAT forbidden",
			line: "chr1\t100\t.\tA\tC\t.\t.\t.\tGT:GQ\t0/1:-2147483648\n",
			want: true,
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "BCF2 forbidden range") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("found=%v, want %v; issues=%+v", found, c.want, issues)
			}
		})
	}
}

// TestValidator_HaploidGTForMtDNA covers VCFv4.5 §1.6.2:
// "Haploid calls, e.g. on Y, male non-pseudoautosomal
// X, or mitochondria, should be indicated by having
// only one allele value." We only check mtDNA — Y is
// ambiguous because of PAR regions.
func TestValidator_HaploidGTForMtDNA(t *testing.T) {
	hdr := "##fileformat=VCFv4.5\n" +
		"##contig=<ID=chrM,length=16569>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
	cases := []struct {
		name, line string
		want       bool
	}{
		{
			// Haploid GT on chrM: ok.
			name: "mtDNA haploid ok",
			line: "chrM\t100\t.\tA\tC\t.\t.\t.\tGT\t1\n",
			want: false,
		},
		{
			// Diploid GT on chrM: error.
			name: "mtDNA diploid forbidden",
			line: "chrM\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\n",
			want: true,
		},
		{
			// Autosome (chr1) diploid GT: no issue
			// (haploid rule only fires on mtDNA).
			name: "autosome diploid ok",
			line: "chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\n",
			want: false,
		},
		{
			// chrY diploid GT: not checked (Y is
			// ambiguous).
			name: "chrY diploid not checked",
			line: "chrY\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\n",
			want: false,
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(hdr + c.line + "\n"))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			rec, err := dec.Record()
			if err != nil {
				t.Fatalf("Record: %v", err)
			}
			issues := NewValidator(ModeStrict).Validate(h, []*Record{rec})
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "haploid GT") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("found=%v, want %v; issues=%+v", found, c.want, issues)
			}
		})
	}
}

// TestValidator_PedigreeReference checks that PEDIGREE
// Father / Mother / Original references resolve to a
// declared sample or another PEDIGREE ID, per
// VCFv4.5 §"PedigreeInDetail".
func TestValidator_PedigreeReference(t *testing.T) {
	hdrOK := "##fileformat=VCFv4.5\n" +
		"##PEDIGREE=<ID=Child,Father=Dad,Mother=Mom,Original=Germline>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

	cases := []struct {
		name string
		hdr  string
		want bool
	}{
		{
			name: "PEDIGREE Father/Mother/Original all declared",
			hdr: "##fileformat=VCFv4.5\n" +
				"##PEDIGREE=<ID=Child,Father=Dad,Mother=Mom,Original=Germline>\n" +
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
				"chr1\t100\t.\tA\tC\t.\t.\t.\n",
			want: false,
		},
		{
			// Father is not declared as a sample.
			// SampleNames is needed for the check to
			// fire (only then is "unknown sample"
			// meaningful).
			name: "PEDIGREE Father undeclared",
			hdr: "##fileformat=VCFv4.5\n" +
				"##PEDIGREE=<ID=Child,Father=Unknown,Mother=Mom>\n" +
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMom\tChild\n" +
				"chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/1\n",
			want: true,
		},
		{
			// PEDIGREE references another PEDIGREE.
			name: "PEDIGREE chain via Original",
			hdr: "##fileformat=VCFv4.5\n" +
				"##PEDIGREE=<ID=Germline>\n" +
				"##PEDIGREE=<ID=Tumour,Original=Germline>\n" +
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGermline\tTumour\n" +
				"chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/1\n",
			want: false,
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dec := NewDecoder(strings.NewReader(c.hdr))
			h, err := dec.Header()
			if err != nil {
				t.Fatalf("Header: %v", err)
			}
			issues := NewValidator(ModeStrict).checkHeader(h)
			found := false
			for _, iss := range issues {
				if strings.Contains(iss.Msg, "references unknown sample") {
					found = true
					break
				}
			}
			if found != c.want {
				t.Errorf("found=%v, want %v; issues=%+v", found, c.want, issues)
			}
		})
	}
	// Unused — silence linter.
	_ = hdrOK
}
