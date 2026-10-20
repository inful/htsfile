package fastq

import "testing"

func TestValidator_Default_EmptyHeader(t *testing.T) {
	v := NewValidator(ModeDefault)
	rec := &Record{Header: "", Seq: []byte("ACGT"), Qual: []byte("!!!!")}
	issues := v.Validate([]*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("got %d issues, want 1: %v", len(issues), issues)
	}
	if issues[0].Field != "header" || issues[0].Severity != SeverityError {
		t.Errorf("issue: %+v", issues[0])
	}
}

func TestValidator_Default_EmptySeq(t *testing.T) {
	v := NewValidator(ModeDefault)
	rec := &Record{Header: "r1", Seq: []byte{}, Qual: []byte{}}
	issues := v.Validate([]*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("got %d issues, want 1: %v", len(issues), issues)
	}
	if issues[0].Field != "sequence" {
		t.Errorf("issue: %+v", issues[0])
	}
}

func TestValidator_Default_SeqQualLengthMismatch(t *testing.T) {
	v := NewValidator(ModeDefault)
	rec := &Record{Header: "r1", Seq: []byte("ACGT"), Qual: []byte("III")}
	issues := v.Validate([]*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("got %d issues, want 1: %v", len(issues), issues)
	}
	if issues[0].Field != "quality" {
		t.Errorf("issue: %+v", issues[0])
	}
}

func TestValidator_Default_QualBelowOffset(t *testing.T) {
	// Phred+33 with a quality byte below '!' (=33) — invalid.
	v := NewValidator(ModeDefault)
	rec := &Record{Header: "r1", Seq: []byte("A"), Qual: []byte(" "), PhredOffset: 33}
	issues := v.Validate([]*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("got %d issues, want 1: %v", len(issues), issues)
	}
	if issues[0].Field != "quality" {
		t.Errorf("issue: %+v", issues[0])
	}
}

func TestValidator_Default_Happy(t *testing.T) {
	v := NewValidator(ModeDefault)
	rec := &Record{Header: "r1", Seq: []byte("ACGT"), Qual: []byte("!!!!"), PhredOffset: 33}
	if issues := v.Validate([]*Record{rec}); len(issues) != 0 {
		t.Errorf("expected no issues, got %v", issues)
	}
}

func TestValidator_Default_MultipleRecords(t *testing.T) {
	v := NewValidator(ModeDefault)
	recs := []*Record{
		{Header: "r1", Seq: []byte("ACGT"), Qual: []byte("!!!!"), PhredOffset: 33},
		{Header: "", Seq: []byte("ACGT"), Qual: []byte("!!!!"), PhredOffset: 33},
		{Header: "r3", Seq: []byte("A"), Qual: []byte("I"), PhredOffset: 33},
	}
	issues := v.Validate(recs)
	if len(issues) != 1 {
		t.Fatalf("got %d issues, want 1: %v", len(issues), issues)
	}
	if issues[0].Line != recs[1].Line {
		t.Errorf("issue line: got %d want %d", issues[0].Line, recs[1].Line)
	}
}

func TestValidator_Strict_Alphabet(t *testing.T) {
	// Strict mode rejects bytes not in the IUPAC nucleotide
	// alphabet (default: ACGTURYSWKMBDHVN, case-insensitive).
	v := NewValidator(ModeStrict)
	rec := &Record{Header: "r1", Seq: []byte("ACGTZ"), Qual: []byte("!!!!!"), PhredOffset: 33}
	issues := v.Validate([]*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("got %d issues, want 1: %v", len(issues), issues)
	}
	if issues[0].Field != "sequence" {
		t.Errorf("issue: %+v", issues[0])
	}
	if issues[0].Severity != SeverityError {
		t.Errorf("severity: got %v want error", issues[0].Severity)
	}
}

func TestValidator_Strict_AmbiguityCodes(t *testing.T) {
	// IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V)
	// and 'N' must pass strict mode.
	v := NewValidator(ModeStrict)
	seq := "ACGTRYSWKMBDHVN"
	rec := &Record{Header: "r1", Seq: []byte(seq), Qual: []byte("!!!!!!!!!!!!!!!"), PhredOffset: 33}
	if issues := v.Validate([]*Record{rec}); len(issues) != 0 {
		t.Errorf("ambiguity codes rejected: %v", issues)
	}
}

func TestValidator_Strict_Lowercase(t *testing.T) {
	v := NewValidator(ModeStrict)
	rec := &Record{Header: "r1", Seq: []byte("acgt"), Qual: []byte("!!!!"), PhredOffset: 33}
	if issues := v.Validate([]*Record{rec}); len(issues) != 0 {
		t.Errorf("lowercase rejected: %v", issues)
	}
}

func TestValidator_Strict_DefaultChecksStillRun(t *testing.T) {
	// Strict mode still surfaces seq/qual length mismatch.
	v := NewValidator(ModeStrict)
	rec := &Record{Header: "r1", Seq: []byte("ACGT"), Qual: []byte("III"), PhredOffset: 33}
	issues := v.Validate([]*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("got %d issues, want 1: %v", len(issues), issues)
	}
	if issues[0].Field != "quality" {
		t.Errorf("issue: %+v", issues[0])
	}
}

func TestValidator_CustomAlphabet(t *testing.T) {
	v := NewValidator(ModeStrict)
	v.Alphabet = "X" // only X is allowed; everything else errors
	rec := &Record{Header: "r1", Seq: []byte("XXXX"), Qual: []byte("!!!!"), PhredOffset: 33}
	if issues := v.Validate([]*Record{rec}); len(issues) != 0 {
		t.Errorf("custom alphabet: %v", issues)
	}
	rec.Seq = []byte("XXXY")
	issues := v.Validate([]*Record{rec})
	if len(issues) != 1 {
		t.Fatalf("custom alphabet: got %d issues want 1: %v", len(issues), issues)
	}
}

func TestSeverity_String(t *testing.T) {
	cases := []struct {
		sev  Severity
		want string
	}{
		{SeverityWarning, "warning"},
		{SeverityError, "error"},
		{Severity(99), "unknown"},
	}
	for _, tc := range cases {
		if got := tc.sev.String(); got != tc.want {
			t.Errorf("String(%v): got %q want %q", tc.sev, got, tc.want)
		}
	}
}
