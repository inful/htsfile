package vcf

import (
	"strings"
	"testing"
)

// reservedProbe builds a minimal header-only VCF with one
// declared INFO or FORMAT line and returns the header-level
// diagnostics the Validator produces for it.
func reservedProbe(t *testing.T, decl string) []ValidationError {
	t.Helper()
	src := "##fileformat=VCFv4.5\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		decl + "\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t100\t.\tA\tC,T\t30\tPASS\t.\tGT:."
	if strings.HasPrefix(decl, "##INFO") {
		src = "##fileformat=VCFv4.5\n" +
			"##contig=<ID=chr1,length=1000>\n" +
			decl + "\n" +
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
			"chr1\t100\t.\tA\tC\t30\tPASS\t.\n"
	}
	dec := NewDecoder(strings.NewReader(src))
	h, err := dec.Header()
	if err != nil {
		t.Fatalf("parse header for %s: %v", decl, err)
	}
	return NewValidator(ModeStrict).Validate(h, nil)
}

func headerIssuesNamed(iss []ValidationError, field string) []ValidationError {
	var out []ValidationError
	for _, e := range iss {
		if e.Field == field {
			out = append(out, e)
		}
	}
	return out
}

// infoDecl / formatDecl build a structured meta-information line.
func infoDecl(id, number, typ string) string {
	return "##INFO=<ID=" + id + ",Number=" + number +
		",Type=" + typ + ",Description=\"probe\">"
}

func formatDecl(id, number, typ string) string {
	return "##FORMAT=<ID=" + id + ",Number=" + number +
		",Type=" + typ + ",Description=\"probe\">"
}

// specReservedInfo is the Reserved INFO keys table, VCFv4.5
// spec lines 460-480, verbatim.
var specReservedInfo = []struct{ id, number, typ string }{
	{"AA", "1", "String"},
	{"AC", "A", "Integer"},
	{"AD", "R", "Integer"},
	{"ADF", "R", "Integer"},
	{"ADR", "R", "Integer"},
	{"AF", "A", "Float"},
	{"AN", "1", "Integer"},
	{"BQ", "1", "Float"},
	{"CIGAR", "A", "String"},
	{"DB", "0", "Flag"},
	{"DP", "1", "Integer"},
	{"END", "1", "Integer"},
	{"H2", "0", "Flag"},
	{"H3", "0", "Flag"},
	{"MQ", "1", "Float"},
	{"MQ0", "1", "Integer"},
	{"NS", "1", "Integer"},
	{"SB", "4", "Integer"},
	{"SOMATIC", "0", "Flag"},
	{"VALIDATED", "0", "Flag"},
	{"1000G", "0", "Flag"},
}

// specReservedFormat is the Reserved genotype keys table,
// VCFv4.5 spec lines 534-596, excluding the M/DPM/ADM
// base-modification families (resolved by MFamily, not by a
// literal table entry).
var specReservedFormat = []struct{ id, number, typ string }{
	{"AD", "R", "Integer"},
	{"ADF", "R", "Integer"},
	{"ADR", "R", "Integer"},
	{"DP", "1", "Integer"},
	{"EC", "A", "Integer"},
	{"LEN", "1", "Integer"},
	{"FT", "1", "String"},
	{"GL", "G", "Float"},
	{"GP", "G", "Float"},
	{"GQ", "1", "Integer"},
	{"GT", "1", "String"},
	{"HQ", "2", "Integer"},
	{"LA", ".", "Integer"},
	{"LAA", ".", "Integer"},
	{"LAD", "LR", "Integer"},
	{"LADF", "LR", "Integer"},
	{"LADR", "LR", "Integer"},
	{"LEC", "LA", "Integer"},
	{"LGL", "LG", "Float"},
	{"LGP", "LG", "Float"},
	{"LPL", "LG", "Integer"},
	{"LPP", "LG", "Integer"},
	{"MQ", "1", "Integer"},
	{"PL", "G", "Integer"},
	{"PP", "G", "Integer"},
	{"PQ", "1", "Integer"},
	{"PS", "1", "Integer"},
	{"PSL", "P", "String"},
	{"PSO", "P", "Integer"},
	{"PSQ", "P", "Integer"},
}

// TestReservedInfoTable_MatchesSpecTable locks T4.25: every
// key in the spec's Reserved INFO keys table must be present
// in reservedInfoKeys with exactly the spec's Number/Type.
func TestReservedInfoTable_MatchesSpecTable(t *testing.T) {
	for _, want := range specReservedInfo {
		got, ok := reservedInfoKeys[want.id]
		if !ok {
			t.Errorf("reservedInfoKeys missing %q (spec Number=%s Type=%s)",
				want.id, want.number, want.typ)
			continue
		}
		if got.Number != want.number || got.Type.String() != want.typ {
			t.Errorf("reservedInfoKeys[%q] = {Number:%q Type:%q}, spec says {Number:%q Type:%q}",
				want.id, got.Number, got.Type, want.number, want.typ)
		}
	}
}

// TestReservedFormatTable_MatchesSpecTable locks T4.25 and
// T4.26 for the FORMAT side.
func TestReservedFormatTable_MatchesSpecTable(t *testing.T) {
	for _, want := range specReservedFormat {
		got, ok := reservedFormatKeys[want.id]
		if !ok {
			t.Errorf("reservedFormatKeys missing %q (spec Number=%s Type=%s)",
				want.id, want.number, want.typ)
			continue
		}
		if got.Number != want.number || got.Type.String() != want.typ {
			t.Errorf("reservedFormatKeys[%q] = {Number:%q Type:%q}, spec says {Number:%q Type:%q}",
				want.id, got.Number, got.Type, want.number, want.typ)
		}
	}
}

// TestReservedKeys_NoFalsePositiveOnSpecDeclarations locks
// T4.26: a header that declares a reserved key exactly as the
// spec's table does must produce no reserved-key diagnostic.
func TestReservedKeys_NoFalsePositiveOnSpecDeclarations(t *testing.T) {
	for _, want := range specReservedInfo {
		field := "header/INFO/" + want.id
		if got := headerIssuesNamed(reservedProbe(t, infoDecl(want.id, want.number, want.typ)), field); len(got) != 0 {
			t.Errorf("false positive on spec-correct %s: %+v", field, got)
		}
	}
	for _, want := range specReservedFormat {
		if want.id == "GT" {
			continue // already declared by the probe fixture
		}
		field := "header/FORMAT/" + want.id
		if got := headerIssuesNamed(reservedProbe(t, formatDecl(want.id, want.number, want.typ)), field); len(got) != 0 {
			t.Errorf("false positive on spec-correct %s: %+v", field, got)
		}
	}
}

// TestReservedKeys_WrongDeclarationWarns locks T4.25: a
// reserved key declared with the wrong Type must be reported.
func TestReservedKeys_WrongDeclarationWarns(t *testing.T) {
	cases := []struct {
		decl  string
		field string
	}{
		{infoDecl("DP", "5", "String"), "header/INFO/DP"},
		{infoDecl("SB", "1", "Float"), "header/INFO/SB"},
		{infoDecl("AC", "1", "Float"), "header/INFO/AC"},
		{infoDecl("H3", "1", "Integer"), "header/INFO/H3"},
		{infoDecl("AF", "R", "Integer"), "header/INFO/AF"},
		{infoDecl("1000G", "1", "Integer"), "header/INFO/1000G"},
		{formatDecl("GQ", "2", "Float"), "header/FORMAT/GQ"},
		{formatDecl("HQ", "G", "Float"), "header/FORMAT/HQ"},
		{formatDecl("PS", "A", "String"), "header/FORMAT/PS"},
		{formatDecl("FT", "2", "Integer"), "header/FORMAT/FT"},
		{formatDecl("LEN", "A", "Integer"), "header/FORMAT/LEN"},
	}
	for _, c := range cases {
		if got := headerIssuesNamed(reservedProbe(t, c.decl), c.field); len(got) == 0 {
			t.Errorf("expected a diagnostic for %s (declared %s)", c.field, c.decl)
		}
	}
}

// TestReservedKeys_RunInDefaultMode locks T4.27: the
// reserved-key check must run in ModeDefault too, not only in
// ModeStrict.
func TestReservedKeys_RunInDefaultMode(t *testing.T) {
	src := "##fileformat=VCFv4.5\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		infoDecl("DP", "5", "String") + "\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tC\t30\tPASS\t.\n"
	dec := NewDecoder(strings.NewReader(src))
	h, err := dec.Header()
	if err != nil {
		t.Fatalf("header: %v", err)
	}
	got := headerIssuesNamed(
		NewValidator(ModeDefault).Validate(h, nil), "header/INFO/DP")
	if len(got) == 0 {
		t.Fatal("ModeDefault emitted no reserved-key diagnostic; " +
			"HeaderChecks still gates checkReservedKeys on ModeStrict")
	}
	if got[0].Severity != SeverityWarning {
		t.Errorf("severity = %v, want SeverityWarning", got[0].Severity)
	}
}
