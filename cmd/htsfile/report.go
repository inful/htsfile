package main

import (
	"fmt"
	"sort"

	"github.com/inful/htsfile/bam"
	"github.com/inful/htsfile/bcf"
	"github.com/inful/htsfile/bed"
	"github.com/inful/htsfile/fasta"
	"github.com/inful/htsfile/fastq"
	"github.com/inful/htsfile/gff"
	"github.com/inful/htsfile/sam"
	"github.com/inful/htsfile/vcf"
)

// Shared issue reporter (extracted from the 5-way duplicated
// runCheck{VCF,BED,SAM,BAM,GFF} functions).
//
// Each per-format package has its own Issue type with its own
// Severity constants. To avoid forcing every package to share a
// common Issue shape (which would couple unrelated packages), the
// reporter works on a normalised reportIssue struct and each
// runCheck* function converts its native Issue slice to
// []reportIssue before calling runIssueReport.

// reportSeverity is the unified severity used by the shared
// reporter. Each per-package SeverityError / SeverityWarning
// constant maps to one of these via the reportFrom* adapters.
type reportSeverity int

const (
	reportWarning reportSeverity = iota
	reportError
)

// String returns the lowercase name of the severity — matches
// the strings used by every per-package Severity.String() so
// the printed output is unchanged.
func (s reportSeverity) String() string {
	switch s {
	case reportWarning:
		return "warning"
	case reportError:
		return "error"
	}
	return "unknown"
}

// reportIssue is the normalised shape consumed by the shared
// reporter. Record is the 1-based record index within the
// alignment section; 0 means "no record index" (e.g. for
// VCF issues that have a record number, BED/SAM/BAM/GFF that
// don't).
type reportIssue struct {
	Line     int
	Record   int
	Field    string
	Severity reportSeverity
	Msg      string
}

// runIssueReport sorts issues, prints them, counts severities,
// and returns an error if any errors were found (or warnings,
// if exitOnWarnings is true). The output format matches the
// pre-refactor per-package inline loop exactly:
//
//	line N        [error] field: message
//	line N rec M   [error] field: message  (when Record > 0)
//
//	<path>: N error(s), M warning(s)
//
// followed by exit code 1 if errors (or warnings with
// --warnings) were found.
func runIssueReport(path string, issues []reportIssue, exitOnWarnings bool) error {
	if len(issues) == 0 {
		return nil
	}
	sort.SliceStable(issues, func(i, j int) bool {
		if issues[i].Line != issues[j].Line {
			return issues[i].Line < issues[j].Line
		}
		if issues[i].Record != issues[j].Record {
			return issues[i].Record < issues[j].Record
		}
		if issues[i].Severity != issues[j].Severity {
			return issues[i].Severity < issues[j].Severity
		}
		return issues[i].Field < issues[j].Field
	})
	var nErr, nWarn int
	for _, iss := range issues {
		switch iss.Severity {
		case reportError:
			nErr++
		case reportWarning:
			nWarn++
		}
		if iss.Record > 0 {
			fmt.Printf("line %d rec %d  [%s] %s: %s\n",
				iss.Line, iss.Record, iss.Severity, iss.Field, iss.Msg)
		} else {
			fmt.Printf("line %d  [%s] %s: %s\n",
				iss.Line, iss.Severity, iss.Field, iss.Msg)
		}
	}
	fmt.Printf("%s: %d error(s), %d warning(s)\n", path, nErr, nWarn)
	if nErr > 0 {
		return fmt.Errorf("%d error(s) found", nErr)
	}
	if exitOnWarnings && nWarn > 0 {
		return fmt.Errorf("%d warning(s) found", nWarn)
	}
	return nil
}

// reportFromVCF converts VCF issues to the shared shape.
func reportFromVCF(in []vcf.ValidationError) []reportIssue {
	out := make([]reportIssue, len(in))
	for i, iss := range in {
		s := reportWarning
		if iss.Severity == vcf.SeverityError {
			s = reportError
		}
		out[i] = reportIssue{
			Line:     iss.Line,
			Record:   iss.Record,
			Field:    iss.Field,
			Severity: s,
			Msg:      iss.Msg,
		}
	}
	return out
}

// reportFromBED converts BED issues to the shared shape.
func reportFromBED(in []bed.Issue) []reportIssue {
	out := make([]reportIssue, len(in))
	for i, iss := range in {
		s := reportWarning
		if iss.Severity == bed.SeverityError {
			s = reportError
		}
		out[i] = reportIssue{
			Line:     iss.Line,
			Field:    iss.Field,
			Severity: s,
			Msg:      iss.Msg,
		}
	}
	return out
}

// reportFromSAM converts SAM issues to the shared shape.
func reportFromSAM(in []sam.Issue) []reportIssue {
	out := make([]reportIssue, len(in))
	for i, iss := range in {
		s := reportWarning
		if iss.Severity == sam.SeverityError {
			s = reportError
		}
		out[i] = reportIssue{
			Line:     iss.Line,
			Field:    iss.Field,
			Severity: s,
			Msg:      iss.Msg,
		}
	}
	return out
}

// reportFromBAM converts BAM issues to the shared shape.
func reportFromBAM(in []bam.Issue) []reportIssue {
	out := make([]reportIssue, len(in))
	for i, iss := range in {
		s := reportWarning
		if iss.Severity == bam.SeverityError {
			s = reportError
		}
		out[i] = reportIssue{
			Line:     iss.Line,
			Field:    iss.Field,
			Severity: s,
			Msg:      iss.Msg,
		}
	}
	return out
}

// reportFromGFF converts GFF issues to the shared shape.
func reportFromGFF(in []gff.Issue) []reportIssue {
	out := make([]reportIssue, len(in))
	for i, iss := range in {
		s := reportWarning
		if iss.Severity == gff.SeverityError {
			s = reportError
		}
		out[i] = reportIssue{
			Line:     iss.Line,
			Field:    iss.Field,
			Severity: s,
			Msg:      iss.Msg,
		}
	}
	return out
}

// reportFromFASTA converts FASTA issues to the shared shape.
func reportFromFASTA(in []fasta.Issue) []reportIssue {
	out := make([]reportIssue, len(in))
	for i, iss := range in {
		s := reportWarning
		if iss.Severity == fasta.SeverityError {
			s = reportError
		}
		out[i] = reportIssue{
			Line:     iss.Line,
			Field:    iss.Field,
			Severity: s,
			Msg:      iss.Msg,
		}
	}
	return out
}

// reportFromFASTQ converts FASTQ issues to the shared shape.
func reportFromFASTQ(in []fastq.Issue) []reportIssue {
	out := make([]reportIssue, len(in))
	for i, iss := range in {
		s := reportWarning
		if iss.Severity == fastq.SeverityError {
			s = reportError
		}
		out[i] = reportIssue{
			Line:     iss.Line,
			Field:    iss.Field,
			Severity: s,
			Msg:      iss.Msg,
		}
	}
	return out
}

// reportFromBCF converts BCF issues to the shared shape.
func reportFromBCF(in []bcf.Issue) []reportIssue {
	out := make([]reportIssue, len(in))
	for i, iss := range in {
		s := reportWarning
		if iss.Severity == bcf.SeverityError {
			s = reportError
		}
		out[i] = reportIssue{
			Line:     iss.Line,
			Field:    iss.Field,
			Severity: s,
			Msg:      iss.Msg,
		}
	}
	return out
}
