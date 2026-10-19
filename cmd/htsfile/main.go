// Command htsfile is a small command-line interface over the
// vcf, bed, sam, bam, gff, fasta, fastq, bcf, tabix, and csi packages.
//
// Usage:
//
//	htsfile check [--strict] [--warnings] file.vcf | file.bed | file.sam | file.bam | file.gff[3] | file.fa[sta] | file.fq | file.bcf | file.[vcf|bcf|gff].tbi | file.[vcf|bcf|gff].csi
//	htsfile roundtrip [--crlf] file.vcf | file.bed | file.sam | file.gff[3] | file.fa[sta] | file.fq
//	htsfile sam2bam < file.sam > file.bam    # convert SAM to BAM
//	htsfile bam2sam < file.bam > file.sam    # convert BAM to SAM
//	htsfile vcf2bcf < file.vcf > file.bcf    # convert VCF to BCF
//	htsfile bcf2vcf < file.bcf > file.vcf    # convert BCF to VCF
//	htsfile vcf query [--index=auto|tabix|csi] file.vcf.gz chr:start-end
//	htsfile vcf index [--format=auto|tabix|csi] [--min-shift=N] [--depth=N] file.vcf.gz
//	htsfile summary file.vcf | file.bed | file.sam | file.bam | file.gff[3] | file.fa[sta] | file.fq | file.bcf | file.[vcf|bcf|gff].tbi | file.[vcf|bcf|gff].csi
//	htsfile version
//
// The format is auto-detected from the file extension (.vcf,
// .bed, .sam, .bam, .gff, .gff3, .fa, .fasta, .fna, .ffn,
// .faa, .fq, .fastq, .bcf, .tbi, .csi). For stdin ("-") or
// unknown extensions, the first non-empty line is sniffed:
// "##fileformat" → VCF;
// "track", "browser", or a tab-separated chrom/start/end
// triple → BED; "@HD" or an 11-column tab-separated record
// with integer FLAG/MAPQ → SAM; the gzip+FEXTRA prefix → BAM
// (also matches BCF's outer wrapper — use the .bcf extension
// to disambiguate); "##gff-version" or a 9-column
// tab-separated feature line → GFF3; a ">" header line → FASTA;
// a "@" header line followed by a sequence line and a "+" line
// → FASTQ. BCF is not auto-sniffed from stdin (it shares its
// outer BGZF container with BAM); use the .bcf extension.
// Tabix (.tbi) and CSI (.csi) are not sniffed — use the
// extension.
//
// Subcommands:
//
//	check     parse and validate the file, exit non-zero on errors.
//	roundtrip parse the file and re-emit it to stdout.
//	summary   print a one-line summary of the file's header.
//	version   print the htsfile version (set via -ldflags at build time).
//
// BAM v1 and BCF v1 are read-only in v1; `roundtrip` is not
// supported for them. Tabix and CSI round-trip is supported
// (read → re-emit) starting with v0.4.0. `check` and `summary`
// work normally for everything.
package main

import (
	"bufio"
	"bytes"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/inful/htsfile/bam"
	"github.com/inful/htsfile/bcf"
	"github.com/inful/htsfile/bed"
	"github.com/inful/htsfile/csi"
	"github.com/inful/htsfile/fasta"
	"github.com/inful/htsfile/fastq"
	"github.com/inful/htsfile/gff"
	"github.com/inful/htsfile/sam"
	"github.com/inful/htsfile/tabix"
	"github.com/inful/htsfile/vcf"
)

// Build-time variables. Populated via:
//
//	go build -ldflags "-X main.version=1.2.3 -X main.commit=$(git rev-parse HEAD) -X main.date=..." .
//
// Default values (when building locally without -ldflags) are intentionally
// informational rather than empty so users see what they're running.
var (
	version = "dev"
	commit  = "none"
	date    = "unknown"
)

// format identifies which parser to use for a file.
type format int

const (
	formatUnknown format = iota
	formatVCF
	formatBED
	formatSAM
	formatBAM
	formatGFF
	formatFASTA
	formatFASTQ
	formatBCF
	formatTabix
	formatCSI
)

func (f format) String() string {
	switch f {
	case formatVCF:
		return "VCF"
	case formatBED:
		return "BED"
	case formatSAM:
		return "SAM"
	case formatBAM:
		return "BAM"
	case formatGFF:
		return "GFF3"
	case formatFASTA:
		return "FASTA"
	case formatFASTQ:
		return "FASTQ"
	case formatBCF:
		return "BCF"
	case formatTabix:
		return "Tabix"
	case formatCSI:
		return "CSI"
	default:
		return "unknown"
	}
}

func main() {
	if len(os.Args) < 2 {
		usage(os.Stderr)
		os.Exit(2)
	}
	cmd := os.Args[1]
	args := os.Args[2:]
	var err error
	switch cmd {
	case "check":
		err = runCheck(args)
	case "roundtrip":
		err = runRoundtrip(args)
	case "sam2bam":
		err = runTranslate(args, runSAMToBAM, "sam2bam", ".sam")
	case "bam2sam":
		err = runTranslate(args, runBAMToSAM, "bam2sam", ".bam")
	case "vcf2bcf":
		err = runTranslate(args, runVCFToBCF, "vcf2bcf", ".vcf")
	case "bcf2vcf":
		err = runTranslate(args, runBCFToVCF, "bcf2vcf", ".bcf")
	case "summary":
		err = runSummary(args)
	case "vcf":
		err = runVCF(args)
	case "version":
		fmt.Printf("htsfile %s (commit %s, built %s)\n", version, commit, date)
		return
	case "-h", "--help", "help":
		usage(os.Stdout)
		return
	default:
		fmt.Fprintf(os.Stderr, "htsfile: unknown command %q\n\n", cmd)
		usage(os.Stderr)
		os.Exit(2)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "htsfile: %v\n", err)
		os.Exit(1)
	}
}

func usage(w io.Writer) {
	_, _ = fmt.Fprintln(w, `Usage: htsfile <command> [flags] <file.vcf|file.bed|file.sam|file.bam|file.gff[3]|file.fa[sta]|file.fq|file.bcf|file.tbi|file.csi|->

Commands:
  check       Validate the file. Exits non-zero if any errors are found.
              --strict     use strict mode (POS monotonic, REF alphabet, ...)
              --warnings   also exit non-zero on warnings
  roundtrip   Re-emit the file to stdout.
              --crlf       use CR+LF line endings
  sam2bam     Convert SAM text to BAM binary (stdout).
  bam2sam     Convert BAM binary to SAM text (stdout).
  vcf2bcf     Convert VCF text to BCF binary (stdout).
  bcf2vcf     Convert BCF binary to VCF text (stdout).
  summary     Print a summary of the file (header / column count / records).
  vcf query    Query an indexed .vcf.gz by region (chrom:start-end).
               --index=auto|tabix|csi   index type (default: auto)
  vcf index    Build a tabix or CSI index next to a BGZF VCF.
               --format=auto|tabix|csi  index format (default: auto → tabix)
               --min-shift=N           CSI binning window exponent (default: 14)
               --depth=N               CSI bin levels (default: 5)

Format is auto-detected from the extension (.vcf, .bed, .sam, .bam, .gff, .gff3, .fa, .fasta, .fq, .bcf, .tbi, .csi). For
stdin ("-") or unknown extensions, the first non-empty line is sniffed.
BAM v1 and BCF v1 are all read-only as text formats, but they
can round-trip through themselves (.bam → .bam, .bcf → .bcf).
The translation subcommands (sam2bam, bam2sam, vcf2bcf, bcf2vcf)
convert between text and binary forms. check and summary work
for everything.`)
}

// openFileOrStdin returns the named file, or os.Stdin if path is "-".
func openFileOrStdin(path string) (io.ReadCloser, error) {
	if path == "-" {
		return io.NopCloser(os.Stdin), nil
	}
	//nolint:gosec // path comes from CLI args; user controls the input
	return os.Open(path)
}

// detectFormat determines the format of path. For known extensions
// (".vcf", ".bed") the extension is trusted. For "-" (stdin) or
// unknown extensions, the file is opened, buffered, and the first
// non-empty line is sniffed. The buffered reader is returned so
// callers can read the data without re-opening (critical for stdin,
// which is single-use). The returned reader is the bufio.Reader
// that did the peek; passing it through to vcf.NewDecoder /
// bed.NewDecoder will re-buffer it, but the peeked data is
// preserved through the double-buffering.
func detectFormat(path string) (format, io.Reader, error) {
	ext := strings.ToLower(filepath.Ext(path))
	switch ext {
	case ".vcf":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatVCF, f, nil
	case ".bed":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatBED, f, nil
	case ".sam":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatSAM, f, nil
	case ".bam":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatBAM, f, nil
	case ".gff", ".gff3":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatGFF, f, nil
	case ".fa", ".fasta", ".fna", ".ffn", ".faa":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatFASTA, f, nil
	case ".fq", ".fastq":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatFASTQ, f, nil
	case ".bcf":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatBCF, f, nil
	case ".tbi":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatTabix, f, nil
	case ".csi":
		f, err := openFileOrStdin(path)
		if err != nil {
			return formatUnknown, nil, err
		}
		return formatCSI, f, nil
	}
	// Unknown extension (or no extension) — open once, buffer, and sniff.
	return sniffFormat(path)
}

// sniffFormat opens the file once, buffers it, and decides the
// format from the first non-empty line. Returns both the detected
// format and the buffered reader; subsequent reads on the returned
// reader see the rest of the data starting from position 0 (the
// Peek does not consume bytes from the caller's perspective once a
// new reader is wrapped around the bufio.Reader).
func sniffFormat(path string) (format, io.Reader, error) {
	f, err := openFileOrStdin(path)
	if err != nil {
		return formatUnknown, nil, err
	}

	br := bufio.NewReader(f)
	// Peek at up to 4 KiB — enough for the first non-empty line in
	// typical BED/VCF headers and for the 4-byte BAM magic.
	peeked, _ := br.Peek(4096)
	// BAM magic check: 0x1f 0x8b 0x08 0x04 (gzip header that every
	// BGZF block starts with — the BAM header itself starts with
	// "BAM\1" but that's *inside* the first BGZF block, so the
	// outer bytes are gzip).
	if len(peeked) >= 4 && bytes.Equal(peeked[:4], bamSniffMagic[:]) {
		return formatBAM, br, nil
	}
	scanner := bufio.NewScanner(bytes.NewReader(peeked))
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, "#") {
			// Comment-style line: VCF comments start with "##",
			// "#CHROM" is the header. Both are VCF.
			return formatVCF, br, nil
		}
		if strings.HasPrefix(line, "##fileformat") {
			return formatVCF, br, nil
		}
		if strings.HasPrefix(line, "track") || strings.HasPrefix(line, "browser") {
			return formatBED, br, nil
		}
		// SAM header line: @HD, @SQ, @RG, @PG, @CO.
		if strings.HasPrefix(line, "@HD") {
			return formatSAM, br, nil
		}
		// GFF3 pragma line: ##gff-version (and other ## pragmas
		// will also be detected below via the 9-column rule).
		if strings.HasPrefix(line, "##gff-version") {
			return formatGFF, br, nil
		}
		// FASTA header line: '>' prefix. The sequence data
		// that follows is multi-line so we only need to look
		// at the first non-blank line.
		if strings.HasPrefix(line, ">") {
			return formatFASTA, br, nil
		}
		// FASTQ sniff is two-stage: the header line starts with
		// '@' but must not be a SAM header tag (@HD/@SQ/@RG/
		// @PG/@CO), and the third non-empty line must start
		// with '+'. We look ahead from the current scanner
		// position; the peek buffer (4 KiB) is enough for the
		// first ~3 lines of any real-world record.
		if strings.HasPrefix(line, "@") && !isSAMHeaderTag(line) {
			if looksLikeFASTQ(scanNonBlankLines(scanner, 4)) {
				return formatFASTQ, br, nil
			}
		}
		// Tab-separated triple: BED.
		fields := strings.Split(line, "\t")
		if len(fields) >= 3 {
			if err := fieldsAreChromStartEnd(fields); err == nil {
				return formatBED, br, nil
			}
		}
		// Tab-separated 11+ field record: SAM alignment. Distinguish
		// from VCF (8 fields) by column count and by the second and
		// fifth fields being integers (FLAG and MAPQ).
		if len(fields) >= 11 && looksLikeSAMRecord(fields) {
			return formatSAM, br, nil
		}
		// 9-column tab-separated line with integer start/end
		// and a single-character strand: GFF3 feature.
		if len(fields) == 9 && looksLikeGFF3Record(fields) {
			return formatGFF, br, nil
		}
		// Default: assume VCF; better error than silent miss.
		return formatVCF, br, nil
	}
	// Empty file: defer to VCF (will fail with a clear error).
	return formatVCF, br, nil
}

// isSAMHeaderTag reports whether line begins with a known
// SAM header tag (@HD, @SQ, @RG, @PG, @CO). The line is
// assumed to be the first non-empty line of a sniffed file.
func isSAMHeaderTag(line string) bool {
	for _, tag := range []string{"@HD", "@SQ", "@RG", "@PG", "@CO"} {
		if strings.HasPrefix(line, tag) {
			return true
		}
	}
	return false
}

// scanNonBlankLines reads up to max lines (skipping blanks)
// from the scanner and returns them as a slice. Returns the
// lines actually read (≤ max). Used by the FASTQ sniff.
func scanNonBlankLines(scanner *bufio.Scanner, maxLines int) []string {
	var out []string
	for scanner.Scan() && len(out) < maxLines {
		s := strings.TrimSpace(scanner.Text())
		if s == "" {
			continue
		}
		out = append(out, s)
	}
	return out
}

// looksLikeFASTQ reports whether lines (the next non-blank
// lines after a '@' header) match the FASTQ shape: a sequence
// line followed by a '+' separator.
func looksLikeFASTQ(lines []string) bool {
	// Need at least 3 lines: @header, sequence, +.
	if len(lines) < 3 {
		return false
	}
	return strings.HasPrefix(lines[2], "+")
}

// looksLikeGFF3Record reports whether fields look like a GFF3
// feature record: exactly 9 tab-separated columns, with the
// 4th and 5th (start, end) being non-negative integers and
// the 7th (strand) being a single character from {+, -, ., ?}.
func looksLikeGFF3Record(fields []string) bool {
	if len(fields) != 9 {
		return false
	}
	if err := parseInt(fields[3]); err != nil {
		return false
	}
	if err := parseInt(fields[4]); err != nil {
		return false
	}
	switch fields[6] {
	case "+", "-", ".", "?":
	default:
		return false
	}
	return true
}

// bamSniffMagic is the first 4 bytes of every BGZF block (and
// therefore of every BAM file): the gzip magic + the deflate
// compression method byte + the FEXTRA flag.
var bamSniffMagic = [4]byte{0x1f, 0x8b, 0x08, 0x04}

// fieldsAreChromStartEnd reports whether the first three fields of
// fields look like chrom / chromStart / chromEnd.
func fieldsAreChromStartEnd(fields []string) error {
	if len(fields) < 3 {
		return errors.New("fewer than 3 fields")
	}
	if fields[1] == "" || fields[2] == "" {
		return errors.New("empty start or end")
	}
	if err := parseInt(fields[1]); err != nil {
		return err
	}
	return parseInt(fields[2])
}

func parseInt(s string) error {
	for _, c := range s {
		if c < '0' || c > '9' {
			return fmt.Errorf("non-numeric: %q", s)
		}
	}
	return nil
}

// looksLikeSAMRecord reports whether fields look like a SAM
// alignment record (11+ tab-separated columns, with column 2
// (FLAG) and column 5 (MAPQ) being non-negative integers).
// Distinguishes SAM from VCF (8 columns) and FASTA (line-based
// headers starting with ">") by structure.
func looksLikeSAMRecord(fields []string) bool {
	if len(fields) < 11 {
		return false
	}
	if err := parseInt(fields[1]); err != nil {
		return false
	}
	if err := parseInt(fields[4]); err != nil {
		return false
	}
	return true
}

// ---------------------------------------------------------------------------
// Shared issue reporter (extracted from the 5-way duplicated
// runCheck{VCF,BED,SAM,BAM,GFF} functions).
//
// Each per-format package has its own Issue type with its own
// Severity constants. To avoid forcing every package to share a
// common Issue shape (which would couple unrelated packages), the
// reporter works on a normalised reportIssue struct and each
// runCheck* function converts its native Issue slice to
// []reportIssue before calling runIssueReport.
// ---------------------------------------------------------------------------

// reportSeverity is the unified severity used by the shared
// reporter. Each per-package SeverityError / SeverityWarning
// constant maps to one of these via the sevFrom* adapters.
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

// reportFromVCF converts a slice of VCF issues to the shared
// reportIssue shape.
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

// runCheck parses and validates the file at args[0].
func runCheck(args []string) error {
	fs := flag.NewFlagSet("check", flag.ExitOnError)
	strict := fs.Bool("strict", false, "use strict validation mode")
	warnings := fs.Bool("warnings", false, "exit non-zero on warnings too")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return errors.New("check: expected exactly one file argument")
	}
	path := fs.Arg(0)
	fmtKind, r, err := detectFormat(path)
	if err != nil {
		return err
	}

	switch fmtKind {
	case formatVCF:
		return runCheckVCF(r, path, *strict, *warnings)
	case formatBED:
		return runCheckBED(r, path, *strict, *warnings)
	case formatSAM:
		return runCheckSAM(r, path, *strict, *warnings)
	case formatBAM:
		return runCheckBAM(r, path, *strict, *warnings)
	case formatGFF:
		return runCheckGFF(r, path, *strict, *warnings)
	case formatFASTA:
		return runCheckFASTA(r, path, *strict, *warnings)
	case formatFASTQ:
		return runCheckFASTQ(r, path, *strict, *warnings)
	case formatBCF:
		return runCheckBCF(r, path, *strict, *warnings)
	case formatTabix:
		return runCheckTabix(r, path)
	case formatCSI:
		return runCheckCSI(r, path)
	default:
		return fmt.Errorf("could not detect format for %s", path)
	}
}

func runCheckVCF(f io.Reader, path string, strict, warnings bool) error {
	mode := vcf.ModeDefault
	if strict {
		mode = vcf.ModeStrict
	}
	dec := vcf.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	var recs []*vcf.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	issues := vcf.NewValidator(mode).Validate(hdr, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	return runIssueReport(path, reportFromVCF(issues), warnings)
}

func runCheckBED(f io.Reader, path string, strict, warnings bool) error {
	mode := bed.ModeDefault
	if strict {
		mode = bed.ModeStrict
	}
	dec := bed.NewDecoder(f)
	track, err := dec.Track()
	if err != nil {
		return fmt.Errorf("parse track: %w", err)
	}
	_ = track
	var recs []*bed.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	issues := bed.NewValidator(mode).Validate(track, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	return runIssueReport(path, reportFromBED(issues), warnings)
}

func runCheckSAM(f io.Reader, path string, strict, warnings bool) error {
	mode := sam.ModeDefault
	if strict {
		mode = sam.ModeStrict
	}
	dec := sam.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	_ = hdr
	var recs []*sam.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	issues := sam.NewValidator(mode).Validate(hdr, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	return runIssueReport(path, reportFromSAM(issues), warnings)
}

// runRoundtrip parses the file at args[0] and re-emits it.
func runRoundtrip(args []string) error {
	fs := flag.NewFlagSet("roundtrip", flag.ExitOnError)
	crlf := fs.Bool("crlf", false, "use CR+LF line endings")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return errors.New("roundtrip: expected exactly one file argument")
	}
	path := fs.Arg(0)
	fmtKind, r, err := detectFormat(path)
	if err != nil {
		return err
	}

	switch fmtKind {
	case formatVCF:
		return runRoundtripVCF(r, *crlf)
	case formatBED:
		return runRoundtripBED(r, *crlf)
	case formatSAM:
		return runRoundtripSAM(r, *crlf)
	case formatBAM:
		return runRoundtripBAM(r)
	case formatGFF:
		return runRoundtripGFF(r, *crlf)
	case formatFASTA:
		return runRoundtripFASTA(r, *crlf)
	case formatFASTQ:
		return runRoundtripFASTQ(r, *crlf)
	case formatBCF:
		return runRoundtripBCF(r)
	case formatTabix:
		return runRoundtripTabix(r)
	case formatCSI:
		return runRoundtripCSI(r)
	default:
		return fmt.Errorf("could not detect format for %s", path)
	}
}

//nolint:dupl // runRoundtripVCF / runRoundtripBED / runRoundtripSAM are
func runRoundtripVCF(f io.Reader, crlf bool) error {
	// intentionally parallel; see runCheckVCF for the rationale.
	dec := vcf.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	var recs []*vcf.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	var enc *vcf.Encoder
	if crlf {
		enc = vcf.NewEncoderCRLF(os.Stdout)
	} else {
		enc = vcf.NewEncoder(os.Stdout)
	}
	if err := enc.EncodeHeader(hdr); err != nil {
		return fmt.Errorf("encode header: %w", err)
	}
	for _, rec := range recs {
		if err := enc.EncodeRecord(rec); err != nil {
			return fmt.Errorf("encode record: %w", err)
		}
	}
	return enc.Flush()
}

func runRoundtripBED(f io.Reader, crlf bool) error {
	dec := bed.NewDecoder(f)
	track, err := dec.Track()
	if err != nil {
		return fmt.Errorf("parse track: %w", err)
	}
	var recs []*bed.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	var enc *bed.Encoder
	if crlf {
		enc = bed.NewEncoderCRLF(os.Stdout)
	} else {
		enc = bed.NewEncoder(os.Stdout)
	}
	if track != nil {
		if err := enc.EncodeTrack(track); err != nil {
			return fmt.Errorf("encode track: %w", err)
		}
	}
	for _, rec := range recs {
		if err := enc.EncodeRecord(rec); err != nil {
			return fmt.Errorf("encode record: %w", err)
		}
	}
	return enc.Flush()
}

//nolint:dupl // See runCheckVCF: same shape, different package.
func runRoundtripSAM(f io.Reader, crlf bool) error {
	dec := sam.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	var recs []*sam.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	var enc *sam.Encoder
	if crlf {
		enc = sam.NewEncoderCRLF(os.Stdout)
	} else {
		enc = sam.NewEncoder(os.Stdout)
	}
	if err := enc.EncodeHeader(hdr); err != nil {
		return fmt.Errorf("encode header: %w", err)
	}
	for _, rec := range recs {
		if err := enc.EncodeRecord(rec); err != nil {
			return fmt.Errorf("encode record: %w", err)
		}
	}
	return enc.Flush()
}

// runTranslate implements the four format-translation
// subcommands (sam2bam, bam2sam, vcf2bcf, bcf2vcf). Each
// takes a single file argument whose extension must match the
// input format, and writes the converted output to stdout.
// The caller passes the per-format conversion function and
// the input extension for the extension check.
func runTranslate(args []string, fn func(io.Reader, io.Writer) error, name, inExt string) error {
	if len(args) != 1 {
		return fmt.Errorf("%s: expected exactly one file argument", name)
	}
	path := args[0]
	ext := strings.ToLower(filepath.Ext(path))
	if ext != inExt {
		return fmt.Errorf("%s: input file %q must have extension %s", name, path, inExt)
	}
	//nolint:gosec // path comes from CLI args; the extension check above is the safety boundary.
	f, err := os.Open(path)
	if err != nil {
		return fmt.Errorf("%s: open: %w", name, err)
	}
	defer func() { _ = f.Close() }()
	return fn(f, os.Stdout)
}

// runSummary prints a one-line summary of the file.
func runSummary(args []string) error {
	fs := flag.NewFlagSet("summary", flag.ExitOnError)
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return errors.New("summary: expected exactly one file argument")
	}
	path := fs.Arg(0)
	fmtKind, r, err := detectFormat(path)
	if err != nil {
		return err
	}

	switch fmtKind {
	case formatVCF:
		return runSummaryVCF(r, path)
	case formatBED:
		return runSummaryBED(r, path)
	case formatSAM:
		return runSummarySAM(r, path)
	case formatBAM:
		return runSummaryBAM(r, path)
	case formatGFF:
		return runSummaryGFF(r, path)
	case formatFASTA:
		return runSummaryFASTA(r, path)
	case formatFASTQ:
		return runSummaryFASTQ(r, path)
	case formatBCF:
		return runSummaryBCF(r, path)
	case formatTabix:
		return runSummaryTabix(r, path)
	case formatCSI:
		return runSummaryCSI(r, path)
	default:
		return fmt.Errorf("could not detect format for %s", path)
	}
}

func runSummaryVCF(f io.Reader, path string) error {
	dec := vcf.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	var nRec int
	for {
		_, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		nRec++
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s\n", path)
	fmt.Fprintf(&b, "  fileformat: %s\n", hdr.FileFormat)
	fmt.Fprintf(&b, "  contigs:    %d\n", len(hdr.Contigs))
	fmt.Fprintf(&b, "  INFO defs:  %d\n", len(hdr.Infos))
	fmt.Fprintf(&b, "  FORMAT defs:%d\n", len(hdr.Formats))
	fmt.Fprintf(&b, "  FILTER defs:%d\n", len(hdr.Filters))
	fmt.Fprintf(&b, "  ALT defs:   %d\n", len(hdr.Alts))
	fmt.Fprintf(&b, "  samples:    %d\n", len(hdr.SampleNames))
	fmt.Fprintf(&b, "  records:    %d\n", nRec)
	_, _ = fmt.Fprint(os.Stdout, b.String())
	return nil
}

func runSummaryBED(f io.Reader, path string) error {
	dec := bed.NewDecoder(f)
	track, err := dec.Track()
	if err != nil {
		return fmt.Errorf("parse track: %w", err)
	}
	columnCounts := make(map[int]int)
	chroms := make(map[string]bool)
	var nRec int
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		columnCounts[rec.ColumnCount()]++
		chroms[rec.Chrom] = true
		nRec++
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s\n", path)
	if track != nil {
		if track.Name != "" {
			fmt.Fprintf(&b, "  track:      %s\n", track.Name)
		}
	}
	bedType := bedTypeFromColumns(columnCounts)
	fmt.Fprintf(&b, "  type:       %s\n", bedType)
	fmt.Fprintf(&b, "  records:    %d\n", nRec)
	fmt.Fprintf(&b, "  chroms:     %d\n", len(chroms))
	_, _ = fmt.Fprint(os.Stdout, b.String())
	return nil
}

//nolint:dupl // See runCheckVCF: same shape as runSummaryBAM and the existing runSummarySAM.
func runSummarySAM(f io.Reader, path string) error {
	dec := sam.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	rnames := make(map[string]int)
	var nRec int
	var nMapped, nUnmapped int
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		rnames[rec.RName]++
		if rec.Flag&sam.FlagUnmap != 0 {
			nUnmapped++
		} else {
			nMapped++
		}
		nRec++
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s\n", path)
	if hdr != nil && hdr.HD != nil {
		if hdr.HD.VN != "" {
			fmt.Fprintf(&b, "  version:    %s\n", hdr.HD.VN)
		}
		if hdr.HD.SO != "" {
			fmt.Fprintf(&b, "  sort order: %s\n", hdr.HD.SO)
		}
		fmt.Fprintf(&b, "  @SQ lines:  %d\n", len(hdr.SQs))
		fmt.Fprintf(&b, "  @RG lines:  %d\n", len(hdr.RGs))
		fmt.Fprintf(&b, "  @PG lines:  %d\n", len(hdr.PGs))
	} else {
		fmt.Fprintf(&b, "  header:     (none)\n")
	}
	fmt.Fprintf(&b, "  records:    %d\n", nRec)
	fmt.Fprintf(&b, "  mapped:     %d\n", nMapped)
	fmt.Fprintf(&b, "  unmapped:   %d\n", nUnmapped)
	_, _ = fmt.Fprint(os.Stdout, b.String())
	return nil
}

func runCheckBAM(f io.Reader, path string, strict, warnings bool) error {
	mode := bam.ModeDefault
	if strict {
		mode = bam.ModeStrict
	}
	dec := bam.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	_ = hdr
	var recs []*bam.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	issues := bam.NewValidator(mode).Validate(hdr, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	return runIssueReport(path, reportFromBAM(issues), warnings)
}

// runRoundtripBAM reads a BAM file and re-emits it to stdout.
// Records stream through verbatim — no sorting — so the
// output is byte-for-byte deterministic given the same input.
func runRoundtripBAM(f io.Reader) error {
	dec := bam.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	enc := bam.NewEncoder(os.Stdout)
	if err := enc.WriteHeader(hdr); err != nil {
		return fmt.Errorf("write header: %w", err)
	}
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		if err := enc.WriteRecord(rec); err != nil {
			return fmt.Errorf("write record: %w", err)
		}
	}
	return enc.Close()
}

func runCheckGFF(f io.Reader, path string, strict, warnings bool) error {
	mode := gff.ModeDefault
	if strict {
		mode = gff.ModeStrict
	}
	dec := gff.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	_ = hdr
	var recs []*gff.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	issues := gff.NewValidator(mode).Validate(hdr, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	return runIssueReport(path, reportFromGFF(issues), warnings)
}

// runRoundtripGFF parses a GFF3 file and re-emits it to stdout.
//
//nolint:dupl // See runCheckVCF: same shape as runRoundtripVCF and the existing runRoundtripSAM.
func runRoundtripGFF(f io.Reader, crlf bool) error {
	dec := gff.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	var recs []*gff.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	var enc *gff.Encoder
	if crlf {
		enc = gff.NewEncoderCRLF(os.Stdout)
	} else {
		enc = gff.NewEncoder(os.Stdout)
	}
	if err := enc.EncodeHeader(hdr); err != nil {
		return fmt.Errorf("encode header: %w", err)
	}
	for _, rec := range recs {
		if err := enc.EncodeRecord(rec); err != nil {
			return fmt.Errorf("encode record: %w", err)
		}
	}
	return enc.Flush()
}

// runSummaryGFF prints a one-line summary of the GFF3 file's
// header.
//
//nolint:dupl // See runCheckVCF: same shape, different package.
func runSummaryGFF(f io.Reader, path string) error {
	dec := gff.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	seqids := make(map[string]bool)
	types := make(map[string]int)
	var nRec int
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		seqids[rec.SeqID] = true
		types[rec.Type]++
		nRec++
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s\n", path)
	if hdr != nil {
		if hdr.Version != 0 {
			fmt.Fprintf(&b, "  version:    %d\n", hdr.Version)
		}
		if hdr.Date != "" {
			fmt.Fprintf(&b, "  date:       %s\n", hdr.Date)
		}
		fmt.Fprintf(&b, "  pragmas:    %d\n", len(hdr.MetaLines))
		fmt.Fprintf(&b, "  sequence-regions: %d\n", len(hdr.SequenceRegions))
	}
	fmt.Fprintf(&b, "  seqids:     %d\n", len(seqids))
	fmt.Fprintf(&b, "  types:      %d unique\n", len(types))
	fmt.Fprintf(&b, "  records:    %d\n", nRec)
	_, _ = fmt.Fprint(os.Stdout, b.String())
	return nil
}

// runSummaryBAM prints a one-line summary of the BAM file's
// header. The format mirrors runSummarySAM.
//
//nolint:dupl // See runCheckVCF: same shape, different package.
func runSummaryBAM(f io.Reader, path string) error {
	dec := bam.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	rnames := make(map[string]int)
	var nRec int
	var nMapped, nUnmapped int
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		rnames[rec.NameRName]++
		if rec.Flag&bam.FlagUnmap != 0 {
			nUnmapped++
		} else {
			nMapped++
		}
		nRec++
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s\n", path)
	if hdr != nil && hdr.HD != nil {
		if hdr.HD.VN != "" {
			fmt.Fprintf(&b, "  version:    %s\n", hdr.HD.VN)
		}
		if hdr.HD.SO != "" {
			fmt.Fprintf(&b, "  sort order: %s\n", hdr.HD.SO)
		}
		fmt.Fprintf(&b, "  @SQ lines:  %d\n", len(hdr.SQs))
		fmt.Fprintf(&b, "  @RG lines:  %d\n", len(hdr.RGs))
		fmt.Fprintf(&b, "  @PG lines:  %d\n", len(hdr.PGs))
	} else {
		fmt.Fprintf(&b, "  header:     (none)\n")
	}
	fmt.Fprintf(&b, "  records:    %d\n", nRec)
	fmt.Fprintf(&b, "  mapped:     %d\n", nMapped)
	fmt.Fprintf(&b, "  unmapped:   %d\n", nUnmapped)
	_, _ = fmt.Fprint(os.Stdout, b.String())
	return nil
}

// bedTypeFromColumns returns a human-readable BED type label
// ("BED3", "BED6-mixed", "BED12", etc.) based on the per-column-
// count record distribution.
func bedTypeFromColumns(counts map[int]int) string {
	if len(counts) == 0 {
		return "empty"
	}
	if len(counts) == 1 {
		for n := range counts {
			return fmt.Sprintf("BED%d", n)
		}
	}
	// Mixed column counts: list the unique values in sorted order.
	keys := make([]int, 0, len(counts))
	for k := range counts {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	parts := make([]string, len(keys))
	for i, k := range keys {
		parts[i] = fmt.Sprintf("BED%d", k)
	}
	return "mixed (" + strings.Join(parts, "+") + ")"
}

// runCheckFASTA parses and validates a FASTA file.
func runCheckFASTA(f io.Reader, path string, strict, warnings bool) error {
	mode := fasta.ModeDefault
	if strict {
		mode = fasta.ModeStrict
	}
	dec := fasta.NewDecoder(f)
	var recs []*fasta.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	issues := fasta.NewValidator(mode).Validate(recs)
	if len(issues) == 0 {
		var totalBases int64
		for _, r := range recs {
			totalBases += r.Length()
		}
		fmt.Printf("%s: clean (%d records, %d bases)\n", path, len(recs), totalBases)
		return nil
	}
	return runIssueReport(path, reportFromFASTA(issues), warnings)
}

// runRoundtripFASTA parses a FASTA file and re-emits it to
// stdout with 80-char sequence wrapping.
func runRoundtripFASTA(f io.Reader, crlf bool) error {
	dec := fasta.NewDecoder(f)
	var recs []*fasta.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	var enc *fasta.Encoder
	if crlf {
		enc = fasta.NewEncoderCRLF(os.Stdout)
	} else {
		enc = fasta.NewEncoder(os.Stdout)
	}
	for _, rec := range recs {
		if err := enc.EncodeRecord(rec); err != nil {
			return fmt.Errorf("encode record: %w", err)
		}
	}
	return enc.Flush()
}

// runSummaryFASTA prints a summary of the FASTA file.
func runSummaryFASTA(f io.Reader, path string) error {
	dec := fasta.NewDecoder(f)
	var nRec int
	var totalBases int64
	ids := make(map[string]int)
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		nRec++
		totalBases += rec.Length()
		if rec.ID != "" {
			ids[rec.ID]++
		}
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s\n", path)
	fmt.Fprintf(&b, "  records:    %d\n", nRec)
	fmt.Fprintf(&b, "  bases:      %d\n", totalBases)
	fmt.Fprintf(&b, "  unique IDs: %d\n", len(ids))
	_, _ = fmt.Fprint(os.Stdout, b.String())
	return nil
}

// runCheckFASTQ parses and validates a FASTQ file. The
// FASTQ format has no header section, so the check is just
// per-record validation against the typed rules.
func runCheckFASTQ(f io.Reader, path string, strict, warnings bool) error {
	mode := fastq.ModeDefault
	if strict {
		mode = fastq.ModeStrict
	}
	dec := fastq.NewDecoder(f)
	var recs []*fastq.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	issues := fastq.NewValidator(mode).Validate(recs)
	if len(issues) == 0 {
		var totalBases int64
		for _, r := range recs {
			totalBases += int64(r.Length())
		}
		fmt.Printf("%s: clean (%d records, %d bases)\n", path, len(recs), totalBases)
		return nil
	}
	return runIssueReport(path, reportFromFASTQ(issues), warnings)
}

// runRoundtripFASTQ parses a FASTQ file and re-emits it to
// stdout. Sequence and quality are written verbatim; Phred
// offset is preserved per-record by re-detecting on the first
// record of the output stream.
func runRoundtripFASTQ(f io.Reader, crlf bool) error {
	dec := fastq.NewDecoder(f)
	var recs []*fastq.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	var enc *fastq.Encoder
	if crlf {
		enc = fastq.NewEncoderCRLF(os.Stdout)
	} else {
		enc = fastq.NewEncoder(os.Stdout)
	}
	for _, rec := range recs {
		if err := enc.EncodeRecord(rec); err != nil {
			return fmt.Errorf("encode record: %w", err)
		}
	}
	return enc.Flush()
}

// runSummaryFASTQ prints a summary of the FASTQ file.
func runSummaryFASTQ(f io.Reader, path string) error {
	dec := fastq.NewDecoder(f)
	var nRec int
	var totalBases int64
	phredOffset := uint8(0)
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		nRec++
		totalBases += int64(rec.Length())
		if phredOffset == 0 {
			phredOffset = rec.PhredOffset
		}
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s\n", path)
	fmt.Fprintf(&b, "  records:    %d\n", nRec)
	fmt.Fprintf(&b, "  bases:      %d\n", totalBases)
	fmt.Fprintf(&b, "  phred+%-3d\tdetected from first record\n", phredOffset)
	_, _ = fmt.Fprint(os.Stdout, b.String())
	return nil
}

// runCheckBCF parses and validates a BCF file. BCF v1 is
// read-only, so the check is structurally identical to the
// other runCheck* family but uses the bcf decoder.
func runCheckBCF(f io.Reader, path string, strict, warnings bool) error {
	mode := bcf.ModeDefault
	if strict {
		mode = bcf.ModeStrict
	}
	dec := bcf.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	_ = hdr
	var recs []*bcf.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			// v1 limitation: FORMAT/sample data is not
			// supported. Surface this with a clear message.
			return fmt.Errorf("parse record: %w", err)
		}
		recs = append(recs, rec)
	}
	issues := bcf.NewValidator(mode).Validate(hdr, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	return runIssueReport(path, reportFromBCF(issues), warnings)
}

// runRoundtripBCF reads a BCF file and re-emits it to stdout.
// Records are written without per-sample FORMAT data (matching
// the v1 BCF reader's scope), so the output is "round-trip
// partial" — header and typed record fields are preserved but
// any per-sample genotype data in the source is dropped.
func runRoundtripBCF(f io.Reader) error {
	dec := bcf.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	enc := bcf.NewEncoder(os.Stdout)
	if err := enc.WriteHeader(hdr); err != nil {
		return fmt.Errorf("write header: %w", err)
	}
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		if err := enc.WriteRecord(rec); err != nil {
			return fmt.Errorf("write record: %w", err)
		}
	}
	return enc.Close()
}

// runSummaryBCF prints a summary of the BCF file's header
// and record count. Per-sample data is not displayed in v1
// (the BCF decoder does not parse the FORMAT column).
func runSummaryBCF(f io.Reader, path string) error {
	dec := bcf.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	var nRec int
	for {
		_, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		nRec++
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s\n", path)
	if hdr != nil {
		if hdr.FileFormat != "" {
			fmt.Fprintf(&b, "  fileformat: %s\n", hdr.FileFormat)
		}
		fmt.Fprintf(&b, "  contigs:    %d\n", len(hdr.Contigs))
		fmt.Fprintf(&b, "  filters:    %d\n", len(hdr.Filters))
		fmt.Fprintf(&b, "  infos:      %d\n", len(hdr.Infos))
		fmt.Fprintf(&b, "  formats:    %d\n", len(hdr.Formats))
	}
	fmt.Fprintf(&b, "  records:    %d\n", nRec)
	_, _ = fmt.Fprint(os.Stdout, b.String())
	return nil
}

// runCheckTabix parses and validates a Tabix index. The
// "validation" is structural (magic, lengths, range checks);
// Tabix has no domain-specific issues like VCF/BED/VCF.
func runCheckTabix(f io.Reader, path string) error {
	idx, err := tabix.NewReader(f).Read()
	if err != nil {
		return fmt.Errorf("parse: %w", err)
	}
	var totalBins, totalChunks int
	for _, ref := range idx.References {
		totalBins += len(ref.Bins)
		for _, chunks := range ref.Bins {
			totalChunks += len(chunks)
		}
	}
	fmt.Printf("%s: clean (%d references, %d bins, %d chunks)\n",
		path, len(idx.References), totalBins, totalChunks)
	return nil
}

// runRoundtripTabix reads a Tabix index and re-emits it to stdout.
// The Writer sorts chunks deterministically by Beg, so the output
// is byte-stable regardless of the input's map iteration order.
func runRoundtripTabix(f io.Reader) error {
	idx, err := tabix.NewReader(f).Read()
	if err != nil {
		return fmt.Errorf("parse: %w", err)
	}
	return tabix.NewWriter(os.Stdout).Write(idx)
}

// runSummaryTabix prints a summary of the Tabix file.
func runSummaryTabix(f io.Reader, path string) error {
	idx, err := tabix.NewReader(f).Read()
	if err != nil {
		return fmt.Errorf("parse: %w", err)
	}
	var totalBins, totalChunks, totalIntv int
	for _, ref := range idx.References {
		totalBins += len(ref.Bins)
		for _, chunks := range ref.Bins {
			totalChunks += len(chunks)
		}
		totalIntv += len(ref.Intervals)
	}
	fmt.Printf("%s\n", path)
	fmt.Printf("  references: %d\n", len(idx.References))
	fmt.Printf("  bins:       %d\n", totalBins)
	fmt.Printf("  chunks:     %d\n", totalChunks)
	fmt.Printf("  intervals:  %d\n", totalIntv)
	return nil
}

// runCheckCSI parses and validates a CSI index.
func runCheckCSI(f io.Reader, path string) error {
	idx, err := csi.NewReader(f).Read()
	if err != nil {
		return fmt.Errorf("parse: %w", err)
	}
	var totalBins, totalChunks int
	for _, ref := range idx.References {
		totalBins += len(ref.Bins)
		for _, chunks := range ref.Bins {
			totalChunks += len(chunks)
		}
	}
	fmt.Printf("%s: clean (min_shift=%d, depth=%d, %d references, %d bins, %d chunks)\n",
		path, idx.MinShift, idx.Depth, len(idx.References), totalBins, totalChunks)
	return nil
}

// runRoundtripCSI reads a CSI index and re-emits it to stdout.
// The Writer sorts chunks deterministically by Beg, so the output
// is byte-stable regardless of the input's map iteration order.
func runRoundtripCSI(f io.Reader) error {
	idx, err := csi.NewReader(f).Read()
	if err != nil {
		return fmt.Errorf("parse: %w", err)
	}
	return csi.NewWriter(os.Stdout).Write(idx)
}

// runSummaryCSI prints a summary of the CSI file.
func runSummaryCSI(f io.Reader, path string) error {
	idx, err := csi.NewReader(f).Read()
	if err != nil {
		return fmt.Errorf("parse: %w", err)
	}
	var totalBins, totalChunks, totalIntv int
	for _, ref := range idx.References {
		totalBins += len(ref.Bins)
		for _, chunks := range ref.Bins {
			totalChunks += len(chunks)
		}
		totalIntv += len(ref.Intervals)
	}
	fmt.Printf("%s\n", path)
	fmt.Printf("  min_shift:  %d  (window = %d bp)\n", idx.MinShift, idx.WindowSize())
	fmt.Printf("  depth:      %d\n", idx.Depth)
	fmt.Printf("  references: %d\n", len(idx.References))
	fmt.Printf("  bins:       %d\n", totalBins)
	fmt.Printf("  chunks:     %d\n", totalChunks)
	fmt.Printf("  intervals:  %d\n", totalIntv)
	return nil
}

// openVCFHeader opens the VCF via the chosen index type
// and returns the parsed header. The returned reader is
// closed before this function returns; callers should
// re-open the index for the main query.
func openVCFHeader(path, indexType string) (*vcf.Header, error) {
	switch indexType {
	case "auto", "tabix":
		r, err := vcf.OpenTabix(path)
		if err != nil {
			return nil, err
		}
		defer func() { _ = r.Close() }()
		return r.Header()
	case "csi":
		r, err := vcf.OpenCSI(path)
		if err != nil {
			return nil, err
		}
		defer func() { _ = r.Close() }()
		return r.Header()
	default:
		return nil, fmt.Errorf("vcf query: unknown index type %q", indexType)
	}
}

// parseVCFRegion parses a samtools-style "chrom:start-end"
// region string into its parts. Both start and end are
// 1-based inclusive (samtools convention); the caller
// converts end to half-open by adding 1.
func parseVCFRegion(s string) (chrom string, start, end int64, err error) {
	idx := strings.LastIndex(s, ":")
	if idx < 0 {
		// Whole-contig query: no range.
		return s, 0, 0, nil
	}
	chrom = s[:idx]
	rangeStr := s[idx+1:]
	parts := strings.SplitN(rangeStr, "-", 2)
	if len(parts) != 2 {
		return "", 0, 0, fmt.Errorf("vcf query: malformed region %q (expected chrom:start-end)", s)
	}
	start, err = strconv.ParseInt(parts[0], 10, 64)
	if err != nil {
		return "", 0, 0, fmt.Errorf("vcf query: bad start in %q: %w", s, err)
	}
	end, err = strconv.ParseInt(parts[1], 10, 64)
	if err != nil {
		return "", 0, 0, fmt.Errorf("vcf query: bad end in %q: %w", s, err)
	}
	if end < start {
		return "", 0, 0, fmt.Errorf("vcf query: end < start in %q", s)
	}
	return chrom, start, end, nil
}

// runVCF dispatches the "vcf" subcommand. Supported
// subcommands are:
//
//   - query — run a region query against an indexed VCF
//   - index — emit a tabix/CSI index next to a VCF
func runVCF(args []string) error {
	if len(args) == 0 {
		return errors.New("vcf: expected subcommand (e.g. \"query\", \"index\")")
	}
	switch args[0] {
	case "query":
		return runVCFQuery(args[1:])
	case "index":
		return runVCFIndex(args[1:])
	default:
		return fmt.Errorf("vcf: unknown subcommand %q", args[0])
	}
}

// clampInt32 narrows an int to int32 with explicit
// clamping at the int32 boundaries so a user-supplied
// `--min-shift=99999` doesn't silently wrap to a tiny
// value at the int32 cast site.
func clampInt32(v int) int32 {
	const maxI32 = int32(^uint32(0) >> 1)
	const minI32 = -maxI32 - 1
	if v > int(maxI32) {
		return maxI32
	}
	if v < int(minI32) {
		return minI32
	}
	return int32(v)
}

// runVCFIndex builds a tabix or CSI index next to a
// BGZF-compressed VCF file. Format selection matches
// `samtools tabix -p vcf` (default: tabix) and
// `bcftools index` (CSI for genomes > 512 MiB).
func runVCFIndex(args []string) error {
	fs := flag.NewFlagSet("vcf index", flag.ExitOnError)
	format := fs.String("format", "auto", "index format: auto, tabix, or csi")
	minShift := fs.Int("min-shift", 14, "CSI min_shift (bin window = 2^N bp; tabix default 14)")
	depth := fs.Int("depth", 5, "CSI depth (number of bin levels; tabix default 5)")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return errors.New("vcf index: expected exactly one file argument")
	}
	path := fs.Arg(0)
	opts := vcf.IndexOpts{
		MinShift: clampInt32(*minShift),
		Depth:    clampInt32(*depth),
	}
	switch *format {
	case "auto":
		// Default to tabix; that is what most VCF
		// tools emit and what the v0.15.1
		// `htsfile vcf query` defaults to.
		opts.Format = vcf.FormatTabix
	case "tabix":
		opts.Format = vcf.FormatTabix
	case "csi":
		opts.Format = vcf.FormatCSI
	default:
		return fmt.Errorf("vcf index: unknown format %q", *format)
	}
	return vcf.BuildIndex(path, opts)
}

// runVCFQuery is the dispatcher for "htsfile vcf query".
// It opens the VCF via either tabix or CSI, runs the
// region query, and writes the matching records to
// stdout in VCF text format.
func runVCFQuery(args []string) error {
	fs := flag.NewFlagSet("vcf query", flag.ExitOnError)
	indexType := fs.String("index", "auto", "index type: auto, tabix, or csi")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 2 {
		return errors.New("vcf query: expected two arguments: <file.vcf.gz> <region>")
	}
	path := fs.Arg(0)
	regionStr := fs.Arg(1)

	chrom, start, end, err := parseVCFRegion(regionStr)
	if err != nil {
		return err
	}
	// If no range was given, fall back to the whole
	// contig. We need to open the reader just to peek at
	// the header and find the contig length.
	if end == 0 {
		hdr, err := openVCFHeader(path, *indexType)
		if err != nil {
			return err
		}
		var found bool
		for _, c := range hdr.Contigs {
			if c.ID == chrom && c.Length > 0 {
				start = 1
				end = c.Length + 1 // half-open
				found = true
				break
			}
		}
		if !found {
			return fmt.Errorf("vcf query: contig %q not declared in header (or length unknown)", chrom)
		}
	} else {
		// The Query API uses half-open [start, end). Convert
		// from samtools-style inclusive end to half-open.
		end++
	}

	// Open the right index. "auto" prefers tabix (the
	// more common format in practice); CSI is used when
	// only .csi exists.
	var recs []*vcf.Record
	switch *indexType {
	case "auto":
		if _, err := os.Stat(path + ".tbi"); err == nil {
			r, err := vcf.OpenTabix(path)
			if err != nil {
				return err
			}
			defer func() { _ = r.Close() }()
			recs, _ = r.Query(chrom, start, end)
		} else if _, err := os.Stat(path + ".csi"); err == nil {
			r, err := vcf.OpenCSI(path)
			if err != nil {
				return err
			}
			defer func() { _ = r.Close() }()
			recs, _ = r.Query(chrom, start, end)
		} else {
			return fmt.Errorf("vcf query: no index found for %s (expected .tbi or .csi)", path)
		}
	case "tabix":
		r, err := vcf.OpenTabix(path)
		if err != nil {
			return err
		}
		defer func() { _ = r.Close() }()
		recs, _ = r.Query(chrom, start, end)
	case "csi":
		r, err := vcf.OpenCSI(path)
		if err != nil {
			return err
		}
		defer func() { _ = r.Close() }()
		recs, _ = r.Query(chrom, start, end)
	default:
		return fmt.Errorf("vcf query: unknown index type %q", *indexType)
	}

	// Emit the matching records in VCF text form.
	enc := vcf.NewEncoder(os.Stdout)
	for _, rec := range recs {
		if err := enc.EncodeRecord(rec); err != nil {
			return err
		}
	}
	return enc.Flush()
}
