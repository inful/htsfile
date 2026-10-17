// Command htsfile is a small command-line interface over the
// vcf, bed, sam, and bam packages.
//
// Usage:
//
//	htsfile check [--strict] [--warnings] file.vcf | file.bed | file.sam | file.bam
//	htsfile roundtrip [--crlf] file.vcf | file.bed | file.sam
//	htsfile summary file.vcf | file.bed | file.sam | file.bam
//	htsfile version
//
// The format is auto-detected from the file extension (.vcf,
// .bed, .sam, .bam). For stdin ("-") or unknown extensions, the
// first non-empty line is sniffed: "##fileformat" → VCF; "track",
// "browser", or a tab-separated chrom/start/end triple → BED;
// "@HD" or an 11-column tab-separated record with integer
// FLAG/MAPQ → SAM; the 4-byte BAM magic at offset 0 → BAM.
//
// Subcommands:
//
//	check     parse and validate the file, exit non-zero on errors.
//	roundtrip parse the file and re-emit it to stdout.
//	summary   print a one-line summary of the file's header.
//	version   print the htsfile version (set via -ldflags at build time).
//
// BAM v1 is read-only: `roundtrip` is not supported for .bam
// files. `check` and `summary` work normally.
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
	"strings"

	"github.com/inful/htsfile/bam"
	"github.com/inful/htsfile/bed"
	"github.com/inful/htsfile/sam"
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
	case "summary":
		err = runSummary(args)
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
	_, _ = fmt.Fprintln(w, `Usage: htsfile <command> [flags] <file.vcf|file.bed|file.sam|file.bam|->

Commands:
  check       Validate the file. Exits non-zero if any errors are found.
              --strict     use strict mode (POS monotonic, REF alphabet, ...)
              --warnings   also exit non-zero on warnings
  roundtrip   Re-emit the file to stdout.
              --crlf       use CR+LF line endings
  summary     Print a summary of the file (header / column count / records).

Format is auto-detected from the extension (.vcf, .bed, .sam, .bam). For
stdin ("-") or unknown extensions, the first non-empty line is sniffed.
BAM v1 is read-only: the roundtrip subcommand does not support .bam files.`)
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
		// Default: assume VCF; better error than silent miss.
		return formatVCF, br, nil
	}
	// Empty file: defer to VCF (will fail with a clear error).
	return formatVCF, br, nil
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
	default:
		return fmt.Errorf("could not detect format for %s", path)
	}
}

//nolint:dupl // runCheckVCF / runCheckBED / runCheckSAM are intentionally parallel
func runCheckVCF(f io.Reader, path string, strict, warnings bool) error {
	// per-format entry points. They share the same control flow (parse,
	// validate, report) but operate on distinct record types from
	// independent packages. Abstracting them behind a common interface
	// would require a generic helper and yield diminishing returns
	// for the third duplicate.
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
	val := vcf.NewValidator(mode)
	issues := val.Validate(hdr, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	sort.SliceStable(issues, func(i, j int) bool {
		if issues[i].Record != issues[j].Record {
			return issues[i].Record < issues[j].Record
		}
		if issues[i].Line != issues[j].Line {
			return issues[i].Line < issues[j].Line
		}
		if issues[i].Severity != issues[j].Severity {
			return issues[i].Severity < issues[j].Severity
		}
		return issues[i].Field < issues[j].Field
	})
	var nErr, nWarn int
	for _, iss := range issues {
		switch iss.Severity {
		case vcf.SeverityError:
			nErr++
		case vcf.SeverityWarning:
			nWarn++
		}
		fmt.Printf("line %d rec %d  [%s] %s: %s\n",
			iss.Line, iss.Record, iss.Severity, iss.Field, iss.Msg)
	}
	fmt.Printf("%s: %d error(s), %d warning(s)\n", path, nErr, nWarn)
	if nErr > 0 {
		return fmt.Errorf("%d error(s) found", nErr)
	}
	if warnings && nWarn > 0 {
		return fmt.Errorf("%d warning(s) found", nWarn)
	}
	return nil
}

//nolint:dupl // See runCheckVCF: same shape, different package.
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
	val := bed.NewValidator(mode)
	issues := val.Validate(track, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	sort.SliceStable(issues, func(i, j int) bool {
		if issues[i].Line != issues[j].Line {
			return issues[i].Line < issues[j].Line
		}
		if issues[i].Severity != issues[j].Severity {
			return issues[i].Severity < issues[j].Severity
		}
		return issues[i].Field < issues[j].Field
	})
	var nErr, nWarn int
	for _, iss := range issues {
		switch iss.Severity {
		case bed.SeverityError:
			nErr++
		case bed.SeverityWarning:
			nWarn++
		}
		fmt.Printf("line %d  [%s] %s: %s\n",
			iss.Line, iss.Severity, iss.Field, iss.Msg)
	}
	fmt.Printf("%s: %d error(s), %d warning(s)\n", path, nErr, nWarn)
	if nErr > 0 {
		return fmt.Errorf("%d error(s) found", nErr)
	}
	if warnings && nWarn > 0 {
		return fmt.Errorf("%d warning(s) found", nWarn)
	}
	return nil
}

//nolint:dupl // See runCheckVCF: same shape, different package.
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
	val := sam.NewValidator(mode)
	issues := val.Validate(hdr, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	sort.SliceStable(issues, func(i, j int) bool {
		if issues[i].Line != issues[j].Line {
			return issues[i].Line < issues[j].Line
		}
		if issues[i].Severity != issues[j].Severity {
			return issues[i].Severity < issues[j].Severity
		}
		return issues[i].Field < issues[j].Field
	})
	var nErr, nWarn int
	for _, iss := range issues {
		switch iss.Severity {
		case sam.SeverityError:
			nErr++
		case sam.SeverityWarning:
			nWarn++
		}
		fmt.Printf("line %d  [%s] %s: %s\n",
			iss.Line, iss.Severity, iss.Field, iss.Msg)
	}
	fmt.Printf("%s: %d error(s), %d warning(s)\n", path, nErr, nWarn)
	if nErr > 0 {
		return fmt.Errorf("%d error(s) found", nErr)
	}
	if warnings && nWarn > 0 {
		return fmt.Errorf("%d warning(s) found", nWarn)
	}
	return nil
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
		return runRoundtripBAM()
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

// runCheckBAM parses and validates a BAM file. The format
// mirrors runCheckSAM but the records come from the bam
// package's binary decoder.
//
//nolint:dupl // See runCheckVCF: same shape, different package.
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
	val := bam.NewValidator(mode)
	issues := val.Validate(hdr, recs)
	if len(issues) == 0 {
		fmt.Printf("%s: clean (%d records)\n", path, len(recs))
		return nil
	}
	sort.SliceStable(issues, func(i, j int) bool {
		if issues[i].Line != issues[j].Line {
			return issues[i].Line < issues[j].Line
		}
		if issues[i].Severity != issues[j].Severity {
			return issues[i].Severity < issues[j].Severity
		}
		return issues[i].Field < issues[j].Field
	})
	var nErr, nWarn int
	for _, iss := range issues {
		switch iss.Severity {
		case bam.SeverityError:
			nErr++
		case bam.SeverityWarning:
			nWarn++
		}
		fmt.Printf("line %d  [%s] %s: %s\n",
			iss.Line, iss.Severity, iss.Field, iss.Msg)
	}
	fmt.Printf("%s: %d error(s), %d warning(s)\n", path, nErr, nWarn)
	if nErr > 0 {
		return fmt.Errorf("%d error(s) found", nErr)
	}
	if warnings && nWarn > 0 {
		return fmt.Errorf("%d warning(s) found", nWarn)
	}
	return nil
}

// runRoundtripBAM returns an error explaining that BAM v1 is
// read-only and the roundtrip subcommand does not support it.
// The error is intentionally clear so users who try the
// subcommand on a .bam file see the reason.
func runRoundtripBAM() error {
	return errors.New("roundtrip: BAM v1 is read-only in this build (no BAM writer)")
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
