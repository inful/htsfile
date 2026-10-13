// Command htsfile is a small command-line interface over the vcf package.
//
// Usage:
//
//	htsfile check [--strict] [--crlf] file.vcf
//	htsfile roundtrip [--crlf] file.vcf
//	htsfile summary file.vcf
//
// Subcommands:
//
//	check     parse and validate the file, exit non-zero on errors.
//	roundtrip parse the file and re-emit it to stdout.
//	summary   print a one-line summary of the file header.
package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"sort"
	"strings"

	"github.com/inful/htsfile/vcf"
)

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
	_, _ = fmt.Fprintln(w, `Usage: htsfile <command> [flags] file.vcf

Commands:
  check       Validate the file. Exits non-zero if any errors are found.
              --strict     use strict mode (POS monotonic, REF alphabet, ...)
              --warnings   also exit non-zero on warnings
              --crlf       write output with CR+LF line endings
  roundtrip   Re-emit the file to stdout.
              --crlf       use CR+LF line endings
  summary     Print a one-line summary of the file's header.`)
}

// openFileOrStdin returns the named file, or os.Stdin if path is "-".
func openFileOrStdin(path string) (io.ReadCloser, error) {
	if path == "-" {
		return io.NopCloser(os.Stdin), nil
	}
	//nolint:gosec // path comes from CLI args; user controls the input
	return os.Open(path)
}

// runCheck parses and validates the file at args[0].
func runCheck(args []string) error {
	fs := flag.NewFlagSet("check", flag.ExitOnError)
	strict := fs.Bool("strict", false, "use strict validation mode")
	warnings := fs.Bool("warnings", false, "exit non-zero on warnings too")
	crlf := fs.Bool("crlf", false, "use CR+LF line endings on output (no effect for check)")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return errors.New("check: expected exactly one file argument")
	}
	_ = crlf // reserved for future use

	f, err := openFileOrStdin(fs.Arg(0))
	if err != nil {
		return err
	}
	defer f.Close() //nolint:errcheck // best-effort close

	mode := vcf.ModeDefault
	if *strict {
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
		fmt.Printf("%s: clean (%d records)\n", fs.Arg(0), len(recs))
		return nil
	}

	// Sort by (record, line, severity, field) for stable output.
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
	fmt.Printf("%s: %d error(s), %d warning(s)\n", fs.Arg(0), nErr, nWarn)
	if nErr > 0 {
		return fmt.Errorf("%d error(s) found", nErr)
	}
	if *warnings && nWarn > 0 {
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

	f, err := openFileOrStdin(fs.Arg(0))
	if err != nil {
		return err
	}
	defer f.Close() //nolint:errcheck // best-effort close

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
	if *crlf {
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

// runSummary prints a one-line summary of the file header.
func runSummary(args []string) error {
	fs := flag.NewFlagSet("summary", flag.ExitOnError)
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return errors.New("summary: expected exactly one file argument")
	}

	f, err := openFileOrStdin(fs.Arg(0))
	if err != nil {
		return err
	}
	defer f.Close() //nolint:errcheck // best-effort close

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
	fmt.Fprintf(&b, "%s\n", fs.Arg(0))
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
