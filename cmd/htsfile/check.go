package main

import (
	"errors"
	"flag"
	"fmt"
	"io"

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

// runCheck parses and validates the file at args[0]. The
// auto-detected format drives the dispatch into the per-format
// runCheck* family.
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

// Per-format check helpers. Each follows the same shape:
//
//	1. build a Decoder for f
//	2. drain the records
//	3. run the per-format Validator
//	4. print a one-line summary on clean files, or hand the
//	   issues to runIssueReport for the unified reporter
//
// The //nolint:dupl annotation acknowledges the structural
// similarity: each format package exposes the same
// Decoder/Validator contract, and the dispatch is the same
// shape. Forcing a generic abstraction earlier would have
// coupled unrelated packages; the dispatch stays here because
// the per-format Issue types differ.

//nolint:dupl // see runCheck.go header for rationale.
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

//nolint:dupl // see runCheck.go header for rationale.
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

//nolint:dupl // see runCheck.go header for rationale.
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

//nolint:dupl // see runCheck.go header for rationale.
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

//nolint:dupl // see runCheck.go header for rationale.
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
