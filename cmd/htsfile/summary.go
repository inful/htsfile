package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"sort"
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

// runSummary prints a one-line-per-block overview of the file
// at args[0]. Dispatch mirrors runCheck.
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

// runSummaryBCF prints a summary of the BCF file's header
// and record count. Per-sample data is not displayed (the BCF
// decoder reads it but the summary focuses on header counts).
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
