package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"os"

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

// runRoundtrip parses the file at args[0] and re-emits it.
// The dispatch mirrors runCheck's; see check.go for the
// rationale on the per-format duplication pattern.
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

//nolint:dupl // See runCheckVCF: same shape, different package.
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

// runRoundtripBCF reads a BCF file and re-emits it to stdout.
// Records are written without per-sample FORMAT data (matching
// the BCF reader's per-sample scope), so the output is
// "round-trip partial" — header and typed record fields are
// preserved but any per-sample genotype data in the source is
// dropped. v0.16.x restored per-field FORMAT encoding with
// END_OF_VECTOR + spec GT-as-int, so this is no longer a v1
// limitation, but the round-trip helper still preserves the
// historical behaviour.
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
