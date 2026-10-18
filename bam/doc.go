// Package bam provides a streaming parser and validator for the
// BAM (Binary Alignment/Map) format used to wrap SAM records
// inside a BGZF-compressed binary container.
//
// The implementation targets the BAM v1 specification as
// documented in https://samtools.github.io/hts-specs/SAMv1.pdf
// (the BAM section, § 4). A BAM file is a BGZF-compressed stream
// containing:
//
//   - 4-byte magic "BAM\1"
//   - A SAM header (length-prefixed UTF-8 text, followed by a
//     count and binary reference-sequence dictionary)
//   - A sequence of variable-length alignment records, each
//     prefixed by a 4-byte block size
//
// Typical usage:
//
//	f, _ := os.Open("sample.bam")
//	defer f.Close()
//	dec := bam.NewDecoder(f)
//	hdr, _ := dec.Header()
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	    // rec.QName, rec.Pos, rec.Flag, rec.Bin, rec.Cigar, ...
//	}
//
// The Decoder is not safe for concurrent use; create one per
// goroutine.
//
// BAM is the binary counterpart of the [sam] package: the SAM
// header text is embedded verbatim in the BAM binary header, and
// each alignment record has the same conceptual fields. This
// package does not share types with [sam] (per the per-package
// self-contained design principle).
package bam
