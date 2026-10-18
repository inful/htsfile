// Package fastq provides a parser, writer, and validator for
// the FASTQ format used to store short-read sequencing data
// (nucleotide reads plus per-base quality scores).
//
// A FASTQ file is a sequence of records. Each record occupies
// exactly four lines:
//
//	@<header>
//	<sequence>
//	+[<repeated header>]
//	<qualities>
//
// The header line begins with `@`. The plus line begins with
// `+` and may optionally repeat the header text after a
// space. The qualities line carries one ASCII character per
// base; the byte minus the Phred offset (33 for Sanger /
// Illumina 1.8+, 64 for older Illumina 1.3–1.7) gives the
// Phred quality score.
//
// v1 supports Sanger Phred+33 and Illumina 1.8+ (also
// Phred+33, distinguished by header conventions). The Phred
// offset is auto-detected from the first record's quality
// bytes: if any byte is below 33 the data is malformed for
// Phred+33; if the minimum is ≥ 64 the data is treated as
// Phred+64. Callers with non-standard data can override via
// [Decoder.SetPhredOffset]. The Solexa Solexa+64 encoding
// is not supported in v1.
//
// The parser is line-oriented; lines longer than
// [DefaultMaxLineLength] (10 MB) produce a [ParseError]
// rather than an OOM. The encoder emits each record as four
// lines with no sequence wrapping.
//
// Typical usage:
//
//	f, _ := os.Open("reads.fq")
//	defer f.Close()
//	dec := fastq.NewDecoder(f)
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	    // rec.Header, rec.Seq, rec.Qual, rec.PhredOffset
//	}
//
// The Decoder is not safe for concurrent use; create one per
// goroutine.
package fastq
