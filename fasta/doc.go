// Package fasta provides a parser, writer, and validator for
// the FASTA format used to store nucleotide or protein
// sequences.
//
// A FASTA file is a sequence of records. Each record starts
// with a header line beginning with `>` followed by an
// identifier and a free-text description. The header is
// followed by one or more lines of sequence data. Sequence
// lines are typically wrapped at 60–80 characters but no
// wrapping is mandatory — the parser accepts any wrapping
// (including unwrapped single-line sequences) and the
// encoder emits with a configurable line width.
//
// The parser is line-oriented and tolerant of mixed case and
// mixed alphabet content (DNA, RNA, protein); strict-mode
// validation enforces a configurable alphabet.
//
// Typical usage:
//
//	f, _ := os.Open("sequences.fa")
//	defer f.Close()
//	dec := fasta.NewDecoder(f)
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	    // rec.ID, rec.Description, rec.Sequence, rec.Raw
//	}
//
// Random access into large FASTA files is supported by the
// sibling [fasta/fai] package, which parses the samtools-style
// .fai index format (5 fields per reference line: name,
// length, offset, lineBases, lineWidth).
//
// The Decoder is not safe for concurrent use; create one per
// goroutine.
package fasta