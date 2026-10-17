// Package gff provides a parser, writer, and validator for the
// GFF3 (Generic Feature Format Version 3) format used to
// describe biological features (genes, transcripts, exons,
// variants, …) on nucleotide sequences.
//
// The implementation targets the GFF3 specification as
// documented at
// https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md.
//
// A GFF3 file consists of an optional set of `##` pragma
// directives (at least `##gff-version 3`), zero or more feature
// records (9 tab-separated columns), and `###` record
// terminators. Feature records use 1-based inclusive
// coordinates (the opposite convention from BED, which uses
// 0-based half-open ranges).
//
// Typical usage:
//
//	f, _ := os.Open("annotations.gff")
//	defer f.Close()
//	dec := gff.NewDecoder(f)
//	hdr, _ := dec.Header()    // required; returns error on bad pragma
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	    // rec.SeqID, rec.Source, rec.Type, rec.Start, rec.End,
//	    // rec.Score, rec.Strand, rec.Phase, rec.Attributes
//	    // rec.Attribute("Parent") returns the first Parent value
//	    // rec.AttributeList("Parent") returns all Parent values
//	}
//
// The Decoder is not safe for concurrent use; create one per
// goroutine.
//
// GFF3 is structurally distinct from the other supported
// formats (VCF, BED, SAM, BAM) and lives in its own package per
// the per-package self-contained design principle.
package gff