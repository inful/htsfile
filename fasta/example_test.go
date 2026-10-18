package fasta_test

import (
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/inful/htsfile/fasta"
)

// This example parses a small multi-record FASTA file and prints
// each record's identifier and sequence length.
func ExampleNewDecoder() {
	src := ">seq1 E. coli K-12\n" +
		"ACGTACGT\n" +
		"ACGTACGT\n" +
		">seq2 some other description\n" +
		"GATTACA\n"

	dec := fasta.NewDecoder(strings.NewReader(src))
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		fmt.Printf("%s\t%d bases\n", rec.ID, rec.Length())
	}
	// Output:
	// seq1	16 bases
	// seq2	7 bases
}

// This example writes a small FASTA document with configurable
// sequence-line wrapping.
func ExampleNewEncoder() {
	var buf strings.Builder
	enc := fasta.NewEncoder(&buf)
	enc.SetLineWidth(10) // wrap each sequence line at 10 chars

	if err := enc.EncodeRecord(&fasta.Record{
		ID:          "seq1",
		Description: "example",
		Sequence:    []byte("ACGTACGTACGTACGT"),
	}); err != nil {
		fmt.Println("error:", err)
		return
	}
	if err := enc.Flush(); err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Print(buf.String())
	// Output:
	// >seq1 example
	// ACGTACGTAC
	// GTACGT
}

// This example shows how the strict-mode validator surfaces
// sequence bytes outside the configured alphabet.
func ExampleValidator() {
	src := ">seq1\n" +
		"ACGTACGT\n" +
		">seq2\n" +
		"ACGZACGT\n" // Z is not a valid IUPAC nucleotide.

	dec := fasta.NewDecoder(strings.NewReader(src))
	var recs []*fasta.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		recs = append(recs, rec)
	}

	for _, iss := range fasta.NewValidator(fasta.ModeStrict).Validate(recs) {
		fmt.Printf("line %d  [%s] %s: %s\n",
			iss.Line, iss.Severity, iss.Field, iss.Msg)
	}
	// Output:
	// line 3  [error] sequence: byte 'Z' at position 3 is not in the configured alphabet
}
