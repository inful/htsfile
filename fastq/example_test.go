package fastq_test

import (
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/inful/htsfile/fastq"
)

// This example parses a small multi-record FASTQ file and
// prints each record's identifier, sequence, and detected
// Phred offset. The decoder auto-detects the offset from
// the FIRST record's quality bytes; subsequent records in a
// real-world stream use the same offset.
func ExampleNewDecoder() {
	src := "@read1\n" +
		"ACGT\n" +
		"+\n" +
		"hhhh\n" + // Phred+64 (hhhh = Q40 in Illumina 1.3–1.7) — first record drives the offset.
		"@read2\n" +
		"GGGG\n" +
		"+\n" +
		"!!!!\n" // Looks like Phred+33, but read1's offset sticks.

	dec := fastq.NewDecoder(strings.NewReader(src))
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		fmt.Printf("%s\t%s\tphred+%-d\n",
			rec.Header, rec.Seq, rec.PhredOffset)
	}
	// Output:
	// read1	ACGT	phred+64
	// read2	GGGG	phred+64
}

// This example writes a small FASTQ document with default
// LF line endings.
func ExampleNewEncoder() {
	var buf strings.Builder
	enc := fastq.NewEncoder(&buf)
	for _, r := range []*fastq.Record{
		{Header: "r1", Seq: []byte("ACGT"), Qual: []byte("!!!!"), PhredOffset: 33},
		{Header: "r2", Seq: []byte("GGGG"), Qual: []byte("JJJJ"), PhredOffset: 33},
	} {
		if err := enc.EncodeRecord(r); err != nil {
			fmt.Println("error:", err)
			return
		}
	}
	if err := enc.Flush(); err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Print(buf.String())
	// Output:
	// @r1
	// ACGT
	// +
	// !!!!
	// @r2
	// GGGG
	// +
	// JJJJ
}

// This example shows how to force a specific Phred offset
// (for example, when parsing data from a known-older
// Illumina run that the auto-detect heuristic should not
// override).
func ExampleDecoder_SetPhredOffset() {
	src := "@r\nACGT\n+\n!!!!\n" // quality bytes look like Phred+33

	dec := fastq.NewDecoder(strings.NewReader(src))
	dec.SetPhredOffset(64) // force Phred+64 anyway

	rec, err := dec.Record()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Printf("forced offset = %d\n", rec.PhredOffset)
	// Output:
	// forced offset = 64
}

// This example shows how the strict-mode validator surfaces
// non-IUPAC sequence bytes.
func ExampleValidator() {
	src := "@r1\nACGTACGT\n+\n!!!!!!!!\n" +
		"@r2\nACGZACGT\n+\n!!!!!!!!\n" // Z is not IUPAC.

	dec := fastq.NewDecoder(strings.NewReader(src))
	var recs []*fastq.Record
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

	for _, iss := range fastq.NewValidator(fastq.ModeStrict).Validate(recs) {
		fmt.Printf("line %d  [%s] %s: %s\n",
			iss.Line, iss.Severity, iss.Field, iss.Msg)
	}
	// Output:
	// line 5  [error] sequence: byte 'Z' at position 3 is not in the configured alphabet
}
