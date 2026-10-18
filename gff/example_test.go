package gff_test

import (
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/inful/htsfile/gff"
)

// This example parses a small GFF3 file with two features and
// prints the identifier, type, and Parent attribute of each.
func ExampleNewDecoder() {
	src := "##gff-version 3\n" +
		"##sequence-region ctgA 1 1497228\n" +
		"ctgA\tEMBL\tgene\t1000\t9000\t.\t+\t.\tID=gene1;Name=ABC1\n" +
		"ctgA\tEMBL\texon\t1000\t2000\t.\t+\t.\tParent=gene1\n"

	dec := gff.NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Printf("gff-version: %d, pragmas: %d\n", hdr.Version, len(hdr.MetaLines))

	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		fmt.Printf("%s %s parent=%s\n",
			rec.Type, rec.Attribute("ID"), rec.Attribute("Parent"))
	}
	// Output:
	// gff-version: 3, pragmas: 2
	// gene gene1 parent=
	// exon  parent=gene1
}

// This example writes a small GFF3 document with default LF
// line endings and round-trips a record's typed columns.
func ExampleNewEncoder() {
	var buf strings.Builder
	enc := gff.NewEncoder(&buf)

	hdr := &gff.Header{
		Version: 3,
		// MetaLines are written verbatim with a "##" prefix
		// added by the encoder; do not include the prefix here.
		MetaLines: []string{
			"gff-version 3",
			"sequence-region ctgA 1 1497228",
		},
	}
	if err := enc.EncodeHeader(hdr); err != nil {
		fmt.Println("error:", err)
		return
	}
	if err := enc.EncodeRecord(&gff.Record{
		SeqID:      "ctgA",
		Source:     "EMBL",
		Type:       "gene",
		Start:      1000,
		End:        9000,
		Strand:     '+',
		Phase:      gff.PhaseUnknown, // emits "."
		Attributes: map[string][]string{"ID": {"gene1"}, "Name": {"ABC1"}},
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
	// ##gff-version 3
	// ##sequence-region ctgA 1 1497228
	// ctgA	EMBL	gene	1000	9000	.	+	.	ID=gene1;Name=ABC1
}

// This example shows how the strict-mode validator surfaces
// semantic problems — here, a `Parent` value that doesn't
// match any `ID` in the document.
func ExampleValidator() {
	src := "##gff-version 3\n" +
		"ctgA\tEMBL\texon\t1000\t2000\t.\t+\t.\tParent=missing\n"

	dec := gff.NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	var recs []*gff.Record
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

	for _, iss := range gff.NewValidator(gff.ModeStrict).Validate(hdr, recs) {
		fmt.Printf("line %d  [%s] %s: %s\n",
			iss.Line, iss.Severity, iss.Field, iss.Msg)
	}
	// Output:
	// line 2  [error] Parent: Parent="missing" not declared as any record's ID
}
