package sam_test

import (
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/inful/htsfile/sam"
)

// This example parses a small SAM file with a header and a few
// records, prints a summary, and re-emits to stdout.
func ExampleNewDecoder() {
	src := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:1000\n" +
		"r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t16\tchr1\t200\t30\t5S5M\t*\t0\t0\tACGTTACGTA\tJJJJJJJJJJ\n"

	dec := sam.NewDecoder(strings.NewReader(src))
	hdr, _ := dec.Header()
	if hdr != nil && hdr.HD != nil {
		fmt.Printf("version: %s, sort: %s\n", hdr.HD.VN, hdr.HD.SO)
	}
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		fmt.Printf("%s pos=%d mapq=%d tags=%d\n",
			rec.QName, rec.Pos, rec.MapQ, len(rec.Tags))
	}
	// Output:
	// version: 1.6, sort: coordinate
	// r1 pos=100 mapq=60 tags=0
	// r2 pos=200 mapq=30 tags=0
}
