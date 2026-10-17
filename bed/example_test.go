package bed_test

import (
	"errors"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/inful/htsfile/bed"
)

// This example parses a small BED6 file, prints the records, and
// re-emits them to stdout to demonstrate round-tripping.
func ExampleNewDecoder() {
	src := "track name=\"peaks\" description=\"MACS2 narrow peaks\" useScore=1\n" +
		"chr1\t100\t200\tpeak_1\t500\t+\n" +
		"chr1\t300\t400\tpeak_2\t900\t-\n" +
		"chr2\t500\t600\tpeak_3\t750\t.\n"

	dec := bed.NewDecoder(strings.NewReader(src))
	track, err := dec.Track()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	if track != nil {
		fmt.Printf("track: %s\n", track.Name)
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
		name := "<unnamed>"
		if rec.Name != nil {
			name = *rec.Name
		}
		fmt.Printf("%s\t%d\t%d\t%s\n", rec.Chrom, rec.ChromStart, rec.ChromEnd, name)
	}
	// Output:
	// track: peaks
	// chr1	100	200	peak_1
	// chr1	300	400	peak_2
	// chr2	500	600	peak_3
}

// This example demonstrates round-tripping: parse a BED6 file with
// a track line, validate it, and re-emit it to stdout.
func ExampleEncoder_EncodeRecord() {
	src := "track name=\"peaks\" useScore=1\n" +
		"chr1\t100\t200\tpeak_1\t500\t+\n"

	dec := bed.NewDecoder(strings.NewReader(src))
	track, _ := dec.Track()
	rec, err := dec.Record()
	if err != nil {
		fmt.Println("error:", err)
		return
	}

	enc := bed.NewEncoder(os.Stdout)
	if track != nil {
		_ = enc.EncodeTrack(track)
	}
	_ = enc.EncodeRecord(rec)
	_ = enc.Flush()
	// Output:
	// track name="peaks" useScore=1
	// chr1	100	200	peak_1	500	+
}
