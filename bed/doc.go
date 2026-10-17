// Package bed provides a parser, writer, and validator for the BED
// (Browser Extensible Data) format used in genomics.
//
// The implementation targets the format described in the UCSC BED FAQ
// (https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and mirrored in
// the samtools/hts-specs BEDv1 document. It supports BED3 through
// BED12 with column-count auto-detection, plus the optional track and
// browser header lines.
//
// BED lines are tab-separated. The first three columns (chrom,
// chromStart, chromEnd) are required; the remaining nine columns are
// optional and may be present in any prefix order. A file that mixes
// BED3 and BED4 records is malformed; the validator surfaces this as
// an error.
//
// Typical usage:
//
//	f, err := os.Open("peaks.bed")
//	if err != nil { log.Fatal(err) }
//	defer f.Close()
//
//	dec := bed.NewDecoder(f)
//	track, _ := dec.Track()   // optional; (nil, nil) if absent
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { log.Fatal(err) }
//	    fmt.Println(rec.Chrom, rec.ChromStart, rec.ChromEnd)
//	}
//
// The Decoder is not safe for concurrent use; create one per goroutine.
package bed
