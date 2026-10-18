// Package sam provides a parser, writer, and validator for the SAM
// (Sequence Alignment/Map) format used in genomics.
//
// The implementation targets the SAM v1 specification as documented
// in https://samtools.github.io/hts-specs/SAMv1.pdf. Each SAM file
// consists of an optional header (lines beginning with "@") followed
// by tab-separated alignment records. Each record has 11 mandatory
// columns (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT,
// TLEN, SEQ, QUAL) plus a variable number of optional fields
// (TAG:TYPE:VALUE).
//
// Header lines:
//
//	@HD   file-level metadata (VN, SO, GO)
//	@SQ   reference sequence dictionary (SN, LN, ...)
//	@RG   read group (ID, SM, LB, ...)
//	@PG   program record (ID, PN, CL, ...)
//	@CO   free-text comment
//
// Tag type codes: A (printable char), C (0-255 int), f (float),
// H (hex), i (int), Z (string), B (array of cCsSiIf).
//
// CIGAR semantic parsing:
//
//	ops, err := sam.ParseCigar(rec.Cigar)
//	// ops is a []sam.CigarUnit{ Len, Op } for each (length, op) pair.
//	// Op is one of CigarM, CigarI, CigarD, CigarN, CigarS, CigarH,
//	// CigarP, CigarEq, CigarX. Per-op helpers:
//
//	sam.ReferenceLength(ops) // total reference span
//	sam.QueryLength(ops)     // total query (SEQ) span
//	sam.AlignmentEnd(rec.Pos, ops) // 1-based rightmost ref position
//
// FormatCigar renders a parsed CIGAR back to its textual form.
//
// Typical usage:
//
//	f, err := os.Open("sample.sam")
//	if err != nil { log.Fatal(err) }
//	defer f.Close()
//
//	dec := sam.NewDecoder(f)
//	hdr, _ := dec.Header()    // optional; nil if no header
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { log.Fatal(err) }
//	    // rec.QName, rec.Flag, rec.Tags, ...
//	}
//
// The Decoder is not safe for concurrent use; create one per
// goroutine.
package sam
