// Package bcf provides a reader for the BCF v2.1 (Binary Call
// Format) format — the binary counterpart of VCF v4.x.
//
// A BCF file is BGZF-compressed and consists of:
//
//   - A 6-byte magic: "BCF\x02\x01"
//   - A length-prefixed VCF-style header (UTF-8 text)
//   - BGZF blocks of typed records
//
// Each BCF record carries the shared per-record fields (CHROM,
// POS, ID, REF, ALT, QUAL, FILTER, INFO) as typed binary data
// keyed by integer indices into the header's @SQ and @FILTER
// tables. The FORMAT column (per-sample genotype data) is
// intentionally out of scope for v1 — see [Decoder] for the
// current field set.
//
// Typical usage:
//
//	f, _ := os.Open("sample.bcf")
//	defer f.Close()
//	dec := bcf.NewDecoder(f)
//	hdr, _ := dec.Header()
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	    // rec.CHROM, rec.Pos, rec.ID, rec.Ref, rec.Alt,
//	    // rec.Qual, rec.Filter, rec.Info
//	}
//
// v1 is read-only. Writer support is a follow-up.
package bcf
