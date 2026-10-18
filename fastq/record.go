package fastq

// Record is one FASTQ record. The four lines of the record
// (header, sequence, plus, qualities) are exposed as Go fields:
//
//   - Header is the text after the leading '@' on the first line.
//   - Seq is the sequence data from the second line. The parser
//     preserves the original byte case and rejects no characters
//     — the IUPAC alphabet check is the validator's job (ModeStrict).
//   - PlusLine is the text after the leading '+' on the third
//     line; many FASTQ files leave this empty but some tools
//     (notably Illumina) repeat the header verbatim.
//   - Qual is the raw ASCII quality bytes from the fourth line.
//     Decode the per-base Phred score as byte - PhredOffset.
//   - PhredOffset records the offset used to decode Qual. It is
//     set by the Decoder's auto-detection on the first record;
//     callers can override via [Decoder.SetPhredOffset].
//
// The zero value is not meaningful; obtain Records from
// [Decoder.Record].
type Record struct {
	// Line is the 1-based line number of the @header line.
	Line int

	// Header is the text after the leading '@', without the
	// trailing newline.
	Header string

	// Seq is the sequence bytes (typically IUPAC nucleotide
	// codes; case is preserved as it appeared on disk).
	Seq []byte

	// PlusLine is the text after the leading '+'. Often
	// empty; some tools (notably Illumina) repeat the
	// header text verbatim.
	PlusLine string

	// Qual is the raw quality ASCII bytes. Decode the
	// per-base Phred score as byte - PhredOffset.
	Qual []byte

	// PhredOffset is the offset used to decode Qual (33 for
	// Sanger / Illumina 1.8+, 64 for older Illumina 1.3–1.7).
	PhredOffset uint8

	// PhredEncoding is the formula the Qual bytes are
	// intended to be decoded with. Defaults to
	// PhredEncodingPhred; callers parsing Solexa+64 data
	// should call [Decoder.SetPhredEncoding] before reading.
	PhredEncoding PhredEncoding
}

// Length returns the number of bases in the record
// (i.e. len(Seq)). Provided for symmetry with [fasta.Record]
// and to avoid having to remember which of len(Seq) and
// len(Qual) is the "right" length — the validator requires
// the two match.
func (r *Record) Length() int {
	if r == nil {
		return 0
	}
	return len(r.Seq)
}
