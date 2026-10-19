package gff

import (
	"errors"
	"fmt"
	"io"
)

// fastaWritten tracks whether EncodeFastaSection has already
// been emitted on a given encoder. A second call would produce
// a second `###` terminator followed by another FASTA block,
// which the spec disallows.
var errFastaAlreadyWritten = errors.New("gff: EncodeFastaSection called twice")

// EncodeFastaSection writes the GFF3 "###" terminator followed
// by the verbatim content of r as the optional FASTA section
// of the document. Per the GFF3 spec §2.3, a FASTA section
// after `###` may contain any number of `>`-prefixed sequence
// records; the encoder does not parse or validate this content
// because it is by definition FASTA, not GFF.
//
// Typical usage (composing with a fasta decoder / file):
//
//	enc := gff.NewEncoder(out)
//	_ = enc.EncodeHeader(hdr)
//	for _, r := range recs { _ = enc.EncodeRecord(r) }
//	f, _ := os.Open("sequences.fa")
//	_ = enc.EncodeFastaSection(f)
//	_ = enc.Flush()
//
// EncodeFastaSection is a no-op (other than the `###` line) if
// the FASTA section is empty. It may be called at most once per
// encoder; a second call returns an error.
func (e *Encoder) EncodeFastaSection(r io.Reader) error {
	if e.fastaWritten {
		return errFastaAlreadyWritten
	}
	e.fastaWritten = true
	if _, err := e.w.WriteString("###"); err != nil {
		return err
	}
	if _, err := e.w.Write(e.lineEnding); err != nil {
		return err
	}
	if r == nil {
		return nil
	}
	written, err := io.Copy(e.w, r)
	if err != nil {
		return fmt.Errorf("gff: fasta section: %w", err)
	}
	_ = written
	return nil
}
