package fastq

// PhredEncoding is the formula the qualities on a record are
// meant to be decoded with. The byte range of Solexa+64 and
// Phred+64 is identical (ASCII 64–104); the distinction is
// only meaningful when a downstream tool converts a quality
// byte to an error probability:
//
//   - Phred:    Q = byte - offset
//               error = 10^(-Q/10)
//   - Solexa:   Q = 10 * log10(1 + 10^((byte - offset) / -10))
//               error = 10^(-Q/10) / (1 + 10^(-Q/10))
//
// The library never applies either formula on the caller's
// behalf; PhredEncoding is recorded on each Record so that
// downstream consumers can apply the right one.
type PhredEncoding uint8

const (
	// PhredEncodingPhred is the Sanger / Phred encoding
	// (default for both Phred+33 and Phred+64).
	PhredEncodingPhred PhredEncoding = iota
	// PhredEncodingSolexa is the Solexa Solexa+64 encoding
	// used by Illumina 1.3–1.5 pipelines. The byte range
	// matches Phred+64.
	PhredEncodingSolexa
)

// String returns "phred" or "solexa".
func (e PhredEncoding) String() string {
	if e == PhredEncodingSolexa {
		return "solexa"
	}
	return "phred"
}

// IsValid reports whether e is one of the two documented
// encodings.
func (e PhredEncoding) IsValid() bool {
	return e <= PhredEncodingSolexa
}