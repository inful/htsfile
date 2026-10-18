package fasta

import "strings"

// Record is one FASTA record. The header (without the leading
// `>` and trailing newline) is split into ID (first
// whitespace-delimited token) and Description (everything
// after the first whitespace, with leading whitespace
// stripped). Sequence is stored as a byte slice with any
// embedded whitespace stripped — the original wrapping (if
// any) is not preserved unless callers configure the encoder
// to emit a specific line width on round-trip.
//
// The zero value is not meaningful; obtain Records from
// [Decoder.Record].
type Record struct {
	// Line is the 1-based line number of the `>` header line.
	Line int

	// ID is the first whitespace-delimited token of the
	// header (e.g. ">NC_000001.1 Escherichia coli ..." has
	// ID="NC_000001.1"). Empty if the header was empty
	// between `>` and the first whitespace.
	ID string

	// Description is the header text after the first
	// whitespace-delimited token (leading whitespace
	// stripped). Empty if the header had no whitespace.
	Description string

	// Sequence is the concatenated sequence data with
	// embedded whitespace removed (preserves the original
	// character case as it appeared on disk). The empty
	// string indicates the record has no sequence data.
	Sequence []byte

	// Raw is the original header line text, without the
	// leading `>` and without the trailing newline. The
	// encoder uses Raw verbatim when set, falling back to a
	// synthesised "ID Description" line otherwise.
	Raw string
}

// Length returns the number of sequence characters (i.e.
// len(Sequence)). Provided for symmetry with the FASTA spec
// and the [fasta/fai] index's "length" field.
func (r *Record) Length() int64 {
	if r == nil {
		return 0
	}
	return int64(len(r.Sequence))
}

// parseHeader splits the header text (already stripped of the
// leading `>` and trailing newline) into ID and Description.
// The first whitespace-delimited token is ID; everything after
// the first run of whitespace is Description with leading
// whitespace trimmed.
func parseHeader(raw string) (id, desc string) {
	raw = strings.TrimRight(raw, "\r")
	i := 0
	for i < len(raw) && raw[i] != ' ' && raw[i] != '\t' {
		i++
	}
	id = raw[:i]
	desc = strings.TrimLeft(raw[i:], " \t")
	return id, desc
}

// stripWhitespace returns a copy of b with all ASCII whitespace
// bytes removed.
func stripWhitespace(b []byte) []byte {
	if !containsWhitespace(b) {
		return b
	}
	out := make([]byte, 0, len(b))
	for _, c := range b {
		if c == ' ' || c == '\t' || c == '\r' || c == '\n' {
			continue
		}
		out = append(out, c)
	}
	return out
}

func containsWhitespace(b []byte) bool {
	for _, c := range b {
		if c == ' ' || c == '\t' || c == '\r' || c == '\n' {
			return true
		}
	}
	return false
}
