package gff

import "strconv"

// Header is the parsed GFF3 header. It carries the version
// line, well-known pragmas (sequence-region, date), and a
// verbatim copy of every `##` pragma line in source order so
// the encoder can preserve them on round-trip.
type Header struct {
	// Version is the integer version reported by the
	// `##gff-version` line (typically 3). Zero if the
	// directive is absent.
	Version int

	// SequenceRegions is the list of ##sequence-region
	// entries in source order. The ##sequence-region pragma
	// is deprecated by the GFF3 spec but commonly emitted;
	// we preserve it for round-trip.
	SequenceRegions []*SequenceRegion

	// Date is the value of ##date (typically "YYYY-MM-DD"),
	// or empty if absent.
	Date string

	// Comments is the list of bare `#` comment lines in
	// source order. The GFF3 spec reserves the `#` prefix
	// for comments (as opposed to `##` pragmas); some
	// producers emit them.
	Comments []string

	// MetaLines is the verbatim list of every `##` pragma
	// line in source order (including the ones parsed into
	// the typed fields above). Used by the Encoder to
	// preserve unknown pragmas on round-trip.
	MetaLines []string
}

// SequenceRegion is a ##sequence-region pragma entry.
type SequenceRegion struct {
	SeqID string
	Start int64
	End   int64
}

// StrUnknown is the rune used to represent "unknown" strand
// values in parsed records.
const StrUnknown = '.'

// PhaseUnknown is the value used to represent "phase not
// applicable / unknown" — see [Record.PhaseUnknown].
const PhaseUnknown = -1

// ScoreUnknown is the sentinel returned by [Record.Score] when
// the source line used "." for the score column.
const ScoreUnknown = -1

// parsePragma parses one `##` pragma line (the leading "##" and
// trailing newline already stripped). It updates hdr with the
// recognised fields and always appends the line to MetaLines
// for round-trip preservation.
func parsePragma(hdr *Header, line string, lineNo int) error {
	hdr.MetaLines = append(hdr.MetaLines, line)
	// Split into keyword + rest (first whitespace).
	kw, rest := splitFirstSpace(line)
	switch kw {
	case "gff-version":
		n, err := strconv.Atoi(rest)
		if err != nil {
			return &ParseError{Line: lineNo, Field: "pragma",
				Msg: "gff-version: non-numeric value " + rest}
		}
		hdr.Version = n
	case "sequence-region":
		// ##sequence-region seqid start end
		parts := splitFields(rest)
		if len(parts) != 3 {
			return &ParseError{Line: lineNo, Field: "pragma",
				Msg: "sequence-region: want 3 fields, got " + strconv.Itoa(len(parts))}
		}
		start, err := strconv.ParseInt(parts[1], 10, 64)
		if err != nil {
			return &ParseError{Line: lineNo, Field: "pragma",
				Msg: "sequence-region: non-numeric start " + parts[1]}
		}
		end, err := strconv.ParseInt(parts[2], 10, 64)
		if err != nil {
			return &ParseError{Line: lineNo, Field: "pragma",
				Msg: "sequence-region: non-numeric end " + parts[2]}
		}
		hdr.SequenceRegions = append(hdr.SequenceRegions, &SequenceRegion{
			SeqID: parts[0],
			Start: start,
			End:   end,
		})
	case "date":
		hdr.Date = rest
	default:
		// Unknown pragma — already appended to MetaLines.
	}
	return nil
}

// splitFirstSpace returns the substring before the first
// whitespace and the remainder. If no whitespace is present, the
// whole string is the keyword and the rest is empty.
func splitFirstSpace(s string) (kw, rest string) {
	for i := 0; i < len(s); i++ {
		if s[i] == ' ' || s[i] == '\t' {
			return s[:i], s[i+1:]
		}
	}
	return s, ""
}

// splitFields splits on runs of whitespace.
func splitFields(s string) []string {
	var out []string
	i := 0
	for i < len(s) {
		for i < len(s) && (s[i] == ' ' || s[i] == '\t') {
			i++
		}
		if i >= len(s) {
			break
		}
		start := i
		for i < len(s) && s[i] != ' ' && s[i] != '\t' {
			i++
		}
		out = append(out, s[start:i])
	}
	return out
}