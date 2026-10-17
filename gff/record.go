package gff

import (
	"strconv"
	"strings"
)

// Record is one GFF3 feature line. The 9 mandatory columns are
// stored as typed fields; the 9th (attributes) is parsed into a
// multi-value map.
//
// GFF3 coordinates are 1-based inclusive — opposite of BED's
// 0-based half-open convention. The parser rejects start > end
// (the spec requires start <= end).
//
// Score, strand, and phase all support the "." sentinel
// meaning "not applicable". [Record.ScoreKnown] reports whether
// the source line had a real score; [Record.Score] returns
// [ScoreUnknown] (a negative sentinel) otherwise. Strand uses
// [StrUnknown] ('.') for unknown and additionally accepts '?'
// per the spec. Phase uses [PhaseUnknown] (-1) for unknown.
//
// The zero value is not meaningful; obtain Records from
// [Decoder.Record].
type Record struct {
	// Line is the 1-based source line number where the record
	// was read.
	Line int

	// SeqID is column 1: the sequence ID (chromosome,
	// contig, …) the feature is on.
	SeqID string

	// Source is column 2: the source of the feature (e.g. the
	// tool or database that produced it). "." means "not
	// available".
	Source string

	// Type is column 3: the feature type (e.g. "gene", "exon",
	// "mRNA", "CDS", "variation"). "." means "not available".
	Type string

	// Start is column 4: 1-based inclusive start coordinate.
	Start int64

	// End is column 5: 1-based inclusive end coordinate. Per
	// the spec, start <= end.
	End int64

	// Score is column 6: a numeric score, or [ScoreUnknown]
	// if the source had "." for this column. Use
	// [Record.ScoreKnown] to disambiguate.
	Score float64

	// ScoreKnown reports whether [Record.Score] was set from
	// a numeric value (true) or from the "." sentinel (false).
	ScoreKnown bool

	// Strand is column 7: one of '+', '-', '.', '?'.
	Strand byte

	// Phase is column 8: 0, 1, 2, or [PhaseUnknown] (-1).
	// Only meaningful for CDS features.
	Phase int

	// Attributes is the parsed 9th column as a multi-value
	// map. Single-value tags (e.g. ID, Name) have a
	// one-element slice; multi-value tags (Parent, Dbxref,
	// Alias, Ontology_term, Note, Target, Gap, Derives_from)
	// have an n-element slice with one entry per value.
	// Use [Record.Attribute] for the first value and
	// [Record.AttributeList] for the full slice.
	Attributes map[string][]string

	// Raw is the original line text (without the leading
	// newline), preserved for round-trip. If empty, the
	// encoder synthesises a line from the structured fields.
	Raw string
}

// Attribute returns the first value of the named attribute, or
// empty if the attribute is absent. Multi-value attributes
// (e.g. Parent with multiple values) return the first value
// only; use [Record.AttributeList] for the full slice.
func (r *Record) Attribute(key string) string {
	if r == nil {
		return ""
	}
	if vs := r.Attributes[key]; len(vs) > 0 {
		return vs[0]
	}
	return ""
}

// AttributeList returns all values for the named attribute, or
// nil if absent.
func (r *Record) AttributeList(key string) []string {
	if r == nil {
		return nil
	}
	return r.Attributes[key]
}

// HasAttribute reports whether key is present in Attributes.
func (r *Record) HasAttribute(key string) bool {
	if r == nil {
		return false
	}
	_, ok := r.Attributes[key]
	return ok
}

// parseAttributes parses the 9th-column attributes string
// ("key=value;key2=value2;...") into a map. Per the GFF3 spec,
// attributes are separated by `;` and individual values within
// a multi-value attribute (Parent, Dbxref, Alias, …) are
// separated by `,`. Percent-encoded characters are preserved
// verbatim — callers can URL-decode if they need to.
// Use [Record.AttributeList] to retrieve the value slices.
func parseAttributes(s string) map[string][]string {
	if s == "" || s == "." {
		return map[string][]string{}
	}
	out := map[string][]string{}
	for _, kv := range strings.Split(s, ";") {
		kv = strings.TrimSpace(kv)
		if kv == "" {
			continue
		}
		eq := strings.IndexByte(kv, '=')
		if eq < 0 {
			// Per spec, attribute tag without '=' is
			// invalid. Surface as an empty value rather
			// than an error so partial files still parse.
			out[kv] = append(out[kv], "")
			continue
		}
		key := kv[:eq]
		val := kv[eq+1:]
		out[key] = append(out[key], strings.Split(val, ",")...)
	}
	return out
}

// formatAttributes serialises an attribute map back into the
// 9th-column "key=value;key=value" string used by GFF3.
// The encoder keeps stable ordering by sorting the keys.
func formatAttributes(attrs map[string][]string) string {
	if len(attrs) == 0 {
		return "."
	}
	keys := make([]string, 0, len(attrs))
	for k := range attrs {
		keys = append(keys, k)
	}
	// Stable alphabetical order for round-trip determinism.
	for i := 1; i < len(keys); i++ {
		for j := i; j > 0 && keys[j-1] > keys[j]; j-- {
			keys[j-1], keys[j] = keys[j], keys[j-1]
		}
	}
	var b strings.Builder
	for i, k := range keys {
		if i > 0 {
			b.WriteByte(';')
		}
		b.WriteString(k)
		b.WriteByte('=')
		for j, v := range attrs[k] {
			if j > 0 {
				b.WriteByte(',')
			}
			b.WriteString(v)
		}
	}
	return b.String()
}

// formatScore formats a score for the 9th-column slot. If
// ScoreKnown is false, returns ".". Floats are formatted with
// the shortest representation that round-trips back to the
// same value, except NaN/Inf which are not representable in
// GFF3 and surface as ".".
func formatScore(score float64, known bool) string {
	if !known {
		return "."
	}
	if score != score { // NaN
		return "."
	}
	if score > 1e15 || score < -1e15 { // Inf approx
		return "."
	}
	return strconv.FormatFloat(score, 'g', -1, 64)
}