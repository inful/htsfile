package bcf

import (
	"sort"
	"strconv"
	"strings"
)

// Header is the parsed BCF header. The BCF header is VCF text
// stored inside a binary container; the parser extracts enough
// structure to resolve the integer indices that appear in
// BCF records (CHROM → contig name, FILTER → filter ID) and
// to type-coerce INFO field values.
//
// Header is per-package self-contained: it intentionally does
// not import the vcf package's Header type, even though the
// two parse the same text format. The ≥3-package threshold for
// extracting shared header-text parsing hasn't been crossed
// yet.
type Header struct {
	// Raw is the verbatim VCF header text (newlines included,
	// no trailing null). The encoder/CLI surface uses this
	// to round-trip and to print a summary.
	Raw string

	// FileFormat is the value of ##fileformat (typically
	// "VCFv4.2"). Empty if absent.
	FileFormat string

	// Contigs maps contig index → name. BCF records store
	// CHROM as an index into this list; the decoder resolves
	// it to the contig name.
	Contigs []string

	// Filters maps filter index → ID. BCF records store
	// FILTER as an index (or special value 0xFFFFFFFF for
	// "missing", 0 for "PASS").
	Filters []string

	// Infos maps INFO ID → definition. The Number / Type
	// fields drive BCF's typed-array value parsing.
	Infos map[string]*InfoDef

	// InfosOrder is the source-order list of INFO IDs.
	// BCF INFO fields are keyed by a position (uint16) into
	// this list, not by name; this slice preserves the
	// declaration order so the decoder can resolve the
	// positional keys back to declarations.
	InfosOrder []string

	// Formats maps FORMAT ID → definition. Reserved for the
	// per-sample data path (out of scope for v1).
	Formats map[string]*FormatDef

	// FormatsOrder is the source-order list of FORMAT IDs,
	// parallel to InfosOrder.
	FormatsOrder []string

	// SampleNames is the list of sample IDs from the #CHROM
	// header line. Empty if no FORMAT column.
	SampleNames []string

	// dictCache memoises the dictionary order derived from Raw.
	dictCache []string
}

// InfoDef describes an INFO field declaration (`##INFO=<ID,…>`).
type InfoDef struct {
	ID          string
	Number      string // "A", "R", "G", ".", or an integer (per VCF §1.4.2)
	Type        string // "Integer", "Float", "Flag", "Character", "String"
	Description string
}

// FormatDef describes a FORMAT field declaration (`##FORMAT=<ID,…>`).
type FormatDef struct {
	ID          string
	Number      string
	Type        string
	Description string
}

// LookupContig returns the contig name at index i, or empty
// if i is out of range.
func (h *Header) LookupContig(i int) string {
	if h == nil || i < 0 || i >= len(h.Contigs) {
		return ""
	}
	return h.Contigs[i]
}

// LookupFilter returns the filter ID at index i, or empty if
// i is out of range. The special value 0 traditionally maps
// to "PASS" in VCF/BCF; the decoder handles that mapping.
func (h *Header) LookupFilter(i int) string {
	// Per BCF §GenotypeEncoding (VCFv4.5 line 2256) and
	// htslib convention: filter index 0 is the implicit
	// PASS, even when explicit ##FILTER lines are declared.
	// Declared filters therefore occupy indices 1..N.
	if i == 0 {
		return "PASS"
	}
	if h == nil || i < 0 || i-1 >= len(h.Filters) {
		return ""
	}
	return h.Filters[i-1]
}

// ParseHeader parses a VCF-format header text into a Header.
// It is intentionally minimal: only the structured ## lines
// (fileformat, contig, FILTER, INFO, FORMAT) and the #CHROM
// sample-name line are parsed. Unrecognised ## lines are
// silently ignored (their content is preserved in Raw but not
// surfaced as typed fields).
//
// Exposed (capitalised) so test helpers can reuse the same
// parsing logic as the decoder.
func ParseHeader(text string) *Header {
	h := &Header{
		Raw:     text,
		Infos:   make(map[string]*InfoDef),
		Formats: make(map[string]*FormatDef),
	}
	for _, line := range strings.Split(text, "\n") {
		if line == "" {
			continue
		}
		switch {
		case strings.HasPrefix(line, "##fileformat="):
			h.FileFormat = strings.TrimPrefix(line, "##fileformat=")
		case strings.HasPrefix(line, "##contig="):
			// ##contig=<ID,…>
			id := extractIDFromAngleField(line, "##contig=")
			if id != "" {
				h.Contigs = append(h.Contigs, id)
			}
		case strings.HasPrefix(line, "##FILTER="):
			// ##FILTER=<ID,…>
			id := extractIDFromAngleField(line, "##FILTER=")
			if id != "" {
				h.Filters = append(h.Filters, id)
			}
		case strings.HasPrefix(line, "##INFO="):
			def := parseStructuredField(line, "##INFO=")
			if def != nil && def.ID != "" {
				h.Infos[def.ID] = def
				h.InfosOrder = append(h.InfosOrder, def.ID)
			}
		case strings.HasPrefix(line, "##FORMAT="):
			def := parseStructuredField(line, "##FORMAT=")
			if def != nil && def.ID != "" {
				h.Formats[def.ID] = &FormatDef{
					ID:          def.ID,
					Number:      def.Number,
					Type:        def.Type,
					Description: def.Description,
				}
				h.FormatsOrder = append(h.FormatsOrder, def.ID)
			}
		case strings.HasPrefix(line, "#CHROM"):
			// #CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT sample1…]
			fields := strings.Split(line, "\t")
			if len(fields) > 8 {
				h.SampleNames = fields[9:]
			}
		}
	}
	return h
}

// extractIDFromAngleField extracts the ID= component from a
// header line of the form `##<TYPE>=<ID=value,…>`.
func extractIDFromAngleField(line, prefix string) string {
	rest := strings.TrimPrefix(line, prefix)
	// rest looks like "<ID=chr1,length=…,…>" or "<ID=PASS,Description=…>".
	// Strip the angle brackets and read the ID= pair.
	if strings.HasPrefix(rest, "<") {
		rest = strings.TrimPrefix(rest, "<")
		if end := strings.Index(rest, ">"); end >= 0 {
			rest = rest[:end]
		}
	}
	for _, kv := range splitStructuredFields(rest, ",") {
		if strings.HasPrefix(kv, "ID=") {
			return strings.TrimPrefix(kv, "ID=")
		}
	}
	return ""
}

// parseStructuredField parses a ##INFO= or ##FORMAT= line
// into an InfoDef / FormatDef.
func parseStructuredField(line, prefix string) *InfoDef {
	rest := strings.TrimPrefix(line, prefix)
	if strings.HasPrefix(rest, "<") {
		rest = strings.TrimPrefix(rest, "<")
		if end := strings.Index(rest, ">"); end >= 0 {
			rest = rest[:end]
		}
	}
	def := &InfoDef{}
	// Comma-aware split: the Description value can contain
	// commas inside double quotes. Treat '"' as toggling a
	// "inside string" state and skip commas in that state.
	for _, kv := range splitStructuredFields(rest, ",") {
		kv = strings.TrimSpace(kv)
		switch {
		case strings.HasPrefix(kv, "ID="):
			def.ID = strings.TrimPrefix(kv, "ID=")
		case strings.HasPrefix(kv, "Number="):
			def.Number = strings.TrimPrefix(kv, "Number=")
		case strings.HasPrefix(kv, "Type="):
			def.Type = strings.TrimPrefix(kv, "Type=")
		case strings.HasPrefix(kv, "Description="):
			// Strip outer quotes if present so callers
			// see the bare description.
			d := strings.TrimPrefix(kv, "Description=")
			d = strings.TrimPrefix(d, "\"")
			if d != "" && d[len(d)-1] == '"' {
				d = d[:len(d)-1]
			}
			def.Description = d
		}
	}
	return def
}

// ParseNumber interprets an INFO/FORMAT Number field. The
// VCF spec allows Number to be "A" (per-alt-allele), "R"
// (per-ref+alt allele), "G" (per-genotype), "." (one per
// type), or a non-negative integer. Returns (n, isVariable)
// where isVariable is true for A/R/G/. modes.
func ParseNumber(s string) (n int, isVariable bool) {
	switch s {
	case "A", "R", "G", ".":
		return 0, true
	case "":
		return 0, false
	}
	v, err := strconv.Atoi(s)
	if err != nil {
		return 0, false
	}
	return v, false
}

// sortedKeys returns the keys of an infos/formats map in
// stable sorted order (used by tests; not part of the
// public API).
//
//nolint:unused // retained for future tests
func sortedKeys[V any](m map[string]V) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// splitStructuredFields splits a key=value comma list
// respecting "..." quoting (commas inside quotes are
// part of the value). Description values like
// "dbSNP membership, build 129" include a comma and
// would be split by strings.Split.
func splitStructuredFields(s, sep string) []string {
	out := []string{}
	cur := []byte{}
	inQuote := false
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c == '"' {
			inQuote = !inQuote
			cur = append(cur, c)
			continue
		}
		if !inQuote && i+len(sep) <= len(s) && s[i:i+len(sep)] == sep {
			out = append(out, string(cur))
			cur = cur[:0]
			i += len(sep) - 1
			continue
		}
		cur = append(cur, c)
	}
	out = append(out, string(cur))
	return out
}

// dictOrder returns the ordered list of FILTER / INFO / FORMAT
// IDs that form the BCF header dictionary, in header-text order.
//
// The spec (§"Dictionary of strings", tex line 2124) defines
// the dictionary as the relative position of each header line:
// "ASP=1, RSPOS=2, dbSNPBuildID=3, and PASS implicitly encoded
// in the first offset PASS=0". The counter is 1-based and spans
// FILTER, INFO and FORMAT lines together; ##contig lines do
// NOT consume an index (contigs are referenced by their own
// 0-based position instead).
//
// Verified against `bcftools view -Ob` output for the spec
// example and for synthetic headers with 0, 1 and 2 contigs.
func (h *Header) dictOrder() []string {
	if h == nil {
		return nil
	}
	if h.dictCache == nil {
		h.dictCache = h.buildDictOrder()
	}
	return h.dictCache
}

// buildDictOrder scans the header text for structured lines.
//
// The index is per *distinct* ID, first occurrence wins: an ID
// used as both an INFO and a FORMAT key shares one index. This
// is what htslib does - for the spec example header it assigns
// NS=1, DP=2, AF=3, AA=4, DB=5, H2=6, q10=7, s50=8, GT=9,
// GQ=10, DP(reuse)=2, HQ=11, so the FORMAT block references DP
// as 2 and HQ as 11.
//
// Verified byte-for-byte against `bcftools view -Ob` output for
// the spec example and for synthetic headers with 0, 1 and 2
// contigs.
func (h *Header) buildDictOrder() []string {
	out := []string{}
	seen := map[string]bool{}
	add := func(id string) {
		if id == "" || seen[id] {
			return
		}
		seen[id] = true
		out = append(out, id)
	}
	if h.Raw != "" {
		for _, line := range strings.Split(h.Raw, "\n") {
			var id string
			switch {
			case strings.HasPrefix(line, "##FILTER="):
				id = extractIDFromAngleField(line, "##FILTER=")
				if id == "PASS" {
					// PASS is implicitly IDX 0 and does not
					// consume a dictionary slot, even when the
					// header declares it explicitly. Verified
					// against `bcftools view -Ob`.
					id = ""
				}
			case strings.HasPrefix(line, "##INFO="):
				id = extractIDFromAngleField(line, "##INFO=")
			case strings.HasPrefix(line, "##FORMAT="):
				id = extractIDFromAngleField(line, "##FORMAT=")
			default:
				continue
			}
			add(id)
		}
		if len(out) > 0 {
			return out
		}
	}
	// Fallback for hand-built headers with no Raw text: use the
	// structured fields in a deterministic order.
	for _, f := range h.Filters {
		if f != "PASS" {
			add(f)
		}
	}
	for _, f := range h.InfosOrder {
		add(f)
	}
	for _, f := range h.FormatsOrder {
		add(f)
	}
	return out
}

// IDXFor returns the 1-based BCF dictionary index of a
// FILTER / INFO / FORMAT ID. "PASS" is implicitly 0.
//
// This is the value the BCF wire format stores for INFO keys,
// FORMAT keys and FILTER offsets.
func (h *Header) IDXFor(id string) int {
	if id == "PASS" {
		return 0
	}
	for i, v := range h.dictOrder() {
		if v == id {
			return i + 1
		}
	}
	return 0
}

// IDAt returns the FILTER / INFO / FORMAT ID with the given
// 1-based dictionary index.
func (h *Header) IDAt(idx int) (string, bool) {
	if idx <= 0 {
		return "", false
	}
	order := h.dictOrder()
	if idx > len(order) {
		return "", false
	}
	return order[idx-1], true
}

// FilterLocalToIDX maps a filter-local index (0 = PASS, 1 =
// the first declared filter) to the BCF dictionary index.
// The two differ because PASS is implicitly 0 and does not
// consume a dictionary slot.
func (h *Header) FilterLocalToIDX(local int) int {
	if local == 0 {
		return 0
	}
	declared := make([]string, 0, len(h.Filters))
	for _, f := range h.Filters {
		if f != "PASS" {
			declared = append(declared, f)
		}
	}
	if local-1 >= len(declared) {
		return 0
	}
	return h.IDXFor(declared[local-1])
}
