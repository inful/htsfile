package vcf

import (
	"fmt"
	"strconv"
	"strings"
)

// FileFormat records the VCF version reported by the leading ##fileformat line,
// for example "##fileformat=VCFv4.5".
//
// The library accepts any v4.x minor version (v4.0 through v4.5 are all
// known-good; future v4.6+ releases will be accepted as long as they
// keep the v4 major). Major versions other than 4 (v3, v5+) are
// rejected. Per-version validation rule differences are not yet
// implemented — the parser and writer treat all v4.x identically,
// preserving the FileFormat verbatim on round-trip.
type FileFormat struct {
	Major int
	Minor int
}

// String returns the canonical "VCFvMAJOR.MINOR" rendering.
func (f FileFormat) String() string {
	return "VCFv" + strconv.Itoa(f.Major) + "." + strconv.Itoa(f.Minor)
}

// ParseFileFormat parses "VCFvMAJOR.MINOR" (case sensitive on "VCF", case
// insensitive on the version digits). The parser enforces major version
// 4 only — VCFv3 and VCFv5+ require a separate code path that this
// library does not implement. Any v4.x minor version is accepted.
func ParseFileFormat(s string) (FileFormat, error) {
	if !strings.HasPrefix(s, "VCFv") && !strings.HasPrefix(s, "vcfv") {
		return FileFormat{}, &ParseError{Field: "fileformat", Msg: "missing 'VCFv' prefix in " + s}
	}
	rest := s[4:]
	dot := strings.Index(rest, ".")
	if dot < 0 {
		return FileFormat{}, &ParseError{Field: "fileformat", Msg: "missing '.' in version " + s}
	}
	major, err := strconv.Atoi(rest[:dot])
	if err != nil {
		return FileFormat{}, &ParseError{Field: "fileformat", Msg: "non-numeric major in " + s}
	}
	minor, err := strconv.Atoi(rest[dot+1:])
	if err != nil {
		return FileFormat{}, &ParseError{Field: "fileformat", Msg: "non-numeric minor in " + s}
	}
	if major != 4 {
		return FileFormat{}, &ParseError{
			Field: "fileformat",
			Msg:   fmt.Sprintf("unsupported VCF major version %d (only v4.x is supported)", major),
		}
	}
	return FileFormat{Major: major, Minor: minor}, nil
}

// Header is the parsed header of a VCF file — the ## meta lines and the
// #CHROM header line. A header must be processed before any [Record] can be
// interpreted (records reference header-declared contigs, INFO/FORMAT keys,
// FILTER codes, etc.).
//
// The maps are populated only for well-known structured lines (##INFO,
// ##FORMAT, ##FILTER, ##ALT). All other ## lines, whether structured or
// unstructured, are kept in [Header.MetaLines] in source order so that
// [Encoder.EncodeHeader] can preserve them verbatim.
type Header struct {
	// FileFormat is the value of the leading ##fileformat line. It is
	// always populated; an empty value indicates the line was missing.
	FileFormat FileFormat

	// FileDate is the value of ##fileDate, or empty if absent.
	FileDate string
	// Source is the value of ##source, or empty if absent.
	Source string
	// Reference is the value of ##reference, or empty if absent.
	Reference string
	// Assembly is the value of ##assembly, or empty if absent.
	Assembly string
	// PedigreeDB is the value of ##pedigreeDB, or empty if absent.
	PedigreeDB string

	// Contigs is the list of ##contig entries in source order.
	Contigs []*Contig

	// Infos maps ID → definition for all ##INFO entries.
	Infos map[string]*InfoDef
	// Filters maps ID → definition for all ##FILTER entries.
	Filters map[string]*FilterDef
	// Formats maps ID → definition for all ##FORMAT entries.
	Formats map[string]*FormatDef
	// Alts maps ID → definition for all ##ALT entries.
	Alts map[string]*AltDef
	// Metas maps ID → definition for all ##META entries. Used by the
	// validator to constrain ##SAMPLE field values.
	Metas map[string]*MetaDef
	// Pedigrees lists all ##PEDIGREE entries in source order.
	Pedigrees []*Pedigree
	// Samples lists all ##SAMPLE entries in source order.
	Samples []*Sample

	// SampleNames lists the sample identifiers from the #CHROM header line.
	// Empty if the file has no FORMAT column / no samples.
	SampleNames []string

	// MetaLines holds every ## meta-information line in source order,
	// including the structured ones (which are also surfaced into the typed
	// maps above). Each entry preserves the original raw text after the
	// leading "##" so that [Encoder.EncodeHeader] can replay it.
	MetaLines []MetaLine
}

// MetaLine is a single ##key=value or ##key=<...> meta line.
//
// Raw is the original text following the leading "##" with its trailing
// newline removed. Key is the value before the first '=' (for structured
// lines the enclosing "<...>" is included on the Value side only when the
// line is structured).
type MetaLine struct {
	Key   string
	Value string
	Raw   string
}

// Contig describes one ##contig=<...> entry.
type Contig struct {
	ID       string
	Length   int64  // 0 if absent
	MD5      string // "" if absent
	URL      string // "" if absent
	Species  string // "" if absent (e.g. "Homo sapiens")
	Taxonomy string // "" if absent (e.g. "NCBI:9606")
	// Other holds any unknown key=value pairs from this contig line.
	Other map[string]string
}

// InfoDef describes one ##INFO=<ID=...,Number=...,Type=...,Description=...>
// entry.
type InfoDef struct {
	ID          string
	Number      Number
	Type        FieldType
	Description string
	// Source is the optional Source="..." field; "" if absent.
	Source string
	// Version is the optional Version="..." field; "" if absent.
	Version string
	// Other holds any unknown key=value pairs from this line.
	Other map[string]string
}

// FormatDef describes one ##FORMAT=<ID=...,Number=...,Type=...,Description=...>
// entry.
type FormatDef struct {
	ID          string
	Number      Number
	Type        FieldType
	Description string
	// Other holds any unknown key=value pairs from this line.
	Other map[string]string
}

// FilterDef describes one ##FILTER=<ID=...,Description=...> entry.
type FilterDef struct {
	ID          string
	Description string
	// Other holds any unknown key=value pairs from this line.
	Other map[string]string
}

// AltDef describes one ##ALT=<ID=...,Description=...> entry.
//
// ID may contain the structured form "TYPE:SUBTYPE" (e.g. "CNV:TR",
// "DUP:TANDEM"). The first segment before ":" is the primary type.
type AltDef struct {
	ID          string
	Description string
	// Other holds any unknown key=value pairs from this line.
	Other map[string]string
}

// Pedigree describes one ##PEDIGREE=<ID=...> entry.
//
// The captured fields cover the common cases listed in the spec
// (Original=, Father=, Mother=, Name_1=…Name_N=). Other custom keys are
// available through Pedigree.Other.
type Pedigree struct {
	ID       string
	Original string
	Father   string
	Mother   string
	Names    []string // ordered list of Name_N values
	Other    map[string]string
}

// Sample describes one ##SAMPLE=<ID=...> entry.
type Sample struct {
	ID          string
	Genomes     map[string]string // Genomes="..." (rarely used)
	MolConvs    map[string]string // MolConvs="..."
	Description string
	DOI         string
	Other       map[string]string
}

// MetaDef describes one ##META=<ID=...,Type=...,Number=...,Values=[...]>
// entry. The Values list, when non-nil, constrains the values that may
// appear on a ##SAMPLE line under this key (per VCFv4.5 §1.4.7
// "Sample field format").
type MetaDef struct {
	ID     string
	Type   FieldType
	Number Number
	// Values is the parsed [a, b, c] list. nil means no Values was
	// declared (the constraint is not enforced).
	Values []string
}

// LookupInfo returns the [InfoDef] for the given ID, or nil if not declared.
// The second return value is false when the ID is not in the header.
func (h *Header) LookupInfo(id string) (*InfoDef, bool) {
	if h == nil || h.Infos == nil {
		return nil, false
	}
	d, ok := h.Infos[id]
	return d, ok
}

// LookupFormat returns the [FormatDef] for the given ID, or nil if not
// declared.
func (h *Header) LookupFormat(id string) (*FormatDef, bool) {
	if h == nil || h.Formats == nil {
		return nil, false
	}
	d, ok := h.Formats[id]
	return d, ok
}

// LookupFilter returns the [FilterDef] for the given ID, or nil if not
// declared.
func (h *Header) LookupFilter(id string) (*FilterDef, bool) {
	if h == nil || h.Filters == nil {
		return nil, false
	}
	d, ok := h.Filters[id]
	return d, ok
}

// LookupAlt returns the [AltDef] for the given ID, or nil if not declared.
func (h *Header) LookupAlt(id string) (*AltDef, bool) {
	if h == nil || h.Alts == nil {
		return nil, false
	}
	d, ok := h.Alts[id]
	return d, ok
}

// LookupContig returns the [Contig] for the given ID, or nil if not declared.
func (h *Header) LookupContig(id string) (*Contig, bool) {
	if h == nil {
		return nil, false
	}
	for _, c := range h.Contigs {
		if c.ID == id {
			return c, true
		}
	}
	return nil, false
}
