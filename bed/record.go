package bed

// Record is one data line of a BED file.
//
// Chrom, ChromStart, and ChromEnd are always present (BED3 minimum).
// Optional fields are pointers (for scalars) or slices (for block
// arrays); nil means the field was not present in the source file.
// Round-tripping distinguishes "field absent" from "field present
// with the zero value" — the latter is preserved exactly.
//
// For a BED12 record, the BlockCount pointer and BlockSizes /
// BlockStarts slices must be consistent: if BlockCount is set, the
// two slices must have exactly BlockCount entries; if BlockCount is
// nil, the slices are read-only with respect to the source column
// count.
//
// The zero value is not meaningful; obtain Records from Decoder.Record.
type Record struct {
	// Line is the 1-based source line number where the record was read.
	// Useful for validator issue messages and CLI tooling. The Decoder
	// sets it; callers constructing Records programmatically may leave
	// it zero (the validator will report "line 0").
	Line int

	// Chrom is the chromosome or scaffold name (column 1). Never empty
	// for a well-formed record.
	Chrom string

	// ChromStart is the 0-based inclusive start position (column 2).
	ChromStart int64

	// ChromEnd is the 0-based exclusive end position (column 3).
	// Per the spec, ChromStart <= ChromEnd (a zero-length region
	// is legal and represents an insertion point).
	ChromEnd int64

	// Name is the display name (column 4). Nil means absent.
	Name *string

	// Score is a value in 0..=1000 (column 5). Nil means absent.
	Score *int

	// Strand is the strand (column 6). Nil means absent.
	Strand *Strand

	// ThickStart is the thickly-drawn start (column 7). Nil means
	// absent; when absent, the spec says it defaults to ChromStart.
	ThickStart *int64

	// ThickEnd is the thickly-drawn end (column 8). Nil means
	// absent; when absent, the spec says it defaults to ChromEnd.
	ThickEnd *int64

	// ItemRgb is the RGB color (column 9). Nil means absent.
	ItemRgb *RGB

	// BlockCount is the number of blocks (column 10). Nil means
	// absent. When present, it must equal len(BlockSizes) and
	// len(BlockStarts).
	BlockCount *int

	// BlockSizes is the per-block size list (column 11). Nil means
	// absent. Each entry must be > 0 and the sum of all entries
	// must satisfy:
	//     blockStarts[len-1] + blockSizes[len-1] == chromEnd - chromStart
	BlockSizes []int64

	// BlockStarts is the per-block relative start (column 12).
	// Nil means absent. The first entry must be 0; all entries
	// must be < chromEnd - chromStart; no two blocks may overlap.
	BlockStarts []int64
}

// ColumnCount returns the BED column count (3..=12) implied by the
// populated fields. Per the spec, the optional fields form a strict
// prefix: a BED5 record has Name AND Score; a BED8 record has Name
// through ThickEnd; etc. ColumnCount walks from the highest optional
// field down to find the highest present one and returns that index
// (1-based).
//
// If the prefix is non-contiguous — e.g. Score is set but Name is
// nil — the column count reflects the highest set field, but the
// Validator will report this as an error. ColumnCount itself does
// not validate; it is used to drive round-trip emission and the
// cross-record consistency check.
func (r *Record) ColumnCount() int {
	switch {
	case r.BlockStarts != nil:
		return 12
	case r.BlockSizes != nil:
		return 11
	case r.BlockCount != nil:
		return 10
	case r.ItemRgb != nil:
		return 9
	case r.ThickEnd != nil:
		return 8
	case r.ThickStart != nil:
		return 7
	case r.Strand != nil:
		return 6
	case r.Score != nil:
		return 5
	case r.Name != nil:
		return 4
	default:
		return 3
	}
}

// Strand is the strand of a BED record. The spec only allows three
// values: '+', '-', and '.' (no strand).
type Strand byte

const (
	// StrandPlus is the forward strand ('+').
	StrandPlus Strand = '+'
	// StrandMinus is the reverse strand ('-').
	StrandMinus Strand = '-'
	// StrandNone means no strand information ('.').
	StrandNone Strand = '.'
)

// String returns the canonical single-character representation.
func (s Strand) String() string {
	return string(byte(s))
}

// ParseStrand returns the Strand for a single character, or an error
// if the character is not one of '+', '-', '.'.
func ParseStrand(c byte) (Strand, error) {
	switch c {
	case '+':
		return StrandPlus, nil
	case '-':
		return StrandMinus, nil
	case '.':
		return StrandNone, nil
	default:
		return 0, &ParseError{Field: "strand", Msg: "strand must be '+', '-', or '.'"}
	}
}

// RGB is a 24-bit color (column 9 of BED). Each component is 0..=255.
type RGB struct {
	R, G, B uint8
}

// Track captures the attributes of a BED "track" header line. BED
// files may have zero or one track line; if present it precedes the
// data lines and declares display metadata (name, description, the
// useScore / itemRgb switches, and any number of other
// implementation-defined attributes).
//
// The Raw field holds the original line text so that round-trip
// preserves every attribute — even ones this package does not
// understand — verbatim.
type Track struct {
	// Name is the value of the name= attribute (required by the
	// UCSC track spec for display, but leniently accepted as empty
	// by this parser).
	Name string

	// Description is the value of the description= attribute.
	Description string

	// Type is the value of the type= attribute. The UCSC
	// track spec defines this as one of bedGraph, bed12,
	// bed12WithBlocks, etc.; this package preserves whatever
	// was declared without interpreting it.
	Type string

	// UseScore is true when the track declares useScore=1.
	UseScore bool

	// ItemRgbOn is true when the track declares itemRgb="On".
	ItemRgbOn bool

	// Visibility is the value of the visibility= attribute
	// (typically "dense", "pack", "full", "hide"). Preserved as-is
	// for round-trip; not interpreted.
	Visibility string

	// Raw is the original track line, including the leading
	// "track" keyword and any attributes this struct does not
	// model. The Encoder emits Raw verbatim when set, falling back
	// to a synthesized line from the named fields otherwise.
	Raw string
}

// HasBlocks reports whether r has all three BED12 block fields
// (BlockCount, BlockSizes, BlockStarts) populated.
func (r *Record) HasBlocks() bool {
	return r.BlockCount != nil && r.BlockSizes != nil && r.BlockStarts != nil
}
