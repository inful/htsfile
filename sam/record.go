package sam

import "fmt"

// Record is one alignment line in a SAM file. The 11 mandatory
// columns are stored as scalar fields; optional fields (tags) are
// stored as a slice of [Tag] in declaration order.
//
// The zero value is not meaningful; obtain Records from
// [Decoder.Record].
//
// The parser is lenient about missing columns: a `*` in any
// mandatory column means "not available" per the SAM spec. Numeric
// columns (FLAG, MAPQ) are 0 for unset; Pos / PNext / TLen are 0
// for unmapped segments; Seq / Qual are empty strings for the `*`
// sentinel.
type Record struct {
	// Line is the 1-based source line number where the record was
	// read. The Decoder sets it; callers constructing Records
	// programmatically may leave it zero.
	Line int

	// QName is column 1: query template NAME. Reads/segments
	// having identical QName come from the same template. `*` means
	// unavailable.
	QName string

	// Flag is column 2: bitwise FLAG (paired, mapped, reverse-strand,
	// etc.). See the SAM spec §1.4 for the bit definitions.
	Flag int

	// RName is column 3: reference sequence NAME of the alignment.
	// Must match an @SQ SN if @SQ lines are present and RName is
	// not `*`. An unmapped segment without coordinate has `*`.
	RName string

	// Pos is column 4: 1-based leftmost mapping POSition of the
	// first matching base. 0 means the segment is unmapped.
	Pos int64

	// MapQ is column 5: mapping quality (-10 log10 Pr{mapping
	// position is wrong}, rounded). 255 means "not available".
	MapQ int

	// Cigar is column 6: the CIGAR string. `*` if unavailable.
	// The library stores it verbatim; for semantic interpretation
	// (per-op reference / query consumption, alignment span) use
	// [ParseCigar].
	Cigar string

	// RNext is column 7: reference sequence name of the primary
	// alignment of the NEXT read in the template. `=` means
	// identical to RName; `*` means unavailable.
	RNext string

	// PNext is column 8: position of the primary alignment of the
	// NEXT read in the template. 0 when unavailable.
	PNext int64

	// TLen is column 9: signed observed template length. 0 for
	// single-segment templates or when unavailable.
	TLen int64

	// Seq is column 10: segment sequence. `*` if not stored. If
	// not `*`, the length must equal the sum of M/I/S/=/X operations
	// in the CIGAR.
	Seq string

	// Qual is column 11: ASCII of Phred-scaled base quality + 33.
	// `*` if not stored. If not `*`, the length must equal the length
	// of Seq.
	Qual string

	// Tags is the ordered list of optional fields (TAG:TYPE:VALUE).
	// Multiple tags with the same two-character key are allowed
	// (their order is preserved).
	Tags []Tag
}

// Tag is one optional field on a SAM record (e.g. NM:i:1, MD:Z:36).
// Key is the two-character tag name, Type is the type code, and
// Value is the raw value as it appears in the file.
//
// Value is stored as a string regardless of Type. Typed accessors
// (Int(), Float(), etc.) are intentionally not provided in v1: the
// SAM tag-value syntax varies subtly by type (B arrays are
// comma-separated, H is hex, Z may contain a colon, etc.) and the
// most common use case is preserving the value verbatim for
// round-trip. Callers that need typed access can parse Value
// themselves; future versions may add typed helpers.
type Tag struct {
	Key   string
	Type  TagType
	Value string
}

// TagType is the single-character type code for a SAM tag. The
// valid codes are defined in the SAMtags spec; this type only
// recognises the seven documented codes and rejects unknown ones
// at parse time.
type TagType byte

const (
	// TagA is `A`: printable character.
	TagA TagType = 'A'
	// TagC is `C`: unsigned 8-bit integer (0-255).
	TagC TagType = 'C'
	// Tagf is `f`: single-precision floating-point number.
	Tagf TagType = 'f'
	// TagH is `H`: hex string (case-insensitive).
	TagH TagType = 'H'
	// Tagi is `i`: signed 32-bit integer.
	Tagi TagType = 'i'
	// TagZ is `Z`: printable string.
	TagZ TagType = 'Z'
	// TagB is `B`: general array (subtypes c, C, s, S, i, I, f).
	TagB TagType = 'B'
)

// String returns the single-character type code.
func (t TagType) String() string { return string(byte(t)) }

// ParseTagType returns the TagType for a single character, or an
// error if the character is not one of the seven documented codes.
func ParseTagType(c byte) (TagType, error) {
	switch c {
	case 'A':
		return TagA, nil
	case 'C':
		return TagC, nil
	case 'f':
		return Tagf, nil
	case 'H':
		return TagH, nil
	case 'i':
		return Tagi, nil
	case 'Z':
		return TagZ, nil
	case 'B':
		return TagB, nil
	default:
		return 0, &ParseError{Field: "tag", Msg: fmt.Sprintf("unknown tag type %q", c)}
	}
}

// IsValid reports whether t is one of the seven documented tag
// types.
func (t TagType) IsValid() bool {
	switch t {
	case TagA, TagC, Tagf, TagH, Tagi, TagZ, TagB:
		return true
	}
	return false
}

// Flag bit values from the SAM spec §1.4. Use with [Record.Flag]
// to query or set properties of the alignment.
const (
	FlagPaired        = 1    // template having multiple segments in sequencing
	FlagProperPair    = 2    // each segment properly aligned according to the aligner
	FlagUnmap         = 4    // segment unmapped
	FlagMateUnmap     = 8    // next segment in template unmapped
	FlagReverse       = 16   // SEQ being reverse complemented
	FlagMateReverse   = 32   // SEQ of next segment reverse complemented
	FlagFirst         = 64   // first segment in template
	FlagLast          = 128  // last segment in template
	FlagSecondary     = 256  // secondary alignment
	FlagQCFail        = 512  // not passing quality controls
	FlagDuplicate     = 1024 // PCR or optical duplicate
	FlagSupplementary = 2048 // supplementary alignment
)
