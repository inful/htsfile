package bam

import (
	"fmt"
)

// Record is one alignment record from a BAM file.
//
// Most scalar fields map 1:1 to their SAM counterparts, but the
// storage shape differs in two important ways:
//
//   - RefID and NextRefID are int32 indices into the header's
//     reference-sequence dictionary (-1 for "no reference", which
//     means RNEXT = "*" in SAM terms), not the sequence names
//     themselves. The decoder sets NameRName / NameRNext from the
//     header so callers don't have to look them up manually.
//   - Cigar is stored as a slice of packed uint32 values, each
//     encoding (op << 4) | len where op is the CIGAR op code
//     (0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 9=B) and len is
//     the operation length. The textual CIGAR string can be
//     reconstructed by [Record.CigarString].
//
// The zero value is not meaningful; obtain Records from
// [Decoder.Record].
type Record struct {
	// Line is the 1-based record number within the alignment
	// section (the first record is 1).
	Line int

	// QName is column 1: query template NAME. `*` means
	// unavailable.
	QName string

	// Flag is column 2: bitwise FLAG. See the SAM spec §1.4 for
	// the bit definitions. Identical bit layout to
	// [sam.FlagPaired] etc.
	Flag uint16

	// RefID is the 0-based index into the BAM header's
	// reference-sequence dictionary (-1 if unmapped). Map to
	// SAM RNAME via [Record.NameRName].
	RefID int32

	// Pos is column 4: 1-based leftmost mapping POSition of the
	// first matching base. 0 means the segment is unmapped. BAM
	// stores 0-based positions internally and the decoder adds 1.
	Pos int64

	// MapQ is column 5: mapping quality. 255 means "not
	// available".
	MapQ uint8

	// Bin is the BAM bin index (a uint16 in [0, 65535]) used
	// by the BAM/BAI indexing scheme to assign a record to one
	// of 37449 + 1 bins. The decoder computes it if it's wrong
	// (some encoders write a stale bin index); the validator
	// checks consistency.
	Bin uint16

	// Cigar is the binary CIGAR representation: each uint32 is
	// (op << 4) | len. See the CIGAR_* constants in this file
	// for the operation codes. Use [Record.CigarString] for the
	// SAM text form.
	Cigar []uint32

	// NextRefID is the reference index of the primary alignment
	// of the NEXT read in the template (-1 if unmapped).
	NextRefID int32

	// NextPos is the 1-based position of the primary alignment
	// of the NEXT read in the template. 0 if unavailable.
	NextPos int64

	// TLen is the signed observed template length. 0 for
	// single-segment templates.
	TLen int64

	// Seq is column 10: segment sequence. "*" if not stored. The
	// length equals the reference length consumed by the CIGAR
	// M/I/S/=/X operations.
	Seq string

	// Qual is column 11: ASCII of Phred-scaled base quality + 33.
	// "*" if not stored. Length equals len(Seq) when set.
	Qual []byte

	// Tags is the ordered list of optional fields. See [Tag] for
	// the value representation. Multiple tags with the same
	// two-character key are allowed.
	Tags []Tag

	// NameRName is the SAM RNAME (sequence name) corresponding
	// to [Record.RefID]. Empty if RefID < 0.
	NameRName string

	// NameRNext is the SAM RNEXT corresponding to
	// [Record.NextRefID]. "*" if NextRefID < 0.
	NameRNext string
}

// Cigar op codes (the high 4 bits of a packed CIGAR uint32).
// Same values as the SAM/BAM specification: M=0, I=1, D=2, N=3,
// S=4, H=5, P=6, ==7, X=8, B=9.
const (
	CigarM  uint32 = 0
	CigarI  uint32 = 1
	CigarD  uint32 = 2
	CigarN  uint32 = 3
	CigarS  uint32 = 4
	CigarH  uint32 = 5
	CigarP  uint32 = 6
	CigarEq uint32 = 7
	CigarX  uint32 = 8
	CigarB  uint32 = 9
)

// CigarOpChar returns the single-character op code for a packed
// CIGAR uint32. Unknown ops (those not listed in the spec) return
// '?'.
func CigarOpChar(p uint32) byte {
	op := p >> 28
	switch op {
	case CigarM:
		return 'M'
	case CigarI:
		return 'I'
	case CigarD:
		return 'D'
	case CigarN:
		return 'N'
	case CigarS:
		return 'S'
	case CigarH:
		return 'H'
	case CigarP:
		return 'P'
	case CigarEq:
		return '='
	case CigarX:
		return 'X'
	case CigarB:
		return 'B'
	default:
		return '?'
	}
}

// CigarOpLen returns the operation length (low 28 bits of the
// packed CIGAR uint32 — BAM uses the high 4 bits for the op
// code and the low 28 bits for the length).
func CigarOpLen(p uint32) uint32 {
	return p & 0x0FFFFFFF
}

// CigarString returns the SAM text form of the CIGAR ("5M2I3M"
// etc.). Returns "*" for an empty CIGAR slice.
func (r *Record) CigarString() string {
	if len(r.Cigar) == 0 {
		return "*"
	}
	buf := make([]byte, 0, len(r.Cigar)*3)
	for _, p := range r.Cigar {
		buf = strconvAppendInt(buf, int64(CigarOpLen(p)))
		buf = append(buf, CigarOpChar(p))
	}
	return string(buf)
}

// strconvAppendInt is a small inline equivalent of
// strconv.AppendInt(buf, n, 10) to avoid pulling strconv into
// the public API surface.
func strconvAppendInt(buf []byte, n int64) []byte {
	if n == 0 {
		return append(buf, '0')
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var tmp [20]byte
	i := len(tmp)
	for n > 0 {
		i--
		tmp[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		tmp[i] = '-'
	}
	return append(buf, tmp[i:]...)
}

// Tag is one optional field on a BAM record (e.g. NM:i:1,
// MD:Z:36). Key is the two-character tag name; Type is the type
// code (see the BAM spec § 1.6 for the type table); Value is the
// value as a string in the same form as SAM would write it.
//
// Type-encoding summary (BAM byte values, not runes):
//
//	A  printable character
//	c  int8
//	C  uint8
//	s  int16 (LE)
//	S  uint16 (LE)
//	i  int32 (LE)
//	I  uint32 (LE)
//	f  float32 (LE)
//	d  float64 (LE)
//	Z  null-terminated string (no terminator stored in Value)
//	H  null-terminated hex string (no terminator stored in Value)
//	B  typed array, subtyped by c/C/S/SiIf
type Tag struct {
	Key   string
	Type  TagType
	Value string
}

// TagType is the single-byte type code for a BAM tag.
type TagType byte

// Valid BAM tag type codes. Naming convention matches the SAM
// package: lowercase = signed type, uppercase = unsigned type.
const (
	TagA TagType = 'A' // printable character
	Tagc TagType = 'c' // int8
	TagC TagType = 'C' // uint8
	Tags TagType = 's' // int16
	TagS TagType = 'S' // uint16
	Tagi TagType = 'i' // int32
	TagI TagType = 'I' // uint32
	Tagf TagType = 'f' // float32
	Tagd TagType = 'd' // float64
	TagZ TagType = 'Z' // null-terminated string
	TagH TagType = 'H' // null-terminated hex
	TagB TagType = 'B' // typed array
)

// IsValid reports whether t is one of the documented BAM tag
// type codes.
func (t TagType) IsValid() bool {
	switch t {
	case TagA, Tagc, TagC, Tags, TagS, Tagi, TagI, Tagf, Tagd, TagZ, TagH, TagB:
		return true
	}
	return false
}

// String returns the single-character type code.
func (t TagType) String() string { return string(byte(t)) }

// ParseTagType returns the TagType for a single byte, or an
// error if the byte is not one of the documented codes.
func ParseTagType(b byte) (TagType, error) {
	t := TagType(b)
	if !t.IsValid() {
		return 0, &ParseError{Field: "tag", Msg: fmt.Sprintf("unknown tag type %q", b)}
	}
	return t, nil
}