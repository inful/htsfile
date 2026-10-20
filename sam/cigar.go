package sam

import (
	"fmt"
	"strconv"
	"strings"
)

// CigarOp is the operation code for one entry in a CIGAR string.
// The nine values defined by SAM v1 §1.4.6 are exposed as named
// constants below. Numeric values match the BAM packed op-code
// convention (M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8) so a
// BAM cigar word (op_len<<4 | op, spec 4.2) can be inspected
// by isolating the low 4 bits and casting to CigarOp.
//
// CigarOp is a single byte; the SAM file format uses the matching
// ASCII character (see Char below).
type CigarOp byte

const (
	// CigarM is M: alignment match (sequence match or mismatch).
	// Consumes query AND reference.
	CigarM CigarOp = 0
	// CigarI is I: insertion to the reference.
	// Consumes query only.
	CigarI CigarOp = 1
	// CigarD is D: deletion from the reference.
	// Consumes reference only.
	CigarD CigarOp = 2
	// CigarN is N: skipped region from the reference (e.g. intron).
	// Consumes reference only.
	CigarN CigarOp = 3
	// CigarS is S: soft clipping. Bases are present in SEQ.
	// Consumes query only.
	CigarS CigarOp = 4
	// CigarH is H: hard clipping. Bases are NOT present in SEQ.
	// Consumes neither query nor reference.
	CigarH CigarOp = 5
	// CigarP is P: padding. Silent deletion from padded reference.
	// Consumes neither query nor reference.
	CigarP CigarOp = 6
	// CigarEq is =: sequence match. Consumes query AND reference.
	// Added in SAM v1.6 to disambiguate from M.
	CigarEq CigarOp = 7
	// CigarX is X: sequence mismatch. Consumes query AND reference.
	// Added in SAM v1.6.
	CigarX CigarOp = 8
)

// String returns the single-character SAM representation of the op
// (e.g. "M", "="). For invalid op values the result is "?".
func (op CigarOp) String() string {
	if c, ok := opChar(op); ok {
		return string(c)
	}
	return "?"
}

// Char returns the ASCII character used in a SAM CIGAR string.
// For invalid op values the result is 0.
func (op CigarOp) Char() byte {
	if c, ok := opChar(op); ok {
		return c
	}
	return 0
}

// IsValid reports whether op is one of the nine documented SAM
// v1 CIGAR operations.
func (op CigarOp) IsValid() bool {
	return op <= CigarX
}

// ConsumesReference reports whether op advances the alignment on
// the reference sequence. True for M, D, N, =, X.
func (op CigarOp) ConsumesReference() bool {
	switch op {
	case CigarM, CigarD, CigarN, CigarEq, CigarX:
		return true
	}
	return false
}

// ConsumesQuery reports whether op advances the alignment on the
// query (SEQ) sequence. True for M, I, S, =, X.
func (op CigarOp) ConsumesQuery() bool {
	switch op {
	case CigarM, CigarI, CigarS, CigarEq, CigarX:
		return true
	}
	return false
}

// opChar maps a CigarOp to its SAM ASCII byte. Returns (0, false)
// for out-of-range op codes.
func opChar(op CigarOp) (byte, bool) {
	switch op {
	case CigarM:
		return 'M', true
	case CigarI:
		return 'I', true
	case CigarD:
		return 'D', true
	case CigarN:
		return 'N', true
	case CigarS:
		return 'S', true
	case CigarH:
		return 'H', true
	case CigarP:
		return 'P', true
	case CigarEq:
		return '=', true
	case CigarX:
		return 'X', true
	}
	return 0, false
}

// parseOpChar maps a single SAM ASCII byte to its CigarOp code.
// Returns (0, false) for any byte that is not one of the nine
// documented CIGAR operations.
func parseOpChar(c byte) (CigarOp, bool) {
	switch c {
	case 'M', 'm':
		return CigarM, true
	case 'I', 'i':
		return CigarI, true
	case 'D', 'd':
		return CigarD, true
	case 'N', 'n':
		return CigarN, true
	case 'S', 's':
		return CigarS, true
	case 'H', 'h':
		return CigarH, true
	case 'P', 'p':
		return CigarP, true
	case '=':
		return CigarEq, true
	case 'X', 'x':
		return CigarX, true
	}
	return 0, false
}

// CigarUnit is one (length, op) pair in a parsed CIGAR string.
// CIGAR strings are sequences of positive integers paired with a
// CIGAR op, e.g. "5M2I3M" parses to [(5,M),(2,I),(3,M)].
type CigarUnit struct {
	Len int32
	Op  CigarOp
}

// ParseCigar parses a SAM CIGAR string into a slice of CigarUnit.
// The string "*" (the SAM "not available" sentinel) returns
// (nil, nil) without error; otherwise an empty string or any
// malformed entry returns an error wrapped in *ParseError with
// the Field set to "cigar".
//
// The parser accepts both upper-case and lower-case ASCII op
// characters and is robust to the common aligner output where 'm'
// (lowercase) appears in extended-CIGAR or clipping contexts.
func ParseCigar(s string) ([]CigarUnit, error) {
	if s == "" {
		return nil, &ParseError{Field: "cigar", Msg: "empty CIGAR string"}
	}
	if s == "*" {
		return nil, nil
	}
	var out []CigarUnit
	i := 0
	for i < len(s) {
		// Read the decimal length. Per spec lengths are
		// non-negative 32-bit integers; we cap at math.MaxInt32.
		j := i
		for j < len(s) && s[j] >= '0' && s[j] <= '9' {
			j++
		}
		if j == i {
			return nil, &ParseError{Field: "cigar",
				Msg: fmt.Sprintf("missing length at offset %d in %q", i, s)}
		}
		n, err := strconv.Atoi(s[i:j])
		if err != nil {
			return nil, &ParseError{Field: "cigar",
				Msg: fmt.Sprintf("invalid length %q in %q: %v", s[i:j], s, err)}
		}
		if n < 0 || n > 1<<31-1 {
			return nil, &ParseError{Field: "cigar",
				Msg: fmt.Sprintf("length %d out of range in %q", n, s)}
		}
		if j >= len(s) {
			return nil, &ParseError{Field: "cigar",
				Msg: fmt.Sprintf("missing op character after %d in %q", n, s)}
		}
		op, ok := parseOpChar(s[j])
		if !ok {
			return nil, &ParseError{Field: "cigar",
				Msg: fmt.Sprintf("unknown op %q at offset %d in %q", s[j], j, s)}
		}
		out = append(out, CigarUnit{Len: int32(n), Op: op}) //nolint:gosec // bounded by the n>1<<31-1 check above
		i = j + 1
	}
	return out, nil
}

// FormatCigar renders a parsed CIGAR back to its textual form.
// A nil or empty input returns "*" (the SAM "not available"
// sentinel).
//
// Unit length 0 is omitted — both because zero-length entries are
// semantically meaningless and because some downstream tools
// reject them. Callers wanting to preserve exact input shape
// should keep their own copy of the original string.
func FormatCigar(ops []CigarUnit) string {
	if len(ops) == 0 {
		return "*"
	}
	var b strings.Builder
	for _, u := range ops {
		if u.Len <= 0 {
			continue
		}
		b.WriteString(strconv.FormatInt(int64(u.Len), 10))
		b.WriteByte(u.Op.Char())
	}
	if b.Len() == 0 {
		return "*"
	}
	return b.String()
}

// ReferenceLength returns the total number of reference bases
// spanned by the CIGAR: the sum of unit lengths for ops that
// consume reference (M, D, N, =, X). Ops that do not consume
// reference (I, S, H, P) contribute 0.
//
// The result equals AlignmentEnd - Pos for a record whose Pos is
// the leftmost reference position. For unmapped records (where
// every op is non-reference-consuming) the result is 0.
func ReferenceLength(ops []CigarUnit) int64 {
	var n int64
	for _, u := range ops {
		if u.Op.ConsumesReference() {
			n += int64(u.Len)
		}
	}
	return n
}

// QueryLength returns the total number of query (SEQ) bases
// spanned by the CIGAR: the sum of unit lengths for ops that
// consume query (M, I, S, =, X). Ops that do not consume query
// (D, N, H, P) contribute 0.
//
// The result equals the length of SEQ for records whose CIGAR is
// well-formed; SAM convention omits hard-clipped bases from SEQ,
// so a "5H10M5H" CIGAR implies SEQ length 10.
func QueryLength(ops []CigarUnit) int64 {
	var n int64
	for _, u := range ops {
		if u.Op.ConsumesQuery() {
			n += int64(u.Len)
		}
	}
	return n
}

// AlignmentEnd returns the 1-based rightmost reference position
// covered by the CIGAR when the alignment starts at pos (1-based,
// per SAM convention). The result is pos + ReferenceLength(ops)
// when ops is non-nil and contains at least one reference-consuming
// op; otherwise the result is pos (which is also the SAM sentinel
// for unmapped segments, where pos is typically 0).
func AlignmentEnd(pos int64, ops []CigarUnit) int64 {
	return pos + ReferenceLength(ops)
}
