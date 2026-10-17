package sam

// Header is the parsed SAM header — the leading "@" lines that
// precede the alignment section. A header is optional in the
// specification; [Decoder.Header] returns (nil, nil) if the file
// has no "@" lines.
//
// Each structured line type (@HD, @SQ, @RG, @PG) is parsed into a
// dedicated struct with well-known fields surfaced as struct
// members and any other key=value pairs kept in the Other map.
// Every header line also has its Raw text preserved so that
// round-trip emission is byte-identical.
//
// Free-text @CO comment lines are stored as plain strings in
// [Header.COs] in source order.
type Header struct {
	// HD is the file-level @HD line. Nil if absent (the spec
	// recommends it, but a parser must accept files without it).
	HD *HD

	// SQs is the list of @SQ lines in source order. Reference
	// sequence names are also exposed via [Header.LookupSQ].
	SQs []*SQ

	// RGs is the list of @RG lines in source order. Read group
	// IDs are also exposed via [Header.LookupRG].
	RGs []*RG

	// PGs is the list of @PG lines in source order. Program IDs
	// are also exposed via [Header.LookupPG].
	PGs []*PG

	// COs is the list of @CO comment lines (raw, including the
	// leading "@CO\t") in source order. Empty if no comments.
	COs []string
}

// HD is the @HD line — file-level header metadata.
//
// VN is the format version (e.g. "1.6"); SO is the sort order
// ("unknown", "unsorted", "queryname", or "coordinate"); GO is
// the grouping ("none", "query", or "reference").
type HD struct {
	VN    string
	SO    string
	GO    string
	Other map[string]string
	// Raw is the original line text (with the leading "@HD\t"
	// and trailing newline removed). The Encoder emits Raw
	// verbatim when set, falling back to a synthesised line
	// from the structured fields otherwise.
	Raw string
}

// SQ is the @SQ line — a reference sequence in the reference
// dictionary. SN is the sequence name (required); LN is the
// sequence length (required). All other fields are optional.
type SQ struct {
	SN    string
	LN    int64
	AH    string
	AN    string
	AS    string
	M5    string
	UR    string
	SP    string
	Other map[string]string
	Raw   string
}

// RG is the @RG line — a read group. ID is the read group
// identifier (required). All other fields are optional.
type RG struct {
	ID    string
	SM    string
	LB    string
	DS    string
	PU    string
	PI    string
	CN    string
	DT    string
	PL    string
	PG    string
	Other map[string]string
	Raw   string
}

// PG is the @PG line — a program record. ID is the program record
// identifier (required); PN is the program name; CL is the command
// line; PP is the previous @PG ID (for chaining).
type PG struct {
	ID    string
	PN    string
	CL    string
	PP    string
	VN    string
	Other map[string]string
	Raw   string
}

// LookupSQ returns the [SQ] with the given SN, or nil if not
// declared. The second return value is false when the SN is not
// in the header.
func (h *Header) LookupSQ(sn string) (*SQ, bool) {
	if h == nil {
		return nil, false
	}
	for _, sq := range h.SQs {
		if sq.SN == sn {
			return sq, true
		}
	}
	return nil, false
}

// LookupRG returns the [RG] with the given ID, or nil if not
// declared. The second return value is false when the ID is not in
// the header.
func (h *Header) LookupRG(id string) (*RG, bool) {
	if h == nil {
		return nil, false
	}
	for _, rg := range h.RGs {
		if rg.ID == id {
			return rg, true
		}
	}
	return nil, false
}

// LookupPG returns the [PG] with the given ID, or nil if not
// declared. The second return value is false when the ID is not in
// the header.
func (h *Header) LookupPG(id string) (*PG, bool) {
	if h == nil {
		return nil, false
	}
	for _, pg := range h.PGs {
		if pg.ID == id {
			return pg, true
		}
	}
	return nil, false
}
