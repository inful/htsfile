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
	// Raw is the complete header text, verbatim: one header
	// line per "\n", with a trailing "\n" after the last
	// line. When it is set, [Encoder.EncodeHeader] emits it
	// as-is instead of re-synthesising the lines from the
	// structured fields.
	//
	// This is what a BAM-to-SAM conversion needs: the BAM
	// header stores the original SAM header text, and
	// re-synthesising the lines would drop the fields this
	// package does not model (M5/UR/SP on @SQ, PL/PU/CN on
	// @RG, VN/CL on @PG).
	Raw string

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
	VN string
	SO string
	// GO is the grouping ("none", "query", or "reference").
	// SS is the sub-sort order (SAMv1 §1.3): zero or more
	// ":"-separated keys describing a secondary sort.
	GO    string
	SS    string
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
	SN string
	LN int64
	AH string
	AN string
	AS string
	DS string
	// M5 is the MD5 checksum of the sequence, lowercase hex.
	M5 string
	// TP is the molecule topology: "linear" or "circular".
	// UCSC's recommended-practice rules key off this.
	TP    string
	UR    string
	SP    string
	Other map[string]string
	Raw   string
}

// RG is the @RG line — a read group. ID is the read group
// identifier (required). All other fields are optional.
type RG struct {
	ID string
	SM string
	// LB is the library; PI is the predicted median insert size;
	// PU is the platform unit; CN is the sequencing center;
	// DT is the run date; PL is the platform; PG is the
	// programs group; BC/FO/KS/PM are the barcode-specific
	// fields (SAMv1 §1.3).
	LB    string
	BC    string
	DS    string
	FO    string
	KS    string
	PU    string
	PI    string
	CN    string
	DT    string
	PL    string
	PG    string
	PM    string
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
	DS    string
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
