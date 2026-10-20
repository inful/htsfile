package bam

import "strconv"

// Header is the parsed BAM header. The binary BAM header carries:
//
//   - A SAM-formatted header text block (Raw)
//   - A reference-sequence dictionary (SQs)
//
// The Decoder additionally parses the embedded SAM text into the
// HD / RGs / PGs / COs fields so callers don't have to re-parse
// it themselves. The parsed fields mirror the [sam.Header]
// structure but are independent types — per the per-package
// self-contained design principle.
type Header struct {
	// Raw is the embedded SAM header text exactly as it appears
	// in the BAM binary header, including embedded newlines but
	// without the trailing null terminator.
	Raw string

	// SQs is the reference-sequence dictionary parsed from the
	// binary portion of the BAM header (after the SAM text). It
	// is also accessible via [Header.LookupSQ].
	SQs []*SQ

	// HD is the file-level @HD line parsed from Raw. Nil if Raw
	// has no @HD line.
	HD *HD

	// RGs is the list of @RG lines parsed from Raw, in source
	// order. Read-group IDs are also accessible via
	// [Header.LookupRG].
	RGs []*RG

	// PGs is the list of @PG lines parsed from Raw, in source
	// order. Program IDs are also accessible via
	// [Header.LookupPG].
	PGs []*PG

	// COs is the list of @CO comment lines parsed from Raw, in
	// source order (raw text including the leading "@CO\t").
	COs []string
}

// HD is the @HD line — file-level header metadata. VN is the
// format version (e.g. "1.6"); SO is the sort order; GO is the
// grouping.
type HD struct {
	VN string
	SO string
	// GO is the grouping; SS is the sub-sort order (SAMv1 §1.3).
	GO    string
	SS    string
	Other map[string]string
	Raw   string
}

// SQ is a reference-sequence entry from the BAM binary reference
// dictionary. Name is the SN (sequence name); Len is the LN
// (sequence length).
type SQ struct {
	Name string
	Len  int64
	// The remaining fields are populated from the binary
	// reference dictionary when present. DS is the
	// description; TP is the topology ("linear" or
	// "circular").
	DS string
	TP string
}

// RG is the @RG line — a read group. ID is the read group
// identifier (required); the other fields are optional.
type RG struct {
	ID    string
	SM    string
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

// PG is the @PG line — a program record. ID is the program
// identifier (required); the other fields are optional.
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

// LookupSQ returns the [SQ] with the given Name, or nil if not
// declared. The second return value is false when Name is not in
// the dictionary.
func (h *Header) LookupSQ(name string) (*SQ, bool) {
	if h == nil {
		return nil, false
	}
	for _, sq := range h.SQs {
		if sq.Name == name {
			return sq, true
		}
	}
	return nil, false
}

// LookupRG returns the [RG] with the given ID, or nil if not
// declared. The second return value is false when the ID is not
// in the header.
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
// declared. The second return value is false when the ID is not
// in the header.
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

// String returns a short one-line representation of the SQ.
// Used by tests and debug output.
func (s *SQ) String() string {
	if s == nil {
		return ""
	}
	return s.Name + ":" + strconv.FormatInt(s.Len, 10)
}
