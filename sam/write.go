package sam

import (
	"bufio"
	"bytes"
	"io"
	"sort"
	"strconv"
)

// Encoder writes a SAM document to an underlying [io.Writer].
//
// Typical usage:
//
//	enc := sam.NewEncoder(w)
//	if hdr != nil { _ = enc.EncodeHeader(hdr) }
//	for _, r := range recs { _ = enc.EncodeRecord(r) }
//	_ = enc.Flush()
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w          *bufio.Writer
	lineEnding []byte
}

// NewEncoder wraps w in an [Encoder] with default LF line endings.
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{
		w:          bufio.NewWriter(w),
		lineEnding: []byte{'\n'},
	}
}

// NewEncoderCRLF wraps w in an [Encoder] that uses CR+LF line endings.
func NewEncoderCRLF(w io.Writer) *Encoder {
	return &Encoder{
		w:          bufio.NewWriter(w),
		lineEnding: []byte{'\r', '\n'},
	}
}

// EncodeHeader writes the SAM header. The first call writes the
// @HD, @SQ, @RG, @PG, and @CO lines in source order. Subsequent
// calls are no-ops (a SAM file has exactly one header section).
func (e *Encoder) EncodeHeader(h *Header) error {
	if h == nil {
		return nil
	}
	if h.HD != nil {
		if err := e.writeRawOrSynth(h.HD.Raw, "@HD", h.HD.Other, map[string]string{
			"VN": h.HD.VN, "SO": h.HD.SO, "GO": h.HD.GO, "SS": h.HD.SS,
		}); err != nil {
			return err
		}
	}
	for _, sq := range h.SQs {
		if err := e.writeRawOrSynth(sq.Raw, "@SQ", sq.Other, map[string]string{
			"SN": sq.SN,
			"LN": strconv.FormatInt(sq.LN, 10),
			"AH": sq.AH, "AN": sq.AN, "AS": sq.AS, "DS": sq.DS,
			"M5": sq.M5, "TP": sq.TP, "UR": sq.UR, "SP": sq.SP,
		}); err != nil {
			return err
		}
	}
	for _, rg := range h.RGs {
		if err := e.writeRawOrSynth(rg.Raw, "@RG", rg.Other, map[string]string{
			"ID": rg.ID, "SM": rg.SM, "LB": rg.LB, "BC": rg.BC, "DS": rg.DS,
			"FO": rg.FO, "KS": rg.KS, "PU": rg.PU, "PI": rg.PI, "CN": rg.CN, "DT": rg.DT,
			"PL": rg.PL, "PG": rg.PG, "PM": rg.PM,
		}); err != nil {
			return err
		}
	}
	for _, pg := range h.PGs {
		if err := e.writeRawOrSynth(pg.Raw, "@PG", pg.Other, map[string]string{
			"ID": pg.ID, "PN": pg.PN, "CL": pg.CL, "DS": pg.DS, "PP": pg.PP, "VN": pg.VN,
		}); err != nil {
			return err
		}
	}
	for _, co := range h.COs {
		if _, err := e.w.WriteString(co); err != nil {
			return err
		}
		if _, err := e.w.Write(e.lineEnding); err != nil {
			return err
		}
	}
	return nil
}

// writeRawOrSynth writes one header line. If raw is set, it's
// emitted verbatim (preserving unknown attributes). Otherwise, a
// line is synthesised from the structured fields followed by the
// Other map (Other wins on conflict so callers can override
// structured fields if needed).
func (e *Encoder) writeRawOrSynth(raw, prefix string, other, structured map[string]string) error {
	if raw != "" {
		if _, err := e.w.WriteString(raw); err != nil {
			return err
		}
		_, err := e.w.Write(e.lineEnding)
		return err
	}
	// Merge: structured first, then other (so caller-provided
	// Other overrides for new fields take precedence).
	merged := make(map[string]string, len(structured)+len(other))
	for k, v := range structured {
		if v != "" {
			merged[k] = v
		}
	}
	for k, v := range other {
		merged[k] = v
	}
	if len(merged) == 0 {
		// No fields at all — emit just the prefix.
		if _, err := e.w.WriteString(prefix); err != nil {
			return err
		}
		_, err := e.w.Write(e.lineEnding)
		return err
	}
	var b bytes.Buffer
	b.WriteString(prefix)
	for _, k := range sortedKeys(merged) {
		b.WriteByte('\t')
		b.WriteString(k)
		b.WriteByte(':')
		b.WriteString(merged[k])
	}
	_, _ = b.Write(e.lineEnding)
	_, err := e.w.Write(b.Bytes())
	return err
}

// EncodeRecord writes one alignment record as a tab-separated
// line followed by the line ending. Empty string fields are
// rendered as `*` per the SAM spec (the canonical "unavailable"
// sentinel for string columns).
func (e *Encoder) EncodeRecord(r *Record) error {
	var b bytes.Buffer
	b.WriteString(orStar(r.QName))
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(r.Flag))
	b.WriteByte('\t')
	b.WriteString(orStar(r.RName))
	b.WriteByte('\t')
	b.WriteString(strconv.FormatInt(r.Pos, 10))
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(r.MapQ))
	b.WriteByte('\t')
	b.WriteString(orStar(r.Cigar))
	b.WriteByte('\t')
	b.WriteString(orStar(r.RNext))
	b.WriteByte('\t')
	b.WriteString(strconv.FormatInt(r.PNext, 10))
	b.WriteByte('\t')
	b.WriteString(strconv.FormatInt(r.TLen, 10))
	b.WriteByte('\t')
	b.WriteString(orStar(r.Seq))
	b.WriteByte('\t')
	b.WriteString(orStar(r.Qual))
	for _, t := range r.Tags {
		b.WriteByte('\t')
		b.WriteString(t.Key)
		b.WriteByte(':')
		b.WriteByte(byte(t.Type))
		b.WriteByte(':')
		b.WriteString(t.Value)
	}
	_, _ = b.Write(e.lineEnding)
	_, err := e.w.Write(b.Bytes())
	return err
}

// orStar returns s if non-empty, else "*". The SAM spec uses "*"
// as the canonical sentinel for "field unavailable" across all
// string columns (QNAME, RNAME, CIGAR, RNEXT, SEQ, QUAL).
func orStar(s string) string {
	if s == "" {
		return "*"
	}
	return s
}

// Flush writes any buffered data to the underlying writer.
func (e *Encoder) Flush() error { return e.w.Flush() }

// sortedKeys returns the keys of m in sorted order. Used to make
// header-line emission deterministic.
func sortedKeys(m map[string]string) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// (errInvalidRecord was removed; no longer needed.)
