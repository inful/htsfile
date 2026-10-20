// Translation helpers for the CLI's sam2bam, bam2sam,
// vcf2bcf, and bcf2vcf subcommands.
//
// These conversions are CLI-level utilities. They translate
// between the project's per-format packages for testing
// purposes — e.g. converting a small SAM text fixture into a
// BAM file (or vice versa) without external tools. They are
// not part of the public API of the per-format packages.
//
// Scope notes:
//   - bam2sam emits SAM with all alignment fields including
//     aux tags. The resulting SAM re-parses cleanly through
//     the existing SAM decoder.
//   - sam2bam emits BAM. Records stream through; the encoder
//     uses default bin assignments (the decoder will read any
//     value, and the validator doesn't check bin).
//   - vcf2bcf and bcf2vcf round-trip header and typed record
//     fields. Per-sample FORMAT data is dropped on vcf2bcf
//     (matching the v1 BCF reader's scope) and synthesised as
//     empty on bcf2vcf.
package main

import (
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/inful/htsfile/bam"
	"github.com/inful/htsfile/bcf"
	"github.com/inful/htsfile/sam"
	"github.com/inful/htsfile/vcf"
)

// ---------------------------------------------------------------------------
// BAM ↔ SAM
// ---------------------------------------------------------------------------

// samToBAMHeader converts a parsed SAM Header to a BAM Header.
// The BAM header carries the verbatim SAM text in hdr.Raw; we
// copy it across and build the SQs list from @SQ lines (with
// any unmapped references appended after).
func samToBAMHeader(h *sam.Header) *bam.Header {
	bh := &bam.Header{}
	if h.HD != nil {
		bh.HD = &bam.HD{VN: h.HD.VN, SO: h.HD.SO, GO: h.HD.GO}
	}
	for _, sq := range h.SQs {
		bh.SQs = append(bh.SQs, &bam.SQ{
			Name: sq.SN, Len: sq.LN,
		})
	}
	for _, rg := range h.RGs {
		bh.RGs = append(bh.RGs, &bam.RG{
			ID: rg.ID, SM: rg.SM, LB: rg.LB, DS: rg.DS,
			PU: rg.PU, PI: rg.PI, CN: rg.CN, DT: rg.DT,
			PL: rg.PL, PG: rg.PG,
		})
	}
	for _, pg := range h.PGs {
		bh.PGs = append(bh.PGs, &bam.PG{
			ID: pg.ID, PN: pg.PN, CL: pg.CL, PP: pg.PP, VN: pg.VN,
		})
	}
	bh.COs = append(bh.COs, h.COs...)
	// Build bh.Raw by re-serialising the structured fields, so
	// the BAM file's SAM-text section actually contains the
	// @HD/@SQ/@RG/@PG/@CO lines. Without this, the BAM file
	// has an empty SAM-text section and the @HD line is lost.
	var sb strings.Builder
	if bh.HD != nil {
		sb.WriteString("@HD")
		if bh.HD.VN != "" {
			sb.WriteString("\tVN:")
			sb.WriteString(bh.HD.VN)
		}
		if bh.HD.SO != "" {
			sb.WriteString("\tSO:")
			sb.WriteString(bh.HD.SO)
		}
		if bh.HD.GO != "" {
			sb.WriteString("\tGO:")
			sb.WriteString(bh.HD.GO)
		}
		sb.WriteByte('\n')
	}
	for _, sq := range bh.SQs {
		sb.WriteString("@SQ\tSN:")
		sb.WriteString(sq.Name)
		sb.WriteString("\tLN:")
		sb.WriteString(strconv.FormatInt(sq.Len, 10))
		sb.WriteByte('\n')
	}
	for _, rg := range bh.RGs {
		sb.WriteString("@RG\tID:")
		sb.WriteString(rg.ID)
		if rg.SM != "" {
			sb.WriteString("\tSM:")
			sb.WriteString(rg.SM)
		}
		if rg.LB != "" {
			sb.WriteString("\tLB:")
			sb.WriteString(rg.LB)
		}
		sb.WriteByte('\n')
	}
	for _, pg := range bh.PGs {
		sb.WriteString("@PG\tID:")
		sb.WriteString(pg.ID)
		if pg.PN != "" {
			sb.WriteString("\tPN:")
			sb.WriteString(pg.PN)
		}
		sb.WriteByte('\n')
	}
	for _, co := range bh.COs {
		sb.WriteString(co)
		sb.WriteByte('\n')
	}
	bh.Raw = sb.String()
	return bh
}

// bamToSAMHeader converts a parsed BAM Header to a SAM Header.
// BAM's hdr.Raw carries the verbatim SAM header text; SAM's
// Header holds the parsed substructures. We populate the
// parsed fields from the BAM header so the resulting SAM
// encoder emits a syntactically equivalent header.
func bamToSAMHeader(h *bam.Header) *sam.Header {
	sh := &sam.Header{
		SQs: bamToSAMSQs(h.SQs),
		RGs: bamToSAMRGs(h.RGs),
		PGs: bamToSAMPGs(h.PGs),
		COs: h.COs,
	}
	if h.HD != nil {
		sh.HD = &sam.HD{VN: h.HD.VN, SO: h.HD.SO, GO: h.HD.GO}
	}
	return sh
}

// bamToSAMSQs converts a []*bam.SQ to []*sam.SQ.
func bamToSAMSQs(in []*bam.SQ) []*sam.SQ {
	if in == nil {
		return nil
	}
	out := make([]*sam.SQ, len(in))
	for i, sq := range in {
		if sq == nil {
			continue
		}
		out[i] = &sam.SQ{SN: sq.Name, LN: sq.Len}
	}
	return out
}

// bamToSAMRGs converts a []*bam.RG to []*sam.RG.
func bamToSAMRGs(in []*bam.RG) []*sam.RG {
	if in == nil {
		return nil
	}
	out := make([]*sam.RG, len(in))
	for i, rg := range in {
		if rg == nil {
			continue
		}
		out[i] = &sam.RG{
			ID: rg.ID, SM: rg.SM, LB: rg.LB, DS: rg.DS,
			PU: rg.PU, PI: rg.PI, CN: rg.CN, DT: rg.DT,
			PL: rg.PL, PG: rg.PG,
		}
	}
	return out
}

// bamToSAMPGs converts a []*bam.PG to []*sam.PG.
func bamToSAMPGs(in []*bam.PG) []*sam.PG {
	if in == nil {
		return nil
	}
	out := make([]*sam.PG, len(in))
	for i, pg := range in {
		if pg == nil {
			continue
		}
		out[i] = &sam.PG{
			ID: pg.ID, PN: pg.PN, CL: pg.CL, PP: pg.PP, VN: pg.VN,
		}
	}
	return out
}

// samToBAMRecord converts a parsed SAM Record to a BAM Record.
// Looks up RefID / NextRefID in the supplied BAM header's SQs.
// Records with unknown contig names produce errors — BAM has
// no way to represent an unknown chromosome name.
func samToBAMRecord(r *sam.Record, hdr *bam.Header) (*bam.Record, error) {
	refID, err := lookupRefID(r.RName, hdr)
	if err != nil {
		return nil, err
	}
	nextRefID, err := lookupRefID(r.RNext, hdr)
	if err != nil {
		return nil, err
	}
	cigar, err := parseCigarString(r.Cigar)
	if err != nil {
		return nil, fmt.Errorf("sam2bam: %w", err)
	}
	// SAM Tag and bam.Tag have identical struct shapes; copy
	// via a typed range loop to translate one to the other.
	bamTags := make([]bam.Tag, len(r.Tags))
	for i, t := range r.Tags {
		bamTags[i] = bam.Tag{
			Key:   t.Key,
			Type:  bam.TagType(t.Type),
			Value: t.Value,
		}
	}
	return &bam.Record{
		QName:     r.QName,
		Flag:      uint16(r.Flag), //nolint:gosec // SAM flag is 12 bits
		RefID:     refID,
		Pos:       r.Pos,
		MapQ:      uint8(r.MapQ), //nolint:gosec // SAM MapQ is 8 bits
		Bin:       computeBin(r.Pos),
		Cigar:     cigar,
		NextRefID: nextRefID,
		NextPos:   r.PNext,
		TLen:      r.TLen,
		Seq:       r.Seq,
		Qual:      []byte(r.Qual),
		Tags:      bamTags,
	}, nil
}

// bamToSAMRecord converts a BAM Record to a SAM Record.
// The BAM Record's NameRName / NameRNext fields already carry
// the contig names cached by the BAM decoder, so no header
// lookup is required here. The packed CIGAR op-codes are
// formatted back to the SAM text form.
func bamToSAMRecord(r *bam.Record) *sam.Record {
	cigar := formatCigarString(r.Cigar)
	samTags := make([]sam.Tag, len(r.Tags))
	for i, t := range r.Tags {
		samTags[i] = sam.Tag{
			Key:   t.Key,
			Type:  sam.TagType(t.Type),
			Value: t.Value,
		}
	}
	return &sam.Record{
		QName: r.QName,
		Flag:  int(r.Flag),
		RName: r.NameRName,
		Pos:   r.Pos,
		MapQ:  int(r.MapQ),
		Cigar: cigar,
		RNext: r.NameRNext,
		PNext: r.NextPos,
		TLen:  r.TLen,
		Seq:   r.Seq,
		Qual:  string(r.Qual),
		Tags:  samTags,
	}
}

// lookupRefID returns the index of contig in hdr.SQs, or -1
// for "*" / empty (unmapped placeholders). An unknown
// contig name returns an error — callers should map
// unknown names to a fresh SQ entry, but for testing
// purposes the round-trip should be lossless.
func lookupRefID(name string, hdr *bam.Header) (int32, error) {
	if name == "" || name == "*" {
		return -1, nil
	}
	for i, sq := range hdr.SQs {
		if sq.Name == name {
			return int32(i), nil
		}
	}
	return -1, fmt.Errorf("unknown reference %q (not in header @SQ list)", name)
}

// computeBin produces a bin index for a 1-based position.
// For Pos <= 0 (unmapped), returns 4680. For Pos > 0, uses
// a simple reg2bins-style level-1 formula (4681 + pos>>26).
// The validator does not check bin so the value is mostly a
// round-trip stability hint.
func computeBin(pos int64) uint16 {
	if pos <= 0 {
		return 4680
	}
	// level-1 bin: 4681 + (pos-1) >> 26, capped at 4688.
	b := 4681 + (pos-1)>>26
	if b > 4688 {
		b = 4688
	}
	return uint16(b) //nolint:gosec // bounded by the check above
}

// formatCigarString converts packed CIGAR op-codes back to
// the SAM text form ("5M2I3M" etc.).
func formatCigarString(packed []uint32) string {
	if len(packed) == 0 {
		return "*"
	}
	var b strings.Builder
	for _, p := range packed {
		op := bam.CigarOpChar(p)
		if op == '?' {
			// Skip unknown ops rather than emit "?".
			op = '?'
		}
		fmt.Fprintf(&b, "%d%c", bam.CigarOpLen(p), op)
	}
	return b.String()
}

// parseCigarString converts a SAM-text CIGAR into the packed
// BAM uint32 form. "" (no CIGAR field) and "*" both map to a
// nil/empty slice; formatCigarString renders those as "*".
// Each op is parsed as <length><char> pairs; an unknown
// char yields an error.
func parseCigarString(s string) ([]uint32, error) {
	if s == "" || s == "*" {
		return nil, nil
	}
	var out []uint32
	for i := 0; i < len(s); {
		// Read the length (one or more digits).
		j := i
		for j < len(s) && s[j] >= '0' && s[j] <= '9' {
			j++
		}
		if j == i {
			return nil, fmt.Errorf("cigar: expected digit at position %d in %q", i, s)
		}
		length, err := strconv.ParseUint(s[i:j], 10, 28)
		if err != nil {
			return nil, fmt.Errorf("cigar: invalid length at position %d: %w", i, err)
		}
		if j >= len(s) {
			return nil, fmt.Errorf("cigar: missing op char after length %d", length)
		}
		op := cigarOpFromChar(s[j])
		if op > 9 {
			return nil, fmt.Errorf("cigar: unknown op %q at position %d", s[j], j)
		}
		// BAM CIGAR packing: op in the high 4 bits, length in
		// the low 28. (See bam.CigarOpChar — it reads op with
		// p >> 28.)
		out = append(out, op<<28|uint32(length))
		i = j + 1
	}
	return out, nil
}

// cigarOpFromChar maps a SAM/BAM CIGAR op character to its
// packed uint32 op code (the high 4 bits). Returns 10 for
// unknown chars (the caller treats this as an error).
func cigarOpFromChar(c byte) uint32 {
	switch c {
	case 'M', 'm':
		return 0
	case 'I', 'i':
		return 1
	case 'D', 'd':
		return 2
	case 'N', 'n':
		return 3
	case 'S', 's':
		return 4
	case 'H', 'h':
		return 5
	case 'P', 'p':
		return 6
	case '=', 'e':
		return 7
	case 'X', 'x':
		return 8
	}
	return 10 // > max valid op (9); triggers error in caller
}

// runTranslateStream implements the standard
// "decode header → convert header → encode header →
// stream records → finalise" pattern shared by all four
// translation subcommands (sam2bam, bam2sam, vcf2bcf,
// bcf2vcf). The four call sites differ only in:
//   - the concrete decoder/encoder types,
//   - the header converter (a no-op for SAM↔BAM but a
//     real type rewrite for VCF↔BCF),
//   - the record converter (which can return an error,
//     notably samToBAMRecord which can fail on mismatched
//     reference IDs),
//   - the encoder's finalise step (Close vs Flush vs
//     "do nothing" — sam.Encoder and vcf.Encoder need
//     Flush, bam.Encoder needs Close, bcf.Encoder is
//     close-on-BGZF-EOF).
//
// All four encoders expose different method names
// (`WriteHeader`/`WriteRecord`/`Close` for bam and bcf,
// `EncodeHeader`/`EncodeRecord`/`Flush` for sam and vcf),
// so the helper is expressed entirely in closures rather
// than via an interface — the per-package naming stays as
// the package authors intended.
func runTranslateStream[DH, EH, DR, ER any](
	readHeader func() (DH, error),
	readRecord func() (DR, error),
	convertHeader func(DH) EH,
	writeHeader func(EH) error,
	convertRecord func(DR, EH) (ER, error),
	writeRecord func(ER) error,
	finalise func() error,
) error {
	hdr, err := readHeader()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	eh := convertHeader(hdr)
	if err := writeHeader(eh); err != nil {
		return fmt.Errorf("write header: %w", err)
	}
	for {
		rec, err := readRecord()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		out, err := convertRecord(rec, eh)
		if err != nil {
			return fmt.Errorf("convert record: %w", err)
		}
		if err := writeRecord(out); err != nil {
			return fmt.Errorf("write record: %w", err)
		}
	}
	if finalise != nil {
		return finalise()
	}
	return nil
}

// runSAMToBAM reads a SAM file from src and writes a BAM
// stream to dst.
func runSAMToBAM(src io.Reader, dst io.Writer) error {
	dec := sam.NewDecoder(src)
	enc := bam.NewEncoder(dst)
	return runTranslateStream(
		dec.Header,
		dec.Record,
		samToBAMHeader,
		enc.WriteHeader,
		func(r *sam.Record, bamHdr *bam.Header) (*bam.Record, error) {
			out, err := samToBAMRecord(r, bamHdr)
			if err != nil {
				return nil, fmt.Errorf("record %q: %w", r.QName, err)
			}
			return out, nil
		},
		enc.WriteRecord,
		func() error {
			if err := enc.Close(); err != nil {
				return fmt.Errorf("close encoder: %w", err)
			}
			return nil
		},
	)
}

// runBAMToSAM reads a BAM stream from src and writes a SAM
// document to dst.
func runBAMToSAM(src io.Reader, dst io.Writer) error {
	dec := bam.NewDecoder(src)
	enc := sam.NewEncoder(dst)
	return runTranslateStream(
		dec.Header,
		dec.Record,
		bamToSAMHeader,
		enc.EncodeHeader,
		func(r *bam.Record, _ *sam.Header) (*sam.Record, error) {
			return bamToSAMRecord(r), nil
		},
		enc.EncodeRecord,
		func() error { return enc.Flush() },
	)
}

// ---------------------------------------------------------------------------
// BCF ↔ VCF
// ---------------------------------------------------------------------------

// vcfToBCFHeader converts a parsed VCF Header to a BCF Header.
// vcf.FileFormat is a struct; BCF uses a string. vcf.Contig.ID
// is the contig name (BCF keeps just a name slice). vcf.InfoDef
// uses typed Number/FieldType; BCF uses stringified versions.
// The translator stringifies the typed values for the BCF
// encoder to consume.
func vcfToBCFHeader(h *vcf.Header) *bcf.Header {
	var sb strings.Builder
	for _, m := range h.MetaLines {
		// m.Raw already includes the leading "##" prefix;
		// only append "\n".
		sb.WriteString(m.Raw)
		sb.WriteByte('\n')
	}
	bh := &bcf.Header{
		FileFormat: h.FileFormat.String(),
		Raw:        sb.String(),
		Infos:      vcfInfoDefsToBCF(h.Infos),
	}
	for _, c := range h.Contigs {
		bh.Contigs = append(bh.Contigs, c.ID)
	}
	bh.SampleNames = append([]string(nil), h.SampleNames...)
	// Build InfosOrder by iterating vcf.Header's MetaLines in
	// source order — each ##INFO line adds its ID at the end.
	// The bcf encoder uses InfosOrder to assign INFO key
	// indices; without it, no INFO bytes are emitted at all.
	for _, m := range h.MetaLines {
		if m.Key != "INFO" {
			continue
		}
		// Value looks like "<ID=DP,Number=1,Type=Integer,Description=\"depth\">".
		const prefix = "<ID="
		if !strings.HasPrefix(m.Value, prefix) {
			continue
		}
		rest := m.Value[len(prefix):]
		end := strings.IndexByte(rest, ',')
		if end < 0 {
			end = strings.IndexByte(rest, '>')
			if end < 0 {
				end = len(rest)
			}
		}
		bh.InfosOrder = append(bh.InfosOrder, rest[:end])
	}
	return bh
}

// _ ensures sam and vcf packages are referenced.
var _ = sam.Header{}
var _ = vcf.Header{}

// vcfInfoDefsToBCF copies a VCF InfoDef map to a BCF InfoDef
// map, stringifying the typed Number / Type fields. Returns
// nil for nil input.
func vcfInfoDefsToBCF(m map[string]*vcf.InfoDef) map[string]*bcf.InfoDef {
	if m == nil {
		return nil
	}
	out := make(map[string]*bcf.InfoDef, len(m))
	for k, v := range m {
		if v == nil {
			out[k] = nil
			continue
		}
		out[k] = &bcf.InfoDef{
			ID:          v.ID,
			Number:      v.Number.String(),
			Type:        v.Type.String(),
			Description: v.Description,
		}
	}
	return out
}

// bcfToVCFHeader converts a parsed BCF Header to a VCF Header.
// The reverse mapping: BCF string-typed Number / Type become
// vcf.Number / vcf.FieldType via their parsers.
func bcfToVCFHeader(h *bcf.Header) *vcf.Header {
	vh := &vcf.Header{
		FileFormat: parseVCFFileFormat(h.FileFormat),
		Infos:      bcfInfosToVCF(h.Infos),
	}
	for _, name := range h.Contigs {
		vh.Contigs = append(vh.Contigs, &vcf.Contig{ID: name})
	}
	vh.SampleNames = append([]string(nil), h.SampleNames...)
	return vh
}

// parseVCFFileFormat parses a "VCFvMAJOR.MINOR" string into a
// vcf.FileFormat struct. Empty string yields the zero value
// (Major=0, Minor=0). Unknown formats yield the zero value.
func parseVCFFileFormat(s string) vcf.FileFormat {
	ff, err := vcf.ParseFileFormat(s)
	if err != nil {
		return vcf.FileFormat{}
	}
	return ff
}

// bcfInfosToVCF copies a BCF InfoDef map to a VCF InfoDef map.
func bcfInfosToVCF(m map[string]*bcf.InfoDef) map[string]*vcf.InfoDef {
	if m == nil {
		return nil
	}
	out := make(map[string]*vcf.InfoDef, len(m))
	for k, v := range m {
		if v == nil {
			out[k] = nil
			continue
		}
		nv := &vcf.InfoDef{
			ID:          v.ID,
			Description: v.Description,
		}
		// Parse the stringified Number / Type back into the
		// typed vcf.Number / vcf.FieldType. Best-effort:
		// unparseable values produce empty fields.
		if v.Number != "" {
			if n, err := vcf.ParseNumber(v.Number); err == nil {
				nv.Number = n
			}
		}
		if v.Type != "" {
			if t, err := vcf.ParseFieldType(v.Type); err == nil {
				nv.Type = t
			}
		}
		out[k] = nv
	}
	return out
}

// vcfToBCFRecord converts a parsed VCF Record to a BCF Record.
// Drops per-sample FORMAT data (matching the v1 BCF reader's
// scope). The vcf.InfoValue carries an unexported typedCore
// (Type, Ints, Floats, Strs, Chars, Missing, Raw); the BCF
// equivalent is exported (Type string, Ints, Floats, Strs, Str,
// Flag). We map Character (Chars) to the single-string Str
// field; Missing / Raw are dropped (BCF uses a different
// sentinel convention).
//
// vcf.Record.ID follows the VCF convention where "." means
// "no ID available". BCF's wire format uses typeMissing for
// this; we normalise both "." and an empty slice to the
// missing-ID form here so the encoder emits typeMissing.
func vcfToBCFRecord(r *vcf.Record) *bcf.Record {
	hasID := false
	for _, id := range r.ID {
		if id != "." && id != "" {
			hasID = true
			break
		}
	}
	var idStr string
	if hasID {
		idStr = strings.Join(r.ID, ",")
	}
	info := make(map[string]bcf.InfoValue, len(r.Info))
	for k, v := range r.Info { //nolint:gocritic // rangeValCopy: InfoValue is a value type by design
		bv := bcf.InfoValue{
			Type:   v.Type.String(),
			Ints:   append([]int64(nil), v.Ints...),
			Floats: append([]float64(nil), v.Floats...),
			Strs:   append([]string(nil), v.Strs...),
			Flag:   v.Flag,
		}
		// Character fields carry a single byte in vcf.InfoValue.Chars.
		if len(v.Chars) > 0 {
			bv.Str = string(v.Chars)
		}
		info[k] = bv
	}
	return &bcf.Record{
		Line:      r.Line,
		CHROM:     r.Chrom,
		Pos:       r.Pos,
		ID:        idStr,
		Ref:       r.Ref,
		Alt:       append([]string(nil), r.Alt...),
		Qual:      vcfQualToBCF(r.Qual),
		QualKnown: r.Qual != nil,
		Filter:    append([]string(nil), r.Filter...),
		Info:      info,
	}
}

// vcfQualToBCF extracts the float64 quality from a VCF Qual
// pointer. Nil means missing.
func vcfQualToBCF(q *float64) float64 {
	if q == nil {
		return 0
	}
	return *q
}

// bcfToVCFRecord converts a parsed BCF Record to a VCF Record.
// Per-sample FORMAT data is not synthesised (BCF v1 doesn't
// carry it). vcf.InfoValue uses unexported typedCore, so we
// set the values via the typedCore path.
func bcfToVCFRecord(r *bcf.Record) *vcf.Record {
	var idList []string
	if r.ID != "" {
		idList = strings.Split(r.ID, ",")
	}
	info := make(map[string]vcf.InfoValue, len(r.Info))
	for k, v := range r.Info {
		t, err := vcf.ParseFieldType(v.Type)
		if err != nil {
			t = vcf.TypeUnknown
		}
		vv := vcf.InfoValue{}
		vv.Type = t
		vv.Ints = append([]int64(nil), v.Ints...)
		vv.Floats = append([]float64(nil), v.Floats...)
		vv.Strs = append([]string(nil), v.Strs...)
		// Character (single character) is stored as a 1-byte
		// Chars slice on the VCF side; treat any TypeChar
		// value with Str as such.
		if t == vcf.TypeCharacter && v.Str != "" {
			vv.Chars = []byte(v.Str)
		}
		vv.Flag = v.Flag
		info[k] = vv
	}
	return &vcf.Record{
		Line:   r.Line,
		Chrom:  r.CHROM,
		Pos:    r.Pos,
		ID:     idList,
		Ref:    r.Ref,
		Alt:    append([]string(nil), r.Alt...),
		Qual:   bcfQualToVCF(r.Qual, r.QualKnown),
		Filter: append([]string(nil), r.Filter...),
		Info:   info,
	}
}

// bcfQualToVCF converts a (qual, qualKnown) pair to a VCF
// pointer. Missing (qualKnown=false) becomes a nil pointer.
func bcfQualToVCF(qual float64, qualKnown bool) *float64 {
	if !qualKnown {
		return nil
	}
	q := qual
	return &q
}

// runVCFToBCF reads a VCF file from src and writes a BCF
// stream to dst.
func runVCFToBCF(src io.Reader, dst io.Writer) error {
	dec := vcf.NewDecoder(src)
	enc := bcf.NewEncoder(dst)
	return runTranslateStream(
		dec.Header,
		dec.Record,
		vcfToBCFHeader,
		enc.WriteHeader,
		func(r *vcf.Record, _ *bcf.Header) (*bcf.Record, error) {
			return vcfToBCFRecord(r), nil
		},
		enc.WriteRecord,
		func() error {
			if err := enc.Close(); err != nil {
				return fmt.Errorf("close encoder: %w", err)
			}
			return nil
		},
	)
}

// runBCFToVCF reads a BCF stream from src and writes a VCF
// document to dst.
func runBCFToVCF(src io.Reader, dst io.Writer) error {
	dec := bcf.NewDecoder(src)
	enc := vcf.NewEncoder(dst)
	return runTranslateStream(
		dec.Header,
		dec.Record,
		bcfToVCFHeader,
		enc.EncodeHeader,
		func(r *bcf.Record, _ *vcf.Header) (*vcf.Record, error) {
			return bcfToVCFRecord(r), nil
		},
		enc.EncodeRecord,
		func() error { return enc.Flush() },
	)
}
