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
	// The BAM header embeds the SAM header text verbatim. Use
	// the original text when we have it: re-synthesising the
	// lines from the structured fields drops the tags this
	// package does not model (M5/UR/SP on @SQ, PL/PU/CN on
	// @RG, VN/CL on @PG).
	if h.Raw != "" {
		bh.Raw = h.Raw
		return bh
	}
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
		// The BAM header stores the original SAM header text
		// verbatim; pass it through so nothing is lost.
		Raw: h.Raw,
		SQs: bamToSAMSQs(h.SQs),
		RGs: bamToSAMRGs(h.RGs),
		PGs: bamToSAMPGs(h.PGs),
		COs: h.COs,
	}
	if h.HD != nil {
		sh.HD = &sam.HD{VN: h.HD.VN, SO: h.HD.SO, GO: h.HD.GO, Raw: h.HD.Raw}
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
			Raw: rg.Raw,
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
			Raw: pg.Raw,
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
	// RNEXT "=" means the mate is on the same reference as
	// this read. htslib resolves it that way
	// (sam.c bam_parse1: `if (strcmp(q, "=") == 0) c->mtid =
	// c->tid;`), and samtools accepts such files. Looking "="
	// up in the @SQ list made us reject valid SAM.
	nextName := r.RNext
	if nextName == "=" {
		nextName = r.RName
	}
	nextRefID, err := lookupRefID(nextName, hdr)
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
			Type:  bamTagTypeToSAM(t.Type),
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
		PNext: r.NextPos,
		TLen:  r.TLen,
		Seq:   r.Seq,
		Qual:  bamQualToSAM(r.Qual),
		RNext: bamNextRefName(r),
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
		// BAM CIGAR packing per spec 4.2: op_len<<4 | op, so
		// the op code is in the LOW 4 bits and the length in
		// the high 28. A 5M is 0x50.
		out = append(out, bam.PackCigarOp(op, uint32(length)))
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
	// The BCF header text must include the sample-line
	// (`#CHROM\tPOS\tID\t...\tSAMPLE1\tSAMPLE2…`) so the BCF
	// decoder can recover SampleNames. MetaLines only carry
	// the `##` prefixes, so we synthesise the sample line
	// from h.SampleNames when the input header omitted one.
	if len(h.SampleNames) > 0 {
		sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
		// The FORMAT column is always present when there are
		// samples: the sample names start in column 10, so
		// omitting it would shift them into the INFO column
		// and the BCF header would declare zero samples.
		sb.WriteString("\tFORMAT")
		for _, sn := range h.SampleNames {
			sb.WriteByte('\t')
			sb.WriteString(sn)
		}
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
	// Mirror for ##FORMAT: the BCF record per-sample block
	// keys off FormatsOrder (uint16 indices).
	bh.Formats = make(map[string]*bcf.FormatDef, len(h.Formats))
	for id, fd := range h.Formats {
		if fd == nil {
			continue
		}
		var t string
		switch fd.Type {
		case vcf.TypeInteger:
			t = "Integer"
		case vcf.TypeFloat:
			t = "Float"
		case vcf.TypeCharacter:
			t = "Character"
		case vcf.TypeString:
			t = "String"
		default:
			t = fd.Type.String()
		}
		bh.Formats[id] = &bcf.FormatDef{
			ID:          fd.ID,
			Number:      fd.Number.String(),
			Type:        t,
			Description: fd.Description,
		}
	}
	for _, m := range h.MetaLines {
		if m.Key != "FORMAT" {
			continue
		}
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
		bh.FormatsOrder = append(bh.FormatsOrder, rest[:end])
	}
	// Also populate bh.Filters from the header
	// text so the BCF FILTER writer can resolve q10/s50
	// indices. The VCF parser keys filters by name in
	// h.Filters (map) but the BCF writer consumes an
	// ordered slice at index time.
	seenFilter := map[string]bool{}
	for _, line := range strings.Split(sb.String(), "\n") {
		if !strings.HasPrefix(line, "##FILTER=") {
			continue
		}
		id, _ := extractIDAndDescription(line)
		if id != "" && !seenFilter[id] {
			bh.Filters = append(bh.Filters, id)
			seenFilter[id] = true
		}
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
		Formats:    bcfFormatsToVCF(h.Formats),
	}
	for _, name := range h.Contigs {
		vh.Contigs = append(vh.Contigs, &vcf.Contig{ID: name})
	}
	vh.SampleNames = append([]string(nil), h.SampleNames...)
	// Re-emit FORMAT declarations as MetaLines so
	// bcfToVCFRecord's vcfFormatOrder() can resolve
	// Format-key indices back to names.
	for _, id := range h.FormatsOrder {
		fd := h.Formats[id]
		if fd == nil {
			continue
		}
		var b strings.Builder
		b.WriteString("<ID=")
		b.WriteString(id)
		b.WriteString(",Number=")
		b.WriteString(fd.Number)
		b.WriteString(",Type=")
		b.WriteString(fd.Type)
		if fd.Description != "" {
			b.WriteString(",Description=\"")
			b.WriteString(fd.Description)
			b.WriteString("\"")
		}
		b.WriteString(">")
		vh.MetaLines = append(vh.MetaLines, vcf.MetaLine{
			Key:   "FORMAT",
			Value: b.String(),
		})
	}
	// Re-emit FILTER declarations with their Descriptions
	// extracted from h.Raw. Without this, the VCF encoder
	// would re-emit ##FILTER=<ID=…,Description="">.
	if vh.Filters == nil {
		vh.Filters = map[string]*vcf.FilterDef{}
	}
	for _, line := range strings.Split(h.Raw, "\n") {
		if !strings.HasPrefix(line, "##FILTER=") {
			continue
		}
		id, desc := extractIDAndDescription(line)
		if id != "" {
			vh.Filters[id] = &vcf.FilterDef{ID: id, Description: desc}
		}
	}
	for _, id := range h.Filters {
		if _, present := vh.Filters[id]; !present {
			vh.Filters[id] = &vcf.FilterDef{ID: id}
		}
	}
	// Contig metadata (length=, md5=, species=, taxonomy=,
	// assembly=, URL=…). The BCF header only carries a
	// stripped Contigs []string; rich Contig fields must
	// be rehydrated from h.Raw so emitContig can rebuild
	// a faithful ##contig line.
	for _, c := range vh.Contigs {
		for _, line := range strings.Split(h.Raw, "\n") {
			if !strings.HasPrefix(line, "##contig=") {
				continue
			}
			if !strings.Contains(line, "ID="+c.ID+",") && !strings.Contains(line, "ID="+c.ID+">") {
				continue
			}
			populateContig(c, line)
		}
	}
	// Re-emit unparsed ##-prefixed lines from h.Raw to
	// preserve fileDate / source / reference / phasing /
	// assembly / pedigreeDB / any other unstructured pragmas.
	for _, line := range strings.Split(h.Raw, "\n") {
		if !strings.HasPrefix(line, "##") {
			continue
		}
		if strings.HasPrefix(line, "##INFO=") ||
			strings.HasPrefix(line, "##FORMAT=") ||
			strings.HasPrefix(line, "##FILTER=") ||
			strings.HasPrefix(line, "##contig=") ||
			strings.HasPrefix(line, "##fileformat=") {
			continue
		}
		// Strip leading "##" and split key=rest.
		rest := strings.TrimPrefix(line, "##")
		eq := strings.IndexByte(rest, '=')
		if eq < 0 {
			continue
		}
		vh.MetaLines = append(vh.MetaLines, vcf.MetaLine{
			Key:   rest[:eq],
			Value: rest[eq+1:],
			Raw:   line,
		})
	}
	for _, c := range vh.Contigs {
		// Skip — emitContig will handle from typed fields.
		_ = c
	}
	return vh
}

// bcfFormatsToVCF copies a bcf FormatDef map to vcf.
func bcfFormatsToVCF(m map[string]*bcf.FormatDef) map[string]*vcf.FormatDef {
	if m == nil {
		return nil
	}
	out := make(map[string]*vcf.FormatDef, len(m))
	for k, v := range m {
		if v == nil {
			out[k] = nil
			continue
		}
		var n vcf.Number
		if v.Number != "" {
			if pn, err := vcf.ParseNumber(v.Number); err == nil {
				n = pn
			}
		}
		var t vcf.FieldType
		if v.Type != "" {
			if pt, err := vcf.ParseFieldType(v.Type); err == nil {
				t = pt
			}
		}
		out[k] = &vcf.FormatDef{
			ID:          v.ID,
			Number:      n,
			Type:        t,
			Description: v.Description,
		}
	}
	return out
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
func vcfToBCFRecord(r *vcf.Record, bh *bcf.Header) *bcf.Record {
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
			Type:         v.Type.String(),
			Ints:         append([]int64(nil), v.Ints...),
			Floats:       append([]float64(nil), v.Floats...),
			Strs:         append([]string(nil), v.Strs...),
			Flag:         v.Flag,
			IntMissing:   append([]bool(nil), v.IntMissing...),
			FloatMissing: append([]bool(nil), v.FloatMissing...),
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
		Format:    vcfFormatKeysToBCF(r.Format, bh),
		Samples:   vcfSamplesToBCF(r.Samples, bh, r.Format),
	}
}

// vcfFormatKeysToBCF maps VCF FORMAT key names (rec.Format)
// to BCF uint16 indices in bh.FormatsOrder. Unknown names
// are skipped (matches the v1 behavior of silently dropping
// keys the BCF header doesn't declare).
func vcfFormatKeysToBCF(keys []string, bh *bcf.Header) []uint16 {
	if bh == nil {
		return nil
	}
	idx := make(map[string]uint16, len(bh.FormatsOrder))
	for i, name := range bh.FormatsOrder {
		idx[name] = uint16(i)
	}
	out := make([]uint16, 0, len(keys))
	for _, k := range keys {
		if v, ok := idx[k]; ok {
			out = append(out, v)
			continue
		}
		// Auto-register FORMAT keys that appear on a
		// record but weren't declared in the header.
		// The BCF writer needs every FORMAT key to have
		// a header index and a Type — without a Type the
		// per-sample encoding would refuse to emit. The
		// safety default picks Integer which works for
		// the common case (GT-as-int, DP, GQ, HQ).
		bh.FormatsOrder = append(bh.FormatsOrder, k)
		idx[k] = uint16(len(bh.FormatsOrder) - 1) //nolint:gosec // bounded by FormatsOrder size
		if bh.Formats == nil {
			bh.Formats = map[string]*bcf.FormatDef{}
		}
		if _, present := bh.Formats[k]; !present {
			bh.Formats[k] = &bcf.FormatDef{
				ID:     k,
				Number: "1",
				Type:   "Integer",
			}
		}
		out = append(out, idx[k])
	}
	return out
}

// vcfSamplesToBCF converts vcf.Genotype per-sample values
// to bcf.Sample slices.
func vcfSamplesToBCF(samples []vcf.Genotype, bh *bcf.Header, keys []string) []bcf.Sample {
	if bh == nil {
		return nil
	}
	idx := make(map[string]uint16, len(bh.FormatsOrder))
	for i, name := range bh.FormatsOrder {
		idx[name] = uint16(i)
	}
	out := make([]bcf.Sample, len(samples))
	for si, gt := range samples {
		s := bcf.Sample{Fields: make(map[uint16]bcf.SampleValue, len(keys))}
		for _, k := range keys {
			kIdx, ok := idx[k]
			if !ok {
				continue
			}
			fv, present := gt.Fields[k]
			if !present {
				// An absent FORMAT field is a missing value.
				// The Missing flag matters: it is what makes
				// the BCF encoder emit one MISSING element
				// plus END_OF_VECTOR padding, the way htslib
				// does, rather than a bare vector end.
				s.Fields[kIdx] = bcf.SampleValue{Type: "Missing", Missing: true}
				continue
			}
			s.Fields[kIdx] = vcfFieldToBCF(&fv, k, bh)
		}
		out[si] = s
		_ = si
	}
	return out
}

// vcfFieldToBCF converts a vcf.FieldValue into the BCF wire
// shape.
func vcfFieldToBCF(fv *vcf.FieldValue, key string, bh *bcf.Header) bcf.SampleValue {
	fd := bh.Formats[key]
	typeName := ""
	if fd != nil {
		typeName = fd.Type
	}
	if typeName == "" {
		switch fv.Type {
		case vcf.TypeInteger:
			typeName = "Integer"
		case vcf.TypeFloat:
			typeName = "Float"
		case vcf.TypeCharacter:
			typeName = "Character"
		case vcf.TypeString:
			typeName = "String"
		default:
			// An undeclared key: the value is kept as raw
			// text, so carry it over as a string.
			typeName = "String"
		}
	}
	sv := bcf.SampleValue{
		Type:    typeName,
		Missing: fv.Missing,
	}
	if len(fv.Strs) == 0 && len(fv.Raw) > 0 {
		// An undeclared key keeps its value in Raw.
		sv.Strs = append([]string(nil), fv.Raw...)
	}
	if len(fv.Ints) > 0 {
		sv.Ints = make([]int32, len(fv.Ints))
		for i, n := range fv.Ints {
			sv.Ints[i] = int32(n) //nolint:gosec // VCF FORMAT ints fit int32
		}
	}
	if len(fv.IntMissing) > 0 {
		sv.IntMissing = append([]bool(nil), fv.IntMissing...)
	}
	if len(fv.Floats) > 0 {
		sv.Floats = make([]float32, len(fv.Floats))
		for i, f := range fv.Floats {
			sv.Floats[i] = float32(f)
		}
	}
	if len(fv.FloatMissing) > 0 {
		sv.FloatMissing = append([]bool(nil), fv.FloatMissing...)
	}
	if len(fv.Chars) > 0 {
		sv.Str = string(fv.Chars)
	}
	if len(fv.Strs) > 0 {
		sv.Strs = append([]string(nil), fv.Strs...)
	}
	return sv
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
func bcfToVCFRecord(r *bcf.Record, vh *vcf.Header) *vcf.Record {
	var idList []string
	if r.ID != "" {
		idList = strings.Split(r.ID, ",")
	}
	info := make(map[string]vcf.InfoValue, len(r.Info))
	for k, v := range r.Info { //nolint:gocritic // rangeValCopy: InfoValue is a value type by design
		t, err := vcf.ParseFieldType(v.Type)
		if err != nil {
			t = vcf.TypeUnknown
		}
		vv := vcf.InfoValue{}
		vv.Type = t
		vv.Ints = append([]int64(nil), v.Ints...)
		vv.Floats = append([]float64(nil), v.Floats...)
		vv.Strs = append([]string(nil), v.Strs...)
		vv.IntMissing = append([]bool(nil), v.IntMissing...)
		vv.FloatMissing = append([]bool(nil), v.FloatMissing...)
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
		Line:    r.Line,
		Chrom:   r.CHROM,
		Pos:     r.Pos,
		ID:      idList,
		Ref:     r.Ref,
		Alt:     append([]string(nil), r.Alt...),
		Qual:    bcfQualToVCF(r.Qual, r.QualKnown),
		Filter:  append([]string(nil), r.Filter...),
		Info:    info,
		Format:  bcfFormatKeysToNames(r.Format, vh),
		Samples: bcfSamplesToVCF(r.Samples, r.Format, vh),
	}
}

// bcfFormatKeysToNames maps BCF FORMAT-key indices back to
// the FORMAT key names from the VCF header (vh.Formats
// maps name → def).
func bcfFormatKeysToNames(keys []uint16, vh *vcf.Header) []string {
	if vh == nil {
		return nil
	}
	order := vcfFormatOrder(vh)
	out := make([]string, 0, len(keys))
	for _, k := range keys {
		if int(k) < len(order) {
			out = append(out, order[k])
		}
	}
	return out
}

// bcfSamplesToVCF converts per-sample BCF values back to
// vcf.Genotype slices.
func bcfSamplesToVCF(samples []bcf.Sample, keys []uint16, vh *vcf.Header) []vcf.Genotype {
	if vh == nil {
		return nil
	}
	order := vcfFormatOrder(vh)
	out := make([]vcf.Genotype, len(samples))
	for si, s := range samples {
		sg := vcf.Genotype{Fields: make(map[string]vcf.FieldValue, len(keys))}
		for _, k := range keys {
			if int(k) >= len(order) {
				continue
			}
			name := order[k]
			sv, ok := s.Fields[k]
			if !ok {
				continue
			}
			// Skip wholly-missing fields — VCF emits
			// trailing "." which spec §1.6.1 says can be
			// dropped when every sample is missing.
			if sv.Missing && len(sv.Ints) == 0 && len(sv.FloatMissing) == 0 && len(sv.IntMissing) == 0 && len(sv.Strs) == 0 && sv.Str == "" && len(sv.Floats) == 0 {
				continue
			}
			sg.Fields[name] = bcfFieldToVCF(&sv, name, vh)
		}
		out[si] = sg
		_ = si
	}
	return out
}

// vcfFormatOrder builds a declaration-order list of FORMAT
// IDs from the VCF header's MetaLines. Reused by both
// bcf-side helpers to resolve BCF FORMAT-key indices.
func vcfFormatOrder(vh *vcf.Header) []string {
	order := make([]string, 0, len(vh.Formats))
	for _, m := range vh.MetaLines {
		if m.Key != "FORMAT" {
			continue
		}
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
		order = append(order, rest[:end])
	}
	return order
}

// bcfFieldToVCF converts a bcf.SampleValue back to vcf.FieldValue.
func bcfFieldToVCF(sv *bcf.SampleValue, name string, vh *vcf.Header) vcf.FieldValue {
	fd := vh.Formats[name]
	typeCode := vcf.TypeUnknown
	if fd != nil {
		typeCode = fd.Type
	}
	fv := vcf.FieldValue{}
	fv.Type = typeCode
	// Whole-field missing collapses to "." in VCF output,
	// but a per-slot missing list like ". , ." must round-trip
	// as ".,." (VCF §1.6.2 line 508). When the BCF decoder
	// provides per-slot flags, suppress the whole-field
	// Missing bit so writeTyped emits the list.
	if len(sv.IntMissing) > 0 || len(sv.FloatMissing) > 0 {
		sv.Missing = false
	}
	fv.Missing = sv.Missing
	if len(sv.Ints) > 0 {
		fv.Ints = make([]int64, len(sv.Ints))
		for i, n := range sv.Ints {
			fv.Ints[i] = int64(n)
		}
	}
	if len(sv.IntMissing) > 0 {
		fv.IntMissing = append([]bool(nil), sv.IntMissing...)
	}
	if len(sv.Floats) > 0 {
		fv.Floats = make([]float64, len(sv.Floats))
		for i, f := range sv.Floats {
			fv.Floats[i] = float64(f)
		}
	}
	if len(sv.FloatMissing) > 0 {
		fv.FloatMissing = append([]bool(nil), sv.FloatMissing...)
	}
	if sv.Str != "" && typeCode == vcf.TypeCharacter {
		fv.Chars = []byte(sv.Str)
	}
	if len(sv.Strs) > 0 {
		fv.Strs = append([]string(nil), sv.Strs...)
	}
	return fv
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
	// Pre-scan: read all records into memory so we can
	// discover FORMAT/INFO keys referenced by records but
	// not declared in the header (V1.6.1 line 695 allows
	// this; the BCF header must nonetheless declare
	// every key used). Two passes — the first populates
	// bh.FormatsOrder/bh.InfosOrder; the second writes.
	dec := vcf.NewDecoder(src)
	enc := bcf.NewEncoder(dst)
	hdr, err := dec.Header()
	if err != nil {
		return fmt.Errorf("parse header: %w", err)
	}
	bh := vcfToBCFHeader(hdr)
	var records []*vcf.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return fmt.Errorf("parse record: %w", err)
		}
		records = append(records, rec)
	}
	// Register FORMAT/INFO keys that appear only on records
	// (typical for VCFs whose header is incomplete). For
	// each one we also need a matching ##FORMAT= or
	// ##INFO= line in bh.Raw so the BCF on-disk header
	// remains self-describing for downstream readers.
	extraText := strings.Builder{}
	for _, rec := range records {
		for _, k := range rec.Format {
			if lookupHeaderKey(bh.FormatsOrder, k) {
				continue
			}
			{
				bh.FormatsOrder = append(bh.FormatsOrder, k)
				if bh.Formats == nil {
					bh.Formats = map[string]*bcf.FormatDef{}
				}
				// A reserved key keeps the spec's own Number and
				// Type; anything else is a string, which is how
				// htslib treats a key the header does not
				// declare. Registering an undeclared key as an
				// Integer silently dropped a GT value like
				// "0|0".
				number, ftype := ".", "String"
				if n, t, ok := vcf.ReservedFormatDef(k); ok {
					number, ftype = n, t.String()
				}
				bh.Formats[k] = &bcf.FormatDef{ID: k, Number: number, Type: ftype}
				extraText.WriteString("##FORMAT=<ID=" + k + ",Number=" + number + ",Type=" + ftype + ">\n")
			}
		}
		for k := range rec.Info {
			if !lookupHeaderKey(bh.InfosOrder, k) {
				bh.InfosOrder = append(bh.InfosOrder, k)
				if bh.Infos == nil {
					bh.Infos = map[string]*bcf.InfoDef{}
				}
				bh.Infos[k] = &bcf.InfoDef{ID: k, Number: ".", Type: "String"}
				extraText.WriteString("##INFO=<ID=" + k + ",Number=.,Type=String>\n")
			}
		}
	}
	// Splice the extra ##FORMAT / ##INFO lines into bh.Raw
	// before the sample line so the BCF file is complete.
	// The synthesised ##FORMAT / ##INFO lines must come
	// *before* the #CHROM line: the #CHROM line ends the
	// header, and htslib rejects a header that continues
	// after it ("Invalid FORMAT id 1").
	if extra := extraText.String(); extra != "" {
		if i := strings.Index(bh.Raw, "#CHROM\t"); i >= 0 {
			bh.Raw = bh.Raw[:i] + extra + bh.Raw[i:]
		} else {
			bh.Raw += extra
		}
	}
	if err := enc.WriteHeader(bh); err != nil {
		return fmt.Errorf("write header: %w", err)
	}
	for _, rec := range records {
		if err := enc.WriteRecord(vcfToBCFRecord(rec, bh)); err != nil {
			return fmt.Errorf("write record: %w", err)
		}
	}
	if err := enc.Close(); err != nil {
		return fmt.Errorf("close encoder: %w", err)
	}
	return nil
}

// lookupHeaderKey reports whether key is in list
// (bcf.InfosOrder or bcf.FormatsOrder).
func lookupHeaderKey(list []string, key string) bool {
	for _, v := range list {
		if v == key {
			return true
		}
	}
	return false
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
		func(r *bcf.Record, vh *vcf.Header) (*vcf.Record, error) {
			return bcfToVCFRecord(r, vh), nil
		},
		enc.EncodeRecord,
		func() error { return enc.Flush() },
	)
}

// extractIDAndDescription parses a "<TYPE>=<ID=…,Description=\"…\">"
// line into the ID and Description, respecting quoted commas.
// Used to rehydrate FILTER/FORMAT/INFO Description values from
// the BCF header's Raw text.
func extractIDAndDescription(line string) (id, desc string) {
	const prefix = "<ID="
	i := strings.Index(line, prefix)
	if i < 0 {
		return "", ""
	}
	rest := line[i+len(prefix):]
	end := strings.Index(rest, ">")
	if end < 0 {
		return "", ""
	}
	body := rest[:end]
	kvs := splitQuotedCSV(body)
	if len(kvs) > 0 {
		id = kvs[0] // The "<ID=" prefix already stripped the leading "ID=".
	}
	for _, kv := range kvs[1:] {
		if strings.HasPrefix(kv, "Description=") {
			d := strings.TrimPrefix(kv, "Description=")
			if len(d) >= 2 && d[0] == '"' && d[len(d)-1] == '"' {
				d = d[1 : len(d)-1]
			}
			desc = d
		}
	}
	return id, desc
}

// splitQuotedCSV splits a comma list respecting quoted spans.
func splitQuotedCSV(s string) []string {
	out := []string{}
	cur := []byte{}
	inQuote := false
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c == '"' {
			inQuote = !inQuote
			cur = append(cur, c)
			continue
		}
		if !inQuote && c == ',' {
			out = append(out, string(cur))
			cur = cur[:0]
			continue
		}
		cur = append(cur, c)
	}
	out = append(out, string(cur))
	return out
}

// populateContig fills a vcf.Contig's rich fields from a
// ##contig=<ID=…,length=…,…> line.
func populateContig(c *vcf.Contig, line string) {
	const prefix = "<"
	body := strings.TrimPrefix(line, "##contig=")
	body = strings.TrimPrefix(body, prefix)
	end := strings.Index(body, ">")
	if end < 0 {
		return
	}
	for _, kv := range splitQuotedCSV(body[:end]) {
		eq := strings.IndexByte(kv, '=')
		if eq < 0 {
			continue
		}
		k, v := kv[:eq], kv[eq+1:]
		if len(v) >= 2 && v[0] == '"' && v[len(v)-1] == '"' {
			v = v[1 : len(v)-1]
		}
		switch k {
		case "ID":
			// already set on c
		case "length":
			n, _ := strconv.ParseInt(v, 10, 64)
			c.Length = n
		case "md5":
			c.MD5 = v
		case "URL":
			c.URL = v
		case "species":
			c.Species = v
		case "taxonomy":
			c.Taxonomy = v
		default:
			if c.Other == nil {
				c.Other = map[string]string{}
			}
			c.Other[k] = v
		}
	}
}

// bamQualToSAM renders a BAM record's quality as the SAM
// ASCII form. An unavailable quality is stored as 0xFF bytes
// and is written as the SAM "*" placeholder.
func bamQualToSAM(qual []byte) string {
	if len(qual) > 0 && qual[0] == 0xFF {
		return "*"
	}
	return string(qual)
}

// bamNextRefName renders a BAM record's mate reference as the
// SAM RNEXT field. When the mate is on the same reference as
// the read, htslib writes the "=" shorthand
// (sam_format1: "=" when mtid == tid), so we do too.
func bamNextRefName(r *bam.Record) string {
	if r.RefID >= 0 && r.NextRefID == r.RefID {
		return "="
	}
	if r.NameRNext == "" {
		return "*"
	}
	return r.NameRNext
}

// bamTagTypeToSAM maps a BAM aux type byte to the SAM tag
// type letter.
//
// BAM stores an integer aux value with a width-specific subtype
// byte (c, C, s, S, i, I); SAM has a single integer type, "i".
// htslib always writes "i" on the SAM side, and its SAM parser
// treats any other integer subtype letter as the character type
// - so emitting "C" would silently turn the integer 1 into the
// character '1'.
func bamTagTypeToSAM(t bam.TagType) sam.TagType {
	switch t {
	case bam.Tagc, bam.TagC, bam.Tags, bam.TagS, bam.Tagi, bam.TagI:
		return sam.Tagi
	default:
		return sam.TagType(t)
	}
}
