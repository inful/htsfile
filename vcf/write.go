package vcf

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"
)

// Encoder writes a VCF document to an underlying [io.Writer].
//
// Typical usage:
//
//	enc := vcf.NewEncoder(w)
//	enc.EncodeHeader(hdr)
//	for _, rec := range records { enc.EncodeRecord(rec) }
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w *bufio.Writer
	// lineEnding is the byte sequence used to terminate each emitted line.
	lineEnding []byte
}

// NewEncoder wraps w in an [Encoder] with default LF line endings.
func NewEncoder(w io.Writer) *Encoder {
	bw, ok := w.(*bufio.Writer)
	if !ok {
		bw = bufio.NewWriterSize(w, 64*1024)
	}
	return &Encoder{w: bw, lineEnding: []byte("\n")}
}

// NewEncoderCRLF returns an Encoder that terminates every emitted line with
// CR+LF instead of LF. Both forms are legal per VCFv4.5 §1.
func NewEncoderCRLF(w io.Writer) *Encoder {
	e := NewEncoder(w)
	e.lineEnding = []byte("\r\n")
	return e
}

// EncodeHeader emits the VCF header: ##fileformat, all meta-information
// lines, and the #CHROM header line.
//
// Structured meta lines are emitted from the typed header fields when those
// fields are populated; the encoder deduplicates against [Header.MetaLines]
// so that round-tripping a parsed-and-written file produces a single
// occurrence of every header line.
//
// If h.FileFormat is the zero value, EncodeHeader errors out — a VCF file
// must declare its version on the first line.
func (e *Encoder) EncodeHeader(h *Header) error {
	if h.FileFormat.Major == 0 && h.FileFormat.Minor == 0 {
		return errors.New("vcf: EncodeHeader: Header.FileFormat is unset")
	}
	// Build a quick lookup of MetaLine keys to detect duplicates.
	metaSeen := map[string]bool{}

	// 1. ##fileformat must be the first line.
	if _, err := e.w.WriteString("##fileformat="); err != nil {
		return err
	}
	if _, err := e.w.WriteString(h.FileFormat.String()); err != nil {
		return err
	}
	if _, err := e.w.Write(e.lineEnding); err != nil {
		return err
	}

	// 2. Walk MetaLines in source order. Skip entries the typed model
	// will produce.
	for _, m := range h.MetaLines {
		if m.Key == "fileformat" {
			continue
		}
		if !e.shouldEmitFromMeta(m, h) {
			continue
		}
		if err := e.writeRawMetaLine(m.Key, m.Value); err != nil {
			return err
		}
		metaSeen[m.Key] = true
	}

	// 3. Emit common unstructured typed fields if non-empty and not
	// already emitted.
	if h.FileDate != "" && !metaSeen["fileDate"] {
		if err := e.writeRawMetaLine("fileDate", h.FileDate); err != nil {
			return err
		}
	}
	if h.Source != "" && !metaSeen["source"] {
		if err := e.writeRawMetaLine("source", h.Source); err != nil {
			return err
		}
	}
	if h.Reference != "" && !metaSeen["reference"] {
		if err := e.writeRawMetaLine("reference", h.Reference); err != nil {
			return err
		}
	}
	if h.Assembly != "" && !metaSeen["assembly"] {
		if err := e.writeRawMetaLine("assembly", h.Assembly); err != nil {
			return err
		}
	}
	if h.PedigreeDB != "" && !metaSeen["pedigreeDB"] {
		if err := e.writeRawMetaLine("pedigreeDB", h.PedigreeDB); err != nil {
			return err
		}
	}

	// 4. ##contig — preserve insertion order.
	for _, c := range h.Contigs {
		if err := e.emitContig(c); err != nil {
			return err
		}
	}

	// 5. ##INFO, ##FILTER, ##FORMAT, ##ALT — sort by ID for deterministic
	// output.
	for _, id := range sortedKeys(h.Infos) {
		if err := e.emitInfo(h.Infos[id]); err != nil {
			return err
		}
	}
	for _, id := range sortedKeys(h.Filters) {
		if err := e.emitFilter(h.Filters[id]); err != nil {
			return err
		}
	}
	for _, id := range sortedKeys(h.Formats) {
		if err := e.emitFormat(h.Formats[id]); err != nil {
			return err
		}
	}
	for _, id := range sortedKeys(h.Alts) {
		if err := e.emitAlt(h.Alts[id]); err != nil {
			return err
		}
	}

	// 6. ##PEDIGREE, ##SAMPLE — preserve order.
	for _, p := range h.Pedigrees {
		if err := e.emitPedigree(p); err != nil {
			return err
		}
	}
	for _, s := range h.Samples {
		if err := e.emitSample(s); err != nil {
			return err
		}
	}

	// 7. #CHROM header line.
	if err := e.emitHeaderLine(h); err != nil {
		return err
	}
	return e.w.Flush()
}

// shouldEmitFromMeta decides whether a MetaLine should be emitted verbatim.
// Lines for typed-model keys (contig/INFO/FORMAT/FILTER/ALT/PEDIGREE/SAMPLE)
// are always skipped here because the typed fields produce them. Common
// short keys (fileDate/source/reference/assembly/pedigreeDB) are skipped if
// the corresponding typed field is populated.
func (e *Encoder) shouldEmitFromMeta(m MetaLine, h *Header) bool {
	switch m.Key {
	case "contig", "INFO", "FILTER", "FORMAT", "ALT", "PEDIGREE", "SAMPLE":
		return false
	case "fileDate":
		return h.FileDate == ""
	case "source":
		return h.Source == ""
	case "reference":
		return h.Reference == ""
	case "assembly":
		return h.Assembly == ""
	case "pedigreeDB":
		return h.PedigreeDB == ""
	}
	return true
}

// writeRawMetaLine writes "##key=value\n" verbatim.
func (e *Encoder) writeRawMetaLine(key, value string) error {
	if _, err := e.w.WriteString("##"); err != nil {
		return err
	}
	if _, err := e.w.WriteString(key); err != nil {
		return err
	}
	if _, err := e.w.WriteString("="); err != nil {
		return err
	}
	if _, err := e.w.WriteString(value); err != nil {
		return err
	}
	_, err := e.w.Write(e.lineEnding)
	return err
}

// emitContig writes a ##contig=<ID=...,length=...,...> line.
func (e *Encoder) emitContig(c *Contig) error {
	var pairs [][2]string
	pairs = append(pairs, [2]string{"ID", quoteIfNeeded(c.ID)})
	if c.Length > 0 {
		pairs = append(pairs, [2]string{"length", strconv.FormatInt(c.Length, 10)})
	}
	if c.MD5 != "" {
		pairs = append(pairs, [2]string{"md5", strconv.Quote(c.MD5)})
	}
	if c.URL != "" {
		pairs = append(pairs, [2]string{"URL", strconv.Quote(c.URL)})
	}
	for k, v := range c.Other {
		pairs = append(pairs, [2]string{k, strconv.Quote(v)})
	}
	return e.writeStructured("contig", pairs)
}

// emitInfo writes a ##INFO=<Id=...,Number=...,Type=...,Description=...,...> line.
func (e *Encoder) emitInfo(d *InfoDef) error {
	pairs := make([][2]string, 0, 4+len(d.Other))
	pairs = append(pairs,
		[2]string{"ID", quoteIfNeeded(d.ID)},
		[2]string{"Number", d.Number.String()},
		[2]string{"Type", d.Type.String()},
		[2]string{"Description", strconv.Quote(d.Description)},
	)
	if d.Source != "" {
		pairs = append(pairs, [2]string{"Source", strconv.Quote(d.Source)})
	}
	if d.Version != "" {
		pairs = append(pairs, [2]string{"Version", strconv.Quote(d.Version)})
	}
	for k, v := range d.Other {
		pairs = append(pairs, [2]string{k, strconv.Quote(v)})
	}
	return e.writeStructured("INFO", pairs)
}

// emitFilter writes a ##FILTER=<ID=...,Description=...> line.
func (e *Encoder) emitFilter(d *FilterDef) error {
	pairs := make([][2]string, 0, 2+len(d.Other))
	pairs = append(pairs,
		[2]string{"ID", quoteIfNeeded(d.ID)},
		[2]string{"Description", strconv.Quote(d.Description)},
	)
	for k, v := range d.Other {
		pairs = append(pairs, [2]string{k, strconv.Quote(v)})
	}
	return e.writeStructured("FILTER", pairs)
}

// emitFormat writes a ##FORMAT=<ID=...,Number=...,Type=...,Description=...> line.
func (e *Encoder) emitFormat(d *FormatDef) error {
	pairs := make([][2]string, 0, 4+len(d.Other))
	pairs = append(pairs,
		[2]string{"ID", quoteIfNeeded(d.ID)},
		[2]string{"Number", d.Number.String()},
		[2]string{"Type", d.Type.String()},
		[2]string{"Description", strconv.Quote(d.Description)},
	)
	for k, v := range d.Other {
		pairs = append(pairs, [2]string{k, strconv.Quote(v)})
	}
	return e.writeStructured("FORMAT", pairs)
}

// emitAlt writes a ##ALT=<ID=...,Description=...> line.
func (e *Encoder) emitAlt(d *AltDef) error {
	pairs := make([][2]string, 0, 2+len(d.Other))
	pairs = append(pairs,
		[2]string{"ID", quoteIfNeeded(d.ID)},
		[2]string{"Description", strconv.Quote(d.Description)},
	)
	for k, v := range d.Other {
		pairs = append(pairs, [2]string{k, strconv.Quote(v)})
	}
	return e.writeStructured("ALT", pairs)
}

// emitPedigree writes a ##PEDIGREE=<ID=...,...> line.
func (e *Encoder) emitPedigree(p *Pedigree) error {
	pairs := [][2]string{
		{"ID", quoteIfNeeded(p.ID)},
	}
	if p.Original != "" {
		pairs = append(pairs, [2]string{"Original", quoteIfNeeded(p.Original)})
	}
	if p.Father != "" {
		pairs = append(pairs, [2]string{"Father", quoteIfNeeded(p.Father)})
	}
	if p.Mother != "" {
		pairs = append(pairs, [2]string{"Mother", quoteIfNeeded(p.Mother)})
	}
	for i, name := range p.Names {
		pairs = append(pairs, [2]string{
			fmt.Sprintf("Name_%d", i+1), quoteIfNeeded(name),
		})
	}
	for k, v := range p.Other {
		pairs = append(pairs, [2]string{k, strconv.Quote(v)})
	}
	return e.writeStructured("PEDIGREE", pairs)
}

// emitSample writes a ##SAMPLE=<ID=...,...> line.
func (e *Encoder) emitSample(s *Sample) error {
	pairs := [][2]string{
		{"ID", quoteIfNeeded(s.ID)},
	}
	if len(s.Genomes) > 0 {
		inner := make([][2]string, 0, len(s.Genomes))
		for k, v := range s.Genomes {
			inner = append(inner, [2]string{k, strconv.Quote(v)})
		}
		sort.SliceStable(inner, func(i, j int) bool { return inner[i][0] < inner[j][0] })
		pairs = append(pairs, [2]string{"Genomes", buildStructuredBody(inner)})
	}
	if len(s.MolConvs) > 0 {
		inner := make([][2]string, 0, len(s.MolConvs))
		for k, v := range s.MolConvs {
			inner = append(inner, [2]string{k, strconv.Quote(v)})
		}
		sort.SliceStable(inner, func(i, j int) bool { return inner[i][0] < inner[j][0] })
		pairs = append(pairs, [2]string{"MolConvs", buildStructuredBody(inner)})
	}
	if s.Description != "" {
		pairs = append(pairs, [2]string{"Description", strconv.Quote(s.Description)})
	}
	if s.DOI != "" {
		pairs = append(pairs, [2]string{"DOI", strconv.Quote(s.DOI)})
	}
	for k, v := range s.Other {
		pairs = append(pairs, [2]string{k, strconv.Quote(v)})
	}
	return e.writeStructured("SAMPLE", pairs)
}

// writeStructured emits "##key=<k1=v1,k2=v2,...>\n".
func (e *Encoder) writeStructured(key string, pairs [][2]string) error {
	if _, err := e.w.WriteString("##"); err != nil {
		return err
	}
	if _, err := e.w.WriteString(key); err != nil {
		return err
	}
	if _, err := e.w.WriteString("=<"); err != nil {
		return err
	}
	if _, err := e.w.WriteString(buildStructuredBody(pairs)); err != nil {
		return err
	}
	if _, err := e.w.WriteString(">"); err != nil {
		return err
	}
	_, err := e.w.Write(e.lineEnding)
	return err
}

// buildStructuredBody joins [k=v, k=v, ...] into "k1=v1,k2=v2,..." text,
// quoting string values via strconv.Quote so that backslashes and inner
// quotes survive the round-trip.
func buildStructuredBody(pairs [][2]string) string {
	var b strings.Builder
	for i, p := range pairs {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(p[0])
		b.WriteByte('=')
		b.WriteString(p[1])
	}
	return b.String()
}

// quoteIfNeeded quotes the value only when it contains characters that
// would otherwise confuse the structured-body parser (commas, equals, quotes,
// newlines, tabs). Most ID values are bare.
func quoteIfNeeded(s string) string {
	if strings.ContainsAny(s, ",=\"\r\n\t") {
		return strconv.Quote(s)
	}
	return s
}

// emitHeaderLine writes the final "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO[\\t<format>\\t<sample>s]" line.
func (e *Encoder) emitHeaderLine(h *Header) error {
	if _, err := e.w.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"); err != nil {
		return err
	}
	if len(h.SampleNames) > 0 {
		if _, err := e.w.WriteString("\tFORMAT"); err != nil {
			return err
		}
		for _, name := range h.SampleNames {
			if _, err := e.w.WriteString("\t"); err != nil {
				return err
			}
			if _, err := e.w.WriteString(name); err != nil {
				return err
			}
		}
	}
	_, err := e.w.Write(e.lineEnding)
	return err
}

// EncodeRecord writes a single VCF data line for rec. The header is not
// required to have been emitted first; that is the caller's responsibility.
func (e *Encoder) EncodeRecord(rec *Record) error {
	if rec == nil {
		return errors.New("vcf: EncodeRecord: nil record")
	}
	var b bytes.Buffer

	encodeFixedColumns(&b, rec)

	// INFO
	if len(rec.Info) == 0 {
		b.WriteByte('.')
	} else {
		// Sort INFO keys for deterministic output.
		keys := make([]string, 0, len(rec.Info))
		for k := range rec.Info {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		writeInfoColumn(&b, keys, rec.Info)
	}

	// FORMAT + sample columns
	if len(rec.Format) == 0 && len(rec.Samples) == 0 {
		_, err := e.w.Write(b.Bytes())
		if err != nil {
			return err
		}
		_, err = e.w.Write(e.lineEnding)
		return err
	}
	b.WriteByte('\t')
	b.WriteString(strings.Join(rec.Format, ":"))

	for _, s := range rec.Samples {
		b.WriteByte('\t')
		encodeSample(&b, rec.Format, s)
	}
	if _, err := e.w.Write(b.Bytes()); err != nil {
		return err
	}
	_, err := e.w.Write(e.lineEnding)
	return err
}

// encodeFixedColumns writes CHROM..FILTER into b. Extracted from
// EncodeRecord to keep the cognitive complexity manageable.
func encodeFixedColumns(b *bytes.Buffer, rec *Record) {
	// CHROM
	b.WriteString(rec.Chrom)
	b.WriteByte('\t')
	// POS
	b.WriteString(strconv.FormatInt(rec.Pos, 10))
	b.WriteByte('\t')
	// ID
	if len(rec.ID) == 0 {
		b.WriteByte('.')
	} else {
		b.WriteString(strings.Join(rec.ID, ";"))
	}
	b.WriteByte('\t')
	// REF — never percent-encode; only ACGTN bases (and the padding base
	// for symbolic SVs) are allowed here.
	b.WriteString(rec.Ref)
	b.WriteByte('\t')
	// ALT — never percent-encode. The spec example for breakends uses
	// literal ':' in position specifiers, and ',' is the only character
	// with structural meaning here. Values inside <...> symbolic alleles
	// are likewise raw.
	if len(rec.Alt) == 0 {
		b.WriteByte('.')
	} else {
		for i, a := range rec.Alt {
			if i > 0 {
				b.WriteByte(',')
			}
			b.WriteString(a)
		}
	}
	b.WriteByte('\t')
	// QUAL
	if rec.Qual == nil {
		b.WriteByte('.')
	} else {
		b.WriteString(formatVCFFloat(*rec.Qual))
	}
	b.WriteByte('\t')
	// FILTER
	if len(rec.Filter) == 0 {
		b.WriteByte('.')
	} else {
		b.WriteString(strings.Join(rec.Filter, ";"))
	}
	b.WriteByte('\t')
}

// writeInfoColumn writes the INFO column contents to b. keys should already
// be in the desired order.
func writeInfoColumn(b *bytes.Buffer, keys []string, info map[string]InfoValue) {
	for i, k := range keys {
		if i > 0 {
			b.WriteByte(';')
		}
		v := info[k]
		b.WriteString(k)
		if isFlagInfo(&v) {
			continue
		}
		b.WriteByte('=')
		if v.Missing {
			b.WriteByte('.')
			continue
		}
		switch v.Type {
		case TypeInteger:
			for j, n := range v.Ints {
				if j > 0 {
					b.WriteByte(',')
				}
				b.WriteString(strconv.FormatInt(n, 10))
			}
		case TypeFloat:
			for j, f := range v.Floats {
				if j > 0 {
					b.WriteByte(',')
				}
				b.WriteString(formatVCFFloat(f))
			}
		case TypeString:
			for j, s := range v.Strs {
				if j > 0 {
					b.WriteByte(',')
				}
				b.WriteString(percentEncode(s))
			}
		case TypeCharacter:
			if len(v.Chars) > 0 {
				b.WriteByte(v.Chars[0])
			}
		default:
			// Unknown type — emit raw if available.
			for j, s := range v.Raw {
				if j > 0 {
					b.WriteByte(',')
				}
				b.WriteString(s)
			}
		}
	}
}

// isFlagInfo reports whether v should be emitted as a key-only Flag token
// (no "=value").
func isFlagInfo(v *InfoValue) bool {
	if v.Type == TypeFlag {
		return v.Flag
	}
	// TypeUnknown may also be present in head; if it has no raw values we
	// can still treat it as flag-like (e.g. "DB" without "=").
	return len(v.Raw) == 0 && v.Type == TypeUnknown
}

// encodeSample emits one sample column for a record. It honours the
// spec §1.6.1 "trailing fields can be dropped" rule: keys past the last
// present FORMAT key are not emitted at all (the absence is implicit).
func encodeSample(b *bytes.Buffer, format []string, s Genotype) {
	// Find the index of the last FORMAT key that has a present value in
	// this sample. Keys beyond it are dropped and must not be emitted.
	lastPresent := -1
	for i, key := range format {
		if _, ok := s.Fields[key]; ok {
			lastPresent = i
		}
	}
	for i, key := range format {
		if i > lastPresent {
			break
		}
		if i > 0 {
			b.WriteByte(':')
		}
		v := s.Fields[key]
		if v.Missing {
			b.WriteByte('.')
			continue
		}
		switch v.Type {
		case TypeInteger:
			for j, n := range v.Ints {
				if j > 0 {
					b.WriteByte(',')
				}
				b.WriteString(strconv.FormatInt(n, 10))
			}
		case TypeFloat:
			for j, f := range v.Floats {
				if j > 0 {
					b.WriteByte(',')
				}
				b.WriteString(formatVCFFloat(f))
			}
		case TypeString:
			for j, s := range v.Strs {
				if j > 0 {
					b.WriteByte(',')
				}
				b.WriteString(percentEncode(s))
			}
		case TypeCharacter:
			if len(v.Chars) > 0 {
				b.WriteByte(v.Chars[0])
			}
		default:
			// Unknown or missing decoded type — fall back to raw.
			if len(v.Raw) > 0 {
				for j, s := range v.Raw {
					if j > 0 {
						b.WriteByte(',')
					}
					b.WriteString(s)
				}
			} else {
				b.WriteByte('.')
			}
		}
	}
}

// formatVCFFloat renders f as the spec permits: "NaN", "Infinity",
// "-Infinity" (case insensitive on parse) or the shortest round-trippable
// decimal.
//
// We use the standard %g format which produces minimal output for
// round-trip purposes.
func formatVCFFloat(f float64) string {
	if math.IsNaN(f) {
		return "NaN"
	}
	if math.IsInf(f, +1) {
		return "Infinity"
	}
	if math.IsInf(f, -1) {
		return "-Infinity"
	}
	// %g uses %e for large exponents or %f otherwise; it is the shortest
	// representation that round-trips losslessly.
	return strconv.FormatFloat(f, 'g', -1, 32)
}

// sortedKeys returns the keys of m sorted ascending. Used for deterministic
// header output.
func sortedKeys[V any](m map[string]V) []string {
	out := make([]string, 0, len(m))
	for k := range m {
		out = append(out, k)
	}
	sort.Strings(out)
	return out
}

// Flush ensures any buffered data is written to the underlying [io.Writer].
// Encoder buffers internally; calling Flush is only required when the
// caller has not yet finished emitting all records.
func (e *Encoder) Flush() error { return e.w.Flush() }
