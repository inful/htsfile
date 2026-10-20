package bam

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/inful/htsfile/internal/bgzf"
)

// BAM magic: the literal ASCII "BAM" followed by 0x01.
var bamMagic = [4]byte{'B', 'A', 'M', 0x01}

// Decoder reads a BAM document incrementally from an underlying
// BGZF-compressed [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("sample.bam")
//	defer f.Close()
//	dec := bam.NewDecoder(f)
//	hdr, _ := dec.Header()    // required; returns error on malformed header
//	for {
//	    rec, err := dec.Record()
//	    if errors.Is(err, io.EOF) { break }
//	    if err != nil { ... }
//	    // rec.QName, rec.Pos, rec.Flag, rec.Cigar, ...
//	}
//
// A Decoder is not safe for concurrent use; create one per
// goroutine.
type Decoder struct {
	// br reads from a bgzf-decoded stream with a small buffer
	// to absorb the frequent short reads inside the
	// alignment-record loop. The bgzf.Reader is wrapped once;
	// a bufio.Reader sits on top of it.
	br *bufio.Reader

	// maxLineLength is preserved to keep the constructor surface
	// symmetric with the [sam] package, but is unused: BAM is
	// framed at the BGZF block level (max 64 KiB each), not by
	// lines, so no line-length cap is meaningful here.
	maxLineLength int

	// Header state.
	hdrRead bool
	hdr     *Header
	hdrErr  error

	// recordLine is the 1-based index of the next alignment
	// record to read.
	recordLine int
}

// DefaultMaxLineLength is the per-line size cap applied by
// [NewDecoder]. For BAM it has no behavioural effect (BAM
// records are framed by BGZF blocks, not lines) but is kept for
// API symmetry with the [sam] and [vcf] packages.
const DefaultMaxLineLength = 10 << 20 // 10 MB

// NewDecoder wraps r in a [Decoder]. r is expected to be a
// BGZF-compressed BAM stream.
func NewDecoder(r io.Reader) *Decoder {
	return NewDecoderWithMaxLine(r, DefaultMaxLineLength)
}

// NewDecoderWithMaxLine is like [NewDecoder] but lets the caller
// specify a max-line cap. For BAM the cap is unused (records are
// BGZF-block framed), but the constructor is provided for API
// symmetry with the [sam] and [vcf] packages.
func NewDecoderWithMaxLine(r io.Reader, maxLineLength int) *Decoder {
	if maxLineLength <= 0 {
		maxLineLength = DefaultMaxLineLength
	}
	gz := bgzf.NewReader(r)
	return &Decoder{
		br:            bufio.NewReader(gz),
		maxLineLength: maxLineLength,
	}
}

// Header returns the file's header. The first call reads the
// magic, embedded SAM text, and reference-sequence dictionary.
// Subsequent calls return the same value.
//
// A missing or malformed magic returns a *ParseError.
func (d *Decoder) Header() (*Header, error) {
	if d.hdrRead {
		return d.hdr, d.hdrErr
	}
	d.hdrRead = true
	// Read the 4-byte magic.
	var magic [4]byte
	if _, err := io.ReadFull(d.br, magic[:]); err != nil {
		if errors.Is(err, io.EOF) || errors.Is(err, io.ErrUnexpectedEOF) {
			d.hdrErr = &ParseError{Field: "magic", Msg: "BAM stream too short for magic"}
			return nil, d.hdrErr
		}
		d.hdrErr = err
		return nil, err
	}
	if magic != bamMagic {
		d.hdrErr = &ParseError{Field: "magic", Msg: fmt.Sprintf("not a BAM file: got %q", magic[:])}
		return nil, d.hdrErr
	}
	// Read l_text + text.
	var lText uint32
	if err := binary.Read(d.br, binary.LittleEndian, &lText); err != nil {
		d.hdrErr = &ParseError{Field: "header", Msg: "failed to read l_text: " + err.Error()}
		return nil, d.hdrErr
	}
	textBytes := make([]byte, lText)
	if _, err := io.ReadFull(d.br, textBytes); err != nil {
		d.hdrErr = &ParseError{Field: "header", Msg: "failed to read header text: " + err.Error()}
		return nil, d.hdrErr
	}
	// Read n_ref + reference dictionary.
	var nRef uint32
	if err := binary.Read(d.br, binary.LittleEndian, &nRef); err != nil {
		d.hdrErr = &ParseError{Field: "header", Msg: "failed to read n_ref: " + err.Error()}
		return nil, d.hdrErr
	}
	sqs := make([]*SQ, 0, nRef)
	for i := uint32(0); i < nRef; i++ {
		var lName uint32
		if err := binary.Read(d.br, binary.LittleEndian, &lName); err != nil {
			d.hdrErr = &ParseError{Field: "ref", Msg: fmt.Sprintf("ref %d: failed to read l_name: %s", i, err.Error())}
			d.hdr = &Header{Raw: string(textBytes), SQs: sqs}
			return d.hdr, d.hdrErr
		}
		nameBytes := make([]byte, lName)
		if _, err := io.ReadFull(d.br, nameBytes); err != nil {
			d.hdrErr = &ParseError{Field: "ref", Msg: fmt.Sprintf("ref %d: failed to read name: %s", i, err.Error())}
			d.hdr = &Header{Raw: string(textBytes), SQs: sqs}
			return d.hdr, d.hdrErr
		}
		// Name is null-terminated; strip the trailing NUL.
		name := string(nameBytes)
		if name != "" && name[len(name)-1] == 0 {
			name = name[:len(name)-1]
		}
		var lRef uint32
		if err := binary.Read(d.br, binary.LittleEndian, &lRef); err != nil {
			d.hdrErr = &ParseError{Field: "ref", Msg: fmt.Sprintf("ref %d: failed to read l_ref: %s", i, err.Error())}
			d.hdr = &Header{Raw: string(textBytes), SQs: sqs}
			return d.hdr, d.hdrErr
		}
		sqs = append(sqs, &SQ{Name: name, Len: int64(lRef)})
	}
	// Parse the embedded SAM header text.
	hdr := &Header{
		Raw: string(textBytes),
		SQs: sqs,
	}
	parseHeaderText(hdr)
	d.hdr = hdr
	return d.hdr, nil
}

// Record returns the next alignment record. When no more records
// are available it returns (nil, [io.EOF]). Calling Record after
// EOF continues to return (nil, [io.EOF]).
//
// If Header has not yet been called, Record calls it implicitly.
func (d *Decoder) Record() (*Record, error) {
	if !d.hdrRead {
		if _, err := d.Header(); err != nil {
			return nil, err
		}
	}
	// Read block_size (4 bytes).
	var blockSize uint32
	if err := binary.Read(d.br, binary.LittleEndian, &blockSize); err != nil {
		if errors.Is(err, io.EOF) || errors.Is(err, io.ErrUnexpectedEOF) {
			return nil, io.EOF
		}
		return nil, &ParseError{Field: "record", Msg: "failed to read block_size: " + err.Error()}
	}
	d.recordLine++
	// Read the fixed-size record core (32 bytes after block_size).
	const coreSize = 32
	if blockSize < coreSize {
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: fmt.Sprintf("block_size %d < core size %d", blockSize, coreSize)}
	}
	core := make([]byte, coreSize)
	if _, err := io.ReadFull(d.br, core); err != nil {
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: "failed to read record core: " + err.Error()}
	}
	refID := int32(binary.LittleEndian.Uint32(core[0:4]))
	pos0 := int32(binary.LittleEndian.Uint32(core[4:8]))
	lReadName := core[8]
	mapq := core[9]
	bin := binary.LittleEndian.Uint16(core[10:12])
	nCigarOp := binary.LittleEndian.Uint16(core[12:14])
	flag := binary.LittleEndian.Uint16(core[14:16])
	lSeq := int32(binary.LittleEndian.Uint32(core[16:20]))
	if lSeq < 0 {
		return nil, &ParseError{Line: d.recordLine, Field: "l_seq",
			Msg: fmt.Sprintf("negative l_seq %d", lSeq)}
	}
	nextRefID := int32(binary.LittleEndian.Uint32(core[20:24]))
	nextPos0 := int32(binary.LittleEndian.Uint32(core[24:28]))
	tlen := int32(binary.LittleEndian.Uint32(core[28:32]))
	// Variable-length remainder: read_name + cigar + seq + qual + aux.
	remSize := int(blockSize) - coreSize
	if remSize < int(lReadName) {
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: fmt.Sprintf("block_size remainder %d < l_read_name %d", remSize, lReadName)}
	}
	rem := make([]byte, remSize)
	if _, err := io.ReadFull(d.br, rem); err != nil {
		return nil, &ParseError{Line: d.recordLine, Field: "record",
			Msg: "failed to read record remainder: " + err.Error()}
	}
	off := 0
	// Bound the variable-length fields against the actual
	// remainder. htslib checks (nCigarOp*4) + lReadName +
	// ((lSeq+1)>>1) + lSeq <= new_l_data; we add a few
	// additional sanity checks (negative or absurd values
	// would otherwise lead to negative seqBytesLen and a
	// make([]byte, -1) panic).
	seqBytesLen := (int(lSeq) + 1) / 2
	want := int(lReadName) + 4*int(nCigarOp) + seqBytesLen + int(lSeq)
	if want < 0 || want > len(rem) {
		return nil, &ParseError{Line: d.recordLine, Field: "remainder",
			Msg: fmt.Sprintf("l_read_name %d + nCigarOp*4 %d + seq %d + qual %d = %d > remainder %d",
				lReadName, 4*nCigarOp, seqBytesLen, lSeq, want, len(rem))}
	}
	// read_name (NUL-terminated but the length is exact).
	readName := string(rem[off : off+int(lReadName)])
	off += int(lReadName)
	if lReadName > 0 && readName[len(readName)-1] == 0 {
		readName = readName[:len(readName)-1]
	}
	// cigar: nCigarOp * 4 bytes (uint32 LE each).
	cigar := make([]uint32, nCigarOp)
	for i := range cigar {
		cigar[i] = binary.LittleEndian.Uint32(rem[off : off+4])
		off += 4
	}
	// seq: (l_seq+1)/2 bytes, 4-bit packed.
	seqLen := int(lSeq)
	seqBytesLen = (seqLen + 1) / 2
	seq := make([]byte, seqBytesLen)
	copy(seq, rem[off:off+seqBytesLen])
	off += seqBytesLen
	// qual: l_seq bytes.
	var qual []byte
	if lSeq > 0 {
		qual = make([]byte, seqLen)
		copy(qual, rem[off:off+seqLen])
		off += seqLen
	}
	// aux: remainder.
	aux := rem[off:]
	tags, err := parseAuxTags(aux, d.recordLine)
	if err != nil {
		return nil, err
	}
	// Convert BAM's 0-based pos to SAM's 1-based POS.
	pos := int64(pos0) + 1
	nextPos := int64(nextPos0) + 1
	if pos0 < 0 {
		pos = 0
	}
	if nextPos0 < 0 {
		nextPos = 0
	}
	// Compute name lookups.
	nameRName := "*"
	nameRNext := "*"
	if refID >= 0 && int(refID) < len(d.hdr.SQs) {
		nameRName = d.hdr.SQs[refID].Name
	}
	if nextRefID >= 0 && int(nextRefID) < len(d.hdr.SQs) {
		nameRNext = d.hdr.SQs[nextRefID].Name
	}
	rec := &Record{
		Line:      d.recordLine,
		QName:     readName,
		Flag:      flag,
		RefID:     refID,
		Pos:       pos,
		MapQ:      mapq,
		Bin:       bin,
		Cigar:     cigar,
		NextRefID: nextRefID,
		NextPos:   nextPos,
		TLen:      int64(tlen),
		Seq:       decodeSeq(seq, seqLen),
		Qual:      qual,
		Tags:      tags,
		NameRName: nameRName,
		NameRNext: nameRNext,
	}
	return rec, nil
}

// decodeSeq decodes the 4-bit packed SEQ bytes back into an ASCII
// string. Bases follow the SAM/BAM mapping:
//
//	=0 A=1 C=2 M=3 G=4 R=5 S=6 V=7 T=8 W=9 Y=10 H=11 K=12 D=13 B=14 N=15
//
// If lSeq == 0 the returned string is "*".
func decodeSeq(packed []byte, lSeq int) string {
	if lSeq == 0 {
		return "*"
	}
	out := make([]byte, lSeq)
	for i := 0; i < lSeq; i++ {
		b := packed[i/2]
		var n byte
		if i%2 == 0 {
			n = b >> 4
		} else {
			n = b & 0x0F
		}
		out[i] = baseChar(n)
	}
	return string(out)
}

// baseChar returns the ASCII base character for a 4-bit BAM code.
// Unknown codes map to 'N'.
func baseChar(n byte) byte {
	switch n {
	case 0:
		return '='
	case 1:
		return 'A'
	case 2:
		return 'C'
	case 3:
		return 'M'
	case 4:
		return 'G'
	case 5:
		return 'R'
	case 6:
		return 'S'
	case 7:
		return 'V'
	case 8:
		return 'T'
	case 9:
		return 'W'
	case 10:
		return 'Y'
	case 11:
		return 'H'
	case 12:
		return 'K'
	case 13:
		return 'D'
	case 14:
		return 'B'
	case 15:
		return 'N'
	default:
		return 'N'
	}
}

// parseAuxTags parses a slice of bytes containing one or more
// TAG:TAG:TYPE:VALUE auxiliary fields. Returns the parsed tags
// in source order.
func parseAuxTags(b []byte, line int) ([]Tag, error) {
	var tags []Tag
	for len(b) > 0 {
		if len(b) < 3 {
			return nil, &ParseError{Line: line, Field: "tag", Msg: "truncated tag header"}
		}
		key := string(b[0:2])
		tt, err := ParseTagType(b[2])
		if err != nil {
			return nil, &ParseError{Line: line, Field: "tag", Msg: err.Error()}
		}
		b = b[3:]
		val, sub, consumed, err := readTagValue(b, tt, line)
		if err != nil {
			return nil, err
		}
		tags = append(tags, Tag{Key: key, Type: tt, Sub: sub, Value: val})
		b = b[consumed:]
	}
	return tags, nil
}

// readTagValue reads one tag value from b and returns its string
// representation plus the number of bytes consumed from b.
func readTagValue(b []byte, tt TagType, line int) (val string, sub TagType, consumed int, err error) {
	switch tt {
	case TagA:
		if len(b) < 1 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "A: missing value"}
		}
		return string(b[0:1]), 0, 1, nil
	case Tagc:
		if len(b) < 1 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "c: missing value"}
		}
		return strconv.FormatInt(int64(int8(b[0])), 10), 0, 1, nil
	case TagC:
		if len(b) < 1 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "C: missing value"}
		}
		return strconv.FormatUint(uint64(b[0]), 10), 0, 1, nil
	case Tags:
		if len(b) < 2 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "s: missing value"}
		}
		return strconv.FormatInt(int64(int16(binary.LittleEndian.Uint16(b[0:2]))), 10), 0, 2, nil
	case TagS:
		if len(b) < 2 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "S: missing value"}
		}
		return strconv.FormatUint(uint64(binary.LittleEndian.Uint16(b[0:2])), 10), 0, 2, nil
	case Tagi:
		if len(b) < 4 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "i: missing value"}
		}
		return strconv.FormatInt(int64(int32(binary.LittleEndian.Uint32(b[0:4]))), 10), 0, 4, nil
	case TagI:
		if len(b) < 4 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "I: missing value"}
		}
		return strconv.FormatUint(uint64(binary.LittleEndian.Uint32(b[0:4])), 10), 0, 4, nil
	case Tagf:
		if len(b) < 4 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "f: missing value"}
		}
		bits := binary.LittleEndian.Uint32(b[0:4])
		return strconv.FormatFloat(float64(math.Float32frombits(bits)), 'g', -1, 32), 0, 4, nil
	case Tagd:
		if len(b) < 8 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "d: missing value"}
		}
		bits := binary.LittleEndian.Uint64(b[0:8])
		return strconv.FormatFloat(math.Float64frombits(bits), 'g', -1, 64), 0, 8, nil
	case TagZ:
		end := bytes.IndexByte(b, 0)
		if end < 0 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "Z: missing null terminator"}
		}
		return string(b[:end]), 0, end + 1, nil
	case TagH:
		end := bytes.IndexByte(b, 0)
		if end < 0 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "H: missing null terminator"}
		}
		return string(b[:end]), 0, end + 1, nil
	case TagB:
		// B array: 1-byte sub_type, uint32 LE count, then elements.
		if len(b) < 5 {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "B: missing sub_type/count"}
		}
		subType := TagType(b[0])
		if !subType.IsValid() || subType == TagZ || subType == TagH || subType == TagB || subType == TagA || subType == Tagd {
			return "", 0, 0, &ParseError{Line: line, Field: "tag",
				Msg: fmt.Sprintf("B: invalid sub_type %q", byte(subType))}
		}
		count := binary.LittleEndian.Uint32(b[1:5])
		eltSize, err := bEltSize(subType)
		if err != nil {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: err.Error()}
		}
		dataLen := int(count) * eltSize
		if len(b) < 5+dataLen {
			return "", 0, 0, &ParseError{Line: line, Field: "tag",
				Msg: fmt.Sprintf("B: data truncated (need %d, have %d)", 5+dataLen, len(b))}
		}
		val, err := formatBArray(b[5:5+dataLen], subType, count)
		if err != nil {
			return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: err.Error()}
		}
		return val, subType, 5 + dataLen, nil
	default:
		return "", 0, 0, &ParseError{Line: line, Field: "tag", Msg: "unknown tag type"}
	}
}

// bEltSize returns the byte size of one element in a B array of
// the given sub-type.
func bEltSize(sub TagType) (int, error) {
	switch sub {
	case Tagc, TagC:
		return 1, nil
	case Tags, TagS:
		return 2, nil
	case Tagi, TagI, Tagf:
		return 4, nil
	default:
		return 0, fmt.Errorf("tag: B sub-type %q unsupported", byte(sub))
	}
}

// formatBArray formats a B-array's data bytes into a
// comma-separated string of element values.
func formatBArray(data []byte, sub TagType, count uint32) (string, error) {
	var parts []string
	switch sub {
	case Tagc:
		for i := 0; i < int(count); i++ {
			parts = append(parts, strconv.FormatInt(int64(int8(data[i])), 10))
		}
	case TagC:
		for i := 0; i < int(count); i++ {
			parts = append(parts, strconv.FormatUint(uint64(data[i]), 10))
		}
	case Tags:
		for i := 0; i < int(count); i++ {
			parts = append(parts, strconv.FormatInt(int64(int16(binary.LittleEndian.Uint16(data[i*2:(i+1)*2]))), 10))
		}
	case TagS:
		for i := 0; i < int(count); i++ {
			parts = append(parts, strconv.FormatUint(uint64(binary.LittleEndian.Uint16(data[i*2:(i+1)*2])), 10))
		}
	case Tagi:
		for i := 0; i < int(count); i++ {
			parts = append(parts, strconv.FormatInt(int64(int32(binary.LittleEndian.Uint32(data[i*4:(i+1)*4]))), 10))
		}
	case TagI:
		for i := 0; i < int(count); i++ {
			parts = append(parts, strconv.FormatUint(uint64(binary.LittleEndian.Uint32(data[i*4:(i+1)*4])), 10))
		}
	case Tagf:
		for i := 0; i < int(count); i++ {
			bits := binary.LittleEndian.Uint32(data[i*4 : (i+1)*4])
			parts = append(parts, strconv.FormatFloat(float64(math.Float32frombits(bits)), 'g', -1, 32))
		}
	default:
		return "", fmt.Errorf("tag: B sub-type %q unsupported", byte(sub))
	}
	return strings.Join(parts, ","), nil
}

// parseHeaderText parses the SAM-formatted header text embedded
// in the BAM binary header. It populates hdr.HD / hdr.RGs /
// hdr.PGs / hdr.COs and leaves hdr.Raw / hdr.SQs untouched.
//
// The parser is a small local copy of [sam]'s header-line rules
// (per the per-package self-contained design principle): each
// non-blank line is "@XX<TAB>key:value..." with XX being a
// two-character type tag. Lines without a leading "@" or with an
// unknown type tag are ignored.
func parseHeaderText(hdr *Header) {
	for _, line := range strings.Split(hdr.Raw, "\n") {
		if line == "" {
			continue
		}
		if line[0] != '@' {
			continue
		}
		if len(line) < 3 {
			continue
		}
		if len(line) > 3 && line[3] != '\t' {
			// Type tag longer than 2 chars — skip silently
			// (the binary spec doesn't enforce a 2-char type
			// for the embedded text, only for the SAM header
			// in the SAM file format itself).
			continue
		}
		tag := line[1:3]
		rest := ""
		if len(line) > 3 {
			rest = strings.TrimPrefix(line[3:], "\t")
		}
		pairs := parseHeaderPairs(rest)
		switch tag {
		case "HD":
			if hdr.HD != nil {
				continue // spec allows exactly one
			}
			hdr.HD = &HD{
				VN:    pairs["VN"],
				SO:    pairs["SO"],
				GO:    pairs["GO"],
				Other: pairs.otherExcept("VN", "SO", "GO"),
				Raw:   line,
			}
		case "RG":
			hdr.RGs = append(hdr.RGs, &RG{
				ID:    pairs["ID"],
				SM:    pairs["SM"],
				LB:    pairs["LB"],
				BC:    pairs["BC"],
				DS:    pairs["DS"],
				FO:    pairs["FO"],
				KS:    pairs["KS"],
				PU:    pairs["PU"],
				PI:    pairs["PI"],
				CN:    pairs["CN"],
				DT:    pairs["DT"],
				PL:    pairs["PL"],
				PG:    pairs["PG"],
				PM:    pairs["PM"],
				Other: pairs.otherExcept("ID", "SM", "LB", "BC", "DS", "FO", "KS", "PU", "PI", "CN", "DT", "PL", "PG", "PM"),
				Raw:   line,
			})
		case "PG":
			hdr.PGs = append(hdr.PGs, &PG{
				ID:    pairs["ID"],
				PN:    pairs["PN"],
				CL:    pairs["CL"],
				DS:    pairs["DS"],
				PP:    pairs["PP"],
				VN:    pairs["VN"],
				Other: pairs.otherExcept("ID", "PN", "CL", "DS", "PP", "VN"),
				Raw:   line,
			})
		case "CO":
			hdr.COs = append(hdr.COs, line)
		case "SQ":
			// @SQ in the embedded text is redundant — the
			// binary reference dictionary is the canonical
			// source. Skip the text version to avoid
			// confusion.
			continue
		default:
			// Unknown @XX line: the embedded text is
			// preserved verbatim in hdr.Raw for any future
			// round-trip use. We don't surface parsed
			// structures for type tags we don't recognise.
		}
	}
}

// headerPairs is a key=value map; the SAM header format is a
// flat tab/space-separated list of "key:value" pairs.
type headerPairs map[string]string

func (m headerPairs) otherExcept(keys ...string) map[string]string {
	excluded := make(map[string]struct{}, len(keys))
	for _, k := range keys {
		excluded[k] = struct{}{}
	}
	out := make(map[string]string)
	for k, v := range m {
		if _, ok := excluded[k]; ok {
			continue
		}
		out[k] = v
	}
	return out
}

func parseHeaderPairs(s string) headerPairs {
	out := make(headerPairs)
	i := 0
	for i < len(s) {
		// skip whitespace
		for i < len(s) && (s[i] == ' ' || s[i] == '\t') {
			i++
		}
		if i >= len(s) {
			break
		}
		kstart := i
		for i < len(s) && s[i] != ':' && s[i] != ' ' && s[i] != '\t' {
			i++
		}
		key := s[kstart:i]
		if i >= len(s) || s[i] != ':' {
			for i < len(s) && s[i] != ' ' && s[i] != '\t' {
				i++
			}
			continue
		}
		i++
		vstart := i
		for i < len(s) && s[i] != ' ' && s[i] != '\t' {
			i++
		}
		out[key] = s[vstart:i]
	}
	return out
}
