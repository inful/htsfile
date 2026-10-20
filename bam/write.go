package bam

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strconv"

	"github.com/inful/htsfile/internal/bgzf"
)

// Encoder writes a BAM document to an underlying [io.Writer].
// The output is BGZF-compressed per the BAM spec.
//
// Typical usage:
//
//	f, _ := os.Create("out.bam")
//	defer f.Close()
//	enc := bam.NewEncoder(f)
//	if err := enc.WriteHeader(hdr); err != nil { ... }
//	for _, rec := range recs { enc.WriteRecord(rec) }
//	enc.Close()
//
// Records may be written in any order (the writer does not
// assume sorted input). The writer flushes BGZF blocks as
// they fill; Close emits the canonical 28-byte EOF marker.
//
// An Encoder is not safe for concurrent use.
type Encoder struct {
	w *bgzf.Writer

	hdrWritten bool
	closed     bool
}

// NewEncoder wraps w in an [Encoder].
func NewEncoder(w io.Writer) *Encoder {
	return &Encoder{w: bgzf.NewWriter(w)}
}

// WriteHeader emits the BAM file header: magic + l_text +
// SAM header text + n_ref + reference dictionary. Must be
// called before any WriteRecord call. A nil header returns
// an error.
func (e *Encoder) WriteHeader(hdr *Header) error {
	if hdr == nil {
		return errors.New("bam: WriteHeader: nil header")
	}
	if e.hdrWritten {
		return errors.New("bam: WriteHeader called twice")
	}
	e.hdrWritten = true

	// Build the BGZF block payload: magic + l_text + text +
	// n_ref + reference dictionary.
	var payload bytes.Buffer
	payload.Write(bamMagic[:])
	var lt [4]byte
	binary.LittleEndian.PutUint32(lt[:], uint32(len(hdr.Raw)))
	payload.Write(lt[:])
	payload.WriteString(hdr.Raw)
	var nr [4]byte
	binary.LittleEndian.PutUint32(nr[:], uint32(len(hdr.SQs)))
	payload.Write(nr[:])
	for _, sq := range hdr.SQs {
		name := sq.Name
		if name == "" {
			name = "*"
		}
		// Null-terminate. The spec says the name "ends with a
		// NUL" but the length field counts the NUL.
		nullName := name + "\x00"
		var ln [4]byte
		binary.LittleEndian.PutUint32(ln[:], uint32(len(nullName)))
		payload.Write(ln[:])
		payload.WriteString(nullName)
		var ll [4]byte
		binary.LittleEndian.PutUint32(ll[:], uint32(sq.Len))
		payload.Write(ll[:])
	}
	if _, err := e.w.Write(payload.Bytes()); err != nil {
		return fmt.Errorf("bam: writing header: %w", err)
	}
	return nil
}

// WriteRecord emits one alignment record. The record's
// RefID, Pos, Flag, MapQ, Bin, Cigar, NextRefID, NextPos,
// TLen, Seq, Qual, and Tags fields are encoded into a single
// BGZF block (one record per block, per the BAM spec). A nil
// record is a no-op.
func (e *Encoder) WriteRecord(r *Record) error {
	if r == nil {
		return nil
	}
	if !e.hdrWritten {
		return errors.New("bam: WriteRecord before WriteHeader")
	}
	if e.closed {
		return errors.New("bam: WriteRecord after Close")
	}
	payload, err := encodeRecord(r)
	if err != nil {
		return fmt.Errorf("bam: encoding record: %w", err)
	}
	if _, err := e.w.Write(payload); err != nil {
		return fmt.Errorf("bam: writing record: %w", err)
	}
	return nil
}

// Close flushes any buffered BGZF data and emits the
// canonical 28-byte EOF marker. A second Close is a no-op.
func (e *Encoder) Close() error {
	if e.closed {
		return nil
	}
	e.closed = true
	return e.w.Close()
}

// encodeRecord serializes one Record into the BGZF block
// payload: 4-byte block_size + 32-byte core + variable-length
// remainder (read_name + cigar + seq + qual + aux tags).
//
// The block_size field counts everything from itself through
// the end of the aux tag section.
func encodeRecord(r *Record) ([]byte, error) {
	// Build the variable-length remainder first so we know
	// block_size when we emit the header.
	var rem bytes.Buffer
	// read_name: NUL-terminated. l_read_name counts the NUL.
	rem.WriteString(r.QName)
	rem.WriteByte(0)
	// cigar: each op is uint32 LE.
	for _, op := range r.Cigar {
		var c [4]byte
		binary.LittleEndian.PutUint32(c[:], op)
		rem.Write(c[:])
	}
	// seq: 4-bit packed. l_seq = 0 means "*"; otherwise
	// pack the ASCII bases into 4-bit codes.
	lSeq := int32(0)
	if r.Seq != "" && r.Seq != "*" {
		lSeq = int32(len(r.Seq))
	}
	if lSeq > 0 {
		packed := encodeSeq(r.Seq)
		rem.Write(packed)
	}
	// qual: ASCII bytes, length = l_seq (if l_seq > 0).
	if lSeq > 0 {
		switch {
		case len(r.Qual) == int(lSeq):
			rem.Write(r.Qual)
		case len(r.Qual) > int(lSeq):
			// Spec: qual length must equal l_seq when seq is
			// present. If the caller gave us a longer qual,
			// truncate.
			rem.Write(r.Qual[:lSeq])
		case len(r.Qual) > 0:
			// Pad with 'I' (Phred 40) when qual is shorter.
			rem.Write(r.Qual)
			rem.Write(bytes.Repeat([]byte{'I'}, int(lSeq)-len(r.Qual)))
		default:
			// Default to 'I' (Phred 40) when no qual supplied.
			rem.Write(bytes.Repeat([]byte{'I'}, int(lSeq)))
		}
	}
	// aux tags.
	if aux, err := encodeAuxTags(r.Tags); err != nil {
		return nil, err
	} else if len(aux) > 0 {
		rem.Write(aux)
	}
	// Core: 32 bytes.
	core := make([]byte, 32)
	binary.LittleEndian.PutUint32(core[0:4], uint32(r.RefID))
	pos0 := uint32(int32(r.Pos) - 1)
	if r.Pos <= 0 {
		pos0 = 0xFFFFFFFF // -1 in int32 (unmapped)
	}
	binary.LittleEndian.PutUint32(core[4:8], pos0)
	// QNAME: spec range [1, 254] including NUL. Refuse to
	// silently truncate; callers that have a longer name
	// should use a truncated convention like htslib does.
	if len(r.QName) > 254-1 {
		return nil, fmt.Errorf("QNAME length %d exceeds BAM max (254-1)", len(r.QName))
	}
	core[8] = uint8(len(r.QName)) + 1 // l_read_name includes NUL
	core[9] = r.MapQ
	binary.LittleEndian.PutUint16(core[10:12], r.Bin)
	// nCigarOp: spec max 65535. Longer CIGARs must be
	// stored in a CG:B,I tag (htslib behaviour). Refusing
	// here matches that, instead of silent truncation.
	if len(r.Cigar) > 0xFFFF {
		return nil, fmt.Errorf("n_cigar_op %d exceeds BAM max (65535); use a CG:B,I tag", len(r.Cigar))
	}
	binary.LittleEndian.PutUint16(core[12:14], uint16(len(r.Cigar)))
	binary.LittleEndian.PutUint16(core[14:16], r.Flag)
	binary.LittleEndian.PutUint32(core[16:20], uint32(lSeq))
	binary.LittleEndian.PutUint32(core[20:24], uint32(r.NextRefID))
	nextPos0 := uint32(int32(r.NextPos) - 1)
	if r.NextPos <= 0 {
		nextPos0 = 0xFFFFFFFF
	}
	binary.LittleEndian.PutUint32(core[24:28], nextPos0)
	binary.LittleEndian.PutUint32(core[28:32], uint32(r.TLen))
	// block_size = 32 (core) + len(rem). The block_size
	// field itself is not counted (it appears in the 4
	// bytes immediately before the core).
	blockSize := uint32(32 + rem.Len())
	if blockSize < 32 {
		return nil, fmt.Errorf("bam: block_size underflow")
	}
	var bs [4]byte
	binary.LittleEndian.PutUint32(bs[:], blockSize)
	out := make([]byte, 0, blockSize)
	out = append(out, bs[:]...)
	out = append(out, core...)
	out = append(out, rem.Bytes()...)
	return out, nil
}

// encodeSeq packs an ASCII sequence into BAM's 4-bit code.
// Unknown bases (anything outside ACGTN ambiguity codes)
// encode as N (15).
func encodeSeq(s string) []byte {
	n := len(s)
	packed := make([]byte, (n+1)/2)
	for i := 0; i < n; i++ {
		code := baseCode(s[i])
		if i%2 == 0 {
			packed[i/2] = code << 4
		} else {
			packed[i/2] |= code & 0x0F
		}
	}
	return packed
}

// baseCode returns the 4-bit BAM code for an ASCII base char.
// Unknown chars map to N (15).
func baseCode(c byte) byte {
	switch c {
	case '=':
		return 0
	case 'A', 'a':
		return 1
	case 'C', 'c':
		return 2
	case 'M', 'm':
		return 3
	case 'G', 'g':
		return 4
	case 'R', 'r':
		return 5
	case 'S', 's':
		return 6
	case 'V', 'v':
		return 7
	case 'T', 't':
		return 8
	case 'W', 'w':
		return 9
	case 'Y', 'y':
		return 10
	case 'H', 'h':
		return 11
	case 'K', 'k':
		return 12
	case 'D', 'd':
		return 13
	case 'B', 'b':
		return 14
	case 'N', 'n':
		return 15
	default:
		return 15
	}
}

// encodeAuxTags serialises a slice of tags into BAM binary form.
// Each tag is laid out as: 2-byte key + 1-byte type + value
// bytes (type-specific). Returns the concatenation of all
// tags' encodings, or an error if a tag's value doesn't
// parse for its declared type.
func encodeAuxTags(tags []Tag) ([]byte, error) {
	var buf bytes.Buffer
	for i, t := range tags {
		if len(t.Key) != 2 {
			return nil, fmt.Errorf("tag[%d]: key %q is not 2 bytes", i, t.Key)
		}
		if !t.Type.IsValid() {
			return nil, fmt.Errorf("tag[%d] %q: invalid type %q", i, t.Key, t.Type)
		}
		buf.WriteString(t.Key)
		buf.WriteByte(byte(t.Type))
		val, err := encodeTagValue(t)
		if err != nil {
			return nil, fmt.Errorf("tag[%d] %q: %w", i, t.Key, err)
		}
		buf.Write(val)
	}
	return buf.Bytes(), nil
}

// encodeTagValue serialises a single tag value. The Value
// field is the SAM-compatible string form; this function
// parses it back to bytes per the tag's declared Type.
func encodeTagValue(t Tag) ([]byte, error) {
	var buf bytes.Buffer
	switch t.Type {
	case TagA:
		if t.Value != "" {
			buf.WriteByte(t.Value[0])
		} else {
			buf.WriteByte(0)
		}
	case Tagc:
		v, err := strconv.ParseInt(t.Value, 10, 8)
		if err != nil {
			return nil, fmt.Errorf("c: %w", err)
		}
		buf.WriteByte(byte(int8(v)))
	case TagC:
		v, err := strconv.ParseUint(t.Value, 10, 8)
		if err != nil {
			return nil, fmt.Errorf("C: %w", err) //nolint:staticcheck // tag-type prefix, not a sentence
		}
		buf.WriteByte(byte(v))
	case Tags:
		v, err := strconv.ParseInt(t.Value, 10, 16)
		if err != nil {
			return nil, fmt.Errorf("s: %w", err)
		}
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(int16(v)))
		buf.Write(b[:])
	case TagS:
		v, err := strconv.ParseUint(t.Value, 10, 16)
		if err != nil {
			return nil, fmt.Errorf("S: %w", err) //nolint:staticcheck // tag-type prefix, not a sentence
		}
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(v))
		buf.Write(b[:])
	case Tagi:
		v, err := strconv.ParseInt(t.Value, 10, 32)
		if err != nil {
			return nil, fmt.Errorf("i: %w", err)
		}
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(int32(v)))
		buf.Write(b[:])
	case TagI:
		v, err := strconv.ParseUint(t.Value, 10, 32)
		if err != nil {
			return nil, fmt.Errorf("I: %w", err) //nolint:staticcheck // tag-type prefix, not a sentence
		}
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(v))
		buf.Write(b[:])
	case Tagf:
		f, err := strconv.ParseFloat(t.Value, 32)
		if err != nil {
			return nil, fmt.Errorf("f: %w", err)
		}
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], math.Float32bits(float32(f)))
		buf.Write(b[:])
	case Tagd:
		f, err := strconv.ParseFloat(t.Value, 64)
		if err != nil {
			return nil, fmt.Errorf("d: %w", err)
		}
		var b [8]byte
		binary.LittleEndian.PutUint64(b[:], math.Float64bits(f))
		buf.Write(b[:])
	case TagZ, TagH:
		buf.WriteString(t.Value)
		buf.WriteByte(0)
	case TagB:
		// B array: 1-byte sub-type + 4-byte count + elements.
		// The sub-type lives in Tag.Sub, not Tag.Type (Tag.Type
		// stays TagB for the wire). The reader refuses unknown
		// sub-types, so we refuse them here too: emitting a
		// B array whose Sub is not one of c, C, s, S, i, I, f
		// would round-trip as a parse error.
		sub := t.Sub
		if !isValidBSubType(sub) {
			return nil, fmt.Errorf("bam: tag %q: B array sub-type %q invalid (must be one of c, C, s, S, i, I, f)", t.Key, byte(sub))
		}
		buf.WriteByte(byte(sub))
		parts := splitComma(t.Value)
		var c [4]byte
		binary.LittleEndian.PutUint32(c[:], uint32(len(parts)))
		buf.Write(c[:])
		for _, p := range parts {
			switch sub {
			case Tagc:
				v, err := strconv.ParseInt(p, 10, 8)
				if err != nil {
					return nil, fmt.Errorf("B:c: %w", err) //nolint:staticcheck // tag-type prefix
				}
				buf.WriteByte(byte(int8(v)))
			case TagC:
				v, err := strconv.ParseUint(p, 10, 8)
				if err != nil {
					return nil, fmt.Errorf("B:C: %w", err)
				}
				buf.WriteByte(byte(v))
			case Tags:
				v, err := strconv.ParseInt(p, 10, 16)
				if err != nil {
					return nil, fmt.Errorf("B:s: %w", err) //nolint:staticcheck // tag-type prefix
				}
				var b [2]byte
				binary.LittleEndian.PutUint16(b[:], uint16(int16(v)))
				buf.Write(b[:])
			case TagS:
				v, err := strconv.ParseUint(p, 10, 16)
				if err != nil {
					return nil, fmt.Errorf("B:S: %w", err)
				}
				var b [2]byte
				binary.LittleEndian.PutUint16(b[:], uint16(v))
				buf.Write(b[:])
			case Tagi:
				v, err := strconv.ParseInt(p, 10, 32)
				if err != nil {
					return nil, fmt.Errorf("B:i: %w", err) //nolint:staticcheck // tag-type prefix
				}
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], uint32(int32(v)))
				buf.Write(b[:])
			case TagI:
				v, err := strconv.ParseUint(p, 10, 32)
				if err != nil {
					return nil, fmt.Errorf("B:I: %w", err)
				}
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], uint32(v))
				buf.Write(b[:])
			case Tagf:
				f, err := strconv.ParseFloat(p, 32)
				if err != nil {
					return nil, fmt.Errorf("B:f: %w", err) //nolint:staticcheck // tag-type prefix
				}
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], math.Float32bits(float32(f)))
				buf.Write(b[:])
			}
		}
	default:
		return nil, fmt.Errorf("unknown tag type %q", t.Type)
	}
	return buf.Bytes(), nil
}

// splitComma splits s on commas and trims surrounding spaces
// from each element. Empty fields are dropped (so ",,"
// yields zero elements, matching how empty arrays are encoded).
func splitComma(s string) []string {
	if s == "" {
		return nil
	}
	var out []string
	for _, p := range bytes.Split([]byte(s), []byte{','}) {
		p = bytes.TrimSpace(p)
		if len(p) > 0 {
			out = append(out, string(p))
		}
	}
	return out
}

func isValidBSubType(b TagType) bool {
	switch b {
	case Tagc, TagC, Tags, TagS, Tagi, TagI, Tagf:
		return true
	}
	return false
}
