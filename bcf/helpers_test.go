package bcf

import (
	"bytes"
	"compress/flate"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"io"
	"math"
	"os"
	"strings"
	"testing"
)

// _ keeps the fmt import used (the panic message above is
// gated behind hdrIdx < 0; the import would otherwise be unused
// on the happy path).
var _ = fmt.Sprintf

// builder is a test helper for constructing synthetic BCF
// streams wrapped in BGZF. We write to a temp file because
// the bgzf package's reader streams through the underlying
// io.Reader; a file is the simplest backing store that works
// with that model.
type builder struct {
	path       string
	f          *os.File
	hdr        *Header  // parsed sampleHeader; used for LookupInfoIndex etc.
	hdrOrder   []string // cached: hdr.InfosOrder
	hdrFilters []string
	hdrContigs []string
}

func newBuilder(t *testing.T) *builder {
	return newBuilderWithHeader(t, sampleHeader)
}

// newBuilderWithHeader builds a builder that writes the given
// header text. Use this when a test needs a header that
// differs from sampleHeader (e.g. for variable-count INFO
// fields, additional ##FILTER declarations, etc.).
func newBuilderWithHeader(t *testing.T, headerText string) *builder {
	t.Helper()
	f, err := os.CreateTemp(t.TempDir(), "bcf-test-*.bcf")
	if err != nil {
		t.Fatal(err)
	}
	hdr := ParseHeader(headerText)
	return &builder{
		path:       f.Name(),
		f:          f,
		hdr:        hdr,
		hdrOrder:   hdr.InfosOrder,
		hdrFilters: hdr.Filters,
		hdrContigs: hdr.Contigs,
	}
}

// writeBGZFBlock writes one BGZF block containing payload
// with valid CRC32 + ISIZE trailer.
func (b *builder) writeBGZFBlock(payload []byte) {
	var cdata bytes.Buffer
	w, _ := flate.NewWriter(&cdata, flate.DefaultCompression)
	_, _ = w.Write(payload)
	_ = w.Close()
	cd := cdata.Bytes()
	bsize := uint16(18 + len(cd) + 8 - 1)
	buf := bytes.Buffer{}
	buf.Write([]byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00,
	})
	var bsizeBytes [2]byte
	binary.LittleEndian.PutUint16(bsizeBytes[:], bsize)
	buf.Write(bsizeBytes[:])
	buf.Write(cd)
	var crcBytes [4]byte
	binary.LittleEndian.PutUint32(crcBytes[:], crc32.ChecksumIEEE(payload))
	buf.Write(crcBytes[:])
	var isBytes [4]byte
	binary.LittleEndian.PutUint32(isBytes[:], uint32(len(payload)))
	buf.Write(isBytes[:])
	_, _ = b.f.Write(buf.Bytes())
}

// writeEOF writes the 28-byte empty BGZF EOF marker.
func (b *builder) writeEOF() {
	_, _ = b.f.Write([]byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00,
	})
}

// finish closes the file and returns its contents.
func (b *builder) finish(t *testing.T) []byte {
	t.Helper()
	b.writeEOF()
	if err := b.f.Close(); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(b.path)
	if err != nil {
		t.Fatal(err)
	}
	return data
}

// writeBCFHeader writes the BCF magic + length-prefixed
// header text wrapped in a single BGZF block.
//
//nolint:unparam // text is parameterized for testability even though all callers pass sampleHeader
func (b *builder) writeBCFHeader(text string) {
	var payload bytes.Buffer
	payload.Write(Magic[:])
	// l_text counts the NUL terminator (spec BCF2 header).
	var l [4]byte
	binary.LittleEndian.PutUint32(l[:], uint32(len(text)+1))
	payload.Write(l[:])
	payload.WriteString(text)
	payload.WriteByte(0)
	b.writeBGZFBlock(payload.Bytes())
}

// testRecord is the in-test struct used to construct BCF
// records. It carries the pre-resolved header indices the
// binary encoder needs (CHROM idx, INFO keys in declaration
// order).
type testRecord struct {
	contigIdx  int
	Pos        int64
	ID         string
	Ref        string
	Alt        []string
	Qual       float64
	QualKnown  bool
	filterIdxs []int
	Info       map[string]InfoValue
	infoDefs   map[string]*InfoDef
	infoOrder  []string
	// formatIdx lists the FORMAT keys present on this record,
	// as indices into the header's FormatsOrder.
	formatIdx []uint16
	// nSample / nFmt are the per-sample block's header counts.
	nSample int
	nFmt    uint8
	// sampleVals holds the per-sample values, keyed by FORMAT
	// index.
	sampleVals []map[uint16]SampleValue
}

// writeRecord encodes one BCF record into a fresh BGZF
// block.
func (b *builder) writeRecord(t *testing.T, r *testRecord) {
	t.Helper()
	var shared bytes.Buffer
	// Spec section 5, BcfSiteEncoding: 24-byte site prefix
	// CHROM, POS, rlen, QUAL, n_info, n_allele, n_sample,
	// n_fmt. n_sample is a uint24 at [20:22] and n_fmt a
	// uint8 at [23].
	nInfo := uint16(len(r.Info))
	nAllele := uint16(1 + len(r.Alt))
	var prefix [24]byte
	binary.LittleEndian.PutUint32(prefix[0:4], uint32(r.contigIdx))
	binary.LittleEndian.PutUint32(prefix[4:8], uint32(int32(r.Pos)-1))
	binary.LittleEndian.PutUint32(prefix[8:12], uint32(len(r.Ref))) // rlen
	qualBits := bcfFloatMissing
	if r.QualKnown {
		qualBits = math.Float32bits(float32(r.Qual))
	}
	binary.LittleEndian.PutUint32(prefix[12:16], qualBits)
	binary.LittleEndian.PutUint16(prefix[16:18], nInfo)
	binary.LittleEndian.PutUint16(prefix[18:20], nAllele)
	binary.LittleEndian.PutUint32(prefix[20:24], uint32(r.nSample)<<24|uint32(r.nFmt))
	shared.Write(prefix[:])

	testWriteTypedString(&shared, r.ID)
	testWriteTypedString(&shared, r.Ref)
	for _, alt := range r.Alt {
		testWriteTypedString(&shared, alt)
	}
	testWriteTypedFilter(&shared, b.hdr, r.filterIdxs)
	for _, id := range r.infoOrder {
		def := r.infoDefs[id]
		val := r.Info[id]
		// INFO key is a typed atomic int holding the header
		// dictionary index (spec §"Dictionary of strings").
		testWriteTypedInt(&shared, int32(b.hdr.IDXFor(id)))
		testWriteTypedInfoValue(&shared, def, val)
	}

	// Per-sample block: one triplet per FORMAT field -
	// fmt_key (typed int dictionary index), fmt_type,
	// fmt_values.
	var indiv bytes.Buffer
	for _, keyIdx := range r.formatIdx {
		name := ""
		if int(keyIdx) < len(b.hdr.FormatsOrder) {
			name = b.hdr.FormatsOrder[keyIdx]
		}
		testWriteTypedInt(&indiv, int32(b.hdr.IDXFor(name)))
		testWriteFormatField(&indiv, b.hdr, r, keyIdx)
	}

	// l_shared + l_indiv framing.
	var lens [8]byte
	binary.LittleEndian.PutUint32(lens[0:4], uint32(shared.Len()))
	binary.LittleEndian.PutUint32(lens[4:8], uint32(indiv.Len()))
	var rec bytes.Buffer
	rec.Write(lens[:])
	rec.Write(shared.Bytes())
	rec.Write(indiv.Bytes())
	b.writeBGZFBlock(rec.Bytes())
}

// testWriteTypedInt emits a typed atomic int32: the descriptor
// byte (size 1, type int32) followed by the little-endian
// value.
func testWriteTypedInt(b *bytes.Buffer, v int32) {
	b.WriteByte((1 << 4) | typeInt32)
	var kb [4]byte
	binary.LittleEndian.PutUint32(kb[:], uint32(v))
	b.Write(kb[:])
}

// testWriteTypedString emits a BCF2 string as a char vector,
// matching the spec (§BcfTypeEncoding, tex line 2276) and
// `bcftools view -Ob` output: the descriptor nibble carries the
// length and there is no separate length prefix.
//
// This deliberately does NOT call the package's own
// writeBCFString helper. Fixtures that mirror the code under
// test cannot catch a change to it - which is exactly how the
// pre-2026-10-22 string encoding shipped broken.
func testWriteTypedString(b *bytes.Buffer, s string) {
	if s == "" {
		b.WriteByte(typeMissing)
		return
	}
	b.WriteByte(byte(len(s))<<4 | typeString)
	b.WriteString(s)
}

func testWriteTypedFilter(b *bytes.Buffer, hdr *Header, idxs []int) {
	if len(idxs) == 0 {
		b.WriteByte(typeMissing)
		return
	}
	// Spec §BcfTypeEncoding: FILTER is a typed *vector* of
	// dictionary offsets, so the descriptor nibble carries the
	// element count. Verified byte-for-byte against
	// `bcftools view -Ob` output.
	size, overflow := testPickBcfSize(len(idxs), false)
	b.WriteByte((uint8(size) << 4) | typeInt32)
	if overflow {
		b.WriteByte((1 << 4) | typeInt8)
		b.WriteByte(byte(len(idxs)))
	}
	// idxs are filter-local (0 = PASS); the wire carries the
	// header dictionary index.
	for i, idx := range idxs {
		var kb [4]byte
		binary.LittleEndian.PutUint32(kb[:], uint32(filterIDXFor(hdr, i, idx)))
		b.Write(kb[:])
	}
}

// filterIDXFor maps a filter-local index to the BCF dictionary
// index using the builder's header.
func filterIDXFor(hdr *Header, _ int, local int) int {
	return hdr.FilterLocalToIDX(local)
}

func testWriteTypedInfoValue(b *bytes.Buffer, def *InfoDef, v InfoValue) {
	if def.Type == "Flag" {
		if v.Flag {
			b.WriteByte(typeInt8)
		} else {
			b.WriteByte(typeMissing)
		}
		return
	}
	defSize, isVar := ParseNumber(def.Number)
	count, hasCount := 0, false
	switch def.Type {
	case "Integer":
		count, hasCount = len(v.Ints), true
	case "Float":
		count, hasCount = len(v.Floats), true
	case "String":
		// Spec/hstlib: one char run with list elements joined
		// by commas; the descriptor nibble is the run's total
		// character count. Matches `bcftools view -Ob`.
		count, hasCount = len(strings.Join(v.Strs, ",")), true
	case "Character":
		count, hasCount = 1, true
	}
	if !hasCount {
		count, isVar = defSize, false
	}
	if !hasCount && defSize == 0 && !isVar {
		count = 1
	}
	size, overflow := testPickBcfSize(count, isVar)
	var lowType uint8
	switch def.Type {
	case "Integer":
		lowType = typeInt32
	case "Float":
		lowType = typeFloat
	case "String", "Character":
		lowType = typeString
	}
	b.WriteByte((uint8(size) << 4) | lowType)
	if overflow {
		b.WriteByte((1 << 4) | typeInt8)
		b.WriteByte(byte(count))
	}
	switch def.Type {
	case "Integer":
		for _, n := range v.Ints {
			var x [4]byte
			binary.LittleEndian.PutUint32(x[:], uint32(int32(n)))
			b.Write(x[:])
		}
	case "Float":
		for _, f := range v.Floats {
			var x [4]byte
			binary.LittleEndian.PutUint32(x[:], float32Bits(float32(f)))
			b.Write(x[:])
		}
	case "String":
		b.WriteString(strings.Join(v.Strs, ","))
	case "Character":
		var l [4]byte
		binary.LittleEndian.PutUint32(l[:], uint32(len(v.Str)))
		b.Write(l[:])
		b.WriteString(v.Str)
	}
}

// testPickBcfSize mirrors pickBcfSize for the test fixture
// builder. Independent copy because tests live in a _test
// file and cannot reuse package-private helpers.
func testPickBcfSize(count int, isVar bool) (size int, overflow bool) {
	if isVar || count >= 15 {
		return 15, true
	}
	return count, false
}

// float64Bits / float32Bits wrap math.Float64bits /
// Float32bits so the helpers don't need a separate math
// import.
func float32Bits(f float32) uint32 { return math.Float32bits(f) }

// readAll reads all records from a BCF byte slice and
// returns the parsed header + record list. Used by the
// individual test bodies in parse_test.go.
func readAll(t *testing.T, data []byte) (*Header, []*Record) {
	t.Helper()
	dec := NewDecoder(bytes.NewReader(data))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	var recs []*Record
	for {
		rec, err := dec.Record()
		if readAllEOF(err) {
			break
		}
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		recs = append(recs, rec)
	}
	return hdr, recs
}

// readAllEOF reports whether err is an io.EOF (the standard
// end-of-stream signal from [Decoder.Record]).
func readAllEOF(err error) bool {
	return errors.Is(err, io.EOF)
}

// testWriteFormatField writes one FORMAT field's fmt_type
// descriptor and fmt_values, following the spec
// (§GenotypeEncoding) independently of the package's own
// encoder. Fixtures that call the code under test cannot catch
// a change to it.
func testWriteFormatField(b *bytes.Buffer, hdr *Header, r *testRecord, keyIdx uint16) {
	def, ok := hdr.FormatAt(int(keyIdx))
	if !ok {
		panic("test bug: unknown FORMAT index")
	}
	// Vector size: max cardinality across samples.
	fieldSize := 0
	for _, sv := range r.sampleVals {
		v, present := sv[keyIdx]
		if !present || v.Missing {
			continue
		}
		n := 0
		switch def.Type {
		case "Integer":
			n = len(v.Ints)
		case "Float":
			n = len(v.Floats)
		case "Character":
			n = 1
		case "String":
			n = len(strings.Join(v.Strs, ","))
		}
		if n > fieldSize {
			fieldSize = n
		}
	}
	var lowType uint8
	switch def.Type {
	case "Integer":
		lowType = typeInt32
	case "Float":
		lowType = typeFloat
	case "Character", "String":
		lowType = typeString
	}
	if fieldSize >= 15 {
		b.WriteByte((15 << 4) | lowType)
		b.WriteByte((1 << 4) | typeInt8)
		b.WriteByte(byte(fieldSize))
	} else {
		b.WriteByte(uint8(fieldSize)<<4 | lowType)
	}
	for _, sv := range r.sampleVals {
		v, present := sv[keyIdx]
		if !present || v.Missing {
			// Missing slot: the per-type sentinel.
			switch def.Type {
			case "Integer":
				var x [4]byte
				binary.LittleEndian.PutUint32(x[:], bcfIntMissingInt32)
				b.Write(x[:])
			case "Float":
				var x [4]byte
				binary.LittleEndian.PutUint32(x[:], bcfFloatMissing)
				b.Write(x[:])
			default:
				b.WriteByte(typeMissing)
			}
			continue
		}
		switch def.Type {
		case "Integer":
			for i := 0; i < fieldSize; i++ {
				var x [4]byte
				if i < len(v.Ints) {
					binary.LittleEndian.PutUint32(x[:], uint32(v.Ints[i]))
				} else {
					binary.LittleEndian.PutUint32(x[:], bcfEndOfVectorInt32)
				}
				b.Write(x[:])
			}
		case "Float":
			for i := 0; i < fieldSize; i++ {
				var x [4]byte
				if i < len(v.Floats) {
					binary.LittleEndian.PutUint32(x[:], math.Float32bits(v.Floats[i]))
				} else {
					binary.LittleEndian.PutUint32(x[:], bcfEndOfVectorFloat32)
				}
				b.Write(x[:])
			}
		case "Character":
			if v.Str != "" {
				b.WriteByte(v.Str[0])
			} else {
				b.WriteByte(typeMissing)
			}
		case "String":
			b.WriteString(strings.Join(v.Strs, ","))
		}
	}
}
