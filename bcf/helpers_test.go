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
	var l [4]byte
	binary.LittleEndian.PutUint32(l[:], uint32(len(text)))
	payload.Write(l[:])
	payload.WriteString(text)
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
}

// writeRecord encodes one BCF record into a fresh BGZF
// block.
func (b *builder) writeRecord(t *testing.T, r *testRecord) {
	t.Helper()
	var shared bytes.Buffer
	// Spec section 5, BcfSiteEncoding: 24-byte site prefix
	// CHROM, POS, rlen, QUAL, n_info, n_allele, n_sample,
	// n_fmt.
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
	// n_sample=0 and n_fmt=0 in the shared prefix; the
	// per-sample block is empty for these fixtures.
	shared.Write(prefix[:])

	testWriteTypedString(&shared, r.ID)
	testWriteTypedString(&shared, r.Ref)
	for _, alt := range r.Alt {
		testWriteTypedString(&shared, alt)
	}
	testWriteTypedFilter(&shared, r.filterIdxs)
	for _, id := range r.infoOrder {
		def := r.infoDefs[id]
		val := r.Info[id]
		hdrIdx, ok := b.hdr.LookupInfoIndex(id)
		if !ok {
			panic(fmt.Sprintf("test bug: %q not in header InfosOrder (have %v)", id, b.hdrOrder))
		}
		// INFO key is a typed atomic int, not a bare uint16.
		shared.WriteByte(typeByteAtomic(typeInt32))
		var kb [4]byte
		binary.LittleEndian.PutUint32(kb[:], uint32(hdrIdx))
		shared.Write(kb[:])
		testWriteTypedInfoValue(&shared, def, val)
	}

	// l_shared + l_indiv framing.
	var lens [8]byte
	binary.LittleEndian.PutUint32(lens[0:4], uint32(shared.Len()))
	binary.LittleEndian.PutUint32(lens[4:8], 0)
	var rec bytes.Buffer
	rec.Write(lens[:])
	rec.Write(shared.Bytes())
	b.writeBGZFBlock(rec.Bytes())
}

func testWriteTypedString(b *bytes.Buffer, s string) {
	if s == "" {
		b.WriteByte(typeMissing)
		return
	}
	b.WriteByte(typeString)
	var l [4]byte
	binary.LittleEndian.PutUint32(l[:], uint32(len(s)))
	b.Write(l[:])
	b.WriteString(s)
}

func testWriteTypedFilter(b *bytes.Buffer, idxs []int) {
	if len(idxs) == 0 {
		b.WriteByte(typeMissing)
		return
	}
	if len(idxs) == 1 {
		b.WriteByte(typeInt8)
		var idx [4]byte
		binary.LittleEndian.PutUint32(idx[:], uint32(idxs[0]))
		b.Write(idx[:])
		return
	}
	b.WriteByte(typeInt32)
	var n [4]byte
	binary.LittleEndian.PutUint32(n[:], uint32(len(idxs)))
	b.Write(n[:])
	for _, i := range idxs {
		var idx [4]byte
		binary.LittleEndian.PutUint32(idx[:], uint32(i))
		b.Write(idx[:])
	}
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
		count, hasCount = len(v.Strs), true
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
		for _, st := range v.Strs {
			var l [4]byte
			binary.LittleEndian.PutUint32(l[:], uint32(len(st)))
			b.Write(l[:])
			b.WriteString(st)
		}
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
