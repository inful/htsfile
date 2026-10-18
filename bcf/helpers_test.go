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
	path        string
	f           *os.File
	hdr         *Header // parsed sampleHeader; used for LookupInfoIndex etc.
	hdrOrder    []string // cached: hdr.InfosOrder
	hdrFilters  []string
	hdrContigs  []string
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
	var payload bytes.Buffer
	var prefix [16]byte
	binary.LittleEndian.PutUint32(prefix[0:4], uint32(r.contigIdx))
	binary.LittleEndian.PutUint32(prefix[4:8], uint32(int32(r.Pos)-1))
	binary.LittleEndian.PutUint32(prefix[8:12], 0)
	binary.LittleEndian.PutUint32(prefix[12:16], uint32(len(r.Info)))
	payload.Write(prefix[:])
	var nFmt [2]byte
	binary.LittleEndian.PutUint16(nFmt[:], 0)
	payload.Write(nFmt[:])
	writeTypedString(&payload, r.ID)
	writeTypedString(&payload, r.Ref)
	for _, alt := range r.Alt {
		writeTypedString(&payload, alt)
	}
	writeTypedFilter(&payload, r.filterIdxs)
	if r.QualKnown {
		payload.WriteByte(typeFloat)
		var bits [8]byte
		binary.LittleEndian.PutUint64(bits[:], float64Bits(r.Qual))
		payload.Write(bits[:])
	} else {
		payload.WriteByte(typeMissing)
	}
	for _, id := range r.infoOrder {
		def := r.infoDefs[id]
		val := r.Info[id]
		// Look up the key index from the header's
		// InfosOrder via the package's LookupInfoIndex
		// helper. BCF stores the key as a positional index
		// into the header's ##INFO declaration list; the
		// test's local infoOrder does NOT match the
		// header's declaration order.
		hdrIdx, ok := b.hdr.LookupInfoIndex(id)
		if !ok {
			panic(fmt.Sprintf("test bug: %q not in header InfosOrder (have %v)", id, b.hdrOrder))
		}
		var keyBytes [2]byte
		binary.LittleEndian.PutUint16(keyBytes[:], uint16(hdrIdx))
		payload.Write(keyBytes[:])
		writeTypedInfoValue(&payload, def, val)
	}
	b.writeBGZFBlock(payload.Bytes())
}

func writeTypedString(b *bytes.Buffer, s string) {
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

func writeTypedFilter(b *bytes.Buffer, idxs []int) {
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

func writeTypedInfoValue(b *bytes.Buffer, def *InfoDef, v InfoValue) {
	if def.Type == "Flag" {
		if v.Flag {
			b.WriteByte(typeInt8)
		} else {
			b.WriteByte(typeMissing)
		}
		return
	}
	count, isVar := ParseNumber(def.Number)
	// For variable-count types (A/R/G/.) or unset counts,
	// determine count from the value's actual length.
	if count == 0 {
		switch def.Type {
		case "Integer":
			count = len(v.Ints)
		case "Float":
			count = len(v.Floats)
		case "String":
			count = len(v.Strs)
		case "Character":
			count = 1
		}
	}
	if isVar {
		b.WriteByte(typeInt32)
		var n [4]byte
		binary.LittleEndian.PutUint32(n[:], uint32(count))
		b.Write(n[:])
	}
	switch def.Type {
	case "Integer":
		b.WriteByte(typeInt32)
		for _, n := range v.Ints {
			var x [4]byte
			binary.LittleEndian.PutUint32(x[:], uint32(int32(n)))
			b.Write(x[:])
		}
	case "Float":
		b.WriteByte(typeFloat)
		for _, f := range v.Floats {
			var x [4]byte
			binary.LittleEndian.PutUint32(x[:], float32Bits(float32(f)))
			b.Write(x[:])
		}
	case "String":
		b.WriteByte(typeString)
		for _, s := range v.Strs {
			var l [4]byte
			binary.LittleEndian.PutUint32(l[:], uint32(len(s)))
			b.Write(l[:])
			b.WriteString(s)
		}
	case "Character":
		b.WriteByte(typeString)
		var l [4]byte
		binary.LittleEndian.PutUint32(l[:], 1)
		b.Write(l[:])
		if len(v.Str) > 0 {
			b.WriteByte(v.Str[0])
		}
	}
}

// float64Bits / float32Bits wrap math.Float64bits /
// Float32bits so the helpers don't need a separate math
// import.
func float64Bits(f float64) uint64 { return math.Float64bits(f) }
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
