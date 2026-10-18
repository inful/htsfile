package bam

import (
	"bytes"
	"compress/flate"
	"encoding/binary"
	"errors"
	"hash/crc32"
	"io"
	"strconv"
	"testing"
)

// builder is a small in-test helper for constructing synthetic
// BAM streams. It writes into a bgzf-wrapped buffer so the
// resulting stream is byte-for-byte a valid BAM file.
type builder struct {
	buf *bytes.Buffer
	bw  *bgzfWriter
}

// bgzfWriter wraps bgzf to expose a Write([]byte) interface
// for tests. We deliberately do not export a real bgzf writer
// yet — the v1 BAM package is read-only and the test helper is
// the only place we need one.
type bgzfWriter struct {
	br io.Writer
}

func newBuilder() *builder {
	var buf bytes.Buffer
	return &builder{buf: &buf, bw: &bgzfWriter{br: &buf}}
}

// Write a single BGZF block containing payload to the underlying writer.
func (w *bgzfWriter) Write(payload []byte) (int, error) {
	var cdata bytes.Buffer
	fw, err := flate.NewWriter(&cdata, flate.DefaultCompression)
	if err != nil {
		return 0, err
	}
	if _, err := fw.Write(payload); err != nil {
		return 0, err
	}
	if err := fw.Close(); err != nil {
		return 0, err
	}
	cd := cdata.Bytes()
	bsize := uint16(26 + len(cd) - 1)
	header := []byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00,
	}
	var bs [2]byte
	binary.LittleEndian.PutUint16(bs[:], bsize)
	header = append(header, bs[:]...)
	header = append(header, cd...)
	var crc [4]byte
	binary.LittleEndian.PutUint32(crc[:], crc32.ChecksumIEEE(payload))
	header = append(header, crc[:]...)
	var is [4]byte
	binary.LittleEndian.PutUint32(is[:], uint32(len(payload)))
	header = append(header, is[:]...)
	_, err = w.br.Write(header)
	if err != nil {
		return 0, err
	}
	return len(payload), nil
}

// writeEOFMarker writes the 28-byte empty BGZF EOF block.
func (w *bgzfWriter) writeEOFMarker() error {
	eof := []byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00,
	}
	_, err := w.br.Write(eof)
	return err
}

// writeBAMHeader writes the BAM magic, SAM header text, and
// reference dictionary to the BGZF-wrapped stream.
func (b *builder) writeBAMHeader(samText string, refs []struct {
	name   string
	length int64
}) error {
	var payload bytes.Buffer
	payload.Write(bamMagic[:])
	var lt [4]byte
	binary.LittleEndian.PutUint32(lt[:], uint32(len(samText)))
	payload.Write(lt[:])
	payload.WriteString(samText)
	var nr [4]byte
	binary.LittleEndian.PutUint32(nr[:], uint32(len(refs)))
	payload.Write(nr[:])
	for _, ref := range refs {
		name := ref.name + "\x00" // null-terminated
		var ln [4]byte
		binary.LittleEndian.PutUint32(ln[:], uint32(len(name)))
		payload.Write(ln[:])
		payload.WriteString(name)
		var ll [4]byte
		binary.LittleEndian.PutUint32(ll[:], uint32(ref.length))
		payload.Write(ll[:])
	}
	if _, err := b.bw.Write(payload.Bytes()); err != nil {
		return err
	}
	return nil
}

// writeRecord appends one BAM alignment record (binary form) to
// the stream.
func (b *builder) writeRecord(r *Record) error {
	// Build the variable-length remainder first.
	var rem bytes.Buffer
	// read_name: include trailing NUL.
	rem.WriteString(r.QName)
	rem.WriteByte(0)
	// cigar: each op is uint32 LE.
	for _, op := range r.Cigar {
		var c [4]byte
		binary.LittleEndian.PutUint32(c[:], op)
		rem.Write(c[:])
	}
	// seq: 4-bit packed. We re-encode r.Seq into packed bytes.
	lSeq := int32(0)
	if r.Seq != "*" {
		lSeq = int32(len(r.Seq))
	}
	if lSeq > 0 {
		packed := make([]byte, (int(lSeq)+1)/2)
		for i := 0; i < int(lSeq); i++ {
			code := testBaseCode(r.Seq[i])
			if i%2 == 0 {
				packed[i/2] = code << 4
			} else {
				packed[i/2] |= code & 0x0F
			}
		}
		rem.Write(packed)
	}
	// qual: ASCII bytes, length = lSeq.
	if lSeq > 0 {
		if len(r.Qual) == int(lSeq) {
			rem.Write(r.Qual)
		} else {
			// pad with 'I' (Phred+33 = 40 = 'I') if the test
			// didn't supply qual bytes.
			rem.Write(bytes.Repeat([]byte{'I'}, int(lSeq)))
		}
	}
	// aux tags.
	aux := testEncodeAuxTags(r.Tags)
	rem.Write(aux)
	// Compute block_size = 32 (core) + len(rem).
	core := make([]byte, 32)
	binary.LittleEndian.PutUint32(core[0:4], uint32(r.RefID))
	pos0 := uint32(int32(r.Pos) - 1)
	if r.Pos <= 0 {
		pos0 = 0xFFFFFFFF // -1 in int32
	}
	binary.LittleEndian.PutUint32(core[4:8], pos0)
	core[8] = uint8(len(r.QName)) + 1 // l_read_name including NUL
	core[9] = r.MapQ
	binary.LittleEndian.PutUint16(core[10:12], r.Bin)
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
	blockSize := uint32(32 + rem.Len())
	var bs [4]byte
	binary.LittleEndian.PutUint32(bs[:], blockSize)
	var payload bytes.Buffer
	payload.Write(bs[:])
	payload.Write(core)
	payload.Write(rem.Bytes())
	_, err := b.bw.Write(payload.Bytes())
	return err
}

// Test-only helpers for tag encoding (testEncodeAuxTags,
// testEncodeTagValue, testBaseCode, testSplitComma) live in
// parse_helpers_test.go.

// small numeric parsers for tag values.
func parseInt64(s string) int64 {
	n := int64(0)
	neg := false
	i := 0
	if len(s) > 0 && s[0] == '-' {
		neg = true
		i = 1
	}
	for ; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			break
		}
		n = n*10 + int64(s[i]-'0')
	}
	if neg {
		n = -n
	}
	return n
}

func parseUint64(s string) uint64 {
	n := uint64(0)
	for i := 0; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			break
		}
		n = n*10 + uint64(s[i]-'0')
	}
	return n
}

func parseFloat(s string, out *float32) {
	v, err := strconv.ParseFloat(s, 32)
	if err != nil {
		*out = 0
		return
	}
	*out = float32(v)
}

// finish writes the EOF marker and returns the assembled stream.
func (b *builder) finish() []byte {
	_ = b.bw.writeEOFMarker()
	return b.buf.Bytes()
}

// cigar packs a BAM CIGAR op as (op << 28) | length. The op is
// in the high 4 bits and the length is in the low 28 bits per
// the BAM spec § 1.4.6.
func cigar(op, length uint32) uint32 {
	return (op << 28) | length
}

// readAll reads all records from a BAM byte slice.
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
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		recs = append(recs, rec)
	}
	return hdr, recs
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

func TestHeader_Magic(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\tSO:coordinate\n", nil); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	_, err := NewDecoder(bytes.NewReader(data)).Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
}

func TestHeader_BadMagic(t *testing.T) {
	data := []byte{'B', 'A', 'M', 0x02, 0, 0, 0, 0}
	_, err := NewDecoder(bytes.NewReader(data)).Header()
	if err == nil {
		t.Fatal("expected error on bad magic, got nil")
	}
}

func TestHeader_EmbeddedSAMText(t *testing.T) {
	samText := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:248956422\n" +
		"@RG\tID:rg1\tSM:sample1\tLB:lib1\n" +
		"@PG\tID:prog1\tPN:htsfile\n" +
		"@CO\tcomment\n"
	b := newBuilder()
	if err := b.writeBAMHeader(samText, []struct {
		name   string
		length int64
	}{{"chr1", 248956422}}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	hdr, err := NewDecoder(bytes.NewReader(data)).Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr.HD == nil {
		t.Fatal("HD missing")
	}
	if hdr.HD.VN != "1.6" || hdr.HD.SO != "coordinate" {
		t.Errorf("HD fields: VN=%q SO=%q", hdr.HD.VN, hdr.HD.SO)
	}
	if len(hdr.SQs) != 1 || hdr.SQs[0].Name != "chr1" || hdr.SQs[0].Len != 248956422 {
		t.Errorf("SQ mismatch: %+v", hdr.SQs)
	}
	if len(hdr.RGs) != 1 || hdr.RGs[0].ID != "rg1" || hdr.RGs[0].SM != "sample1" {
		t.Errorf("RG mismatch: %+v", hdr.RGs)
	}
	if len(hdr.PGs) != 1 || hdr.PGs[0].ID != "prog1" || hdr.PGs[0].PN != "htsfile" {
		t.Errorf("PG mismatch: %+v", hdr.PGs)
	}
	if len(hdr.COs) != 1 {
		t.Errorf("CO mismatch: %+v", hdr.COs)
	}
	if hdr.Raw != samText {
		t.Errorf("Raw mismatch")
	}
}

func TestRecord_BasicFields(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000}}); err != nil {
		t.Fatal(err)
	}
	rec := &Record{
		QName:     "read1",
		Flag:      0,
		RefID:     0,
		Pos:       100,
		MapQ:      60,
		Bin:       4681,
		Cigar:     []uint32{cigar(CigarM, 5)}, // 5M
		NextRefID: -1,
		NextPos:   0,
		TLen:      0,
		Seq:       "ACGTA",
		Qual:      []byte("IIIII"),
		Tags:      nil,
	}
	if err := b.writeRecord(rec); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	_, recs := readAll(t, data)
	if len(recs) != 1 {
		t.Fatalf("expected 1 record, got %d", len(recs))
	}
	got := recs[0]
	if got.QName != "read1" {
		t.Errorf("QName: %q", got.QName)
	}
	if got.Pos != 100 {
		t.Errorf("Pos: %d", got.Pos)
	}
	if got.MapQ != 60 {
		t.Errorf("MapQ: %d", got.MapQ)
	}
	if got.Seq != "ACGTA" {
		t.Errorf("Seq: %q", got.Seq)
	}
	if string(got.Qual) != "IIIII" {
		t.Errorf("Qual: %q", got.Qual)
	}
	if got.CigarString() != "5M" {
		t.Errorf("CigarString: %q", got.CigarString())
	}
	if got.NameRName != "chr1" {
		t.Errorf("NameRName: %q", got.NameRName)
	}
}

func TestRecord_Unmapped(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n", nil); err != nil {
		t.Fatal(err)
	}
	rec := &Record{
		QName:     "unmapped",
		Flag:      FlagUnmap,
		RefID:     -1,
		Pos:       0,
		MapQ:      0,
		Bin:       4680,
		Seq:       "*",
		Qual:      nil,
		NextRefID: -1,
		NextPos:   0,
		TLen:      0,
	}
	if err := b.writeRecord(rec); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	_, recs := readAll(t, data)
	if len(recs) != 1 {
		t.Fatalf("expected 1 record, got %d", len(recs))
	}
	if recs[0].Flag&FlagUnmap == 0 {
		t.Errorf("FlagUnmap bit not set")
	}
	if recs[0].RefID != -1 {
		t.Errorf("RefID: %d", recs[0].RefID)
	}
	if recs[0].Seq != "*" {
		t.Errorf("Seq: %q", recs[0].Seq)
	}
	if recs[0].NameRName != "*" {
		t.Errorf("NameRName: %q (want *)", recs[0].NameRName)
	}
}

func TestRecord_PairedRead(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n", []struct {
		name   string
		length int64
	}{{"chr1", 100000}}); err != nil {
		t.Fatal(err)
	}
	// First mate.
	if err := b.writeRecord(&Record{
		QName: "pair1",
		Flag:  FlagPaired | FlagProperPair | FlagFirst,
		RefID: 0, Pos: 100,
		MapQ:      60,
		Bin:       4681,
		Cigar:     []uint32{cigar(CigarM, 5)}, // 5M
		NextRefID: 0, NextPos: 200,
		TLen: 100,
		Seq:  "ACGTA",
		Qual: []byte("IIIII"),
	}); err != nil {
		t.Fatal(err)
	}
	// Second mate.
	if err := b.writeRecord(&Record{
		QName: "pair1",
		Flag:  FlagPaired | FlagProperPair | FlagLast,
		RefID: 0, Pos: 200,
		MapQ:      60,
		Bin:       4681,
		Cigar:     []uint32{cigar(CigarM, 5)}, // 5M
		NextRefID: 0, NextPos: 100,
		TLen: -100,
		Seq:  "TCGTA",
		Qual: []byte("IIIII"),
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	_, recs := readAll(t, data)
	if len(recs) != 2 {
		t.Fatalf("expected 2 records, got %d", len(recs))
	}
	if recs[0].Flag&FlagFirst == 0 {
		t.Errorf("first mate missing First flag")
	}
	if recs[1].Flag&FlagLast == 0 {
		t.Errorf("second mate missing Last flag")
	}
	if recs[0].TLen != 100 {
		t.Errorf("first mate TLen: %d", recs[0].TLen)
	}
	if recs[1].TLen != -100 {
		t.Errorf("second mate TLen: %d", recs[1].TLen)
	}
}

func TestRecord_LongSequence(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000000}}); err != nil {
		t.Fatal(err)
	}
	// 200-base sequence.
	seq := ""
	for i := 0; i < 200; i++ {
		seq += string("ACGT"[i%4])
	}
	qual := bytes.Repeat([]byte{'I'}, 200)
	cigar := []uint32{cigar(CigarM, 200)} // 200M
	if err := b.writeRecord(&Record{
		QName: "long",
		RefID: 0, Pos: 500,
		MapQ: 60, Bin: 4681,
		Cigar:     cigar,
		Seq:       seq,
		Qual:      qual,
		NextRefID: -1,
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	_, recs := readAll(t, data)
	if len(recs) != 1 {
		t.Fatalf("expected 1 record, got %d", len(recs))
	}
	if recs[0].Seq != seq {
		t.Errorf("Seq round-trip failed (len=%d want=%d)", len(recs[0].Seq), len(seq))
	}
	if len(recs[0].Qual) != 200 {
		t.Errorf("Qual length: %d", len(recs[0].Qual))
	}
}

func TestRecord_Tags(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000}}); err != nil {
		t.Fatal(err)
	}
	tags := []Tag{
		{Key: "NM", Type: Tagi, Value: "1"},
		{Key: "MD", Type: TagZ, Value: "5"},
		{Key: "XA", Type: TagA, Value: "A"},
		{Key: "XB", Type: Tagf, Value: "1.5"},
		{Key: "XC", Type: TagH, Value: "1A2B"},
		{Key: "XD", Type: Tagd, Value: "3.14"},
		{Key: "XE", Type: Tagc, Value: "-1"},
		{Key: "XF", Type: TagC, Value: "200"},
		{Key: "XG", Type: Tags, Value: "-100"},
		{Key: "XH", Type: TagS, Value: "500"},
		{Key: "XI", Type: TagI, Value: "1234567"},
	}
	if err := b.writeRecord(&Record{
		QName: "tagged", RefID: 0, Pos: 100, MapQ: 60, Bin: 4681,
		Cigar: []uint32{cigar(CigarM, 5)},
		Seq:   "ACGTA", Qual: []byte("IIIII"),
		Tags: tags,
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	_, recs := readAll(t, data)
	if len(recs) != 1 {
		t.Fatalf("expected 1 record, got %d", len(recs))
	}
	got := recs[0].Tags
	if len(got) != len(tags) {
		t.Fatalf("tag count: got %d want %d", len(got), len(tags))
	}
	for i, tg := range tags {
		if got[i].Key != tg.Key {
			t.Errorf("tag[%d].Key: %q want %q", i, got[i].Key, tg.Key)
		}
		if got[i].Type != tg.Type {
			t.Errorf("tag[%d].Type: %q want %q", i, got[i].Type, tg.Type)
		}
	}
	// Spot-check a couple of values.
	if got[0].Value != "1" {
		t.Errorf("NM value: %q want 1", got[0].Value)
	}
	if got[1].Value != "5" {
		t.Errorf("MD value: %q want 5", got[1].Value)
	}
	if got[2].Value != "A" {
		t.Errorf("XA value: %q want A", got[2].Value)
	}
	if got[3].Value != "1.5" {
		t.Errorf("XB value: %q want 1.5", got[3].Value)
	}
	if got[7].Value != "200" {
		t.Errorf("XF value: %q want 200", got[7].Value)
	}
}

func TestValidator_CleanRecords(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000}}); err != nil {
		t.Fatal(err)
	}
	if err := b.writeRecord(&Record{
		QName: "clean", RefID: 0, Pos: 100, MapQ: 60, Bin: 4681,
		Cigar: []uint32{cigar(CigarM, 5)}, Seq: "ACGTA", Qual: []byte("IIIII"),
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	hdr, recs := readAll(t, data)
	val := NewValidator(ModeDefault)
	issues := val.Validate(hdr, recs)
	if len(issues) != 0 {
		t.Errorf("expected no issues, got %v", issues)
	}
}

func TestValidator_BadRefID(t *testing.T) {
	b := newBuilder()
	if err := b.writeBAMHeader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n", []struct {
		name   string
		length int64
	}{{"chr1", 1000}}); err != nil {
		t.Fatal(err)
	}
	if err := b.writeRecord(&Record{
		QName: "bad", RefID: 5, Pos: 100, MapQ: 60, Bin: 4681,
		Cigar: []uint32{cigar(CigarM, 5)}, Seq: "ACGTA", Qual: []byte("IIIII"),
	}); err != nil {
		t.Fatal(err)
	}
	data := b.finish()
	hdr, recs := readAll(t, data)
	val := NewValidator(ModeDefault)
	issues := val.Validate(hdr, recs)
	if len(issues) == 0 {
		t.Errorf("expected refID error")
	}
}

func TestCigarString(t *testing.T) {
	tests := []struct {
		cigar []uint32
		want  string
	}{
		{nil, "*"},
		{[]uint32{}, "*"},
		{[]uint32{cigar(CigarM, 5)}, "5M"},
		{[]uint32{cigar(CigarM, 5), cigar(CigarI, 2), cigar(CigarM, 3)}, "5M2I3M"},
		{[]uint32{cigar(CigarS, 3), cigar(CigarI, 4), cigar(CigarD, 5)}, "3S4I5D"},
	}
	for _, tc := range tests {
		got := (&Record{Cigar: tc.cigar}).CigarString()
		if got != tc.want {
			t.Errorf("CigarString(%v): got %q want %q", tc.cigar, got, tc.want)
		}
	}
}
