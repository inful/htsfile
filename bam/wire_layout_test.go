package bam

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"strings"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
)

// B-array tags must round-trip the sub-type byte. Before
// v0.20.x the Tag struct had no Sub field, the writer
// emitted the parent TagB letter ('B') as the sub-type
// (which the spec reader rejects), and the reader had no way
// to convey the sub-type back to callers.
func TestBArraySubTypeRoundTrip(t *testing.T) {
	hdr := &Header{
		HD:  &HD{VN: "1.6"},
		SQs: []*SQ{{Name: "chr1", Len: 1000}},
	}
	rec := &Record{
		QName: "r1", RefID: 0, Pos: 1, MapQ: 60, Bin: 4681,
		Cigar: []uint32{PackCigarOp(CigarM, 5)},
		Seq:   "ACGTA", Qual: []byte("IIIII"),
		Tags: []Tag{
			{Key: "XI", Type: TagB, Sub: Tagc, Value: "1,2,3"},
			{Key: "XJ", Type: TagB, Sub: Tagi, Value: "100,-100,32768"},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(rec); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := enc.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	dec := NewDecoder(bytes.NewReader(buf.Bytes()))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	got, err := dec.Record()
	if err != nil {
		t.Fatalf("Record: %v", err)
	}
	if got.Tags[0].Type != TagB {
		t.Errorf("Tags[0].Type = %c, want B", byte(got.Tags[0].Type))
	}
	if got.Tags[0].Sub != Tagc {
		t.Errorf("Tags[0].Sub = %c, want c", byte(got.Tags[0].Sub))
	}
	if got.Tags[0].Value != "1,2,3" {
		t.Errorf("Tags[0].Value = %q, want 1,2,3", got.Tags[0].Value)
	}
	if got.Tags[1].Sub != Tagi {
		t.Errorf("Tags[1].Sub = %c, want i", byte(got.Tags[1].Sub))
	}
}

func TestBArrayInvalidSubTypeRejected(t *testing.T) {
	hdr := &Header{HD: &HD{VN: "1.6"}, SQs: []*SQ{{Name: "chr1", Len: 1000}}}
	rec := &Record{
		QName: "r1", RefID: 0, Pos: 1, MapQ: 60, Bin: 4681,
		Cigar: []uint32{PackCigarOp(CigarM, 5)},
		Seq:   "ACGTA", Qual: []byte("IIIII"),
		Tags: []Tag{
			{Key: "XB", Type: TagB, Sub: TagA, Value: "1"},
		},
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(rec); err == nil {
		t.Errorf("WriteRecord: expected error for invalid B sub-type")
	} else if !strings.Contains(err.Error(), "sub-type") {
		t.Errorf("error %q does not mention sub-type", err)
	}
}

func TestWriter_RejectsOversizedQName(t *testing.T) {
	hdr := &Header{HD: &HD{VN: "1.6"}, SQs: []*SQ{{Name: "chr1", Len: 1000}}}
	rec := &Record{
		QName: strings.Repeat("a", 300),
		RefID: 0, Pos: 1, MapQ: 60, Bin: 4681,
		Cigar: []uint32{PackCigarOp(CigarM, 5)},
		Seq:   "ACGTA", Qual: []byte("IIIII"),
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(rec); err == nil {
		t.Errorf("WriteRecord: expected error for QNAME > 254")
	} else if !strings.Contains(err.Error(), "QNAME") {
		t.Errorf("error %q does not mention QNAME", err)
	}
}

func TestWriter_RejectsOversizedCigarCount(t *testing.T) {
	hdr := &Header{HD: &HD{VN: "1.6"}, SQs: []*SQ{{Name: "chr1", Len: 1000}}}
	rec := &Record{
		QName: "r1", RefID: 0, Pos: 1, MapQ: 60, Bin: 4681,
		Cigar: makeCigar(70000),
		Seq:   "A", Qual: []byte("I"),
	}
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	if err := enc.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := enc.WriteRecord(rec); err == nil {
		t.Errorf("WriteRecord: expected error for n_cigar_op > 65535")
	} else if !strings.Contains(err.Error(), "cigar") && !strings.Contains(err.Error(), "65535") {
		t.Errorf("error %q does not mention cigar/65535", err)
	}
}

func TestParser_BoundsCheckOnTruncatedRemainder(t *testing.T) {
	data, err := makeTruncatedBAM()
	if err != nil {
		t.Fatalf("makeTruncatedBAM: %v", err)
	}
	dec := NewDecoder(bytes.NewReader(data))
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	_, recErr := dec.Record()
	if recErr == nil {
		t.Errorf("expected ParseError on inflated l_seq, got nil")
	} else if !strings.Contains(recErr.Error(), "remainder") &&
		!strings.Contains(recErr.Error(), "l_read_name") {
		t.Errorf("error %q does not mention the bound violation", recErr)
	}
}

func makeCigar(n int) []uint32 {
	out := make([]uint32, n)
	for i := range out {
		out[i] = PackCigarOp(CigarM, 1)
	}
	return out
}

// makeTruncatedBAM returns BGZF-wrapped bytes containing a
// real BAM header followed by a single record whose core
// advertises l_seq=200 but the block carries only a handful
// of payload bytes. The reader must return a ParseError,
// not a panic.
func makeTruncatedBAM() ([]byte, error) {
	var buf bytes.Buffer
	enc := NewEncoder(&buf)
	hdr := &Header{
		HD:  &HD{VN: "1.6"},
		SQs: []*SQ{{Name: "chr1", Len: 1000}},
	}
	if err := enc.WriteHeader(hdr); err != nil {
		return nil, err
	}
	if err := enc.Close(); err != nil {
		return nil, err
	}
	hdrBytes := buf.Bytes()
	buf.Reset()

	// Hand-craft the truncated record.
	var core [32]byte
	copy(core[0:4], bamMagic[:])
	core[8] = 3 // l_read_name
	binary.LittleEndian.PutUint16(core[12:14], 1)
	binary.LittleEndian.PutUint32(core[16:20], 200) // LIE
	qname := []byte{'r', '1', 0}
	var cb [4]byte
	binary.LittleEndian.PutUint32(cb[:], PackCigarOp(CigarM, 5))
	seq := []byte{0x10}
	qual := []byte{'I'}
	inner := append([]byte{}, core[:]...)
	inner = append(inner, qname...)
	inner = append(inner, cb[:]...)
	inner = append(inner, seq...)
	inner = append(inner, qual...)
	bsize := uint32(len(inner)) // 32 + 9 = 41
	var sb [4]byte
	binary.LittleEndian.PutUint32(sb[:], bsize)
	rec := append(sb[:], inner...)

	// Re-emit: write a real BGZF stream that contains the
	// header bytes (the BGZF wrapping the encoder emitted
	// already) followed by a fresh BGZF block carrying the
	// truncated record. Easiest: re-parse the header bytes,
	// copy their first inner block, then add a second block.
	// Even simpler: write both as raw bytes into a single
	// concat of the encoded header stream and a new record
	// stream.
	out, err := appendRecordAfterHeader(hdrBytes, rec)
	if err != nil {
		return nil, err
	}
	return out, nil
}

// appendRecordAfterHeader appends a single raw record
// (already encoded as [bsize][core][rem] bytes) as a new
// BGZF block after the header's BGZF stream. The
// encoder's Close emits a 28-byte BGZF EOF marker at the
// end of its output; we strip that marker so the new block
// is read by the decoder's record loop, not skipped as
// end-of-stream.
func appendRecordAfterHeader(hdrBytes, rec []byte) ([]byte, error) {
	const bgzfEOFMagic = 28 // bytes in the canonical BGZF EOF block
	if len(hdrBytes) < bgzfEOFMagic {
		return nil, fmt.Errorf("header stream too short: %d", len(hdrBytes))
	}
	truncated := hdrBytes[:len(hdrBytes)-bgzfEOFMagic]
	var out bytes.Buffer
	out.Write(truncated)
	bw := bgzf.NewWriter(&out)
	if _, err := bw.Write(rec); err != nil {
		return nil, err
	}
	if err := bw.Close(); err != nil {
		return nil, err
	}
	return out.Bytes(), nil
}
