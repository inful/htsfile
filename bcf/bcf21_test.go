package bcf

import (
	"bytes"
	"io"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
)

// A BCF file with the 2.1 minor version must still be readable.
//
// htslib 1.24 writes 2.2 only and rejects a 2.1 file outright
// ("could not parse header"), and no 2.1 file can be produced
// with the tools available here, so the 2.1 path is exercised
// by patching the minor-version byte of a file we build.
//
// The record layout is the same in both: the spec's list of
// changes between BCFv2.1 and BCFv2.2 (tex lines 2818-2825) is
// about the header's optional IDX field, the END_OF_VECTOR
// byte, and where PASS sits in the dictionary - none of which
// changes the record layout.

// TestRead_BCF21Magic reads a file whose minor version has
// been patched to 1.
func TestRead_BCF21Magic(t *testing.T) {
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	defs := map[string]*InfoDef{
		"DP": {ID: "DP", Number: "1", Type: "Integer"},
	}
	b.writeRecord(t, &testRecord{
		contigIdx: 0, Pos: 100, Ref: "A",
		Info: map[string]InfoValue{
			"DP": {Type: "Integer", Ints: []int64{42}},
		},
		infoDefs:  defs,
		infoOrder: []string{"DP"},
	})
	data := b.finish(t)

	patched := patchMinorVersion(t, data, 1)
	if patched == nil {
		t.Fatal("could not patch the minor version")
	}

	dec := NewDecoder(bytes.NewReader(patched))
	hdr, err := dec.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if dec.Version() != 1 {
		t.Errorf("Version: got %d want 1", dec.Version())
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
	_ = hdr
	if len(recs) != 1 {
		t.Fatalf("records: got %d want 1", len(recs))
	}
	if recs[0].Pos != 100 {
		t.Errorf("Pos: got %d want 100", recs[0].Pos)
	}
}

// TestWrite_AlwaysBCF22 checks that we never write a 2.1 file:
// the minor version byte is 2 unless a caller asks for
// something else.
func TestWrite_AlwaysBCF22(t *testing.T) {
	b := newBuilder(t)
	b.writeBCFHeader(sampleHeader)
	data := b.finish(t)

	raw := decompress(t, data)
	if len(raw) < 5 {
		t.Fatalf("short file: %d bytes", len(raw))
	}
	if string(raw[:5]) != "BCF\x02\x02" {
		t.Errorf("magic: got %q want %q", raw[:5], "BCF\\x02\\x02")
	}
}

// decompress returns the uncompressed contents of a BGZF
// stream.
func decompress(t *testing.T, data []byte) []byte {
	t.Helper()
	raw, err := io.ReadAll(bgzf.NewReader(bytes.NewReader(data)))
	if err != nil {
		t.Fatalf("decompress: %v", err)
	}
	return raw
}

// patchMinorVersion returns the BGZF stream with the BCF
// minor-version byte set to minor.
func patchMinorVersion(t *testing.T, data []byte, minor byte) []byte {
	t.Helper()
	raw := decompress(t, data)
	if len(raw) < 5 || string(raw[:3]) != "BCF" {
		t.Fatalf("not a BCF stream: %q", raw[:min(5, len(raw))])
	}
	raw[4] = minor
	var out bytes.Buffer
	zw := bgzf.NewWriter(&out)
	if _, err := zw.Write(raw); err != nil {
		t.Fatalf("write: %v", err)
	}
	if err := zw.Close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	return out.Bytes()
}
