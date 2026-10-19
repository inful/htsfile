package vcf

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/inful/htsfile/internal/bgzf"
	"github.com/inful/htsfile/tabix"
)

// makeIndexedVCF writes a BGZF-compressed VCF file and a
// tabix index that covers the entire file in one chunk.
// Returns the VCF path.
func makeIndexedVCF(t *testing.T) (vcfPath string) {
	t.Helper()
	dir := t.TempDir()
	vcfPath = filepath.Join(dir, "sample.vcf.gz")
	tbiPath := vcfPath + ".tbi"
	_ = tbiPath

	// Step 1: BGZF-compress the VCF content.
	const content = "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=10000>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t.\t.\t.\n" +
		"chr1\t300\t.\tC\tT\t.\t.\t.\n" +
		"chr1\t500\t.\tG\tA\t.\t.\t.\n"

	f, err := os.Create(vcfPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("create vcf: %v", err)
	}
	bw := bgzf.NewWriter(f)
	if _, err := bw.Write([]byte(content)); err != nil {
		t.Fatalf("bw.Write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bw.Close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("f.Close: %v", err)
	}

	// Step 2: compute the BGZF virtual offset for each
	// line. Tabix chunk Beg/End are virtual offsets.
	// For a one-chunk-per-record index, we use a single
	// chunk that covers the whole file (offset 0 to EOF).
	//
	// The EOF marker of a BGZF file is at virtual offset
	// (block_count << 16) | 0. For a file with one block,
	// that's 0x10000 | 0 = 65536.
	idx := &tabix.Index{
		References: []*tabix.Reference{
			{
				Bins: map[uint32][]tabix.Chunk{
					0: {{Beg: 0, End: 0x10000}},
				},
				Intervals: []tabix.VirtualOffset{0, 0x10000},
			},
		},
	}

	// Step 3: write the index file.
	tbiFile, err := os.Create(tbiPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("create tbi: %v", err)
	}
	if err := tabix.NewWriter(tbiFile).Write(idx); err != nil {
		_ = tbiFile.Close()
		t.Fatalf("write tbi: %v", err)
	}
	if err := tbiFile.Close(); err != nil {
		t.Fatalf("close tbi: %v", err)
	}

	return vcfPath
}

func TestTabixReader_OpenAndHeader(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	hdr, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if len(hdr.Contigs) != 1 || hdr.Contigs[0].ID != "chr1" {
		t.Errorf("Contigs: %+v, want one chr1", hdr.Contigs)
	}
}

func TestTabixReader_QueryRegion(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	// Query region that overlaps pos 100 and 300 but not
	// 500. With one chunk covering the whole file, we'll
	// fetch all three records and filter by POS.
	recs, err := r.Query("chr1", 200, 400)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 1 {
		t.Fatalf("Query: got %d records, want 1", len(recs))
	}
	if recs[0].Pos != 300 {
		t.Errorf("Query: got POS %d, want 300", recs[0].Pos)
	}
}

func TestTabixReader_QueryWholeContig(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	recs, err := r.Query("chr1", 1, 10000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 3 {
		t.Errorf("Query: got %d records, want 3", len(recs))
	}
}

func TestTabixReader_QueryUnknownContig(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	recs, err := r.Query("chrX", 1, 100)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	if len(recs) != 0 {
		t.Errorf("Query: got %d records, want 0", len(recs))
	}
}

func TestTabixReader_QueryInvalidRegion(t *testing.T) {
	vcfPath := makeIndexedVCF(t)
	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Query("chr1", 0, 100); err == nil {
		t.Error("Query(0, 100): expected error for 0-based start")
	}
	if _, err := r.Query("chr1", 100, 100); err == nil {
		t.Error("Query(100, 100): expected error for end <= start")
	}
}

func TestTabixReader_MissingIndex(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "no-index.vcf")
	if err := os.WriteFile(path, []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"), 0o600); err != nil {
		t.Fatalf("write: %v", err)
	}
	if _, err := OpenTabix(path); err == nil {
		t.Error("OpenTabix: expected error for missing .tbi")
	}
}

func TestTabixReader_ViaBGZFStream(t *testing.T) {
	// Confirm the BGZF stream round-trips through the
	// internal bgzf.Reader. The first two bytes must be
	// the gzip magic (1F 8B).
	vcfPath := makeIndexedVCF(t)
	f, err := os.Open(vcfPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer func() { _ = f.Close() }()
	hdr := make([]byte, 2)
	if _, err := f.Read(hdr); err != nil {
		t.Fatalf("read: %v", err)
	}
	if hdr[0] != 0x1f || hdr[1] != 0x8b {
		t.Fatalf("not BGZF: got %x %x", hdr[0], hdr[1])
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		t.Fatalf("seek: %v", err)
	}
	br := bgzf.NewReader(f)
	var buf bytes.Buffer
	if _, err := buf.ReadFrom(br); err != nil {
		t.Fatalf("bgzf read: %v", err)
	}
	if !bytes.Contains(buf.Bytes(), []byte("chr1\t100")) {
		t.Errorf("BGZF content missing chr1:100: %q", buf.String())
	}
}

// makeMultiBlockVCF writes a BGZF-compressed VCF with
// enough records that the decompressed payload spans
// multiple BGZF blocks. Each record carries a long INFO
// string of incompressible random bases so 200 records
// × ~700 bytes (>140 KiB) crosses the 64 KiB block
// boundary. Returns the VCF path.
func makeMultiBlockVCF(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "big.vcf.gz")

	// Random-base INFO so the deflate stream cannot
	// collapse the payload. The bgzf.Writer forces a
	// block split when its buffer reaches 65535 bytes
	// (maxBGZFBlockSize), so uncompressible > 70 KiB
	// yields multiple blocks deterministically.
	pad := randomBases(600)
	var buf bytes.Buffer
	buf.WriteString("##fileformat=VCFv4.2\n")
	buf.WriteString("##contig=<ID=chr1,length=200000>\n")
	buf.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for i := 1; i <= 200; i++ {
		fmt.Fprintf(&buf, "chr1\t%d\t.\tA\tG\t.\t.\tP=%s\n", i*100, pad)
	}

	f, err := os.Create(vcfPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("create vcf: %v", err)
	}
	bw := bgzf.NewWriter(f)
	if _, err := bw.Write(buf.Bytes()); err != nil {
		t.Fatalf("bw.Write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bw.Close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("f.Close: %v", err)
	}
	return vcfPath
}

// randomBases returns n bytes chosen from {A,C,G,T,N}
// using a deterministic LCG so test fixtures are
// reproducible. The output is incompressible for deflate,
// which forces BGZF block boundaries.
func randomBases(n int) string {
	const alphabet = "ACGTN"
	b := make([]byte, n)
	state := uint64(0xdeadbeefcafe1234)
	for i := 0; i < n; i++ {
		state = state*6364136223846793005 + 1442695040888963407
		b[i] = alphabet[state%5]
	}
	return string(b)
}

// TestTabixReader_ChunkEndStop verifies that decodeChunk
// stops reading when the BGZF virtual offset reaches the
// chunk's End, instead of decoding until EOF. The fixture
// spans multiple BGZF blocks; we discover the file
// offset of the second block and build an index whose
// chunk End falls exactly at that block boundary. The
// Query must not return records that lie in the second
// block.
func TestTabixReader_ChunkEndStop(t *testing.T) {
	vcfPath := makeMultiBlockVCF(t)

	// OpenTabix requires a .tbi. Write a placeholder
	// index so we can open the file; we overwrite it
	// below with the chunk-boundary-aware index.
	tbiPath := vcfPath + ".tbi"
	tbiFile, err := os.Create(tbiPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("placeholder tbi: %v", err)
	}
	placeholder := &tabix.Index{
		References: []*tabix.Reference{{
			Bins:      map[uint32][]tabix.Chunk{0: {{Beg: 0, End: tabix.VirtualOffset(0xFFFFFFFFFFFFFFFF)}}},
			Intervals: []tabix.VirtualOffset{0, tabix.VirtualOffset(0xFFFFFFFFFFFFFFFF)},
		}},
	}
	if err := tabix.NewWriter(tbiFile).Write(placeholder); err != nil {
		_ = tbiFile.Close()
		t.Fatalf("placeholder write: %v", err)
	}
	_ = tbiFile.Close()

	r, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix: %v", err)
	}
	defer func() { _ = r.Close() }()

	if err := r.seek.Seek(0); err != nil {
		t.Fatalf("seek 0: %v", err)
	}

	// Walk the BGZF stream block-by-block and record
	// each block's start (compressed file offset).
	// Read records along the way; the first record of
	// each block is the "boundary record" whose offset
	// marks the chunk End we will set below.
	type boundary struct {
		blockOff int64
		firstPos int64
	}
	var boundaries []boundary
	dec := NewDecoder(r.seek)
	if _, err := dec.Header(); err != nil {
		t.Fatalf("Header: %v", err)
	}
	var prevBlockOff int64
	for {
		// Capture the BGZF virtual offset BEFORE
		// reading the next record. If its high 48
		// bits (block offset) differ from the
		// previous iteration, this record starts in
		// a new block — and that block's start
		// offset is the chunk End we want.
		prevOff := r.seek.VirtualOffset()
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Record: %v", err)
		}
		curOff := r.seek.VirtualOffset()
		blockOff := int64(prevOff >> 16)
		if blockOff != prevBlockOff {
			boundaries = append(boundaries, boundary{
				blockOff: blockOff,
				firstPos: rec.Pos,
			})
			prevBlockOff = blockOff
		}
		_ = curOff
	}
	if len(boundaries) < 1 {
		t.Skip("fixture did not span multiple BGZF blocks")
	}

	// Use the second-block file offset (after the
	// boundary) as chunk End. The chunk End is a BGZF
	// virtual offset; we set the low 16 bits to 0
	// because the chunk boundary is at the block
	// start.
	block2End := uint64(boundaries[0].blockOff&0x0000FFFFFFFFFFFF) << 16

	// Build the index: a single chunk that spans from
	// the start of the file to just before block 2.
	idx := &tabix.Index{
		References: []*tabix.Reference{
			{
				Bins: map[uint32][]tabix.Chunk{
					0: {{Beg: 0, End: tabix.VirtualOffset(block2End)}},
				},
				Intervals: []tabix.VirtualOffset{0, tabix.VirtualOffset(block2End)},
			},
		},
	}
	tbiFile, err = os.Create(tbiPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("recreate tbi: %v", err)
	}
	if err := tabix.NewWriter(tbiFile).Write(idx); err != nil {
		_ = tbiFile.Close()
		t.Fatalf("write tbi: %v", err)
	}
	_ = tbiFile.Close()

	// Re-open and query. With chunk End at the block-2
	// boundary, decodeChunk must not read into block 2.
	// We expect the same number of records we observed
	// in block 1.
	r2, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix #2: %v", err)
	}
	defer func() { _ = r2.Close() }()
	if _, err := r2.Header(); err != nil {
		t.Fatalf("Header #2: %v", err)
	}
	recs, err := r2.Query("chr1", 1, 200000)
	if err != nil {
		t.Fatalf("Query (chunked): %v", err)
	}
	if len(recs) >= 200 {
		t.Fatalf("Query (chunked): got %d records, want < 200 (chunk End=%d should have stopped at block boundary)", len(recs), block2End)
	}
	// All returned records must be < boundaries[0].firstPos
	// (they live entirely inside block 1).
	for i, rec := range recs {
		if rec.Pos >= boundaries[0].firstPos {
			t.Errorf("rec[%d].Pos=%d >= block-2 firstPos=%d; decodeChunk read into block 2",
				i, rec.Pos, boundaries[0].firstPos)
		}
	}

	// Now overwrite the index with one whose End is
	// "the whole file" (high sentinel) and verify the
	// Query yields all 200 records.
	idxFull := &tabix.Index{
		References: []*tabix.Reference{
			{
				Bins: map[uint32][]tabix.Chunk{
					0: {{Beg: 0, End: tabix.VirtualOffset(0xFFFFFFFFFFFFFFFF)}},
				},
				Intervals: []tabix.VirtualOffset{0, tabix.VirtualOffset(0xFFFFFFFFFFFFFFFF)},
			},
		},
	}
	tbiFile, err = os.Create(tbiPath) //nolint:gosec // test file in t.TempDir
	if err != nil {
		t.Fatalf("recreate tbi full: %v", err)
	}
	if err := tabix.NewWriter(tbiFile).Write(idxFull); err != nil {
		_ = tbiFile.Close()
		t.Fatalf("write tbi full: %v", err)
	}
	_ = tbiFile.Close()

	r3, err := OpenTabix(vcfPath)
	if err != nil {
		t.Fatalf("OpenTabix #3: %v", err)
	}
	defer func() { _ = r3.Close() }()
	if _, err := r3.Header(); err != nil {
		t.Fatalf("Header #3: %v", err)
	}
	fullRecs, err := r3.Query("chr1", 1, 200000)
	if err != nil {
		t.Fatalf("Query (full): %v", err)
	}
	if len(fullRecs) != 200 {
		t.Fatalf("Query (full): got %d records, want 200", len(fullRecs))
	}
	if len(recs) >= len(fullRecs) {
		t.Fatalf("chunked query (%d records) should be smaller than full query (%d)", len(recs), len(fullRecs))
	}
}
