package vcf

import (
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/inful/htsfile/csi"
	"github.com/inful/htsfile/internal/bgzf"
	"github.com/inful/htsfile/tabix"
)

// buildBlockedVCF writes a BGZF-compressed VCF and
// returns the path plus the BGZF virtual offsets of each
// non-EOF block. The fixture is sized so the writer
// produces several blocks; callers use the offsets to
// build an index that mirrors samtools tabix output
// (one chunk per block, End = start of next block).
//
// numRecords is fixed at 5000 across the benchmarks
// because varying it just changes the fixture size
// without exercising different code paths.
func buildBlockedVCF(b *testing.B) (vcfPath string, blockOffs []uint64) {
	const numRecords = 5000
	b.Helper()
	dir := b.TempDir()
	vcfPath = filepath.Join(dir, "bench.vcf.gz")

	pad := randomBases(600)
	var content []byte
	content = append(content, []byte("##fileformat=VCFv4.2\n")...)
	content = append(content, []byte("##contig=<ID=chr1,length=1000000>\n")...)
	content = append(content, []byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")...)
	for i := 1; i <= numRecords; i++ {
		line := fmt.Sprintf("chr1\t%d\t.\tA\tG\t.\t.\tP=%s\n", i*100, pad)
		content = append(content, line...)
	}

	f, err := os.Create(vcfPath) //nolint:gosec // benchmark fixture in b.TempDir
	if err != nil {
		b.Fatal(err)
	}
	bw := bgzf.NewWriter(f)
	if _, err := bw.Write(content); err != nil {
		b.Fatal(err)
	}
	if err := bw.Close(); err != nil {
		b.Fatal(err)
	}
	if err := f.Close(); err != nil {
		b.Fatal(err)
	}

	// Probe each block's start offset with a
	// SeekReader. We don't know the offsets yet, so
	// read forward one byte at a time and track when
	// the BGZF virtual offset's high 48 bits change.
	probeFile, err := os.Open(vcfPath) //nolint:gosec // benchmark fixture in b.TempDir
	if err != nil {
		b.Fatal(err)
	}
	defer func() { _ = probeFile.Close() }()
	sr, err := bgzf.NewSeekReader(probeFile)
	if err != nil {
		b.Fatal(err)
	}
	defer func() { _ = sr.Close() }()
	if err := sr.Seek(0); err != nil {
		b.Fatal(err)
	}

	buf := make([]byte, 1)
	prevHigh := uint64(0)
	blockOffs = append(blockOffs, 0)
	for {
		n, err := sr.Read(buf)
		if n == 0 {
			break
		}
		if err != nil && !errors.Is(err, io.EOF) {
			b.Fatal(err)
		}
		cur := sr.VirtualOffset()
		high := cur >> 16
		if high != prevHigh {
			blockOffs = append(blockOffs, high<<16)
			prevHigh = high
		}
		if errors.Is(err, io.EOF) {
			break
		}
	}
	// Drop the BGZF EOF marker block (28-byte empty
	// trailer at the very end of the file). It would
	// otherwise create a "chunk" pointing at empty
	// bytes, which makes the chunk-End stop bail out
	// with "empty document" when the decoder tries
	// to parse a header from it.
	if len(blockOffs) >= 2 {
		last := blockOffs[len(blockOffs)-1]
		stat, err := os.Stat(vcfPath) //nolint:gosec // benchmark fixture in b.TempDir
		if err == nil {
			fileSize := stat.Size()
			if fileSize > 28 {
				if (last >> 16) >= uint64(fileSize-28) {
					blockOffs = blockOffs[:len(blockOffs)-1]
				}
			}
		}
	}
	return vcfPath, blockOffs
}

// buildTabixForBlockedVCF writes a tabix .tbi whose
// chunks are aligned to the BGZF block boundaries of
// vcfPath. One chunk per block; End = start of next
// block (or MaxUint64 for the last chunk). Mirrors
// samtools' default output.
func buildTabixForBlockedVCF(b *testing.B, vcfPath string, blockOffs []uint64) {
	b.Helper()
	tbiPath := vcfPath + ".tbi"
	tbiFile, err := os.Create(tbiPath) //nolint:gosec // benchmark fixture in b.TempDir
	if err != nil {
		b.Fatal(err)
	}
	chunks := make([]tabix.Chunk, 0, len(blockOffs))
	for i, off := range blockOffs {
		end := uint64(math.MaxUint64)
		if i+1 < len(blockOffs) {
			end = blockOffs[i+1]
		}
		chunks = append(chunks, tabix.Chunk{
			Beg: tabix.VirtualOffset(off),
			End: tabix.VirtualOffset(end),
		})
	}
	idx := &tabix.Index{
		References: []*tabix.Reference{{
			Bins:      map[uint32][]tabix.Chunk{0: chunks},
			Intervals: []tabix.VirtualOffset{tabix.VirtualOffset(0)},
		}},
	}
	if err := tabix.NewWriter(tbiFile).Write(idx); err != nil {
		_ = tbiFile.Close()
		b.Fatal(err)
	}
	_ = tbiFile.Close()
}

// buildTabixWithEnd writes a single-chunk index whose
// End is `chunkEnd`. Used to demonstrate the chunk-End
// stop's effect: the reader breaks out of decodeChunk
// as soon as the BGZF virtual offset reaches chunkEnd.
func buildTabixWithEnd(b *testing.B, vcfPath string, chunkEnd uint64) {
	b.Helper()
	tbiPath := vcfPath + ".tbi"
	tbiFile, err := os.Create(tbiPath) //nolint:gosec // benchmark fixture in b.TempDir
	if err != nil {
		b.Fatal(err)
	}
	idx := &tabix.Index{
		References: []*tabix.Reference{{
			Bins:      map[uint32][]tabix.Chunk{0: {{Beg: 0, End: tabix.VirtualOffset(chunkEnd)}}},
			Intervals: []tabix.VirtualOffset{tabix.VirtualOffset(0)},
		}},
	}
	if err := tabix.NewWriter(tbiFile).Write(idx); err != nil {
		_ = tbiFile.Close()
		b.Fatal(err)
	}
	_ = tbiFile.Close()
}

// buildCSIForBlockedVCF writes a CSI .csi sibling for
// vcfPath using the same chunk list as the tabix index.
// The CSI uses (min_shift=14, depth=5), which produces
// the same bin IDs as tabix for the positions used in
// the fixture.
func buildCSIForBlockedVCF(b *testing.B, vcfPath string, blockOffs []uint64) {
	b.Helper()
	csiPath := vcfPath + ".csi"
	csiFile, err := os.Create(csiPath) //nolint:gosec // benchmark fixture in b.TempDir
	if err != nil {
		b.Fatal(err)
	}
	chunks := make([]csi.Chunk, 0, len(blockOffs))
	for i, off := range blockOffs {
		end := uint64(math.MaxUint64)
		if i+1 < len(blockOffs) {
			end = blockOffs[i+1]
		}
		chunks = append(chunks, csi.Chunk{
			Beg: off,
			End: end,
		})
	}
	idx := &csi.Index{
		MinShift: 14,
		Depth:    5,
		References: []*csi.Reference{{
			Bins:      map[uint32][]csi.Chunk{0: chunks},
			Intervals: []csi.VirtualOffset{0},
		}},
	}
	if err := csi.NewWriter(csiFile).Write(idx); err != nil {
		_ = csiFile.Close()
		b.Fatal(err)
	}
	_ = csiFile.Close()
}

// queryOnceTabix opens a fresh tabix reader, parses
// the header, runs Query, and closes. Mirrors the
// "htsfile vcf query" CLI usage.
func queryOnceTabix(path, chrom string, start, end int64) ([]*Record, error) {
	r, err := OpenTabix(path)
	if err != nil {
		return nil, err
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		return nil, err
	}
	return r.Query(chrom, start, end)
}

// queryOnceCSI is the symmetric CSI helper.
func queryOnceCSI(path, chrom string, start, end int64) ([]*Record, error) {
	r, err := OpenCSI(path)
	if err != nil {
		return nil, err
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		return nil, err
	}
	return r.Query(chrom, start, end)
}

// queryWarmTabix opens the reader once and reuses it
// across iterations. Measures per-query work without
// the BGZF/tabix open overhead.
func queryWarmTabix(b *testing.B, path string, start, end int64) {
	r, err := OpenTabix(path)
	if err != nil {
		b.Fatal(err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := r.Query("chr1", start, end); err != nil {
			b.Fatal(err)
		}
	}
}

// queryWarmCSI is the CSI twin of queryWarmTabix.
func queryWarmCSI(b *testing.B, path string, start, end int64) {
	r, err := OpenCSI(path)
	if err != nil {
		b.Fatal(err)
	}
	defer func() { _ = r.Close() }()
	if _, err := r.Header(); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := r.Query("chr1", start, end); err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkTabixQuery_SmallRegion queries a 10 bp window
// on a 5,000-record multi-block fixture. Only one
// record overlaps; the chunk-End stop is what bounds
// the per-chunk decode work.
func BenchmarkTabixQuery_SmallRegion(b *testing.B) {
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	recs, err := queryOnceTabix(vcfPath, "chr1", 10000, 10010)
	if err != nil {
		b.Fatal(err)
	}
	if len(recs) == 0 {
		b.Fatal("expected at least one record for chr1:10000-10010")
	}
	b.SetBytes(int64(len(recs)))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := queryOnceTabix(vcfPath, "chr1", 10000, 10010); err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkTabixQuery_WholeContig queries the entire
// contig on the same fixture. All records are returned;
// the chunk-End stop does not skip anything because
// every chunk's End is reached.
func BenchmarkTabixQuery_WholeContig(b *testing.B) {
	const n = 5000
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	b.SetBytes(int64(n))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := queryOnceTabix(vcfPath, "chr1", 1, 1000000); err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkCSIQuery_SmallRegion mirrors the tabix
// small-region benchmark on the CSI side.
func BenchmarkCSIQuery_SmallRegion(b *testing.B) {
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	buildCSIForBlockedVCF(b, vcfPath, bOffs)
	recs, err := queryOnceCSI(vcfPath, "chr1", 10000, 10010)
	if err != nil {
		b.Fatal(err)
	}
	if len(recs) == 0 {
		b.Fatal("expected at least one record for chr1:10000-10010")
	}
	b.SetBytes(int64(len(recs)))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := queryOnceCSI(vcfPath, "chr1", 10000, 10010); err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkCSIQuery_WholeContig mirrors the tabix
// whole-contig benchmark on the CSI side.
func BenchmarkCSIQuery_WholeContig(b *testing.B) {
	const n = 5000
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	buildCSIForBlockedVCF(b, vcfPath, bOffs)
	b.SetBytes(int64(n))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := queryOnceCSI(vcfPath, "chr1", 1, 1000000); err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkTabixQuery_Warm_SmallRegion opens the reader
// once and reuses it across iterations. Isolates the
// per-query work from the BGZF/tabix open cost.
func BenchmarkTabixQuery_Warm_SmallRegion(b *testing.B) {
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	queryWarmTabix(b, vcfPath, 10000, 10010)
}

// BenchmarkTabixQuery_Warm_WholeContig is the warm
// variant of BenchmarkTabixQuery_WholeContig.
func BenchmarkTabixQuery_Warm_WholeContig(b *testing.B) {
	const n = 5000
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	b.SetBytes(int64(n))
	queryWarmTabix(b, vcfPath, 1, 1000000)
}

// BenchmarkCSIQuery_Warm_SmallRegion is the warm CSI
// small-region benchmark.
func BenchmarkCSIQuery_Warm_SmallRegion(b *testing.B) {
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	buildCSIForBlockedVCF(b, vcfPath, bOffs)
	queryWarmCSI(b, vcfPath, 10000, 10010)
}

// BenchmarkCSIQuery_Warm_WholeContig is the warm CSI
// whole-contig benchmark.
func BenchmarkCSIQuery_Warm_WholeContig(b *testing.B) {
	const n = 5000
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	buildCSIForBlockedVCF(b, vcfPath, bOffs)
	b.SetBytes(int64(n))
	queryWarmCSI(b, vcfPath, 1, 1000000)
}

// BenchmarkBGZFSeekReader_VirtualOffset measures the
// overhead of calling SeekReader.VirtualOffset on every
// record. Mirrors the new decodeChunk check.
func BenchmarkBGZFSeekReader_VirtualOffset(b *testing.B) {
	vcfPath, _ := buildBlockedVCF(b)
	f, err := os.Open(vcfPath) //nolint:gosec // benchmark fixture in b.TempDir
	if err != nil {
		b.Fatal(err)
	}
	defer func() { _ = f.Close() }()
	sr, err := bgzf.NewSeekReader(f)
	if err != nil {
		b.Fatal(err)
	}
	defer func() { _ = sr.Close() }()
	if err := sr.Seek(0); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = sr.VirtualOffset()
	}
}

// BenchmarkTabixQuery_ChunkEndStop compares two query
// shapes on the same fixture:
//
//   - "no-stop" uses a single chunk with End =
//     MaxUint64, so the reader never short-circuits
//     and reads every record in every block.
//   - "per-block" uses the per-block chunk index (End =
//     next block start), which is what samtools emits
//     by default. The chunk-End stop fires at every
//     block boundary.
//
// Both queries fetch the same records; the per-block
// variant should be no slower, and on narrow queries
// may even be slightly faster because the decoder
// skips records whose BGZF position has crossed End.
func BenchmarkTabixQuery_ChunkEndStop(b *testing.B) {
	vcfPath, bOffs := buildBlockedVCF(b)
	buildTabixWithEnd(b, vcfPath, math.MaxUint64)
	b.Run("no-stop", func(b *testing.B) {
		queryWarmTabix(b, vcfPath, 1, 1000000)
	})
	buildTabixForBlockedVCF(b, vcfPath, bOffs)
	b.Run("per-block", func(b *testing.B) {
		queryWarmTabix(b, vcfPath, 1, 1000000)
	})
}
