// Package fasta is a small read-only FASTA reference
// reader. It supports plain and BGZF-compressed FASTA
// files with an optional samtools-faidx-style .fai
// index.
//
// The .fai format (samtools faidx output) is a
// tab-separated, one-record-per-line index:
//
//	name<TAB>length<TAB>offset<TAB>linebases<TAB>linewidth
//
// where `offset` is the byte offset of the first base
// (just past the header line), `linebases` is the number
// of bases per line (typically 60, 70, or 80), and
// `linewidth` is the number of bytes per line including
// the line terminator (linebases + 1 for LF, +2 for
// CRLF). For BGZF files, .fai offsets are virtual
// (compressed-stream positions) and require a seekable
// BGZF reader; this package ignores .fai for BGZF and
// falls back to streaming.
//
// If no .fai is found at "<path>.fai", Open streams the
// FASTA once at open time to build an in-memory index.
// For plain FASTA, file access then uses a seekable
// file handle. For BGZF, the whole stream is decompressed
// into memory and accessed via a *bytes.Reader.
package fasta

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/inful/htsfile/internal/bgzf"
)

// Reference is a read-only handle on a FASTA file
// (plain or BGZF). It supports random access to a
// single base by contig and 1-based position, plus
// efficient multi-base range reads. The handle is not
// safe for concurrent use.
type Reference struct {
	path string

	// rs is the seekable source for reads. For plain
	// FASTA it's the underlying *os.File. For BGZF
	// (decompressed into memory at Open) it's a
	// *bytes.Reader wrapping the decompressed bytes.
	rs io.ReadSeeker
	// bgzfData is non-nil only when rs wraps decompressed
	// BGZF content; we hold the slice so the GC doesn't
	// collect it from under the reader.
	bgzfData []byte

	// fai maps contig name to its index entry. nil until
	// initIndex runs.
	fai map[string]*faiEntry
	// names lists contig names in the order they appear
	// in the file (or in the .fai).
	names []string
	// loaded records whether initIndex has run.

	// Cursor state for the scanned path. When callers
	// issue BaseAt/BaseRange calls in order on the same
	// contig, we resume scanning from the cursor instead
	// of restarting at the first base of the contig.
	// This makes sequential VCF validation O(N) instead
	// of O(N^2). Reset on chrom change, on backward
	// seeks, and on .fai-path lookups.
	cursorCh        string
	cursorOffset    int64
	cursorBasesSeen int64
	loaded          bool
}

// faiEntry is one row of a samtools faidx .fai file.
type faiEntry struct {
	Length        int64   // total bases
	Offset        int64   // byte offset of first base (past '>' line)
	LineBases     int     // bases per line (used by .fai path)
	LineWidth     int     // bytes per line including terminator (used by .fai path)
	lineBasesSeen int     // running tally used while building the in-memory index
	lineOffsets   []int64 // byte offset of each sequence line (scanned path only)
}

// Open opens a FASTA file at path. Plain FASTA is
// accessed via a file handle. BGZF (gzip-compressed)
// FASTA is detected by magic and decompressed into
// memory at open time; the caller pays a memory cost
// equal to the decompressed size. If a sibling .fai
// index exists at "<path>.fai", it is loaded directly
// for plain FASTA; for BGZF the .fai is ignored (its
// virtual offsets require a seekable BGZF reader).
//
// The caller must call Close.
func Open(path string) (*Reference, error) {
	f, err := os.Open(path) //nolint:gosec // path is caller-supplied
	if err != nil {
		return nil, err
	}

	// Peek the first few bytes to detect BGZF. The
	// BGZF magic is the gzip header with FEXTRA flag
	// (0x04) set in byte 3.
	hdr := make([]byte, 4)
	n, _ := io.ReadFull(f, hdr)
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		_ = f.Close()
		return nil, err
	}
	isBGZF := n == 4 && hdr[0] == 0x1f && hdr[1] == 0x8b && hdr[2] == 0x08 && hdr[3]&0x04 != 0

	if isBGZF {
		// Decompress the entire BGZF stream into memory.
		// BGZF .fai offsets are virtual; without a
		// seekable BGZF reader we can't use them, so
		// the streamed scan index is used.
		br := bgzf.NewReader(f)
		data, derr := io.ReadAll(br)
		if derr != nil {
			_ = f.Close()
			return nil, fmt.Errorf("fasta: BGZF decompression of %s: %w", path, derr)
		}
		_ = f.Close()
		r := &Reference{path: path, rs: bytes.NewReader(data), bgzfData: data}
		if err := r.initIndex(); err != nil {
			return nil, fmt.Errorf("fasta: build index for %s: %w", path, err)
		}
		return r, nil
	}

	r := &Reference{path: path, rs: f}

	// Try .fai first (plain FASTA only).
	idx, ierr := readFAI(path + ".fai")
	if ierr == nil && idx != nil {
		r.fai = idx
		r.loaded = true
		for k := range idx {
			r.names = append(r.names, k)
		}
		return r, nil
	}

	// Fallback: stream the file once to build the index.
	if err := r.initIndex(); err != nil {
		_ = f.Close()
		return nil, fmt.Errorf("fasta: build index for %s: %w", path, err)
	}
	return r, nil
}

// Close releases the underlying file. The Reference
// must not be used after Close. For BGZF references
// (which hold the decompressed bytes in memory), Close
// is a no-op.
func (r *Reference) Close() error {
	if f, ok := r.rs.(*os.File); ok && f != nil {
		err := f.Close()
		r.rs = nil
		return err
	}
	r.rs = nil
	return nil
}

// Contigs returns the contig names in file order.
func (r *Reference) Contigs() []string {
	out := make([]string, len(r.names))
	copy(out, r.names)
	return out
}

// Length returns the number of bases in the named contig.
// Returns an error if the contig is unknown.
func (r *Reference) Length(chrom string) (int64, error) {
	if !r.loaded {
		if err := r.initIndex(); err != nil {
			return 0, err
		}
	}
	e, ok := r.fai[chrom]
	if !ok {
		return 0, fmt.Errorf("fasta: contig %q not found", chrom)
	}
	return e.Length, nil
}

// BaseAt returns the 1-based base at the given position
// on chrom. Returns an error if chrom is unknown, pos is
// out of range, or the file cannot be read.
//
// Two access paths:
//   - If the index was loaded from a sibling .fai, we
//     use the simple offset+skipLines*linewidth+col
//     formula. .fai is the canonical samtools format
//     and lines are assumed contiguous and free of
//     internal whitespace / digits.
//   - If the index was built by streaming the file, we
//     seek to the first-base offset and scan forward,
//     line by line, counting bases while skipping
//     spaces, tabs, asterisks, and digits. This is
//     robust against comment lines and internal
//     formatting. Cost is O(target-pos / lineBases).
func (r *Reference) BaseAt(chrom string, pos int) (byte, error) {
	if pos < 1 {
		return 0, fmt.Errorf("fasta: position %d is 0-based; use 1-based", pos)
	}
	if !r.loaded {
		if err := r.initIndex(); err != nil {
			return 0, err
		}
	}
	e, ok := r.fai[chrom]
	if !ok {
		return 0, fmt.Errorf("fasta: contig %q not found", chrom)
	}
	if int64(pos) > e.Length {
		return 0, fmt.Errorf("fasta: %s:%d past contig length %d", chrom, pos, e.Length)
	}
	zero := int64(pos - 1)
	if len(e.lineOffsets) == 0 {
		// .fai path: assume contiguous lines. This
		// seeks arbitrarily so the cursor is invalidated.
		r.cursorCh = ""
		skipLines := zero / int64(e.LineBases)
		col := int(zero % int64(e.LineBases))
		off := e.Offset + skipLines*int64(e.LineWidth) + int64(col)
		if _, err := r.rs.Seek(off, io.SeekStart); err != nil {
			return 0, err
		}
		buf := make([]byte, 1)
		if _, err := io.ReadFull(r.rs, buf); err != nil {
			return 0, err
		}
		b := buf[0]
		if b == '\n' || b == '\r' {
			return 0, fmt.Errorf("fasta: unexpected newline at %s:%d (file truncated?)", chrom, pos)
		}
		return b, nil
	}
	// Scanned path: resume from the cursor when the next
	// call targets the same contig at a forward position.
	// Otherwise restart at the first-base offset.
	startOffset := e.Offset
	startSeen := int64(0)
	if r.cursorCh == chrom && zero >= r.cursorBasesSeen {
		startOffset = r.cursorOffset
		startSeen = r.cursorBasesSeen
	} else {
		r.cursorCh = ""
	}
	// Read line by line without bufio so the file offset
	// is exact and the cursor can be advanced precisely.
	if _, err := r.rs.Seek(startOffset, io.SeekStart); err != nil {
		return 0, err
	}
	seen := startSeen
	for {
		line, rawLen, off, err := readLineTracked(r.rs)
		if line == "" && errors.Is(err, io.EOF) {
			break
		}
		if err != nil && !errors.Is(err, io.EOF) {
			return 0, err
		}
		stripped := strings.TrimRight(line, "\r\n")
		if stripped == "" || strings.HasPrefix(stripped, ";") || strings.HasPrefix(stripped, ">") {
			if errors.Is(err, io.EOF) {
				break
			}
			if strings.HasPrefix(stripped, ">") {
				return 0, fmt.Errorf("fasta: %s:%d not found (file may be truncated)", chrom, pos)
			}
			continue
		}
		for i := 0; i < len(stripped); i++ {
			c := stripped[i]
			if c == ' ' || c == '\t' || c == '*' {
				continue
			}
			if c >= '0' && c <= '9' {
				continue
			}
			if seen == zero {
				// Advance cursor to the byte right
				// after the one we just returned.
				r.cursorCh = chrom
				r.cursorOffset = off - int64(rawLen) + int64(i) + 1
				r.cursorBasesSeen = seen + 1
				return c, nil
			}
			seen++
		}
		if errors.Is(err, io.EOF) {
			break
		}
	}
	r.cursorCh = ""
	return 0, fmt.Errorf("fasta: %s:%d past contig length %d", chrom, pos, e.Length)
}

// readLineTracked reads one line (up to '\n' or EOF) from
// rs and returns the line (without trailing newline), the
// raw line length on disk (line + terminator), and the
// file offset of the byte right after the line's
// terminator (or, if EOF was hit, right after the last
// byte returned).
//
// Reads are buffered via a 4 KiB scratch buffer so a
// line that's longer than the buffer still terminates
// at '\n' or EOF.
func readLineTracked(rs io.ReadSeeker) (line string, rawLen int, off int64, err error) {
	const chunkSize = 4096
	var buf []byte
	scratch := make([]byte, chunkSize)
	eof := false
	for {
		n, rerr := rs.Read(scratch)
		if n > 0 {
			for i := 0; i < n; i++ {
				buf = append(buf, scratch[i])
				if scratch[i] == '\n' {
					trimmed := buf
					if len(trimmed) > 0 && trimmed[len(trimmed)-1] == '\n' {
						trimmed = trimmed[:len(trimmed)-1]
					}
					if len(trimmed) > 0 && trimmed[len(trimmed)-1] == '\r' {
						trimmed = trimmed[:len(trimmed)-1]
					}
					cur, _ := rs.Seek(0, io.SeekCurrent)
					// Read may have returned more bytes than
					// we needed for this line; seek to the
					// actual position past the terminator so
					// the next call to Read picks up at the
					// right place.
					curPos := cur - int64(n-i-1)
					if _, err := rs.Seek(curPos, io.SeekStart); err != nil {
						return "", 0, 0, err
					}
					return string(trimmed), len(buf), curPos, nil
				}
			}
		}
		if errors.Is(err, io.EOF) {
			eof = true
			break
		}
		if rerr != nil {
			return "", 0, 0, rerr
		}
		if n == 0 {
			break
		}
	}
	cur, _ := rs.Seek(0, io.SeekCurrent)
	trimmed := buf
	if len(trimmed) > 0 && trimmed[len(trimmed)-1] == '\r' {
		trimmed = trimmed[:len(trimmed)-1]
	}
	if eof {
		return string(trimmed), len(buf), cur, io.EOF
	}
	return string(trimmed), len(buf), cur, nil
}

// BaseRange returns the bases at 1-based positions start
// through end-1 (inclusive of start, exclusive of end) as
// a byte slice of length end-start. Returns an error if
// chrom is unknown, the range is out of bounds, or the
// file cannot be read.
//
// BaseRange is more efficient than calling [BaseAt] in a
// loop because it does a single forward scan instead of
// N independent seeks + scans. The two paths use the same
// correctness machinery as BaseAt: skip spaces, tabs,
// asterisks, and digits in the scanned path; trust the
// .fai layout in the indexed path.
func (r *Reference) BaseRange(chrom string, start, end int) ([]byte, error) {
	if start < 1 {
		return nil, fmt.Errorf("fasta: start %d is 0-based; use 1-based", start)
	}
	if end <= start {
		return nil, fmt.Errorf("fasta: range end %d must be > start %d", end, start)
	}
	if !r.loaded {
		if err := r.initIndex(); err != nil {
			return nil, err
		}
	}
	e, ok := r.fai[chrom]
	if !ok {
		return nil, fmt.Errorf("fasta: contig %q not found", chrom)
	}
	if int64(end-1) > e.Length {
		return nil, fmt.Errorf("fasta: %s:%d past contig length %d", chrom, end-1, e.Length)
	}
	wantLen := end - start
	if len(e.lineOffsets) == 0 {
		// .fai path: each base is one byte on disk.
		off := e.Offset + int64(start-1)
		if _, err := r.rs.Seek(off, io.SeekStart); err != nil {
			return nil, err
		}
		buf := make([]byte, wantLen)
		n, err := io.ReadFull(r.rs, buf)
		if err != nil {
			return nil, err
		}
		if n != wantLen {
			return nil, fmt.Errorf("fasta: short read for %s:%d-%d", chrom, start, end-1)
		}
		for _, b := range buf {
			if b == '\n' || b == '\r' {
				return nil, fmt.Errorf("fasta: unexpected newline in %s:%d-%d (file truncated?)", chrom, start, end-1)
			}
		}
		return buf, nil
	}
	// Scanned path: resume from the cursor when the next
	// call targets the same contig at a forward position.
	startOffset := e.Offset
	startSeen := int64(0)
	if r.cursorCh == chrom && int64(start-1) >= r.cursorBasesSeen {
		startOffset = r.cursorOffset
		startSeen = r.cursorBasesSeen
	} else {
		r.cursorCh = ""
	}
	if _, err := r.rs.Seek(startOffset, io.SeekStart); err != nil {
		return nil, err
	}
	seen := startSeen
	out := make([]byte, 0, wantLen)
	targetStart := int64(start - 1)
	targetEnd := int64(end - 1)
	for {
		line, rawLen, off, err := readLineTracked(r.rs)
		if line == "" && errors.Is(err, io.EOF) {
			break
		}
		if err != nil && !errors.Is(err, io.EOF) {
			return nil, err
		}
		stripped := strings.TrimRight(line, "\r\n")
		if stripped == "" || strings.HasPrefix(stripped, ";") || strings.HasPrefix(stripped, ">") {
			if errors.Is(err, io.EOF) {
				break
			}
			if strings.HasPrefix(stripped, ">") {
				return nil, fmt.Errorf("fasta: %s:%d-%d not found (file may be truncated)", chrom, start, end-1)
			}
			continue
		}
		for i := 0; i < len(stripped); i++ {
			c := stripped[i]
			if c == ' ' || c == '\t' || c == '*' {
				continue
			}
			if c >= '0' && c <= '9' {
				continue
			}
			if seen >= targetStart && seen < targetEnd {
				out = append(out, c)
			}
			seen++
			if seen >= targetEnd {
				r.cursorCh = chrom
				r.cursorOffset = off - int64(rawLen) + int64(i) + 1
				r.cursorBasesSeen = seen
				return out, nil
			}
		}
		if errors.Is(err, io.EOF) {
			break
		}
	}
	r.cursorCh = ""
	return nil, fmt.Errorf("fasta: %s:%d-%d past contig length %d", chrom, start, end-1, e.Length)
}

// stripIgnoredChars returns the sequence line with the
// characters that the FASTA format explicitly says are
// ignored (spaces, tabs, asterisks, and digits) removed.
// Used by [initIndex] to compute the true base count of
// each line so the offset arithmetic stays correct for
// files written with human-readable spacing. Digits are
// also stripped because Wikipedia notes they "are not
// allowed but are used in some databases to indicate
// the position in the sequence" — counting them would
// misreport the base total.
func stripIgnoredChars(line string) string {
	out := make([]byte, 0, len(line))
	for i := 0; i < len(line); i++ {
		c := line[i]
		if c == ' ' || c == '\t' || c == '*' {
			continue
		}
		if c >= '0' && c <= '9' {
			continue
		}
		out = append(out, c)
	}
	return string(out)
}

// initIndex streams the file once and builds an in-memory
// index of every contig: name, length, offset of first
// base, line width in bases and bytes. Each contig
// records its own line width because different contigs
// in a single FASTA may legitimately differ (the last
// line of a contig is often shorter).
func (r *Reference) initIndex() error {
	if _, err := r.rs.Seek(0, io.SeekStart); err != nil {
		return err
	}
	r.fai = make(map[string]*faiEntry)
	br := bufio.NewReaderSize(r.rs, 64*1024)

	var cur *faiEntry
	// off tracks the byte offset just past the '>'
	// header line — where the first base sits.
	off := int64(0)

	for {
		lineStart := off
		line, err := br.ReadString('\n')
		if line != "" {
			off += int64(len(line))
		}
		if err != nil && !errors.Is(err, io.EOF) {
			return err
		}
		// The line as stored on disk (with its line
		// terminator) has total bytes (off - lineStart).
		// We use that as lineWidth directly — it equals
		// bases-per-line plus 1 (LF) for Unix endings,
		// or plus 2 (CRLF). ReadString splits on '\n'
		// only, so any '\r' is included in `line`.
		stripped := strings.TrimRight(line, "\r\n")
		if stripped == "" {
			if errors.Is(err, io.EOF) {
				break
			}
			continue
		}
		if stripped[0] == '>' {
			if cur != nil {
				cur.Length = int64(cur.lineBasesSeen)
			}
			name := strings.TrimPrefix(stripped, ">")
			if sp := strings.IndexAny(name, " \t"); sp >= 0 {
				name = name[:sp]
			}
			cur = &faiEntry{Offset: off}
			r.fai[name] = cur
			r.names = append(r.names, name)
			if errors.Is(err, io.EOF) {
				break
			}
			continue
		}
		if cur == nil {
			return fmt.Errorf("fasta: sequence data before any header at offset %d", lineStart)
		}
		// Deprecated Pearson `;` comment lines appear
		// after the header and were meant to be ignored.
		// Wikipedia (FASTA format): "one or more comments,
		// distinguished by a semi-colon at the beginning
		// of the line, may occur after the header". We
		// skip them so they don't corrupt the offset
		// arithmetic.
		if strings.HasPrefix(stripped, ";") {
			if errors.Is(err, io.EOF) {
				break
			}
			continue
		}
		// Internal whitespace (spaces, tabs) is allowed
		// inside sequence lines per the format ("Anything
		// other than a valid character would be ignored,
		// including spaces, tabulators, asterisks, etc").
		// Strip them so len(seqLine) is the base count,
		// not the byte count.
		seqLine := stripIgnoredChars(stripped)
		// Per-contig line width: use the first line we
		// see for this contig. The last line of a contig
		// is usually shorter, but the formula
		// (skipLines = zero / lineBases, col = zero %
		// lineBases) still produces correct results
		// when lineBases is the full-width line count.
		if cur.LineBases == 0 {
			cur.LineBases = len(seqLine)
			cur.LineWidth = int(off - lineStart)
		}
		cur.lineBasesSeen += len(seqLine)
		cur.Length = int64(cur.lineBasesSeen)
		// Record the file offset of each sequence line so
		// BaseAt can seek directly to the right line even
		// when comment lines or blank lines appear in
		// between (which would break the simple
		// offset+skipLines*linewidth formula).
		cur.lineOffsets = append(cur.lineOffsets, lineStart)
		if errors.Is(err, io.EOF) {
			break
		}
	}
	r.loaded = true
	return nil
}

// readFAI parses a samtools faidx .fai file. Each line
// is `name<TAB>length<TAB>offset<TAB>linebases<TAB>linewidth`.
// Returns nil if the file cannot be opened (caller should
// fall back to scanning the FASTA directly).
func readFAI(path string) (map[string]*faiEntry, error) {
	f, err := os.Open(path) //nolint:gosec // path is sibling of caller-supplied
	if err != nil {
		return nil, err
	}
	defer func() { _ = f.Close() }()
	idx := make(map[string]*faiEntry)
	br := bufio.NewReader(f)
	for {
		line, err := br.ReadString('\n')
		if err != nil && !errors.Is(err, io.EOF) {
			return nil, err
		}
		trimmed := strings.TrimRight(line, "\r\n")
		if trimmed == "" {
			if errors.Is(err, io.EOF) {
				break
			}
			continue
		}
		parts := strings.Split(trimmed, "\t")
		if len(parts) != 5 {
			return nil, fmt.Errorf("fasta: malformed .fai line %q (expected 5 tab-separated fields)", trimmed)
		}
		length, err1 := strconv.ParseInt(parts[1], 10, 64)
		offset, err2 := strconv.ParseInt(parts[2], 10, 64)
		lineBases, err3 := strconv.Atoi(parts[3])
		lineWidth, err4 := strconv.Atoi(parts[4])
		if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
			return nil, errors.New("fasta: malformed numeric fields in .fai")
		}
		idx[parts[0]] = &faiEntry{
			Length:    length,
			Offset:    offset,
			LineBases: lineBases,
			LineWidth: lineWidth,
		}
		if errors.Is(err, io.EOF) {
			break
		}
	}
	return idx, nil
}
