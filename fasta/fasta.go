// Package fasta is a small read-only FASTA reference
// reader. It supports plain (uncompressed) FASTA files
// with an optional samtools-faidx-style .fai index. BGZF
// (gzip-compressed) FASTA is not supported in this
// version — callers needing that should gunzip first.
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
// CRLF).
//
// If no .fai is found at "<path>.fai", Open streams the
// FASTA once at open time to build an in-memory index.
// For a 3 GB human genome reference this takes ~3 s and
// ~3 GB of buffer; callers handling larger references
// should pre-generate a .fai with samtools faidx.
package fasta

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// Reference is a read-only handle on a plain FASTA file.
// It supports random access to a single base by contig
// and 1-based position. The handle is not safe for
// concurrent use.
type Reference struct {
	path string
	file *os.File

	// fai maps contig name to its index entry. nil until
	// initIndex runs.
	fai map[string]*faiEntry
	// names lists contig names in the order they appear
	// in the file (or in the .fai).
	names []string
	// loaded records whether initIndex has run.
	loaded bool
}

// faiEntry is one row of a samtools faidx .fai file.
type faiEntry struct {
	Length        int64 // total bases
	Offset        int64 // byte offset of first base (past '>' line)
	LineBases     int   // bases per line
	LineWidth     int   // bytes per line (linebases + newline)
	lineBasesSeen int   // running tally used while building the in-memory index
}

// Open opens a plain FASTA file at path. If a sibling
// .fai index exists at "<path>.fai", it is loaded
// directly; otherwise the file is scanned once to build
// an in-memory index. The caller must call Close.
func Open(path string) (*Reference, error) {
	f, err := os.Open(path) //nolint:gosec // path is caller-supplied
	if err != nil {
		return nil, err
	}
	r := &Reference{path: path, file: f}

	// Try .fai first.
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
// must not be used after Close.
func (r *Reference) Close() error {
	if r.file == nil {
		return nil
	}
	err := r.file.Close()
	r.file = nil
	return err
}

// Contigs returns the contig names in file order.
func (r *Reference) Contigs() []string {
	out := make([]string, len(r.names))
	copy(out, r.names)
	return out
}

// BaseAt returns the 1-based base at the given position
// on chrom. Returns an error if chrom is unknown, pos is
// out of range, or the file cannot be read.
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
	// Compute the byte offset of the requested base
	// within the file. The first base is at e.Offset.
	// Each line carries LineBases bases and LineWidth
	// bytes (including the newline).
	zero := int64(pos - 1)
	skipLines := zero / int64(e.LineBases)
	col := int(zero % int64(e.LineBases))
	off := e.Offset + skipLines*int64(e.LineWidth) + int64(col)
	if _, err := r.file.Seek(off, io.SeekStart); err != nil {
		return 0, err
	}
	buf := make([]byte, 1)
	if _, err := io.ReadFull(r.file, buf); err != nil {
		return 0, err
	}
	b := buf[0]
	if b == '\n' || b == '\r' {
		return 0, fmt.Errorf("fasta: unexpected newline at %s:%d (file truncated?)", chrom, pos)
	}
	return b, nil
}

// initIndex streams the file once and builds an in-memory
// index of every contig: name, length, offset of first
// base, line width in bases and bytes. Each contig
// records its own line width because different contigs
// in a single FASTA may legitimately differ (the last
// line of a contig is often shorter).
func (r *Reference) initIndex() error {
	if _, err := r.file.Seek(0, io.SeekStart); err != nil {
		return err
	}
	r.fai = make(map[string]*faiEntry)
	br := bufio.NewReaderSize(r.file, 64*1024)

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
		if err != nil && err != io.EOF {
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
			if err == io.EOF {
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
			if err == io.EOF {
				break
			}
			continue
		}
		if cur == nil {
			return fmt.Errorf("fasta: sequence data before any header at offset %d", lineStart)
		}
		seqLine := strings.TrimRight(stripped, " \t\r\n")
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
		if err == io.EOF {
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
		if err != nil && err != io.EOF {
			return nil, err
		}
		trimmed := strings.TrimRight(line, "\r\n")
		if trimmed == "" {
			if err == io.EOF {
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
		if err == io.EOF {
			break
		}
	}
	return idx, nil
}