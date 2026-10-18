// Package fai reads the .fai (FASTA Index) format used by
// samtools and htslib to enable random access into FASTA
// files.
//
// The .fai format is line-oriented with one tab-separated
// record per reference sequence:
//
//	name<TAB>length<TAB>offset<TAB>lineBases<TAB>lineWidth
//
// Where:
//
//   - name is the FASTA record ID (the first
//     whitespace-delimited token of the `>` header).
//   - length is the number of bases in the record.
//   - offset is the byte offset of the first base in the FASTA
//     file (i.e. the position immediately after the header
//     line's trailing newline).
//   - lineBases is the number of bases per wrapped sequence
//     line.
//   - lineWidth is the byte width of one wrapped sequence
//     line including its newline character.
//
// For records that are not wrapped (lineBases == lineWidth),
// the field has no effect; for wrapped records, the pair
// tells the reader how to skip sequence newlines during
// random-access seeks.
//
// This package is read-only. Index generation is a separate
// concern (multi-pass, atomic-rename) that the v1 release
// does not implement.
package fai

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Reader parses a .fai file from an underlying [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("sequences.fa.fai")
//	defer f.Close()
//	idx, err := fai.NewReader(f).Read()
//	if err != nil { ... }
//	for _, ref := range idx.References {
//	    // ref.Name, ref.Length, ref.Offset, ref.LineBases, ref.LineWidth
//	}
//
// A Reader is not safe for concurrent use.
type Reader struct {
	scanner *bufio.Scanner
}

// NewReader wraps r in a [Reader] that reads one .fai record
// per line.
func NewReader(r io.Reader) *Reader {
	return &Reader{scanner: bufio.NewScanner(r)}
}

// Index is the parsed .fai file. References are stored in
// file order.
type Index struct {
	References []*Reference
}

// Reference is one .fai record.
type Reference struct {
	// Name is the FASTA record ID.
	Name string

	// Length is the number of bases in the record.
	Length int64

	// Offset is the byte offset (in the .fa file) of the
	// first base — i.e. the position immediately after the
	// header line's trailing newline.
	Offset int64

	// LineBases is the number of bases per wrapped sequence
	// line. For records with no wrapping (a single sequence
	// line), LineBases equals LineWidth.
	LineBases int

	// LineWidth is the byte width of one wrapped sequence
	// line, including the trailing newline. For example,
	// 80-char sequence wrapping with a single LF gives
	// LineWidth = 81.
	LineWidth int
}

// Read parses the entire .fai file from r.
func (r *Reader) Read() (*Index, error) {
	idx := &Index{}
	lineNo := 0
	for r.scanner.Scan() {
		lineNo++
		raw := r.scanner.Bytes()
		// Strip trailing CR (CRLF tolerance).
		if n := len(raw); n > 0 && raw[n-1] == '\r' {
			raw = raw[:n-1]
		}
		if len(raw) == 0 {
			continue
		}
		fields := strings.Split(string(raw), "\t")
		if len(fields) != 5 {
			return nil, fmt.Errorf("fai: line %d: want 5 fields, got %d", lineNo, len(fields))
		}
		length, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("fai: line %d: length: %w", lineNo, err)
		}
		offset, err := strconv.ParseInt(fields[2], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("fai: line %d: offset: %w", lineNo, err)
		}
		lineBases, err := strconv.Atoi(fields[3])
		if err != nil {
			return nil, fmt.Errorf("fai: line %d: lineBases: %w", lineNo, err)
		}
		lineWidth, err := strconv.Atoi(fields[4])
		if err != nil {
			return nil, fmt.Errorf("fai: line %d: lineWidth: %w", lineNo, err)
		}
		idx.References = append(idx.References, &Reference{
			Name:      fields[0],
			Length:    length,
			Offset:    offset,
			LineBases: lineBases,
			LineWidth: lineWidth,
		})
	}
	if err := r.scanner.Err(); err != nil {
		if errors.Is(err, bufio.ErrTooLong) {
			return nil, fmt.Errorf("fai: line exceeds bufio.Scanner max (default 64 KiB)")
		}
		return nil, err
	}
	return idx, nil
}