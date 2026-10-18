// Package bai provides a reader for the BAM Index (BAI) format
// used to enable random access into BAM files.
//
// The BAI format is a two-level index over a BGZF-compressed BAM
// stream:
//
//   - A bin index per reference, using the samtools 6-level
//     binning scheme (37449 bins total covering different genomic
//     intervals at increasing resolution). Each bin records the
//     smallest and largest BGZF virtual offsets that contain
//     alignments overlapping that bin.
//   - A linear (interval) index per reference, partitioning the
//     reference into 16384-bp windows and recording the BGZF
//     virtual offset of the first alignment overlapping (or
//     starting in) each window. Used as a coarse-grained
//     "minimum offset" filter before the bin index is consulted.
//
// BGZF virtual offsets are uint64 values encoded as
// (compressed_file_offset << 16) | (uncompressed_block_offset).
// They are the natural addressing unit for BGZF block-level
// random access.
//
// v1 scope: read-only parsing of the BAI file format plus
// region queries over the parsed index. [Reg2Bins] returns
// the candidate bin numbers for a half-open range; the
// [Reference.Query] method applies that against an already-
// parsed [Reference] and returns a [Query] describing the
// BGZF regions the caller must scan to find alignments
// overlapping the range. [Reference.FindOffset] is a
// position-only variant for the common "give me everything
// from this offset onwards" use case.
package bai

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
)

// Magic is the four-byte BAI file signature ("BAI\1").
var Magic = [4]byte{'B', 'A', 'I', 0x01}

// Index is a parsed BAI index. References are stored in the
// same order as the BAM header's reference-sequence dictionary.
type Index struct {
	// References has one entry per reference sequence in the
	// indexed BAM. len(References) equals n_ref in the BAM
	// header.
	References []*Reference
}

// Reference is the per-reference index. Bins are stored in
// ascending bin-number order; Intervals is the linear index
// (one entry per 16384-bp window).
type Reference struct {
	// Bins is the per-bin chunk index. Each Bin records the
	// smallest and largest BGZF virtual offset of alignments
	// assigned to that bin.
	Bins []*Bin

	// Intervals is the linear index: one BGZF virtual offset
	// per 16384-bp window across the reference. Intervals[i]
	// is the offset of the first alignment that starts in
	// window i (or any later window if no alignment starts in
	// window i).
	Intervals []uint64

	// NoCoor is the count of unplaced reads (records with no
	// coordinate). Optional in older files; zero if absent.
	NoCoor uint64
}

// Bin is one entry in the bin index. Number is the bin number
// (0..=37448); Lo and Hi are BGZF virtual offsets bounding
// the alignments assigned to this bin.
type Bin struct {
	Number uint32
	Lo     uint64
	Hi     uint64
}

// Reader parses a BAI index from an underlying [io.Reader].
//
// Typical usage:
//
//	f, _ := os.Open("sample.bam.bai")
//	defer f.Close()
//	idx, err := bai.NewReader(f).Read()
//	if err != nil { ... }
//	for _, ref := range idx.References {
//	    for _, b := range ref.Bins { ... }
//	}
//
// A Reader is not safe for concurrent use.
type Reader struct {
	r io.Reader
}

// NewReader wraps r in a [Reader]. r is expected to be a
// BAI-formatted stream.
func NewReader(r io.Reader) *Reader {
	return &Reader{r: r}
}

// Read parses the entire BAI file from the underlying reader
// and returns the resulting [Index]. The full file is loaded
// into memory; BAI files are typically small (a few hundred
// bytes per reference per bin + 16 KB of intervals per
// reference).
func (r *Reader) Read() (*Index, error) {
	// Magic.
	var magic [4]byte
	if _, err := io.ReadFull(r.r, magic[:]); err != nil {
		if errors.Is(err, io.EOF) {
			return nil, io.ErrUnexpectedEOF
		}
		return nil, err
	}
	if magic != Magic {
		return nil, fmt.Errorf("bai: bad magic: got %q", magic[:])
	}
	// n_ref.
	var nRef uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nRef); err != nil {
		return nil, fmt.Errorf("bai: read n_ref: %w", err)
	}
	idx := &Index{References: make([]*Reference, 0, nRef)}
	for i := uint32(0); i < nRef; i++ {
		ref, err := r.readReference()
		if err != nil {
			return nil, fmt.Errorf("bai: reference %d: %w", i, err)
		}
		idx.References = append(idx.References, ref)
	}
	return idx, nil
}

func (r *Reader) readReference() (*Reference, error) {
	// n_bin + bin entries.
	var nBin uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nBin); err != nil {
		return nil, fmt.Errorf("read n_bin: %w", err)
	}
	bins := make([]*Bin, nBin)
	for i := uint32(0); i < nBin; i++ {
		var b Bin
		if err := binary.Read(r.r, binary.LittleEndian, &b.Number); err != nil {
			return nil, fmt.Errorf("read bin[%d].number: %w", i, err)
		}
		if err := binary.Read(r.r, binary.LittleEndian, &b.Lo); err != nil {
			return nil, fmt.Errorf("read bin[%d].lo: %w", i, err)
		}
		if err := binary.Read(r.r, binary.LittleEndian, &b.Hi); err != nil {
			return nil, fmt.Errorf("read bin[%d].hi: %w", i, err)
		}
		bins[i] = &b
	}
	// n_intv + intervals.
	var nIntv uint32
	if err := binary.Read(r.r, binary.LittleEndian, &nIntv); err != nil {
		return nil, fmt.Errorf("read n_intv: %w", err)
	}
	intervals := make([]uint64, nIntv)
	for i := uint32(0); i < nIntv; i++ {
		if err := binary.Read(r.r, binary.LittleEndian, &intervals[i]); err != nil {
			return nil, fmt.Errorf("read interval[%d]: %w", i, err)
		}
	}
	// n_no_coor (optional in older files). We probe for it by
	// attempting to read 8 bytes; if EOF, NoCoor is 0.
	ref := &Reference{Bins: bins, Intervals: intervals}
	var noCoor uint64
	if err := binary.Read(r.r, binary.LittleEndian, &noCoor); err != nil {
		if !errors.Is(err, io.EOF) && !errors.Is(err, io.ErrUnexpectedEOF) {
			return nil, fmt.Errorf("read n_no_coor: %w", err)
		}
		// EOF here means the file ended without n_no_coor —
		// older format.
	} else {
		ref.NoCoor = noCoor
	}
	return ref, nil
}

// Reg2Bins returns the set of BAI bin numbers that overlap the
// half-open genomic range [beg, end). beg and end are 0-based
// positions in the reference's coordinate space.
//
// The bin assignment follows the samtools 6-level scheme
// (37449 bins total). This helper is the inverse of the
// bin-assignment step used by BAM encoders when writing BAI
// files; it's exported so future region-query code can ask
// "which bins contain region X?" and look them up in an
// already-parsed [Index].
//
// Based on samtools/htslib's reg2bins (sam.h).
func Reg2Bins(beg, end int64) []uint32 {
	if end <= beg {
		return nil
	}
	end--
	var bins []uint32
	for k := 0; k < 6; k++ {
		switch k {
		case 0:
			bins = append(bins, 0)
		case 1:
			for b := 1 + (beg >> 26); b <= 1+(end>>26); b++ {
				bins = append(bins, uint32(b))
			}
		case 2:
			for b := 9 + (beg >> 23); b <= 9+(end>>23); b++ {
				bins = append(bins, uint32(b))
			}
		case 3:
			for b := 73 + (beg >> 20); b <= 73+(end>>20); b++ {
				bins = append(bins, uint32(b))
			}
		case 4:
			for b := 585 + (beg >> 17); b <= 585+(end>>17); b++ {
				bins = append(bins, uint32(b))
			}
		case 5:
			for b := 4681 + (beg >> 14); b <= 4681+(end>>14); b++ {
				bins = append(bins, uint32(b))
			}
		}
	}
	return bins
}
