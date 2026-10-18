// Package csi provides a reader for the CSI (Coordinate-Sorted
// Index) v1 format — the successor to Tabix used by
// samtools/htslib ≥1.10. CSI is the generic index for sorted,
// bgzipped, position-sorted text files (VCF, GFF, BED, SAM,
// BCF, BAM).
//
// The CSI format is similar to Tabix but adds:
//   - An explicit version marker (4-byte magic "CSI\x01")
//   - Tunable bin sizes via min_shift (the linear-index
//     window size is 2^min_shift bp; default 14 → 16 KiB
//     windows matching Tabix)
//   - Tunable depth (number of bin levels; default 5)
//
// Format:
//
//	magic    : "CSI\x01"             (4 bytes)
//	min_shift: int32 little-endian
//	depth    : int32 little-endian
//	n_ref    : uint32 little-endian
//	For each ref:
//	    n_bin      : uint32
//	    For each bin:
//	        bin      : uint32
//	        n_chunk  : uint32
//	        For each chunk:
//	            cnk_beg : uint64 BGZF virtual offset
//	            cnk_end : uint64 BGZF virtual offset
//	    n_intv     : uint32
//	    For each interval:
//	        ioff      : uint64 BGZF virtual offset
//
// CSI does NOT have an optional n_no_coor field; the bin and
// linear index together fully describe the indexed positions.
//
// v1 is read-only. CSI generation is out of scope per the
// ROADMAP's "Index generation" out-of-scope note.
package csi