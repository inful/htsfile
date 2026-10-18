// Package tabix provides a reader for the Tabix index format
// (file extension .tbi), the general-purpose index for
// sorted, bgzipped, position-sorted text files (VCF, GFF,
// BED, SAM, …).
//
// The Tabix format is line-oriented and per-reference:
//
//	magic    : "TBI\x01"            (4 bytes)
//	n_ref    : uint32 little-endian
//	For each ref:
//	    n_bin      : uint32
//	    For each bin:
//	        bin        : uint32 (samtools 5-level bin number)
//	        n_chunk    : uint32
//	        For each chunk:
//	            cnk_beg : uint64 BGZF virtual offset
//	            cnk_end : uint64 BGZF virtual offset
//	    n_intv     : uint32
//	    For each interval:
//	        ioff      : uint64 BGZF virtual offset
//	    [optional n_no_coor (uint64) — older files]
//
// Bin numbers follow the samtools 5-level binning scheme
// (levels 0-4680). Virtual offsets are 64-bit BGZF addresses:
// the high 48 bits are the compressed-block offset; the low
// 16 bits are the offset within the decompressed block.
//
// v1 is read-only. Tabix generation is out of scope per
// the ROADMAP's "Index generation" out-of-scope note.
package tabix
