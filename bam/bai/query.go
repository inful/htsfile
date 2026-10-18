package bai

// Query describes the set of BGZF regions of a BAM file that
// may contain alignments overlapping a half-open genomic
// range. The caller uses the chunks to seek the BAM stream
// and the MinOffset as a starting point for the linear scan.
type Query struct {
	// Chunks are the (Lo, Hi) BGZF virtual offset pairs
	// across the candidate bins. Lo is the smallest offset
	// of any alignment in the bin; Hi is the largest
	// offset. Callers use these as (start, end) seek
	// hints into the BGZF stream.
	Chunks []Chunk

	// MinOffset is the minimum BGZF virtual offset of any
	// alignment that starts at or after the query range,
	// derived from the linear (interval) index. Callers
	// can seek here as a coarse-grained lower bound.
	MinOffset uint64
}

// Chunk is one (Lo, Hi) BGZF virtual offset pair.
type Chunk struct {
	Lo uint64
	Hi uint64
}

// Query returns the BGZF regions of r that may contain
// alignments overlapping the half-open range [beg, end).
// beg and end are 0-based positions in the reference's
// coordinate space. An empty result (no chunks, MinOffset=0)
// means the index has no entries for this range — either
// the reference is empty or the region is past the indexed
// alignments.
//
// The result is an over-approximation: a chunk may contain
// alignments that do not actually overlap the query range,
// because BAI groups alignments by bin rather than by exact
// position. The caller must scan each chunk and filter to
// the records that genuinely overlap.
//
// The function is O(B + W) where B is the number of
// candidate bins (≤ 6 + log2(region length) per the bin
// scheme) and W is the number of 16384-bp windows in the
// range, capped at the reference's interval count.
func (r *Reference) Query(beg, end int64) Query {
	if r == nil || end <= beg {
		return Query{}
	}
	candidates := Reg2Bins(beg, end)
	// Bin lookup by Number. Bins are stored in ascending
	// Number order (per the writer) so a linear scan is
	// fine here, but the candidate list is short so the
	// constant overhead of a map doesn't pay for itself
	// either — we use a sort.Search-style index walk.
	binIdx := make(map[uint32]*Bin, len(r.Bins))
	for _, b := range r.Bins {
		binIdx[b.Number] = b
	}
	var q Query
	seen := make(map[Chunk]bool, len(candidates))
	for _, num := range candidates {
		b, ok := binIdx[num]
		if !ok {
			continue
		}
		c := Chunk{Lo: b.Lo, Hi: b.Hi}
		if seen[c] {
			continue
		}
		seen[c] = true
		q.Chunks = append(q.Chunks, c)
	}
	// Linear-index min offset. Per the BAI semantics
	// (mirroring htslib's bam_intv_to_offset), this is the
	// smallest non-zero interval offset across all windows
	// from 0 to end/16384 inclusive — i.e. the offset of
	// the earliest alignment whose position could be at or
	// before the query's end. Alignments starting earlier
	// may extend (via CIGAR) into the query range, so we
	// have to scan from the very beginning of the
	// reference. Zero entries mean "no alignment starts
	// in this window" and are skipped.
	if n := int64(len(r.Intervals)); n > 0 {
		last := end >> 14
		if last >= n {
			last = n - 1
		}
		for i := int64(0); i <= last; i++ {
			if off := r.Intervals[i]; off != 0 {
				q.MinOffset = off
				break
			}
		}
	}
	return q
}

// FindOffset returns the smallest BGZF virtual offset of
// any alignment in r that could overlap position pos
// (0-based). Returns 0 if the index has no such alignment
// (no intervals, or all intervals are zero).
//
// Equivalent to Query(pos, pos+1).MinOffset but slightly
// cheaper (no bin lookup).
func (r *Reference) FindOffset(pos int64) uint64 {
	if r == nil || pos < 0 {
		return 0
	}
	n := int64(len(r.Intervals))
	if n == 0 {
		return 0
	}
	last := pos >> 14
	if last >= n {
		return 0
	}
	for i := int64(0); i <= last; i++ {
		if off := r.Intervals[i]; off != 0 {
			return off
		}
	}
	return 0
}