package bed

import (
	"fmt"
	"sort"
)

// Severity classifies a validator [Issue] by impact.
type Severity int

const (
	// SeverityWarning is a suspicious condition that does not make
	// the data unusable but may indicate a mistake.
	SeverityWarning Severity = iota
	// SeverityError is a definite problem: malformed data, a
	// violated invariant, or a spec violation.
	SeverityError
)

// String returns the lowercase name of the severity.
func (s Severity) String() string {
	switch s {
	case SeverityWarning:
		return "warning"
	case SeverityError:
		return "error"
	default:
		return "unknown"
	}
}

// Issue is one problem found by the [Validator].
type Issue struct {
	Line     int
	Field    string
	Severity Severity
	Msg      string
}

// Mode controls the strictness of the [Validator].
type Mode int

const (
	// ModeDefault checks structural well-formedness: column-count
	// consistency, coordinate sanity, optional-field invariants,
	// and BED12 block invariants.
	ModeDefault Mode = iota
	// ModeStrict adds sort-order checks: records must be sorted by
	// (chrom, chromStart).
	ModeStrict
)

// Validator runs semantic checks against a parsed BED document.
type Validator struct {
	Mode Mode
}

// NewValidator returns a Validator with the given mode.
func NewValidator(m Mode) *Validator {
	return &Validator{Mode: m}
}

// Validate runs all applicable checks and returns the issues found.
// It does not short-circuit; every issue is reported.
func (v *Validator) Validate(track *Track, recs []*Record) []Issue {
	var issues []Issue
	// Default checks.
	issues = append(issues, v.checkColumnConsistency(recs)...)
	issues = append(issues, v.checkRecords(recs)...)
	// Strict additions.
	if v.Mode == ModeStrict {
		issues = append(issues, v.checkSortOrder(recs)...)
	}
	return issues
}

// checkColumnConsistency verifies that all records have the same
// column count (ColumnCount() returns the same value for all).
func (v *Validator) checkColumnConsistency(recs []*Record) []Issue {
	var issues []Issue
	if len(recs) == 0 {
		return issues
	}
	want := recs[0].ColumnCount()
	for i, r := range recs[1:] {
		if r.ColumnCount() != want {
			issues = append(issues, Issue{
				Line:     r.Line,
				Severity: SeverityError,
				Msg: fmt.Sprintf("column count %d differs from first record's %d",
					r.ColumnCount(), want),
			})
		}
		_ = i
	}
	return issues
}

// checkRecords runs per-record checks.
func (v *Validator) checkRecords(recs []*Record) []Issue {
	issues := make([]Issue, 0, len(recs))
	for _, r := range recs {
		issues = append(issues, v.checkOne(r)...)
	}
	return issues
}

func (v *Validator) checkOne(r *Record) []Issue {
	var issues []Issue
	// chromStart >= 0.
	if r.ChromStart < 0 {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "chromStart",
			Severity: SeverityError,
			Msg:      fmt.Sprintf("chromStart must be >= 0; got %d", r.ChromStart),
		})
	}
	// chromStart <= chromEnd.
	if r.ChromStart > r.ChromEnd {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "chromEnd",
			Severity: SeverityError,
			Msg:      fmt.Sprintf("chromStart (%d) > chromEnd (%d)", r.ChromStart, r.ChromEnd),
		})
	}
	// Score 0..=1000.
	if r.Score != nil {
		if *r.Score < 0 || *r.Score > 1000 {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "score",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("score must be in [0, 1000]; got %d", *r.Score),
			})
		}
	}
	// Contiguous prefix check.
	if r.hasGapBeforeColumn(r.ColumnCount()) {
		issues = append(issues, Issue{
			Line:     r.Line,
			Severity: SeverityError,
			Msg:      fmt.Sprintf("non-contiguous optional-field prefix: %d columns but a lower field is unset", r.ColumnCount()),
		})
	}
	// BED12 block invariants.
	if r.HasBlocks() {
		issues = append(issues, v.checkBlocks(r)...)
	}
	return issues
}

// hasGapBeforeColumn reports whether r has column n populated but a
// lower-numbered optional column unset (e.g. Score set but Name nil).
func (r *Record) hasGapBeforeColumn(n int) bool {
	switch n {
	case 12:
		if r.BlockStarts == nil {
			return true
		}
		fallthrough
	case 11:
		if r.BlockSizes == nil {
			return true
		}
		fallthrough
	case 10:
		if r.BlockCount == nil {
			return true
		}
		fallthrough
	case 9:
		if r.ItemRgb == nil {
			return true
		}
		fallthrough
	case 8:
		if r.ThickEnd == nil {
			return true
		}
		fallthrough
	case 7:
		if r.ThickStart == nil {
			return true
		}
		fallthrough
	case 6:
		if r.Strand == nil {
			return true
		}
		fallthrough
	case 5:
		if r.Score == nil {
			return true
		}
		fallthrough
	case 4:
		return r.Name == nil
	default:
		return false
	}
}

func (v *Validator) checkBlocks(r *Record) []Issue {
	var issues []Issue
	bc := *r.BlockCount
	if len(r.BlockSizes) != bc {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "blockSizes",
			Severity: SeverityError,
			Msg:      fmt.Sprintf("blockSizes has %d entries; blockCount is %d", len(r.BlockSizes), bc),
		})
	}
	if len(r.BlockStarts) != bc {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "blockStarts",
			Severity: SeverityError,
			Msg:      fmt.Sprintf("blockStarts has %d entries; blockCount is %d", len(r.BlockStarts), bc),
		})
	}
	if bc == 0 {
		return issues
	}
	// first blockStart must be 0
	if r.BlockStarts[0] != 0 {
		issues = append(issues, Issue{
			Line:     r.Line,
			Field:    "blockStarts",
			Severity: SeverityError,
			Msg:      fmt.Sprintf("first blockStart must be 0; got %d", r.BlockStarts[0]),
		})
	}
	// each blockSize > 0
	for i, sz := range r.BlockSizes {
		if sz <= 0 {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "blockSizes",
				Severity: SeverityError,
				Msg:      fmt.Sprintf("blockSize[%d] = %d; must be > 0", i, sz),
			})
		}
	}
	// no overlapping blocks
	for i := 1; i < len(r.BlockStarts); i++ {
		prev := r.BlockStarts[i-1] + r.BlockSizes[i-1]
		cur := r.BlockStarts[i]
		if cur < prev {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "blockStarts",
				Severity: SeverityError,
				Msg: fmt.Sprintf("block[%d] starts at %d which overlaps block[%d] ending at %d",
					i, cur, i-1, prev),
			})
		}
	}
	// sum invariant: lastStart + lastSize == chromEnd - chromStart
	if len(r.BlockSizes) > 0 {
		regionLen := r.ChromEnd - r.ChromStart
		lastStart := r.BlockStarts[len(r.BlockStarts)-1]
		lastSize := r.BlockSizes[len(r.BlockSizes)-1]
		if lastStart+lastSize != regionLen {
			issues = append(issues, Issue{
				Line:     r.Line,
				Field:    "blockStarts",
				Severity: SeverityError,
				Msg: fmt.Sprintf("last block (start=%d, size=%d) ends at %d; expected chromEnd-chromStart=%d",
					lastStart, lastSize, lastStart+lastSize, regionLen),
			})
		}
	}
	return issues
}

// checkSortOrder verifies records are sorted by (chrom, chromStart).
func (v *Validator) checkSortOrder(recs []*Record) []Issue {
	var issues []Issue
	sorted := make([]*Record, len(recs))
	copy(sorted, recs)
	sort.SliceStable(sorted, func(i, j int) bool {
		if sorted[i].Chrom != sorted[j].Chrom {
			return sorted[i].Chrom < sorted[j].Chrom
		}
		return sorted[i].ChromStart < sorted[j].ChromStart
	})
	for i := range recs {
		if recs[i] != sorted[i] {
			issues = append(issues, Issue{
				Line:     recs[i].Line,
				Severity: SeverityWarning,
				Msg:      "records are not sorted by (chrom, chromStart)",
			})
			// Report just the first out-of-order record; the rest
			// would be noise.
			break
		}
	}
	return issues
}
