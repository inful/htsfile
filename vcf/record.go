package vcf

// Record is one data line of a VCF file — the 8 fixed fields plus, when
// present, a FORMAT column and one column per sample.
//
// The zero value is not meaningful; obtain Records from [Decoder.Record].
type Record struct {
	// Chrom is the chromosome / contig name (column 1). It must match a
	// declared ##contig ID when the header is present.
	Chrom string

	// Pos is the 1-based reference position (column 2).
	//
	// Telomeres are represented by Pos = 0 or Pos = N+1, where N is the
	// length of the matching contig.
	Pos int64

	// ID is the semicolon-separated identifier list (column 3). Use "." or
	// an empty slice to indicate that no ID is available.
	ID []string

	// Ref is the reference allele(s) (column 4). Each base must be one of
	// A, C, G, T, N (case insensitive).
	Ref string

	// Alt is the list of alternate allele strings (column 5), one entry
	// per comma-separated value. Each entry is either a nucleotide string,
	// the "*" symbol, a symbolic <ID>, the <* > unspecified allele, or a
	// breakend replacement string (t]p], t[p[, ]p]t, [p[t).
	Alt []string

	// Qual is the Phred-scaled quality score (column 6). A nil pointer
	// represents the MISSING value.
	Qual *float64

	// Filter is the list of filter codes (column 7). A nil slice represents
	// PASS (the single token "PASS"). An empty slice followed by a "." entry
	// would be unusual; PASS is the common case.
	//
	// To represent "filters have not been applied", set Filter to nil — by
	// the spec's writing convention, the literal "." should be encoded as a
	// single-entry slice containing "." so that round-trip preserves intent.
	Filter []string

	// Info is the per-record INFO field map (column 8). The presence of a
	// key implies a non-missing value (Flag = true for [TypeFlag] entries).
	// To record a MISSING value for a particular INFO key, set
	// Info[key].Missing = true.
	Info map[string]InfoValue

	// Format is the ordered list of FORMAT column keys (column 9). It is
	// nil for records that have no sample columns.
	Format []string

	// Samples is the per-sample data, in the same order as
	// [Header.SampleNames]. It is nil when no FORMAT column is present.
	Samples []Genotype
}

// Genotype is the per-sample column block for one record.
//
// Fields maps the FORMAT key to its typed value for this sample. Keys
// present in [Record.Format] but absent from Fields are treated as
// "trailing fields dropped" by the spec and are valid.
//
// Samples may have heterogeneous ploidy; downstream code can inspect
// Genotype.Fields["GT"] for the GT string and re-parse it.
type Genotype struct {
	Fields map[string]FieldValue
}
