// Package vcf implements a parser, writer, and validator for the Variant Call
// Format version 4.5 (VCFv4.5) as specified in
// https://github.com/samtools/hts-specs/blob/master/VCFv4.5.tex.
//
// The package exposes three primary types:
//
//   - [Decoder] for streaming VCF input into a [Header] and a sequence of
//     [Record] values.
//   - [Encoder] for serialising a [Header] and [Record] values back into a
//     VCF-conforming text stream.
//   - [Validator] for cross-checking that records agree with the metadata
//     declared in the header (cardinality, type, key reference) and, when
//     configured for [Strict] mode, additional rules mandated by the spec.
//
// The implementation is intentionally dependency-light: only the Go standard
// library is used.
//
// # Streaming
//
// A [Decoder] reads one line at a time from its underlying [io.Reader]; the
// whole file is never buffered in memory. [Encoder] similarly writes a single
// line per call to [Encoder.EncodeRecord].
//
// # Round-trip
//
// Reading a VCF with [Decoder] and writing the parsed [Header] and
// [Record] values back with [Encoder] is guaranteed to produce a stream that
// parses again into equivalent values (same field content, same logical
// metadata). The textual representation may differ in field ordering inside
// structured meta-information lines (the spec explicitly allows this), but
// field content is preserved.
//
// # Coverage
//
// Currently supported (per the chosen scope):
//
//   - File header: ##fileformat, ##fileDate, ##source, ##reference, ##contig,
//     ##INFO, ##FILTER, ##FORMAT, ##ALT, ##assembly, ##pedigreeDB,
//     ##PEDIGREE, ##META, ##SAMPLE, and other unstructured ##key=value lines.
//   - 8 fixed columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) and the
//     optional FORMAT column with one data field per sample.
//   - Percent-encoding of reserved characters in values.
//   - Reserved INFO keys (AA, AC, AD, ADF, ADR, AF, AN, BQ, CIGAR, DB, DP,
//     END, H2, H3, MQ, MQ0, SB, SOMATIC, VALIDATED, 1000G) and SV keys
//     (IMPRECISE, NOVEL, SVLEN, CIPOS, CIEND, HOMLEN, HOMSEQ, BKPTID, MEINFO,
//     METRANS, DGVID, DBVARID, DBRIPID, MATEID, PARID, EVENT, EVENTTYPE,
//     CILEN, CN, CICN, SVCLAIM, RN, RUS, RUL, RUC, RB, CIRUC, CIRB, RUB).
//     v4.5-specific reserved keys (tandem repeats, base modifications) are
//     looked up by [IsReservedInfoKey]; the validator checks declared
//     Number/Type against the spec table in [ModeStrict] mode.
//   - Reserved FORMAT keys (GT, DP, EC, FT, GL, GP, GQ, HQ, LEN, MQ, PL, PP,
//     PQ, PS, AD, ADF, ADR) and SV keys (CN, CICN, CNQ, CNL, CNP, NQ, HAP,
//     AHAP). v4.5 local-allele FORMAT keys (LAA, LAD, LADF, LADR, LEC, LGL,
//     LGP, LPL, LPP) are looked up by [IsReservedFormatKey]; LAA is also
//     structurally validated to follow GT in the FORMAT column.
//   - Symbolic structural variant alleles: DEL, INS, DUP, INV, CNV (and the
//     recommended subtypes CNV:TR, DUP:TANDEM, DEL:ME, INS:ME) plus the
//     unspecified-allele <*>.
//   - Breakend replacement strings (t]p], t[p[, ]p]t, [p[t), including
//     contig-anchored and circular-contig forms.
//   - The phasing field PS (phase set identifier) is acknowledged.
//
// Currently out of scope: full structural validation of the
// tandem-repeat and base-modification record values
// (e.g. enforcing that RB is a repeat motif present in
// REF, or that M0A's value count matches the per-base
// cardinality). The library recognises the keys (see
// [IsReservedInfoKey], [IsBaseModKey]) and warns on
// mismatched header metadata; the value-level grammar is
// the user's responsibility.
//
// # Safety
//
// The package never panics on malformed input; all errors are returned as
// [ParseError] / [ValidationError] values with line and (where relevant)
// record context.
package vcf
