# Changelog

All notable changes to `htsfile` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **FASTA format support.** A new `fasta` package
  (`github.com/inful/htsfile/fasta`) implements parsing,
  writing, and validation for the FASTA sequence format.
  Streaming `Decoder`, round-trippable `Encoder` with
  configurable sequence-line wrapping (default 80 chars), and
  a per-record `Validator` with `ModeDefault` (non-empty ID,
  non-empty sequence) and `ModeStrict` (alphabet check
  against a configurable character set, default IUPAC
  nucleotide codes). The parser accepts any sequence-line
  wrapping (including unwrapped single-line records) and
  preserves the original character case; the encoder re-wraps
  at the configured width.

  The companion `fasta/fai` package reads `.fai` index files
  (read-only, per the open-decision recommendation): 5
  tab-separated fields per reference (name, length, offset,
  lineBases, lineWidth). Index generation is out of scope
  for v1 — it requires a multi-pass atomic-rename pipeline
  that's a separate effort.

  The `cmd/htsfile` CLI auto-dispatches by extension (`.fa`,
  `.fasta`, `.fna`, `.ffn`, `.faa`) and sniffs the first
  non-empty line for stdin / unknown extensions (a `>` header
  line indicates FASTA). `check`, `roundtrip`, and `summary`
  all work for FASTA.

- **CLI duplication refactor — shared issue reporter.** The
  per-format `runCheck{VCF,BED,SAM,BAM,GFF}` functions had
  grown into five near-duplicates of the same sort/print/count
  loop, each carrying a `//nolint:dupl` with rationale. FASTA
  was the ROADMAP-pinned trigger to refactor this. The
  duplication is now replaced by a single `runIssueReport`
  helper plus per-package adapters (`reportFrom{VCF,BED,SAM,
  BAM,GFF,FASTA}`) that convert native `Issue` types to a
  normalised `reportIssue` struct. The five `runCheck*`
  functions shrink to their essential ~20 lines (parse +
  validate + delegate). FASTA's `runCheckFASTA` is the first
  caller that never needed a `//nolint:dupl`. This is the
  refactor the ROADMAP flagged at the ≥3-packages threshold;
  it now lands for the sixth format.

- **GFF3 format support.** A new `gff` package
  (`github.com/inful/htsfile/gff`) implements parsing, writing,
  and validation for the Generic Feature Format Version 3.
  Streaming `Decoder`, round-trippable `Encoder`, and a
  per-record `Validator` with `ModeDefault` and `ModeStrict`.
  Supports the full 9-column feature record (SeqID, Source,
  Type, Start, End, Score, Strand, Phase, Attributes) plus the
  `##` pragma directives (`##gff-version`, `##sequence-region`,
  `##date`, and arbitrary tool-specific pragmas preserved
  verbatim for round-trip) and the `###` record terminator.
  GFF3 coordinates are 1-based inclusive (the opposite of
  BED's 0-based half-open); the parser rejects `start < 1` and
  `end < start`. Strand accepts `+`, `-`, `.`, `?`; phase
  accepts `0`, `1`, `2`, `.`. Multi-value attributes
  (`Parent`, `Dbxref`, `Alias`, `Note`, `Target`, `Gap`,
  `Derives_from`, `Ontology_term`) are exposed as
  `map[string][]string` via the `AttributeList(key)` accessor;
  `Attribute(key)` returns the first value as a convenience.
  The validator's strict mode cross-references `Parent` IDs
  against the document's `ID` set. FASTA emission after `###`
  is intentionally out of scope for v1 (the `###` line just
  ends the feature section). The `cmd/htsfile` CLI
  auto-dispatches by extension (`.gff`, `.gff3`) and sniffs
  `##gff-version` or a 9-column tab-separated feature line
  with single-character strand for stdin / unknown
  extensions.

  This is the smallest of the v1 packages — similar size to
  BED. Round-trip preserves all `##` pragmas in source order
  plus the bare `#` comments; multi-value attributes are
  serialised back as `key=v1,v2,v3`.

  CLI duplication for the four per-format `runCheck*` /
  `runSummary*` entry points has now hit the ≥3-packages
  threshold (VCF/BED/SAM + BAM + GFF = five per-format
  duplicates). Each duplicate carries a `//nolint:dupl` with
  rationale. The refactor to extract a shared
  "issue-reporting" helper is queued for FASTA, where the
  refactor's payoff will be largest; the current pattern is
  tolerable but acknowledged as a known wart.

- **BAM v1 format support.** A new `bam` package
  (`github.com/inful/htsfile/bam`) reads BGZF-compressed BAM
  files: streaming `Decoder` plus a per-record `Validator` with
  `ModeDefault` and `ModeStrict`. Supports the full binary
  framing (4-byte magic `BAM\1`, length-prefixed SAM header
  text, binary reference-sequence dictionary) and every BAM
  alignment-record field — refID, pos (0-based → 1-based POS),
  l_read_name, mapq, bin, n_cigar_op, flag, l_seq, next_refID,
  next_pos, tlen, read_name, packed CIGAR, 4-bit packed SEQ,
  qual, and the full aux-tag set (A, c, C, s, S, i, I, f, d, Z,
  H, B). SEQ is decoded back to ASCII; CIGAR is exposed as a
  packed `[]uint32` and as a SAM-text `CigarString()`; tags are
  parsed for all 12 documented BAM type codes with values
  exposed as SAM-compatible strings. BAM v1 is **read-only**
  in this build — there is no BAM writer. The companion
  `bam/bai` package reads BAI index files (read-only; bin
  index + 16384-bp linear intervals + `n_no_coor`) and exposes
  the samtools 6-level bin-assignment helper `Reg2Bins` for
  future region-query code.

  The `internal/bgzf` package is a small zero-dependency
  implementation of the BGZF block format: every block is a
  valid gzip stream with the FEXTRA flag set, carrying the
  `BC` subfield (BSIZE-1, uint16 LE), constrained to 64 KiB;
  each block is terminated with a 4-byte CRC32 + 4-byte ISIZE
  trailer. The reader decompresses blocks on the fly, validates
  CRC32 and ISIZE, and recognises the 28-byte empty EOF marker.

  The `cmd/htsfile` CLI auto-dispatches by extension (`.vcf`,
  `.bed`, `.sam`, `.bam`) and sniffs the first 4 bytes for
  stdin / unknown extensions — the gzip+FEXTRA prefix indicates
  a BGZF-compressed BAM. `check` and `summary` work for BAM;
  `roundtrip` is not supported (BAM v1 is read-only) and returns
  a clear "BAM v1 is read-only" error rather than silently
  emitting broken output.

  CLI duplication for the four per-format `runCheck*` /
  `runSummary*` entry points has hit the ≥3-packages threshold
  the ROADMAP flagged; the four duplicates each carry a
  `//nolint:dupl` with rationale. The refactor to extract a
  shared "issue-reporting" helper is queued for the next format
  landing (GFF3), where the savings will be more substantial.

- **SAM v1 format support.** A new `sam` package
  (`github.com/inful/htsfile/sam`) mirrors the vcf/ and bed/
  API: streaming `Decoder`, round-trippable `Encoder`, and a
  per-record `Validator` with `ModeDefault` and `ModeStrict`.
  Supports the full SAM header line set (`@HD`, `@SQ`, `@RG`,
  `@PG`, `@CO`, plus unknown `@XX` lines preserved verbatim for
  round-trip) and 11-column alignment records with optional
  `TAG:TYPE:VALUE` fields. Tag values are stored as raw strings
  in v1; typed accessors are deferred until a real use case
  surfaces (per `pref-similar-shape-not-duplication`). The
  `cmd/htsfile` CLI auto-dispatches by extension (`.vcf`, `.bed`,
  `.sam`) and sniffs the first line for stdin or unknown
  extensions.

  This is the largest of the v1 packages — ~3× the BED v1 work.
  CIGAR string *semantic* parsing is intentionally out of scope
  (CIGAR is stored as a string); index-file support (BAI) is
  deferred to BAM in Tier 2.

- **VCF v4.1–v4.5 acceptance.** The parser, writer, and validator
  now accept any v4.x minor version of the VCF specification.
  Previously documented as "v4.5 only", the underlying
  `ParseFileFormat` has always accepted any v4.x; this change
  makes the behaviour explicit, fixes the misleading doc comment,
  and adds explicit test coverage for v4.0, v4.1, v4.2, v4.3, and
  v4.5. Major versions other than 4 (v3, v5+) are still rejected.

  Per-version validation rule differences (which reserved keys
  are recognised, which `Number=` values are valid, etc.) are
  **not** yet implemented — the library treats all v4.x
  identically and preserves the `FileFormat` verbatim on
  round-trip. Adding per-version rules would require reading the
  v4.1–v4.4 specs to enumerate the differences, and is deferred
  until users report concrete per-version bugs.

- **BED format support.** A new `bed` package
  (`github.com/inful/htsfile/bed`) mirrors the VCF API:
  streaming `Decoder`, round-trippable `Encoder`, and a
  per-record `Validator` with `ModeDefault` and `ModeStrict`.
  Supports BED3 through BED12 with auto-detected column counts,
  optional track lines (preserved verbatim for round-trip), and
  the standard BED12 block invariants (blockCount ↔ array
  length, first blockStart is 0, sum invariant, no overlaps, all
  blockSizes > 0). Browser lines (`browser position ...`) are
  skipped during parse and not preserved. The `cmd/htsfile` CLI
  auto-dispatches by extension (`.vcf` / `.bed`) and sniffs the
  first line for stdin or unknown extensions.

  This is a sibling format to VCF; it does not share any code
  with the VCF package. The new format is intentionally narrow
  — BEDGraph, narrowPeak, broadPeak, and GFF/GTF are out of
  scope for this release.

- **Test coverage on structured-meta parse paths.** Pinned the
  duplicate INFO / FILTER / FORMAT ID detection,
  Flag-with-non-zero-Number rejection, contig length parsing,
  missing-ID error reporting for every meta-line type, and the
  positive PEDIGREE `Name_*` and SAMPLE `Genomes`/`MolConvs`
  nested-body paths. The VCF package coverage rose from 81.7% to
  85.2%.

### Changed

- **Bounded per-line read size in all three Decoders.** Each
  format package (`vcf/`, `bed/`, `sam/`) now configures its line
  reader with a `MaxLineLength` cap (default 10 MB). Lines
  longer than the cap produce a clear `*ParseError` rather than
  unbounded allocation. This is a real DoS-surface fix for
  library use cases (e.g. parsing untrusted VCF/BED/SAM input);
  the previous `bufio.Reader.ReadBytes` would happily allocate
  a 10 GB byte slice for a malicious 10 GB line. The internal
  reader is now `bufio.Scanner` with the buffer cap pinned via
  `make([]byte, 0, maxLineLength)` so the cap is enforced (a
  subtle Scanner-API gotcha: passing `make([]byte, 64*1024, 100)`
  to `Buffer` would silently grow to 64 KB regardless of the
  max). The new constructors are `NewDecoderWithMaxLine(r,
  maxLineLength)`; the existing `NewDecoder(r)` uses the
  default. Behaviour is unchanged for any real-world file.

- **Internal refactor — value-type coercion.** `InfoValue` and
  `FieldValue` now embed an unexported `typedCore` so the four
  near-identical `FieldType` switches (decode/encode ×
  INFO/FORMAT) collapse to two: `decodeTyped` in `parse.go`
  and `writeTyped` in `write.go`. Type-coercion rules are now
  defined in one place; a future change (e.g. how `NaN` is
  formatted, how percent-encoding round-trips a new edge case)
  is fixed once instead of four times.

  **API impact:** because Go does not allow promoted fields in
  struct literals, callers that construct `InfoValue` /
  `FieldValue` directly via `{Type: X, Strs: Y}` must switch to
  `{typedCore: typedCore{...}}`. As `typedCore` is unexported,
  external callers use field-by-field assignment:

  ```go
  v := vcf.InfoValue{}
  v.Type  = vcf.TypeInteger
  v.Ints  = []int64{42}
  ```

  Reads (`v.Ints`, `v.HasValue()`, `v.Int()`, `v.Float()`,
  `v.String()`) are unchanged — field promotion handles them
  transparently. The godoc example in `ExampleNewEncoder` was
  updated to reflect the new construction style. Behaviour is
  byte-for-byte identical: the corpus tests and the CLI's
  `roundtrip` subcommand both pass unchanged.

- **SAM parser errors on malformed header type tags and duplicate
  `@HD`.** A line like `@H	random` (1-character type tag) was
  previously silently stored in `Header.COs`; the parser now
  requires the type tag to be exactly two characters. A second
  `@HD` line in the same file (spec allows exactly one) now
  errors with a clear message; the first `@HD` is still
  returned alongside the error so callers can inspect what
  parsed.

## [0.0.1] — 2026-10-13

### Added

- Initial release of the VCF v4.5 parser, writer, and validator.
- Streaming `Decoder` (line-by-line) for `io.Reader` input.
- Round-trippable `Encoder` for LF and CR+LF output.
- `Validator` with two modes:
  - `ModeDefault` — well-formed + semantically consistent (header-declared
    Number/Type checks, declared keys, declared FILTER codes, sample-count
    matching, dropped-trailing-FORMAT semantics).
  - `ModeStrict` — adds spec-compliance rules (POS monotonic within CHROM,
    REF must contain only ACGTN, ALT ≠ REF, fileformat major = 4).
- `cmd/htsfile` CLI: `check`, `roundtrip`, `summary`, `version`.
- Coverage of the core v4.5 grammar:
  - All fixed columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
  - FORMAT column and per-sample data
  - Header meta lines (##fileformat, ##fileDate, ##source, ##reference,
    ##contig, ##INFO, ##FILTER, ##FORMAT, ##ALT, ##PEDIGREE, ##SAMPLE,
    ##phasing, ##assembly, ##pedigreeDB)
  - Reserved INFO keys (DP, AF, AN, NS, DB, H2, …) including SV keys
    (SVLEN, CIPOS, CIEND, IMPRECISE, EVENT, MATEID, SVCLAIM, …)
  - Reserved FORMAT keys (GT, GQ, PL, GL, GP, AD, DP, FT, PS, …)
  - Symbolic structural-variant alleles (DEL, INS, DUP, INV, CNV and
    subtypes CNV:TR, DUP:TANDEM, DEL:ME, INS:ME)
  - Breakend replacement strings (`t]p]`, `t[p[`, `]p]t`, `[p[t`)
- Standard library only; no runtime dependencies.
- CI: GitHub Actions matrix (Go 1.22 / 1.24 / stable × linux / macos) plus
  golangci-lint v2.
- Release: GoReleaser v2 producing linux/darwin × amd64/arm64 archives.

### Out of scope (deferred)

- Base modification fields (`M*`).
- Tandem-repeat fields (`RN` / `RUS` / `RUL` / `RUC` / `RB` / `CIRUC` / `CIRB`
  / `RUB`).
- Local-allele encoding (`LAA` / `LAD` / `LGL` / `LPL` / `LGP` / `LPP`).
- BCF (binary counterpart).

[Unreleased]: https://github.com/inful/htsfile/compare/v0.0.1...HEAD
[0.0.1]: https://github.com/inful/htsfile/releases/tag/v0.0.1