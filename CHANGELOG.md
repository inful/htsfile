# Changelog

All notable changes to `htsfile` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] — 2026-10-18

This release closes the "close symmetry" Phase 1 described
in the pinned context: every format the library reads now has
a Writer, and the CLI roundtrip subcommand works for all of
them.

Six new commits:

- feat(bgzf): add BGZF block writer (commit `a7a2895`)
- feat(bam/bai): add BAI index writer (commit `c478b9c`)
- feat(bam): add BAM record encoder (commit `43bbe5d`)
- feat(bcf): add BCF record encoder (commit `3af191d`)
- feat(cli): enable roundtrip for BAM and BCF (commit `6098edf`)
- docs: document the writers + CLI round-trip in CHANGELOG, README, and ROADMAP

### Added

- **BGZF block writer** (`internal/bgzf.Writer`). Symmetric
  counterpart of the existing Reader. Buffers up to 64 KiB
  per block (the spec maximum), flushes on full or on Close.
  Emits the canonical 28-byte empty EOF marker. The CDATA in
  the EOF marker is the empty-deflate-stream sequence (`03 00`)
  so the marker decodes cleanly through any flate reader. The
  header layout is fixed at 18 bytes: 4 (magic + CM + FLG) + 4
  (MTIME) + 2 (XFL) + 2 (XLEN) + 2 (SI1 + SI2) + 2 (SLEN) + 2
  (BSIZE). Foundational for both BAM and BCF writers (both use
  BGZF compression).

- **BAI index writer** (`bam/bai.Writer`). Symmetric
  counterpart of the existing Reader. Takes a pre-built *Index
  and serializes; bins are emitted in ascending Number order
  for deterministic output.

- **BAM record encoder** (`bam.Encoder`). Writes BGZF-
  compressed BAM. Surfaces the `bam.Encoder` type with
  `NewEncoder`, `WriteHeader`, `WriteRecord`, `Close`. Records
  round-trip cleanly including all 11 non-B-array aux tag
  types (B arrays require a sub-type field on `Tag` that
  the struct doesn't currently carry — documented
  limitation).

- **BCF record encoder** (`bcf.Encoder`). Writes BGZF-
  compressed BCF v2.1. Records are written WITHOUT per-sample
  FORMAT data (matching the v1 BCF reader's scope). n_fmt is
  always 0; the decoder skips the FORMAT block when n_fmt == 0,
  so the produced stream round-trips cleanly through the
  existing decoder. Surfaces the `bcf.Encoder` type with
  `NewEncoder`, `WriteHeader`, `WriteRecord`, `Close`.

- **CLI round-trip for BAM and BCF.** The `htsfile roundtrip`
  subcommand now works for `.bam` and `.bcf` files (in addition
  to the previously-supported formats). BAM round-trip is
  full; BCF round-trip is partial — header and typed record
  fields preserve, per-sample FORMAT data is dropped.

### ROADMAP §Open decisions

The "Index generation — Tabix + CSI writers" item landed in
v0.4.0. The new v0.5.0 work is the symmetric counterpart of
v0.4.0: every format the library reads now has a Writer. The
BCF writer still omits per-sample FORMAT data (matching v1's
decoder scope); per-sample binary FORMAT encoding is the only
remaining reader/writer asymmetry and is the natural v2
follow-up if/when needed.

[0.5.0]: https://github.com/inful/htsfile/compare/v0.4.0...v0.5.0
[Unreleased]: https://github.com/inful/htsfile/compare/v0.5.0...HEAD

## [0.4.0] — 2026-10-18

This release closes the Tabix + CSI halves of the ROADMAP
§Open decisions "Index generation" item. Four new commits:

- feat(tabix): add Tabix index writer (commit `22ac357`)
- feat(csi): add CSI index writer (commit `996a347`)
- feat(cli): enable roundtrip for Tabix and CSI indexes (commit `9b1c2e3`)
- docs: document Tabix + CSI writers (commit `c6731bc`)

### Added

- **Tabix and CSI index writers.** Each package now has a
  symmetric `Writer` (NewWriter + Write(*Index)) alongside its
  Reader. The Writer takes an already-built `*Index` and
  serializes it to the underlying `io.Writer` in the standard
  Tabix / CSI file format. Bin numbers are emitted in
  ascending order and chunks within each bin are sorted by
  `Beg` for deterministic output (byte-for-byte round-trip
  through the existing Reader).

  This closes the Tabix + CSI halves of the ROADMAP §Open
  decisions item **"Index generation (BCF + Tabix + CSI
  writers)"**. The BCF half is still deferred — the BCF
  writer requires per-sample FORMAT data support, which v1's
  BCF reader explicitly deferred for the same reason (see
  the v0.2.0 release notes).

  The `htsfile roundtrip` subcommand now works for `.tbi` and
  `.csi` files: it reads the file and re-emits it to stdout.
  BAM v1 and BCF v1 remain read-only (`roundtrip` returns a
  clear "read-only" error for those).

[0.4.0]: https://github.com/inful/htsfile/compare/v0.3.0...v0.4.0
[Unreleased]: https://github.com/inful/htsfile/compare/v0.4.0...HEAD

## [0.3.0] — 2026-10-18

This release builds on the v0.2.0 baseline by adding the
FASTQ format (the first Tier 4 item) and bringing the
codebase into gofmt compliance.

- **FASTQ format support.** A new `fastq` package
  (`github.com/inful/htsfile/fastq`) implements parsing,
  writing, and validation for the FASTQ short-read
  sequencing format. Streaming `Decoder` plus a per-record
  `Validator` with `ModeDefault` (non-empty header,
  non-empty sequence, sequence and quality same length,
  every quality byte ≥ PhredOffset) and `ModeStrict` (the
  default checks plus IUPAC nucleotide alphabet validation
  for the sequence). Round-trippable `Encoder` with LF or
  CR+LF output.

  Records carry `Header`, `Seq`, `PlusLine`, `Qual`, and a
  per-record `PhredOffset` (33 for Sanger / Illumina 1.8+,
  64 for older Illumina 1.3–1.7). The Decoder auto-detects
  the offset from the first record's quality bytes
  (`SetPhredOffset` disables detection for callers with
  non-standard data). Solexa Solexa+64 is not supported in
  v1.

  The `cmd/htsfile` CLI auto-dispatches by extension
  (`.fq`, `.fastq`) and sniffs the first three lines for
  stdin / unknown extensions (a `@`-prefixed header line
  followed by a sequence line and a `+` line). All three
  subcommands — `check`, `roundtrip`, `summary` — work
  for FASTQ.

### Style

- **gofmt across the codebase.** Applied `gofmt -w` to the
  existing source tree. Pre-commit, 44 files were
  gofmt-dirty (mostly missing trailing newlines; some
  column-alignment drift in struct literals). Both the
  CI `gofmt -l` check and `golangci-lint run ./...` now
  pass on a clean checkout.

[0.3.0]: https://github.com/inful/htsfile/compare/v0.2.0...v0.3.0
[Unreleased]: https://github.com/inful/htsfile/compare/v0.3.0...HEAD

## [0.2.0] — 2026-10-18

This release builds on the v0.1.0 baseline by adding three
read-only companion packages: BCF v2.1 (binary VCF), Tabix v1
(generic position index for sorted bgzipped text files), and
CSI v1 (coordinate-sorted index). All three are read-only in
this build per the open-decision scope — index generation is
a v2 follow-up.

### Added

- **BCF v2.1 format support.** A new `bcf` package
  (`github.com/inful/htsfile/bcf`) reads BGZF-compressed
  Binary Call Format files. Streaming `Decoder` plus a per-
  record `Validator` with `ModeDefault` and `ModeStrict`.
  Decodes the BCF magic + length-prefixed VCF-style header
  text, the binary reference dictionary, and every
  alignment-record field — refID, pos (0-based → 1-based POS),
  the typed ID/REF/ALT, FILTER (resolved by header),
  QUAL (float64), and the typed INFO array (Integer, Float,
  Flag, Character, String) with header-declared type
  coercion. INFO keys are resolved positionally from the
  header's `##INFO` declarations (the parser preserves
  `InfosOrder` for this). `LookupContig`, `LookupFilter`,
  `LookupInfoIndex`, and `LookupFormatIndex` provide direct
  ID → index access for callers that need it without
  iterating the order slice.

  **v1 is read-only.** The per-sample FORMAT column is not
  parsed (the binary FORMAT/sample format is complex enough
  to warrant a dedicated follow-up; v1 surfaces a clear
  ParseError rather than guess). The companion BCF *writer*
  is deferred to v2 per the open-decision scope in the
  ROADMAP — building it now would roughly double the package
  size, and the round-trip story is incomplete without sample
  data. The `cmd/htsfile` CLI auto-dispatches by extension
  (`.bcf`) and surfaces a clear "BCF v1 is read-only" error
  for `roundtrip` on `.bcf` files. `check` and `summary`
  work for BCF.

  BCF and BAM share the same outer BGZF wrapper, so the CLI
  content-sniff heuristic can't disambiguate them on stdin;
  users should pass `.bcf` files via the file path (which
  triggers the extension dispatch) or feed BCF on stdin
  knowing it will be misclassified as BAM and produce a
  parsing error.

- **Tabix v1 index reader.** A new `tabix` package
  (`github.com/inful/htsfile/tabix`) reads the general-purpose
  position index for sorted bgzipped text files (file
  extension `.tbi`). Decodes the 4-byte magic (`TBI\x01`),
  the per-reference bin index (samtools 5-level binning with
  bin numbers 0-4680), and the linear index of BGZF virtual
  offsets. Convenience helpers `BinFromRange(beg, end)` and
  `LinearIndexOffset(pos)` cover the common caller-side
  workflow of "which chunks might contain my region?".
  **Read-only in v1.** The optional `n_no_coor` field is
  always read; EOF acts as the absence signal for single-ref
  legacy files. Multi-reference legacy files without
  `n_no_coor` (pre-samtools-0.1.18) will misalign — those
  should use CSI v1 instead.

- **CSI v1 index reader.** A new `csi` package
  (`github.com/inful/htsfile/csi`) reads the Coordinate-
  Sorted Index format (file extension `.csi`) used by
  samtools ≥ 1.10. Adds an explicit version marker (`CSI\x01`)
  plus tunable `min_shift` (linear-index window size exponent)
  and `depth` (number of bin levels). The format is otherwise
  similar to Tabix. Unlike Tabix, CSI has no optional
  `n_no_coor` field — the explicit version byte is enough to
  unambiguously parse the binary format. **Read-only in v1.**
  CSI is the recommended choice for new code (use Tabix only
  when reading pre-existing `.tbi` files from samtools
  0.1.18-1.9).

- **CLI integration for Tabix and CSI.** The `cmd/htsfile`
  CLI auto-dispatches by extension (`.tbi`, `.csi`) and
  surfaces a clear "Tabix is read-only" / "CSI is read-only"
  error for `roundtrip`. `check` and `summary` work normally
  for both formats. (Tabix and CSI are not auto-sniffed from
  stdin — use the file extension to disambiguate.)

### Test coverage notes

`bcf/` rose from 65.6% to 79.0% with the addition of tests
for the previously-untested `InfoValue.String()` (each type
branch), `Severity.String()` (each branch), `sortedKeys`
(the reserved-for-future-tests helper), `skipFmtBlock` (the
always-error placeholder for FORMAT/sample data), and
`skipTypedInfoValue` (each type-byte branch plus two
truncated-stream error paths). Plus four small wins on
previously-partial paths: `Validator.checkOne` strict-mode
INFO-not-declared, `Header.LookupFormatIndex` nil receiver,
`Decoder.Header` cached-return on second call, and
`Decoder.Record` implicit-Header invocation.

[0.2.0]: https://github.com/inful/htsfile/compare/v0.1.0...v0.2.0
[Unreleased]: https://github.com/inful/htsfile/compare/v0.2.0...HEAD

## [0.1.0] — 2026-10-18

This release builds on the v0.0.1 VCF-only baseline with five
new format packages (SAM v1, BED v1, BAM v1, GFF3, FASTA),
VCF v4.x acceptance, the CLI duplication refactor that
eliminates the five-way issue-reporting duplication, internal
value-type coercion, bounded per-line reads, and SAM
header-tag validation.

### Added

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

- **BAM v1 format support.** A new `bam` package
  (`github.com/inful/htsfile/bam`) reads BGZF-compressed BAM
  files: streaming `Decoder` plus a per-record `Validator`
  with `ModeDefault` and `ModeStrict`. Supports the full
  binary framing (4-byte magic `BAM\1`, length-prefixed SAM
  header text, binary reference-sequence dictionary) and every
  BAM alignment-record field — refID, pos (0-based → 1-based
  POS), l_read_name, mapq, bin, n_cigar_op, flag, l_seq,
  next_refID, next_pos, tlen, read_name, packed CIGAR,
  4-bit packed SEQ, qual, and the full aux-tag set (A, c, C,
  s, S, i, I, f, d, Z, H, B). SEQ is decoded back to ASCII;
  CIGAR is exposed as a packed `[]uint32` and as a SAM-text
  `CigarString()`; tags are parsed for all 12 documented BAM
  type codes with values exposed as SAM-compatible strings.
  BAM v1 is **read-only** in this build — there is no BAM
  writer. The companion `bam/bai` package reads BAI index
  files (read-only; bin index + 16384-bp linear intervals +
  `n_no_coor`) and exposes the samtools 6-level
  bin-assignment helper `Reg2Bins` for future region-query
  code.

  The `internal/bgzf` package is a small zero-dependency
  implementation of the BGZF block format: every block is a
  valid gzip stream with the FEXTRA flag set, carrying the
  `BC` subfield (BSIZE-1, uint16 LE), constrained to 64 KiB;
  each block is terminated with a 4-byte CRC32 + 4-byte ISIZE
  trailer. The reader decompresses blocks on the fly, validates
  CRC32 and ISIZE, and recognises the 28-byte empty EOF marker.

  The `cmd/htsfile` CLI auto-dispatches by extension (`.vcf`,
  `.bed`, `.sam`, `.bam`) and sniffs the first 4 bytes for
  stdin / unknown extensions — the gzip+FEXTRA prefix
  indicates a BGZF-compressed BAM. `check` and `summary` work
  for BAM; `roundtrip` is not supported (BAM v1 is read-only)
  and returns a clear "BAM v1 is read-only" error rather than
  silently emitting broken output.

- **GFF3 format support.** A new `gff` package
  (`github.com/inful/htsfile/gff`) implements parsing,
  writing, and validation for the Generic Feature Format
  Version 3. Streaming `Decoder`, round-trippable `Encoder`,
  and a per-record `Validator` with `ModeDefault` and
  `ModeStrict`. Supports the full 9-column feature record
  (SeqID, Source, Type, Start, End, Score, Strand, Phase,
  Attributes) plus the `##` pragma directives
  (`##gff-version`, `##sequence-region`, `##date`, and
  arbitrary tool-specific pragmas preserved verbatim for
  round-trip) and the `###` record terminator. GFF3
  coordinates are 1-based inclusive (the opposite of BED's
  0-based half-open); the parser rejects `start < 1` and
  `end < start`. Strand accepts `+`, `-`, `.`, `?`; phase
  accepts `0`, `1`, `2`, `.`. Multi-value attributes
  (`Parent`, `Dbxref`, `Alias`, `Note`, `Target`, `Gap`,
  `Derives_from`, `Ontology_term`) are exposed as
  `map[string][]string` via the `AttributeList(key)` accessor;
  `Attribute(key)` returns the first value as a convenience.
  The validator's strict mode cross-references `Parent` IDs
  against the document's `ID` set. FASTA emission after `###`
  is intentionally out of scope for v1. The
  `cmd/htsfile` CLI auto-dispatches by extension (`.gff`,
  `.gff3`) and sniffs `##gff-version` or a 9-column
  tab-separated feature line with single-character strand for
  stdin / unknown extensions.

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
  was the ROADMAP-pinned trigger to refactor. The duplication
  is now replaced by a single `runIssueReport` helper plus
  per-package adapters (`reportFrom{VCF,BED,SAM,BAM,GFF,
  FASTA}`) that convert native `Issue` types to a normalised
  `reportIssue` struct. The five `runCheck*` functions shrink
  to their essential ~20 lines (parse + validate + delegate).
  FASTA's `runCheckFASTA` is the first caller that never
  needed a `//nolint:dupl`. This is the refactor the ROADMAP
  flagged at the ≥3-packages threshold; it now lands for the
  sixth format.

### Changed

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

- **Bounded per-line read size in all three Decoders.** Each
  format package (`vcf/`, `bed/`, `sam`) now configures its line
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

### Test coverage notes

- **BAM `checkOne` strict-mode branches.** Reached 100%
  coverage on the validator's per-record check with the
  addition of four tests for the previously-untested
  strict-mode branches (QUAL length mismatch, strict-mode
  POS check, strict-mode duplicate-tag warning, strict-mode
  next_refID out-of-bounds). The `bam/` package overall
  coverage rose from 76.2% to 78.8%.

- **VCF structured-meta parse paths.** Pinned the
  duplicate INFO / FILTER / FORMAT ID detection,
  Flag-with-non-zero-Number rejection, contig length parsing,
  missing-ID error reporting for every meta-line type, and the
  positive PEDIGREE `Name_*` and SAMPLE `Genomes`/`MolConvs`
  nested-body paths. The VCF package coverage rose from 81.7%
  to 85.2%.

[0.1.0]: https://github.com/inful/htsfile/compare/v0.0.1...v0.1.0

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

[0.0.1]: https://github.com/inful/htsfile/releases/tag/v0.0.1