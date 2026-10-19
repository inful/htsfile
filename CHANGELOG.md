# Changelog

All notable changes to `htsfile` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- *(nothing yet)*

## [0.15.8] — 2026-10-19

Patch release closing four smaller spec gaps surfaced
by the v0.15.7 audit.

### Added — vcf

- **`Validator.checkIntegerBCFRange`** (T4.14). Per
  `§1.4.2`, the eight Integer values from $-2^{31}$ to
  $-2^{31}+7$ are reserved by BCF2's type encoding and
  therefore disallowed in VCF too. Catching them at the
  VCF layer prevents silent corruption on a VCF→BCF→VCF
  round-trip. Walks `rec.Info` and `rec.Samples[*].Fields`.
- **mtDNA haploid GT check** (T4.15). Per `§1.6.2`,
  mitochondrial GT must have exactly one allele. Y is
  skipped because PAR regions legitimately carry diploid
  calls. 4 sub-tests.
- **PEDIGREE Father/Mother/Original resolution** (T4.16).
  Header-only ModeStrict check: PEDIGREE references must
  resolve to a declared sample or another PEDIGREE ID.

### Changed — vcf

- **`parseRecord`** now enforces the FORMAT key regex
  `^([A-Za-z_][0-9A-Za-z_.]*|1000G)$` per-record (T4.13,
  `§1.6.2`). The header-side check already enforced it
  for declared keys.

### Tests

- `TestValidator_IntegerBCFRange` with 5 sub-tests.
- `TestValidator_HaploidGTForMtDNA` with 4 sub-tests.
- `TestValidator_PedigreeReference` with 3 sub-tests.
- `TestDecoder_Errors` gained 2 new FORMAT-key-regex
  cases.

`go test -race -count=1 ./...` is green; `golangci-lint run ./...` reports 0 issues; `make fmt-check` is clean.

## [0.15.7] — 2026-10-19

Patch release closing two record-level spec gaps
surfaced by a fresh VCF v4.5 / BCF 2.2 audit, plus the
audit itself in `SPEC_COMPLIANCE.md`.

### Added — vcf

- **`Validator.CIRUC/CIRB bounds-span-zero`** (T4.12).
  The new check rejects records where a CIRUC or CIRB
  lower bound is positive or upper bound is negative,
  per `§1.6.1.7`.
- **`SPEC_COMPLIANCE.md` "Spec audit" section**:
  full section-by-section status against VCF v4.5 / BCF
  2.2, with ✅ / ❓ markers and line numbers in the
  .tex source, and a prioritised list of remaining
  gaps.

### Changed — vcf

- **`parseRecord`** now rejects duplicate IDs within a
  single record's ID column and duplicate keys within a
  single record's FORMAT column (T4.11, §1.6.1 line
  401 and §1.6.2 line 469).

### Tests

- `TestDecoder_Errors` gained two entries: "duplicate
  ID within record" and "duplicate FORMAT key within
  record".
- New `TestValidator_CIRUCBIRBBoundsSpanZero` with 5
  sub-tests (valid CIRUC, lower/upper-bound violations
  for each of CIRUC and CIRB).

`go test -race -count=1 ./...` is green; `golangci-lint run ./...` reports 0 issues; `make fmt-check` is clean.

## [0.15.6] — 2026-10-19

Patch release adding breakend mate position round-trip
validation.

### Added — vcf

- **`Validator.checkBreakendMatePositions`** (`599f6f0`,
  SPEC_COMPLIANCE T4.10). The two records in a
  breakend pair must round-trip — each mate's breakend
  remote position must equal the partner's local
  (Chrom, POS). Fires only when both mates carry
  breakend ALTs (paired SNVs sharing a MATEID are
  skipped). Both directions (A→B and B→A) are
  checked.

### Changed — vcf

- **`mateSnapshot`** gained a `chrom`, `pos`, and
  `breakends` field so the round-trip check has the
  data it needs from the first mate seen.

### Fixed — vcf

- `testdata/sv_example.vcf`: the `delbp1`/`delbp2`
  breakend pair had a real round-trip error (record
  1's breakend pointed to chrA:5 but the partner was
  at chrA:2). Fixed to a proper round-trip
  (`T[chrA:2[` + `]chrA:2]A`).

### Tests — vcf

- `TestValidator_BreakendMatePositionRoundTrip` with
  6 sub-tests: balanced round-trip, wrong partner
  POS, wrong partner remote, wrong remote contig,
  missing MATEID, non-breakend ALT.

`go test -race -count=1 ./...` is green; `golangci-lint run ./...` reports 0 issues; `make fmt-check` is clean.

## [0.15.5] — 2026-10-19

## [0.15.5] — 2026-10-19

Patch release adding §1.7 gVCF non-variant block
validation.

### Added — vcf

- **`Validator.checkGVCFBlock`** (`bca688d`,
  SPEC_COMPLIANCE T4.9). A non-variant block record —
  `ALT="<*>"` (parser-normalised from `<NON_REF>`) plus
  an INFO END — now gets four checks:
  1. `<*>`-ALT must declare END (Strict → error;
     default → warning).
  2. END ≥ POS (block cannot go backwards).
  3. Cross-record overlap: when the previous record
     was a block, the next record's POS must be >
     previous block's END.
  4. With a FASTA reference attached, REF must match
     the reference over the block range `[POS, END]`.

- **`prevEnd` field** on the streaming Validator
  tracks the previous block's END so the overlap
  check (rule 3) can fire without a separate pass.

### Tests

- `TestValidator_GVCFBlockRules` covers rules 1-2
  across modes.
- `TestValidator_GVCFBlockOverlap` covers rule 3.
- `TestValidator_GVCFBlockRefMatchesReference`
  covers rule 4 with a real FASTA reference.

### Docs

- `SPEC_COMPLIANCE.md` T4.9 documented; quick index
  §1.7 → `T1.6, T4.9`.

`go test -race -count=1 ./...` is green; `golangci-lint run ./...` reports 0 issues; `make fmt-check` is clean.

## [0.15.4] — 2026-10-19

Patch release adding `htsfile vcf index` and the
index-builder API.

### Added — vcf

- **`BuildIndex(path, opts)`** (`7b3d11e`). Reads a
  BGZF-compressed VCF and writes either a tabix
  `.tbi` (default) or a CSI `.csi` index next to it.
  Algorithm: the vcf.Decoder buffers blocks ahead,
  so per-record offsets inside a block are
  unreliable; we use the BGZF block boundary as the
  natural chunk boundary and emit one chunk per
  (ref, bin) per block. End is set to the start of
  the next block, or to the BGZF virtual offset of
  the EOF marker for the last block, which is what
  samtools / bcftools emit by default and what
  makes the chunk-End stop actually fire on
  single-block files. The linear index records the
  BGZF offset of the first record in each
  `(2^min_shift)`-bp window. 6 tests cover the
  happy path, both formats, options validation,
  small-region round-trip, EOF handling, and bad
  paths.

- **`IndexOpts` / `IndexFormat` types** alongside
  `BuildIndex`. `FormatTabix` is the default;
  `FormatCSI` takes `MinShift` (bin window
  exponent) and `Depth` (bin levels).

- **`TabixReader.Query` / `CSIReader.Query` CHROM
  + dedup fix** (`7b3d11e`). Two correctness
  changes needed for the index writer to round-trip:
  1. Chunks are aligned to BGZF blocks, which can
     hold records from multiple CHROMs. Queries now
     filter by `rec.Chrom == chrom` in the inner
     loop, in addition to the existing POS range
     filter.
  2. Multiple bins can cover the same record (a
     record sits in every bin its POS falls into),
     so the same record can appear in multiple
     chunks. Queries now dedup emitted records via a
     `(chromIndex, pos)` set.

### Added — cmd/htsfile

- **`htsfile vcf index`** (`7b3d11e`). Wraps
  `vcf.BuildIndex`. Flags:
  `--format=auto|tabix|csi` (default: auto →
  tabix), `--min-shift=N` (CSI bin window
  exponent, default 14), `--depth=N` (CSI bin
  levels, default 5). 4 new tests cover both
  formats, the round-trip with `vcf query`, and
  argument validation.

`go test -race -count=1 ./...` is green; `golangci-lint run ./...` reports 0 issues.

## [0.15.3] — 2026-10-19

## [0.15.3] — 2026-10-19

Patch release adding benchmark coverage for the
indexed-read paths.

### Added — internal/bgzf, vcf

- **Benchmarks** (`fe3a81e`):
  - `internal/bgzf/BenchmarkReader_Decompress`:
    raw BGZF block decompression throughput on a
    single 64 KiB block of incompressible data
    (~1.3 GB/s on M1 Pro).
  - `internal/bgzf/BenchmarkSeekReader_RandomSeek`:
    SeekReader.Seek + Read cost on a two-block
    stream.
  - `vcf/BenchmarkTabixQuery_*` /
    `BenchmarkCSIQuery_*`: cold (open per call,
    mirrors `htsfile vcf query`) and warm
    (open once) variants for both small-region
    and whole-contig queries, plus
    `BenchmarkBGZFSeekReader_VirtualOffset` for
    the per-record cost the new decodeChunk
    check adds.
  - `vcf/BenchmarkTabixQuery_ChunkEndStop`:
    A/B between a single-chunk index with
    End=MaxUint64 and the per-block index
    samtools emits. The per-block variant is
    1-2 % faster on the multi-block fixture,
    confirming the chunk-End stop is doing
    useful work for records that straddle
    block boundaries.

`go test -race -count=1 ./...` is green; `golangci-lint run ./...` reports 0 issues.

## [0.15.2] — 2026-10-19

Patch release adding the chunk-End stop
optimisation on top of v0.15.1.

### Changed — vcf, internal/bgzf

- **`decodeChunk` stops at the chunk's End BGZF
  virtual offset instead of decoding until EOF**
  (`2aa03d7`). For wide queries this stops
  fetching bytes past the last overlapping
  record. New `bgzf.Reader.VirtualOffset()`
  tracks the BGZF virtual offset of the next
  byte to return (high 48 bits = block start in
  the underlying file; low 16 bits = within-
  decompressed-block cursor), backed by a
  `countingReader` that captures the file
  position at NewReader time. `SeekReader.VirtualOffset`
  forwards it. `TabixReader.decodeChunk` and
  `CSIReader.decodeChunk` accept `chunkEndOff`
  and break the read loop when the BGZF position
  reaches it. `TestTabixReader_ChunkEndStop`
  builds a multi-block fixture and rebuilds the
  index with chunk End at the block-2 boundary
  to prove the optimisation works.

## [0.15.1] — 2026-10-19

Patch release adding the `htsfile vcf query` CLI
subcommand and a chunk-boundary stop optimisation
on top of v0.15.0.

### Added — cmd/htsfile

- **`htsfile vcf query`** (`b8f4c91`). Reads a
  tabix- or CSI-indexed BGZF-compressed VCF and writes
  the records overlapping a region to stdout in VCF
  text form. Region is samtools-style `chrom:start-end`
  (1-based inclusive). Whole-contig query (`chrom` with
  no `:start-end`) reads the contig length from the
  header. `--index=auto|tabix|csi` chooses the index
  type; the default prefers `.tbi` and falls back to
  `.csi`. 5 new tests cover region match, whole
  contig, malformed region, missing index, and the
  region parser.

### Changed — vcf, internal/bgzf

- **`decodeChunk` stops at the chunk's End BGZF
  virtual offset instead of decoding until EOF**
  (`e1a52bf`). Tabix/CSI chunks carry an explicit End
  byte offset that marks the last byte of the last
  record in the chunk; the reader now respects it.
  Two pieces were needed:

  1. **`bgzf.Reader.VirtualOffset()`** tracks the
     absolute BGZF virtual offset of the next byte to
     return, by wrapping the underlying io.Reader in
     a `countingReader` that captures the file
     position at NewReader time. `SeekReader.VirtualOffset`
     forwards it.
  2. **`TabixReader.decodeChunk` and
     `CSIReader.decodeChunk`** accept `chunkEndOff`
     and break the read loop when the BGZF position
     reaches it. Callers (`Query`) pass each chunk's
     `End` value.

  One new test (`TestTabixReader_ChunkEndStop`)
  builds a multi-block fixture, discovers the second
  block's BGZF virtual offset via the new
  `VirtualOffset()`, rebuilds the index with a chunk
  End at that boundary, and confirms Query returns
  fewer records than the full-file query.

## [0.15.0] — 2026-10-19

This release adds tabix- and CSI-indexed random access
to BGZF-compressed VCF, plus CIGAR consistency
validation. It is the single largest user-facing
feature since v0.10.0.

### Added — internal/bgzf

- **`SeekReader`** (`a3f1c08`). Random-access read of a
  BGZF-compressed stream backed by `*os.File` via
  virtual offsets (block << 16 | byte-within-block).
  Each `Seek` closes the previous decompressor and
  creates a new one positioned at the right block
  boundary; the caller can then `Read` from there.
  Used by the indexed VCF readers.

### Added — vcf

- **`vcf.OpenTabix(path)`** (`a3f1c08`). Opens a
  tabix-indexed BGZF-compressed VCF and returns a
  `*TabixReader`. The sibling `.tbi` must exist.
  - `Header()` returns the parsed header (cached).
  - `Query(chrom, start, end)` returns records
    overlapping the half-open 1-based region
    `[start, end)` in O(n) per overlap, not per
    file. Unknown contigs return `nil, nil`.
- **`vcf.OpenCSI(path)`** (`a3f1c08`). CSI sibling of
  `OpenTabix` for `.csi` indexes. Same `Header` and
  `Query` surface.
- **CIGAR consistency check** (`a3f1c08`). When
  `INFO/CIGAR` is present, the validator verifies
  the sum of REF-consuming op lengths equals
  `len(REF)` and the sum of ALT-consuming op
  lengths equals `len(ALT)`. Malformed CIGAR
  strings (missing op, unknown op) are rejected.
  7 new sub-tests cover SNV, indel, insertion,
  malformed and unknown-op cases.

9 new tests for `TabixReader` and `CSIReader` cover
open, header, region query, unknown contig, invalid
region, missing index, and BGZF-stream sanity.

## [0.14.0] — 2026-10-18

This release makes the `fasta` package substantially
faster and adds BGZF support. The vcf validator has
no API or behaviour changes.

### Added — fasta

- **BGZF support** (`8c9a44a`). `Open` now detects the
  BGZF magic (`0x1f 0x8b 0x08` with FEXTRA flag) and
  decompresses the stream into memory before indexing.
  The `.fai` index is ignored for BGZF because its
  virtual offsets require a seekable BGZF reader
  which this package does not currently provide. The
  caller pays a memory cost equal to the decompressed
  size. All other `Reference` APIs work identically.
- **`Reference` now uses `io.ReadSeeker`** instead of
  `*os.File` directly. This unifies the plain-FASTA
  and BGZF code paths.

### Changed — fasta

- **Cursor optimisation** (`8c9a44a`). When `BaseAt`
  or `BaseRange` is called sequentially on the same
  contig with increasing positions, the scanned-path
  scan resumes from a recorded file offset rather than
  restarting at the first base of the contig. This
  makes sequential VCF validation O(N) instead of
  O(N²). Cursor resets on chrom change, on backward
  seeks, and on `.fai`-path lookups.
- **Plain-FASTA scan uses `readLineTracked`** instead
  of `bufio.NewReader` so the file offset of the next
  byte is exact after each line. (The previous
  bufio-based approach made the offset depend on
  `Buffered()`, which is unreliable for cursor
  bookkeeping.)

### Not implemented

- **IUPAC ALT reduction cross-check**: the parser
  doesn't track which ALTs were originally IUPAC codes
  (declared via `##ALT`), so we cannot distinguish a
  reduced-IUPAC base from a concrete call. The full-REF
  match added in v0.11.0 already catches the case
  where the IUPAC reduction is inconsistent with the
  reference for REF; the same coverage applies to ALT
  substitutions where the ALT base replaces the REF
  base.

3 new tests cover the cursor optimisation
(sequential calls, backward seek, BGZF-decompressed
content).

## [0.13.0] — 2026-10-18

This release adds a reference-aware check for
breakend remote positions (§1.6.5.1). For every
breakend ALT in a record, the remote contig and
position are parsed and the position is verified
against the reference contig length.

### Added — vcf validator

- **`checkBreakendPositions`** (`f3a0a8c`). Fires
  only when a FASTA reference is attached. For every
  ALT matching the breakend shape `t]p]`, `t[p[`,
  `]p]t`, or `[p[t`, parse the remote position
  (`chr:pos`) and verify `pos <= ref.Length(chr)`.
  - Unknown remote contig → SeverityWarning in
    ModeDefault, SeverityError in ModeStrict.
  - Telomere breakends (no position) → skip.
  - Non-breakend ALTs → skip.
- **`parseBreakendPosition`** (exported inside
  the package). Handles the four breakend forms
  plus the optional strand-bracket trailing the
  position.

9 new sub-tests cover the four breakend forms,
within-contig, past-contig, cross-contig, telomere,
unknown-remote-contig, and non-breakend cases.

## [0.12.0] — 2026-10-18

This release adds two reference-aware bounds checks
that fire only when a FASTA reference is attached:
`END` and `SVLEN` must stay within the contig. Both
checks are skipped for mitochondrial chromosomes
(matching the existing `checkEND` exception).

### Added — fasta

- **`(*Reference).Length(chrom) (int64, error)`**
  (`7d4bbaa`). Returns the number of bases in the
  named contig, or an error if the contig is unknown.
  Used by the validator for bounds checks.

### Added — vcf validator

- **`checkENDBounds`** (`7d4bbaa`). When a reference
  is attached and the record declares `INFO END`,
  verify `END <= ref.Length` and `END >= POS`.
  Fires in both modes; replaces the need for a
  declared `##contig=<length>`.
- **`checkSVLENBounds`** (`7d4bbaa`). When a
  reference is attached and the record declares
  `INFO SVLEN`, verify `POS + |SVLEN| - 1 <=
  ref.Length` for each entry. Negative SVLEN is
  caught by the existing positive-value check.

13 new sub-tests cover the within-contig, past-
contig, before-POS, mtDNA, and unknown-contig cases.

## [0.11.0] — 2026-10-18

This release replaces the single-base REF padding
check (v0.10.0) with a full-REF match against the
reference, plus adds an efficient multi-base lookup
to the `fasta` package. The new check strictly
subsumes the old one and adds three new error classes:
wrong contig, wrong POS, and undeclared IUPAC codes.

### Added — fasta

- **`(*Reference).BaseRange(chrom, start, end)`**
  (`1ba3d28`). Half-open range lookup returning the
  bases at 1-based positions `start` through `end-1`.
  Uses a single forward scan instead of N independent
  BaseAt calls. Three new tests cover normal use,
  out-of-range, and across-comment-line.

### Changed — vcf validator

- **`Validator.WithReference` now performs full-REF
  matching** (`1ba3d28`). For every byte of REF,
  the validator compares it against the corresponding
  reference base at `POS+i`. One SeverityError per
  mismatch; lookup failure degrades to SeverityWarning
  in ModeStrict, skipped in ModeDefault. This catches:
  - Wrong contig (REF[0] doesn't match ref[POS]).
  - Wrong POS (off-by-one or larger).
  - Transcription errors anywhere in REF.
  - Undeclared IUPAC codes that didn't reduce correctly.
  - Multi-base substitution REF mismatches.

  The previous single-base padding check
  (`checkPaddingBase`) is removed; the new check
  strictly subsumes it.

### Removed

- **`checkPaddingBase` and `needsPadding`** in
  `vcf/validate.go`. The full-REF check is the single
  reference-aware REF validator.

## [0.10.0] — 2026-10-18

This release adds a top-level `fasta` package and uses
it to give the VCF validator a reference-aware padding
check (§1.6.1.1). It is a small, additive feature
release: no public API in `vcf` changed beyond the new
`Validator.WithReference` method.

### Added — new package

- **`fasta` package** (`e6f9b91`). A small
  read-only FASTA reference reader. Supports plain
  FASTA with optional samtools-faidx-style `.fai`
  random access. When no `.fai` is found, the package
  streams the file once at open time to build an
  in-memory index. The package exposes `Open`,
  `(*Reference).BaseAt`, `(*Reference).Contigs`, and
  `(*Reference).Close`. Five unit tests cover scanned,
  indexed, unknown-contig, out-of-range, and
  variable-width cases. BGZF compression is out of
  scope for this release.

### Added — vcf validator

- **`Validator.WithReference`** (`e6f9b91`). When a
  `*fasta.Reference` is attached, the validator
  enforces the §1.6.1.1 REF/ALT padding-base rule:
  for any record where REF or any ALT is symbolic, or
  any ALT has a different length from REF (a true
  indel), `REF[0]` must equal the reference base at
  `POS`. Reference lookup failures (unknown contig,
  position out of range) downgrade to a SeverityWarning
  in `ModeStrict` and are skipped in `ModeDefault`.
  Six new sub-tests cover right/wrong padding for
  insertion, deletion, symbolic ALT, SNV (no
  padding required), and unknown contig.

## [0.9.1] — 2026-10-18

This patch release closes the last open item from the
v0.9.0 spec-compliance plan: T1.3 (`##META` Values
enforcement). It is a strict-mode warning that fires
when a `##SAMPLE` field whose key matches a `##META` ID
that declared a `Values=[...]` list uses a value outside
that list. v0.9.0 marked T1.3 as deferred; this release
adds the check.

### Added — header

- **`##META` Values enforcement** (`9e7dc91`). The
  parser now recognises `##META` as a structured line
  with its own typed schema (`MetaDef{ID, Type,
  Number, Values}`); the structured-body splitter
  tracks bracket depth so that commas inside
  `[Values=[a,b,c]]` don't split fields. The validator
  fires a `SeverityWarning` in ModeStrict when a
  `##SAMPLE` field value is outside the declared list.
  ModeDefault skips the check to avoid surprising
  callers.

## [0.9.0] — 2026-10-18

This release finishes the VCF v4.5 spec-compliance work
started in v0.7.0 and v0.8.0. All 17 remaining tasks in
`vcf/SPEC_COMPLIANCE.md` are done; the only deferred item
is T1.3 (`##META` Values enforcement), which is
low-impact and can land later if needed.

### Added — parser hardening

- **`FILTER="0"` is reserved** (`d6f8384`).
  Per §1.6.1 the literal `0` cannot appear as a
  filter value; both the bare-zero and mixed-zero
  forms are rejected at parse time.

- **`##ALT` first-level type is one of
  `DEL | INS | DUP | INV | CNV`** (`ac1d4b5`).
  Per §1.4.4 the structural-variant first level must
  be one of those five tokens (or a bare IUPAC code
  — see below). Subtypes via `:` are accepted.

- **Contig ID must not be a reserved symbolic allele**
  (`e56fa7a`). `##contig=<ID=DEL,…>` and similar
  reserved symbolic-allele names are rejected per
  §1.4.7.

- **`<NON_REF>` is an alias for `<*>`** (`154f733`).
  The ALT parser normalises the gVCF alias form to the
  canonical `<*>` so downstream code only sees one
  form.

- **IUPAC ambiguity-code reduction in ALT** (`f1cb3a0`).
  Per §1.4.4 the parser reduces `R, Y, S, W, K, M, B,
  D, H, V` in ALT to their first-alphabetically
  constituent base when the corresponding `##ALT`
  declaration exists. An undeclared IUPAC code in ALT
  is rejected by the validator.

### Added — per-record checks

- **LAA required when local-allele fields present**
  (`8fb631c`). `LAD / LADF / LADR / LEC / LGL / LGP /
  LPL / LPP` without a non-empty LAA fires a
  SeverityError in ModeStrict.

- **PS and PSL mutual exclusion per sample**
  (`5e9c7c1`). Per §1.6.2 a sample carrying both `PS`
  and `PSL` (with non-empty values) fires a
  SeverityError.

- **`CIPOS / CIEND / CILEN` cardinality = `2 × len(ALT)`**
  (`b633a2b`). Per §1.6.1.7 each confidence-interval
  field must contain exactly `2 * len(ALT)` entries
  when present.

- **`MEINFO / METRANS` cardinality = `4 × len(ALT)`**
  (`dc22cff`). Each mobile-element field must contain
  exactly `4 * len(ALT)` entries.

- **`CICN` cardinality = `2 × len(ALT)`** (`76737f9`).

- **`<CNV> / <DEL> / <DUP>` alleles share `SVLEN` when
  `FORMAT CN` is declared** (`df189f1`). Per the
  FORMAT-keys section of §1.6.1.7 the SVLEN of every
  CNV-family ALT must match.

- **`<CNV:TR>` requires `RUS` or `RUL`** (`0c2a109`).
  Per the tandem-repeats section of §1.6.1.7 a record
  carrying `<CNV:TR>` without `RUS` and without `RUL`
  is invalid.

- **`FT` per-sample no whitespace or `;`** (`87d350f`).
  Per §1.6.2 each FT token must be free of whitespace
  and semicolons.

- **`RUS / RUL / RUC / RB / CIRUC / CIRB` length =
  sum of `RN`; `RUS` length matches `RUL`** (`924d61a`).
  Both rules from §1.6.1.7.

### Added — cross-record checks (streaming Validator)

- **CHROM contiguous block** (`f1e9858`). Per §1.6.1
  every CHROM block must be contiguous; the streaming
  Validator tracks every CHROM seen and fires when a
  reorder is detected.

- **ID uniqueness across records** (`d90c7da`). Per
  §1.6.1 no identifier should appear in more than
  one data record.

- **Breakend mate PSL / PSO consistency** (`d1cdb32`).
  Per §1.6.2 PSO the PSL and PSO of two records forming
  a breakend mate pair must match. The check keys by
  both the record's own IDs and its MATEID so either
  mate may arrive first.

### Added — feature

- **`vcf.Open(path)`** (`c0427d3`). A new public
  helper that opens a path and returns a streaming
  `*Decoder`. The first two bytes are sniffed for the
  gzip magic; on match the reader is wrapped in
  `internal/bgzf.NewReader`, otherwise the raw stream
  is read directly. Both plain `.vcf` and compressed
  `.vcf.gz` use the same call site.

## [0.8.0] — 2026-10-18

This release closes out the v4.5 spec-compliance work begun in
v0.7.0: 15 commits since v0.7.0, all in the `vcf` package. The
checks are additions only — no breaking changes to the public
API. Existing files that don't trigger the new checks continue
to validate as before.

### Added — VCF v4.5 spec compliance

- **POS / REF must fit within `##contig` length**
  (`e0709ea`). POS and POS+len(REF)-1 must not exceed the
  declared contig length (per §1.6.1.2). Records on an
  undeclared contig or one without a length are unaffected.

- **Symbolic ALT `<ID>` must be declared in `##ALT`**
  (`eee560f`). Records using `<DEL>`, `<DUP>`, etc. must have a
  matching `##ALT=<ID=...>` line. Breakend forms and `<*>`
  are accepted as the spec defines them rather than via `##ALT`.

- **Declared END must match the computed value**
  (`2bdfd3c`). When a record carries END, the value must equal
  the computed END (POS+len(REF)-1, plus SVLEN/LEN adjustments
  for symbolic ALT alleles).

- **ID column rejects `,` and `=`** (`ba45d9e`). The parser
  already rejected whitespace and semicolon in column 3;
  comma and equals are added so IDs round-trip cleanly through
  BCF v2.2.

- **Reserved FORMAT keys PL / GL / GP** (`9933e54`). All three
  use `Number=G`. PL must be Integer; GL and GP must be Float.
  `IsReservedFormatKey` now returns true for all three.

- **GT must be first FORMAT field** (`1913e5a`). Per §1.6.2,
  GT must be the first key in the FORMAT column when declared.
  Records without GT are not flagged.

- **PL (Number=G) and AD (Number=R) cardinality regression
  tests** (`e65561c`). The existing `Number.Cardinality` path
  already returns the right count for both: PL is
  `genotypeCount(alleleCount, ploidy)`, AD is `altCount+1`.
  Locks in behaviour with tests.

- **Ploidy consistency across samples in a record**
  (`cf102ca`). Per §1.6.2, all samples in a record must share
  the same ploidy. Mixed ploidies make Number=G / Number=P
  cardinality ambiguous.

- **Mitochondrial chromosome exemptions** (`7e072fb`). Per
  §1.6.1.6, mtDNA uses a circular coordinate system. POS can
  legitimately exceed the contig length and END can wrap past
  the origin. The POS-length and END-match checks skip
  mtDNA. Recognised names: MT, Mt, mt, chrM, chrm, M, chrMT,
  chrMt, MTDNA, mtdna.

- **`##contig` species / taxonomy sub-fields** (`02f1bdd`).
  Previously these landed in the catch-all `c.Other`; now they
  are parsed into typed fields on the `Contig` struct
  (`Species`, `Taxonomy`) and round-trip through the encoder.

## [0.7.0] — 2026-10-18

This release rounds out the VCF v4.5 spec compliance surface
with 14 per-record and cross-record checks, the streaming
`Validator` API, full BCF 2.2 read/write, semantic CIGAR
parsing, `.fai` index generation, and the FASTA section in
GFF3 encoders. 28 commits since v0.6.0.

### Added — VCF spec compliance

- **Streaming `Validator` API** (`Validator.HeaderChecks`,
  `Next`, `Finish`, `Reset`). The previous batch `Validate`
  still works; the streaming API lets callers pipe records
  through the validator without buffering them. The state
  (cross-record `prevChrom`/`prevPos`, PS phase-style map)
  lives on the `Validator` itself. `fbc068a`.

- **GT structural checks** (VCFv4.5 §1.6.1.7). Mixed
  `|` and `/` separators inside a single GT fire a
  `SeverityError`; allele indices outside `[0, len(ALT)]`
  fire an error. `e0c9730`.

- **END is a computed field** per VCF v4.5 §1.6.1.6. The
  parser computes `END = POS + len(REF) - 1` when END is
  not declared; an explicitly declared END that disagrees
  with the computation fires an error in `ModeStrict`.
  `244e456`.

- **PSL / PSO / PSQ recognition** with `Number=P` cardinality
  (1 + per-sample ploidy). Reserved FORMAT keys declared in
  the header now match the spec's `Number=P` rule. `cb01bf4`.

- **M-family value range** (`MFamily`). `M` must be in
  `[0, 1]`; `DPM` and `ADM` must be non-negative integers.
  `2252f20`.

- **Non-printable byte rejection** in any header or data
  field. Bytes in `U+0000..U+0008`, `U+000B..U+000C`, and
  `U+000E..U+001F` are rejected at parse time. `1947130`.

- **Integer range check** for BCF-unsafe values: integers
  outside `[-2^31, -2^31 + 7]` fire a warning
  (the spec says these cannot round-trip through BCF).
  `1947130`.

- **SVLEN sign rule** per §1.6.1.7: positive for INS, negative
  for DEL, either sign for INV/DUP/CNV. `1947130`.

- **CIPOS / CIEND / CILEN must span 0** (the spec's
  confidence interval must include the breakpoint).
  `1947130`.

- **##contig name regex** (SAM-compatible, §1.4.7). Contig
  IDs that fail the regex are rejected at parse time.
  `7c1b302`.

- **Header-line trailing-tab rejection.** A `##` line whose
  final character is a tab (other than inside a quoted
  string) is rejected. `7c1b302`.

- **Strict Float format regex.** Rejects `1.`, `1e`, etc.;
  accepts `INF`, `-INF`, `NAN`, `NaN` case-insensitively.
  `65301fb`.

- **Sample-name whitespace check.** Sample names must not
  contain whitespace. `65301fb`.

- **SVCLAIM per-allele rules** per §1.6.1.7. Each symbolic
  ALLELE's SVCLAIM value must be in the allowed set for that
  allele. `65301fb`.

- **EVENTTYPE reserved-value check** (`IsReservedEventType`).
  The 15 reserved values from §1.6.1.7 (DEL, DEL:ME, INS,
  INS:ME, DUP, DUP:TANDEM, DUP:DISPERSED, INV, TRA,
  TRA:BALANCED, TRA:UNBALANCED, CHROMOTHRIPSIS,
  CHROMOPLEXY, BFB, DOUBLEMINUTE) are accepted; unreserved
  values fire a `SeverityWarning` (default mode) or
  `SeverityError` (`ModeStrict`). `eb4eed3`.

- **LAA value-range check.** Each LAA value must be a 1-based
  index in `[1, len(ALT)]`. `0658c19`.

- **LAA distinctness check.** LAA values must be distinct
  integers per the spec's "n distinct integers" rule.
  Duplicate values fire a `SeverityError`. `162c941`.

- **Cross-record PS phasing-style consistency.** Records
  that share a `PS` value (and the same sample + chromosome)
  must use the same GT separator style (phased `|` vs
  unphased `/`). The check is restricted to records that
  declare an explicit PS — the spec is silent on the implicit
  default phase set, and the spec example itself mixes
  separators on records without PS. State held in
  `Validator.phaseStyle`, cleared on chromosome change.
  `0a93ec3`.

- **##contig MD5 format validation.** The MD5 must be
  32 lowercase hex characters per §1.4.7. The actual digest
  is not recomputed (the sequence is not in the header);
  malformed MD5 strings produce a `ParseError` at parse
  time. `5996115`.

- **IUPAC ambiguity code reduction in REF.** REF now reduces
  each IUPAC code (R, Y, S, W, K, M, B, D, H, V) to one of
  its constituent bases during parsing (e.g. R → A, Y → C,
  K → G). Case is preserved. This satisfies the spec's
  rule that REF must be only A/C/G/T/N. `cbfce74`.

### Added — Other packages

- **GFF3 `###` FASTA section** (`Encoder.EncodeFastaSection`).
  GFF3 allows an optional FASTA block after the `###`
  directive; the encoder previously stopped at the feature
  section. The new method writes `###` then streams the
  verbatim content of an `io.Reader` (typically a FASTA
  file). One-shot: a second call returns an error.

- **Solexa Solexa+64 encoding** (`fastq.PhredEncoding`,
  `Decoder.SetPhredEncoding`, `Record.PhredEncoding`). The
  Solexa+64 byte range is identical to Phred+64; the
  difference is the per-byte error-probability formula.
  Callers with known Solexa data set the encoding flag; the
  library records it on each `Record` and never applies the
  formula itself. Default remains Phred+33 / Phred+64.

- **v4.5 reserved-key table** (`vcf.IsReservedInfoKey`,
  `IsReservedFormatKey`, `IsBaseModKey`, and a new
  `Validator.checkReservedKeys` step in `ModeStrict`). The
  package now knows the expected `Number` / `Type` for the
  v4.5-reserved INFO keys (M* family, RN/RUS/RUL/RUC/RB,
  CIRUC/CIRB/RUB, LAA/LAD) and FORMAT keys (LAA/LAD/LADF/
  LADR/LEC/LGL/LGP/LPL/LPP), and warns when a header
  declaration disagrees with the spec. Structural value
  validation (e.g. RB's value being a repeat motif present
  in REF) remains the caller's responsibility.

- **BAI region queries** (`bam/bai.Reference.Query`,
  `Reference.FindOffset`). Given a half-open range `[beg,
  end)`, `Query` returns the deduplicated BGZF `(Lo, Hi)`
  virtual-offset pairs from the candidate bins (looked up
  via `Reg2Bins`) plus the linear-index minimum offset per
  htslib's `bam_intv_to_offset` semantics. `FindOffset` is
  the position-only shortcut.

- **BCF per-sample FORMAT** (v2 scope). `bcf.Record` gains
  `Format []uint16` and `Samples []bcf.Sample`. The decoder
  reads the FORMAT list and per-sample typed values
  (cardinality driven by the FORMAT key's declared Number);
  the encoder mirrors. Records with `Format` nil / empty
  are byte-identical to v1 output (n_fmt = 0).

- **BCF 2.2 read/write**. The decoder already accepted
  the 2.2 magic; the encoder now does too. New API:
  `Encoder.SetBCFVersion(v uint8)` (1 for 2.1, 2 for 2.2)
  must be called before `WriteHeader`. The decoder exposes
  `Decoder.Version() uint8` so callers can detect 2.2 input.

  BCF 2.2 introduces a long-string type tag (`0xF` +
  `uint32` length) for strings of length 15 or more. The
  encoder emits `typeLong` when the configured version is
  2 and any string in the record meets the length threshold
  (ID, REF, ALT, or per-sample FORMAT string). The decoder
  accepts either `typeString` (7) or `typeLong` (0xF) for
  these fields.

- **`.fai` index generation** (`fai.Build`, `fai.Writer`).
  Build reads a FASTA document and returns an `*Index` in source
  order; `Writer.Write` serialises an `*Index` in canonical
  `.fai` format. The two together let callers produce `.fai`
  files without a separate atomic-rename dance. Build handles
  LF and CR-LF line endings and captures the byte offset of
  each record's first base along with `lineBases` / `lineWidth`
  from the first sequence line of the record.

- **Semantic CIGAR parsing** (`sam.ParseCigar`,
  `sam.FormatCigar`, `sam.CigarUnit`, `sam.CigarOp`). The
  library previously stored CIGAR strings verbatim and left
  interpretation to callers; the new parser decodes a CIGAR
  string into `[]CigarUnit{ Len, Op }`, where `Op` is one of
  nine named constants (`CigarM`, `CigarI`, `CigarD`, `CigarN`,
  `CigarS`, `CigarH`, `CigarP`, `CigarEq`, `CigarX`). Per-op
  helpers expose the SAM semantics: `ConsumesReference`,
  `ConsumesQuery`, `ReferenceLength`, `QueryLength`,
  `AlignmentEnd`. The numeric op-code matches the BAM packed
  format (M=0 … X=8) so a BAM CIGAR word can be inspected by
  isolating `op>>28`. Lower-case op chars are accepted because
  aligners sometimes emit them in extended-CIGAR contexts.

## [0.6.0] — 2026-10-18

This release adds the four text↔binary translation subcommands
to the CLI. Round-trip through every format combination tested
in `cmd/htsfile/translate_test.go` is clean. The BCF↔VCF side
drops per-sample FORMAT data (matching the v1 BCF reader scope).

### Added

- **CLI format-translation subcommands.** Four new
  subcommands — `sam2bam`, `bam2sam`, `vcf2bcf`, `bcf2vcf` —
  convert between text and binary forms of the alignment and
  variant formats. Each takes a single file argument whose
  extension is checked and writes the converted output to
  stdout. Wired through the existing CLI dispatcher in
  `cmd/htsfile/main.go`.

  ```sh
  htsfile sam2bam reads.sam > reads.bam
  htsfile bam2sam aligned.bam > aligned.sam
  htsfile vcf2bcf variants.vcf > variants.bcf
  htsfile bcf2vcf variants.bcf > variants.vcf
  ```

- **Translation helpers** in `cmd/htsfile/translate.go`:
  per-format header and record converters (e.g.
  `samToBAMHeader`, `bamToSAMRecord`, `vcfToBCFHeader`,
  `bcfToVCFRecord`), packed-CIGAR pack/unpack, bin
  computation. The BCF↔VCF translator drops per-sample FORMAT
  data (matching the v1 BCF reader's scope); VCF "." IDs
  are normalised to the BCF typeMissing form.

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
[0.7.0]: https://github.com/inful/htsfile/compare/v0.6.0...v0.7.0
[0.8.0]: https://github.com/inful/htsfile/compare/v0.7.0...v0.8.0
[0.9.0]: https://github.com/inful/htsfile/compare/v0.8.0...v0.9.0
[0.9.1]: https://github.com/inful/htsfile/compare/v0.9.0...v0.9.1
[Unreleased]: https://github.com/inful/htsfile/compare/v0.9.1...HEAD

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