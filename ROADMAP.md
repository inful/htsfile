# Roadmap

`htsfile` is a Go library for parsing, writing, and validating
bioinformatics file formats. As of v0.0.x it ships three formats
([VCF v4.5](https://samtools.github.io/hts-specs/VCFv4.5.pdf),
[BED v1](https://samtools.github.io/hts-specs/BEDv1.pdf), and
[SAM v1](https://samtools.github.io/hts-specs/SAMv1.pdf)). This
document tracks the planned expansion to additional formats from
the [samtools/hts-specs](https://samtools.github.io/hts-specs/) family
and adjacent bioinformatics formats.

> **Where to continue from this point in a new session:**
> Tiers 1–3 are complete. Tier 4 status: FASTQ shipped in
> v0.3.0; Tabix + CSI writers + BCF record encoder shipped
> in v0.4.0 + v0.5.0 (every format the library reads now has a
> Writer; CLI roundtrip works for all of them). CRAM v3 is
> still at the design-survey stage (see `cram/SPEC.md`); the
> codec gap blocks implementation. htsget/refget and
> crypt4gh remain pending (htsget/refget's architectural
> shape — separate directory vs separate Go module — is
> explicitly deferred per the ROADMAP §Open decisions). The roadmap is up to date as of the latest commit
> on `main`. See the [CHANGELOG](CHANGELOG.md) for what just
> landed and `git log --oneline` for the commit history.

## Design principles

These are the constraints that guide every expansion. They are
deliberately conservative — we have learned (see commit `a9fae1f`
reverting the `applyStructured` refactor) that premature abstraction
is a net cost.

- **Per-package self-contained.** Each format is its own Go package
  (`vcf/`, `bed/`, `sam/`, future `gff/`, `fasta/`, `bam/`, …). No
  cross-package imports of internal helpers.
- **No shared infrastructure until ≥3 packages need it.** When a
  second package reaches for a helper that the first one has
  already written, the helper is duplicated. When a third package
  does the same, we extract. Two is not a pattern.
- **CLI auto-dispatch by extension.** The existing
  `cmd/htsfile` command routes by file extension (`.vcf`, `.bed`,
  `.sam`) and sniffs the first line for stdin / unknown
  extensions. New formats plug in here with no new top-level
  subcommands.
- **Streaming, round-trippable, validating.** Every new format
  package exposes a streaming `Decoder`, a round-trippable
  `Encoder`, and a `Validator` with `ModeDefault` and
  `ModeStrict`. This is the contract the existing packages
  already implement.
- **Zero runtime dependencies.** Standard library only. This is a
  project-level rule, not a per-format rule.
- **TDD-first.** Tests for the public surface are written before
  the implementation. Coverage targets stay at the existing level
  (85%+ across all packages).
- **Bounded-line reads.** Every Decoder caps per-line reads at
  `DefaultMaxLineLength` (10 MB) via `bufio.Scanner` to prevent
  DoS via hostile input. A new format inherits this pattern.

## Scope decisions

These are the explicit decisions made about what is and isn't in
this project's scope. Each is reversible, but reversing one is a
deliberate decision.

- **In scope:** BAM, GFF3, FASTA, FASTQ, BCF v2.1, Tabix, CSI, BAI
  (read-side). See tier list below.
- **Out of scope:**
  - **CRAM v3** — sophisticated codec system; deferred until BAM
    and the BGZF/BAI infrastructure are stable.
  - **htsget / refget** — network protocols with a different shape
    (server + client). If added, likely as a separate module
    (`htsnet/`) rather than a format package.
  - **crypt4gh** — encryption wrapper; cross-cutting concern that
    doesn't fit the per-format package model.
  - **bigBed / bigWig** — binary indexed formats; their own
    ecosystem. Not in this roadmap.
  - **Index generation** — we read BAI/CSI/Tabix but do not write
    them. Writing is a different problem (multi-pass, atomic
    rename) and a future effort if needed.
  - **FASTA/FASTQ multi-variant handling beyond Sanger / Illumina
    conventions** — FASTQ has several incompatible variants
    (Illumina 1.3–1.7 binned quality, Sanger Phred+33, Solexa
    Solexa+64, …). We will support the two most common encodings
    explicitly and document the rest as future work.
  - **Per-version VCF validation** (v4.1–v4.5) — the library
    accepts and round-trips any v4.x `##fileformat`, but does
    **not** implement per-version validation rule differences
    (which reserved keys are recognised, which `Number=`
    values are valid). Add only if a user reports a concrete
    per-version bug.

## Tiers

Each tier is a logical commit set, not a hard release boundary.
Tiers are sequenced so that work in tier N is unblocked by tier
N−1 and unblocks tier N+1. Where useful, a tier is broken into
multiple phases (types → decoder → encoder → validator → CLI → docs)
in the style established for the BED v1 work.

### Tier 1 — SAM v1 — **DONE** (commit `3484a4f`)

- [x] `sam/` package: `Record` (11 mandatory columns + variable
      tag map), `Header` (HD / SQ / RG / PG / CO lines + optional
      fields), `Decoder`, `Encoder`, `Validator`
- [x] CLI integration via extension (`.sam` → SAM) and content
      sniff (`@HD` header line → SAM; 11-column tab-separated
      with integer FLAG/MAPQ for headerless files)
- [x] SAMtags companion: predefined optional tag types (`A`, `C`,
      `f`, `H`, `i`, `Z`, `B`) — kept as raw `Value string` in v1;
      typed accessors deferred until needed
- [x] Test corpus: small synthetic SAM file (5 records covering
      BED-style paired/stranded/secondary cases)
- [x] Public godoc example
- [x] README + CHANGELOG
- [x] Post-merge fixes: malformed header type tags error
      (`@H\trandom` no longer silently stored in `COs`);
      duplicate `@HD` lines error; dead code removed
      (commit `60a300a`)
- [x] DoS hardening: bounded per-line read at
      `DefaultMaxLineLength` (10 MB) across all three packages
      (commit `cd6cb08`)

### Tier 2 — Companion and adjacent formats

- [x] **BAM** — binary SAM with BGZF compression. Reader first;
      BAI index reader for random access. *Completed.* Includes:
      - `internal/bgzf/` — zero-dependency BGZF block decoder
        (BSIZE, CRC32, ISIZE validation; EOF marker)
      - `bam/` — streaming `Decoder`, per-record `Validator`
        with `ModeDefault` and `ModeStrict`, all 12 BAM tag
        type codes, full aux-tag parsing including B-arrays,
        4-bit packed SEQ decoding, CIGAR exposed as
        `[]uint32` and as `CigarString()`. Read-only — no BAM
        writer yet.
      - `bam/bai/` — BAI file parser (bin index + linear
        intervals + `n_no_coor`); `Reg2Bins` helper for the
        samtools 6-level bin scheme; region queries deferred.
      - CLI: `check` and `summary` for `.bam`; `roundtrip`
        returns "read-only" error.
- [x] **GFF3** — feature annotation. 1-based inclusive coordinates
      (different convention from BED's 0-based half-open). 9
      mandatory columns + free-form `attributes` column. Strand and
      phase are typed values. *Completed.* Includes:
      - `gff/` — streaming `Decoder`, round-trippable `Encoder`,
        per-record `Validator` (`ModeDefault` + `ModeStrict`
        with Parent ID cross-reference check); `##gff-version`
        + `##sequence-region` + `##date` parsed; arbitrary `##`
        pragmas preserved verbatim; `###` terminator ends the
        feature section. Multi-value attributes exposed as
        `map[string][]string` with `Attribute` /
        `AttributeList` accessors.
      - CLI: `check`, `roundtrip`, `summary` for `.gff` /
        `.gff3`; content sniff `##gff-version` or 9-column
        feature line.
      - FASTA emission after `###` is intentionally out of
        scope for v1.
- [x] **FASTA** — sequence foundation. Two record types (header
      starting with `>` + multi-line sequence). Optional `.fai`
      index file for random access. *Completed.* Includes:
      - `fasta/` — streaming `Decoder`, round-trippable `Encoder`
        with configurable sequence-line wrapping (default 80),
        per-record `Validator` (`ModeDefault` + `ModeStrict`
        alphabet check).
      - `fasta/fai/` — `.fai` index parser (5 tab-separated
        fields per reference: name, length, offset, lineBases,
        lineWidth). Read-only; generation out of scope.
      - CLI: `check`, `roundtrip`, `summary` for `.fa` /
        `.fasta` / `.fna` / `.ffn` / `.faa`; content sniff
        `>` header line.

  Also in this commit: the ≥3-packages CLI duplication
  refactor the ROADMAP flagged. The five per-format `runCheck*`
  duplicates have been replaced by a single shared
  `runIssueReport` helper plus per-package adapters.

### Tier 3 — VCF ecosystem and indexing

- [x] **VCF v4.1–v4.5 acceptance (any v4.x minor).** Parser
      reads `##fileformat=VCFvX.Y` and stores the version on the
      `Header`; the writer round-trips it verbatim. Major
      versions other than 4 are rejected. Per-version
      validation rule differences (which reserved keys are
      recognised, which `Number=` values are valid) are **not**
      yet implemented — the library treats all v4.x identically
      and asks users to report per-version bugs as they surface.
      (commit `2b48959`)
- [x] **BCF v2.1 (read-only sub-tier)** — binary VCF. Companion
      to VCF. **Reader-only in v1** per the open-decision scope
      ("per the ROADMAP's Index generation out-of-scope note").
      Includes:
      - `bcf/` — streaming `Decoder`, per-record `Validator`
        (`ModeDefault` + `ModeStrict`), all 12 BAM-equivalent
        INFO type codes, and a small but useful per-record
        API (`LookupContig`, `LookupFilter`, `InfosOrder`).
      - `ParseHeader` is exposed so test helpers can reuse the
        VCF-text parser.
      - `LookupInfoIndex` / `LookupFormatIndex` helpers fix the
        off-by-N bug class for future writers.
      - **No BCF writer** in v1 — round-trip would require a
        ~2× larger package (writer alone) and is incomplete
        without per-sample data. Deferred to v2.
      - **No per-sample FORMAT column** parsing in v1 — too
        complex for an initial sub-tier; surfaced as a clear
        ParseError rather than guessed at.
      - CLI: `check`, `summary` for `.bcf`; `roundtrip`
        returns "BCF v1 is read-only" error.
- [x] **Tabix (read-only sub-tier 3B)** — generic position
      index for sorted bgzipped text files. **Reader-only in
      v1.** Includes:
      - `tabix/` — `Reader` + `Index` + per-reference bin /
        linear index. Bin numbers follow the samtools 5-level
        scheme; helper `BinFromRange(beg, end)` returns the
        bins that overlap a region.
      - Magic `"TBI\x01"`, BGZF virtual offsets stored as
        `uint64` (high 48 = compressed-block offset; low 16 =
        offset within the decompressed block).
      - The optional `n_no_coor` field is always read; EOF is
        treated as the absence signal. Multi-reference legacy
        files (pre-samtools-0.1.18) without `n_no_coor` will
        misalign — use CSI for those.
      - CLI: `check`, `summary` for `.tbi`; `roundtrip`
        returns "Tabix is read-only" error.
- [x] **CSI v1 (read-only sub-tier 3B)** — Coordinate-Sorted
      Index. **Reader-only in v1.** Includes:
      - `csi/` — `Reader` + `Index` + per-reference bin /
        linear index. Tunable `min_shift` (window size exponent)
        and `depth` (bin levels); defaults match Tabix (14 / 5).
      - Magic `"CSI\x01"`, explicit version byte disambiguates
        format variants — no `n_no_coor` ambiguity.
      - Helper `Index.WindowSize()` returns `1 << MinShift` bp.
      - CLI: `check`, `summary` for `.csi`; `roundtrip`
        returns "CSI is read-only" error.

Both Tabix and CSI share the same 5-level binning math and
BGZF virtual-offset concept; the only meaningful divergence is
the format-detection story (CSI has a version byte, Tabix has
an optional trailer field). v1 inlines `reg2bins` in both
packages rather than extracting a shared `internal/genomic`
package — the bin math is small and the duplication is
intentional (a shared helper would be a third cross-package
import, which the per-package self-contained rule discourages).

If/when a shared `index/` package is needed for writers (Tabix
+ CSI + BAI in three places), it gets extracted then. Not before.

### Tier 4 — Research / optional

- [x] **FASTQ** — Sanger Phred+33 + Illumina 1.8+ (auto-detect
      Phred+64 via max-byte heuristic). Parser, writer,
      validator, CLI integration. Done in commit `f41a671`
      (v0.3.0).

- [ ] **CRAM v3** — design phase only at this point. CRAM v3.1
      requires codecs (LZMA, rANS 4x8, rANS Nx16, adaptive
      arithmetic, fqzcomp, name tokeniser) that have no Go
      stdlib equivalents. Implementing them from scratch is
      ~3500–5000 LOC of subtle codec code, and the project's
      "no runtime dependencies" rule precludes using
      samtools/htscodecs. Detailed scope survey in
      `cram/SPEC.md` (commit `3ac192b`). Resumable once one of:
      (a) stdlib adds the codecs, (b) the no-deps rule is
      relaxed for a separate `htscram/` module with its own
      go.mod, or (c) real-world CRAM testdata becomes available
      so codec correctness can be validated against known-good
      files.

- [ ] **htsget / refget** — separate module if at all.

- [ ] **crypt4gh** — encryption; cross-cutting, may live in
      `internal/crypt/`.

## Cross-cutting

These are concerns that affect all tiers. They are tracked here
rather than per-tier to avoid duplication.

- **Test corpus sourcing.** Public datasets:
  - SAM/BAM: 1000 Genomes, gnomAD
  - GFF3: Ensembl, GENCODE
  - FASTA/FASTQ: NCBI
  - VCF: 1000 Genomes (already a small corpus in `vcf/testdata/`)
  
  For each new format, build a curated `testdata/` with one or two
  small representative files plus a larger file for streaming
  performance tests.

- **Public API stability.** New packages start at v0.0.x. We
  follow semver and document API-impact notes in CHANGELOG under
  `[Unreleased]`.

- **CI matrix.** Existing matrix (Go 1.22/1.24/stable × linux /
  macos) stays. New packages inherit the same lint, vet, test,
  and coverage checks.

- **CHANGELOG discipline.** Every tier's worth of work lands as
  one atomic commit with a Conventional Commits message and a
  CHANGELOG entry under `[Unreleased]`. Mirror the SAM v1 commit
  (`3484a4f`) as the template.

- **CLI duplication.** The five per-format `runCheck*`
  functions used to be near-duplicates with `//nolint:dupl`.
  The ≥3-packages refactor flagged in this ROADMAP has now
  landed: a single shared `runIssueReport` helper plus
  per-package adapters replaces the duplicated sort/print/
  count loop. The `runRoundtrip*` and `runSummary*` families
  are still parallel (different stats per format make a
  shared abstraction less clear), but each new format keeps
  the same pattern and the per-format lines are short.

## What we have today

| Package        | Status | Coverage | Test surface |
|----------------|--------|----------|--------------|
| `vcf/`         | v4.5   | 85.3%    | Decoder (with `NewDecoderWithMaxLine`), Encoder, Validator, errors, types, corpus round-trip |
| `bed/`         | v1     | 90.5%    | Decoder (with `NewDecoderWithMaxLine`), Encoder, Validator, errors, types, corpus round-trip |
| `sam/`         | v1     | 83.6%    | Decoder (with `NewDecoderWithMaxLine`), Encoder, Validator, errors, types, corpus round-trip |
| `bam/`         | v1     | 76.2%    | Decoder (with `NewDecoderWithMaxLine`), Validator, errors, types, B-array round-trip, CLI smoke |
| `bam/bai/`     | n/a    | 84.1%    | Index parser + `Reg2Bins` helper |
| `internal/bgzf/` | n/a  | 77.5%    | Block decoder (BSIZE, CRC32, ISIZE, EOF marker) |
| `gff/`         | v3     | 75.7%    | Decoder (with `NewDecoderWithMaxLine`), Encoder, Validator, errors, types, multi-value attrs, CLI smoke |
| `fasta/`       | n/a    | 75.6%    | Decoder (with `NewDecoderWithMaxLine`), Encoder, Validator, errors, types, round-trip, CLI smoke |
| `fasta/fai/`   | n/a    | 80.6%    | Index parser (5 fields: name, length, offset, lineBases, lineWidth) |
| `bcf/`         | v2.1   | 56.4%    | Decoder (BGZF-wrapped), Validator (`ModeDefault`/`ModeStrict`), errors, types, typed INFO, CLI smoke; **read-only** — no writer |
| `tabix/`       | v1     | 87.4%    | Tabix index reader; bin + linear index + chunk lookup; `BinFromRange` helper; CLI smoke |
| `csi/`         | v1     | 87.2%    | CSI index reader (version-aware via `CSI\x01` magic); tunable `min_shift`/`depth`; CLI smoke |
| `cmd/htsfile/` | n/a    | n/a      | CLI: `check`, `roundtrip`, `summary`, `version`; auto-dispatch by `.vcf`/`.bed`/`.sam`/`.bam`/`.gff`/`.gff3`/`.fa`/`.fasta`/`.fna`/`.ffn`/`.faa`/`.bcf`/`.tbi`/`.csi` extension + content sniff for stdin; shared `runIssueReport` helper replaces the 5-way `runCheck*` duplication |

All three format packages follow the same shape (one `Decoder`, one
`Encoder`, one `Validator`, one `Record`/`Header` type, one
`Errors` file, one `*_test.go` per source file plus a corpus
round-trip test and a public godoc example). New packages should
land in the same shape and inherit the bounded-line-read pattern
(`DefaultMaxLineLength` / `NewDecoderWithMaxLine`).

## Open decisions

These are the unresolved questions for the next tier of work.
Each has a recommendation; the recommendation should be confirmed
before starting implementation.

### For Tier 4 — FASTQ / SAM/VCF converters / CRAM v3 / crypt4gh

FASTQ shipped in v0.3.0. CRAM v3 is at the design-survey
stage (see `cram/SPEC.md`); the codec gap is the blocker.
htsget/refget and crypt4gh remain pending.

- [x] **Index generation — Tabix + CSI writers.** Done in
      commits `22ac357` (tabix) and `996a347` (csi). Both
      packages now expose `Writer` + `Write(*Index)`. The CLI
      `roundtrip` subcommand works for `.tbi` and `.csi`
      files. Closes the Tabix + CSI halves of this item.
- [x] **Index generation — BCF writer (record-side).**
      The v0.5.0 release added the BCF record encoder; it
      writes records without per-sample FORMAT data (matching
      v1's BCF reader scope). The remaining asymmetry is the
      per-sample binary FORMAT encoding; deferred unless a
      real-world use case surfaces.
      follow-up. Tabix + CSI writers shipped (see the Tier 4
      item above). BCF writer still deferred for per-sample
      FORMAT data.
- [x] **CLI format translation — sam2bam / bam2sam / vcf2bcf /
      bcf2vcf.** Landed in the v0.6.0 section of the
      CHANGELOG as a CLI-only feature. Translators live in
      `cmd/htsfile/translate.go`; conversion is text ↔ binary
      with the per-sample FORMAT caveat above. Round-trips
      cleanly through every format combination tested.
- [ ] **htsget / refget** — separate package from the per-format
      work, scope TBD (see §Open decisions for the two
      architectural readings).
- [ ] **VCF v4.5 spec compliance — remaining checks.** The
      high-value batch landed in this release (POS ≤ contig
      length, symbolic ALT declared in `##ALT`, END matches
      computed value, ID column character set). Remaining:

  - **PL / GL / GP cardinality = `ploidy*(ploidy-1)/2 + ploidy`**
    (the "Number=G" formula). Needs per-sample ploidy from
    GT. Medium scope.
  - **AD cardinality = len(REF) + len(ALT)** for
    `Number=A`. Same plumbing as above.
  - **`Number=G` reserved keys must have `Type=Float` or
    `Type=Integer`** (PL, GL, GP). The reserved-key table
    knows `Number=G` but does not yet validate the Type
    against it.
  - **GT must be the first FORMAT field** if present
    (recommended, not required).
  - **Ploidy consistency across a record** — when GTs
    have different ploidies, `Number=G` and `Number=P`
    cardinality must use the maximum ploidy per record.
  - **`##contig` species / taxonomy / assembly** field
    validation.
  - **Mitochondrial DNA (`chrM`) special-case rules** from
    §1.6.1.6 (mtDNA uses a circular-aware coordinate system).
  - **Compressed `.vcf.gz` reader pipeline** (BGZF wrapper
    + tabix/CSI region access). All the building blocks
    exist (`internal/bgzf`, `tabix`, `csi`); no wired-up
    `vcf.Open(path)` that returns a streaming decoder.

  Tracked here as Tier 3 follow-up; not blockers for any
  released feature. The full list of checks shipped in
  v0.7.0 lives in the `v0.7.0` CHANGELOG entry.
- [ ] **CRAM v3** — design-survey stage only. Codecs 3–8 are
      not in Go stdlib and the no-runtime-deps rule precludes
      samtools/htscodecs. Resume when stdlib adds the codecs,
      a separate `htscram/` module is acceptable, or
      real-world testdata becomes available. See `cram/SPEC.md`.
- [ ] **htsget / refget** — separate package from the per-format
      work. See §Tier 4 above for the two architectural readings
      (separate directory vs separate Go module); decision
      deferred to first user pull for the protocols.
- [ ] **crypt4gh** — encryption; cross-cutting.

## Anti-patterns we have learned to avoid

These are codified in commit history and in the Brain preferences
(`pref-mid-task-reassessment`, `pref-similar-shape-not-duplication`).
Future expansion work should respect them.

- **Don't refactor for the sake of it.** Similar shape is not
  duplicated content. Before starting a structural refactor,
  count the identical lines per site vs the unique lines per
  site. Only refactor when the shared body is substantial (>50%
  identical). See commit history for the `applyStructured`
  revert (`a9fae1f` was a successful refactor; the earlier
  attempt was reverted before this session).
- **Don't extract shared abstractions until ≥3 packages need
  them.** Two is not a pattern.
- **Don't bundle orthogonal changes in a single commit.** Each
  tier lands as one atomic commit, but a tier is not "everything
  I touched this month" — it's the smallest meaningful unit.
- **Don't introduce runtime dependencies.** The project is
  stdlib-only by design.
- **Don't ship without tests.** New packages follow TDD: tests
  for the public surface are written first.
- **Don't hide errors in default fallbacks.** When a parser
  defaults to a format on ambiguous input (e.g. the content
  sniff defaulting to VCF), document the heuristic explicitly
  and surface the real parse error to the user, not a
  misleading one.
