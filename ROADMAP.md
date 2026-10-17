# Roadmap

`htsfile` is a Go library for parsing, writing, and validating
bioinformatics file formats. As of v0.0.x it ships three formats
([VCF v4.5](https://samtools.github.io/hts-specs/VCFv4.5.pdf),
[BED v1](https://samtools.github.io/hts-specs/BEDv1.pdf), and
[SAM v1](https://samtools.github.io/hts-specs/SAMv1.pdf)). This
document tracks the planned expansion to additional formats from
the [samtools/hts-specs](https://samtools.github.io/hts-specs/) family
and adjacent bioinformatics formats.

> **Where to continue from this point in a new session:** start with
> **Tier 2 (GFF3)**. BAM Tier 2 just landed (commit replacing this
> note). The roadmap is up to date as of the latest commit on
> `main`. See the [CHANGELOG](CHANGELOG.md) for what just
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
- [ ] **GFF3** — feature annotation. 1-based inclusive coordinates
      (different convention from BED's 0-based half-open). 9
      mandatory columns + free-form `attributes` column. Strand and
      phase are typed values. *Next up.*
- [ ] **FASTA** — sequence foundation. Two record types (header
      starting with `>` + multi-line sequence). Optional `.fai`
      index file for random access.

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
- [ ] **BCF v2.1** — binary VCF. Companion to VCF; useful for
      round-tripping large VCFs without the text overhead.
- [ ] **Tabix** — generic index for sorted, bgzipped position
      files. Most useful once BAM lands.
- [ ] **CSI v1** — successor to BAI. Same problem class as
      Tabix; may share infrastructure.

If/when a shared `index/` package is needed (Tabix + CSI + BAI in
three places), it gets extracted then. Not before.

### Tier 4 — Research / optional

- [ ] **CRAM v3** — depends on Tier 2 (BAM) being stable.
- [ ] **htsget / refget** — separate module if at all.
- [ ] **crypt4gh** — encryption; cross-cutting, may live in
      `internal/crypt/`.
- [ ] **FASTQ** — multiple variants; needs design discussion
      before implementation. Will start with Sanger Phred+33
      (the dominant modern convention) and Illumina 1.8+ (also
      Phred+33 but different header conventions).

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

- **CLI duplication.** `runCheck{VCF,BED,SAM,BAM}`,
  `runRoundtrip{VCF,BED,SAM}`, and `runSummary{VCF,BED,SAM,BAM}`
  are intentionally near-duplicates and carry `//nolint:dupl`
  with rationale. The duplication is structural across
  per-format entry points. BAM landed before the threshold
  refactor was due (the ROADMAP originally pinned GFF3 as the
  trigger), so the BAM entries joined the existing pattern;
  the four-way refactor (extract an issue-reporting helper)
  is queued for the next format landing.

## What we have today

| Package        | Status | Coverage | Test surface |
|----------------|--------|----------|--------------|
| `vcf/`         | v4.5   | 85.3%    | Decoder (with `NewDecoderWithMaxLine`), Encoder, Validator, errors, types, corpus round-trip |
| `bed/`         | v1     | 90.5%    | Decoder (with `NewDecoderWithMaxLine`), Encoder, Validator, errors, types, corpus round-trip |
| `sam/`         | v1     | 83.6%    | Decoder (with `NewDecoderWithMaxLine`), Encoder, Validator, errors, types, corpus round-trip |
| `bam/`         | v1     | 76.2%    | Decoder (with `NewDecoderWithMaxLine`), Validator, errors, types, B-array round-trip, CLI smoke |
| `bam/bai/`     | n/a    | 84.1%    | Index parser + `Reg2Bins` helper |
| `internal/bgzf/` | n/a  | 77.5%    | Block decoder (BSIZE, CRC32, ISIZE, EOF marker) |
| `cmd/htsfile/` | n/a    | n/a      | CLI: `check`, `roundtrip`, `summary`, `version`; auto-dispatch by `.vcf`/`.bed`/`.sam`/`.bam` extension + content sniff for stdin |

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

### For Tier 2 — GFF3 and FASTA

- [ ] **GFF3 attributes column.** Free-form `key=value; key2=value2`
      pairs. The v1 parser will store the raw text plus
      well-known keys (ID, Name, Parent, Alias, Target, Gap, etc.)
      in a map. Open: do we also expose typed accessors (e.g.
      `rec.Attributes["Parent"]` returns a `[]string` for
      multi-value keys)?

- [ ] **FASTA index scope.** Read-only, write, or both?
      *Recommendation:* read-only for v1 (consistent with
      ROADMAP's "Index generation" out-of-scope note). If
      streaming-only access turns out to be enough, this can
      stay read-only forever.

- [ ] **FASTQ quality encoding.** Multiple incompatible
      encodings (Phred+33 Sanger/Illumina 1.8+, Phred+64 Solexa,
      Phred+64 Illumina 1.3–1.7). v1 should at minimum auto-detect
      and parse Phred+33 (the dominant modern convention); an
      `Encoding` field on the Record would let callers override.

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
