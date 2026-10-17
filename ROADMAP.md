# Roadmap

`htsfile` is a Go library for parsing, writing, and validating
bioinformatics file formats. As of v0.0.x it ships two formats
([VCF v4.5](https://samtools.github.io/hts-specs/VCFv4.5.pdf) and
[BED v1](https://samtools.github.io/hts-specs/BEDv1.pdf)). This
document tracks the planned expansion to additional formats from
the [samtools/hts-specs](https://samtools.github.io/hts-specs/) family
and adjacent bioinformatics formats.

## Design principles

These are the constraints that guide every expansion. They are
deliberately conservative — we have learned (see "lessons from the
applyStructured refactor" in this repo's commit history) that
premature abstraction is a net cost.

- **Per-package self-contained.** Each format is its own Go package
  (`vcf/`, `bed/`, future `sam/`, `gff/`, `fasta/`, …). No cross-package
  imports of internal helpers.
- **No shared infrastructure until ≥3 packages need it.** When a
  second package reaches for a helper that the first one has
  already written, the helper is duplicated. When a third package
  does the same, we extract. Two is not a pattern.
- **CLI auto-dispatch by extension.** The existing
  `cmd/htsfile` command already routes by file extension and
  sniffs the first line for stdin / unknown extensions. New
  formats plug in here with no new top-level subcommands.
- **Streaming, round-trippable, validating.** Every new format
  package exposes a streaming `Decoder`, a round-trippable
  `Encoder`, and a `Validator` with `ModeDefault` and
  `ModeStrict`. This is the contract the vcf/ and bed/ packages
  already implement.
- **Zero runtime dependencies.** Standard library only. This is a
  project-level rule, not a per-format rule.
- **TDD-first.** Tests for the public surface are written before
  the implementation. Coverage targets stay at the existing level
  (85%+ on the vcf/ and bed/ packages).

## Scope decisions

These are the explicit decisions made about what is and isn't in
this project's scope. Each is reversible, but reversing one is a
deliberate decision.

- **In scope:** SAM, BAM, GFF3, FASTA, FASTQ, VCF v4.1–v4.5 (currently
  v4.5-only), BCF v2.1, Tabix, CSI, BAI (read-side). See tier list
  below.
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

## Tiers

Each tier is a logical commit set, not a hard release boundary.
Tiers are sequenced so that work in tier N is unblocked by tier
N−1 and unblocks tier N+1. Where useful, a tier is broken into
multiple phases (types → decoder → encoder → validator → CLI → docs)
in the style established for the BED v1 work.

### Tier 1 — SAM v1 (next)

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

### Tier 2 — Companion and adjacent formats

- [ ] **BAM** — binary SAM with BGZF compression. Reader first;
      BAI index reader for random access. New CLI subcommand or
      `--binary` flag.
- [ ] **GFF3** — feature annotation. 1-based inclusive coordinates
      (different convention from BED's 0-based half-open). 9
      mandatory columns + free-form `attributes` column. Strand and
      phase are typed values.
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
  CHANGELOG entry under `[Unreleased]`. Mirror the BED v1 commit
  (`f5a5530`) as the template.

## Open decisions to resolve before Tier 1 starts

These are not blockers, but answering them up front saves churn
mid-implementation.

- [ ] **Header parsing for SAM** — SAM headers have a free-form
      key=value syntax on `@HD` lines and a structured format on
      `@SQ`. Decide whether to model the structured fields
      explicitly (`SQ.Length`, `SQ.MD5`, `SQ.Species`, …) or
      store the raw `Other` map for unknown fields.
- [ ] **BAM CLI surface** — separate `htsfile sam view` and
      `htsfile bam view` subcommands, or a single `htsfile view`
      with format detection? Recommendation: keep separate
      subcommands if the read/write semantics diverge
      meaningfully; otherwise rely on extension dispatch.
- [ ] **FASTA index scope** — do we read `.fai` files, write
      them, or both? Recommendation: read-only for v1, since
      index generation is a separate effort.
- [ ] **VCF v4.1–v4.5 strictness** — does "support v4.1" mean
      "accept the file and produce identical output" or
      "validate against the v4.1 rules"? Recommendation: accept
      + emit, but validate per the strictest applicable version
      unless the user passes a flag.

## What we have today

| Package | Status | Coverage | Test surface |
|---------|--------|----------|--------------|
| `vcf/`   | v4.5   | 85.2%    | Decoder, Encoder, Validator, errors, types, corpus round-trip |
| `bed/`   | v1     | 100% on data; corpus round-trip | Decoder, Encoder, Validator, errors, types |
| `cmd/htsfile/` | n/a | n/a    | CLI: `check`, `roundtrip`, `summary`, `version`; auto-dispatch by extension |

Both packages follow the same shape (one `Decoder`, one
`Encoder`, one `Validator`, one `Record` type, one `Errors` file,
one `*_test.go` per source file plus a corpus round-trip test and
a public godoc example). New packages should land in the same
shape.

## Anti-patterns we have learned to avoid

These are codified in commit history and in the Brain preferences
(`pref-mid-task-reassessment`, `pref-similar-shape-not-duplication`).
Future expansion work should respect them.

- **Don't refactor for the sake of it.** Similar shape is not
  duplicated content. Before starting a structural refactor,
  count the identical lines per site vs the unique lines per
  site. Only refactor when the shared body is substantial (>50%
  identical). See commit history for the `applyStructured`
  revert.
- **Don't extract shared abstractions until ≥3 packages need
  them.** Two is not a pattern.
- **Don't bundle orthogonal changes in a single commit.** Each
  tier lands as one atomic commit, but a tier is not "everything
  I touched this month" — it's the smallest meaningful unit.
- **Don't introduce runtime dependencies.** The project is
  stdlib-only by design.
- **Don't ship without tests.** New packages follow TDD: tests
  for the public surface are written first.
