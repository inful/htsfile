# Changelog

All notable changes to `htsfile` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **Internal refactor — value-type coercion.** `InfoValue` and `FieldValue`
  now embed an unexported `typedCore` so the four near-identical
  `FieldType` switches (decode/encode × INFO/FORMAT) collapse to two:
  `decodeTyped` in `parse.go` and `writeTyped` in `write.go`. Type-coercion
  rules are now defined in one place; a future change (e.g. how `NaN` is
  formatted, how percent-encoding round-trips a new edge case) is fixed
  once instead of four times.

  **API impact:** because Go does not allow promoted fields in struct
  literals, callers that construct `InfoValue` / `FieldValue` directly via
  `{Type: X, Strs: Y}` must switch to `{typedCore: typedCore{...}}`. As
  `typedCore` is unexported, external callers use field-by-field
  assignment:

  ```go
  v := vcf.InfoValue{}
  v.Type  = vcf.TypeInteger
  v.Ints  = []int64{42}
  ```

  Reads (`v.Ints`, `v.HasValue()`, `v.Int()`, `v.Float()`, `v.String()`)
  are unchanged — field promotion handles them transparently. The
  godoc example in `ExampleNewEncoder` was updated to reflect the new
  construction style. Behaviour is byte-for-byte identical: the
  `TestEncoder_RoundTripSpecExample` / `TestEncoder_RoundTripSVExample`
  corpus tests and the CLI's `roundtrip` subcommand both pass unchanged.

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