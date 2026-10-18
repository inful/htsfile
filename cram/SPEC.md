# CRAM v3 design — scope survey

This document surveys what a CRAM v3 implementation would entail
in this project, given the project's "standard library only" rule
(see [ROADMAP.md §Design principles](../ROADMAP.md) and
[README.md §Spec-compliant](../README.md)).

It is intentionally *not* an implementation plan. The goal is to
quantify the trade-offs and surface the codec-scope decision before
committing to code.

## What CRAM v3 actually requires

CRAM v3.1 (and v3.0) is two things stacked:

1. **The CRAM format** — a container-and-block structure with
   reference-based record encoding. Container header, block
   headers (method, content-type, content-id, compressed size,
   uncompressed size, CRC32), slice headers, alignment-record
   encoding using CRAM's tag-typed integer encodings (itf8,
   ltf8, htf8), reference-based sequence and quality
   reconstruction.

2. **The CRAM block codecs** — eight compression methods that
   can appear on any block (raw, gzip, bzip2, LZMA, rANS 4x8,
   rANS Nx16, adaptive arithmetic, fqzcomp, name tokeniser).
   Codecs 3–8 all post-date the project's existing codec
   coverage (BGZF = raw deflate + a 2-byte sub-field).

### Codec inventory vs Go stdlib

| Codec ID | Name        | First in | In Go stdlib?  | Decoder effort (LOC) |
|---------:|-------------|----------|----------------|----------------------|
| 0        | raw         | v2.0     | n/a            | 0                    |
| 1        | gzip (raw deflate) | v2.0 | `compress/flate` | 0 (already used) |
| 2        | bzip2       | v2.0     | `compress/bzip2` (Go ≥1.21) | 0 |
| 3        | LZMA        | v3.0     | **no**         | ~800–1200 (range coder + LZ state) |
| 4        | rANS 4x8    | v3.0     | **no**         | ~600–900 (4-state ANS renormalisation) |
| 5        | rANS Nx16   | v3.1     | **no**         | ~1200–1800 (transforms + pack + RLE) |
| 6        | adaptive arith | v3.1   | **no**         | ~500–700 |
| 7        | fqzcomp     | v3.1     | **no**         | ~500–800 (quality-specific) |
| 8        | name tokeniser | v3.1   | **no**         | ~400–700 |

Precedent: `internal/bgzf/` is 260 lines and BGZF is essentially
"raw deflate with a 2-byte sub-field". LZMA and rANS are each
roughly **5× the complexity** of BGZF and would live in
`internal/cramcodecs/` (analogous to `internal/bgzf/`).

The reference implementation lives in samtools/htscodecs
(separate from htslib precisely because these codecs need their
own test surface and history). Implementing them from scratch in
this project, with no spec-comparable test corpus, is high-risk.

## Format-only effort (no codec gaps)

Even if every codec were in stdlib, the CRAM format itself is
substantial:

- File definition (26-byte header, magic `CRAM`, version, file-id)
- Container header (variable length, length-itf8-encoded fields)
- Block header (method + content-type + content-id + sizes + CRC32)
- Slice header (ref-id, record-counter, record-count, ...)
- Alignment-record encoding with the CRAM tag-typed integer
  encodings (itf8, ltf8, htf8) and reference-based sequence /
  quality reconstruction (the "external blocks" pipeline)
- CRC32 over compressed payloads
- EOF container (required)

Estimate: ~1800–2400 LOC of format code + ~1000 LOC of tests.

## Scope options

### A. CRAM v2.1 read-only (codecs 0/1/2 only) — minimal, shippable

- Format-level support for v2.1 + v3.x files.
- Codec support: raw, gzip, bzip2 (all in stdlib).
- v3 files using LZMA / rANS / etc. fail with a clear
  `ParseError{Field: "block.method", Msg: "codec 3 (LZMA) is not supported in this build"}`.
- Region queries deferred (like BAM).
- No CRAM writer — same open-decision scope as BAM/BCF/Tabix/CSI.

Effort: ~2800–3400 LOC total (format + tests + CLI integration +
codec-adapter plumbing). Single-session work is feasible but
tight; probably wants the same "spec → design → commit" cadence
as the other Tier 3/4 packages.

### B. CRAM v3 with stdlib codecs only (same as A)

This is the same scope as A; the v3 vs v2.1 difference is in the
file-header version byte and one or two slice-header fields. v3
files using LZMA / rANS still fail with the same clear error.

### C. CRAM v3 with full codec coverage (codecs 0–8)

Adds ~3500–5000 LOC of codec implementations in
`internal/cramcodecs/`. Each codec needs its own decoder test
suite — the only meaningful test surface is real-world CRAM files
produced by samtools/htslib, which we do not currently have in
`testdata/`. This is multi-session work and the highest-risk
scope: codec correctness is subtle, and a single off-by-one in
the rANS renormalisation step produces silent data corruption.

### D. CRAM design doc only (no implementation)

This document. Useful if the scope decision is "we want CRAM
eventually but not this milestone".

### E. Defer CRAM entirely

Mark CRAM as deferred (no implementation) and update the
ROADMAP to record why: codecs 3–8 are not in stdlib and the
project's "no runtime dependencies" rule precludes using
htscodecs. Revisit only when (a) stdlib adds them, or (b) the
project relaxes the no-deps rule for CRAM, or (c) the work
shifts to a separate module (e.g. `htscram/`) with its own
dependency policy.

## Recommendation

**A or B** — implement the format + stdlib codecs, defer the
LZMA / rANS codecs to a future "v2 follow-up". Rationale:

1. The format work is the bulk of the value regardless of which
   codecs we ship. Once the format is in, adding codec X is a
   self-contained `internal/cramcodecs/<x>.go` task.
2. Bzip2 + raw deflate cover the "fast" and "normal" CRAM
   profiles. Files using LZMA / rANS are the "archive" profile,
   typically 1–5% of files in the wild depending on the data
   type. A clear error for those is better than no CRAM support
   at all.
3. The 3500–5000 LOC of codec implementations would be the
   largest single addition to the project to date — bigger than
   the entire BAM+BAI+BGZF work combined. It would also be the
   only code in the project that has no spec-comparable test
   corpus (we don't ship real-world CRAM testdata).
4. The "no deps" rule is a *project* rule, not a per-package
   rule. If we ever wanted to do C, we could revisit by adding
   CRAM to a separate `htscram/` module with its own go.mod and
   a `htscodecs` dependency. That's a different milestone.

The decision between A and B is mostly cosmetic — both produce
the same code. **B** is the right label because the spec talks
about v3, and we should label the package `cram/v3` (or just
`cram`) and document the codec gap explicitly in the package
doc-comment.

## Open questions (need user input before implementation)

1. **Scope: A/B, C, D, or E?** (see above)
2. **Test corpus**: do we have access to small real CRAM files
   for testdata/? Without them, the test suite relies entirely
   on synthetic byte streams (hard to exercise reference-based
   reconstruction realistically). For scope A/B, the codec
   subset we support (0/1/2) is small enough that synthetic
   tests cover the format adequately; for scope C, real
   testdata is mandatory.
3. **Reference handling**: the reference-based sequence /
   quality reconstruction requires *external* reference lookup.
   The simplest contract is "the caller provides a reference
   fetcher function"; a more ambitious contract is to embed
   MD5s and resolve via a URI / CRAM-side auxiliary file. The
   ROADMAP does not commit to either. Scope A/B can ship with
   the caller-provided-fetcher shape (matching how
   [htslib](https://github.com/samtools/htslib) handles it).
4. **`.crai` index**: read-only, similar to `.bai` / `.csi`.
   Out of scope for the format commit; flagged for a follow-up.

## What scope A/B would look like as a commit set

Following the per-tier-commit pattern (cf. ROADMAP §Cross-cutting
→ CHANGELOG discipline):

1. `feat(cram): add CRAM v3 read-only parser, validator, and CLI
   integration` — format parsing, decoder, validator, CLI
   extension dispatch + sniff, codec 0/1/2 adapters, clear error
   for codec 3–8.

3. `docs(changelog): add CRAM v3 entry to Unreleased` — entries
   for Added + a note in Style/Changed for the codec limitation.

No separate "tests" commit — the feature commit carries the
public-surface tests in the same change, per the project's
preference for tests-inside-feature-commit.

## What scope A/B would NOT cover (flagged for future work)

- CRAM writer (out of scope per the ROADMAP Index-generation
  open decision).
- Codecs 3–8 (LZMA, rANS 4x8, rANS Nx16, adaptive arith, fqzcomp,
  name tokeniser) — require codec implementations outside
  stdlib.
- `.crai` index reader — small follow-up, same shape as the BAI
  / CSI readers.
- Region queries — same deferral as BAM.
- Lossy CRAM (quality binning, read-name tokenisation) — out of
  scope.
- v4.0 draft features — out of scope; v3 is the current
  standard.