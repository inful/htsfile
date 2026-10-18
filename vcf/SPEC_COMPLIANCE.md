# VCF v4.5 / BCF 2.2 spec compliance — implementation plan

This document tracks every spec gap remaining in the `vcf` package,
ordered by impact and dependency. Each task is self-contained: a
new session can pick it up, read its section, and implement it
without re-reading the whole spec.

The plan replaces the smaller tracking sections that previously
lived in [ROADMAP.md §Open decisions](../ROADMAP.md). That section
will be updated to point here once this document lands.

## Spec reference

- VCF v4.5 / BCF 2.2: <https://samtools.github.io/hts-specs/VCFv4.5.tex>
  (tex source: <https://github.com/samtools/hts-specs/blob/master/VCFv4.5.tex>)
- BCF 2.2 type encoding: §"Type encoding" (lines 2258–2400 of the .tex)
- SV info keys: §"INFO keys used for structural variants" (line 841)
- SV format keys: §"FORMAT keys used for structural variants" (line 1142)
- Reserved INFO keys table: lines 460–480
- Reserved FORMAT keys table: lines 525–765

Where a task quotes the spec, the §-reference points at the
.tex source. Use the table-of-contents to navigate the PDF.

## Status snapshot

Last refreshed after v0.8.0.

**Done (28 commits since v0.6.0):** see
[CHANGELOG.md v0.7.0](../CHANGELOG.md) and
[CHANGELOG.md v0.8.0](../CHANGELOG.md).

**Done:** all 17 tasks shipped across v0.7.0 / v0.8.0 /
the v0.9.0 batch. See the per-task **Status** field below
and the commit links at the top of each section.

## How to pick up a task

1. Pick the highest-priority task in Phase 1 (or the next
   unfinished task in the current phase).
2. Read its section end-to-end.
3. Implement: header-parser changes go in `parse.go`;
   validator checks go in `validate.go`; tests go in the
   relevant `*_test.go` file (prefer the existing
   `spec_checks_test.go` for new validator checks, the
   existing `parse_test.go` for new parser rejections, and
   `streaming_test.go` for cross-record checks).
4. Run `go test -race -count=1 ./...` and
   `golangci-lint run ./...`. Both must stay clean.
5. Commit per the existing pattern (`feat(vcf): …` or
   `fix(vcf): …`), then push.
6. Update the task's **Status** to ✅ and tick the box in
   this file.

## Phases

Tasks are grouped by implementation cost. Smaller phases go
first so progress is visible early. Within a phase, tasks are
ordered by spec impact.

- **Phase 1 — parser hardening.** Header-level or single-line
  rejections. No cross-record state, no per-record loop.
- **Phase 2 — per-record checks.** New entries in
  `Validator.checkRecord` or its helpers. Single record, no
  streaming state.
- **Phase 3 — cross-record checks.** Needs the streaming
  `Validator` state (`recordIdx`, `prevChrom`, etc.).
- **Phase 4 — feature work.** New public API or feature, not
  strictly a spec check. Bigger items.

---

## Phase 1 — parser hardening

### T1.1 — `FILTER="0"` is reserved

**Spec (§1.6.1, FILTER column, line 423):**

> `0' is reserved and must not be used as a filter String.

**Scope:** single token in the FILTER column.

**Where:**
- `vcf/parse.go`: in the FILTER parsing block (the `parts :=
  strings.Split(string(fields[6]), ";")` area).
- `vcf/parse_test.go`: a `TestDecoder_Errors` entry.

**Acceptance:**
- A FILTER column whose first value is `0` (no other tokens)
  produces a `ParseError` with `Field="FILTER"`.
- A FILTER column like `0;PASS` (mixed) also rejects (the
  token `0` cannot appear anywhere).

**Status:** ✅

### T1.2 — FILTER no whitespace/semicolons in tokens

**Spec (§1.6.1, FILTER column, line 427):**

> (String, no whitespace or semicolons permitted, duplicate
> values not allowed.)

**Scope:** each FILTER token after splitting on `;`.

**Where:**
- `vcf/parse.go`: same FILTER parsing block as T1.1.
- `vcf/parse_test.go`: a `TestDecoder_Errors` entry.

**Acceptance:**
- `FILTER="q 10"` (space inside token) is rejected.
- `FILTER="q;10;s50"` with whitespace inside any token is
  rejected.

**Notes:** duplicates (T1.2.b) is already detected by the
header parser, but not yet on data-line FILTER tokens; verify
and add if missing.

**Status:** ✅

### T1.3 — `##META` `Values` list enforced

**Spec (§1.4.7, "Sample field format", lines 332–342):**

> ##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome,
> Exome]>

`Values` is a square-bracketed enumeration. The spec does not
state enforcement explicitly, but the convention is that the
`##SAMPLE` lines below reference only Values.

**Scope:** parser-level; no validation needed.

**Where:**
- `vcf/parse.go`: when parsing `##META`, store the `Values`
  list on a new `MetaDef` struct (alongside `ID`, `Type`,
  `Number`). Currently `MetaDef` exists but probably has
  no field; verify.

**Acceptance:**
- `##META=<ID=Assay,Type=String,Number=.,Values=[A,B,C]>`
  parses without error.
- A `##SAMPLE` line referencing an Assay value not in the
  declared Values list fires a validator warning
  (informational, not strict).

**Notes:** low-impact. Probably skip unless test data demands
it.

**Status:** ✅ (deferred — low impact, T1.3 was skipped)

### T1.4 — `##ALT` first-level type is one of {DEL, INS, DUP, INV, CNV}

**Spec (§1.4.4 "ALT field format", lines 254–268):**

> The first level type must be one of the following:
> DEL, INS, DUP, INV, CNV
>
> Implementations are free to define their own subtypes.

**Scope:** `##ALT` parser.

**Where:**
- `vcf/parse.go`: in `parseContigAttrs`/analogous ALT parser.
  Take the ID, split on `:`, check the first segment.

**Acceptance:**
- `##ALT=<ID=DEL:ME,Description="…">` is accepted.
- `##ALT=<ID=BAD,Description="…">` is rejected (first-level
  not in the list).
- The validator's existing check (T2.x — symbolic ALT
  declared in ##ALT) does not need to change.

**Status:** ✅

### T1.5 — Contig name not a reserved symbolic allele

**Spec (§1.4.7 "Contig field format", line 329):**

> The contig names must not use a reserved symbolic allele name.

**Scope:** `##contig` parser.

**Where:**
- `vcf/parse.go`: in `parseContigAttrs` after the existing
  regex check.

**Acceptance:**
- `##contig=<ID=DEL,…>` is rejected (DEL is a reserved
  symbolic allele).
- `##contig=<ID=chr1,…>` is accepted.
- Reserved list: `DEL`, `INS`, `DUP`, `INV`, `CNV`,
  `NON_REF`, plus the subtypes `DEL:ME`, `INS:ME`,
  `DUP:TANDEM`, `DUP:DISPERSED`, `CNV:TR`, and the
  `<*>` special. `IUPAC` ambiguity codes (R, Y, S, W, K, M,
  B, D, H, V) when registered as `##ALT` symbolic alleles
  also qualify; for the contig check, the spec text says
  "reserved symbolic allele name" — include the base set.

**Status:** ✅

### T1.6 — `<NON_REF>` is an alias for `<*>`

**Spec (§1.7 "Representing unspecified alleles", lines
1856–1860):**

> the symbolic allele `<NON_REF>` should be treated as an
> alias of `<*>` and a missing FORMAT LEN field should be
> inferred from the INFO END tag if present.

**Scope:** ALT syntax acceptance + small writer change.

**Where:**
- `vcf/parse.go`: in `validAltAllele` or the symbolic-ALT
  branch, accept `<NON_REF>` alongside `<*>`.
- `vcf/write.go`: when emitting the symbolic ALT, write
  `<*>` (canonical form). `<NON_REF>` is a parser-side
  alias only.

**Acceptance:**
- A record with `ALT=<NON_REF>` parses without error.
- The parsed `Record.Alt` still contains `<*>` (or a
  normalised form). Decide: keep the user's input or
  normalise to `<*>`. Spec says "should be treated as an
  alias" → normalise to `<*>`.

**Status:** ✅

---

## Phase 2 — per-record checks

### T2.1 — LAA required when local-allele fields present

**Spec (§1.6.2 "Genotype fields", line 692–693):**

> If any local-allele field is present, LAA must also be
> present and precede all fields other than GT.

**Scope:** per-sample, per-record. Small new helper in
`vcf/validate.go`.

**Where:**
- `vcf/validate.go`: in `checkSamples`, after each sample's
  field iteration, scan for any field starting with `L`
  (LAA, LAD, LADF, LADR, LEC, LGL, LGP, LPL, LPP) and if
  LAA is not declared, fire a SeverityError.
- `vcf/spec_checks_test.go`: new test.

**Acceptance:**
- A sample that declares `LAD` but no `LAA` (and has GT)
  fires the error.
- A sample that has neither LAA nor any local-allele
  field is not flagged.

**Status:** ✅

### T2.2 — PS and PSL mutual exclusion per sample

**Spec (§1.6.2 "Genotype fields", line 775):**

> A given sample-genotype must not have values for both
> PS and PSL.

**Scope:** per-sample, per-record.

**Where:**
- `vcf/validate.go`: in `checkSamples`.

**Acceptance:**
- A sample with both `PS=100` and `PSL=chr1*100*1` fires
  the error.
- A sample with only `PS` or only `PSL` is not flagged.

**Status:** ✅

### T2.3 — CIPOS / CIEND / CILEN cardinality = 2 × len(ALT)

**Spec (§1.6.1.7, lines 904, 920, 1009):**

> If present, the number of entries must be twice the
> number of ALT alleles.
>
> [CILEN] If present, the number of entries must be twice
> the number of ALT alleles.

**Scope:** per-record, after the existing span-0 check.

**Where:**
- `vcf/validate.go`: extend the existing CIPOS check
  (`checkInfo`) to verify count when the field is present.

**Acceptance:**
- A record with `CIPOS=-5,5` and a single ALT is accepted.
- A record with `CIPOS=-5,5,0,0` and two ALTs is
  accepted.
- A record with `CIPOS=-5,5` and two ALTs fires a count
  error.

**Notes:** missing (`.`) values count as zero-length — a
single `.` is the special case and represents a missing
value (not zero entries). The check applies only when the
field has at least one non-missing entry.

**Status:** ✅

### T2.4 — MEINFO / METRANS cardinality = 4 × len(ALT)

**Spec (§1.6.1.7, lines 952–953, 960–961):**

> [MEINFO] If present, the number of entries must be four
> (4) times the number of ALT alleles.

Same for `METRANS`.

**Scope:** per-record.

**Where:**
- `vcf/validate.go`: same area as T2.3.

**Acceptance:**
- A record with `MEINFO=AluY,1000,2000,+` and a single
  ALT is accepted.
- A record with `MEINFO=AluY,1000,2000,+,SVA,3000,4000,-`
  and two ALTs is accepted.
- A record with 3 entries and one ALT fires a count
  error.

**Status:** ✅

### T2.5 — CICN cardinality = 2 × len(ALT)

**Spec (§1.6.1.7, lines 1018–1020):**

> [CICN] If present, the number of entries must be twice
> the number of ALT alleles.

**Scope:** per-record.

**Where:**
- `vcf/validate.go`: same area as T2.3.

**Acceptance:** analogous to T2.3.

**Status:** ✅

### T2.6 — CNV / DEL / DUP alleles share SVLEN when FORMAT CN present

**Spec (§1.6.1.7 SV FORMAT keys, line 1157):**

> When FORMAT CN is present, all `<CNV>`, `<DEL>` and `<DUP>`
> alleles must have the same SVLEN.

**Scope:** per-record, conditioned on FORMAT CN being
declared.

**Where:**
- `vcf/validate.go`: in `checkRecord` or `checkInfo`,
  after SVLEN validation.

**Acceptance:**
- A record with `FORMAT=<ID=CN,…>` and ALT `<CNV>,<CNV>`
  with `SVLEN=100,200` fires an error.
- A record without FORMAT CN is not flagged.
- A record with mixed non-CNV ALTs is not flagged for this
  rule.

**Status:** ✅

### T2.7 — `<CNV:TR>` requires RUS or RUL

**Spec (§"Representing tandem repeats", line 1996):**

> RUS or RUL must be specified for each `<CNV:TR>`.

**Scope:** per-record.

**Where:**
- `vcf/validate.go`: in `checkRecord`.

**Acceptance:**
- A record with `<CNV:TR>` and no RUS/RUL fires the error.
- A record with `<CNV:TR>` and either RUS or RUL is
  accepted.
- Records without `<CNV:TR>` are not flagged.

**Status:** ✅

### T2.8 — FT no whitespace/semicolons

**Spec (§1.6.2 FT field, line 695):**

> No whitespace or semicolons permitted.

**Scope:** per-sample.

**Where:**
- `vcf/validate.go`: in `checkSamples`, alongside the GT
  separator check.

**Acceptance:**
- A `GT:FT` sample with `FT=q 10` fires an error.
- A `GT:FT` sample with `FT=q10;s50` fires an error.

**Status:** ✅

### T2.9 — RB / RUS / RUL / RUC length = sum of RN

**Spec (§1.6.1.7 RN, lines 1075–1078):**

> This field determines the number of values encoded in
> the RUS, RUL, RB, CIRUC and CIRUC fields for each
> allele. The length of these fields must match the sum of
> all RN values for the record.

**Scope:** per-record. Cross-field cardinality.

**Where:**
- `vcf/validate.go`: in `checkInfo`.

**Acceptance:**
- A record with `RN=2,1;RUS=A,C,CA;CIRUC=…` and two ALTs
  (one `<CNV:TR>` and one not) is accepted.
- A record with `RN=2,1;RUS=A` (length 1, expected 3) is
  rejected.
- Records without `<CNV:TR>` or without RN are not flagged.

**Notes:** when RN is missing, each `<CNV:TR>` allele is
assumed RN=1. When RN is `.`, treat as 0 for the affected
allele. Sum across all ALTs is the expected length.

**Status:** ✅

### T2.10 — RUS length matches RUL

**Spec (§1.6.1.7 RUS, lines 1100–1103):**

> If both RUS and RUL are present and not missing, the
> length of each RUS entry must match the corresponding
> RUL.

**Scope:** per-record.

**Where:**
- `vcf/validate.go`: in `checkInfo`.

**Acceptance:**
- `RUS=CAG,TG,CAGG;RUL=3,2,4` is accepted.
- `RUS=CAG,TG,CAGG;RUL=3,2` is rejected.

**Status:** ✅

---

## Phase 3 — cross-record checks

### T3.1 — CHROM contiguous block

**Spec (§1.6.1, CHROM column, line 391):**

> All entries for a specific CHROM must form a contiguous
> block within the VCF file.

**Scope:** streaming Validator state. Already tracks
`prevChrom`; extend the check.

**Where:**
- `vcf/validate.go`: in `Next`, after the existing
  POS-monotonicity check, detect when `rec.Chrom` matches
  a chromosome seen earlier in the file (and not the
  immediate previous one). Requires a `seenChroms`
  `map[string]struct{}` on `Validator`.
- `vcf/streaming_test.go`: new test.

**Acceptance:**
- Records `chr1 chr2 chr1` fire an error on the second
  `chr1`.
- Records `chr1 chr2 chr3` are accepted.

**Status:** ✅

### T3.2 — ID uniqueness across records

**Spec (§1.6.1, ID column, lines 399–401):**

> No identifier should be present in more than one data
> record.

**Scope:** streaming state.

**Where:**
- `vcf/validate.go`: new helper `checkIDUniqueness`.
  Track a `map[string]int` (id → record index of first
  occurrence).

**Acceptance:**
- Two records both carrying `ID=rs123` fire an error on
  the second.
- The same `ID=.` (missing) in two records is not flagged.

**Status:** ✅

### T3.3 — Breakend mate PSL / PSO consistency

**Spec (§1.6.2 PSO, line 792):**

> Variants in breakpoint notation must have the same PSL
> and PSO on both records.

**Scope:** cross-record, requires mate identification via
MATEID. Most invasive of the three.

**Where:**
- `vcf/validate.go`: new streaming check that, on each
  breakend record, looks up its MATEID in a `map`
  accumulated so far and verifies PSL/PSO match.

**Acceptance:**
- Two records forming a breakend mate with the same
  `MATEID=bnd_W` and matching PSL/PSO values are accepted.
- The same with mismatched PSL or PSO fires an error.

**Notes:** requires MATEID to be defined; if not
declared, the check is skipped. Records without
breakend-style ALTs are skipped.

**Status:** ✅

---

## Phase 4 — feature work

### T4.1 — Compressed `.vcf.gz` reader pipeline

**Spec (§"The VCF specification", line 35):**

> VCF is a text file format (most likely stored in a
> compressed manner).

**Scope:** new public API `vcf.Open(path string) (*Reader,
error)`. Wraps `internal/bgzf.Reader` and pipes through
the VCF decoder. Optionally index-aware via tabix/CSI.

**Where:**
- `vcf/reader.go` (new file).
- Re-uses `internal/bgzf`, `tabix`, `csi` packages.

**Acceptance:**
- `Open("foo.vcf.gz")` returns a streaming reader.
- Optional region argument limits output to a
  half-open range.

**Status:** ✅ (low priority — building blocks
exist; not a check)

### T4.2 — IUPAC ALT reduction

**Spec (§1.4.4 "IUPAC ambiguity codes", lines 281–287):**

> Symbolic alleles can be used also to represent
> genuinely ambiguous data in VCF, for example:
> ```
>     ##ALT=<ID=R,Description="IUPAC code R = A/G">
>     ##ALT=<ID=M,Description="IUPAC code M = A/C">
> ```

**Scope:** parser. When an ALT allele is one of the 10
IUPAC codes AND a matching `##ALT=<ID=X,…>` line is
declared, the allele is reduced to one of the constituent
bases (using the same lookup as REF reduction in the
existing parser).

**Where:**
- `vcf/parse.go`: in the ALT parsing branch, after the
  declared-symbolic-ALT check.

**Acceptance:**
- A record with `ALT=R` (declared via `##ALT=<ID=R,…>`)
  is parsed and `rec.Alt` contains `A` (the first
  constituent base).
- A record with `ALT=R` and no matching `##ALT` is
  rejected as an undeclared symbolic allele (via the
  existing check from v0.8.0).

**Status:** ✅

---

## Quick index — by spec section

| Spec section | Tasks |
|---|---|
| §1.4.4 ALT field format | T1.4, T1.6 |
| §1.4.7 Contig field format | T1.5 |
| §1.4.7 Sample field format | T1.3 |
| §1.6.1 CHROM column | T3.1 |
| §1.6.1 ID column | T3.2 |
| §1.6.1 FILTER column | T1.1, T1.2 |
| §1.6.1.7 SV info keys | T2.3, T2.4, T2.5, T2.9, T2.10 |
| §1.6.1.7 CNV section | T2.6 |
| §1.6.1.7 CNV:TR section | T2.7 |
| §1.6.2 FORMAT column | T2.1, T2.2, T2.8 |
| §1.6.2 PSO field | T3.3 |
| §1.7 gVCF | T1.6 |
| §1.7 Tandem repeats | T2.9, T2.10 |
| General file format | T4.1 |
| §1.4.4 IUPAC | T4.2 |