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

**Status:** ✅

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

### T4.3 — Reference-aware REF/ALT padding check (§1.6.1.1)

**Spec (§1.6.1.1, lines 406–407):**

> For simple insertions and deletions in which either
> the REF or one of the ALT alleles would otherwise be
> null/empty, the REF and ALT Strings must include the
> base before the variant (which must be reflected in
> the POS field), unless the variant occurs at position
> 1 on the contig in which case it must include the base
> after the variant.

**Scope:** per-record, conditioned on a FASTA reference
being attached. The check fires only when padding is
required: any ALT is symbolic, or any ALT has a
different length from REF (true indel). SNVs and
same-length substitutions are exempt — the spec allows
padding to be omitted for "complex substitutions or
other variants where all alleles have at least one base
represented".

**Where:**
- `fasta/fasta.go` (new package, v0.10.0): plain
  FASTA reader with optional `.fai` random access.
- `vcf/validate.go`: `Validator.WithReference` and
  `checkPaddingBase` look up `ref[POS]` and compare to
  `REF[0]`.

**Acceptance:**
- A record with `POS=3`, `REF="G"`, `ALT="GX"` on a
  reference where `ref[3]="G"` is accepted.
- A record with the same POS/ALT but `REF="A"` fires
  SeverityError.
- An SNV with `REF="A"`, `ALT="T"` (same length) is
  accepted regardless of `ref[POS]`.
- Unknown contig or position-out-of-range degrades to
  SeverityWarning (skipped in ModeDefault).

**Notes:** plain FASTA only in this release. BGZF (the
common `.fa.gz` form) is not supported; callers should
gunzip first. The `.fai` format follows the samtools
faidx output: 5-column TSV
(`name`, `length`, `offset`, `linebases`, `linewidth`).
The offset is the byte position of the first base
(just past the header line).

**Status:** ✅ (v0.10.0)

### T4.4 — Full REF match against reference (§1.6.1.1, supersedes v0.10.0)

**Spec (§1.6.1.1):**

> the REF and ALT Strings must include the base
> before the variant (which must be reflected in
> the POS field)

**Scope:** per-record, conditioned on a FASTA
reference being attached. For every byte `i` of
REF, the validator compares `REF[i]` against
`ref[POS+i]`. One SeverityError per mismatch.

**Where:**
- `fasta/fasta.go`: new `(*Reference).BaseRange` for
  efficient multi-base lookup.
- `vcf/validate.go`: `checkRefMatchesReference`
  replaces the v0.10.0 single-base
  `checkPaddingBase`. The call site in `checkRecord`
  is a one-liner; no `needsPadding` gate (the full
  check fires for any non-empty REF).

**Acceptance:**
- SNV `REF="G" ALT="A"` at `POS=3` against
  `ref[3]="G"` passes.
- SNV `REF="T" ALT="A"` at `POS=3` against
  `ref[3]="G"` fires (REF mismatch). The v0.10.0
  padding check would NOT have caught this; the new
  full check does.
- Insertion at `POS=3` `REF="G" ALT="GX"` against
  `ref[3]="G"` passes.
- Multi-base substitution at `POS=3` `REF="GT"
  ALT="AC"` against `ref[3..5]="GT"` passes.
- Symbolic `<DEL>` at `POS=2` with `REF="AG"` against
  `ref[2]="A"`, `ref[3]="G"` passes.
- Unknown contig or out-of-range position degrades
  to SeverityWarning (skipped in ModeDefault).

**Notes:** `BaseRange` is one forward scan instead of
N independent seeks; it correctly skips comment
lines and internal whitespace in the scanned path.

**Status:** ✅ (v0.11.0)

### T4.5 — Reference-aware END / SVLEN bounds (§1.6.1.6, §1.6.1.7)

**Spec (§1.6.1.6):**

> END is a computed field that, when present, must
> be set to the maximum end reference position
> (1-based).

**Spec (§1.6.1.7):**

> SVLEN must be specified for symbolic structural
> variant alleles.

**Scope:** per-record, conditioned on a FASTA
reference being attached. Without a reference, the
existing `checkEND` already enforces value
consistency (END == computed from REF / SVLEN),
and the existing positive-value check handles
negative SVLEN. The new bounds checks verify that
END and SVLEN stay within the contig even when no
`##contig=<length>` declaration is present.

**Where:**
- `fasta/fasta.go`: new `(*Reference).Length(chrom)`
  returning the contig length (or an error).
- `vcf/validate.go`: `checkENDBounds` and
  `checkSVLENBounds` wired into `checkRecord`
  after the existing `checkEND`. Both are no-ops
  without a reference; both skip mitochondrial
  chromosomes (matching the existing exception).

**Acceptance:**
- `POS=3, END=6, SVLEN=4` against `ref.Length=16`
  passes.
- `POS=3, END=100, SVLEN=98` against `ref.Length=16`
  fires `END exceeds contig length`.
- `POS=3, END=2, SVLEN=-1` (which the existing
  positive-value check rejects anyway) fires
  `END before POS`.
- `POS=3, SVLEN=100` against `ref.Length=16` fires
  `SVLEN implies end position past contig`.
- `POS=1, SVLEN=16` against `ref.Length=16` passes
  (end equals last base).
- `chrM` (mitochondrial) skips both checks even
  when bounds would be violated.
- Unknown contig: lookup error is silent; check is
  skipped.

**Notes:** bounds checks use absolute SVLEN per the
spec's "negative SVLEN should be treated as positive"
rule. The 0 sentinel in the Ints slice (parser
convention for missing) is filtered out via
`isMissing`.

**Status:** ✅ (v0.12.0)

### T4.6 — Breakend remote-position bounds (§1.6.5.1)

**Spec (§1.6.5.1, "Specifying complex rearrangements
with breakends"):**

> The position p of the mate breakend, indicated by
> a string of the form "chr:pos".

The four valid breakend shapes are `t]p]`, `t[p[`,
`]p]t`, `[p[t`. The position is bracketed by two
matching brackets (`[` or `]`).

**Scope:** per-record, conditioned on a FASTA
reference being attached. For every breakend ALT,
parse the remote position and verify it lies within
the remote contig's length. Without a reference,
this check is skipped.

**Where:**
- `vcf/validate.go`: `parseBreakendPosition` parses
  the four breakend forms plus the optional
  trailing strand-bracket. `checkBreakendPositions`
  fires `SeverityError` when `pos > ref.Length` for
  the remote contig; unknown remote contigs degrade
  to SeverityWarning in ModeDefault, SeverityError in
  ModeStrict. Telomere breakends (no position) and
  non-breakend ALTs are skipped.

**Acceptance:**
- `A]chr1:10]` (t]p] form, pos=10) against
  `ref.Length(chr1)=16` passes.
- `A]chr1:100]` (pos past contig) fires
  `breakend position exceeds contig length`.
- `A]chrX:100]` (unknown remote contig) fires
  warning (default) / error (strict).
- `]chr1]A` (telomere, no position) is skipped.
- `T` (non-breakend ALT) is skipped.

**Notes:** breakend ALTs are written as literal
strings (e.g. `A]chr1:10]`), not as symbolic IDs
in `<>` form. The local base character is part of
the replacement string per the spec, not a
declared ALT type.

**Status:** ✅ (v0.13.0)

### T4.7 — fasta performance + BGZF (v0.14.0, no vcf change)

This is a perf/feature release in the `fasta` package;
no new vcf checks land here.

**Scope:** the `fasta` package now supports BGZF
(`.fa.gz`) references and an internal cursor that
makes sequential per-record lookups O(N) instead of
O(N²).

**Where:**
- `fasta/fasta.go`: `Open` detects the BGZF magic
  (`0x1f 0x8b 0x08` + FEXTRA flag) and decompresses
  into memory. `Reference` carries `cursorCh`,
  `cursorOffset`, and `cursorBasesSeen` for the
  forward-scan fast path.

**Notes:** the `.fai` index is ignored for BGZF
because its virtual offsets require a seekable BGZF
reader, which this package does not currently provide.
For BGZF references, Open pays a memory cost equal to
the decompressed size.

The IUPAC ALT reduction cross-check originally listed
here as a candidate v0.14.0 item was dropped: the
parser does not record which ALTs were originally
IUPAC codes, so we cannot distinguish a reduced-IUPAC
base from a concrete call. The full-REF match
(v0.11.0) already covers the case where the
reduction is inconsistent with the reference.

**Status:** ✅ (v0.14.0)

### T4.8 — Tabix/CSI indexed VCF reads + CIGAR validation (v0.15.0)

**Spec (§1.6.1.1, "CIGAR" in INFO):**

> a cigar string used in the original alignment

The CIGAR string describes how to align each ALT to
the reference. The sum of REF-consuming op lengths
(M, =, X, D, N) must equal `len(REF)`, and the sum of
ALT-consuming op lengths (M, =, X, I, S, H) must
equal `len(ALT)`. Padding (P) consumes nothing.

**Scope:** per-record, always fires (no reference
required). CIGAR is checked independently of
reference attachment.

**Where:**
- `vcf/validate.go`: `checkCIGAR` extracted from
  `checkInfo` (cognitive-complexity). `parseCIGAR`
  helper.
- `internal/bgzf/seek.go`: `SeekReader` for random
  access via virtual offsets.
- `vcf/tabix.go` + `vcf/csi.go`: `OpenTabix`,
  `OpenCSI`, `Query`.
- `vcf/tabix_test.go` + `vcf/csi_test.go`: 9 tests.

**Acceptance:**
- `REF=ACGT ALT=A CIGAR=1M3D`: 4 reference bases
  consumed, 1 query; passes.
- `REF=A ALT=ACGT CIGAR=1M3I`: 1 reference, 4 query;
  passes.
- `REF=A ALT=T CIGAR=2M`: 2 reference but REF has
  length 1; error.
- `CIGAR=1` (missing op): malformed; error.
- `CIGAR=1Z` (unknown op): malformed; error.

**Notes:** for v0.15.0 there is no CLI wrapper; the
package-level API is the only surface. Indexed VCF
reads use the new BGZF `SeekReader`; plain
(non-indexed) reads are unchanged.

**Status:** ✅ (v0.15.0)

### T4.9 — gVCF non-variant block rules (§1.7)

**Spec (§1.7 gVCF, "Genomic Variant Call Format"):**

> A VCF file may be provided in a gVCF format... a
> block record uses the <NON_REF> symbolic allele (which is
> an alias for <*>) to represent a stretch of
> reference-matching positions from POS through END.

Per the spec a "non-variant block" record is identified by:

- ALT contains exactly `<*>` (or its `<NON_REF>` alias,
  which the parser already normalises to `<*>` per T1.6).
- The header declares `##INFO=<ID=END,...>`, and the
  record carries an END value.

Rules to enforce:

- ALT is `<*>` (single-element ALT) ⇒ END must be
  declared in INFO. Without END the block is
  malformed.
- END must be ≥ POS (the block cannot go
  backwards).
- A non-variant block must not overlap with the
  previous record on the same CHROM. That is, the
  next record's POS must be > previous record's
  effective "last position" — `END` for a previous
  block, `POS` for a previous variant record.
- (Reference-attached only) The REF column of a
  non-variant block must match the reference
  sequence over `[POS, END]`.

**Scope:** per-record (single) for rules 1-3;
reference-aware for rule 4.

**Where:**
- `vcf/validate.go`: `checkGVCFBlock` per-record
  checks (rules 1-3) plus the cross-record
  monotonic-overlap check that uses the existing
  `prevPos` machinery with a new `prevEnd` field.
- `vcf/validate.go`: `checkGVCFBlockRef` for rule 4,
  analogous to the existing `checkRefMatchesReference`
  but reading `[POS, END]` instead of `[POS, len(REF)]`.
- `vcf/spec_checks_test.go`: new sub-tests for each
  rule.

**Acceptance:**
- `chr1 100 . A <NON_REF> . . END=200`: block from
  100..200; passes.
- `chr1 100 . A <NON_REF> . .` (no END): error
  (ModeStrict) / warning (ModeDefault).
- `chr1 100 . A <NON_REF> . . END=50`: END < POS;
  error.
- Two records on the same contig with
  `chr1 100 . A <NON_REF> . . END=200` followed by
  `chr1 150 . A G ...`: overlap; error.
- With a reference attached, REF column must
  match over the block range; mismatch error.

**Status:** ✅ (v0.15.5)

### T4.10 — Breakend mate position round-trip (§1.6.5.1)

**Spec (§1.6.5.1, "Specifying complex rearrangements
with breakends"):**

> The position p of the mate breakend, indicated by
> a string of the form "chr:pos".

In addition to T4.6's per-record bound check (the
remote position must lie within the remote contig),
the **two records in a breakend pair must point to
each other**. That is, if record A is at (chrA, posA)
with breakend ALT `t]chrB:posB]`, then record B (its
MATEID mate) must be at exactly (chrB, posB) — and
symmetrically B's breakend ALT must point back to
A's (chrA, posA).

**Scope:** cross-record, requires MATEID.

**Where:**
- `vcf/validate.go`: extend `mateSnapshot` with
  the local position and parsed breakend remote
  positions; new `checkBreakendMatePositions`
  compares the mate's breakends to the partner's
  local position.
- `vcf/spec_checks_test.go`: new
  `TestValidator_BreakendMatePositionRoundTrip`.

**Acceptance:**
- A pair `chr1 100 . N ]chr2:200]N . . MATEID=bnd1`
  followed by `chr2 200 . N ]chr1:100]N . .
  MATEID=bnd1`: mates point to each other; passes.
- Same with one mate's POS wrong (`chr2 201 . N
  ]chr1:100]N ...`): error in both modes.
- Same with one mate's breakend remote position
  wrong (`chr2 200 . N ]chr1:999]N ...`): error.
- Same with the remote contig wrong (`chr2 200 . N
  ]chr3:100]N ...`): error.
- Records without breakend ALTs (e.g. SNVs with a
  MATEID linking paired SNVs) are skipped — the
  position round-trip only applies to breakends.

**Status:** ✅ (v0.15.6)

## Spec audit

This section is a one-shot audit of every section of the
VCF v4.5 / BCF 2.2 specification, with status against
the `vcf` (and supporting `bcf`, `fasta`, `tabix`,
`csi`) packages.

A ✅ means the rule is fully implemented in the parser
or validator and exercised by tests. A ❓ means a known
gap. Numbers in parentheses are the line numbers in the
.tex source at
`<https://github.com/samtools/hts-specs/blob/master/VCFv4.5.tex>`.

### §1.4.1 Character encoding

- ✅ UTF-8 declared (line 73).
- ✅ Non-printable bytes U+0000–U+0008, U+000B–U+000C,
  U+000E–U+001F rejected (T1.1 parser hardening).
- ✅ LF and CR+LF line separators accepted (line 78).
- ❓ Percent encoding (`%3A`, `%3B`, …) for special chars
  inside INFO/FORMAT values is decoded by the parser
  but not enforced on input — we accept literal `:`
  inside an INFO value too.

### §1.4.2 Data types

- ✅ Integer / Float / Flag / Character / String types
  (line 97).
- ❓ **BCF Integer range $-2^{31}$ to $-2^{31}+7$ is
  disallowed** (line 100). Not enforced — VCF and BCF
  share the same parser path, and we accept the
  disallowed range.
- ❓ **Float format** must match one of the two
  regexes; `parseVCFFloat` uses `strconv.ParseFloat`
  which is more permissive than the spec regex.
- ❓ **Float INF/INFINITY/NAN parsing** matches the spec
  intent, but should be tied to a dedicated test
  (currently exercised indirectly).

### §1.4.3 Meta-information lines

- ✅ `##fileformat` must be first and unique.
- ✅ `##INFO`, `##FILTER`, `##FORMAT` parsing.
- ✅ `##ALT` first-level type (T1.4).
- ✅ `##contig` name validation (T1.5).
- ✅ `##META` `Values` list (T1.3).
- ❓ `##pedigreeDB` is parsed but not validated as a URL.
- ❓ `##assembly` is parsed but not validated as a URL.
- ❓ `##PEDIGREE` parses `<ID=…,Original=…>` etc. but
  does not cross-check that `Father`/`Mother` reference
  declared sample IDs.

### §1.4.4 Header line

- ✅ Eight fixed columns present.
- ❓ **Duplicate sample IDs in the header line** (line
  365: "Duplicate sample IDs are not allowed") is
  parsed but not rejected.

### §1.6.1 Fixed fields

- ✅ `CHROM` no whitespace (line 384).
- ✅ `POS` numeric, 1-based, telomere positions 0/N+1
  accepted.
- ✅ `POS` sorted within CHROM (T3.1, ModeStrict).
- ✅ `ID` no whitespace/semicolons/equals (T2.x).
- ✅ `ID` unique across records (T3.2).
- ❓ **Duplicate IDs within a single record's ID
  field** (line 401: "duplicate values not allowed") is
  not enforced.
- ✅ `REF` bases A/C/G/T/N (T2.x in ModeStrict).
- ✅ `REF` IUPAC reduction (T4.2).
- ✅ `ALT` rules (T1.4, T1.6, T2.x).
- ✅ `QUAL` is Phred-scaled.
- ✅ `FILTER` token rules (T1.1, T1.2).
- ❓ **INFO/FORMAT key regex**
  `^([A-Za-z_][0-9A-Za-z_.]*|1000G)$` is enforced for
  declared keys (T1.6.1 parser hardening covered it),
  but not for record-level keys in the parser path.

### §1.6.1.1 Reserved INFO keys

- ✅ All standard keys (AA, AC, AD, ADF, ADR, AF, AN, BQ,
  CIGAR, DB, DP, END, H2, H3, MQ, MQ0, NS, SB, SOMATIC,
  VALIDATED, 1000G).
- ✅ END computed value check (existing T2.x).
- ✅ END/SVLEN bounds (T4.5).

### §1.6.2 Genotype fields

- ✅ `GT` first FORMAT field if present (existing check).
- ✅ `LAA` precedes all non-GT fields.
- ✅ Trailing FORMAT fields dropped.
- ✅ PS / PSL / PSO / PSQ (T2.2, T3.3).
- ❓ **Duplicate FORMAT keys within a record**
  (line 469: "duplicate keys are not allowed") is not
  enforced.
- ❓ **GT allele values 0..N** are parsed but not range-
  checked; a record with ALT=`A` and `GT=2` would parse.
- ❓ **Haploid GT must have one allele** (line 553) is
  not checked.
- ❓ **LA** reserved "for future use" — we silently accept
  it without complaint.
- ❓ **LR, LG** cardinalities on local-allele equivalents
  are checked through the matching field's `Number` but
  not via a dedicated LR/LG count check.

### §1.6.5 SV INFO keys

- ✅ IMPRECISE, NOVEL, END, SVTYPE, SVLEN, CIPOS, CIEND,
  HOMLEN, HOMSEQ, BKPTID, MEINFO (T2.4), METRANS (T2.4),
  CILEN (T2.3), CN, CICN (T2.5), SVCLAIM (T2.7),
  MATEID, PARID, EVENT, EVENTTYPE, RN, RUS, RUL, RUC, RB
  (T2.9, T2.10), CIRUC, CIRB, RUB.

### §1.6.5 FORMAT keys for SV

- ✅ CN, CICN, CNQ, CNL, CNP, NQ, HAP, AHAP.

### §1.6.5.1 Breakends

- ✅ Four breakend forms parsed.
- ✅ Remote position bounds (T4.6).
- ✅ PSL/PSO consistency (T3.3).
- ✅ Mate position round-trip (T4.10).
- ❓ **PARID** declared in the spec but our index only
  handles MATEID. PARID is the partner's local ID — the
  same handling as MATEID would apply.

### §1.7 gVCF

- ✅ T4.9 (gVCF non-variant block rules).

### §1.7.1 CNV

- ✅ T2.6 (CNV/DEL/DUP alleles share SVLEN when FORMAT CN).
- ✅ T2.7 (CNV:TR requires RUS or RUL).

### §1.7.2 Tandem repeats

- ✅ T2.9 (RB/RUS/RUL/RUC length = sum of RN).
- ✅ T2.10 (RUS length matches RUL).
- ❓ **CIRB confidence interval bounds span zero**
  (line 1119) is not checked.

### §5 BCF2

The `bcf` package already implements most of the BCF2
binary type encoding (§BcfTypeEncoding). The v0.15.7
audit wording ("type encoding not implemented") was
written when only the v1 text-based reader existed; by
v0.15.7 the v2 reader/writer had landed.

Implemented in the v2 reader/writer:
- BGZF container.
- 6-byte magic + length-prefixed VCF-style header.
- BCF 2.1 (default) and BCF 2.2 (long-string `0xF`
  type tag) magic accept.
- CHROM (int32), POS (int32), ID/REF/ALT
  (typed-string), FILTER (typed), QUAL (typed float64),
  INFO (typed per declaration), per-sample FORMAT
  (typed values).
- Type-byte vocabulary: `typeMissing=0`, `typeInt8=1`,
  `typeInt16=2`, `typeInt32=3`, `typeFloat=5`,
  `typeString=7`, `typeLong=0xF` (BCF 2.2 long strings).
- Per-spec missing-value sentinels for FILTER and length
  prefixes (`bcfFilterMissing = bcfInt32Missing =
  0xFFFFFFFF`).

Remaining spec-conformance gaps:
1. ✅ **Atomic type byte inline size** (line 2271).
   Closed by T4.21, v0.15.11. The writer now emits
   `(1<<4) | type` for atomic Integer/Float; the
   reader masks the high nibble to accept both forms.
2. ✅ **Float MISSING encoding** (line 2333). Closed
   by T4.24, v0.16.1. Per-type MISSING for a Float
   slot is the signaling NaN `0x7F800001`; for an
   Integer slot it is `0x80000000`. The writer
   emits these sentinels when a sample's value is
   Missing; the reader detects them and surfaces
   Missing=true on the per-sample slot.
3. ✅ **END_OF_VECTOR** (line 2304–2305, 2325–2326).
   Closed by T4.23, v0.16.0. The encoder pads shorter
   entries with `0x81` (int8) or `0x80000001` (int32)
   or `0x7F800002` (float32); the reader strips them.
4. ✅ **Per-field FORMAT encoding** (line 2236–2242).
   Closed by T4.23, v0.16.0. Each FORMAT key now has
   exactly one type byte for the whole N-sample vector.
5. ✅ **Integer reserved range on writer side** (line
   2306). Closed by T4.22, v0.15.11.

§5 BCF2 is now spec-compliant.

### Spec rules exercised by tests but currently
**outside SPEC_COMPLIANCE.md's task taxonomy**

The following checks live in `validate.go` but aren't
called out as their own SPEC_COMPLIANCE task:

- `checkCardinality` for INFO/FORMAT (line 696).
- `checkPloidy` and `checkHaploidGenotype` (existing
  ploidy-consistency rules).
- `checkMitochondrialChromExemptions` for mtDNA.
- `checkFTNoWhitespaceOrSemicolons` (T2.8).
- `checkPSLUnphased` and `checkPSLTrio`.

### Highest-impact remaining gaps

In priority order (most spec-rule-fixing value):

1. ✅ **FORMAT column duplicate keys within a record**
   (closed by T4.11, v0.15.7).
2. ✅ **ID column duplicate values within a record**
   (closed by T4.11, v0.15.7).
3. ✅ **FORMAT key regex enforcement** (closed by T4.13,
   v0.15.8).
4. ✅ **CIRB bounds span zero** (closed by T4.12, v0.15.7).
5. ✅ **BCF2 forbidden Integer range** (closed by T4.14,
   v0.15.8).
6. ✅ **mtDNA haploid GT length** (closed by T4.15,
   v0.15.8).
7. ✅ **PEDIGREE Father/Mother/Original resolution**
   (closed by T4.16, v0.15.8).

Remaining:

1. ✅ **Per-field FORMAT encoding + END_OF_VECTOR**
   (§5, line 2208–2255, 2304, 2325). Closed by T4.23,
   v0.16.0. The encoder now emits the per-field block
   (one type byte per field, N values per field padded
   with END_OF_VECTOR); the reader strips the
   padding so each sample sees only its own
   cardinality.
2. ✅ **Inline size in type byte** (§BcfTypeEncoding,
   line 2271). Closed by T4.21, v0.15.11. The
   writer now emits `(1<<4) | type` for atomic
   Integer and Float; the reader masks the high
   nibble to accept both the bare and inline-size
   form.
3. ✅ **Integer reserved range on writer**
   (§1.4.2, line 2306). Closed by T4.22,
   v0.15.11. The writer refuses to emit values
   in `[-2^31, -2^31+7]`.

§5 BCF2 is now substantially spec-compliant.

### T4.24 — Per-type MISSING sentinels (§BcfTypeEncoding, line 2304, 2333)

**Spec (§BcfTypeEncoding, line 2304):**

> For integer types, `0x80000000` is interpreted as a
> MISSING value (line 2304).

**Spec (line 2333):**

> MISSING & 0x7F800001

The previous writer (v0.16.0) encoded "missing" floats
or integers via a size-0 byte (`typeMissing`) — spec-
compliant but not the per-type MISSING form. Some
external tools emit the per-type MISSING sentinel and
expect readers to recognise it.

**Where:**
- `bcf/end_of_vector.go`: new constants
  `bcfFloatMissing = 0x7F800001` and
  `bcfIntMissingInt32 = 0x80000000` keep the spec's
  reserved sentinels next to the existing
  END_OF_VECTOR constants.
- `bcf/format.go`:
  - the writer emits per-type MISSING sentinels for
    Integer/Float slots, keeping the field's value
    count uniform across samples;
  - the reader detects either size-0 typeMissing
    (whole-field collapse) or per-type MISSING
    (one sample slot only) and surfaces Missing=true
    on the corresponding sample.
- Tests: `TestFormatFloatMissingSignalingNaN`,
  `TestFormatIntMissingSentinel`.

**Acceptance:**
- Mixed missing/present Float field (one sample is
  Missing, one has floats): the missing sample's slot
  reads back as Missing=true with no Floats; the
  present sample's slot reads back as Floats=[0.1,
  0.2].
- Mixed missing/present Integer field: same.

**Status:** ✅ (v0.16.1) (§5 line 2208–2255, 2304, 2325)

**Spec (§5, line 2208–2255):**

> "Genotype fields are encoded not by sample as in VCF
> but rather by field, with a vector of values for each
> sample following each field."

**Spec (§BcfTypeEncoding, line 2304–2305):**

> For integer types, the values `0x81`, `0x8001`,
> `0x80000001` are END_OF_VECTOR indicators. For
> float, `0x7F800002`.

The v0.15.x FORMAT block was per-sample (one type byte
per sample per key) which works for self-roundtrip but
is not the spec's compact form. **T4.23 restructures the
block to be per-field**: for each FORMAT key, the
encoder emits exactly one type byte (size in high nibble,
type in low) plus N typed values, padded with
END_OF_VECTOR sentinels for shorter entries. The reader
strips the padding so each sample sees only its own
cardinality.

GT is stored in the per-field block as a typed int8
vector (the spec's wire form at line 2463): each allele
encodes as `(allele+1)<<1 | phased`. Mixed-ploidy GTs
(e.g. chrX-style haploid + diploid in one record) are
padded with the int8 END_OF_VECTOR (`0x81`) for the
shorter entries.

Where:
- `bcf/gt.go`: `encodeGT(gt string) []uint8` and
  `decodeGT(ints []uint8) (string, error)` implement the
  spec's int-vector encoding and its inverse. The phasing
  rule attaches the phased bit to the trailing allele of
  each pair.
- `bcf/end_of_vector.go`: sentinel constants
  `bcfEndOfVectorInt8 = 0x81`, `bcfEndOfVectorInt32 =
  0x80000001`, `bcfEndOfVectorFloat32 = 0x7F800002`.
- `bcf/format.go`: per-field reader
  `(*Decoder).readFormatBlock` and writer
  `writeFormatBlock`. `fieldMaxSize` computes the
  uniform vector size; `writeFormatSampleValue`
  pads shorter entries with the appropriate
  END_OF_VECTOR.
- Tests: `TestEncodeGT`, `TestDecodeGT`,
  `TestGTRoundTrip`, `TestFormatPerFieldMixedPloidy`,
  `TestFormatPerFieldAllMissing`.

Acceptance:
- "0/1" → int vector `[0x02, 0x04]` and back.
- "0|1" → `[0x03, 0x05]`.
- "./." → `[0x00, 0x00]`.
- Mixed-ploidy (`S1=0`, `S2=0/1`): the encoder pads S1
  with END_OF_VECTOR, the decoder strips it, both samples
  recover their own GT cardinality.
- All-missing fields collapse to a single `typeMissing`
  byte for the whole field.

Status: ✅ (v0.16.0)

### T4.17 — LAA cross-check (§1.6.2)

The cross-check ensures LA-field values are
interpreted as indices into LAA, not the full ALT
list. With `LAA=[2,4]` (2 elements) and `GT=0/1`, an
LGL field with `len(LGL)==6` is wrong: under the
refactored rules the cardinality for an LA-format
`Number=G` field is `C(len(LAA)+ploidy−1, ploidy)`,
not `C(len(ALT)+ploidy, ploidy)`.

**Where:**
- `vcf/validate.go`:
  - new `wantedCardinality` helper picks `len(LAA)`
    when the key is a `localAlleleKeys` entry and
    the LA count is non-zero; falls back to
    `len(rec.Alt)` when LAA is absent.
  - For `Number=G` LA fields the helper subtracts 1
    from `cardinalAlt` because the conventional
    formula adds REF back in, which doesn't apply
    when LAA is the scope.
  - `laaLenFor` looks at `Ints` / `Strs` / `Raw` to
    handle the spec's `LAA Type=String` default.
- `vcf/spec_checks_test.go`:
  `TestValidator_LAACrossCheck` with 4 sub-tests
  covering LEC with two values, LEC with four
  values (error), LGL with three values (match),
  LGL with five values (error).

**Acceptance:**
- `LAA=[2,4]`, `LEC=10,20`: cardinality 2 → OK.
- `LAA=[2,4]`, `LEC=10,20,30,40`: cardinality 2 →
  error.
- `LAA=[2,4]`, ploidy 2, `LGL=0.1,0.2,0.3`:
  cardinality 3 → OK.
- `LAA=[2,4]`, ploidy 2, `LGL=<5 vals>`: cardinality
  3 → error.

**Status:** ✅ (v0.15.9)

### T4.11 — Per-record ID/FORMAT duplicate checks (§1.6.1, §1.6.2)

**Spec (§1.6.1 ID column, line 401):**

> No identifier should be present in more than one
> data record.
> (String, no whitespace or semicolons permitted,
> duplicate values not allowed.)

**Spec (§1.6.2 FORMAT keys, line 469):**

> The first key must always be the genotype (GT) if it
> is present.
> […] duplicate keys are not allowed.

The header-level checks (T2.x parser hardening) reject
duplicate `##INFO`, `##FILTER`, and `##FORMAT` IDs.
The cross-record ID uniqueness check (T3.2) rejects the
same ID on two separate records. The remaining two cases
— duplicate IDs within a single record's ID field, and
duplicate keys within a single record's FORMAT column —
were not enforced; this task closes them.

**Where:**
- `vcf/parse.go`: `parseRecord` walks the ID field
  with a `seenID` map and the FORMAT column with a
  `seenFmt` map, both producing `ParseError` on a
  duplicate.
- `vcf/parse_test.go`: two new entries in
  `TestDecoder_Errors` cover "duplicate ID within
  record" and "duplicate FORMAT key within record".

**Acceptance:**
- A record with `ID=rs1;rs1` is rejected.
- A record with `FORMAT=GT:GT` is rejected.
- A record with `ID=rs1;rs2` (unique within the record)
  parses fine.

**Status:** ✅ (v0.15.7)

### T4.12 — CIRUC/CIRB bounds span zero (§1.6.1.7)

**Spec (§1.6.1.7 CIRUC/CIRB, line 1119):**

> The lower bound must be less than or equal to zero
> and the upper bound must be greater than or equal to
> zero.

**Scope:** per-record.

**Where:**
- `vcf/validate.go`: `checkRNDeterminedFields` walks
  each CIRUC/CIRB pair and emits a `SeverityError`
  per violation.
- `vcf/spec_checks_test.go`:
  `TestValidator_CIRUCBIRBBoundsSpanZero` with 5
  sub-tests (valid CIRUC, lower/upper bound violations
  for each of CIRUC and CIRB).

**Acceptance:**
- CIRUC=-5,5 valid.
- CIRUC=1,5 (lower positive): error.
- CIRUC=-5,-1 (upper negative): error.
- CIRB=2,5 (lower positive): error.
- CIRB=-5,-2 (upper negative): error.

**Status:** ✅ (v0.15.7)

### T4.25 — Reserved INFO/FORMAT key tables incomplete (§1.6.1.4, §1.6.2)

Found by the v0.19.0 audit.

**Spec (Reserved INFO keys table, line 460–480):**

```
AA 1 String | AC A Integer | AD R Integer | ADF R Integer
ADR R Integer | AF A Float | AN 1 Integer | BQ 1 Float
CIGAR A String | DB 0 Flag | DP 1 Integer | END 1 Integer
H2 0 Flag | H3 0 Flag | MQ 1 Float | MQ0 1 Integer
NS 1 Integer | SB 4 Integer | SOMATIC 0 Flag
VALIDATED 0 Flag | 1000G 0 Flag
```

**Spec (Reserved genotype keys table, line 534–596):**

```
AD R Integer | ADF R Integer | ADR R Integer | DP 1 Integer
EC A Integer | LEN 1 Integer | FT 1 String | GL G Float
GP G Float | GQ 1 Integer | GT 1 String | HQ 2 Integer
LA . Integer | LAA . Integer | LAD LR Integer | LADF LR Integer
LADR LR Integer | LEC LA Integer | LGL LG Float | LGP LG Float
LPL LG Integer | LPP LG Integer | MQ 1 Integer | PL G Integer
PP G Integer | PQ 1 Integer | PS 1 Integer | PSL P String
PSO P Integer | PSQ P Integer
(+ the M / DPM / ADM base-modification families, already
 covered via MFamily)
```

**Gap:** `vcf/reserved.go` `reservedInfoKeys` carries only
the v4.5-*new* subset (RN, RUS, RUL, RUC, RB, CIRUC, CIRB,
RUB, LAA, LAD) and `reservedFormatKeys` carries only the
local-allele + phase-set + GL/GP/PL families. Every key in
the two spec tables above that is not in those subsets is
absent, so `Validator.checkReservedKeys` silently skips
it. A file declaring, say,
`##INFO=<ID=DP,Number=5,Type=String,...>` or
`##FORMAT=<ID=HQ,Number=G,Type=Float,...>` produces no
diagnostic in either `ModeDefault` or `ModeStrict`.

**Acceptance:**
- All 21 INFO keys and all 30 FORMAT keys listed in the two
  spec tables are present in `reservedInfoKeys` /
  `reservedFormatKeys` with the spec's Number and Type.
- A header declaring any of them with a wrong Number or Type
  emits a `header/INFO/<key>` or `header/FORMAT/<key>`
  diagnostic.
- A header declaring them correctly emits nothing.
- The base-modification families keep resolving through
  `MFamily` rather than literal table entries.

**Where:**
- `vcf/reserved.go`: `reservedInfoKeys` gains all 21 spec-
  table keys (grouped under a "Reserved INFO keys table"
  comment with the spec's line range); `reservedFormatKeys`
  gains all 18 missing non-base-mod spec-table keys.
- `vcf/reserved_test.go`: `TestIsReservedInfoKey` /
  `TestIsReservedFormatKey` — both previously asserted
  that DP/AF/AC/NS and GT/DP were **not** reserved, which
  pinned the gap in place. Rewritten to assert the spec
  tables, with genuinely-arbitrary keys (MN, RM, FOO, ZZZ)
  as the negative cases.
- `vcf/reserved_table_test.go`: new. `specReservedInfo` /
  `specReservedFormat` transcribe the two spec tables
  verbatim; `TestReservedInfoTable_MatchesSpecTable` and
  `TestReservedFormatTable_MatchesSpecTable` diff the
  tables against them, `TestReservedKeys_WrongDeclarationWarns`
  proves 11 wrong declarations are caught.

**Status:** ✅ (v0.20.0)

### T4.26 — Nine reserved FORMAT keys have the wrong Number/Type (§1.6.2)

Found by the v0.19.0 audit. This is a **false-positive**
bug: the validator emits a diagnostic on files that follow
the spec exactly.

The audit originally listed three keys; diffing the whole
spec table against the implementation during the fix found
nine. All nine were wrong:

| Key | Spec (lines 547-555) | Was | Now |
|---|---|---|---|
| LAA  | `.` Integer  | `.` String   | `.` Integer |
| LAD  | `LR` Integer | `.` String   | `LR` Integer |
| LADF | `LR` Integer | `.` String   | `LR` Integer |
| LADR | `LR` Integer | `.` String   | `LR` Integer |
| LEC  | `LA` Integer | `.` Integer  | `LA` Integer |
| LGL  | `LG` Float   | `G` String   | `LG` Float |
| LGP  | `LG` Float   | `G` String   | `LG` Float |
| LPL  | `LG` Integer | `P` String   | `LG` Integer |
| LPP  | `LG` Integer | `G` Integer  | `LG` Integer |

The common thread is that all nine are local-allele keys
whose Number should be one of `LA`/`LR`/`LG`; the original
table used `.` / `G` / `P` instead. LPL and LPP were also
wrong on Type — they are the phred-scaled likelihood and
posterior, so Integer, not Float (Float is LGL/LGP).

Note `laaLenFor` and `checkLAA` already tolerate either
wire form (the LAA values parse as either `Ints` or `Strs`),
so the *parsing* path is fine; only the header-metadata
cross-check is wrong. Related: the T4.17 comment says
"the spec's `LAA Type=String` default", which is
incorrect — v4.5 says Integer.

**Where:**
- `vcf/reserved.go`: all nine entries corrected, each with
  an inline comment naming the spec line it comes from.
- `vcf/reserved_table_test.go`:
  `TestReservedKeys_NoFalsePositiveOnSpecDeclarations`
  walks `specReservedInfo` and `specReservedFormat` and
  asserts a spec-correct declaration produces no
  `header/INFO/<key>` or `header/FORMAT/<key>` diagnostic.

**Acceptance:**
- All nine entries match the spec table (locked by
  `TestReservedFormatTable_MatchesSpecTable`).
- A file with any spec-correct declaration validates clean.

**Status:** ✅ (v0.20.0)

### T4.27 — `checkReservedKeys` doc comment contradicts its gating (§1.6.1.4)

Found by the v0.19.0 audit.

`Validator.checkReservedKeys` carries this doc comment:

> The check is non-fatal: the spec lists these as
> reservations, not strict requirements, so we surface a
> warning the caller can promote to an error by switching
> to ModeStrict.

But `HeaderChecks` gates the entire call on
`if v.Mode == ModeStrict` (vcf/validate.go:234), so in
`ModeDefault` there is *no* reserved-key check at all —
there is nothing to "promote". The comment describes a
default-mode check that does not exist.

Shipped the first option: the gate moved off the
diagnostic. `checkReservedKeys` now runs in **both** modes
and still reports `SeverityWarning`. A reserved key
declared with the wrong Number or Type is a header defect
regardless of mode, and the spec frames these keys as
reserved "albeit optional" — visible, not fatal. This is
also the larger user win: the check is now reachable
without asking callers to opt into strict mode.

**Where:**
- `vcf/validate.go`: `HeaderChecks` drops the
  `if v.Mode == ModeStrict` guard and preallocates `issues`
  to `len(h.Infos)+len(h.Formats)` (the bound the two
  contributing checks can hit).
- `vcf/reserved_table_test.go`:
  `TestReservedKeys_RunInDefaultMode` asserts a wrong
  `##INFO=<ID=DP,Number=5,Type=String,...>` produces a
  `SeverityWarning` under `ModeDefault`.

**Acceptance:**
- Doc comment matches the implemented gating.
- The check is reachable in `ModeDefault`.

**Status:** ✅ (v0.20.0)

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
| §1.6.1 Reserved INFO keys table | T4.25 |
| §1.6.1.7 SV info keys | T2.3, T2.4, T2.5, T2.9, T2.10 |
| §1.6.1.7 CNV section | T2.6 |
| §1.6.1.7 CNV:TR section | T2.7 |
| §1.6.2 FORMAT column | T2.1, T2.2, T2.8, T4.25, T4.26, T4.27 |
| §1.6.2 PSO field | T3.3 |
| §1.7 gVCF | T1.6, T4.9 |
| §1.6.5.1 Breakend round-trip | T4.10 |
| §1.7 Tandem repeats | T2.9, T2.10 |
| General file format | T4.1 |
| §1.4.4 IUPAC | T4.2 |
---

## Phase 5 — conformance verified against the spec's own worked examples

The v0.19.x-v0.21.0 audits were a **static read of the code
against the spec prose**. That method demonstrably misses
whole classes of defect: on the same pass it found a
catastrophic BAM CIGAR bit-layout inversion, a Tabix reader
missing 8 header fields, and a CSI reader missing 3 - none
of which the prose-level walk surfaced, and all of which
the existing round-trip tests passed.

The independent check is the spec's own **worked example**
(§1.6.1, lines 43-68): a complete 24-line VCF with header,
4 INFO/FILTER/FORMAT declaration groups and 5 records
covering phased/unphased GT, a Flags INFO pair, `Number=A`
Float cardinality, a multi-value missing `.,.`, an ALT=`.`
record, and a 2-allele indel with differing FORMAT widths.

Run against `htsfile check --strict`, `roundtrip`, and
`vcf2bcf`/`bcf2vcf`, the example surfaces the following.

### T5.1 - `.,.` list-of-missing-values corrupted to `0,0`

**Spec (§1.6.2, line 508):**

> If a field contains a list of missing values, it can be
> represented either as a single MISSING value (`.') or
> as a list of missing values (e.g. `.,.,.' if the field
> was Number=3).

**Defect:** a `Number=2` FORMAT field written `.,.` reads
back and re-emits as `0,0`. Minimal repro:

```
in:  20  14370  .  G  A  29  PASS  .  GT:HQ  0|0:.,.
out: 20  14370  .  G  A  29  PASS  .  GT:HQ  0|0:0,0
```

Not sample-count dependent - reproduces on a single sample.
The spec's own example record carries `1/1:43:5:.,.`.

**Fix (v0.22.0):** added per-slot `IntMissing []bool` /
`FloatMissing []bool` to `typedCore` in `vcf/value.go`.
`vcf.decodeTyped` records the missing marker for every
"."-parsed slot; `vcf.writeTyped` re-emits "." for those
slots even when `v.Missing` is unset. Mirrored into the
BCF layer (`bcf.SampleValue.IntMissing`, `FloatMissing`)
so per-slot round-trips work both ways.

**Status:** shipped in v0.22.0

### T5.2 - BCF round-trip drops FORMAT and all sample columns

**Defect:** `vcf2bcf` then `bcf2vcf` on the spec example
returns a header with 8 tab-separated fields instead of 12:
the FORMAT column header and all three sample IDs are
gone, and every record loses its genotype data.

**Fix (v0.22.0):** the VCF->BCF translator now carries
`Format []uint16` and `Samples []bcf.Sample` (vs. the
old "FORMAT dropped on vcf2bcf" scope-cut). Symmetrically
`bcfToVCFRecord` rebuilds `Format []string` and `Samples
[]vcf.Genotype` keyed by the FORMAT declaration order
(preserved across the BCF wire via `FormatsOrder`).
FORMAT keys present on a record but absent from the
header are auto-registered as `Number=1,Type=Integer`
(a conservative default matching GT/DP/GQ/HQ) so
incomplete VCFs round-trip without dropping fields.
Trailing wholly-missing FORMAT fields are also dropped
on output per §1.6.1 "trailing fields can be dropped".

**Status:** shipped in v0.22.0

### T5.3 - BCF writer output is rejected by our own BCF reader

**Defect:** `htsfile check` on a `.bcf` produced by
`htsfile vcf2bcf` fails with

```
parse record: line 0 INFO: Float field has type byte 4
```

The typed-INFO encoding written by `bcf.writeRecord` is
not what `bcf.parseAux` expects - a writer/reader
disagreement inside one package, independent of any
external tool.

**Fix (v0.22.0):** the BCF record layer was completely
mis-implemented against the spec's `BcfSiteEncoding` and
`BcfTypeEncoding`. The actual spec layout - every record
prefixed by `l_shared` + `l_indiv` uint32 length prefixes,
then CHROM/POS/rlen/QUAL/n_info/n_allele/n_sample/n_fmt,
then a typed ID+REF+ALT+FILTER+INFO block using the
`(size<<4)|type` byte with size=15 overflow followed by
a typed Integer count - was rebuilt in `bcf/encodeRecord`,
`bcf.Encoder.WriteRecord`, `bcf.Decoder.Record`,
`bcf.writeTypedInt` / `bcf.readTypedInt`. INFO keys are
typed ints (not bare uint16) per the spec, and the
variable-length type byte plus overflow count are now
emitted/parsed consistently.

**Status:** shipped in v0.22.0

### T5.4 - BCF header Description strings corrupted

**Defect:** `Description="Allele Frequency"` becomes
`Description="\""` after a BCF round-trip. Descriptions
containing a comma (`"dbSNP membership, build 129"`) are
additionally truncated at the comma to
`Description="dbSNP membership`.

**Fix (v0.22.0):** two bugs in `bcf/header.go`:

1. `parseStructuredField` split structured-line fields
   with `strings.Split(rest, ",")` - a naive split that
   didn't respect double-quoted sub-spans. So
   `Description="dbSNP membership, build 129"` was cut at
   the comma. Replaced with a quote-aware
   `splitStructuredFields` that toggles a `inQuote` flag
   on `"` and only splits on unquoted commas.
2. `extractIDFromAngleField` (used for ##contig and
   ##FILTER) had the same naive split. Also routed
   through `splitStructuredFields`.

The BCF->VCF translator also rehydrates FILTER and
contig Description metadata from the BCF `hdr.Raw` text
(via `extractIDAndDescription` and `populateContig` in
`cmd/htsfile/translate.go`) so a full VCF<->BCF round-trip
preserves Descriptions and `##contig=length,md5,species,...`.

**Status:** shipped in v0.22.0

### Verified clean against the spec example

- `htsfile check --strict` parses and validates the §1.6.1
  example with **zero** diagnostics.
- Single-value and multi-value `Number=2` FORMAT fields
  round-trip intact when populated (`0|0:51,51` is
  preserved).
- Phased (`0|0`, `1|2`), unphased (`1/1`), and
  mixed-phase GT all round-trip intact.
- INFO Flag fields (`DB`, `H2`) round-trip intact.
- `Number=A` Float cardinality (`AF=0.333,0.667` against
  2 ALTs) validates and round-trips.
- ALT=`.` (no-variant) record round-trips intact.
- `md5` on `##contig` gains quotes on output. This is
  **more** spec-compliant than the spec's own example:
  §1.4.3 line 126 requires optional structured-line field
  values to be quoted. Not a defect.
- INFO key ordering is not preserved. **Spec-allowed:**
  §1.4.3 line 133 - "Implementations must not rely on the
  order of the fields within structured lines and are not
  required to preserve field ordering."

## Phase 6 — conformance verified against an external oracle

Phases 4 and 5 were verified against the spec's own worked
examples and against our own reader/writer pair. That is not
enough: a writer and a reader that agree with each other and
with nobody else pass every round-trip test in the suite.

So the reference implementation was installed as an oracle
(htslib 1.24.1 / samtools 1.24 / bcftools 1.24, via Homebrew)
and our output was fed to it. **Every defect below was
invisible to the 600+ green tests in this repository** and
was found in minutes by a tool that does not share a line of
code with us.

The v0.22.0 BCF output was unreadable by *any* modern htslib
tool:

```
[E::bcf_hdr_read] Invalid BCF2 magic string: only BCFv2.2 is supported
[W::bcf_sr_add_hreader] No BGZF EOF marker; file 'x.bcf' may be truncated
[W::bcf_record_check] Bad BCF record at 20:14370: Invalid REF/ALT type 9 (unknown)
[E::bcf_record_check] Bad BCF record at 20:14370 - shared section malformed or too short
[W::bcf_record_check] Bad BCF record at 20:17330: Invalid FORMAT id 12
[E::bcf_format_gt_v2] Unexpected type 7
```

A mutation experiment confirmed the blind spot. Swapping
POS and rlen in the writer *and* the reader together — a
spec-breaking change that both sides agree on — was caught by
**1 of 633 test functions**.

### The BCF defects

### T6.1 - BCF magic was 6 bytes with `l_text` at offset 6

**Spec (§BCF2 header, tex lines 2086-2091):**

> magic char[3] — File format identifier, `BCF`.
> major_version uint8 — Version of the specification.
> minor_version uint8 — Version of the specification.
> l_text uint32_t — Length of the 'text' field.

**Defect.** We wrote `BCF\x02\x01\x00` (6 bytes) and put
`l_text` at offset 6. The magic is 5 bytes: `BCF` plus the
major and minor version. BCF 2.2 is `BCF\x02\x02`.

**Fix.** `Magic = [5]byte{'B','C','F',0x02,0x02}`; `l_text`
at offset 5. The 2.1 magic is kept for *reading* only.

**Status:** shipped in v0.23.0

### T6.2 - BGZF EOF marker had the wrong OS byte

**Defect.** htslib's `bgzf_check_EOF` is a byte-exact
`memcmp` against the 28-byte EOF block. We wrote `0x00` in
the OS-byte position; the marker requires `0xff`.

**Fix.** The EOF block is now byte-identical to htslib's.

**Status:** shipped in v0.23.0

### T6.3 - `n_sample` / `n_fmt` were packed in the wrong byte order

**Spec (§BcfSiteEncoding, tex lines 2191-2196):**

> n_sample uint24 - Number of samples in this record.
> n_fmt uint8 - Number of FORMAT fields.

**Defect.** We packed `n_sample<<8 | n_fmt`. The spec puts
`n_sample` in bytes [20:22] as a uint24 and `n_fmt` in byte
[23].

**Fix.** `n_sample | n_fmt<<24`. The reader also cross-checks
`n_sample` against the header's sample count.

**Status:** shipped in v0.23.0

### T6.4 - strings were not char vectors

**Spec (§BcfTypeEncoding, tex line 2276):** a string is a
character vector whose length is carried in the descriptor's
high nibble.

**Defect.** We wrote type `0x07` followed by a `uint32`
length and then the bytes. htslib read the length as the
vector size and desynchronised the stream.

**Fix.** `writeBCFString` emits the descriptor with the
character count in the high nibble and the raw characters
after it.

**Status:** shipped in v0.23.0

### T6.5 - FILTER wrote a bare type code

**Defect.** FILTER wrote `0x01` (a bare type code) instead of
a `(size<<4)|type` descriptor, which declares a vector of
size 0.

**Fix.** `pickBcfSize` + `typeByteSize`; `readTypedFilter`
reads the descriptor.

**Status:** shipped in v0.23.0

### T6.6 - multi-element String INFO lost its commas

**Defect.** We wrote the elements as separate strings.
htslib keeps the commas inside the character run:
`CSQ=alpha,beta,gamma` becomes one 16-byte run.

**Fix.** The writer comma-joins; the reader splits on commas.

**Status:** shipped in v0.23.0

### T6.7 - FORMAT String size was the element count, not the character count

**Defect.** A FORMAT String field's size nibble must carry
the character count, not the number of elements.

**Fix.** `fieldMaxSize` counts the joined characters.

**Status:** shipped in v0.23.0

### T6.8 - the FORMAT overflow length was not a typed integer

**Defect.** For a vector of 15 or more elements the
descriptor's nibble is 15 and the true length follows as a
typed Integer. We wrote a bare byte.

**Fix.** The overflow length is a typed `int8`.

**Status:** shipped in v0.23.0

### T6.9 - the per-sample block was not the spec's triplet layout

**Spec (§GenotypeEncoding, tex line 2233):**

> Suppose there are i genotype fields in a specific record.
> Each i is encoded by a triplet

**Defect.** We wrote every `fmt_key` first and then every
`(fmt_type, fmt_values)` pair. htslib rejected the record
with *"individuals section malformed or too short"*.

**Fix.** `encodeIndivid` writes `fmt_key`, `fmt_type`,
`fmt_values` per field, interleaved.

**Status:** shipped in v0.23.0

### T6.10 - GT was written as a string

**Defect.** The VCF header declares `GT` as
`Number=1,Type=String`, and we wrote it that way. BCF
requires an **integer allele vector**; htslib's
`bcf_format_gt_v2` rejects any other type
(*"Unexpected type 7"*).

**Fix.** `bcf/gt_wire.go` transcribes htslib's encoding
(htslib 1.24 `vcf.c`, `vcf_parse_format`'s GT branch): each
allele is `(allele+1)<<1 | phased`, a missing allele is 0,
and the phased bit of an allele describes the separator
*before* it. The first allele's bit is the "prefixed
phasing" flag. Verified byte-for-byte against
`bcftools view -Ob` for `0|1`, `0/1`, `1|0`, `1/0`, `./.`,
`0`, `0|1/2` and `0/1|2`.

**Status:** shipped in v0.23.0

### T6.11 - keys used positional indices instead of the header dictionary index

**Spec (§"Dictionary of strings", tex line 2124):**

> ASP=1, RSPOS=2, dbSNPBuildID=3, and PASS implicitly
> encoded in the first offset PASS=0

**Defect.** We wrote the position of a key in `InfosOrder` /
`FormatsOrder` and the position of a filter in `Filters`.
htslib's `BCF_DT_ID` is a *single* dictionary shared by
FILTER, INFO and FORMAT, so `DP` used as both an INFO and a
FORMAT key has one index.

**Fix.** `Header.IDXFor` / `Header.IDAt` implement the
dictionary: 1-based, per *distinct* ID (first occurrence
wins), in header-text order; `##contig` lines do not consume
an index; an explicit `##FILTER=<ID=PASS>` line does not
either (PASS is implicitly 0). Verified against
`bcftools view -Ob` for the spec example and for synthetic
headers with 0, 1 and 2 contigs.

**Status:** shipped in v0.23.0

### T6.12 - the header text was not NUL-terminated

**Spec (§BCF2 header, tex lines 2086-2091):** `l_text` is the
length of the text field *including* the terminating NUL.

**Defect.** We wrote `l_text = len(text)` and no terminator.

**Fix.** The writer emits `len(text)+1` and a NUL; the
reader strips it.

**Status:** shipped in v0.23.0

### T6.13 - a missing ID was written as typeMissing

**Spec (§BcfSiteEncoding, tex line 2196):**

> ID - typed string - Variant identifier; 0x07 for a missing
> value

**Defect.** We wrote `0x00`, which htslib rejects with
*"Invalid ID type 0 (null)"*.

**Fix.** A missing ID is `0x07` (also the descriptor for a
zero-length character vector).

**Status:** shipped in v0.23.0

### T6.14 - integers were always int32

**Spec (§BcfTypeEncoding, tex line 2266):** *"It is up to the
encoder to determine the appropriate ranged value to use"* -
8, 16 or 32 bit values are all legal.

**Defect.** Not a defect, but htslib picks the narrowest width
that fits and our int32-only output was 4x larger.

**Fix.** `intWidthFor` picks int8, then int16, then int32.
Our output is now byte-identical to htslib's.

**Status:** shipped in v0.23.0

### T6.15 - a present Flag wrote `0x01`

**Defect.** htslib writes a single `0x00` byte for a present
Flag; the INFO entry's presence carries the information.

**Fix.** Writer and reader both follow htslib.

**Status:** shipped in v0.23.0

### T6.16 - an absent FORMAT value padded with MISSING instead of END_OF_VECTOR

**Defect.** For a fixed-size field with no value for a sample
we wrote a MISSING sentinel in every slot. htslib writes one
MISSING element and then pads with END_OF_VECTOR.

**Fix.** Matched. An explicit `.,.` still writes one MISSING
element per slot, which is a different thing.

**Status:** shipped in v0.23.0

### The SAM/BAM defects

The same oracle sweep found three defects in the BAM path.

### T6.17 - BAM wrote ASCII quality bytes instead of raw Phred values

**Defect.** BAM stores the base quality as a raw Phred score
(0-93). We wrote the ASCII `+33` form, so a Phred 40 became
the byte `73` on the wire; samtools read it as Phred 73 and
printed `j` instead of `I`.

**Fix.** `encodeQual` converts on write; the decoder converts
back. `0xFF` (an unavailable quality) is preserved and
renders as `*`.

**Status:** shipped in v0.23.0

### T6.18 - `RNEXT "="` was rejected

**Spec (SAMv1 §1.4.1):** `RNEXT` may be `=`, meaning the mate
is on the same reference as the read. htslib
(`sam.c`, `bam_parse1`) resolves it as `mtid = tid`.

**Defect.** We looked `=` up in the `@SQ` list and failed with
*"unknown reference "=" (not in header @SQ list)"*.

**Fix.** The SAM-to-BAM translator resolves `=` to the
record's own `RNAME`. The BAM-to-SAM translator writes `=`
when the mate is on the same reference, as htslib does.

**Status:** shipped in v0.23.0

### T6.19 - the BAM-to-SAM header was re-synthesised

**Defect.** Converting a BAM back to SAM rebuilt the header
lines from the parsed structures, dropping every tag the
package does not model (`M5`, `UR`, `SP` on `@SQ`; `PL`,
`PU`, `CN` on `@RG`; `VN`, `CL` on `@PG`).

**Fix.** The BAM header stores the SAM header text verbatim,
so the original text is now carried through the conversion
(`sam.Header.Raw`).

**Status:** shipped in v0.23.0

### T6.20 - the BAM header text lacked a trailing newline

**Defect.** htslib warned *"Missing trailing newline on SAM
header. Possibly truncated"*.

**Fix.** The header text ends with a newline, as htslib's
does.

**Status:** shipped in v0.23.0

### Verified against the oracle

After the fixes above:

- `bcftools view`, `bcftools query`, `bcftools index`,
  `bcftools sort` and `htsfile` all read our BCF with
  **zero** output on stderr.
- The BCF records we write are **byte-identical** to
  `bcftools view -Ob` output for the spec's own worked
  example - all 5 records, both the shared and the
  per-sample block.
- `samtools view`, `samtools quickcheck` and `samtools
  idxstats` all read our BAM with zero output on stderr, and
  the records match `samtools view` of a samtools-produced
  BAM exactly.
- `sam -> bam -> sam` is byte-identical to the input SAM.

### The second oracle pass

The v0.23.0 release was verified against one input file: the
spec's own worked example. Asking "would you consider this ready
for release?" prompted the obvious question that had not been
asked - **can we read htslib's files?** - and a second pass that
built a fixture covering the types the spec example never
reaches.

That pass found **8 more defects**, including one that silently
corrupted data:

> `DBL=.,.` was written to the BCF as `DBL=0,0`.

Every test in the suite passed while that happened, because the
writer and the reader agreed with each other: the reader read
the zero back and the round trip looked fine.

### T6.21 - a missing INFO value was written as zero

**Spec (§BcfTypeEncoding, tex line 2304):** the values
`0x80000000`, `0x8000` and `0x80` are the per-type MISSING
sentinels for 32-, 16- and 8-bit integers.

**Defect.** A missing INFO element was written as the number 0.
htslib writes the per-type sentinel. `DBL=.,.` became
`DBL=0,0`; `VR=` became `VR=0`.

**Fix.** `bcf.InfoValue` gained `IntMissing` / `FloatMissing`
parallel slices (mirroring the `vcf` side added in v0.22.0), the
writer emits the sentinel, and the reader recognises it. An
INFO key with an empty value list is one missing element, which
is what htslib writes.

**Status:** shipped in v0.23.1

### T6.22 - a missing or empty INFO String was written as typeMissing

**Defect.** An empty String INFO was written as `0x00`
(typeMissing). htslib writes `0x07` - a zero-length char
vector, the same form the spec gives for a missing ID.

**Fix.** The empty char vector.

**Status:** shipped in v0.23.1

### T6.23 - an empty Character INFO was written as typeMissing

**Defect.** Same as T6.22, for the Character type.

**Fix.** The empty char vector.

**Status:** shipped in v0.23.1

### T6.24 - a FORMAT field collapsed to the whole-field typeMissing form

**Defect.** When every sample of a FORMAT field was missing, we
wrote a single `0x00` for the whole field. htslib always writes
the field with a per-sample value: the per-type MISSING
sentinel for numbers, and the literal `.` text for strings and
characters.

**Fix.** The per-sample model.

**Status:** shipped in v0.23.1

### T6.25 - a FORMAT String dropped the commas between its elements

**Defect.** A multi-element FORMAT String was written as
`aabb` - the element boundaries were lost, so the value could
not be read back. htslib copies the field text verbatim,
commas included (`aa,bb`).

**Fix.** The raw field text, NUL padded to the field size.

**Status:** shipped in v0.23.1

### T6.26 - a FORMAT Character was one byte per sample

**Defect.** htslib stores a FORMAT Character as the character
plus a NUL, so two bytes per sample and a size nibble of 2. We
wrote one byte.

**Fix.** The character plus NUL padding.

**Status:** shipped in v0.23.1

### T6.27 - a BAM integer tag became an invalid SAM type

**Defect.** BAM stores an integer aux value with a
width-specific subtype byte (c, C, s, S, i, I). Converting a BAM
back to SAM, we emitted the subtype letter - `NM:C:1`. SAM has no
`C` type: the valid set is A, i, f, Z, H, B. htslib reads
`NM:C:1` as the *character* `1`, so the value silently changed
from the integer 1 to the character '1'.

**Fix.** All six integer subtypes map to the SAM letter `i`,
which is what htslib writes.

**Status:** shipped in v0.23.1

### T6.28 - a VCF with samples but no `##FORMAT` declarations was rejected

**Defect.** Three separate problems: the BCF header's `#CHROM`
line omitted the FORMAT column when the header declared no
FORMAT keys, so the sample names landed in the wrong column and
the BCF header declared zero samples; the synthesised
`##FORMAT` lines were appended *after* the `#CHROM` line, which
htslib rejects; and an undeclared FORMAT key was registered as
`Number=1,Type=Integer`, which silently dropped a `GT` value
like `0|0`.

htslib accepts such a file - a VCF with samples must declare
GT, but htslib tolerates its absence.

**Fix.** The FORMAT column is always written when there are
samples; the synthesised header lines go before the `#CHROM`
line; an undeclared key is registered as a String (htslib's
rule), and a reserved key keeps the spec's own Number and Type.

**Status:** shipped in v0.23.1

### The gVCF END rule

The spec's own gVCF example (§1.7, tex lines 1878-1886)
triggered five diagnostics from `check --strict`. Two of them
are typos in the spec's example - a semicolon where the FORMAT
separator should be a colon in three rows, and a missing `PL`
in one FORMAT column - and the rest were our END rule.

The rule is that a gVCF reference block's END is `POS+LEN-1`.
The spec's example uses `POS+LEN` in three of its seven rows and
`POS+LEN-1` in the rest, so the validator now accepts either
form for a `<*>` block (`vcf/validate_record.go`, `checkEND`).

### The guard

`bcf/spec_wire_test.go` holds the reference bytes as literal
hex strings and asserts that our encoder emits exactly them.
A fixture built by the code under test cannot catch a change
to it; these bytes come from `bcftools view -Ob`.

The guard is mutation-tested: each of the six defects above
that changes the record layout (T6.3, T6.4, T6.9, T6.10 and
the two halves of T6.11) was reintroduced and the guard
failed for every one of them.

`cmd/htsfile/type_coverage_test.go` pins the type-coverage
fixture to the bytes `bcftools view -Ob` produces for it, so
T6.21 - T6.26 cannot regress.

`vcf/testdata/type_coverage.vcf` exercises String (single and
multi-element), Character, Float, Number=R / G / . variable
cardinality, a String longer than the 15-byte descriptor
nibble, and a missing value in every position.

## Phase 7 — the tabix / CSI index, verified against the oracle

The tabix and CSI index code was written by reading the
format descriptions, not by consulting an implementation, and
it had never been run against htslib's. Doing so found **nine
defects**, all of them invisible to the suite.

### T7.1 - the readers could not open htslib's index at all

htslib writes `.tbi` and `.csi` files BGZF-compressed. Our
readers read the raw bytes, so every index htslib produced
was rejected:

```
tabix: not a Tabix file: got "\x1f\x8b\b\x04"
```

**Fix.** `internal/bgzf.InflateIfBGZF` sniffs the gzip magic
and inflates the index transparently; a raw index is still
accepted (`tabix/tabix.go`, `csi/csi.go`).

### T7.2 - the writers wrote an uncompressed index

htslib rejects an uncompressed `.csi` outright:

```
[E::index_load] Invalid index header for x.vcf.gz
```

**Fix.** Both writers now frame the index with BGZF, the way
htslib does (`tabix/write.go`, `csi/write.go`).

### T7.3 - the `.tbi` had no reference names

The index carried no name table, so a region query could not
resolve a reference and htslib's `tabix` fell back to
scanning the file.

**Fix.** The names come from the header's contigs and go into
the `.tbi`'s name table and the `.csi`'s aux field
(`vcf/index_writer.go`).

### T7.4 - the `.tbi` had no VCF preset

`format`, `col_seq`, `col_beg`, `col_end` and `meta` were all
zero. htslib's preset is `format=2` (VCF), `col_seq=1`,
`col_beg=2`, `col_end=0`, `meta='#'`.

**Fix.** The writer sets the preset (`vcf/index_writer.go`).

### T7.5 - the `.csi` had the wrong depth and no aux

A CSI records the number of bin *levels*, which is one more
than the tabix depth: htslib writes 6 for the tabix depth of
5. Ours wrote 5, and the aux field - which carries the tabix
preset and the reference names - was empty.

**Fix.** `Depth: b.depth + 1`, and the aux field holds the
tabix preset followed by the name table, byte for byte as
htslib writes it.

### T7.6 - a CSI was written with the tabix bin numbering

The bin numbers depend on the number of levels: a tabix index
files a record under bin 4681+, a CSI under 37449+. We used
the tabix numbering for both, so every CSI bin was in the
wrong place.

**Fix.** The builder picks the leaf bin with the format's own
level count (`csi.BinFromRangeWith`).

### T7.7 - a CSI carried a linear index it should not have

A CSI has no linear index: the per-bin `loffset` replaces it.
htslib frees the linear index once the loffsets are computed.
Our `.csi` wrote one anyway, so htslib's reader mis-parsed
the file from that point on - which is why its region queries
came back wrong.

**Fix.** The CSI writer emits no linear index, and writes
`n_no_coor` once at the end of the file rather than per
reference. The reader was changed to match.

### T7.8 - the chunks pointed at the wrong offsets

Two bugs stacked here.

The BGZF reader reported a stale virtual offset once a block
was exhausted: the offset stayed at the end of the block just
read instead of moving to the start of the next one. The
reader also read blocks through a `bufio.Reader`, so the
counting reader's position no longer marked the block
boundary, and every offset after the first block was wrong.

**Fix.** The blocks are read straight from the counting
reader, and the next block is loaded as soon as the current
one is drained, so the virtual offset of the next byte to
return is always right (`internal/bgzf/bgzf.go`).

The second bug was in the index builder: it took the record's
position from the VCF decoder's reader, but the decoder reads
lines through a `bufio.Scanner` that reads ahead, so the
position was past the start of the next record.

**Fix.** `Decoder.SetExactReads` makes the decoder read lines
one byte at a time, so the reader's position stays exactly at
the end of the line just read. The index builder and the
chunk decoder use it.

### T7.9 - the chunking did not follow htslib's convention

htslib's chunk runs from the start of the first record in the
bin to the start of the first record in the next bin, and the
linear index entry for a window is the start of the first
record at or after that window, with an empty window
inheriting the entry of the window after it. Ours used
block-granular chunks and a linear index that was mostly
zeroes, so a query that started mid-file scanned from the
beginning of the file - correct, but far too slow, and not
the index htslib writes.

**Fix.** The builder follows htslib's rules. Our `.tbi` is now
byte-identical to htslib's apart from the order the bins are
written in, which is an artefact of htslib's hash table and
carries no meaning.

### The guard

`vcf/index_oracle_test.go` builds an index over the spec
example with our code and with htslib's `tabix`, and asserts
that the two hold the same bins, the same chunks, the same
linear index, the same per-bin loffsets and the same aux
bytes - and that htslib's own `tabix` answers the same region
queries with either index. The test is skipped when htslib's
tools are not installed.

### Status

Shipped in v0.24.0.

### The remaining items

Three items were left over from the audit of v0.23.1. None of
them turned out to be a defect in the code, except one small
one; what they needed was a decision, a fixture and a guard.

#### The spec's gVCF example

The spec's own gVCF example (VCFv4.5 §1.7, tex lines
1862-1874) has a typo: five of its seven rows separate the
last two FORMAT values with a semicolon where a colon belongs,
so `PL;LEN` arrives as the single value `900;14`, which is not
an Integer. The LEN value never arrives, and the END check
fails with it, because a reference block's END is POS+LEN-1.

The diagnostics are right and the example is wrong. The
corrected example is shipped as
`vcf/testdata/gvcf_example.vcf` and must validate clean; the
example exactly as printed must produce exactly the nine
issues its typos cause (five FORMAT values that are not
Integers, and the four END checks that fail because the LEN
value never arrived). Both are pinned by
`vcf/gvcf_example_test.go`.

#### A BCF 2.1 file

htslib 1.24 writes BCF 2.2 only and rejects a 2.1 file
outright, and no 2.1 file can be produced with the tools
available here. The spec's list of changes between BCFv2.1 and
BCFv2.2 (tex lines 2818-2825) is about the header's optional
IDX field, the END_OF_VECTOR byte and where PASS sits in the
dictionary - none of which changes the record layout - so a 2.1
file is a 2.2 file with a different minor-version byte.

We read one: `bcf/bcf21_test.go` patches the minor-version byte
of a file it builds and reads it back. That test found a small
defect: the decoder reported the *major* version (always 2)
rather than the minor one, because it read the wrong byte of the
magic. The spec fixes magic[3] = major_version and magic[4] =
minor_version, and htslib writes the literal `BCF\x02\x02`, so
the two bytes are indistinguishable in a 2.2 file - which is how
the mistake survived. Fixed in `bcf/parse.go`.

#### The `##dictionary=` pragma

The spec's BCF2 header section mentions a string-dictionary
pragma (tex line 2100):

> BCF2 supports encoding values in a dictionary of strings. The
> string map is provided by the keyword
> `##dictionary=S0,S1,...,SN` as a comma-separate ordered list of
> strings.

The "Dictionary of strings" section that follows (tex lines
2103-2143) then describes the *header* dictionary: the position
of each FILTER / INFO / FORMAT key in the header, which is what
a BCF2 record's typed integers refer to. The `s0`, `s1` in that
section's example are a visual guide - the spec says so
explicitly - not the wire format. That encoding is the one we
implement (T6.11).

So in VCFv4.5 the pragma has no wire encoding of its own: the
header dictionary does that job. A header line carrying it is
still valid, so it must be accepted and preserved rather than
rejected. It is: the VCF parser takes it, the encoder replays
it, and it survives a round trip through BCF (htslib's
`bcftools` keeps it too). Pinned by
`vcf/dictionary_test.go`.
