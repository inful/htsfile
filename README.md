# htsfile

[![CI](https://github.com/inful/htsfile/actions/workflows/ci.yml/badge.svg)](https://github.com/inful/htsfile/actions/workflows/ci.yml)
[![Go Report Card](https://goreportcard.com/badge/github.com/inful/htsfile)](https://goreportcard.com/report/github.com/inful/htsfile)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`htsfile` is a Go library for parsing, writing, and validating bioinformatics
file formats. It currently targets nine specifications from
[samtools/hts-specs](https://github.com/samtools/hts-specs) and
the Sequence Ontology:

- **VCF v4.5** — Variant Call Format (accepts any v4.x minor)
- **BED v1** — Browser Extensible Data (UCSC-derived, BED3 through BED12)
- **SAM v1** — Sequence Alignment/Map (header + 11-column records + optional tags)
- **BAM v1** — Binary SAM (BGZF-compressed) with read-only BAI index reader
- **GFF3** — Generic Feature Format Version 3 (9 columns + free-form attributes)
- **FASTA** — Sequence records (`>` header + multi-line sequence); read-only `.fai` index reader
- **FASTQ** — Short-read sequencing records (`@` header + sequence + `+` + qualities), with Phred+33 / Phred+64 auto-detection
- **BCF v2.1** — Binary VCF (BGZF-compressed). **Read-only** in v1.
- **Tabix v1** — Generic position index for sorted bgzipped files. **Read-only** in v1.
- **CSI v1** — Coordinate-sorted index (successor to Tabix). **Read-only** in v1.

See [ROADMAP.md](ROADMAP.md) for the planned expansion (CRAM v3,
htsget/refget, crypt4gh, etc.).

The implementation aims to be:

- **Spec-compliant** — covers the v4.5 VCF grammar for core records, INFO/FORMAT
  metadata, symbolic structural variants, and breakend notation; and the BED v1
  grammar for BED3 through BED12 with auto-detected column counts, track lines,
  and block invariants. Base modifications (`M*` fields), tandem repeats
  (`RUS/RUL/RUC/RB/...`), and local allele encoding (`LAA/LAD/LPL/LGL/...`) are
  deferred to a follow-up release for VCF. BEDGraph, narrowPeak, broadPeak, and
  GFF/GTF are out of scope for v1; BED-derivatives with column-count extension
  may be added later.
- **Streaming** — `Decoder` reads headers and records incrementally from any
  `io.Reader`. No buffering of the whole file.
- **Round-trippable** — reading a VCF or BED file and writing it back produces
  output that re-parses into the same logical content. Track lines and structured
  meta-information fields are preserved verbatim when possible.
- **Validating** — `Validator` runs header-declared cardinality (`Number=1/A/R/G/...`)
  and type (`Integer/Float/Flag/Character/String`) checks against every record
  and surfaces every issue it finds, not just the first. For BED, the validator
  additionally enforces column-count consistency, coordinate sanity, the score
  range `[0, 1000]`, strand value, RGB range, and BED12 block invariants.

## Install

```sh
go get github.com/inful/htsfile
```

## Usage

### Parsing

```go
import "github.com/inful/htsfile/vcf"

f, err := os.Open("sample.vcf")
if err != nil { log.Fatal(err) }
defer f.Close()

dec := vcf.NewDecoder(f)
hdr, err := dec.Header()
if err != nil { log.Fatal(err) }

for {
    rec, err := dec.Record()
    if err == io.EOF { break }
    if err != nil { log.Fatal(err) }

    fmt.Println(rec.Chrom, rec.Pos, rec.Ref, rec.Alt, rec.Info["DP"])
}
```

### Writing

```go
enc := vcf.NewEncoder(os.Stdout)
if err := enc.EncodeHeader(hdr); err != nil { log.Fatal(err) }
for _, rec := range records {
    if err := enc.EncodeRecord(rec); err != nil { log.Fatal(err) }
}
```

### Validating

```go
v := &vcf.Validator{Mode: vcf.Default}   // or vcf.Strict
issues := v.Validate(hdr, records)
for _, iss := range issues {
    fmt.Printf("line %d: [%s] %s\n", iss.Line, iss.Severity, iss.Message)
}
```

### BED Parsing

The `bed` package mirrors the VCF streaming API:

```go
dec := bed.NewDecoder(r)
track, _ := dec.Track()    // optional; nil if absent
for {
    rec, err := dec.Record()
    if errors.Is(err, io.EOF) { break }
    // rec.Chrom, rec.ChromStart, rec.ChromEnd, ...
    // Optional fields (Name, Score, Strand, ItemRgb, ...) are pointers.
}
```

The decoder auto-detects BED3 through BED12 by column count. Track
lines are parsed and round-tripped; browser lines are skipped. See
`bed.ExampleNewDecoder` for a full runnable example.

### BAM Parsing

The `bam` package reads BGZF-compressed BAM files:

```go
dec := bam.NewDecoder(f)
hdr, _ := dec.Header()    // required; BGZF + binary framing
for {
    rec, err := dec.Record()
    if errors.Is(err, io.EOF) { break }
    // rec.QName, rec.Pos, rec.Flag, rec.MapQ, rec.Cigar, rec.Seq,
    // rec.Qual, rec.Tags, ...
    // rec.CigarString() returns the SAM text form ("5M2I3M" etc.).
}
```

The decoder parses all 11 BAM core fields plus the variable-length
remainder (read name, packed SEQ, qual, optional tags). SEQ is decoded
from 4-bit packed form back to ASCII; tags are parsed for all
documented BAM type codes (A, c, C, s, S, i, I, f, d, Z, H, B) with
values exposed as SAM-compatible strings. BAI files are read by the
sibling `bam/bai` package; see `bai.ExampleNewReader` if you want
random access. BAM v1 is read-only in this build.

### GFF3 Parsing

The `gff` package reads and writes GFF3 (Generic Feature Format
Version 3):

```go
dec := gff.NewDecoder(r)
hdr, _ := dec.Header()    // required; returns Version + pragmas
for {
    rec, err := dec.Record()
    if errors.Is(err, io.EOF) { break }
    // rec.SeqID, rec.Source, rec.Type,
    // rec.Start, rec.End,       // 1-based inclusive
    // rec.Score, rec.ScoreKnown,
    // rec.Strand, rec.Phase,
    // rec.Attributes map[string][]string
    // rec.Attribute("Parent")      // first value
    // rec.AttributeList("Parent") // all values
}
```

GFF3 coordinates are 1-based inclusive (the opposite of BED).
Multi-value attributes (`Parent`, `Dbxref`, `Alias`, `Note`, …) use
comma separation within a single attribute; multiple attributes are
separated by `;`. The parser rejects `start < 1` and `end < start`,
and the validator's `ModeStrict` cross-references `Parent` IDs
against the document's `ID` set. The encoder is round-trippable:
`htsfile roundtrip foo.gff > out.gff` re-emits a clean file.

### FASTA Parsing

The `fasta` package reads and writes FASTA sequence records:

```go
dec := fasta.NewDecoder(r)
for {
    rec, err := dec.Record()
    if errors.Is(err, io.EOF) { break }
    // rec.ID, rec.Description, rec.Sequence, rec.Raw
    // rec.Length() returns the base count
}
```

The parser accepts any sequence-line wrapping (or no wrapping at
all) and stores the sequence as a byte slice with embedded
whitespace stripped. The encoder re-wraps at 80 chars by default
(override with `enc.SetLineWidth(n)`). The validator's `ModeStrict`
checks every byte against a configurable alphabet (default:
IUPAC nucleotide codes plus `*` and `-`). The sibling `fasta/fai`
package reads `.fai` index files for random-access into large
FASTA files; index generation is out of scope for v1.

### BCF Parsing

The `bcf` package reads BCF v2.1 (the binary counterpart of VCF):

```go
dec := bcf.NewDecoder(f)
hdr, _ := dec.Header()   // required; returns Version + contigs + filters + infos
for {
    rec, err := dec.Record()
    if errors.Is(err, io.EOF) { break }
    // rec.CHROM (resolved name), rec.Pos (1-based),
    // rec.ID, rec.Ref, rec.Alt, rec.Qual + rec.QualKnown,
    // rec.Filter (resolved IDs), rec.Info
}
```

BCF v1 is read-only — the package parses the shared per-record
fields (CHROM, POS, ID, REF, ALT, FILTER, QUAL, INFO) but skips
the per-sample FORMAT column. A BCF writer (with full round-trip
support) is queued for v2. The validator's `ModeStrict` cross-
references undeclared INFO keys.

### Tabix / CSI Index Parsing

The `tabix` and `csi` packages read the two general-purpose
position indexes used by samtools/htslib:

```go
// Tabix (magic "TBI\x01")
idx, _ := tabix.NewReader(f).Read()
for _, ref := range idx.References {
    // ref.Bins (bin number → chunks), ref.Intervals (linear index),
    // ref.NoCoor (unmapped count)
    _ = tabix.BinFromRange(1000, 5000) // bin numbers covering [1000, 5000)
}

// CSI (magic "CSI\x01"; version-aware with min_shift/depth)
idx, _ := csi.NewReader(f).Read()
// idx.MinShift (window size exponent), idx.Depth (bin levels),
// idx.WindowSize() = 1 << MinShift bp
```

Both packages expose `NewReader`/`Read` plus a symmetric
`NewWriter`/`Write(*Index)` (added in v0.4.0). The CLI's
`roundtrip` subcommand works for `.tbi` and `.csi` files.
The Tabix format has an
optional `n_no_coor` field per reference that creates an
ambiguity for multi-reference legacy files (the decoder reads
it unconditionally and treats EOF as the absence signal). CSI
adds an explicit version byte and avoids that ambiguity; if you
need to read pre-samtools-1.10 multi-reference Tabix files,
CSI is the safer choice.

## CLI

A small command-line tool ships in `cmd/htsfile`. Format is
auto-detected from the file extension (`.vcf`, `.bed`, `.sam`,
`.bam`, `.gff`, `.gff3`, `.fa`, `.fasta`, `.fna`, `.ffn`,
`.faa`, `.fq`, `.fastq`, `.bcf`, `.tbi`, `.csi`); for stdin
(`-`) or unknown extensions, the first non-empty line is sniffed
(`##fileformat` → VCF; `track`, `browser`, or a tab-separated
`chrom\tstart\tend` triple → BED; `@HD` or an 11-column
tab-separated record with integer FLAG/MAPQ → SAM; the
gzip+FEXTRA prefix → BAM (this also matches BCF's outer
wrapper — use the `.bcf` extension to disambiguate; BCF is
not auto-sniffed from stdin); `##gff-version` or a 9-column
tab-separated feature line → GFF3; a `>` header line →
FASTA; a `@` header line followed by a sequence line and a
`+` line → FASTQ. Tabix (`.tbi`) and CSI (`.csi`) are not
sniffed — use the extension.

```sh
# Validate a file
htsfile check sample.vcf                  # default mode (warns on issues)
htsfile check --strict sample.vcf         # strict spec-compliance checks
htsfile check --strict --warnings file.vcf   # exit non-zero on warnings too
htsfile check peaks.bed                   # BED is auto-detected by extension
htsfile check reads.sam                   # SAM is auto-detected by extension
htsfile check aligned.bam                 # BAM is auto-detected by extension
htsfile check genes.gff3                  # GFF3 is auto-detected by extension
htsfile check sequences.fa                # FASTA is auto-detected by extension
htsfile check reads.fq                    # FASTQ is auto-detected by extension
htsfile check sample.bcf                  # BCF is auto-detected by extension
htsfile check sample.bcf.tbi              # Tabix index (.tbi) by extension
htsfile check sample.bcf.csi              # CSI index (.csi) by extension

# Re-emit a file (verifies round-trip)
htsfile roundtrip sample.vcf > out.vcf
htsfile roundtrip --crlf sample.vcf       # CR+LF line endings
htsfile roundtrip peaks.bed > out.bed     # BED round-trip
htsfile roundtrip reads.sam > out.sam     # SAM round-trip
htsfile roundtrip genes.gff3 > out.gff3   # GFF3 round-trip
htsfile roundtrip sequences.fa > out.fa  # FASTA round-trip (80-char wrapping)
# roundtrip does not support .bam or .bcf files (both are
# read-only in v1). .bam and .bcf round-trip is supported
# as of v0.5.0. .tbi and .csi round-trip is supported.

# Print a header summary
htsfile summary sample.vcf
htsfile summary peaks.bed                 # BED summary: track, type, records, chroms
htsfile summary reads.sam                 # SAM summary: version, @SQ/@RG/@PG counts, records
htsfile summary aligned.bam               # BAM summary: version, @SQ/@RG/@PG counts, records
htsfile summary genes.gff3                # GFF3 summary: version, pragmas, seqids, types, records
htsfile summary sequences.fa              # FASTA summary: records, bases, unique IDs
htsfile summary sample.bcf                # BCF summary: fileformat, contigs, filters, records
htsfile summary sample.bcf.tbi            # Tabix summary: references, bins, chunks, intervals
htsfile summary sample.bcf.csi            # CSI summary: min_shift, depth, references, bins, chunks
```

## Status

### vcf/ package (VCF v4.5)

| Feature                                | Status |
|----------------------------------------|:------:|
| File header (`##fileformat` etc.)      |   ✅   |
| 8 fixed columns (CHROM–INFO)           |   ✅   |
| FORMAT + sample columns + GT           |   ✅   |
| Reserved INFO/FORMAT keys (DP/AF/GT…)  |   ✅   |
| Symbolic SVs (DEL/INS/DUP/INV/CNV)     |   ✅   |
| Breakend notation (`t]p]` etc.)        |   ✅   |
| Reserved SV INFO keys (SVLEN, CIPOS…)  |   ✅   |
| Phasing (PS / PQ)                      |   ✅   |
| PSL / PSO / PSQ lists                  |   ✅   |
| Base modifications (M\* fields)         |   🔜   |
| Tandem repeats (RN/RUS/RUC/RB…)        |   🔜   |
| Local alleles (LAA / LAD / LPL…)       |   🔜   |
| VCF v4.1–v4.5 acceptance (any v4.x)     |   ✅   |
| BCF v2.1 (binary counterpart)           |   ❌   |

### bed/ package (BED v1)

| Feature                                | Status |
|----------------------------------------|:------:|
| BED3–BED12 with column-count auto-detect |  ✅   |
| Track lines (parsed and round-tripped) |   ✅   |
| Browser lines (skipped, not preserved)  |   ✅   |
| BED12 block invariants (blockCount, sum, overlap) | ✅ |
| `ModeDefault` and `ModeStrict` validators | ✅ |
| BEDGraph / narrowPeak / broadPeak      |   ❌   |
| bigBed (binary indexed)                |   ❌   |

### sam/ package (SAM v1)

| Feature                                | Status |
|----------------------------------------|:------:|
| `@HD` / `@SQ` / `@RG` / `@PG` / `@CO` header lines | ✅ |
| 11-column alignment records (QNAME–QUAL) | ✅ |
| Optional tag fields (TAG:TYPE:VALUE) | ✅ |
| `ModeDefault` and `ModeStrict` validators | ✅ |
| Lenient parsing of malformed lines | ✅ |
| CIGAR string parsing (semantic) | ❌ |
| Read/write of unknown `@XX` header lines | ✅ |

### bam/ package (BAM v1)

| Feature                                | Status |
|----------------------------------------|:------:|
| BGZF block decoder (internal/bgzf)     |   ✅   |
| `BAM\1` magic + reference dictionary   |   ✅   |
| Embedded SAM header text (HD/SQ/RG/PG/CO) | ✅   |
| 11-column alignment records + tags     |   ✅   |
| 4-bit packed SEQ / CIGAR / qual        |   ✅   |
| All 12 BAM tag type codes (A, c, C, s, S, i, I, f, d, Z, H, B) | ✅ |
| `ModeDefault` and `ModeStrict` validators | ✅ |
| BAM **writer** (round-trip)             |   ❌   |
| BAI index reader (`bam/bai`)            |   ✅   |
| BAI-based region queries                |   🔜   |

### gff/ package (GFF3)

| Feature                                | Status |
|----------------------------------------|:------:|
| `##gff-version` + `##sequence-region` + `##date` pragmas | ✅ |
| All other `##` pragmas preserved verbatim for round-trip | ✅ |
| 9-column feature records (SeqID/Source/Type/Start/End/Score/Strand/Phase/Attributes) | ✅ |
| 1-based inclusive coordinates (start ≥ 1, end ≥ start) | ✅ |
| Strand: +, -, ., ?                      |   ✅   |
| Phase: 0, 1, 2, .                      |   ✅   |
| Multi-value attributes (`Parent`, `Dbxref`, …) | ✅ |
| `###` record terminator (ends feature section) | ✅ |
| Round-trippable `Encoder`              |   ✅   |
| `ModeDefault` and `ModeStrict` validators (strict = Parent ID cross-references) | ✅ |
| FASTA emission after `###`             |   ❌   |

### fasta/ package (FASTA)

| Feature                                | Status |
|----------------------------------------|:------:|
| `>` header + multi-line sequence records | ✅   |
| Sequence bytes preserved verbatim (whitespace stripped) | ✅ |
| Configurable sequence-line wrapping on encode (default 80) | ✅ |
| `ModeDefault` and `ModeStrict` validators (strict = alphabet check) | ✅ |
| Header parsing (ID + Description split on whitespace) | ✅ |
| `.fai` index reader (`fasta/fai`, read-only) | ✅ |
| `.fai` index generation                |   ❌   |
| Multi-sequence FASTA (multi-FASTA, .mfa) | ❌  |

### fastq/ package (FASTQ v1)

| Feature                                | Status |
|----------------------------------------|:------:|
| `@` header + sequence + `+` + qualities records | ✅ |
| Quality bytes preserved verbatim        | ✅ |
| Round-trippable encoder (LF / CR+LF)    | ✅ |
| `ModeDefault` and `ModeStrict` validators (strict = IUPAC alphabet) | ✅ |
| Phred+33 / Phred+64 auto-detection (first record) | ✅ |
| Manual offset override via `SetPhredOffset` | ✅ |
| Per-record `PhredOffset` carried on Record | ✅ |
| Solexa Solexa+64 encoding support    |   ❌  |

### bcf/ package (BCF v2.1)

| Feature                                | Status |
|----------------------------------------|:------:|
| BGZF block decoder (internal/bgzf)     |   ✅   |
| `BCF\x02\x01` magic + reference dictionary |   ✅ |
| Embedded VCF-style header text (fileformat, contig, FILTER, INFO, FORMAT) | ✅ |
| 9-column feature records equivalent: CHROM/POS/ID/REF/ALT/FILTER/QUAL/INFO | ✅ |
| 4-bit packed SEQ / packed CIGAR / qual |   ✅   |
| All 5 BCF INFO type codes (Integer, Float, Flag, Character, String) | ✅ |
| `ModeDefault` and `ModeStrict` validators | ✅ |
| **BCF writer (record-side)**           |   ✅   |
| Per-sample FORMAT column               |   ❌   |
| BCF 2.2 (latest version) support       |   ❌   |

### tabix/ package (Tabix v1)

| Feature                                | Status |
|----------------------------------------|:------:|
| `TBI\x01` magic + per-reference bin and linear index |   ✅   |
| Bin number → chunks resolution (`BinFromRange`, `LookupChunks`) | ✅ |
| Multiple references with optional `n_no_coor` | ✅ (modern format; legacy multi-ref without n_no_coor produces a parse error) | ✅ |
| **Tabix writer (generation)**          |   ✅   |

### csi/ package (CSI v1)

| Feature                                | Status |
|----------------------------------------|:------:|
| `CSI\x01` magic + tunable `min_shift`/`depth` |   ✅   |
| Bin number → chunks resolution (`BinFromRange`, `LookupChunks`) | ✅ |
| Multiple references                    | ✅   |
| No `n_no_coor` ambiguity (version byte disambiguates) | ✅ |
| **CSI writer (generation)**            |   ✅   |

### Planned (see [ROADMAP.md](ROADMAP.md))

- **Tier 4**: CRAM v3, htsget/refget, crypt4gh. FASTQ shipped
  in v0.3.0; Tabix + CSI writers shipped in v0.4.0. BCF
  writer still deferred (per-sample FORMAT data).

## Releases

Tagged releases are published via [GoReleaser](https://goreleaser.com/):

- linux/darwin × amd64/arm64 tar.gz archives
- SHA-256 checksums
- Auto-generated release notes from [Conventional Commits](https://www.conventionalcommits.org/)
- The Go module is tagged as `github.com/inful/htsfile/v0.0.1`, `v0.1.0`, etc.

To cut a release locally:

```sh
git tag -s v0.1.0 -m "Release v0.1.0"
git push origin v0.1.0
# GitHub Actions runs goreleaser and publishes the release.
```

## Testing

```sh
make test          # unit + golden tests, race detector
make cover         # coverage report
make lint          # golangci-lint
```

## Specification

Conforms to [VCFv4.5 and BCFv2.2](https://github.com/samtools/hts-specs/blob/master/VCFv4.5.tex).

## License

MIT — see [LICENSE](LICENSE).