# htsfile

[![CI](https://github.com/inful/htsfile/actions/workflows/ci.yml/badge.svg)](https://github.com/inful/htsfile/actions/workflows/ci.yml)
[![Go Report Card](https://goreportcard.com/badge/github.com/inful/htsfile)](https://goreportcard.com/report/github.com/inful/htsfile)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`htsfile` is a Go library for parsing, writing, and validating bioinformatics
file formats. It currently targets four specifications from
[samtools/hts-specs](https://github.com/samtools/hts-specs):

- **VCF v4.5** — Variant Call Format (accepts any v4.x minor)
- **BED v1** — Browser Extensible Data (UCSC-derived, BED3 through BED12)
- **SAM v1** — Sequence Alignment/Map (header + 11-column records + optional tags)
- **BAM v1** — Binary SAM (BGZF-compressed) with read-only BAI index reader

See [ROADMAP.md](ROADMAP.md) for the planned expansion (BAM, GFF3, FASTA,
FASTQ, BCF v2.1, Tabix/CSI).

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

## CLI

A small command-line tool ships in `cmd/htsfile`. Format is
auto-detected from the file extension (`.vcf`, `.bed`, `.sam`,
`.bam`); for stdin (`-`) or unknown extensions, the first
non-empty line is sniffed (`##fileformat` → VCF; `track`,
`browser`, or a tab-separated `chrom\tstart\tend` triple → BED;
`@HD` or an 11-column tab-separated record with integer
FLAG/MAPQ → SAM; the gzip+FEXTRA prefix → BAM).

```sh
# Validate a file
htsfile check sample.vcf                  # default mode (warns on issues)
htsfile check --strict sample.vcf         # strict spec-compliance checks
htsfile check --strict --warnings file.vcf   # exit non-zero on warnings too
htsfile check peaks.bed                   # BED is auto-detected by extension
htsfile check reads.sam                   # SAM is auto-detected by extension
htsfile check aligned.bam                 # BAM is auto-detected by extension

# Re-emit a file (verifies round-trip)
htsfile roundtrip sample.vcf > out.vcf
htsfile roundtrip --crlf sample.vcf       # CR+LF line endings
htsfile roundtrip peaks.bed > out.bed     # BED round-trip
htsfile roundtrip reads.sam > out.sam     # SAM round-trip
# roundtrip does not support .bam files (BAM v1 is read-only).

# Print a header summary
htsfile summary sample.vcf
htsfile summary peaks.bed                 # BED summary: track, type, records, chroms
htsfile summary reads.sam                 # SAM summary: version, @SQ/@RG/@PG counts, records
htsfile summary aligned.bam               # BAM summary: version, @SQ/@RG/@PG counts, records
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

### Planned (see [ROADMAP.md](ROADMAP.md))

- **Tier 2 (next)**: GFF3, FASTA
- **Tier 3**: BCF v2.1, Tabix, CSI
- **Tier 4**: CRAM v3, htsget/refget, crypt4gh, FASTQ

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