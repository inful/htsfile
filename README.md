# htsfile

[![CI](https://github.com/inful/htsfile/actions/workflows/ci.yml/badge.svg)](https://github.com/inful/htsfile/actions/workflows/ci.yml)
[![Go Report Card](https://goreportcard.com/badge/github.com/inful/htsfile)](https://goreportcard.com/report/github.com/inful/htsfile)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`htsfile` is a Go library for parsing, writing, and validating bioinformatics
file formats. It currently targets the **VCF v4.5** specification
([samtools/hts-specs](https://github.com/samtools/hts-specs)).

The implementation aims to be:

- **Spec-compliant** — covers the v4.5 grammar for core records, INFO/FORMAT
  metadata, symbolic structural variants, and breakend notation. Base
  modifications (`M*` fields), tandem repeats (`RUS/RUL/RUC/RB/...`), and local
  allele encoding (`LAA/LAD/LPL/LGL/...`) are deferred to a follow-up release.
- **Streaming** — `Decoder` reads headers and records incrementally from any
  `io.Reader`. No buffering of the whole file.
- **Round-trippable** — reading a VCF and writing it back produces output that
  round-trips through `Decode` again. Field ordering inside structured
  meta-information lines is preserved when known.
- **Validating** — `Validator` runs header-declared cardinality (`Number=1/A/R/G/...`)
  and type (`Integer/Float/Flag/Character/String`) checks against every record
  and surfaces every issue it finds, not just the first.

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

## CLI

A small command-line tool ships in `cmd/htsfile`:

```sh
# Validate a file
htsfile check sample.vcf                  # default mode (warns on issues)
htsfile check --strict sample.vcf         # strict spec-compliance checks
htsfile check --strict --warnings file.vcf   # exit non-zero on warnings too

# Re-emit a file (verifies round-trip)
htsfile roundtrip sample.vcf > out.vcf
htsfile roundtrip --crlf sample.vcf       # CR+LF line endings

# Print a header summary
htsfile summary sample.vcf
```

## Status

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
| BCF (binary counterpart)               |   ❌   |

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