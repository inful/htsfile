#!/usr/bin/env python3
"""Mutation-test the BCF wire-format guard.

Each mutation reintroduces one of the defects that were found by
the external oracle (bcftools/htslib) and that the package's own
round-trip tests could not see. The guard
bcf.TestEncoderMatchesReferenceBytes must FAIL for every one of
them.

Run:  python3 tools/mutate_bcf_wire.py

The script restores the tree afterwards, so it is safe to run on
a dirty checkout - but run it with a clean tree anyway, so that
"restored" is unambiguous.
"""

import io
import subprocess
import sys

MUTATIONS = [
    (
        "n_sample/n_fmt packing swapped",
        "bcf/write.go",
        "binary.LittleEndian.PutUint32(prefix[20:24], nSample|uint32(nFmt)<<24)",
        "binary.LittleEndian.PutUint32(prefix[20:24], nSample<<8|uint32(nFmt))",
    ),
    (
        "string written as type 7 + uint32 length (not a char vector)",
        "bcf/write.go",
        "\twriteBCFString(w, s)\n}",
        "\tw.WriteByte(typeString)\n"
        "\tvar lb [4]byte\n"
        "\tbinary.LittleEndian.PutUint32(lb[:], uint32(len(s)))\n"
        "\tw.Write(lb[:])\n"
        "\tw.WriteString(s)\n}",
    ),
    (
        "GT written as a string instead of an allele vector",
        "bcf/format.go",
        "\tif def.ID == \"GT\" {\n\t\treturn typeInt8\n\t}",
        "\tif def.ID == \"GT\" {\n\t\treturn typeString\n\t}",
    ),
    (
        "all fmt_keys written first instead of one triplet per field",
        "bcf/write.go",
        "\t\twriteTypedInt(&buf, int32(hdr.IDXFor(name)))\n"
        "\t\tif err := writeFormatField(&buf, r, hdr, version, keyIdx); err != nil {\n"
        "\t\t\treturn nil, err\n\t\t}\n\t}",
        "\t\t_ = name\n"
        "\t\twriteTypedInt(&buf, int32(hdr.IDXFor(name)))\n\t}\n"
        "\tfor _, keyIdx := range r.Format {\n"
        "\t\tif err := writeFormatField(&buf, r, hdr, version, keyIdx); err != nil {\n"
        "\t\t\treturn nil, err\n\t\t}\n\t}",
    ),
    (
        "FORMAT key written as the FORMAT index instead of the header dictionary index",
        "bcf/write.go",
        "writeTypedInt(&buf, int32(hdr.IDXFor(name)))",
        "writeTypedInt(&buf, int32(keyIdx))",
    ),
    (
        "INFO key written as the INFO index instead of the header dictionary index",
        "bcf/write.go",
        "writeTypedInt(w, int32(hdr.IDXFor(key)))",
        "writeTypedInt(w, int32(keyIdx))",
    ),
]


def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def main():
    failures = []
    for name, path, old, new in MUTATIONS:
        src = io.open(path, encoding="utf-8").read()
        if old not in src:
            print(f"SKIP  {name}: pattern not found in {path}")
            failures.append(name)
            continue
        io.open(path, "w", encoding="utf-8").write(src.replace(old, new, 1))
        r = run(
            "go test -count=1 -run 'TestEncoderMatchesReferenceBytes' ./bcf/ 2>&1"
        )
        run("gofmt -w bcf/")
        caught = r.returncode != 0
        print(f"{'CAUGHT' if caught else 'MISSED'}  {name}")
        if not caught:
            failures.append(name)
        io.open(path, "w", encoding="utf-8").write(src)

    # The tree must be back to green.
    r = run("go test -count=1 ./bcf/ 2>&1")
    if r.returncode != 0:
        print("RESTORE FAILED - tree is not green:")
        print(r.stdout[-2000:])
        failures.append("restore")
    else:
        print("OK     tree restored and green")

    if failures:
        print(f"\n{len(failures)} mutation(s) NOT caught: {failures}")
        return 1
    print("\nAll mutations caught.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
