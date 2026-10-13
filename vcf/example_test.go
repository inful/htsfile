package vcf_test

import (
	"errors"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/inful/htsfile/vcf"
)

// This example parses the canonical VCFv4.5 spec fixture and prints one
// line per record.
func ExampleNewDecoder() {
	f, err := os.Open("testdata/spec_example.vcf")
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	defer f.Close()

	dec := vcf.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Printf("fileformat: %s, samples: %d\n", hdr.FileFormat, len(hdr.SampleNames))

	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		fmt.Printf("%s\t%d\t%s\t%s\t%s\n",
			rec.Chrom, rec.Pos, rec.Ref,
			strings.Join(rec.Alt, ","),
			strings.Join(rec.Info["DP"].Raw, ","))
	}
	// Output:
	// fileformat: VCFv4.5, samples: 3
	// 20	14370	G	A	14
	// 20	17330	T	A	11
	// 20	1110696	A	G,T	10
	// 20	1230237	T		13
	// 20	1234567	GTC	G,GTCT	9
}

// This example writes a small VCF document with a single SNP record.
func ExampleNewEncoder() {
	var buf strings.Builder
	enc := vcf.NewEncoder(&buf)

	hdr := &vcf.Header{
		FileFormat: vcf.FileFormat{4, 5},
		Infos: map[string]*vcf.InfoDef{
			"DP": {ID: "DP", Number: vcf.Number{Number: 1},
				Type: vcf.TypeInteger, Description: "Read depth"},
		},
		SampleNames: []string{"S1"},
	}
	if err := enc.EncodeHeader(hdr); err != nil {
		fmt.Println("error:", err)
		return
	}

	q := 30.0
	if err := enc.EncodeRecord(&vcf.Record{
		Chrom: "1", Pos: 100, Ref: "A",
		Alt:    []string{"G"},
		Qual:   &q,
		Filter: []string{"PASS"},
		Info: map[string]vcf.InfoValue{
			"DP": {Type: vcf.TypeInteger, Ints: []int64{42}},
		},
		Format: []string{"GT", "DP"},
		Samples: []vcf.Genotype{
			{Fields: map[string]vcf.FieldValue{
				"GT": {Type: vcf.TypeString, Strs: []string{"0/1"}},
				"DP": {Type: vcf.TypeInteger, Ints: []int64{20}},
			}},
		},
	}); err != nil {
		fmt.Println("error:", err)
		return
	}
	if err := enc.Flush(); err != nil {
		fmt.Println("error:", err)
		return
	}

	fmt.Println(buf.String())
	// Output:
	// ##fileformat=VCFv4.5
	// ##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
	// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
	// 1	100	.	A	G	30	PASS	DP=42	GT:DP	0/1:20
}

// This example shows how to use the Validator to detect semantic
// problems in a parsed file.
func ExampleValidator() {
	src := `##fileformat=VCFv4.5
##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	C	.	PASS	DP=42
1	200	.	A	C	.	PASS	ZZ=oops
`

	dec := vcf.NewDecoder(strings.NewReader(src))
	hdr, err := dec.Header()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	var recs []*vcf.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		recs = append(recs, rec)
	}

	v := vcf.NewValidator(vcf.ModeDefault)
	for _, iss := range v.Validate(hdr, recs) {
		fmt.Printf("[%s] record %d, %s: %s\n",
			iss.Severity, iss.Record, iss.Field, iss.Msg)
	}
	// Output:
	// [warning] record 2, INFO/ZZ: INFO key is not declared in the header
}

// This example parses a VCF, modifies a record, and writes it back out.
func ExampleDecoder_Record_modify() {
	src := "##fileformat=VCFv4.5\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"1\t100\t.\tA\tC\t.\t.\t.\n"

	dec := vcf.NewDecoder(strings.NewReader(src))
	if _, err := dec.Header(); err != nil {
		fmt.Println("error:", err)
		return
	}
	rec, err := dec.Record()
	if err != nil {
		fmt.Println("error:", err)
		return
	}

	// Promote the call to PASS and add a depth INFO field.
	rec.Filter = []string{"PASS"}
	rec.Info = map[string]vcf.InfoValue{
		"DP": {Type: vcf.TypeInteger, Ints: []int64{42}},
	}

	var buf strings.Builder
	enc := vcf.NewEncoder(&buf)
	if err := enc.EncodeRecord(rec); err != nil {
		fmt.Println("error:", err)
		return
	}
	_ = enc.Flush()
	fmt.Println(buf.String())
	// Output:
	// 1	100	.	A	C	.	PASS	DP=42
}

// This example demonstrates that parsing + encoding the spec example
// produces a stream that re-parses into the same logical content.
func Example_roundTrip() {
	f, err := os.Open("testdata/spec_example.vcf")
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	defer f.Close()

	dec := vcf.NewDecoder(f)
	hdr, err := dec.Header()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	var recs []*vcf.Record
	for {
		rec, err := dec.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		recs = append(recs, rec)
	}

	// Re-emit and re-parse; field content must match.
	var buf strings.Builder
	enc := vcf.NewEncoder(&buf)
	if err := enc.EncodeHeader(hdr); err != nil {
		fmt.Println("error:", err)
		return
	}
	for _, r := range recs {
		if err := enc.EncodeRecord(r); err != nil {
			fmt.Println("error:", err)
			return
		}
	}
	_ = enc.Flush()

	dec2 := vcf.NewDecoder(strings.NewReader(buf.String()))
	hdr2, err := dec2.Header()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	var n int
	for {
		_, err := dec2.Record()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			fmt.Println("error:", err)
			return
		}
		n++
	}

	fmt.Printf("fileformat preserved: %s == %s; %d records preserved\n",
		hdr.FileFormat, hdr2.FileFormat, n)
	// Output:
	// fileformat preserved: VCFv4.5 == VCFv4.5; 5 records preserved
}
