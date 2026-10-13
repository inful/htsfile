package vcf

import (
	"bytes"
	"io"
	"os"
	"testing"
)

// BenchmarkParseSpecExample measures parsing throughput on the canonical
// spec example fixture.
func BenchmarkParseSpecExample(b *testing.B) {
	data, err := os.ReadFile("testdata/spec_example.vcf")
	if err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		dec := NewDecoder(bytes.NewReader(data))
		if _, err := dec.Header(); err != nil {
			b.Fatal(err)
		}
		for {
			_, err := dec.Record()
			if err == io.EOF {
				break
			}
			if err != nil {
				b.Fatal(err)
			}
		}
	}
}

// BenchmarkRoundTripSpecExample measures the full parse + encode cycle.
func BenchmarkRoundTripSpecExample(b *testing.B) {
	data, err := os.ReadFile("testdata/spec_example.vcf")
	if err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		dec := NewDecoder(bytes.NewReader(data))
		hdr, err := dec.Header()
		if err != nil {
			b.Fatal(err)
		}
		var recs []*Record
		for {
			rec, err := dec.Record()
			if err == io.EOF {
				break
			}
			if err != nil {
				b.Fatal(err)
			}
			recs = append(recs, rec)
		}
		var buf bytes.Buffer
		enc := NewEncoder(&buf)
		if err := enc.EncodeHeader(hdr); err != nil {
			b.Fatal(err)
		}
		for _, r := range recs {
			if err := enc.EncodeRecord(r); err != nil {
				b.Fatal(err)
			}
		}
	}
}

// BenchmarkValidateSpecExample measures semantic validation throughput.
func BenchmarkValidateSpecExample(b *testing.B) {
	data, err := os.ReadFile("testdata/spec_example.vcf")
	if err != nil {
		b.Fatal(err)
	}
	dec := NewDecoder(bytes.NewReader(data))
	hdr, err := dec.Header()
	if err != nil {
		b.Fatal(err)
	}
	var recs []*Record
	for {
		rec, err := dec.Record()
		if err == io.EOF {
			break
		}
		if err != nil {
			b.Fatal(err)
		}
		recs = append(recs, rec)
	}
	v := NewValidator(ModeDefault)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = v.Validate(hdr, recs)
	}
}