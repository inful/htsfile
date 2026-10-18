package fai

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

func TestWrite_Basic(t *testing.T) {
	idx := &Index{References: []*Reference{
		{Name: "seq1", Length: 100, Offset: 6, LineBases: 80, LineWidth: 81},
		{Name: "seq2", Length: 200, Offset: 200, LineBases: 100, LineWidth: 101},
	}}
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(idx); err != nil {
		t.Fatalf("Write: %v", err)
	}
	want := "seq1\t100\t6\t80\t81\nseq2\t200\t200\t100\t101\n"
	if buf.String() != want {
		t.Errorf("got %q, want %q", buf.String(), want)
	}
}

func TestWrite_Empty(t *testing.T) {
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(&Index{}); err != nil {
		t.Fatalf("Write empty: %v", err)
	}
	if buf.Len() != 0 {
		t.Errorf("empty write produced %q", buf.String())
	}
	if err := NewWriter(&buf).Write(nil); err != nil {
		t.Fatalf("Write nil: %v", err)
	}
}

func TestWrite_RejectsTabInName(t *testing.T) {
	idx := &Index{References: []*Reference{
		{Name: "bad\tname", Length: 1, Offset: 0, LineBases: 1, LineWidth: 2},
	}}
	if err := NewWriter(io.Discard).Write(idx); err == nil {
		t.Error("expected error on tab in name")
	}
}

func TestWrite_RejectsNewlineInName(t *testing.T) {
	idx := &Index{References: []*Reference{
		{Name: "bad\nname", Length: 1, Offset: 0, LineBases: 1, LineWidth: 2},
	}}
	if err := NewWriter(io.Discard).Write(idx); err == nil {
		t.Error("expected error on newline in name")
	}
}

func TestWrite_RoundTripsThroughReader(t *testing.T) {
	idx := &Index{References: []*Reference{
		{Name: "chr1", Length: 248956422, Offset: 6, LineBases: 80, LineWidth: 81},
		{Name: "chr2", Length: 5000000, Offset: 320000000, LineBases: 80, LineWidth: 81},
	}}
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(idx); err != nil {
		t.Fatalf("Write: %v", err)
	}
	parsed, err := NewReader(&buf).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(parsed.References) != len(idx.References) {
		t.Fatalf("References: got %d, want %d", len(parsed.References), len(idx.References))
	}
	for i, r := range idx.References {
		got := parsed.References[i]
		if got.Name != r.Name || got.Length != r.Length || got.Offset != r.Offset ||
			got.LineBases != r.LineBases || got.LineWidth != r.LineWidth {
			t.Errorf("ref %d: got %+v, want %+v", i, got, r)
		}
	}
}

// Build tests

func TestBuild_SingleLineNoWrap(t *testing.T) {
	// >seq1\nACGTACGT\n  — header is 6 bytes ("X>seq1\n"),
	// first base at offset 6, sequence is one line of 8 bases
	// with lineBases=8 and lineWidth=9 (8 + LF).
	src := ">seq1\nACGTACGT\n"
	idx, err := Build(strings.NewReader(src))
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	if len(idx.References) != 1 {
		t.Fatalf("References: %d", len(idx.References))
	}
	r := idx.References[0]
	if r.Name != "seq1" {
		t.Errorf("Name: %q", r.Name)
	}
	if r.Length != 8 {
		t.Errorf("Length: %d", r.Length)
	}
	if r.Offset != 6 {
		t.Errorf("Offset: %d", r.Offset)
	}
	if r.LineBases != 8 {
		t.Errorf("LineBases: %d", r.LineBases)
	}
	if r.LineWidth != 9 {
		t.Errorf("LineWidth: %d", r.LineWidth)
	}
}

func TestBuild_WrappedSequence(t *testing.T) {
	// Two records, second is wrapped at 60.
	src := ">a\nACGT\n>b\n" + strings.Repeat("A", 60) + "\n" + strings.Repeat("C", 5) + "\n"
	idx, err := Build(strings.NewReader(src))
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	if len(idx.References) != 2 {
		t.Fatalf("References: %d", len(idx.References))
	}
	a := idx.References[0]
	if a.Name != "a" || a.Length != 4 || a.Offset != 3 || a.LineBases != 4 || a.LineWidth != 5 {
		t.Errorf("a: %+v", a)
	}
	b := idx.References[1]
	if b.Name != "b" || b.Length != 65 || b.LineBases != 60 || b.LineWidth != 61 {
		t.Errorf("b: %+v", b)
	}
	// 'b' header is ">b\n" (3 bytes) at offset 8, so first
	// base is at 11.
	if b.Offset != 11 {
		t.Errorf("b.Offset: got %d, want 11", b.Offset)
	}
}

func TestBuild_CRLF(t *testing.T) {
	// Windows line endings: ">seq1\r\nACGT\r\n" — first base
	// at offset 7, lineBases=4, lineWidth=6 (4 + CRLF).
	src := ">seq1\r\nACGT\r\n"
	idx, err := Build(strings.NewReader(src))
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	r := idx.References[0]
	if r.Offset != 7 {
		t.Errorf("Offset: %d, want 7", r.Offset)
	}
	if r.LineBases != 4 {
		t.Errorf("LineBases: %d, want 4", r.LineBases)
	}
	if r.LineWidth != 6 {
		t.Errorf("LineWidth: %d, want 6 (CRLF)", r.LineWidth)
	}
}

func TestBuild_BlankLinesBetween(t *testing.T) {
	src := ">a\nAAAA\n\n>b\nCCCC\n"
	idx, err := Build(strings.NewReader(src))
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	if len(idx.References) != 2 {
		t.Fatalf("References: %d", len(idx.References))
	}
	if idx.References[0].Length != 4 || idx.References[1].Length != 4 {
		t.Errorf("Lengths: %d, %d", idx.References[0].Length, idx.References[1].Length)
	}
}

func TestBuild_HeaderWithDescription(t *testing.T) {
	// ID is the first whitespace-delimited token; description
	// is everything after and must not be captured in Name.
	src := ">chr1 Homo sapiens chromosome 1\nACGT\n"
	idx, err := Build(strings.NewReader(src))
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	if idx.References[0].Name != "chr1" {
		t.Errorf("Name: %q, want chr1", idx.References[0].Name)
	}
}

func TestBuild_TabInHeader(t *testing.T) {
	// Tab is also a valid header ID/description separator.
	src := ">chr1\tbuild37\nACGT\n"
	idx, err := Build(strings.NewReader(src))
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	if idx.References[0].Name != "chr1" {
		t.Errorf("Name: %q, want chr1", idx.References[0].Name)
	}
}

func TestBuild_EmptyInput(t *testing.T) {
	_, err := Build(strings.NewReader(""))
	if !errors.Is(err, errMissingGT) {
		t.Errorf("empty input: got err %v, want errMissingGT", err)
	}
}

func TestBuild_MissingHeader(t *testing.T) {
	// First non-empty line is not a '>' header.
	_, err := Build(strings.NewReader("ACGTACGT\n"))
	if !errors.Is(err, errMissingGT) {
		t.Errorf("missing header: got err %v, want errMissingGT", err)
	}
}

func TestBuild_EmptySequence(t *testing.T) {
	// Header followed by another header (no sequence between).
	src := ">a\n>b\nACGT\n"
	idx, err := Build(strings.NewReader(src))
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	if len(idx.References) != 2 {
		t.Fatalf("References: %d", len(idx.References))
	}
	if idx.References[0].Length != 0 {
		t.Errorf("empty record length: %d, want 0", idx.References[0].Length)
	}
	// Offset/lineBases/lineWidth for the empty record stay at
	// the zero value — there is no first sequence line to
	// measure.
	if idx.References[0].Offset != 0 || idx.References[0].LineBases != 0 || idx.References[0].LineWidth != 0 {
		t.Errorf("empty record: %+v", idx.References[0])
	}
}

func TestBuildAndWrite_RoundTrip(t *testing.T) {
	src := ">chr1 description one\n" +
		strings.Repeat("A", 60) + "\n" +
		strings.Repeat("C", 60) + "\n" +
		">chr2\n" +
		strings.Repeat("G", 30) + "\n"
	idx, err := Build(strings.NewReader(src))
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	var buf bytes.Buffer
	if err := NewWriter(&buf).Write(idx); err != nil {
		t.Fatalf("Write: %v", err)
	}
	parsed, err := NewReader(&buf).Read()
	if err != nil {
		t.Fatalf("Read back: %v", err)
	}
	if len(parsed.References) != 2 {
		t.Fatalf("References: %d", len(parsed.References))
	}
	r0 := parsed.References[0]
	if r0.Name != "chr1" || r0.Length != 120 || r0.LineBases != 60 || r0.LineWidth != 61 {
		t.Errorf("chr1 round-trip: %+v", r0)
	}
	r1 := parsed.References[1]
	if r1.Name != "chr2" || r1.Length != 30 || r1.LineBases != 30 || r1.LineWidth != 31 {
		t.Errorf("chr2 round-trip: %+v", r1)
	}
}