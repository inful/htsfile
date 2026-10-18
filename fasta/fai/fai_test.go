package fai

import (
	"errors"
	"strings"
	"testing"
)

func TestRead_Basic(t *testing.T) {
	src := "seq1\t100\t100\t80\t81\n" +
		"seq2\t200\t500\t100\t101\n"
	idx, err := NewReader(strings.NewReader(src)).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(idx.References) != 2 {
		t.Fatalf("References: %d", len(idx.References))
	}
	r1 := idx.References[0]
	if r1.Name != "seq1" || r1.Length != 100 || r1.Offset != 100 {
		t.Errorf("r1: %+v", r1)
	}
	if r1.LineBases != 80 || r1.LineWidth != 81 {
		t.Errorf("r1 line bases/width: %d/%d", r1.LineBases, r1.LineWidth)
	}
}

func TestRead_Empty(t *testing.T) {
	idx, err := NewReader(strings.NewReader("")).Read()
	if err != nil {
		t.Fatal(err)
	}
	if len(idx.References) != 0 {
		t.Errorf("References: %d", len(idx.References))
	}
}

func TestRead_BadFieldCount(t *testing.T) {
	_, err := NewReader(strings.NewReader("seq1\t100\n")).Read()
	if err == nil {
		t.Fatal("expected error on wrong field count")
	}
}

func TestRead_NonNumeric(t *testing.T) {
	_, err := NewReader(strings.NewReader("seq1\tabc\t0\t80\t81\n")).Read()
	if err == nil {
		t.Fatal("expected error on non-numeric length")
	}
}

func TestRead_ScanErr(t *testing.T) {
	// Force a scanner error by closing the underlying reader
	// mid-stream.
	_, err := NewReader(&erroringReader{}).Read()
	if err == nil {
		t.Fatal("expected error from underlying reader")
	}
}

type erroringReader struct{}

func (e *erroringReader) Read(p []byte) (int, error) {
	return 0, errors.New("simulated reader error")
}
