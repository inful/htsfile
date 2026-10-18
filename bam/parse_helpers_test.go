package bam

import (
	"bytes"
	"encoding/binary"
	"math"
	"strconv"
)

// Test-only helpers for constructing synthetic BAM streams
// for the parser tests. These are kept in parse_helpers_test.go
// (with test- prefixed names) to avoid colliding with the
// production encodeAuxTags / baseCode / encodeTagValue /
// splitComma in write.go.
//
// The canonical encoding functions live in write.go and are
// what production callers use. These test-prefixed copies
// exist solely to keep the parser's existing fixture-building
// helpers working without rewriting every test.

func testBaseCode(c byte) byte {
	switch c {
	case '=':
		return 0
	case 'A', 'a':
		return 1
	case 'C', 'c':
		return 2
	case 'M', 'm':
		return 3
	case 'G', 'g':
		return 4
	case 'R', 'r':
		return 5
	case 'S', 's':
		return 6
	case 'V', 'v':
		return 7
	case 'T', 't':
		return 8
	case 'W', 'w':
		return 9
	case 'Y', 'y':
		return 10
	case 'H', 'h':
		return 11
	case 'K', 'k':
		return 12
	case 'D', 'd':
		return 13
	case 'B', 'b':
		return 14
	case 'N', 'n':
		return 15
	default:
		return 15
	}
}

func testEncodeAuxTags(tags []Tag) []byte {
	var buf bytes.Buffer
	for _, t := range tags {
		buf.WriteString(t.Key)
		buf.WriteByte(byte(t.Type))
		buf.Write(testEncodeTagValue(t))
	}
	return buf.Bytes()
}

func testEncodeTagValue(t Tag) []byte {
	var buf bytes.Buffer
	switch t.Type {
	case TagA:
		if len(t.Value) > 0 {
			buf.WriteByte(t.Value[0])
		} else {
			buf.WriteByte(0)
		}
	case Tagc:
		v := parseInt64(t.Value)
		buf.WriteByte(byte(int8(v)))
	case TagC:
		v := parseUint64(t.Value)
		buf.WriteByte(byte(v))
	case Tags:
		v := parseInt64(t.Value)
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(int16(v)))
		buf.Write(b[:])
	case TagS:
		v := parseUint64(t.Value)
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(v))
		buf.Write(b[:])
	case Tagi:
		v := parseInt64(t.Value)
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(int32(v)))
		buf.Write(b[:])
	case TagI:
		v := parseUint64(t.Value)
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(v))
		buf.Write(b[:])
	case Tagf:
		var b [4]byte
		var f float32
		parseFloat(t.Value, &f)
		binary.LittleEndian.PutUint32(b[:], math.Float32bits(f))
		buf.Write(b[:])
	case Tagd:
		var b [8]byte
		f, _ := strconv.ParseFloat(t.Value, 64)
		binary.LittleEndian.PutUint64(b[:], math.Float64bits(f))
		buf.Write(b[:])
	case TagZ, TagH:
		buf.WriteString(t.Value)
		buf.WriteByte(0)
	case TagB:
		sub := t.Type
		buf.WriteByte(byte(sub))
		parts := testSplitComma(t.Value)
		var c [4]byte
		binary.LittleEndian.PutUint32(c[:], uint32(len(parts)))
		buf.Write(c[:])
		for _, p := range parts {
			switch sub {
			case Tagc:
				v := parseInt64(p)
				buf.WriteByte(byte(int8(v)))
			case TagC:
				v := parseUint64(p)
				buf.WriteByte(byte(v))
			case Tags:
				v := parseInt64(p)
				var b [2]byte
				binary.LittleEndian.PutUint16(b[:], uint16(int16(v)))
				buf.Write(b[:])
			case TagS:
				v := parseUint64(p)
				var b [2]byte
				binary.LittleEndian.PutUint16(b[:], uint16(v))
				buf.Write(b[:])
			case Tagi:
				v := parseInt64(p)
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], uint32(int32(v)))
				buf.Write(b[:])
			case TagI:
				v := parseUint64(p)
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], uint32(v))
				buf.Write(b[:])
			case Tagf:
				var f float32
				parseFloat(p, &f)
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], math.Float32bits(f))
				buf.Write(b[:])
			}
		}
	}
	return buf.Bytes()
}

func testSplitComma(s string) []string {
	if s == "" {
		return nil
	}
	var out []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == ',' {
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	out = append(out, s[start:])
	return out
}
