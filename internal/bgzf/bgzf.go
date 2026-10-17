// Package bgzf implements a streaming reader for the BGZF block
// format used to wrap BAM (and other) data inside a gzip container.
//
// BGZF (Blocked GNU Zip Format) is a constrained variant of gzip
// defined by the samtools project:
//
//   - Every block is a valid gzip stream with the FEXTRA flag set,
//     carrying a single extra subfield "BC" (SI1=0x42='B',
//     SI2=0x43='C', SLEN=2) whose value is a uint16 little-endian
//     "BSIZE": the total block size minus 1.
//   - Each BGZF block is at most 64 KiB (BSIZE <= 0xFFFF), which
//     gives BGZF its random-access property: a block index can
//     seek to any block boundary by jumping to the compressed
//     offset.
//   - The CDATA portion (compressed bytes inside a block) is a
//     raw deflate stream; we use compress/flate to decode it.
//   - Each block ends with a 4-byte CRC32 of the uncompressed
//     payload and a 4-byte ISIZE (uncompressed size mod 2^32).
//
// A BGZF "file" is just a concatenation of these blocks followed
// by a 28-byte empty EOF marker (BSIZE=27, ISIZE=0, CRC32=0).
//
// The package is zero-dependency: only stdlib compress/flate and
// hash/crc32 are used.
package bgzf

import (
	"bufio"
	"bytes"
	"compress/flate"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"io"
)

// BGZF block layout constants.
const (
	// headerSize is the fixed size of the gzip header with the
	// BGZF extra subfield. The header is:
	//   1f 8b 08 04 <4-byte MTIME> 00 ff <XLEN=6> BC 02 <BSIZE>
	headerSize = 18

	// maxBSIZE is the largest legal block size (in bytes). The
	// spec requires BSIZE <= 0xFFFF; the highest two bits are
	// reserved and must be zero.
	maxBSIZE = 0xFFFF

	// trailerSize is the size of the CRC32 + ISIZE trailer at
	// the end of each block.
	trailerSize = 8

	// emptyEOFFenceBSIZE is the BSIZE of the canonical empty
	// BGZF EOF marker block (the 28-byte "trailer" used by
	// samtools). BSIZE = 27 → total = 28 → header 18 + CDATA 2
	// (an empty deflate stream) + 8 trailer.
	emptyEOFFenceBSIZE = 27

	// ioBufSize is the internal buffer used by the read-ahead
	// bufio.Reader. 64 KiB is one full BGZF block.
	ioBufSize = 64 * 1024
)

// Errors returned by the reader.
var (
	// ErrBadMagic is returned when the first two bytes of a
	// block are not the gzip magic (0x1f 0x8b) or the
	// compression method byte is not deflate (0x08).
	ErrBadMagic = errors.New("bgzf: bad block magic")

	// ErrBadExtra is returned when the FEXTRA flag is unset or
	// the BGZF extra subfield is missing or malformed.
	ErrBadExtra = errors.New("bgzf: bad extra subfield")

	// ErrBadBSIZE is returned when BSIZE is outside the legal
	// range (must be <= 0xFFFF) or the block's claimed CDATA
	// length does not match the BSIZE math.
	ErrBadBSIZE = errors.New("bgzf: bad BSIZE")

	// ErrBadCRC is returned when the trailing CRC32 of a block
	// does not match the decompressed payload's CRC32.
	ErrBadCRC = errors.New("bgzf: CRC32 mismatch")

	// ErrBadISIZE is returned when the trailing ISIZE of a
	// block does not match the decompressed payload's length.
	ErrBadISIZE = errors.New("bgzf: ISIZE mismatch")
)

// Reader is a streaming BGZF block reader. It buffers one block
// at a time: a Read fills the internal buffer with the
// decompressed payload of the next block, then serves bytes from
// that buffer until empty. EOF is signalled by reading past the
// final block's payload (typically the 28-byte EOF marker is the
// last block seen, with an empty payload).
//
// A Reader is not safe for concurrent use.
type Reader struct {
	// br reads the compressed stream. We use a 64 KiB buffer
	// because each BGZF block is at most that size, so a single
	// block always fits in one fill.
	br *bufio.Reader

	// buf holds the decompressed payload of the current block
	// (or the next block that has been read-ahead).
	buf []byte

	// offset is the offset of the next byte to return from buf.
	offset int

	// err is the sticky error from the last Read attempt; once
	// non-nil, all future Reads return it.
	err error

	// endOfStream is set once the EOF marker has been consumed,
	// so subsequent Reads return io.EOF without trying to read
	// another block.
	endOfStream bool
}

// NewReader returns a Reader that decompresses BGZF data from r.
// The Reader buffers one block at a time; no Read on the
// underlying r happens until the first call to Read.
func NewReader(r io.Reader) *Reader {
	return &Reader{br: bufio.NewReaderSize(r, ioBufSize)}
}

// Read decompresses the next BGZF block as needed and copies
// bytes from it into p. It returns the number of bytes copied
// and io.EOF when the BGZF stream has been fully consumed.
//
// Read never returns (0, nil) per the io.Reader contract; a zero
// return is always accompanied by an error (typically io.EOF).
func (r *Reader) Read(p []byte) (int, error) {
	if r.err != nil {
		return 0, r.err
	}
	if len(p) == 0 {
		return 0, nil
	}
	// Drain the current buffer first.
	if r.offset < len(r.buf) {
		n := copy(p, r.buf[r.offset:])
		r.offset += n
		return n, nil
	}
	// Buffer is empty; load the next block.
	if r.endOfStream {
		r.err = io.EOF
		return 0, io.EOF
	}
	if err := r.nextBlock(); err != nil {
		r.err = err
		return 0, err
	}
	// nextBlock loaded a fresh payload; copy from it.
	n := copy(p, r.buf[r.offset:])
	r.offset += n
	return n, nil
}

// nextBlock reads one BGZF block from the underlying reader,
// decompresses its CDATA, and validates CRC32 + ISIZE. The EOF
// marker is recognised and recorded; once seen, the next call
// to Read returns io.EOF.
func (r *Reader) nextBlock() error {
	hdr := make([]byte, headerSize)
	n, err := io.ReadFull(r.br, hdr)
	if err != nil {
		if errors.Is(err, io.EOF) {
			// io.ReadFull only returns plain io.EOF when zero
			// bytes were read — i.e. the stream is empty (no
			// blocks at all). For "truncated mid-header" cases
			// it returns io.ErrUnexpectedEOF directly, which we
			// pass through below.
			if n == 0 {
				return io.EOF
			}
			return io.ErrUnexpectedEOF
		}
		return err
	}
	// Validate gzip magic + compression method.
	if hdr[0] != 0x1f || hdr[1] != 0x8b {
		return fmt.Errorf("%w: got %02x %02x", ErrBadMagic, hdr[0], hdr[1])
	}
	if hdr[2] != 0x08 {
		return fmt.Errorf("%w: compression method 0x%02x not deflate", ErrBadMagic, hdr[2])
	}
	// FLG must have FEXTRA (0x04) set for BGZF.
	if hdr[3]&0x04 == 0 {
		return fmt.Errorf("%w: FEXTRA flag not set", ErrBadExtra)
	}
	// XLEN (little-endian uint16 at offset 10).
	xlen := binary.LittleEndian.Uint16(hdr[10:12])
	if xlen < 6 {
		return fmt.Errorf("%w: XLEN=%d too small for BGZF subfield", ErrBadExtra, xlen)
	}
	// BGZF extra subfield: SI1='B', SI2='C', SLEN=2, BSIZE (uint16 LE).
	if hdr[12] != 0x42 || hdr[13] != 0x43 {
		return fmt.Errorf("%w: SI1=%02x SI2=%02x not 'BC'", ErrBadExtra, hdr[12], hdr[13])
	}
	subLen := binary.LittleEndian.Uint16(hdr[14:16])
	if subLen != 2 {
		return fmt.Errorf("%w: BGZF SLEN=%d, want 2", ErrBadExtra, subLen)
	}
	bsize := binary.LittleEndian.Uint16(hdr[16:18])
	if bsize > maxBSIZE {
		return fmt.Errorf("%w: BSIZE=%d > %d", ErrBadBSIZE, bsize, maxBSIZE)
	}
	// CDATA length: total block size = BSIZE + 1, of which
	// headerSize + trailerSize are non-CDATA. So CDATA = BSIZE + 1
	// - headerSize - trailerSize.
	cdataLen := int(bsize) + 1 - headerSize - trailerSize
	if cdataLen < 0 {
		return fmt.Errorf("%w: BSIZE=%d too small for header+trailer", ErrBadBSIZE, bsize)
	}

	// Read CDATA + trailer.
	tail := make([]byte, cdataLen+trailerSize)
	if _, err := io.ReadFull(r.br, tail); err != nil {
		if errors.Is(err, io.EOF) {
			return io.ErrUnexpectedEOF
		}
		return err
	}
	cdata := tail[:cdataLen]
	gotCRC := binary.LittleEndian.Uint32(tail[cdataLen : cdataLen+4])
	gotISIZE := binary.LittleEndian.Uint32(tail[cdataLen+4 : cdataLen+8])

	// Decompress CDATA (raw deflate, no gzip header).
	fr := flate.NewReader(bytes.NewReader(cdata))
	decompressed, err := io.ReadAll(fr)
	_ = fr.Close()
	if err != nil {
		return fmt.Errorf("bgzf: deflate: %w", err)
	}

	// Validate CRC32 + ISIZE.
	if uint32(len(decompressed)) != gotISIZE {
		return fmt.Errorf("%w: have %d bytes, want %d", ErrBadISIZE, len(decompressed), gotISIZE)
	}
	if crc32.ChecksumIEEE(decompressed) != gotCRC {
		return fmt.Errorf("%w: have %08x, want %08x", ErrBadCRC, crc32.ChecksumIEEE(decompressed), gotCRC)
	}

	// The 28-byte empty EOF marker (BSIZE = 27, ISIZE = 0) marks
	// the end of the BGZF stream. Record that we've seen it;
	// the next Read returns io.EOF.
	if bsize == emptyEOFFenceBSIZE && len(decompressed) == 0 {
		r.endOfStream = true
		r.buf = nil
		r.offset = 0
		return nil
	}

	r.buf = decompressed
	r.offset = 0
	return nil
}