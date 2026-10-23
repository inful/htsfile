package bgzf

import (
	"bytes"
	"compress/flate"
	"encoding/binary"
	"errors"
	"hash/crc32"
	"io"
)

// Writer compresses payloads into BGZF blocks written to an
// underlying [io.Writer]. The output is a sequence of valid
// BGZF blocks (each ≤ 64 KiB; BSIZE ≤ 0xFFFF), followed by the
// canonical 28-byte empty EOF marker. Bytes produced by this
// Writer round-trip cleanly through [Reader].
//
// Typical usage:
//
//	f, _ := os.Create("out.bgz")
//	defer f.Close()
//	w := bgzf.NewWriter(f)
//	w.Write([]byte("payload 1"))
//	w.Write([]byte("payload 2"))
//	w.Close()
//
// The Writer buffers at most [maxBSIZE] bytes per block;
// Write calls flush the current buffer once it crosses the
// block boundary. Close flushes any pending data, then emits
// the canonical EOF marker.
//
// A Writer is not safe for concurrent use; create one per
// goroutine.
type Writer struct {
	w io.Writer

	// blockSize is the maximum uncompressed payload per block
	// (≤ 64 KiB). The compressed payload plus header + trailer
	// fits within the spec's 64 KiB block budget for any
	// realistic deflate ratio.
	blockSize int

	// buf holds the uncompressed payload of the current block.
	// It is flushed to the underlying writer when (a) it fills
	// past blockSize on a Write call, or (b) Close is called.
	buf bytes.Buffer

	// closed reports whether Close has been called. A second
	// Close is a no-op; a Write after Close returns an error.
	closed bool
}

// NewWriter wraps w in a [Writer] using the spec-maximum
// block size (64 KiB). Block size is fixed at construction;
// callers that need a specific block size should construct
// the Writer directly.
func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w, blockSize: maxBSIZE}
}

// Write appends p to the current block. If the resulting
// block would exceed the configured size, the current block
// is flushed first and p starts a fresh block. The interface
// is [io.Writer] so callers can use a Writer wherever an
// io.Writer is expected; the only difference is that errors
// from the underlying io.Writer are reported via this Write
// return rather than a separate Flush call.
func (wr *Writer) Write(p []byte) (int, error) {
	if wr == nil {
		return 0, errors.New("bgzf: Write on nil Writer")
	}
	if wr.closed {
		return 0, errors.New("bgzf: Write after Close")
	}
	written := 0
	for len(p) > 0 {
		// Fill the current block up to blockSize.
		room := wr.blockSize - wr.buf.Len()
		n := len(p)
		if n > room {
			n = room
		}
		wr.buf.Write(p[:n])
		p = p[n:]
		written += n
		// If the block is full, flush it.
		if wr.buf.Len() >= wr.blockSize {
			if err := wr.flushBlock(); err != nil {
				return written, err
			}
		}
	}
	return written, nil
}

// Close flushes any pending data and emits the canonical
// 28-byte empty EOF marker. A second Close is a no-op.
// Close must be called to produce a valid BGZF stream —
// without it, downstream readers will hang waiting for the
// EOF marker.
func (wr *Writer) Close() error {
	if wr == nil {
		return errors.New("bgzf: Close on nil Writer")
	}
	if wr.closed {
		return nil
	}
	wr.closed = true
	if wr.buf.Len() > 0 {
		if err := wr.flushBlock(); err != nil {
			return err
		}
	}
	return wr.writeEOF()
}

// flushBlock compresses and emits the current block's
// payload. After this call wr.buf is empty.
func (wr *Writer) flushBlock() error {
	payload := wr.buf.Bytes()
	wr.buf.Reset()

	// Compress payload as raw deflate.
	var cbuf bytes.Buffer
	cbuf.Grow(len(payload) + 64) // rough upper bound
	fw, err := flate.NewWriter(&cbuf, flate.DefaultCompression)
	if err != nil {
		return err
	}
	if _, err := fw.Write(payload); err != nil {
		_ = fw.Close()
		return err
	}
	if err := fw.Close(); err != nil {
		return err
	}
	cdata := cbuf.Bytes()

	// Compute CRC32 + ISIZE of the uncompressed payload.
	crc := crc32.ChecksumIEEE(payload)
	isize := uint32(len(payload))

	// Build the block: header (18) + CDATA + trailer (8).
	// Total block size = headerSize + len(cdata) + trailerSize.
	// BSIZE = total - 1.
	totalSize := headerSize + len(cdata) + trailerSize
	bsize := uint16(totalSize - 1)
	if totalSize > maxBSIZE+1 {
		// Deflate should make this unreachable for any
		// realistic input (compressed size < uncompressed size
		// on average; only pathological inputs inflate).
		return errors.New("bgzf: compressed block exceeds 64 KiB")
	}

	// Build the gzip + BGZF header. The BGZF header is a
	// standard gzip header with FEXTRA set; the FEXTRA field
	// (XLEN bytes) starts at offset 12 — after MTIME (4 bytes)
	// and XFL (2 bytes). BGZF always writes XFL=0.
	var block bytes.Buffer
	block.Grow(totalSize)
	block.WriteByte(0x1f) // gzip magic 1
	block.WriteByte(0x8b) // gzip magic 2
	block.WriteByte(0x08) // deflate
	block.WriteByte(0x04) // FLG = FEXTRA
	// MTIME (4 bytes, zero).
	block.Write([]byte{0, 0, 0, 0})
	// XFL (2 bytes, zero — BGZF does not use this).
	block.WriteByte(0)
	block.WriteByte(0)
	// XLEN (uint16 LE) = 6 (the size of the BGZF subfield).
	var xlen [2]byte
	binary.LittleEndian.PutUint16(xlen[:], 6)
	block.Write(xlen[:])
	// BGZF subfield: SI1='B', SI2='C', SLEN=2, BSIZE uint16 LE.
	block.WriteByte('B')
	block.WriteByte('C')
	block.WriteByte(2)
	block.WriteByte(0) // SLEN is uint16 LE
	var bsizeBytes [2]byte
	binary.LittleEndian.PutUint16(bsizeBytes[:], bsize)
	block.Write(bsizeBytes[:])
	// CDATA.
	block.Write(cdata)
	// Trailer: CRC32 uint32 LE + ISIZE uint32 LE.
	var crcBytes [4]byte
	binary.LittleEndian.PutUint32(crcBytes[:], crc)
	block.Write(crcBytes[:])
	var isizeBytes [4]byte
	binary.LittleEndian.PutUint32(isizeBytes[:], isize)
	block.Write(isizeBytes[:])

	if _, err := wr.w.Write(block.Bytes()); err != nil {
		return err
	}
	return nil
}

// writeEOF emits the canonical 28-byte empty EOF marker.
// The marker is itself a valid BGZF block with BSIZE=27 and
// an empty decompressed payload (CRC32=0, ISIZE=0). The
// CDATA section contains the canonical empty-deflate-stream
// sequence `03 00` so downstream flate decoders produce an
// empty payload when this block is read.
func (wr *Writer) writeEOF() error {
	var eof [28]byte
	eof[0] = 0x1f
	eof[1] = 0x8b
	eof[2] = 0x08
	eof[3] = 0x04 // FEXTRA
	// MTIME = 0 (bytes 4-7).
	// XFL = 0 (byte 8).
	eof[8] = 0
	// OS = 0xff "unknown" (byte 9). This byte is load-bearing:
	// htslib's bgzf_check_EOF() does a byte-exact memcmp against
	// its hardcoded 28-byte marker, so writing 0x00 here makes
	// samtools/bcftools report "No BGZF EOF marker; file may be
	// truncated" on a perfectly complete file.
	eof[9] = 0xff
	// XLEN = 6 (uint16 LE → bytes 10-11).
	eof[10] = 6
	eof[11] = 0
	// BGZF subfield: SI1='B', SI2='C', SLEN=2, BSIZE=27.
	eof[12] = 'B'
	eof[13] = 'C'
	eof[14] = 2
	eof[15] = 0
	// BSIZE = 27 (uint16 LE → bytes 16-17 → 0x1B 0x00).
	eof[16] = 27
	eof[17] = 0
	// CDATA: 2 bytes of empty deflate stream (03 00).
	eof[18] = 0x03
	eof[19] = 0x00
	// CRC32 = 0, ISIZE = 0 — the trailing 8 bytes are already 0.
	_, err := wr.w.Write(eof[:])
	return err
}
