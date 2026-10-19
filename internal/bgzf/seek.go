package bgzf

import (
	"fmt"
	"io"
)

// SeekReader provides random access to a BGZF-compressed
// stream backed by an *os.File. It supports Seek to a
// BGZF virtual offset (the high 48 bits are the byte
// offset of the block in the underlying file; the low
// 16 bits are the byte offset within the decompressed
// block payload).
//
// A SeekReader is not safe for concurrent use.
type SeekReader struct {
	// file is the underlying compressed stream. We keep
	// the file handle open for the lifetime of the
	// SeekReader so Seek can re-position it cheaply.
	file io.ReadSeekCloser

	// rs is the current BGZF reader. Each Seek creates a
	// new BGZF reader positioned at the right block
	// boundary; the previous one is discarded.
	rs io.Reader

	// blockOff is the byte offset in the underlying file
	// of the block that rs is currently positioned at.
	blockOff int64
}

// NewSeekReader returns a SeekReader positioned at the
// start of the BGZF stream. Callers should defer Close
// to release the file handle.
func NewSeekReader(file io.ReadSeekCloser) (*SeekReader, error) {
	if _, err := file.Seek(0, io.SeekStart); err != nil {
		return nil, fmt.Errorf("bgzf: seek to start: %w", err)
	}
	return &SeekReader{file: file, rs: NewReader(file), blockOff: 0}, nil
}

// Close releases the underlying file.
func (s *SeekReader) Close() error {
	if s.file == nil {
		return nil
	}
	err := s.file.Close()
	s.file = nil
	return err
}

// Seek positions the reader at the given virtual offset.
// Subsequent Read calls return bytes starting from that
// offset.
func (s *SeekReader) Seek(virtualOff uint64) error {
	blockOff := int64(virtualOff >> 16)
	byteOff := int(virtualOff & 0xFFFF)
	if _, err := s.file.Seek(blockOff, io.SeekStart); err != nil {
		return fmt.Errorf("bgzf: seek to block %d: %w", blockOff, err)
	}
	// Discard the previous BGZF reader and create a fresh
	// one at the new block offset.
	s.rs = NewReader(s.file)
	s.blockOff = blockOff
	// Skip byteOff bytes within the new block's decompressed
	// payload. BGZF blocks are at most 64 KiB so this loop
	// always terminates quickly.
	if byteOff > 0 {
		skipped := 0
		for skipped < byteOff {
			todo := byteOff - skipped
			buf := make([]byte, todo)
			n, err := s.rs.Read(buf)
			skipped += n
			if err != nil {
				return fmt.Errorf("bgzf: skip %d bytes within block: %w", byteOff, err)
			}
		}
	}
	return nil
}

// Read decompresses bytes from the current position and
// writes them to p. Returns io.EOF when the BGZF stream
// ends.
func (s *SeekReader) Read(p []byte) (int, error) {
	return s.rs.Read(p)
}