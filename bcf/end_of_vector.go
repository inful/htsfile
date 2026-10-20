package bcf

// Per-type and END_OF_VECTOR sentinels for the BCF2 mixed-
// length vector encoding (§BcfTypeEncoding, VCFv4.5 §5).
//
// When a variable-length vector is padded out, real values
// are distinguished from padding by exact match against the
// appropriate END_OF_VECTOR sentinel. When a single slot
// within a typed vector is "missing", the writer emits the
// per-type MISSING sentinel (signaling NaN for floats, the
// matching MSB-set value for integers); the reader detects
// these sentinels and surfaces them as Missing=true on the
// corresponding sample.
//
// Eight values per integer width just above the matching
// MISSING sentinel are reserved by the BCF2 encoding
// (line 2306).
const (
	// bcfEndOfVectorInt8 is the END_OF_VECTOR sentinel
	// for int8 GT vectors. The int8 width is the one
	// actually emitted today (line 2463) since int8
	// covers up to ~63 alleles.
	bcfEndOfVectorInt8 uint8 = 0x81

	// bcfEndOfVectorInt32 is the END_OF_VECTOR sentinel
	// for int32 FORMAT integer vectors.
	bcfEndOfVectorInt32 uint32 = 0x80000001

	// bcfEndOfVectorFloat32 is the END_OF_VECTOR sentinel
	// for float32 FORMAT float vectors: signaling NaN
	// (line 2325–2326).
	bcfEndOfVectorFloat32 uint32 = 0x7F800002

	// bcfFloatMissing is the per-type MISSING sentinel for
	// float32 values (line 2333). It is the signaling NaN
	// 0x7F800001 — distinguished from the END_OF_VECTOR
	// sentinel 0x7F800002 by the last two payload bits.
	bcfFloatMissing uint32 = 0x7F800001

	// bcfIntMissingInt32 is the per-type MISSING sentinel
	// for int32 values (line 2304). It is the first of
	// the eight reserved values just above MinInt32.
	// T4.22 already forbids the eight reserved values
	// from the encoder; this constant mirrors the spec's
	// per-type MISSING form on the reader side.
	bcfIntMissingInt32 uint32 = 0x80000000
)
