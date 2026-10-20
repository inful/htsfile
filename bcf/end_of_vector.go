package bcf

// END_OF_VECTOR sentinels for the BCF2 mixed-length
// vector encoding (§BcfTypeEncoding, VCFv4.5 §5
// line 2304–2305 for integers and line 2325–2326 for
// floats). When a variable-length vector is padded out
// with END_OF_VECTOR values, real values are
// distinguished from padding by exact match against the
// appropriate sentinel.
//
// The eight values per integer width just above the
// matching MISSING sentinel are reserved by the BCF2
// encoding (line 2306).

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
	// for float32 FORMAT float vectors: signaling NaN.
	bcfEndOfVectorFloat32 uint32 = 0x7F800002
)
