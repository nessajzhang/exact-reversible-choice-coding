# Paper 2 channel-error boundary audit

The deterministic audit covered 64 FullContext-selected generated codewords and all 55,552 single-base edit operations.
The noiseless payload decoder accepted 12,093 edits with a wrong payload; CRC-16/CCITT-FALSE reduced this to 0, and CRC plus the optional source-specific canonical verifier left 0 in this sample.

Insertions and deletions are rejected by the fixed-length interface. This analysis is an error-detection composition audit, not a DNA-channel model, correction experiment, dropout benchmark, consensus analysis or file-recovery result.
