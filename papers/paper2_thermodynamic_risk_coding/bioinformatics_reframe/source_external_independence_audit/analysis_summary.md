# Source--external fixed-length Hamming-separation audit

Status: `COMPLETE_EXHAUSTIVE_SOURCE_EXTERNAL_AUDIT`

- GCall vs external_codec_eligible_2053: 24,631,894 exact pair comparisons; minimum Hamming distance 55/108 (maximum identity 49.074%); duplicates=0, reverse complements=0, pairs at >=70% identity=0.
- GCall vs external_locked_codebook_2048: 24,571,904 exact pair comparisons; minimum Hamming distance 55/108 (maximum identity 49.074%); duplicates=0, reverse complements=0, pairs at >=70% identity=0.
- GCfix vs external_codec_eligible_2053: 24,623,682 exact pair comparisons; minimum Hamming distance 56/108 (maximum identity 48.148%); duplicates=0, reverse complements=0, pairs at >=70% identity=0.
- GCfix vs external_locked_codebook_2048: 24,563,712 exact pair comparisons; minimum Hamming distance 56/108 (maximum identity 48.148%); duplicates=0, reverse complements=0, pairs at >=70% identity=0.

The fixed 0F and 0R-prime sequences are assay context, not variable content. This audit establishes only the reported fixed-length Hamming separation; it does not test edit-distance neighbourhoods, design genealogy or cluster-level independence and does not create biological replication or prospective codec-output evidence.
