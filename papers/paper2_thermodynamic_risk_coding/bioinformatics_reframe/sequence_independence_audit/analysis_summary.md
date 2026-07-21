# PCR fixed-length Hamming-separation audit

Status: `COMPLETE_EXHAUSTIVE_CROSS_POOL_AUDIT`

All `143,904,012` oriented GCall--GCfix variable-region pairs were compared exactly. The closest pair differed at `54` of `108` positions (maximum identity `50.000%`). There were `0` exact duplicates, `0` exact cross-pool reverse complements, and `0` pairs with at least 90% identity.

The sequence-level analysis unit is the 108-nt synthetic variable region. The 41 shared adapter nucleotides are declared assay context and were excluded. This audit establishes the reported fixed-length Hamming separation only; it does not test edit-distance neighbourhoods, design genealogy or cluster-level independence and does not create biological replication or prospective validation.
