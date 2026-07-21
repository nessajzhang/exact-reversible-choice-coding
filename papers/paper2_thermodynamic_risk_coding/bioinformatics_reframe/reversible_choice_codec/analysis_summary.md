# Reversible Choice Codec Validation

Status: `COMPLETE_REPRODUCIBLE_ASSAY_CALIBRATED_CHOICE_CODEC`

## Exact formal result

The 108-nt GC 49--59 / homopolymer <= 3 language contains exactly `21520867790325216400381593480072139809549357542204871770858223072` sequences and supports a total 213-bit domain (1.972222 bits/nt).
An independent bottom-up recurrence returned the same total (`21520867790325216400381593480072139809549357542204871770858223072`), and an explicit reduced-length enumeration completed 45128/45128 all-rank checks with 0 failures.
Reserving `r` low-order selector bits gives every payload `2^r` exact candidates; the selected sequence reranks to the same physical rank and decoding discards only those `r` bits.  The rate cost is therefore exactly `r` bits per sequence.

## Held-out public experimental libraries

Each direction uses a source-only frozen model and an outcome-blind 2,048-sequence target codebook satisfying the same GC/homopolymer constraints.  One candidate fiber, not one candidate sequence or hash seed, is the paired analysis unit.

- GCall_to_GCfix, r=2: paired experimental gain 0.2503 outcome SD; raw gain 0.00234903, 95% fiber bootstrap [0.00190028, 0.00281840]; 512/512 payloads round-tripped.
- GCall_to_GCfix, r=4: paired experimental gain 0.3398 outcome SD; raw gain 0.00318967, 95% fiber bootstrap [0.00255258, 0.00385161]; 128/128 payloads round-tripped.
- GCfix_to_GCall, r=2: paired experimental gain 0.2744 outcome SD; raw gain 0.00213546, 95% fiber bootstrap [0.00166214, 0.00258502]; 512/512 payloads round-tripped.
- GCfix_to_GCall, r=4: paired experimental gain 0.3812 outcome SD; raw gain 0.00296666, 95% fiber bootstrap [0.00230855, 0.00360908]; 128/128 payloads round-tripped.

These are retrospective exact codebooks assembled from previously assayed sequences. They test the selector on held-out measured outcomes but are not a prospective synthesis or PCR experiment on newly generated codewords.

Across 32 prespecified outcome-blind library/hash mappings:
- GCall_to_GCfix, r=2: P5 standardized gain range [0.1988, 0.2668], positive in 1.000 of mappings; P5-minus-P2 range [0.0124, 0.0993], positive in 1.000.
- GCall_to_GCfix, r=4: P5 standardized gain range [0.2505, 0.3850], positive in 1.000 of mappings; P5-minus-P2 range [0.0338, 0.1951], positive in 1.000.
- GCfix_to_GCall, r=2: P5 standardized gain range [0.2456, 0.3169], positive in 1.000 of mappings; P5-minus-P2 range [0.0406, 0.1189], positive in 1.000.
- GCfix_to_GCall, r=4: P5 standardized gain range [0.3481, 0.4896], positive in 1.000 of mappings; P5-minus-P2 range [0.0400, 0.1903], positive in 1.000.
These seed ranges are algorithmic sensitivity analyses, not additional biological replicates.

## Generated-language rate--utility audit

- GCall selector, r=2: retained payload 211 bits (1.953704 bits/nt); own-model gain 0.00239274; opposite-pool-model gain 0.00250810; selected round-trip failures 0.
- GCall selector, r=4: retained payload 209 bits (1.935185 bits/nt); own-model gain 0.00364235; opposite-pool-model gain 0.00378689; selected round-trip failures 0.
- GCfix selector, r=2: retained payload 211 bits (1.953704 bits/nt); own-model gain 0.00259213; opposite-pool-model gain 0.00231263; selected round-trip failures 0.
- GCfix selector, r=4: retained payload 209 bits (1.935185 bits/nt); own-model gain 0.00395302; opposite-pool-model gain 0.00348049; selected round-trip failures 0.

Threshold rejection had failure fractions up to 0.867188 at the tested caps, whereas fixed-choice encoding remained total on every sampled payload.
The minimum generated-within-source-min/max fraction across active P5 features was 0.997927. 2 low-overlap invariant feature row(s) had a zero P5 coefficient; feature-level overlap is nevertheless reported because generated-sequence scores remain computational predictions.

## Statistical and claim contract

- Public confidence intervals use 2,000 paired bootstrap resamples of candidate fibers.
- Generated-score intervals use payload-fiber resampling and do not represent biological replication.
- Source-only fitting, regularization and thresholds use no target-pool outcomes.
- P0 composition, P2 assay-context and P5 complete-context selectors use identical fibers.
- No result claims mechanism, causal PCR improvement, emitted-sequence wet-lab validation, NUPACK, recovery, full weighted-language capacity or storage-system superiority.

## Reproduction

```bash
$PYTHON analysis_tools/validate_paper2_reversible_choice_codec.py --jobs 12 --bootstrap 2000 --generated-payloads 512
```
