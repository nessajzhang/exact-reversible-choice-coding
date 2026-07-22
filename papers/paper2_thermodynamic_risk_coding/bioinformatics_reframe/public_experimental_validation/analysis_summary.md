# Public Experimental Data Validation

Status: `COMPLETE_REPRODUCIBLE_PUBLIC_DATA_AUDIT`
Route decision: `PUBLIC_DATA_GATE_NARROWS_TO_VALIDATION_AND_DOMAIN_MISMATCH`

## Data actually analysed

- DT4DDS: 24,472 public sequence-level records, analysed separately as 12,472 Genscript 102-nt references and 12,000 Twist 108-nt references across 0a, 0b, 2d, 4d and 7d measurements.
- Gimpel 2025 PCR: 23,992 filtered public 108-nt variable sequences (GCall 11,998; GCfix 11,994), with 240 author-defined worst-2% labels in each pool. The publication design size was 12,000 per pool, but two GCall and six GCfix references were excluded by the authors' processing rule and are not silently restored here.
- Source-data cross-check: all 11,998 GCall `eff` and `x0` values matched the independent publication Source Data member exactly.
- Full PCR amplicons were reconstructed as publication-defined 0F (20 nt) + variable region (108 nt) + 0R-prime (21 nt), total 149 nt.

## Prespecified primary comparisons

DT4DDS primary metric is the change in out-of-fold Spearman correlation for `D3_full - D1_structural_rules`:

### mean_log2_cpm_plus1: `no incremental association`

- Genscript_GCall: Δ=-0.000325, 95% paired bootstrap [-0.001464, 0.000879]
  Standardized weighted-pair coefficient: median -0.0557053; negative in 25/25 outer fits.
- Twist_GCall: Δ=-0.001381, 95% paired bootstrap [-0.003593, 0.000810]
  Standardized weighted-pair coefficient: median -0.00536627; negative in 25/25 outer fits.

### day7_vs_day0_log2fc_plus1: `not replicated`

- Genscript_GCall: Δ=0.000092, 95% paired bootstrap [-0.000413, 0.000599]
  Standardized weighted-pair coefficient: median -0.0112965; negative in 25/25 outer fits.
- Twist_GCall: Δ=-0.000430, 95% paired bootstrap [-0.001015, 0.000167]
  Standardized weighted-pair coefficient: median -0.00715958; negative in 21/25 outer fits.

### PCR worst-2% labels: `pool-specific performance with risk-direction discordance`

Primary metric is average precision; prevalence is approximately 2% in each processed pool.

- GCall: Δ=0.010861, 95% paired bootstrap [-0.001343, 0.030219]
  Standardized variable-region weighted-pair coefficient: median -0.000948473; negative in 25/25 outer fits (low-efficiency label is 1, so the intended risk direction was positive).
- GCfix: Δ=0.001284, 95% paired bootstrap [-0.017994, 0.018448]
  Standardized variable-region weighted-pair coefficient: median -0.0555961; negative in 25/25 outer fits (low-efficiency label is 1, so the intended risk direction was positive).
- GCall_to_GCfix: Δ=0.002292, 95% paired bootstrap [-0.001966, 0.006357]
- GCfix_to_GCall: Δ=-0.001988, 95% paired bootstrap [-0.007605, 0.002762]

### PCR continuous efficiency: `replicated predictive increment with risk-direction discordance`

The same `P3-P2` comparison was prespecified for the authors' continuous relative-efficiency estimate. Higher efficiency is favorable, so an adverse risk proxy was expected to have a negative standardized coefficient.

- GCall: Δ=0.011278, 95% paired bootstrap [0.007978, 0.014704]
  Standardized variable-region weighted-pair coefficient: median 0.000344073; negative in 0/25 outer fits (prespecified adverse direction: negative).
- GCfix: Δ=0.030236, 95% paired bootstrap [0.024896, 0.035700]
  Standardized variable-region weighted-pair coefficient: median 0.000566213; negative in 0/25 outer fits (prespecified adverse direction: negative).

The complete baseline (`P1-P0`), full-amplicon sensitivity, continuous-efficiency, binary-DT4DDS gate and secondary metrics are retained in the TSV outputs whether favorable or unfavorable.

## Event gates

- DT4DDS / Genscript_GCall / missing_any: 357/12472 events; `nested_5x5_secondary`.
- DT4DDS / Genscript_GCall / missing_day7: 100/12472 events; `nested_5x5_secondary`.
- DT4DDS / Twist_GCall / missing_any: 28/12000 events; `descriptive_only_below_30_events`.
- DT4DDS / Twist_GCall / missing_day7: 11/12000 events; `descriptive_only_below_30_events`.
- Gimpel2025_PCR / GCall / worst_2_percent_label: 240/11998 events; `nested_5x5_confirmatory`.
- Gimpel2025_PCR / GCfix / worst_2_percent_label: 240/11994 events; `nested_5x5_confirmatory`.

## Allowed interpretation

The results quantify whether the frozen Paper 2 string-level proxy carries incremental information for public, sequence-level experimental endpoints after prespecified baselines. For PCR continuous efficiency, the increment is reproducible but its coefficient is opposite to the intended risk direction; this is evidence of domain mismatch, not molecular-risk validation. Pool heterogeneity, null results and confidence intervals are part of the result.

## Forbidden interpretation

These public-data analyses do not show that Paper 2 emitted oligos were synthesized, aged, amplified or recovered; do not isolate a causal molecular mechanism; do not validate NUPACK; do not establish storage-density or full 110-nt language capacity; and do not import evidence from an unpublished sister paper. DT4DDS read counts are composite post-workflow detectability measurements, not isolated aging or synthesis outcomes. Gimpel labels concern their published random pools and adapter/PCR workflow, not Paper 2 codewords.

## Reproducibility and manuscript freeze

- Nested preprocessing and regularization tuning were confined to training folds.
- Every primary delta uses identical outer predictions and 2,000 paired sequence-level bootstrap resamples.
- Source hashes, release commits, exact adapter provenance, environment, all model metrics, coefficients and output hashes are saved beside this report.
- `main.tex` remained frozen at SHA-256 `db06c885d761468bbff885e0db09caf85437dd4ac5481dad90a379f76de37b26` during this run.
