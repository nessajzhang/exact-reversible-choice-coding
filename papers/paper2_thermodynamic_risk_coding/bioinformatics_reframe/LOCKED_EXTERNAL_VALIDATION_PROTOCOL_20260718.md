# Locked external-laboratory validation protocol

Lock date: 2026-07-18 (Asia/Shanghai)

Purpose: add a genuinely separate, external-laboratory test of frozen
GCall/GCfix-trained selectors without changing the codec, feature sets,
candidate budget or primary effect definition after inspecting selector
performance on this pool.

This is a secondary validation. It does not convert the public sequences into
prospectively emitted Paper 2 codewords.

## Frozen public input

- Publication and repository: Gimpel et al. 2025, DOI
  `10.1038/s41467-025-64221-4`; repository commit
  `af62c57f9a90ecdfdd0f1623441e82bdb7e082c1`.
- Sequence/prediction table:
  `external-validation-predictions.csv`, SHA-256
  `2195ede4ce1a7172840c3091553e73a394c777a1d61dbad292b0ba556b984fb2`.
- Independent sequence table:
  `sequence_data_anonymized_no_duplicates.csv`, SHA-256
  `53e9e4bdd4b740061b1650d3b4c8adcb3012312af20bbfbe5b75748865017864`.
- The two tables must join one-to-one on `seq_id_anonymized` with identical
  108-nt sequence strings.

## Cohort and outcome

1. Include only rows with `has_insertedmotif == False` (the 10,000 originally
   random sequences).
2. Require a unique 108-nt A/C/G/T sequence and a non-missing `eff_Taq` value.
3. `eff_Taq` is the primary measured endpoint because it is the external-
   laboratory workflow using the GCall/GCfix KAPA/Taq conditions. `eff_Q5` is
   a workflow-shift sensitivity endpoint and cannot support same-assay claims.
4. A descriptive low-efficiency event is defined by the article repository's
   empirical bottom-2% procedure on eligible random sequences. Continuous
   efficiency remains the primary endpoint.
5. Missing-outcome rows are excluded before outcome-blind library selection;
   no imputation is allowed.

## Frozen candidate library and fibers

- Apply the unchanged base-code eligibility rule: length 108 nt, whole-sequence
  GC count 49--59 inclusive, and maximum homopolymer length at most 3.
- Sort eligible sequences by
  `sha256("paper2-choice-external-v1|external_Taq|" + sequence_sha256)`, then by
  `sequence_sha256`; take the first 2,048 sequences.
- Assign base ranks 0--2,047 in that order and contiguous non-overlapping
  fibers of width `2**r`.
- `r=2` is primary; `r=4` is secondary. `r=0` is the no-choice reference. No
  alternative key, library size, fiber permutation or selector-bit value may
  replace these after results are computed.
- The same library and fibers are used for every selector.

## Frozen selectors

1. `P0_composition`, `P2_assay_context` and `P5_combined_context` are fitted
   separately on all GCall and all GCfix source sequences using the existing
   fixed feature definitions, scaling, ridge-alpha grid, five-fold source-only
   selection and seeds. External outcomes are never used in fitting, scaling,
   alpha selection or feature construction.
2. Both GCall-trained and GCfix-trained models are reported. No post-result
   winner selection or outcome-tuned ensemble is permitted.
3. Published `GCall 2perc` and `GCfix 2perc` 1D-CNN probabilities are external
   SOTA comparators. The selected candidate minimizes the published probability
   of the positive (low-efficiency) class; ties use the smallest base rank.
4. The no-choice expectation is the within-fiber mean. A deterministic
   first-candidate selector is reported only as a null/tie-break check, not as
   a biological method.

## Frozen estimands and inference

For each selector and `r`, report:

- selected-minus-fiber-mean measured efficiency for every fiber;
- mean, 95% interval, median, IQR and fractions below/equal to/above zero;
- standardized mean gain using the fixed 2,048-sequence library population SD;
- selected low-efficiency event rate and paired event-rate change;
- payload count, candidate width, retained payload bits and exact round-trip
  audit status.

For P5, use a two-stage bootstrap that resamples the source sequences, repeats
source scaling/alpha selection/model fitting, and resamples target fibers. The
published CNN outputs are fixed upstream predictions, so their uncertainty
interval is conditional on those predictions and resamples target fibers only.
Use a within-fiber random-choice randomization test for each directional mean
gain. Mapping-seed sensitivity is algorithmic sensitivity, not replication and
must not narrow biological intervals.

## Interpretation boundary

A positive result supports transfer of a frozen selector to measured sequences
generated and amplified in an external laboratory under matched assay
conditions. It does not establish a prospective codec-output experiment,
mechanism, causal effect, NUPACK validity, material storage recovery, complete
110-nt language capacity or superiority over end-to-end DNA-storage systems.

