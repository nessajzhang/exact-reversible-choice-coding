# Assay-Calibrated Selection Utility

Status: `COMPLETE_REPRODUCIBLE_SELECTION_AUDIT`

## Question

At an identical retained fraction of each public sequence pool, does training-fold-only calibration select sequences with more favourable measured outcomes than a frozen low-weighted-pair rule or a matched structural baseline?

## Primary 25% retention comparisons

- DT4DDS / repeated_nested_oof / Genscript_GCall / day7_vs_day0_log2fc_plus1 / calibration_vs_frozen_proxy: delta selected mean = 0.00325081, 95% paired bootstrap [-0.00941718, 0.0157497]; delta adverse-event rate = -0.000641437.
- DT4DDS / repeated_nested_oof / Genscript_GCall / day7_vs_day0_log2fc_plus1 / weighted_increment: delta selected mean = -0.000570526, 95% paired bootstrap [-0.00269522, 0.00148819]; delta adverse-event rate = 0.
- DT4DDS / repeated_nested_oof / Genscript_GCall / mean_log2_cpm_plus1 / calibration_vs_frozen_proxy: delta selected mean = -0.00494925, 95% paired bootstrap [-0.0510799, 0.0407988]; delta adverse-event rate = 0.00224503.
- DT4DDS / repeated_nested_oof / Genscript_GCall / mean_log2_cpm_plus1 / weighted_increment: delta selected mean = 0.00502479, 95% paired bootstrap [-0.00985122, 0.0198647]; delta adverse-event rate = 0.
- DT4DDS / repeated_nested_oof / Twist_GCall / day7_vs_day0_log2fc_plus1 / calibration_vs_frozen_proxy: delta selected mean = 0.000263568, 95% paired bootstrap [-0.00578611, 0.00678277]; delta adverse-event rate = -0.00133333.
- DT4DDS / repeated_nested_oof / Twist_GCall / day7_vs_day0_log2fc_plus1 / weighted_increment: delta selected mean = 0.00010369, 95% paired bootstrap [-0.0010965, 0.00132365]; delta adverse-event rate = 0.
- DT4DDS / repeated_nested_oof / Twist_GCall / mean_log2_cpm_plus1 / calibration_vs_frozen_proxy: delta selected mean = -0.00553878, 95% paired bootstrap [-0.0240421, 0.0128913]; delta adverse-event rate = 0.000666667.
- DT4DDS / repeated_nested_oof / Twist_GCall / mean_log2_cpm_plus1 / weighted_increment: delta selected mean = -0.00277003, 95% paired bootstrap [-0.00885769, 0.00229747]; delta adverse-event rate = 0.000333333.
- Gimpel2025_PCR / repeated_nested_oof / GCall / relative_efficiency / calibration_vs_frozen_proxy: delta selected mean = 0.00310245, 95% paired bootstrap [0.00281135, 0.00342994]; delta adverse-event rate = -0.0196667.
- Gimpel2025_PCR / repeated_nested_oof / GCall / relative_efficiency / full_structural_context_increment: delta selected mean = 0.000248054, 95% paired bootstrap [6.56579e-05, 0.000422875]; delta adverse-event rate = -0.000333333.
- Gimpel2025_PCR / repeated_nested_oof / GCall / relative_efficiency / full_weighted_increment: delta selected mean = 0.000517082, 95% paired bootstrap [0.00036414, 0.000695393]; delta adverse-event rate = -0.00466667.
- Gimpel2025_PCR / repeated_nested_oof / GCall / relative_efficiency / weighted_increment: delta selected mean = -4.15099e-06, 95% paired bootstrap [-2.12795e-05, 1.09521e-05]; delta adverse-event rate = 0.
- Gimpel2025_PCR / repeated_nested_oof / GCfix / relative_efficiency / calibration_vs_frozen_proxy: delta selected mean = 0.00296237, 95% paired bootstrap [0.00261933, 0.00332214]; delta adverse-event rate = -0.0226742.
- Gimpel2025_PCR / repeated_nested_oof / GCfix / relative_efficiency / full_structural_context_increment: delta selected mean = 0.000172653, 95% paired bootstrap [4.65384e-05, 0.000307302]; delta adverse-event rate = -0.000333444.
- Gimpel2025_PCR / repeated_nested_oof / GCfix / relative_efficiency / full_weighted_increment: delta selected mean = 0.000405451, 95% paired bootstrap [0.000297197, 0.00052689]; delta adverse-event rate = -0.00166722.
- Gimpel2025_PCR / repeated_nested_oof / GCfix / relative_efficiency / weighted_increment: delta selected mean = 1.74667e-05, 95% paired bootstrap [-7.23494e-07, 3.53741e-05]; delta adverse-event rate = 0.
- Gimpel2025_PCR / source_only_transfer / GCall_to_GCfix / relative_efficiency / calibration_vs_frozen_proxy: delta selected mean = 0.00295584, 95% paired bootstrap [0.00261674, 0.00329983]; delta adverse-event rate = -0.0233411.
- Gimpel2025_PCR / source_only_transfer / GCall_to_GCfix / relative_efficiency / full_structural_context_increment: delta selected mean = 0.000152259, 95% paired bootstrap [2.7369e-05, 0.000276491]; delta adverse-event rate = -0.000333444.
- Gimpel2025_PCR / source_only_transfer / GCall_to_GCfix / relative_efficiency / full_weighted_increment: delta selected mean = 0.000423329, 95% paired bootstrap [0.000308651, 0.000541973]; delta adverse-event rate = -0.00200067.
- Gimpel2025_PCR / source_only_transfer / GCall_to_GCfix / relative_efficiency / weighted_increment: delta selected mean = 1.19446e-05, 95% paired bootstrap [-3.09154e-06, 2.84438e-05]; delta adverse-event rate = 0.
- Gimpel2025_PCR / source_only_transfer / GCfix_to_GCall / relative_efficiency / calibration_vs_frozen_proxy: delta selected mean = 0.00270252, 95% paired bootstrap [0.00239832, 0.00300073]; delta adverse-event rate = -0.0203333.
- Gimpel2025_PCR / source_only_transfer / GCfix_to_GCall / relative_efficiency / full_structural_context_increment: delta selected mean = 0.000293842, 95% paired bootstrap [0.00015174, 0.000431019]; delta adverse-event rate = -0.001.
- Gimpel2025_PCR / source_only_transfer / GCfix_to_GCall / relative_efficiency / full_weighted_increment: delta selected mean = 0.000399199, 95% paired bootstrap [0.000239492, 0.000560026]; delta adverse-event rate = -0.004.
- Gimpel2025_PCR / source_only_transfer / GCfix_to_GCall / relative_efficiency / weighted_increment: delta selected mean = -2.20373e-05, 95% paired bootstrap [-6.84644e-05, 1.03338e-05]; delta adverse-event rate = 0.000333333.

## Interpretation contract

- Positive within-pool selection utility supports assay-specific retrospective calibration, not a universal molecular-risk direction.
- A null weighted increment means the weighted-pair scalar adds no utility beyond candidate count, longest stem and declared assay context in that comparison.
- Source-only transfer is reported separately and determines whether a calibration can move between sequence pools without target tuning.
- The finite-pool retention penalty is computed as -log2(selected_n/N) from the achieved retained fraction and is reported per nucleotide only as finite-pool accounting. It is not full-language capacity or achieved storage density.
- Public sequences were not emitted by the Paper 2 codec. No synthesis, PCR, sequencing or recovery experiment on Paper 2 oligos is claimed.

## Statistical contract

- The independent analysis unit is one public reference sequence.
- Within-pool predictions are means over five repeated five-fold outer passes; scaling and ridge tuning occur inside training folds.
- Transfer models are fitted and tuned on the complete source pool without target-pool outcomes.
- Intervals use 2,000 paired sequence-level bootstrap resamples conditional on the frozen predictions and fixed selected sets.
- Retention fractions of 10%, 25% and 50% are all retained in the output; 25% is the primary descriptive operating point.
