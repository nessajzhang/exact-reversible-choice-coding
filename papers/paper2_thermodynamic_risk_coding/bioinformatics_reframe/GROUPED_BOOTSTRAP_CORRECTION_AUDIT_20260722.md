# Grouped-bootstrap correction audit

Status: `PASS_CORRECTED_GROUPED_SOURCE_BOOTSTRAP`

## Corrected contract

- Unit sampled at the source stage: one original source-sequence ID.
- Bootstrap multiplicity: explicit duplicate rows, equivalent to integer frequency weights for the ridge objective.
- Cross-validation: five folds are assigned once to original source IDs; every bootstrap copy of an ID remains in that fold.
- Cross-fold duplicate leakage: zero by construction and checked in tests.
- Target stage: one fixed candidate fiber.
- Replicates: 2,000 source-plus-target replicates for each of 48 estimands.
- Reference check: fast grouped ridge predictions agree with a `PredefinedSplit` implementation to below `2e-10` in the frozen test.

The complete run was:

```text
$PYTHON analysis_tools/validate_paper2_sota_external_uncertainty.py --jobs 10 --bootstrap 2000 --two-stage 2000 --randomizations 100000
```

## Old-versus-corrected comparison

Relative to the 2026-07-20 baseline, the maximum absolute changes were:

| Quantity | Maximum absolute change |
|---|---:|
| Bootstrap mean | 0.000241376247 |
| Lower 2.5% endpoint | 0.000158108787 |
| Upper 97.5% endpoint | 0.000368791848 |
| Fraction above zero | 0.0160 |

No estimand changed its interval sign classification: all previously positive matched-condition intervals remained positive, all reported negative Q5 FullContext/paired intervals remained negative, and the Q5 AssayContext intervals that crossed zero continued to cross zero. The corrected TSV, all 2,000-replicate rows, environment contract, regenerated figure source data and figures replace the baseline values.

## Inference boundary

The correction addresses source-bootstrap grouping and weighting. The intervals remain conditional on published sequence-level outcome estimates and on the final feature families, model class, endpoints, filters, codebook and mapping plans. They do not propagate technical-replicate, batch, normalization or source-assay measurement uncertainty, nor the complete earlier exploratory model-selection process.
