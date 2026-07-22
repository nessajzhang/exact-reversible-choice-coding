# Paper 2 Manuscript Consistency Audit

Date: 2026-07-22 (Asia/Shanghai)

Status: `PASS_SUBMISSION_TECHNICAL_OUTPUT_AND_BOUNDARY_CONSISTENCY`

- Visible structured abstract count: 149 words, including headings and placeholders (maximum recommended: 150).
- FullContext within-pool and source-only transfer Spearman values match the frozen assay table.
- Public FullContext choice-fiber gains and exact round trips match both transfer directions at `r=2` and `r=4`.
- Analysis-plan-locked external KAPA gains match both frozen source models and widths on 2,048 outcome-blind sequences drawn from already-public measurements.
- Corrected grouped-source plus target-fiber bootstrap intervals are positive in all eight cross-pool/KAPA FullContext settings; resampled duplicate source rows retain their original fold, equivalent to integer weighting.
- FullContext-minus-released-1D-CNN continuous-efficiency intervals are positive in all eight cross-pool/KAPA comparisons after the corrected two-stage propagation.
- Endpoint-matched FullContext-minus-CNN low-efficiency intervals include zero in all eight cross-pool/KAPA comparisons; no cross-endpoint dominance is asserted.
- The Q5 workflow-shift FullContext, FullContext-minus-AssayContext and FullContext-minus-CNN two-stage intervals are negative in all four source-width settings.
- Figure 2d places eight positive matched-assay and four negative exploratory-Q5 FullContext intervals on one zero-effect axis; its source-data rows are checked directly.
- Cross-pool/KAPA FullContext random-choice tests retain 100,000 draws, zero exceedances and the declared minimum attainable Monte Carlo P formula.
- FullContext and FullContext-minus-AssayContext gains are positive in all 32 outcome-blind public-codebook mappings.
- FullContext gain, FullContext-minus-AssayContext and FullContext-minus-CNN are positive in all 32 post hoc external KAPA mappings for both source models and widths; mappings are not treated as replicates.
- The data-partition/evidence hierarchy, r=0..6 rate-benefit, negative-fiber distributions, feature stratification, score-margin abstention and fixed-precision selector audits are complete and manifest-locked.
- Fixed-precision score keys changed no AssayContext, FullContext or released-CNN primary choice; exploratory Composition differences are disclosed.
- Exact 108-nt total, modulo-N rank dispersion, 213-bit domain, 45,128 reduced all-rank checks and 65,024 generated candidate round trips are internally consistent.
- Generated first-base and all-position composition audits exclude the former lexicographic-prefix bias; eight generated mapping namespaces remain positive with zero round-trip failures.
- Exhaustive 143,904,012-pair fixed-length Hamming audit retains zero exact duplicates or reverse complements, minimum distance 54/108 and no pair at 70% identity; no edit-distance, genealogy or cluster-independence claim is made.
- Exhaustive source-external audits retain zero duplicates, reverse complements or pairs at 70% identity across both the 2,048 lock and all 2,053 eligible sequences.
- Threshold-rejection failure fractions and zero choice-codec failures match the frozen table.
- Runtime prose and SI table match the current fixed-input, machine-specific benchmark; all timed round trips passed.
- The deterministic single-edit audit covers 55,552 substitutions, insertions and deletions over 64 emitted codewords; 12,093 wrong-payload acceptances without a checksum fell to zero after CRC-16/CCITT-FALSE in this finite composition audit, which is not a channel or recovery claim.
- Reported intervals condition on already-public sequence-level measurements and do not propagate technical-replicate, batch, normalization or assay-measurement uncertainty.
- All 19 active citation keys exist in the paper-local BibTeX file and match the separately verified reference set.
- Required public/generated, computational/material and observational/causal boundaries are present; forbidden positive-claim patterns are absent.
- LaTeX compilation-log audit: main.log=SKIP (log absent from release archive; compiled-log checks skipped); supplementary_codec_evidence.log=SKIP (log absent from release archive; compiled-log checks skipped). Compiled PDFs remain required even when release-excluded logs are skipped.
- The compiled manuscript, seven-page OUP preflight and Supplementary pages are covered by the rendered-page audit in `PDF_VISUAL_QC_20260722.md`.

## Verified manifests

- `public_experimental_validation`: 10 entries; manifest SHA-256 `35fec335cb2a0d28b7a940b91ccbf96c05f67b7b465e6d63edd956c3c4139e3d`.
- `assay_calibrated_selection`: 9 entries; manifest SHA-256 `c11ec35c45188162b621d65881465d378ddf26d283edda2e03314249f5efb2fc`.
- `reversible_choice_codec`: 20 entries; manifest SHA-256 `d9d7faa34165feddc46657d1091ea347afec8436fa9d8ecf303d475f24a5c455`.
- `sequence_independence_audit`: 4 entries; manifest SHA-256 `d55d5c5ec5739de0527a10807968cd98b2efa91284c7857b00541483b84ccf19`.
- `runtime_benchmark`: 4 entries; manifest SHA-256 `468f591c79e12c891f98cf4dd303a70413ecd121842113090f8272a9b3a2b7ea`.
- `sota_and_external_validation`: 32 entries; manifest SHA-256 `c6d89b8463fd15afb1d4c7ed4437e408e81da93c37d8990f4a44b27660ddbe0c`.
- `source_external_independence_audit`: 5 entries; manifest SHA-256 `6f2ae722e19adf8b1f6191a38b9e5df521a39d1c99109f37671eabca3dafc28b`.
- `external_mapping_sensitivity`: 6 entries; manifest SHA-256 `0ccc3ff132f02681aa20b6c672424f6a6603de5d74198eb72216a38356ae45fe`.
- `major_revision_diagnostics`: 13 entries; manifest SHA-256 `63b70de86b2a9250ef864d9301110cb3ec8056f26c04a11951d87dc6a14c5cec`.
- `channel_error_boundary`: 7 entries; manifest SHA-256 `53612d8081a7fcf944b67c8557502d3d54368621eb4c88493016a4ce53938c1e`.

## Current artifact hashes

- `main.tex`: `5cae1a6370ea76df9a56e98b6f88c75abfc5c46add098bdbadf4cf4b90856daa`
- `supplementary_codec_evidence.tex`: `eaa3b16b42b8eb31b7bc13a12e5ef3dbf7c40c47a96b7c4f2e5eea2d77f49869`
- `build/main.pdf`: `74b8c746497c96ae37cddae52dedf08a0867dd5393b4c0817d3e0de6961e4c5c`
- `build/supplementary_codec_evidence.pdf`: `7f607d790c5d7da16dbbaf5a615988e7a2816a6e3cc909115857c1a484854873`
- `figures/assay_calibrated_selection.pdf`: `edaec0fbd4c23809066738a2e70e59bcd80c687b75e7ade20f30a8a4ccf655b5`
- `figures/reversible_choice_codec.pdf`: `4fc98d1ac184d25ff791a0d128aa80509dbdacd5b68f8a21224d39e60f004168`
- `figures/paired_two_stage_sensitivity.pdf`: `3f7e5563dc601e6bbfb4fd4e3f98723113decf5655c328ec27c0021183cb29ac`

## Submission boundary

Scientific and technical consistency is a narrow PASS, not an acceptance guarantee. The public repository and BSD-3-Clause licence are machine-verifiable. The immutable archive DOI remains blocked by Zenodo sign-in and author-approved creator metadata; journal upload also requires author-controlled identity/contact/funding/conflict/CRediT fields and all-author verification of the actual AI-use disclosure and final wording. The generated codec still has no prospective emitted-codeword wet-lab validation.
