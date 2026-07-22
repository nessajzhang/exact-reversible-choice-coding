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
- `external_mapping_sensitivity`: 6 entries; manifest SHA-256 `2f7ba4e956f3bf02de5cc1b9c49bf9edf0549ebb5190141574846f0ed8efa207`.
- `major_revision_diagnostics`: 13 entries; manifest SHA-256 `5d376c268070fe233f28205cbe5f7adbb83c2fdae276292618735225435d2839`.
- `channel_error_boundary`: 7 entries; manifest SHA-256 `53612d8081a7fcf944b67c8557502d3d54368621eb4c88493016a4ce53938c1e`.

## Current artifact hashes

- `main.tex`: `db06c885d761468bbff885e0db09caf85437dd4ac5481dad90a379f76de37b26`
- `supplementary_codec_evidence.tex`: `d9c027abf38eea4bed118ffb4c13573ed5e1684f21f3a3d614ed09602e4450de`
- `build/main.pdf`: `095ba48e4b0d6d86560af4f679dd43743fd7c98a2f18b3d9f04292cd28025dbb`
- `build/supplementary_codec_evidence.pdf`: `0246c7586b568e85304d713d0c8d47e9891027e3edc83de2807a2bd75f3b50be`
- `figures/assay_calibrated_selection.pdf`: `75156a3361366d095a34402670c6437f72321a6a2a4660c202981e842996898e`
- `figures/reversible_choice_codec.pdf`: `fb2efd7598e8289876f606a767e9b51684228a8994adb3cbafe187652ca4a6d4`
- `figures/paired_two_stage_sensitivity.pdf`: `887f5b35198844a4dd481658fb9ae576c28d0c698fd36ae9137f8694d478700e`

## Submission boundary

Scientific and technical consistency is a narrow PASS, not an acceptance guarantee. The public repository and BSD-3-Clause licence are machine-verifiable. A persistent archive URL or identifier for the manuscript-matched software and test-data snapshot remains outstanding because the author deferred that external publication action; a DOI is one possible identifier, not the only route. Journal upload also requires author-controlled identity/contact/funding/conflict/CRediT fields and author verification of any policy-required declarations and final wording. The generated codec still has no prospective emitted-codeword wet-lab validation.
