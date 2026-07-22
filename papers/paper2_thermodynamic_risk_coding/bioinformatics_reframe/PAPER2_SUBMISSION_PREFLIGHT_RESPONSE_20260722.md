# Response to final pre-submission review

Manuscript: *Exact reversible choice coding with retrospective assay-guided selection in PCR-tested oligonucleotide libraries*

Status: `TECHNICAL_REVISION_COMPLETE; AUTHOR_AND_ARCHIVE_FIELDS_REQUIRED_BEFORE_SUBMISSION`

We thank the reviewer for the final pre-submission assessment. We implemented every automatable required correction without changing frozen scientific values or deleting adverse results. The Q5 reversal, matched-condition negative fibers, endpoint-dependent released-1D-CNN comparisons and absence of prospective assay for generated candidates remain explicit.

## 1. Release consistency checker

**Addressed.** Compiled-log inspection now has a three-state contract. An omitted log in a clean release is `SKIP`, a present clean log is `PASS` and a present log containing a declared LaTeX error is `FAIL`. Missing logs no longer produce `FileNotFoundError`; compiled PDFs remain mandatory. Three unit tests cover all states, and the expanded full suite contains 21 tests. The release builder now runs the checker before writing the root SHA-256 manifest, so the packaged `SKIP` audit is locked and an idempotent checker rerun does not invalidate the manifest.

**Locations:** `analysis_tools/verify_paper2_bioinformatics_choice_consistency.py`; `tests/test_paper2_consistency_log_contract.py`; `verify_integrity.sh`; `reproduce_analysis.sh`; package `README.md`.

## 2. Q5 reversal in the main figure

**Addressed.** Figure 2d now shows eight positive matched-assay FullContext intervals and four negative secondary exploratory Q5 intervals on a shared zero-effect axis. Values come from the existing grouped-source plus target-fiber two-stage output; none was refitted or changed for plotting. The source-data table contains all 12 rows, and the consistency checker verifies the sign and row count directly.

**Locations:** Figure 2d and caption in `main.tex`; `analysis_tools/plot_paper2_reversible_choice_codec.py`; `figures/reversible_choice_codec_source_data.tsv`; `MANUSCRIPT_CONSISTENCY_AUDIT.md`.

## 3. Abstract estimand and negative fibers

**Addressed.** The 149-word structured abstract now names the primary estimand as measured selected-minus-fiber-mean gain and retains the 17.2--30.9% negative-fiber fraction, endpoint-dependent CNN comparison, exploratory Q5 reversal and unassayed status of generated candidates.

**Location:** structured abstract, Results paragraph in `main.tex`.

## 4. Absolute-effect interpretation

**Addressed.** Results state that absolute FullContext gains were modest (`0.0021--0.0041` on the published relative-efficiency scale) and direct interpretation to standardized effects, low-efficiency event risk and prospective synthesis-level validation.

**Location:** Results, `main.tex`.

## 5. CRC scope

**Addressed.** The main text now states that zero CRC silent misdecodes in the finite deterministic enumeration does not guarantee detection under arbitrary corruptions. The analysis remains labelled an error-detection composition audit, not correction, recovery or a stochastic DNA-channel experiment.

**Locations:** Algorithm, `main.tex`; channel-error section in `supplementary_codec_evidence.tex`; `bioinformatics_reframe/channel_error_boundary/`.

## 6. Edit-distance and k-mer audit

**Not added; scope retained explicitly.** This was recommended rather than required. The manuscript reports only what was actually audited: exact and declared Hamming-near overlap for fixed-length sequences. It explicitly disclaims edit-distance neighbourhood, design-genealogy and cluster-independence analysis. Adding an unplanned metric at this final stage was not necessary to support the bounded Hamming-separation claim.

**Locations:** System and methods and Discussion in `main.tex`; sequence-independence sections in the Supplementary Information.

## 7. Public repository and licence

**Addressed.** The public source repository is https://github.com/nessajzhang/exact-reversible-choice-coding and project code uses the BSD-3-Clause licence. The release contains split integrity/analysis/figure/TeX workflows, frozen derived outputs, source data for every main figure, tests, manifests and a seven-page OUP preflight. The release builder now includes the Figure 1 TSV at extraction time rather than relying on a later figure-reproduction run to create it. Third-party public data and model assets retain their original terms.

## 8. Persistent software archive

**Technical metadata complete; external deposition deferred by the author.** `ZENODO_DEPOSIT_METADATA_20260722.md` supplies an optional Zenodo route. The journal asks for an archived submitted version and test data with a stable URL; a DOI is not the only possible persistent identifier. No deposit, tag or identifier was fabricated or published in this closeout.

## 9. AI provenance and final wording

**Machine-side audit complete; all-author verification remains required.** The provenance record lists the actual assistance with wording, code, tests, literature metadata, figures and release auditing, and states that no experimental observations were generated. A factual disclosure draft is present in the cover letter. Only the authors can confirm whether it covers every generative-AI tool used and approve final wording.

**Locations:** `AI_ASSISTANCE_PROVENANCE_RECORD_20260722.md`; `AUTHOR_LED_FINAL_TEXT_VERIFICATION_20260722.md`; `BIOINFORMATICS_COVER_LETTER_20260722.md`.

## 10. Authorship and declarations

**Author action required.** The review manuscript remains anonymous by design. No author names, order, affiliations, ORCIDs, contact, CRediT roles, funding, conflicts or acknowledgements were inferred. A completion matrix records every missing field.

**Location:** `BIOINFORMATICS_AUTHOR_CONTROLLED_FIELDS_20260722.md`.

## Technical verification

- main review PDF: 9 pages;
- Supplementary Information: 17 pages;
- OUP preflight: 7 pages without font, margin or spacing manipulation;
- structured abstract: 149 words;
- exact generated-candidate round trips: 65,024/65,024;
- test suite: 21/21 pass;
- main/Supplement/OUP compilation with BibTeX: pass;
- rendered-page visual QC: pass;
- clean-release missing-log contract: `SKIP`/`SKIP`, followed by overall consistency pass and an unchanged root manifest after checker rerun;
- root and nested SHA-256 verification: pass.

The technical revision is complete. The manuscript should not be labelled formally submission-ready until the author identity/declaration record and a persistent archive URL or identifier for the manuscript-matched software/test-data snapshot are supplied. Any policy-required disclosure wording must be verified by the authors rather than inferred by software.

## Final-review technical corrections

The abstract's decoder sentence now says that exact inversion recovers the payload from any candidate, rather than saying that it recovers a fiber member. The reproduction bootstrap validates portable Python 3.12.10 and all six pinned distribution versions before permitting canonical-output writes; file-only integrity checks no longer trigger this bootstrap. The release builder disables Python bytecode generation and rejects interpreter, editor and operating-system caches before creating `SHA256SUMS.txt`. Five environment and release-hygiene tests cover these contracts. No tag or external archive was created because that external publication action remains author-controlled.
