# Exact reversible choice coding pre-submission completion report

Date: 2026-07-22

Overall status: `SCIENTIFIC_AND_TECHNICAL_CLOSEOUT_PASS; PUBLIC_REPOSITORY_AND_LICENSE_LIVE; DOI_AND_HUMAN_AUTHOR_FIELDS_PENDING`

## Required work completed

### Grouped/weighted source bootstrap

- The two-stage source bootstrap now samples original source-sequence IDs with replacement and represents multiplicity as duplicate-row integer weights.
- Every copy of an original ID remains in its fixed validation fold through `PredefinedSplit`; duplicate copies cannot enter both training and validation.
- The fast implementation is tested against a scikit-learn reference implementation.
- The complete 2,000-replicate grouped-source plus target-fiber bootstrap was rerun.
- Across 48 estimands, the maximum absolute change in an interval endpoint was `0.000368791848`; no interval changed sign class.
- Evidence: `GROUPED_BOOTSTRAP_CORRECTION_AUDIT_20260722.md`, `sota_and_external_validation/two_stage_bootstrap_replicates.tsv`, `tests/test_paper2_grouped_bootstrap.py`.

### QC path and rendered-page audit

- Active manuscript/release checks use `bioinformatics_reframe/PDF_VISUAL_QC_20260722.md`.
- The review manuscript, Supplementary Information and OUP preflight were rendered and inspected; final page counts are 9, 17 and 7, respectively.
- No visible clipping or overlap was found; the OUP preflight meets the seven-page limit without changing template font, margins or spacing.
- Figure 2d now places the eight positive matched-assay intervals and four negative exploratory-Q5 intervals on one zero-effect axis; the underlying 12 source-data rows are checked directly by the manuscript consistency audit.

### Measurement and public-data boundaries

- Main text, captions and Supplementary Information state that intervals condition on published sequence-level outcomes and do not propagate technical-replicate, batch, normalization or assay-measurement uncertainty.
- External KAPA is described as an analysis-plan-locked evaluation of already-public measurements from an external laboratory reported in the same source publication.
- The text explicitly disclaims prospective preregistration, third-party replication and prospective assay of codec-emitted sequences.

### Current prior art

- The manuscript and prior-art matrix now include Explorer, HEDGES, Composite Hedges, FrameD, the 2026 standardized codec benchmark and Gungnir, alongside enumerative coding, constrained coding, shaping/list selection and sequence-design tools.
- The novelty claim is limited to the combined fixed-rate, full-domain, scorer-replaceable and scorer-independent-decoding interface plus identical-fiber retrospective PCR evaluation.
- Error-correction novelty is explicitly disclaimed.
- Evidence: `REFERENCE_VERIFICATION_REPORT_20260722.md` and `references.bib`.

### Public code and licence

- Public repository: https://github.com/nessajzhang/exact-reversible-choice-coding
- Licence detected by GitHub: BSD-3-Clause.
- The repository contains manifested code, tests, frozen derived outputs, source-data tables, figures, review PDFs and a seven-page OUP preflight. The final repository commit is the public `main` head produced by the release synchronization; no self-referential commit hash is embedded inside the commit itself.
- A clean-tree audit passed all root and nested manifests, 15 deterministic/grouped-bootstrap/channel/log-contract tests, frozen-data analysis reproduction and all three figures.
- The consistency checker now gives `PASS`, `SKIP` or `FAIL` for compiled-log inspection. In a clean release with TeX logs intentionally absent, both log checks return `SKIP` while the required compiled PDFs and all scientific consistency checks still pass. The builder runs this release-context audit before writing the root manifest, making a checker rerun manifest-idempotent.

### AI provenance and author-led verification

- `AI_ASSISTANCE_PROVENANCE_RECORD_20260722.md` records the actual Codex-assisted tasks and states that no experimental data were generated.
- `AUTHOR_LED_FINAL_TEXT_VERIFICATION_20260722.md` separates completed machine checks from declarations that only the human authors can make.
- The cover letter contains a factual disclosure draft, explicitly marked for all-author verification.
- Human author signoff remains unchecked; Codex has not claimed to certify author consent, authorship, funding, conflicts or final approval.

## Strongly recommended work completed

### Single-edit, silent-misdecode and CRC composition audit

- Sixty-four emitted codewords were audited under every single substitution, insertion and deletion: 55,552 edit operations in total.
- Fixed length rejected every insertion and deletion.
- The base decoder silently accepted 12,093 edits as a wrong payload; CRC-16/CCITT-FALSE reduced this count to zero in this finite deterministic sample.
- The optional canonical verifier alone did not eliminate silent wrong-payload decoding.
- The manuscript labels this as a bounded error-detection composition audit, not an error-correcting code, stochastic channel model or recovery experiment.
- Evidence: `channel_error_boundary/`, `analysis_tools/audit_paper2_channel_error_boundary.py` and `tests/test_paper2_channel_error_boundary.py`.

### Negative-fiber fraction in the abstract

- The 149-word structured abstract reports that 17.2--30.9% of matched-condition fibers had negative gain and defines the primary estimand as measured selected-minus-fiber-mean gain.
- It continues to report endpoint-dependent CNN comparisons, exploratory Q5 reversal and the absence of synthesis or assay for generated candidates.

### Finite CRC and absolute-effect boundaries

- The main text reports that absolute FullContext gains were modest (`0.0021--0.0041` on the published relative-efficiency scale) and directs practical interpretation to standardized effects, event risk and prospective validation.
- The CRC statement is explicitly finite-sample: zero silent misdecodes in the enumerated audit does not guarantee detection under arbitrary corruptions.

## Remaining author-controlled blockers

### Immutable archive DOI

Zenodo redirected the deposit workflow to sign-in. A DOI cannot be minted without an authenticated deposit and a truthful public creator list. `ZENODO_DEPOSIT_METADATA_20260722.md` contains all technical metadata, but exact creator names, order, affiliations and ORCIDs must be supplied and approved by the authors. No GitHub username or local identity was used as a substitute.

### Submission identity and declarations

The anonymous manuscript intentionally retains an empty author block and withheld Contact. Before journal upload, the authors must provide and approve author order, affiliations, corresponding email, CRediT, funding, conflicts, acknowledgements and the actual AI-use disclosure. These are not scientific-analysis blockers, but they are formal submission blockers.
