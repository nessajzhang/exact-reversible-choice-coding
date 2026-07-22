# AI-assistance provenance record

Record date: 2026-07-22

Status: `MACHINE_SIDE_PROVENANCE_COMPLETE; AUTHOR_CONFIRMATION_REQUIRED_BEFORE_SUBMISSION`

## Tool and scope

OpenAI Codex, based on the GPT-5 model family, was used during author-directed manuscript revision. The assistance covered:

- inspection of manuscript, Supplementary Information, code, frozen outputs and rendered PDFs;
- suggestions and direct edits to scientific wording, structure, captions, response material and reproducibility documentation;
- review and modification of analysis code, including the grouped-source bootstrap correction and deterministic single-edit/CRC boundary audit;
- generation of unit tests, fixed numerical test vectors and consistency checks;
- literature discovery and bibliographic field checking against publisher, Crossref and indexed records;
- figure regeneration, release packaging, repository preparation and manifest checking.
- revision of the release consistency checker so omitted TeX logs are reported as `SKIP` rather than causing a missing-file traceback;
- integration of the already-frozen negative Q5 two-stage intervals into the main Figure 2 without changing their values.

## What the tool did not supply

- It did not generate experimental observations, PCR measurements, synthesized oligonucleotides or wet-laboratory results.
- It did not perform prospective validation of codec-emitted sequences.
- It did not determine authorship, author order, affiliations, funding, conflicts, contributions or correspondence details.
- It did not replace author responsibility for mathematical correctness, statistical interpretation, code, references, data licensing or final claims.

## Verification trail

- Numerical claims are tied to versioned scripts and machine-readable frozen outputs rather than accepted from model text.
- The corrected grouped bootstrap was run for 2,000 source-plus-target replicates and compared with the prior frozen baseline.
- Deterministic-selection, grouped-bootstrap and channel-boundary unit tests are executable from the release.
- Added references are recorded in `REFERENCE_VERIFICATION_REPORT_20260722.md` with DOI and field checks.
- Manuscript/Supplement consistency, hashes, compilation and rendered-page checks are separate release gates.
- Twenty-one deterministic, grouped-bootstrap, channel-boundary, release-log, frozen-environment and release-hygiene tests pass in the final technical package.

## Disclosure text for author approval

> During author-directed revision, OpenAI Codex (GPT-5 model family) assisted with code review and modification, test generation, literature-metadata checking, internal-consistency auditing, presentation review and language suggestions. All reported numerical results were produced by the versioned analysis scripts; the tool generated no experimental data. The authors retain responsibility for the mathematical arguments, code, references, analyses and final wording.

The quoted text is a factual draft. All authors must confirm that it accurately describes their complete use of generative-AI tools before it is copied into the cover letter and the journal-designated disclosure location.
