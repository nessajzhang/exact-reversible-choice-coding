# Paper 2 Python 3.12.10 portability audit

Audit date: 2026-07-23

Status: `LOCAL_CANONICAL_MIGRATION_PASS; CLEAN_LINUX_WORKFLOW_CONFIGURED; AUTHOR_CONTROLLED_SUBMISSION_FIELDS_OPEN`

## Scope

This audit closes the release-bootstrap failure caused by requiring Python 3.12.13 through a binary-only `uv` path. It changes the canonical interpreter and release entry points, not the frozen scientific claims. Q5 reversal, negative-fiber fractions, endpoint-dependent CNN comparisons and the absence of prospective assays for generated candidates remain unchanged.

## Environment repair

- Canonical Python: 3.12.10.
- Pinned distributions: joblib 1.5.3, Matplotlib 3.11.0, NumPy 2.5.0, pandas 3.0.3, scikit-learn 1.9.0 and SciPy 1.18.0.
- Bootstrap path: `uv python install 3.12.10`, followed by an isolated `uv venv` and the exact `requirements.txt`.
- Tested bootstrap: `uv 0.10.0` in a new temporary environment; the interpreter and all six distributions satisfied the exact contract.
- Default environment location: `$HOME/.cache/paper2-thermodynamic-risk-coding/py31210`, overridable with `PAPER2_UV_ENV`.

## Entry-point separation

- `verify_integrity.sh --files-only` uses a basic Python interpreter plus `sha256sum` or `shasum`; it does not create or require the canonical scientific environment.
- `verify_integrity.sh --contracts-only` runs the 21 non-mutating tests with an injected compatible interpreter or an isolated canonical fallback.
- `reproduce_analysis.sh`, `reproduce_figures.sh` and canonical-output modes of `run_full_audit.sh` continue to enforce the exact interpreter and package contract.
- The cache-hygiene integration test injects its running interpreter and tests `PYTHONDONTWRITEBYTECODE` directly; it no longer attempts to download a source-only Python patch.

## Local verification evidence

- File-only verification passed while `PYTHON` and `PAPER2_UV_ENV` pointed to deliberately unavailable locations; no canonical environment was created.
- All 21 deterministic-selection, grouped-bootstrap, channel-boundary, log-contract, environment-contract and release-hygiene tests passed in the isolated Python 3.12.10 environment.
- The full public-input analysis completed with the exact Python 3.12.10 contract.
- Seventy-eight scientific CSV/TSV artifacts from the nine non-runtime analysis directories were compared with the previous public commit after normalizing CRLF/LF line endings: 78 identical, 0 different.
- The three manuscript figures in PDF, PNG and SVG form were byte-identical to the previous public commit: 9 identical, 0 different.
- All 10 finalized scientific-output manifests passed and locked 110 files.
- The machine-specific runtime benchmark was rerun under Python 3.12.10; its timing and memory values were updated in the Supplementary Information and remain explicitly non-biological.
- The manuscript and Supplementary Information compiled to 9 and 17 pages; the OUP two-column preflight remained 7 pages. Both TeX logs passed the declared consistency checks, and rendered pages showed no new clipping or overlap.

Machine-readable summary: `PYTHON_31210_MIGRATION_COMPARISON_20260723.tsv`.

## Clean-Linux verification

`.github/workflows/paper2-linux-portability.yml` defines an Ubuntu 24.04 audit with pinned `uv 0.10.0`. Its default job checks the root and nested manifests, runs all 21 tests, reproduces frozen analyses and figures, compiles with LaTeX/BibTeX, runs the consistency checker, inspects PDF metadata/fonts and renders every PDF page. A manually dispatched second job runs the complete public-input audit. GitHub Actions retains the text logs, compiled PDFs and rendered-page evidence.

The workflow result and immutable run URL must be inserted here only after the public commit has been pushed and the remote run has completed. A configured workflow is not described as a passed Linux run before that event.

## Remaining author-controlled actions

No author name, affiliation, correspondence address, ORCID, CRediT statement, funding statement, conflict declaration, acknowledgement, release tag or persistent archive identifier was inferred. Those actions require approved creator metadata and, for an external deposit, an authenticated publication action. The AI-assistance record and cover-letter wording remain drafts requiring all-author factual verification; they were not inserted into the anonymous manuscript.
