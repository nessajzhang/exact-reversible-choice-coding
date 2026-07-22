# Paper 2: assay-calibrated reversible choice coding

## Current manuscript

Working title: **Exact reversible choice coding with retrospective assay-guided selection in PCR-tested oligonucleotide libraries**

Primary route: a Bioinformatics Original Paper built around an algorithmic interface between exact constrained coding and a source-trained public experimental assay score. The scientific contribution is the reversible choice/fiber construction and its retrospective evaluation in PCR-tested libraries; the work is not positioned as a new software product or as prospective validation of codec-emitted sequences.

Current scientific status: `MINOR_REVISION_TECHNICAL_CLOSEOUT; PYTHON_31210_PORTABILITY_FIXED; PUBLIC_REPOSITORY_AND_LICENSE_READY; SUBMISSION_BLOCKED_BY_AUTHOR_METADATA_AND_PERSISTENT_ARCHIVE_GATE`. The formal construction, public-data analysis, matched baselines, negative transfer result, figures, manuscript and Supplementary Information are claim-bounded. The public repository is `https://github.com/nessajzhang/exact-reversible-choice-coding` and project code is BSD-3-Clause licensed. Before upload, authors must supply and approve the submission metadata and provide a persistent archived snapshot of the submitted software and test data. A DOI is one possible identifier, not the only acceptable form.

## Central result

An exact base codec assigns each payload `2^r` legal sequence representatives. A frozen assay model selects one representative, but decoding uses only exact rerank, inverse permutation and an `r`-bit shift. The score is not needed by the decoder, and the rate cost is exactly `r` bits per oligonucleotide.

The current evidence has two deliberately separate layers:

1. **Measured public layer.** Source-only models trained on one public PCR pool rank a Hamming-separated target pool. Outcome-blind 2,048-sequence codebooks show positive measured selection utility under matched KAPA conditions. At `r=2`, cross-pool standardized gains are 0.250 and 0.274; external-KAPA gains span 0.229--0.390. Comparisons with released 1D-CNN outputs are endpoint-dependent, 17.2--30.9% of matched-condition fibers are negative, and a secondary exploratory Q5 analysis reverses direction.
2. **Generated exact-code layer.** The 108-nt GC 49--59 / homopolymer <=3 language contains exactly `21520867790325216400381593480072139809549357542204871770858223072` words. A fixed 213-bit domain carries 1.972222 bits/nt; `r=2` and `r=4` retain 211 and 209 bits. All 65,024 generated candidate checks round-trip exactly.

The public assayed sequences were not prospectively emitted by the generated codec, and generated sequences have model scores but no measured PCR outcomes.

## Primary files

- `main.tex` and `build/main.pdf`: current main manuscript.
- `supplementary_codec_evidence.tex` and `build/supplementary_codec_evidence.pdf`: current Supplementary Information.
- `figures/assay_calibrated_selection.*`: public calibration, fixed-length Hamming-separation and negative-boundary figure.
- `figures/reversible_choice_codec.*`: choice method, measured effects, sensitivity, rate--score and threshold-baseline figure.
- `bioinformatics_reframe/reversible_choice_codec/`: exact language, finite public codebooks, generated candidates, baselines, seeds, provenance and SHA-256 manifest.
- `bioinformatics_reframe/sequence_independence_audit/`: exhaustive 143,904,012-pair audit and manifest.
- `bioinformatics_reframe/assay_calibrated_selection/`: frozen model, selection and coefficient outputs.
- `bioinformatics_reframe/runtime_benchmark/`: frozen manuscript timing snapshot, environment and manifest; one-command machine-local reruns go to `runtime_benchmark/latest_local/` and do not overwrite the snapshot.
- `bioinformatics_reframe/00_source_inventory.md` through `08_priority_workplan.md`: literature, story, claim and submission audit trail.
- `bioinformatics_reframe/BIOINFORMATICS_LIVE_SUBMISSION_AUDIT_20260717.md`: dated journal-scope, JCR/CAS and LLM-policy audit.
- `bioinformatics_reframe/REFERENCE_VERIFICATION_REPORT_20260722.md`: cited-reference field and status verification.
- `bioinformatics_reframe/MANUSCRIPT_CONSISTENCY_AUDIT.md`: current claim/number/table/figure/manifest boundary PASS.
- `bioinformatics_reframe/PAPER2_SUBMISSION_PREFLIGHT_RESPONSE_20260722.md`: final reviewer-response closeout, including resolved technical items and author-controlled gates.
- `bioinformatics_reframe/REPRODUCIBLE_EXPORT_AUDIT_20260718.md`: byte-identical double-export evidence for figures and compiled PDFs.

## Reproduction entry points

From the repository root:

```bash
./papers/paper2_thermodynamic_risk_coding/reproduce_bioinformatics_choice_study.sh
```

The compatibility command delegates to `run_full_audit.sh`. The workflow is also split into independently callable stages:

```bash
./papers/paper2_thermodynamic_risk_coding/verify_integrity.sh --files-only
./papers/paper2_thermodynamic_risk_coding/verify_integrity.sh --contracts-only
./papers/paper2_thermodynamic_risk_coding/reproduce_analysis.sh
./papers/paper2_thermodynamic_risk_coding/reproduce_figures.sh
./papers/paper2_thermodynamic_risk_coding/compile_manuscript.sh
./papers/paper2_thermodynamic_risk_coding/run_full_audit.sh --skip-tex
```

`--files-only` verifies the 10 nested manifests and, in a built release, the root SHA-256 manifest and transient-artifact policy using only a basic Python interpreter and standard shell hash tools. It neither installs nor requires the canonical scientific environment. `--contracts-only` runs the 21 non-mutating contract tests with an already compatible scientific Python or an isolated canonical fallback. Only analysis and figure reproduction require the exact frozen environment.

`--skip-tex` completes scientific analysis, figure and integrity checks without treating absent LaTeX/BibTeX as an analysis failure. `--analysis-only`, `--figures-only` and `--from-frozen` provide narrower runs. The full public-input run validates source hashes, rebuilds the analyses, rechecks fixed-length Hamming separation, regenerates figures, records a machine-local timing rerun, compiles when both `latexmk` and `bibtex` are available, and then checks manuscript consistency. A fixed `SOURCE_DATE_EPOCH` removes creation-time-only PDF/figure timestamp drift; it does not freeze machine timing.

The manuscript consistency checker reports compiled-log checks as `PASS`, `SKIP` or `FAIL`. Release archives intentionally omit TeX logs, so an absent `build/main.log` or Supplement log is a documented `SKIP`; a present clean log is `PASS`, and a present log containing a declared LaTeX error is `FAIL`. Compiled PDFs remain required release artifacts.

Useful environment variables:

```bash
PYTHON=/path/to/python JOBS=12 SEQ_JOBS=10 \
  ./papers/paper2_thermodynamic_risk_coding/reproduce_bioinformatics_choice_study.sh
```

Set `PAPER2_SKIP_UPSTREAM=1` only for a fast rebuild from the frozen upstream analysis tables already present in the repository. The primary choice script still revalidates source hashes, refits source-only selectors and checks the frozen prediction ledger. `PAPER2_LOCAL_RUNTIME_DIR` may redirect the non-frozen local timing output.

The portable canonical run uses Python 3.12.10 with NumPy 2.5.0, pandas 3.0.3, SciPy 1.18.0, scikit-learn 1.9.0, joblib 1.5.3 and Matplotlib 3.11.0. Exact versions, commands and seeds are stored in the output `environment_and_seeds.json` files. The runner accepts `PYTHON` for canonical output only when Python and all six distribution versions match this frozen contract exactly; an importable but version-mismatched environment is rejected. Otherwise it explicitly asks `uv` to install Python 3.12.10 and creates or reuses the isolated pinned environment at `${PAPER2_UV_ENV:-$HOME/.cache/paper2-thermodynamic-risk-coding/py31210}`. This path was tested with `uv 0.10.0`; it does not modify a system or Anaconda base environment.

The repository workflow `.github/workflows/paper2-linux-portability.yml` runs the root/nested integrity checks, all 21 contract tests, frozen-output reproduction, LaTeX/BibTeX compilation, consistency checking and PDF inspection on Ubuntu 24.04. A manually dispatched second job runs the complete public-input analysis. CI logs and rendered PDF-QC artifacts are retained by GitHub Actions.

## Claim boundaries

Do not convert any of the following into positive evidence:

- no NUPACK result;
- no synthesis, material or emitted-codeword wet-lab result;
- no empirical sequence or file recovery;
- no storage-density or end-to-end system superiority;
- no practical full 110-nt weighted-language capacity;
- the 174-bit construction is only a strict zero-score 110-nt sublanguage;
- the 4,096-byte test is noiseless deterministic serialization only;
- DT4DDS is observational and its controlling weighted-feature increment is null;
- ViennaRNA, Primer3 and seqfold are computational screens only;
- no unpublished sister-paper conclusion is evidence for this manuscript.

## Remaining author-controlled gates

Before upload, authors must supply and verify names, order, affiliations, corresponding email, ORCID identifiers, contributions, funding, conflicts and acknowledgements. Any journal-required declarations remain subject to author verification. The public repository and BSD-3-Clause licence are in place; the author has deferred creation of a persistent archived snapshot in Zenodo, Figshare, Software Heritage or an equivalent service. A Zenodo DOI would require authenticated deposition and truthful creator metadata, whereas another approved archive may issue a different persistent identifier. None of these author-controlled fields is inferred in the manuscript.
