#!/usr/bin/env python3
"""Build a clean, portable Paper 2 public-release repository tree."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
REFRAME = PAPER / "bioinformatics_reframe"
DEFAULT_OUTPUT = PAPER / "public_release" / "exact-reversible-choice-coding"

ANALYSIS_FILES = (
    "analyze_paper2_external_mapping_sensitivity.py",
    "analyze_paper2_major_revision_diagnostics.py",
    "audit_paper2_channel_error_boundary.py",
    "audit_paper2_dt4dds_sequence_detectability.py",
    "audit_paper2_pcr_sequence_independence.py",
    "audit_paper2_source_external_independence.py",
    "benchmark_paper2_reversible_choice_codec.py",
    "build_paper2_public_release.py",
    "download_paper2_public_inputs.py",
    "finalize_paper2_output_manifests.py",
    "hairpin_risk_features.py",
    "make_paper2_assay_calibrated_selection_figure.py",
    "make_submission_figures.py",
    "paper2_deterministic_selection.py",
    "plot_paper2_reversible_choice_codec.py",
    "validate_paper2_assay_calibrated_selection.py",
    "validate_paper2_public_experimental_data.py",
    "validate_paper2_reversible_choice_codec.py",
    "validate_paper2_sota_external_uncertainty.py",
    "verify_paper2_environment.py",
    "verify_paper2_bioinformatics_choice_consistency.py",
    "verify_paper2_output_manifests.py",
)

PUBLIC_INPUT_RUNTIME_FILES = (
    "analysis_tools/download_paper2_public_inputs.py",
    "analysis_tools/audit_paper2_dt4dds_sequence_detectability.py",
    "papers/paper2_thermodynamic_risk_coding/reproduce_from_public_inputs.sh",
)

PAPER_FILES = (
    ".python-version",
    ".gitignore",
    "LICENSE",
    "THIRD_PARTY_NOTICES.md",
    "VERSION",
    "requirements.txt",
    "README.md",
    "main.tex",
    "supplementary_codec_evidence.tex",
    "references.bib",
    "paper2_repro_common.sh",
    "reproduce_analysis.sh",
    "reproduce_bioinformatics_choice_study.sh",
    "reproduce_figures.sh",
    "reproduce_from_frozen_outputs.sh",
    "reproduce_from_public_inputs.sh",
    "compile_manuscript.sh",
    "run_full_audit.sh",
    "verify_contracts.sh",
    "verify_file_integrity.sh",
    "verify_integrity.sh",
)

CI_FILES = (
    "paper2-linux-portability.yml",
)

FIGURE_FILES = (
    "assay_calibrated_selection.pdf",
    "assay_calibrated_selection.png",
    "assay_calibrated_selection.svg",
    "paired_two_stage_sensitivity.pdf",
    "paired_two_stage_sensitivity.png",
    "paired_two_stage_sensitivity.svg",
    "paired_two_stage_sensitivity_source_data.tsv",
    "reversible_choice_codec.pdf",
    "reversible_choice_codec.png",
    "reversible_choice_codec.svg",
    "reversible_choice_codec_source_data.tsv",
)

DATA_FILES = (
    "paper2_assay_calibrated_selection_figure_source_data.tsv",
)

BUILD_FILES = (
    "main.pdf",
    "supplementary_codec_evidence.pdf",
)

OUP_PREFLIGHT_FILES = (
    "main_oup_preflight.tex",
    "build/main_oup_preflight.pdf",
)

OUTPUT_DIRS = (
    "public_experimental_validation",
    "assay_calibrated_selection",
    "sequence_independence_audit",
    "source_external_independence_audit",
    "reversible_choice_codec",
    "sota_and_external_validation",
    "external_mapping_sensitivity",
    "major_revision_diagnostics",
    "channel_error_boundary",
    "runtime_benchmark",
)

RELEASE_DOCUMENTS = (
    "GROUPED_BOOTSTRAP_CORRECTION_AUDIT_20260722.md",
    "REFERENCE_VERIFICATION_REPORT_20260722.md",
    "AI_ASSISTANCE_PROVENANCE_RECORD_20260722.md",
    "AUTHOR_LED_FINAL_TEXT_VERIFICATION_20260722.md",
    "BIOINFORMATICS_AUTHOR_CONTROLLED_FIELDS_20260722.md",
    "BIOINFORMATICS_COVER_LETTER_20260722.md",
    "LOCKED_EXTERNAL_VALIDATION_PROTOCOL_20260718.md",
    "MANUSCRIPT_CONSISTENCY_AUDIT.md",
    "OUP_PAGE_PREFLIGHT_20260722.md",
    "PAPER2_SUBMISSION_PREFLIGHT_CHANGE_INDEX_20260722.tsv",
    "PAPER2_SUBMISSION_PREFLIGHT_RESPONSE_20260722.md",
    "PDF_VISUAL_QC_20260722.md",
    "PRE_SUBMISSION_COMPLETION_REPORT_20260722.md",
    "PYTHON_31210_MIGRATION_COMPARISON_20260723.tsv",
    "PYTHON_31210_PORTABILITY_AUDIT_20260723.md",
    "ZENODO_DEPOSIT_METADATA_20260722.md",
)

RELEASE_EVIDENCE_DIRS = (
    "linux_portability_evidence_20260723",
)

TEXT_SUFFIXES = {
    ".bib",
    ".cfg",
    ".csv",
    ".ini",
    ".json",
    ".log",
    ".md",
    ".py",
    ".sh",
    ".tex",
    ".toml",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}
TEXT_NAMES = {".gitignore", "LICENSE", "VERSION", "requirements.txt"}
POSIX_USER_PREFIXES = ("/" + "Users" + "/", "/" + "home" + "/")


def copy_file(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def copy_tree(source: Path, destination: Path) -> None:
    if not source.is_dir():
        raise FileNotFoundError(source)
    shutil.copytree(
        source,
        destination,
        ignore=shutil.ignore_patterns("latest_local", "__pycache__", ".DS_Store"),
    )


def is_text(path: Path) -> bool:
    return path.name in TEXT_NAMES or path.suffix.lower() in TEXT_SUFFIXES


def sanitize_text(path: Path) -> None:
    if not is_text(path):
        return
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return
    text = text.replace(str(ROOT), "$PROJECT_ROOT")
    for prefix in POSIX_USER_PREFIXES:
        text = re.sub(re.escape(prefix) + r"[^/\s\"']+", "$HOME", text)
    text = re.sub(r"[A-Za-z]:\\\\Users\\\\[^\\\s\"']+", "$HOME", text)
    path.write_text(text, encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_root_manifest(output: Path) -> None:
    rows = []
    for path in sorted(output.rglob("*")):
        if not path.is_file() or path.name == "SHA256SUMS.txt" or ".git" in path.parts:
            continue
        rows.append(f"{sha256(path)}  {path.relative_to(output).as_posix()}")
    (output / "SHA256SUMS.txt").write_text("\n".join(rows) + "\n", encoding="utf-8")


def check_no_local_identity(output: Path) -> None:
    forbidden = (
        *(re.compile(re.escape(prefix)) for prefix in POSIX_USER_PREFIXES),
        re.compile(r"[A-Za-z]:\\\\Users\\\\"),
        re.compile(r"\b" + "shi" + "ki" + r"\b", re.I),
        re.compile("Live" + "Sync"),
        re.compile("DNA" + "存储编码顶刊论文项目"),
    )
    leaks: list[str] = []
    for path in sorted(output.rglob("*")):
        if not path.is_file() or not is_text(path):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if any(pattern.search(text) for pattern in forbidden):
            leaks.append(path.relative_to(output).as_posix())
    if leaks:
        raise RuntimeError("local identity/path leakage remains: " + ", ".join(leaks))


def check_no_transient_artifacts(output: Path) -> None:
    """Reject interpreter, editor and operating-system cache artifacts."""

    transient_names = {".DS_Store", "Thumbs.db"}
    transient_suffixes = {".pyc", ".pyo"}
    offenders = [
        path.relative_to(output).as_posix()
        for path in sorted(output.rglob("*"))
        if (
            "__pycache__" in path.relative_to(output).parts
            or path.name in transient_names
            or path.suffix.lower() in transient_suffixes
        )
    ]
    if offenders:
        raise RuntimeError(
            "transient build artifacts entered release: " + ", ".join(offenders)
        )


def check_public_input_runtime_complete(output: Path) -> None:
    """Fail the release build if a public-input entrypoint dependency is absent."""

    missing = [name for name in PUBLIC_INPUT_RUNTIME_FILES if not (output / name).is_file()]
    if missing:
        raise RuntimeError(
            "public-input runtime dependency is missing: " + ", ".join(missing)
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    allowed_parent = (PAPER / "public_release").resolve()
    if allowed_parent not in output.parents:
        raise ValueError(f"output must be below {allowed_parent}")
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)

    for name in ANALYSIS_FILES:
        copy_file(ROOT / "analysis_tools" / name, output / "analysis_tools" / name)

    release_paper = output / "papers" / "paper2_thermodynamic_risk_coding"
    for name in PAPER_FILES:
        copy_file(PAPER / name, release_paper / name)
    for name in CI_FILES:
        copy_file(PAPER / "ci" / name, output / ".github" / "workflows" / name)
    copy_tree(PAPER / "tests", release_paper / "tests")
    copy_tree(PAPER / "examples", release_paper / "examples")
    for name in FIGURE_FILES:
        copy_file(PAPER / "figures" / name, release_paper / "figures" / name)
    for name in DATA_FILES:
        copy_file(ROOT / "data" / name, output / "data" / name)
    for name in BUILD_FILES:
        copy_file(PAPER / "build" / name, release_paper / "build" / name)
    for name in OUP_PREFLIGHT_FILES:
        copy_file(PAPER / "oup_preflight" / name, release_paper / "oup_preflight" / name)
    for name in OUTPUT_DIRS:
        copy_tree(REFRAME / name, release_paper / "bioinformatics_reframe" / name)
    for name in RELEASE_DOCUMENTS:
        copy_file(REFRAME / name, release_paper / "bioinformatics_reframe" / name)
    for name in RELEASE_EVIDENCE_DIRS:
        copy_tree(REFRAME / name, release_paper / "bioinformatics_reframe" / name)

    for name in (
        ".python-version",
        "LICENSE",
        "THIRD_PARTY_NOTICES.md",
        "VERSION",
        "requirements.txt",
    ):
        copy_file(PAPER / name, output / name)
    copy_file(PAPER / "README.md", output / "README.md")
    copy_file(PAPER / ".gitignore", output / ".gitignore")

    for path in sorted(output.rglob("*")):
        if path.is_file():
            sanitize_text(path)

    # Sanitization changes portable provenance files, so rebuild nested locks.
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(output / "analysis_tools"))
    import finalize_paper2_output_manifests as finalizer

    release_reframe = release_paper / "bioinformatics_reframe"
    for name in finalizer.OUTPUT_DIRS:
        finalizer.write_manifest(release_reframe / name)

    # Generate the release-context audit before the root manifest.  Review
    # archives intentionally omit TeX logs, so the two log checks are SKIP;
    # running the checker again must not invalidate SHA256SUMS.txt.
    subprocess.run(
        [
            sys.executable,
            str(output / "analysis_tools" / "verify_paper2_bioinformatics_choice_consistency.py"),
        ],
        cwd=output,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        check=True,
    )
    check_public_input_runtime_complete(output)
    check_no_local_identity(output)
    check_no_transient_artifacts(output)
    write_root_manifest(output)

    print(f"Public release tree: {output}")
    print(f"Files: {sum(1 for path in output.rglob('*') if path.is_file())}")
    print(f"Bytes: {sum(path.stat().st_size for path in output.rglob('*') if path.is_file())}")
    print(f"Root manifest: {sha256(output / 'SHA256SUMS.txt')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
