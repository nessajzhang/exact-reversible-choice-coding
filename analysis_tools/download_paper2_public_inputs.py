#!/usr/bin/env python3
"""Acquire and verify the public inputs used by Paper 2.

The repository does not redistribute third-party source files.  This helper
checks out the two public upstream repositories at fixed commits, downloads
the publication supplementary files from the publisher, derives the frozen
DT4DDS sequence-level table, and copies the required PCR-bias files into the
flat, hash-locked input directory consumed by the analyses.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PCR_REPO = ROOT / "external_data" / "paper2_gimpel2025_pcr_bias"
DT_REPO = ROOT / "external_data" / "paper2_dt4dds_sequence_recovery"
SUPP_DIR = ROOT / "external_data" / "paper2_gimpel2025_supplement"
FLAT_DIR = ROOT / "external_data" / "paper2_gimpel2025_public_release_af62c57"

PCR_URL = "https://github.com/BorgwardtLab/PCR-bias.git"
PCR_COMMIT = "af62c57f9a90ecdfdd0f1623441e82bdb7e082c1"
DT_URL = "https://github.com/fml-ethz/dt4dds_notebooks.git"
DT_COMMIT = "126e6da5c41f4e5de072b7a1a0934068b743de6c"

SUPPLEMENTS = {
    "41467_2025_64221_MOESM1_ESM.pdf": (
        "https://static-content.springer.com/esm/"
        "art%3A10.1038%2Fs41467-025-64221-4/MediaObjects/"
        "41467_2025_64221_MOESM1_ESM.pdf",
        "a0ebb08bd6c6da6abc32afd91a388b646d579b87409474efbfe647c945c9063c",
    ),
    "41467_2025_64221_MOESM4_ESM.zip": (
        "https://static-content.springer.com/esm/"
        "art%3A10.1038%2Fs41467-025-64221-4/MediaObjects/"
        "41467_2025_64221_MOESM4_ESM.zip",
        "3fdc4f6dbe21f17e4741e13794d91283b7cf4ee2caa4af4d1cbc71adb53dcd68",
    ),
}

FLAT_FILES = {
    "External_GCall2GCfix_2perc_regression_plus_probs.csv": (
        "CNN/results_revision/External_GCall2GCfix_2perc_regression_plus_probs.csv",
        "b694f853322114eb7bf0dbd249c027520e91d5b1d9db93c5b198b011892ac9de",
    ),
    "External_GCfix2GCall_2perc_regression_plus_probs.csv": (
        "CNN/results_revision/External_GCfix2GCall_2perc_regression_plus_probs.csv",
        "7bfc2440f3963505a8885cd4cd18b44fe0c99fffa854956793f083f5f6acbcff",
    ),
    "external-validation-predictions.csv": (
        "analysis/data/machine_learning_results/external-validation-predictions.csv",
        "2195ede4ce1a7172840c3091553e73a394c777a1d61dbad292b0ba556b984fb2",
    ),
    "sequence_data_anonymized_no_duplicates.csv": (
        "analysis/41_external_validation_generation/sequence_data_anonymized_no_duplicates.csv",
        "53e9e4bdd4b740061b1650d3b4c8adcb3012312af20bbfbe5b75748865017864",
    ),
    "GCall_seqprops.csv": (
        "analysis/data/internal_datasets/GCall/seqprops.csv",
        "154ea0a1cf48b368a03037f1681fccf434bb3d3bd66ff4f3197eef55fceff864",
    ),
    "GCfix_seqprops.csv": (
        "analysis/data/internal_datasets/GCfix/seqprops.csv",
        "6463d9d91a96cf33e67362df5e9f07985df61e08f6b49a09f2992da0fa9d48dd",
    ),
    "GCall_params.csv": (
        "analysis/data/internal_datasets/GCall/params.csv",
        "376eb0db5c0c391cc9c5941795516ca26c0dbe6ff289b1e8e4b0981828db7918",
    ),
    "GCfix_params.csv": (
        "analysis/data/internal_datasets/GCfix/params.csv",
        "bd8b345ca06d21af630656b725951c93bd166afb1f23f67b8f04f36334cd3565",
    ),
    "45_external_validation_filtering_run_analysis.ipynb": (
        "analysis/45_external_validation_filtering/run_analysis.ipynb",
        "455fd443ba799224d6cf9185644ce2fbfdc395a4c6cc20aa9c87e296cdeefa16",
    ),
}

DT_DERIVED = {
    ROOT / "data" / "paper2_dt4dds_sequence_detectability.csv":
        "32ed334ae29154c7e515d6d36409d3b7d7399d6f865265460c2e5c0b7f64dba3",
    ROOT / "data" / "paper2_dt4dds_sequence_detectability_summary.csv":
        "78939aa80b9ef33fc72c8c7112ace7eb8dcc4dc4c9f2151c8db2374b83964fb2",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(*args: str) -> str:
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()


def require_hash(path: Path, expected: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    observed = sha256(path)
    if observed != expected:
        raise ValueError(f"SHA-256 mismatch for {path}: {observed} != {expected}")


def checkout(url: str, destination: Path, commit: str) -> None:
    if not (destination / ".git").is_dir():
        if destination.exists() and any(destination.iterdir()):
            raise ValueError(f"refusing to replace non-empty path: {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(
            ["git", "clone", "--filter=blob:none", "--no-checkout", url, str(destination)],
            cwd=ROOT,
        )
    subprocess.check_call(["git", "-C", str(destination), "fetch", "--depth", "1", "origin", commit])
    subprocess.check_call(["git", "-C", str(destination), "checkout", "--detach", commit])
    observed = run("git", "-C", str(destination), "rev-parse", "HEAD")
    if observed != commit:
        raise ValueError(f"commit mismatch for {destination}: {observed} != {commit}")


def download(url: str, destination: Path, expected: str) -> None:
    if destination.is_file() and sha256(destination) == expected:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "Paper2-reproducibility/1.0"})
        with urllib.request.urlopen(request) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output)
        require_hash(temporary, expected)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def prepare_flat_release() -> None:
    FLAT_DIR.mkdir(parents=True, exist_ok=True)
    for local_name, (upstream_path, expected) in FLAT_FILES.items():
        source = PCR_REPO / upstream_path
        destination = FLAT_DIR / local_name
        require_hash(source, expected)
        if not destination.is_file() or sha256(destination) != expected:
            shutil.copy2(source, destination)
        require_hash(destination, expected)


def derive_dt4dds() -> None:
    if all(path.is_file() and sha256(path) == expected for path, expected in DT_DERIVED.items()):
        return
    subprocess.check_call(
        [sys.executable, "analysis_tools/audit_paper2_dt4dds_sequence_detectability.py", "--permutations", "500"],
        cwd=ROOT,
    )
    for path, expected in DT_DERIVED.items():
        require_hash(path, expected)


def verify_all() -> None:
    if run("git", "-C", str(PCR_REPO), "rev-parse", "HEAD") != PCR_COMMIT:
        raise ValueError("PCR-bias checkout is not at the frozen commit")
    if run("git", "-C", str(DT_REPO), "rev-parse", "HEAD") != DT_COMMIT:
        raise ValueError("DT4DDS checkout is not at the frozen commit")
    for name, (_url, expected) in SUPPLEMENTS.items():
        require_hash(SUPP_DIR / name, expected)
    for name, (_upstream, expected) in FLAT_FILES.items():
        require_hash(FLAT_DIR / name, expected)
    for path, expected in DT_DERIVED.items():
        require_hash(path, expected)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if not args.check_only:
        checkout(PCR_URL, PCR_REPO, PCR_COMMIT)
        checkout(DT_URL, DT_REPO, DT_COMMIT)
        for name, (url, expected) in SUPPLEMENTS.items():
            download(url, SUPP_DIR / name, expected)
        prepare_flat_release()
        derive_dt4dds()
    verify_all()
    print("PASS_PAPER2_PUBLIC_INPUT_CONTRACT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
