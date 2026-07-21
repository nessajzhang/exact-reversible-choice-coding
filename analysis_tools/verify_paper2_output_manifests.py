#!/usr/bin/env python3
"""Verify the frozen Paper 2 scientific-output manifests."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFRAME = ROOT / "papers" / "paper2_thermodynamic_risk_coding" / "bioinformatics_reframe"
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_manifested_path(manifest: Path, row: dict[str, str]) -> Path:
    if "file" in row:
        return manifest.parent / row["file"]
    if "relative_path" in row:
        return manifest.parent / row["relative_path"]
    raw = Path(row["path"])
    if raw.is_absolute():
        return raw
    project_relative = ROOT / raw
    return project_relative if project_relative.exists() else manifest.parent / raw


def main() -> int:
    checked_files = 0
    for directory_name in OUTPUT_DIRS:
        manifest = REFRAME / directory_name / "sha256_manifest.tsv"
        if not manifest.is_file():
            raise FileNotFoundError(f"missing manifest: {manifest.relative_to(ROOT)}")
        with manifest.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            path = resolve_manifested_path(manifest, row)
            if not path.is_file():
                raise FileNotFoundError(f"missing manifested file: {path}")
            observed_bytes = path.stat().st_size
            if observed_bytes != int(row["bytes"]):
                raise RuntimeError(f"size mismatch: {path}")
            if sha256(path) != row["sha256"]:
                raise RuntimeError(f"SHA-256 mismatch: {path}")
            checked_files += 1
    print(f"Paper 2 scientific-output integrity passed: {len(OUTPUT_DIRS)} manifests, {checked_files} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
