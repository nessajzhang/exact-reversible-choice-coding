#!/usr/bin/env python3
"""Regenerate portable manifests for the frozen Paper 2 output directories."""

from __future__ import annotations

import argparse
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


def manifestable_files(directory: Path) -> list[Path]:
    return sorted(
        path
        for path in directory.rglob("*")
        if path.is_file()
        and path.name != "sha256_manifest.tsv"
        and "latest_local" not in path.relative_to(directory).parts
    )


def write_manifest(directory: Path) -> tuple[int, str]:
    rows = []
    for path in manifestable_files(directory):
        rows.append(
            {
                "file": path.relative_to(directory).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest = directory / "sha256_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("file", "bytes", "sha256"), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows), sha256(manifest)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=REFRAME,
        help="bioinformatics_reframe directory containing the frozen output folders",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    for name in OUTPUT_DIRS:
        directory = root / name
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        count, digest = write_manifest(directory)
        print(f"{name}: {count} files; manifest {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
