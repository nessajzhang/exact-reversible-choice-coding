#!/usr/bin/env python3
"""Verify the exact Python and package versions used for canonical Paper 2 outputs."""

from __future__ import annotations

import argparse
import json
import sys
from importlib import metadata
from typing import Mapping, Sequence


CANONICAL_PYTHON = (3, 12, 10)
CANONICAL_DISTRIBUTIONS = {
    "joblib": "1.5.3",
    "matplotlib": "3.11.0",
    "numpy": "2.5.0",
    "pandas": "3.0.3",
    "scikit-learn": "1.9.0",
    "scipy": "1.18.0",
}


def compare_environment(
    python_version: Sequence[int], distributions: Mapping[str, str | None]
) -> list[str]:
    """Return deterministic mismatch messages for a proposed environment."""

    mismatches: list[str] = []
    actual_python = tuple(int(value) for value in python_version[:3])
    if actual_python != CANONICAL_PYTHON:
        mismatches.append(
            "python: expected "
            + ".".join(map(str, CANONICAL_PYTHON))
            + ", found "
            + ".".join(map(str, actual_python))
        )
    for name, expected in CANONICAL_DISTRIBUTIONS.items():
        actual = distributions.get(name)
        if actual != expected:
            mismatches.append(f"{name}: expected {expected}, found {actual or 'missing'}")
    return mismatches


def current_environment() -> tuple[tuple[int, int, int], dict[str, str | None]]:
    versions: dict[str, str | None] = {}
    for name in CANONICAL_DISTRIBUTIONS:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return tuple(sys.version_info[:3]), versions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    python_version, distributions = current_environment()
    mismatches = compare_environment(python_version, distributions)
    report = {
        "status": "PASS" if not mismatches else "FAIL",
        "canonical_python": ".".join(map(str, CANONICAL_PYTHON)),
        "actual_python": ".".join(map(str, python_version)),
        "canonical_distributions": CANONICAL_DISTRIBUTIONS,
        "actual_distributions": distributions,
        "mismatches": mismatches,
    }
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    elif not args.quiet or mismatches:
        if mismatches:
            print("Paper 2 frozen-environment check failed:", file=sys.stderr)
            for mismatch in mismatches:
                print(f"- {mismatch}", file=sys.stderr)
        else:
            print("PASS_PAPER2_FROZEN_ENVIRONMENT")
    return 0 if not mismatches else 2


if __name__ == "__main__":
    raise SystemExit(main())
