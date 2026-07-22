#!/usr/bin/env python3

from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analysis_tools"))

import verify_paper2_environment as environment  # noqa: E402


class FrozenEnvironmentContractTests(unittest.TestCase):
    def test_exact_environment_passes(self) -> None:
        self.assertEqual(
            environment.compare_environment(
                environment.CANONICAL_PYTHON,
                environment.CANONICAL_DISTRIBUTIONS,
            ),
            [],
        )

    def test_python_patch_version_is_part_of_contract(self) -> None:
        mismatches = environment.compare_environment(
            (3, 12, 13), environment.CANONICAL_DISTRIBUTIONS
        )
        self.assertEqual(
            mismatches,
            ["python: expected 3.12.10, found 3.12.13"],
        )

    def test_missing_and_wrong_distribution_versions_fail(self) -> None:
        versions = dict(environment.CANONICAL_DISTRIBUTIONS)
        versions["numpy"] = "2.4.0"
        versions["scipy"] = None
        mismatches = environment.compare_environment(
            environment.CANONICAL_PYTHON, versions
        )
        self.assertEqual(
            mismatches,
            [
                "numpy: expected 2.5.0, found 2.4.0",
                "scipy: expected 1.18.0, found missing",
            ],
        )


if __name__ == "__main__":
    unittest.main()
