#!/usr/bin/env python3
"""Deterministic tests for the original-sequence grouped bootstrap CV contract."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analysis_tools"))

import validate_paper2_sota_external_uncertainty as sota  # noqa: E402


class GroupedBootstrapContractTests(unittest.TestCase):
    def test_duplicate_source_ids_never_cross_folds(self) -> None:
        groups = np.asarray([0, 0, 1, 2, 2, 2, 3, 4, 5, 5, 6, 7, 8, 9])
        test_fold = sota.fixed_group_test_fold(
            groups,
            original_group_count=10,
            alpha_seed=20260718,
        )
        self.assertEqual(set(test_fold.tolist()), set(range(5)))
        for group in np.unique(groups):
            self.assertEqual(len(np.unique(test_fold[groups == group])), 1)

    def test_grouped_fast_and_sklearn_reference_match(self) -> None:
        rng = np.random.default_rng(20260722)
        original_n = 80
        features = rng.normal(size=(original_n, 4))
        outcome = features @ np.asarray([0.4, -0.2, 0.0, 0.1]) + rng.normal(
            scale=0.2, size=original_n
        )
        sampled = rng.integers(0, original_n, size=original_n)
        bootstrap_x = features[sampled]
        bootstrap_y = outcome[sampled]
        fast = sota.fit_fast_ridge_cv(
            bootstrap_x,
            bootstrap_y,
            alpha_seed=12345,
            source_group_ids=sampled,
            original_group_count=original_n,
        )
        reference = sota.fit_bootstrap_model(
            bootstrap_x,
            bootstrap_y,
            alpha_seed=12345,
            source_group_ids=sampled,
            original_group_count=original_n,
        )
        probe = rng.normal(size=(17, 4))
        self.assertEqual(fast.alpha, reference.best_params_["model__alpha"])
        self.assertLess(
            float(np.max(np.abs(fast.predict(probe) - reference.predict(probe)))),
            2e-10,
        )

    def test_out_of_domain_group_id_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            sota.fixed_group_test_fold(
                np.asarray([0, 1, 5]),
                original_group_count=5,
                alpha_seed=1,
            )


if __name__ == "__main__":
    unittest.main()
