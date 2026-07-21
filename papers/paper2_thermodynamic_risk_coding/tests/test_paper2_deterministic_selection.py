#!/usr/bin/env python3
"""Cross-platform contract tests for Paper 2 deterministic selection."""

from __future__ import annotations

import csv
import sys
import unittest
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "analysis_tools"))

import paper2_deterministic_selection as deterministic  # noqa: E402


class DeterministicSelectionContractTest(unittest.TestCase):
    def test_frozen_round_to_even_vectors(self) -> None:
        vector_path = Path(__file__).with_name("data") / "paper2_score_quantization_vectors.tsv"
        with vector_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        scaled = np.asarray([float(row["scaled_binary64"]) for row in rows], dtype=np.float64)
        expected = np.asarray([int(row["expected_int64"]) for row in rows], dtype=np.int64)
        np.testing.assert_array_equal(deterministic.round_to_even_int64(scaled), expected)

    def test_score_scaling_and_negative_scores(self) -> None:
        scores = np.asarray([-2.0e-12, -1.0e-12, 0.0, 1.0e-12, 2.0e-12])
        np.testing.assert_array_equal(
            deterministic.quantized_score_keys(scores),
            np.asarray([-2, -1, 0, 1, 2], dtype=np.int64),
        )

    def test_smallest_choice_index_breaks_equal_key_ties(self) -> None:
        self.assertEqual(deterministic.argmax_smallest([1.0, 1.0, -1.0]), 0)
        self.assertEqual(deterministic.argmin_smallest([-1.0, -1.0, 1.0]), 0)

    def test_nonfinite_values_are_rejected(self) -> None:
        for value in (np.nan, np.inf, -np.inf):
            with self.subTest(value=value), self.assertRaises(ValueError):
                deterministic.quantized_score_keys([value])

    def test_signed_int64_range_is_checked_after_rounding(self) -> None:
        self.assertEqual(
            int(deterministic.round_to_even_int64([float(-(1 << 63))])[0]),
            -(1 << 63),
        )
        with self.assertRaises(OverflowError):
            deterministic.round_to_even_int64([float(1 << 63)])

    def test_reference_contract_metadata_is_explicit(self) -> None:
        self.assertEqual(deterministic.SCORE_SCALE, 10**12)
        self.assertIn("ties-to-even", deterministic.ROUNDING_CONTRACT)


if __name__ == "__main__":
    unittest.main()
