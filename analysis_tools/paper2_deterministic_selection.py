#!/usr/bin/env python3
"""Deterministic score-key and tie-breaking helpers for Paper 2.

The payload codec never consumes these scores.  They are used only by the
encoder and by the optional canonical-output verifier.  Floating predictions
are converted to signed integer keys at a declared precision before ranking;
the first extremum along each fiber is therefore the smallest choice index.
"""

from __future__ import annotations

from typing import Any

import numpy as np


SCORE_DECIMAL_PLACES = 12
SCORE_SCALE = 10**SCORE_DECIMAL_PLACES
INT64_MIN = -(1 << 63)
INT64_MAX = (1 << 63) - 1
ROUNDING_CONTRACT = "IEEE-754 binary64 round-to-nearest, ties-to-even via numpy.rint"


def round_to_even_int64(scaled_values: Any) -> np.ndarray:
    """Apply the canonical RTE rule and return range-checked int64 keys.

    Inputs are interpreted as IEEE-754 binary64 values that have already been
    scaled.  ``numpy.rint`` supplies round-to-nearest, ties-to-even.  Conversion
    through Python integers makes the signed-int64 range check exact, including
    the asymmetric endpoints, before constructing the NumPy array.
    """

    scaled = np.asarray(scaled_values, dtype=np.float64)
    if not np.isfinite(scaled).all():
        raise ValueError("scaled selector scores must be finite float64 values")
    rounded = np.rint(scaled)
    python_keys = [int(value) for value in rounded.ravel()]
    if any(value < INT64_MIN or value > INT64_MAX for value in python_keys):
        raise OverflowError("rounded selector key exceeds the signed int64 range")
    return np.asarray(python_keys, dtype=np.int64).reshape(rounded.shape)


def quantized_score_keys(values: Any) -> np.ndarray:
    """Return scores under the canonical binary64/RTE/int64 contract."""

    scores = np.asarray(values, dtype=np.float64)
    if not np.isfinite(scores).all():
        raise ValueError("selector scores must be finite float64 values")
    with np.errstate(over="ignore", invalid="ignore"):
        scaled = scores * np.float64(SCORE_SCALE)
    if not np.isfinite(scaled).all():
        raise OverflowError("scaled selector score is outside the finite binary64 range")
    return round_to_even_int64(scaled)


def argmax_smallest(values: Any, axis: int = -1) -> np.ndarray | int:
    """Select the largest quantized key, breaking ties by smallest index."""

    keys = quantized_score_keys(values)
    result = np.argmax(keys, axis=axis)
    return int(result) if np.ndim(result) == 0 else result


def argmin_smallest(values: Any, axis: int = -1) -> np.ndarray | int:
    """Select the smallest quantized key, breaking ties by smallest index."""

    keys = quantized_score_keys(values)
    result = np.argmin(keys, axis=axis)
    return int(result) if np.ndim(result) == 0 else result


def choice_stability_audit(values: Any, mode: str) -> dict[str, float | int]:
    """Compare raw-float and declared-key choices for a two-dimensional fiber array."""

    scores = np.asarray(values, dtype=np.float64)
    if scores.ndim != 2:
        raise ValueError("choice stability audit requires a two-dimensional fiber array")
    keys = quantized_score_keys(scores)
    if mode == "max":
        raw_choice = np.argmax(scores, axis=1)
        key_choice = np.argmax(keys, axis=1)
        ordered = np.sort(scores, axis=1)
        raw_margin = ordered[:, -1] - ordered[:, -2] if scores.shape[1] > 1 else np.full(len(scores), np.inf)
        extreme = keys.max(axis=1, keepdims=True)
    elif mode == "min":
        raw_choice = np.argmin(scores, axis=1)
        key_choice = np.argmin(keys, axis=1)
        ordered = np.sort(scores, axis=1)
        raw_margin = ordered[:, 1] - ordered[:, 0] if scores.shape[1] > 1 else np.full(len(scores), np.inf)
        extreme = keys.min(axis=1, keepdims=True)
    else:
        raise ValueError(f"unsupported choice mode: {mode}")
    tie_count = (keys == extreme).sum(axis=1)
    finite_margin = raw_margin[np.isfinite(raw_margin)]
    return {
        "fibers": int(len(scores)),
        "choices_changed_by_quantization": int((raw_choice != key_choice).sum()),
        "quantized_extremum_tie_fibers": int((tie_count > 1).sum()),
        "minimum_raw_top_gap": float(finite_margin.min()) if len(finite_margin) else float("inf"),
        "median_raw_top_gap": float(np.median(finite_margin)) if len(finite_margin) else float("inf"),
    }
