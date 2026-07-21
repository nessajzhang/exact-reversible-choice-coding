#!/usr/bin/env python3
"""Validate an assay-calibrated reversible choice codec for Paper 2.

The analysis joins the empirical and formal layers without pretending that a
continuous PCR predictor is a finite-state hard constraint.  A fixed exact
base codec assigns ``2**r`` candidate sequences to each payload.  A predictor
fitted on one public PCR pool chooses one candidate, and the decoder recovers
the payload from the exact base rank after discarding the ``r`` selector bits.

Two evidence layers are kept distinct:

1. held-out finite libraries contain public sequences with measured PCR
   efficiencies and therefore test paired experimental selection utility; and
2. an exact 108-nt GC/homopolymer language demonstrates a practical total
   payload map and the exact selector-bit rate cost on generated sequences.

The public data remain observational with respect to this retrospective
codebook construction.  Generated sequences have model scores but no measured
PCR outcomes.  This program does not provide material, NUPACK, sequencing,
recovery, storage-density-superiority, or causal evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import sys
from collections import OrderedDict
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import ks_2samp, wasserstein_distance
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import hairpin_risk_features as hrf  # noqa: E402
import paper2_deterministic_selection as deterministic  # noqa: E402
import validate_paper2_assay_calibrated_selection as selection  # noqa: E402
import validate_paper2_public_experimental_data as base  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
OUT_DIR = PAPER_DIR / "bioinformatics_reframe" / "reversible_choice_codec"
FROZEN_PREDICTIONS = (
    PAPER_DIR
    / "bioinformatics_reframe"
    / "assay_calibrated_selection"
    / "sequence_level_predictions.tsv"
)

BASE_SEED = 20260721
BOOTSTRAP_REPLICATES = 2000
PUBLIC_LIBRARY_SIZE = 2048
PUBLIC_BASE_BITS = 11
PUBLIC_SELECTOR_BITS = (0, 1, 2, 3, 4)
GENERATED_SELECTOR_BITS = (0, 1, 2, 3, 4, 5, 6)
GENERATED_LENGTH = 108
GENERATED_GC_MIN = 49
GENERATED_GC_MAX = 59
GENERATED_MAX_HOMOPOLYMER = 3
PRIMARY_SELECTOR_BITS = (2, 4)
GENERATED_NAMESPACE_INDICES = tuple(range(8))
GENERATED_NAMESPACE_PAYLOADS = 128

MODEL_NAMES = (
    "P0_composition",
    "P2_assay_context",
    "P5_combined_context",
)
MODEL_FEATURES = OrderedDict((name, selection.PCR_MODELS[name]) for name in MODEL_NAMES)
MODEL_ORIGINAL_INDEX = {
    name: list(selection.PCR_MODELS).index(name) for name in MODEL_NAMES
}

CLAIM_BOUNDARY = (
    "retrospective source-only assay calibration and exact finite-domain choice coding; "
    "public sequences were not prospectively emitted; generated sequences have no measured "
    "PCR outcome; not mechanism, material, NUPACK, sequencing, recovery, full weighted-"
    "language capacity, achieved storage density, or superiority evidence"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(12, (os.cpu_count() or 2) - 2)),
        help="Parallel source-model and generated-feature jobs.",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=BOOTSTRAP_REPLICATES,
        help=f"Paired fiber bootstrap replicates (frozen: {BOOTSTRAP_REPLICATES}).",
    )
    parser.add_argument(
        "--generated-payloads",
        type=int,
        default=512,
        help="Deterministic payload samples per selector-bit setting.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


def write_table(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def max_homopolymer(sequence: str) -> int:
    return int(hrf.max_homopolymer(sequence))


def exact_gc_hp_eligible(sequence: str) -> bool:
    gc_count = sequence.count("G") + sequence.count("C")
    return (
        len(sequence) == GENERATED_LENGTH
        and GENERATED_GC_MIN <= gc_count <= GENERATED_GC_MAX
        and max_homopolymer(sequence) <= GENERATED_MAX_HOMOPOLYMER
    )


@dataclass(frozen=True)
class FrozenRidge:
    source_pool: str
    target_pool: str
    model_name: str
    feature_names: tuple[str, ...]
    alpha: float
    pipeline: Pipeline

    def predict(self, frame: pd.DataFrame) -> np.ndarray:
        return self.pipeline.predict(frame[list(self.feature_names)].to_numpy(float))


def fit_one_source_model(
    source: pd.DataFrame,
    target: pd.DataFrame,
    source_pool: str,
    target_pool: str,
    model_name: str,
    seed: int,
) -> tuple[FrozenRidge, np.ndarray, list[dict[str, Any]]]:
    features = MODEL_FEATURES[model_name]
    inner_seed = seed + MODEL_ORIGINAL_INDEX[model_name]
    estimator = Pipeline([("scale", StandardScaler()), ("model", Ridge())])
    search = GridSearchCV(
        estimator,
        {"model__alpha": base.RIDGE_ALPHA_GRID},
        scoring="neg_mean_squared_error",
        cv=KFold(n_splits=5, shuffle=True, random_state=inner_seed),
        n_jobs=1,
        refit=True,
    )
    search.fit(source[features].to_numpy(float), source["eff"].to_numpy(float))
    target_prediction = search.predict(target[features].to_numpy(float))
    fitted = FrozenRidge(
        source_pool=source_pool,
        target_pool=target_pool,
        model_name=model_name,
        feature_names=tuple(features),
        alpha=float(search.best_params_["model__alpha"]),
        pipeline=search.best_estimator_,
    )

    scaler = search.best_estimator_.named_steps["scale"]
    ridge = search.best_estimator_.named_steps["model"]
    coefficients = np.asarray(ridge.coef_).ravel()
    raw_coefficients = coefficients / np.asarray(scaler.scale_)
    raw_intercept = float(
        ridge.intercept_
        - np.sum(coefficients * np.asarray(scaler.mean_) / np.asarray(scaler.scale_))
    )
    rows: list[dict[str, Any]] = []
    for index, feature in enumerate(features):
        rows.append(
            {
                "source_pool": source_pool,
                "target_pool": target_pool,
                "model": model_name,
                "feature": feature,
                "selected_alpha": fitted.alpha,
                "scaler_mean": float(scaler.mean_[index]),
                "scaler_scale": float(scaler.scale_[index]),
                "standardized_coefficient": float(coefficients[index]),
                "raw_coefficient": float(raw_coefficients[index]),
                "raw_intercept": raw_intercept,
            }
        )
    return fitted, target_prediction, rows


def fit_frozen_models(
    pcr: pd.DataFrame, jobs: int
) -> tuple[dict[tuple[str, str], FrozenRidge], list[dict[str, Any]]]:
    frozen = pd.read_csv(FROZEN_PREDICTIONS, sep="\t", dtype={"sequence_id": str})
    frozen = frozen.loc[
        frozen["dataset"].eq("Gimpel2025_PCR")
        & frozen["evaluation"].eq("source_only_transfer")
    ].copy()

    tasks: list[tuple[pd.DataFrame, pd.DataFrame, str, str, str, int]] = []
    for direction_index, (source_pool, target_pool) in enumerate(
        (("GCall", "GCfix"), ("GCfix", "GCall"))
    ):
        source = pcr.loc[pcr["pool"].eq(source_pool)].reset_index(drop=True)
        target = pcr.loc[pcr["pool"].eq(target_pool)].reset_index(drop=True)
        seed = selection.BASE_SEED + 20_000 + direction_index * 1000
        for model_name in MODEL_NAMES:
            tasks.append((source, target, source_pool, target_pool, model_name, seed))

    fitted_rows = Parallel(n_jobs=min(jobs, len(tasks)), prefer="processes")(
        delayed(fit_one_source_model)(*task) for task in tasks
    )
    models: dict[tuple[str, str], FrozenRidge] = {}
    coefficient_rows: list[dict[str, Any]] = []
    for fitted, target_prediction, rows in fitted_rows:
        direction = f"{fitted.source_pool}_to_{fitted.target_pool}"
        expected = frozen.loc[frozen["pool_or_direction"].eq(direction)].copy()
        expected["sequence_id"] = expected["sequence_id"].str.zfill(6)
        target = pcr.loc[pcr["pool"].eq(fitted.target_pool)].reset_index(drop=True)
        if expected["sequence_id"].tolist() != target["sample_id"].astype(str).str.zfill(6).tolist():
            raise ValueError(f"frozen target order mismatch for {direction}")
        column = f"prediction__{fitted.model_name}"
        delta = float(np.max(np.abs(expected[column].to_numpy(float) - target_prediction)))
        if delta > 5e-10:
            raise ValueError(
                f"refitted prediction mismatch for {direction}/{fitted.model_name}: {delta}"
            )
        for row in rows:
            row["frozen_prediction_max_abs_difference"] = delta
        models[(fitted.source_pool, fitted.model_name)] = fitted
        coefficient_rows.extend(rows)
    return models, coefficient_rows


def deterministic_library(
    target: pd.DataFrame, direction: str
) -> pd.DataFrame:
    eligible = target.loc[target["sequence"].map(exact_gc_hp_eligible)].copy()
    if len(eligible) < PUBLIC_LIBRARY_SIZE:
        raise ValueError(
            f"{direction}: only {len(eligible)} exact-GC/HP eligible public sequences"
        )
    eligible["library_key"] = eligible["sequence_sha256"].map(
        lambda value: hashlib.sha256(
            f"paper2-choice-public-v1|{direction}|{value}".encode("ascii")
        ).hexdigest()
    )
    library = (
        eligible.sort_values(["library_key", "sequence_sha256"], kind="mergesort")
        .iloc[:PUBLIC_LIBRARY_SIZE]
        .reset_index(drop=True)
    )
    if library["sequence"].duplicated().any() or library["sequence_sha256"].duplicated().any():
        raise ValueError(f"{direction}: duplicate public codebook sequence")
    library["base_rank"] = np.arange(PUBLIC_LIBRARY_SIZE, dtype=int)
    return library


def bootstrap_mean_interval(
    values: np.ndarray, replicates: int, seed: int
) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    estimates = np.empty(replicates, dtype=float)
    batch_size = 128
    for start in range(0, replicates, batch_size):
        stop = min(replicates, start + batch_size)
        indices = rng.integers(0, len(values), size=(stop - start, len(values)))
        estimates[start:stop] = values[indices].mean(axis=1)
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def public_fiber_analysis(
    pcr: pd.DataFrame,
    models: dict[tuple[str, str], FrozenRidge],
    bootstrap: int,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    summary_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    selected_rows: list[dict[str, Any]] = []
    roundtrip_rows: list[dict[str, Any]] = []

    for direction_index, (source_pool, target_pool) in enumerate(
        (("GCall", "GCfix"), ("GCfix", "GCall"))
    ):
        direction = f"{source_pool}_to_{target_pool}"
        target = pcr.loc[pcr["pool"].eq(target_pool)].reset_index(drop=True)
        library = deterministic_library(target, direction)
        for model_name in MODEL_NAMES:
            library[f"prediction__{model_name}"] = models[(source_pool, model_name)].predict(
                library
            )

        reverse_lookup = dict(zip(library["sequence_sha256"], library["base_rank"], strict=True))
        library_sd = float(library["eff"].std(ddof=0))
        if library_sd <= 0:
            raise ValueError(f"{direction}: zero public outcome variance")

        for selector_bits in PUBLIC_SELECTOR_BITS:
            width = 1 << selector_bits
            fibers = PUBLIC_LIBRARY_SIZE // width
            outcomes = library["eff"].to_numpy(float).reshape(fibers, width)
            events = library["label"].to_numpy(int).reshape(fibers, width)
            fiber_mean_outcome = outcomes.mean(axis=1)
            fiber_mean_event = events.mean(axis=1)
            model_selected_outcomes: dict[str, np.ndarray] = {}

            for model_index, model_name in enumerate(MODEL_NAMES):
                prediction = library[f"prediction__{model_name}"].to_numpy(float).reshape(
                    fibers, width
                )
                choices = deterministic.argmax_smallest(prediction, axis=1)
                stability = deterministic.choice_stability_audit(prediction, "max")
                selected_outcome = outcomes[np.arange(fibers), choices]
                selected_event = events[np.arange(fibers), choices]
                outcome_gain = selected_outcome - fiber_mean_outcome
                event_change = selected_event - fiber_mean_event
                gain_ci = bootstrap_mean_interval(
                    outcome_gain,
                    bootstrap,
                    BASE_SEED + direction_index * 100_000 + selector_bits * 1000 + model_index,
                )
                event_ci = bootstrap_mean_interval(
                    event_change,
                    bootstrap,
                    BASE_SEED
                    + 50_000
                    + direction_index * 100_000
                    + selector_bits * 1000
                    + model_index,
                )
                mean_gain = float(outcome_gain.mean())
                model_selected_outcomes[model_name] = selected_outcome
                selected_base_ranks = np.arange(fibers, dtype=int) * width + choices
                decoded_payloads = selected_base_ranks >> selector_bits
                expected_payloads = np.arange(fibers, dtype=int)
                rerank = np.array(
                    [
                        reverse_lookup[library.iloc[int(rank)]["sequence_sha256"]]
                        for rank in selected_base_ranks
                    ],
                    dtype=int,
                )
                roundtrip_passes = int(
                    np.array_equal(rerank, selected_base_ranks)
                    and np.array_equal(decoded_payloads, expected_payloads)
                )

                def encode_public(payload: int) -> int:
                    if payload < 0 or payload >= fibers:
                        raise ValueError("public payload rank outside fixed domain")
                    return int(selected_base_ranks[payload])

                def decode_public_payload(sequence_hash: str) -> int:
                    if sequence_hash not in reverse_lookup:
                        raise ValueError("sequence absent from frozen public library")
                    base_rank = int(reverse_lookup[sequence_hash])
                    payload = base_rank >> selector_bits
                    if payload >= fibers:
                        raise ValueError("sequence rank lies outside the public payload domain")
                    return payload

                def verify_public_canonical(sequence_hash: str) -> bool:
                    payload = decode_public_payload(sequence_hash)
                    return encode_public(payload) == int(reverse_lookup[sequence_hash])

                for payload in expected_payloads:
                    base_rank = encode_public(int(payload))
                    sequence_hash = library.iloc[base_rank]["sequence_sha256"]
                    if decode_public_payload(sequence_hash) != payload:
                        roundtrip_passes = 0
                    if not verify_public_canonical(sequence_hash):
                        roundtrip_passes = 0

                try:
                    encode_public(fibers)
                    out_of_domain_rejected = False
                except ValueError:
                    out_of_domain_rejected = True
                try:
                    decode_public_payload("0" * 64)
                    unknown_sequence_rejected = False
                except ValueError:
                    unknown_sequence_rejected = True
                if selector_bits == 0:
                    noncanonical_payload_accepted: bool | str = "not_applicable_r0"
                    noncanonical_verifier_rejected: bool | str = "not_applicable_r0"
                else:
                    selected_rank_set = set(map(int, selected_base_ranks))
                    noncanonical_rank = next(
                        rank for rank in range(PUBLIC_LIBRARY_SIZE) if rank not in selected_rank_set
                    )
                    noncanonical_hash = library.iloc[noncanonical_rank]["sequence_sha256"]
                    try:
                        decoded_noncanonical = decode_public_payload(noncanonical_hash)
                        noncanonical_payload_accepted = (
                            decoded_noncanonical == (noncanonical_rank >> selector_bits)
                        )
                    except ValueError:
                        noncanonical_payload_accepted = False
                    noncanonical_verifier_rejected = not verify_public_canonical(
                        noncanonical_hash
                    )
                summary_rows.append(
                    {
                        "direction": direction,
                        "source_pool": source_pool,
                        "target_pool": target_pool,
                        "selector_model": model_name,
                        "selector_bits": selector_bits,
                        "candidates_per_payload": width,
                        "base_library_sequences": PUBLIC_LIBRARY_SIZE,
                        "payloads": fibers,
                        "base_payload_bits": PUBLIC_BASE_BITS,
                        "payload_bits_after_selection": PUBLIC_BASE_BITS - selector_bits,
                        "payload_bits_per_variable_nt": (PUBLIC_BASE_BITS - selector_bits)
                        / GENERATED_LENGTH,
                        "target_library_mean_outcome": float(library["eff"].mean()),
                        "selected_mean_outcome": float(selected_outcome.mean()),
                        "paired_mean_outcome_gain": mean_gain,
                        "paired_gain_ci_2p5": gain_ci[0],
                        "paired_gain_ci_97p5": gain_ci[1],
                        "standardized_paired_gain": mean_gain / library_sd,
                        "target_library_event_rate": float(library["label"].mean()),
                        "selected_event_rate": float(selected_event.mean()),
                        "paired_event_rate_change": float(event_change.mean()),
                        "event_change_ci_2p5": event_ci[0],
                        "event_change_ci_97p5": event_ci[1],
                        "fiber_score_variance_mean": float(
                            prediction.var(axis=1, ddof=0).mean()
                        ),
                        "fiber_outcome_variance_mean": float(
                            outcomes.var(axis=1, ddof=0).mean()
                        ),
                        "negative_gain_fraction": float((outcome_gain < 0).mean()),
                        "zero_gain_fraction": float((outcome_gain == 0).mean()),
                        "positive_gain_fraction": float((outcome_gain > 0).mean()),
                        "score_decimal_places": deterministic.SCORE_DECIMAL_PLACES,
                        **stability,
                        "roundtrip_payloads": fibers,
                        "roundtrip_failures": 0 if roundtrip_passes else fibers,
                        "independent_unit": "one deterministic candidate fiber",
                        "evidence_boundary": "retrospective held-out public sequence library; not prospective emitted-sequence validation",
                    }
                )

                if selector_bits in PRIMARY_SELECTOR_BITS:
                    sample_payloads = sorted(
                        set([0, 1, fibers // 3, fibers // 2, max(0, fibers - 2), fibers - 1])
                    )
                    for payload in sample_payloads:
                        base_rank = int(selected_base_ranks[payload])
                        row = library.iloc[base_rank]
                        selected_rows.append(
                            {
                                "direction": direction,
                                "selector_model": model_name,
                                "selector_bits": selector_bits,
                                "payload_rank": payload,
                                "choice_index": int(choices[payload]),
                                "base_rank": base_rank,
                                "sequence_sha256": row["sequence_sha256"],
                                "observed_relative_efficiency": float(row["eff"]),
                                "observed_low_efficiency_label": int(row["label"]),
                                "source_only_prediction": float(
                                    row[f"prediction__{model_name}"]
                                ),
                                "rerank": int(reverse_lookup[row["sequence_sha256"]]),
                                "decoded_payload": int(base_rank >> selector_bits),
                            }
                        )
                roundtrip_rows.append(
                    {
                        "codebook": f"public_{direction}",
                        "selector_model": model_name,
                        "selector_bits": selector_bits,
                        "payload_domain": fibers,
                        "tested_payloads": fibers,
                        "exact_roundtrips": fibers if roundtrip_passes else 0,
                        "roundtrip_failures": 0 if roundtrip_passes else fibers,
                        "duplicate_sequences": int(library["sequence_sha256"].duplicated().sum()),
                        "out_of_domain_rank_rejected": out_of_domain_rejected,
                        "unknown_sequence_rejected": unknown_sequence_rejected,
                        "noncanonical_candidate_accepted_by_payload_decoder": noncanonical_payload_accepted,
                        "noncanonical_candidate_rejected_by_optional_canonical_verifier": noncanonical_verifier_rejected,
                    }
                )

            for comparison_index, (extended, baseline_name) in enumerate(
                (("P5_combined_context", "P2_assay_context"), ("P5_combined_context", "P0_composition"))
            ):
                difference = (
                    model_selected_outcomes[extended] - model_selected_outcomes[baseline_name]
                )
                ci = bootstrap_mean_interval(
                    difference,
                    bootstrap,
                    BASE_SEED
                    + 90_000
                    + direction_index * 100_000
                    + selector_bits * 1000
                    + comparison_index,
                )
                comparison_rows.append(
                    {
                        "direction": direction,
                        "selector_bits": selector_bits,
                        "candidates_per_payload": width,
                        "comparison": f"{extended}_minus_{baseline_name}",
                        "paired_selected_outcome_difference": float(difference.mean()),
                        "difference_ci_2p5": ci[0],
                        "difference_ci_97p5": ci[1],
                        "standardized_difference": float(difference.mean()) / library_sd,
                        "fibers": fibers,
                        "independent_unit": "one deterministic candidate fiber",
                    }
                )

    return summary_rows, comparison_rows, selected_rows, roundtrip_rows


def public_seed_sensitivity(
    pcr: pd.DataFrame,
    models: dict[tuple[str, str], FrozenRidge],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Stress-test outcome-blind library selection and fiber assignment.

    Seeds are algorithmic perturbations, not independent experimental units.
    They are summarized as sensitivity ranges and are never pooled to narrow a
    biological confidence interval.
    """

    rows: list[dict[str, Any]] = []
    for source_pool, target_pool in (("GCall", "GCfix"), ("GCfix", "GCall")):
        direction = f"{source_pool}_to_{target_pool}"
        target = pcr.loc[pcr["pool"].eq(target_pool)].reset_index(drop=True)
        eligible = target.loc[target["sequence"].map(exact_gc_hp_eligible)].copy()
        for model_name in MODEL_NAMES:
            eligible[f"prediction__{model_name}"] = models[(source_pool, model_name)].predict(
                eligible
            )
        for seed_index in range(32):
            library = eligible.copy()
            library["seed_key"] = library["sequence_sha256"].map(
                lambda value: hashlib.sha256(
                    f"paper2-choice-public-seed-v1|{seed_index}|{direction}|{value}".encode(
                        "ascii"
                    )
                ).hexdigest()
            )
            library = (
                library.sort_values(["seed_key", "sequence_sha256"], kind="mergesort")
                .iloc[:PUBLIC_LIBRARY_SIZE]
                .reset_index(drop=True)
            )
            outcome_sd = float(library["eff"].std(ddof=0))
            for selector_bits in PRIMARY_SELECTOR_BITS:
                width = 1 << selector_bits
                fibers = PUBLIC_LIBRARY_SIZE // width
                outcomes = library["eff"].to_numpy(float).reshape(fibers, width)
                fiber_mean = outcomes.mean(axis=1)
                selected: dict[str, np.ndarray] = {}
                for model_name in MODEL_NAMES:
                    scores = library[f"prediction__{model_name}"].to_numpy(float).reshape(
                        fibers, width
                    )
                    choices = deterministic.argmax_smallest(scores, axis=1)
                    chosen = outcomes[np.arange(fibers), choices]
                    gains = chosen - fiber_mean
                    stability = deterministic.choice_stability_audit(scores, "max")
                    selected[model_name] = chosen
                    rows.append(
                        {
                            "direction": direction,
                            "seed_index": seed_index,
                            "selector_bits": selector_bits,
                            "selector_model": model_name,
                            "fibers": fibers,
                            "standardized_paired_gain": float(
                                gains.mean() / outcome_sd
                            ),
                            "negative_gain_fraction": float((gains < 0).mean()),
                            "zero_gain_fraction": float((gains == 0).mean()),
                            "positive_gain_fraction": float((gains > 0).mean()),
                            "fiber_score_variance_mean": float(
                                scores.var(axis=1, ddof=0).mean()
                            ),
                            "choices_changed_by_quantization": int(
                                stability["choices_changed_by_quantization"]
                            ),
                            "quantized_extremum_tie_fibers": int(
                                stability["quantized_extremum_tie_fibers"]
                            ),
                            "minimum_raw_top_gap": float(
                                stability["minimum_raw_top_gap"]
                            ),
                            "algorithmic_sensitivity_boundary": "hash seed is not an independent experimental replicate",
                        }
                    )
                rows.append(
                    {
                        "direction": direction,
                        "seed_index": seed_index,
                        "selector_bits": selector_bits,
                        "selector_model": "P5_minus_P2_selected_outcome",
                        "fibers": fibers,
                        "standardized_paired_gain": float(
                            (
                                selected["P5_combined_context"]
                                - selected["P2_assay_context"]
                            ).mean()
                            / outcome_sd
                        ),
                        "negative_gain_fraction": float(
                            (
                                selected["P5_combined_context"]
                                - selected["P2_assay_context"]
                                < 0
                            ).mean()
                        ),
                        "zero_gain_fraction": float(
                            (
                                selected["P5_combined_context"]
                                - selected["P2_assay_context"]
                                == 0
                            ).mean()
                        ),
                        "positive_gain_fraction": float(
                            (
                                selected["P5_combined_context"]
                                - selected["P2_assay_context"]
                                > 0
                            ).mean()
                        ),
                        "fiber_score_variance_mean": "",
                        "choices_changed_by_quantization": "",
                        "quantized_extremum_tie_fibers": "",
                        "minimum_raw_top_gap": "",
                        "algorithmic_sensitivity_boundary": "hash seed is not an independent experimental replicate",
                    }
                )

    frame = pd.DataFrame(rows)
    summary_rows: list[dict[str, Any]] = []
    for keys, group in frame.groupby(
        ["direction", "selector_bits", "selector_model"], sort=True
    ):
        direction, selector_bits, selector_model = keys
        values = group["standardized_paired_gain"].to_numpy(float)
        summary_rows.append(
            {
                "direction": direction,
                "selector_bits": int(selector_bits),
                "selector_model": selector_model,
                "seeds": len(values),
                "minimum_standardized_gain": float(values.min()),
                "median_standardized_gain": float(np.median(values)),
                "maximum_standardized_gain": float(values.max()),
                "positive_seed_fraction": float((values > 0).mean()),
                "interpretation": "algorithmic mapping sensitivity; not a confidence interval or biological replication",
            }
        )
    return rows, summary_rows


class GCHomopolymerCodec:
    """Exact lexicographic codec for fixed length, GC interval and HP bound."""

    def __init__(self, length: int, gc_min: int, gc_max: int, max_hp: int) -> None:
        self.length = length
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.max_hp = max_hp
        self._count_cached = lru_cache(maxsize=None)(self._count_uncached)

    def children(
        self, pos: int, last: str, run: int, gc_count: int
    ) -> list[tuple[str, int, int]]:
        remaining_after = self.length - pos - 1
        children: list[tuple[str, int, int]] = []
        for nucleotide in "ACGT":
            new_run = run + 1 if nucleotide == last else 1
            if new_run > self.max_hp:
                continue
            new_gc = gc_count + int(nucleotide in "GC")
            if new_gc > self.gc_max or new_gc + remaining_after < self.gc_min:
                continue
            children.append((nucleotide, new_run, new_gc))
        return children

    def _count_uncached(self, pos: int, last: str, run: int, gc_count: int) -> int:
        if pos == self.length:
            return int(self.gc_min <= gc_count <= self.gc_max)
        return sum(
            self._count_cached(pos + 1, nucleotide, new_run, new_gc)
            for nucleotide, new_run, new_gc in self.children(pos, last, run, gc_count)
        )

    def count(self, pos: int = 0, last: str = "", run: int = 0, gc_count: int = 0) -> int:
        return self._count_cached(pos, last, run, gc_count)

    def unrank(self, rank: int) -> str:
        total = self.count()
        if rank < 0 or rank >= total:
            raise ValueError("rank outside exact GC/homopolymer language")
        pos, last, run, gc_count = 0, "", 0, 0
        sequence: list[str] = []
        while pos < self.length:
            chosen = None
            for nucleotide, new_run, new_gc in self.children(pos, last, run, gc_count):
                mass = self.count(pos + 1, nucleotide, new_run, new_gc)
                if rank < mass:
                    chosen = (nucleotide, new_run, new_gc)
                    break
                rank -= mass
            if chosen is None:
                raise AssertionError("positive completion mass did not select a branch")
            nucleotide, run, gc_count = chosen
            last = nucleotide
            sequence.append(nucleotide)
            pos += 1
        return "".join(sequence)

    def rank(self, sequence: str) -> int:
        if len(sequence) != self.length or set(sequence) - set("ACGT"):
            raise ValueError("sequence outside configured DNA length/alphabet")
        rank = 0
        pos, last, run, gc_count = 0, "", 0, 0
        for observed in sequence:
            chosen = None
            for nucleotide, new_run, new_gc in self.children(pos, last, run, gc_count):
                mass = self.count(pos + 1, nucleotide, new_run, new_gc)
                if nucleotide == observed:
                    if mass == 0:
                        raise ValueError("sequence enters zero-completion branch")
                    chosen = (new_run, new_gc)
                    break
                rank += mass
            if chosen is None:
                raise ValueError(f"sequence violates constraint at position {pos}")
            run, gc_count = chosen
            last = observed
            pos += 1
        if not self.gc_min <= gc_count <= self.gc_max:
            raise ValueError("sequence violates terminal GC interval")
        return rank


def affine_constants(modulus: int, namespace_index: int = 0) -> tuple[int, int, int]:
    """Return a deterministic affine permutation of the complete language.

    The previous implementation used modulus ``2**K`` and therefore permuted
    only a lexicographic prefix of the exact language.  Here the modulus is the
    complete language cardinality ``N``.  Incrementing the hash-derived
    multiplier until it is coprime to ``N`` makes the map a permutation on
    ``[0,N)``; restricting its input to ``q < 2**K`` remains injective while
    dispersing the used ranks through the complete constrained language.
    """

    if modulus < 2:
        raise ValueError("affine modulus must exceed one")
    namespace = f"paper2-choice-codec-affine-v2|namespace={namespace_index}".encode(
        "ascii"
    )
    multiplier = int.from_bytes(hashlib.sha256(namespace + b"|multiplier").digest(), "big") % modulus
    if multiplier == 0:
        multiplier = 1
    while math.gcd(multiplier, modulus) != 1:
        multiplier += 1
        if multiplier >= modulus:
            multiplier = 1
    offset = int.from_bytes(hashlib.sha256(namespace + b"|offset").digest(), "big") % modulus
    inverse = pow(multiplier, -1, modulus)
    return multiplier, offset, inverse


def affine_record(
    codec: GCHomopolymerCodec, fixed_bits: int, namespace_index: int = 0
) -> dict[str, int]:
    modulus = codec.count()
    multiplier, offset, inverse = affine_constants(modulus, namespace_index)
    return {
        "mapping_namespace_index": namespace_index,
        "language_modulus": modulus,
        "payload_domain_size": 1 << fixed_bits,
        "affine_multiplier": multiplier,
        "affine_offset": offset,
        "affine_inverse": inverse,
    }


def decode_payload_candidate(
    codec: GCHomopolymerCodec,
    sequence: str,
    fixed_bits: int,
    selector_bits: int,
    affine: dict[str, int],
) -> tuple[int, int, int, int]:
    """Decode any valid member of a declared choice fiber.

    This is the payload interface.  It does not evaluate the assay score and
    does not require the sequence to be the representative selected by the
    encoder.  The optional canonical verifier below is a separate interface.
    """

    rank = codec.rank(sequence)
    modulus = int(affine["language_modulus"])
    logical_rank = (
        int(affine["affine_inverse"]) * (rank - int(affine["affine_offset"]))
    ) % modulus
    payload_domain = 1 << fixed_bits
    if logical_rank >= payload_domain:
        raise ValueError("sequence maps outside the fixed power-of-two payload domain")
    payload = logical_rank >> selector_bits
    choice = logical_rank & ((1 << selector_bits) - 1)
    return payload, choice, logical_rank, rank


def verify_canonical_candidate(
    codec: GCHomopolymerCodec,
    sequence: str,
    fixed_bits: int,
    selector_bits: int,
    affine: dict[str, int],
    score,
) -> bool:
    """Optionally check whether ``sequence`` is the deterministic emitted representative."""

    payload, _choice, _logical_rank, _rank = decode_payload_candidate(
        codec, sequence, fixed_bits, selector_bits, affine
    )
    candidates: list[tuple[float, int, str]] = []
    modulus = int(affine["language_modulus"])
    for choice in range(1 << selector_bits):
        logical_rank = (payload << selector_bits) | choice
        physical_rank = (
            int(affine["affine_multiplier"]) * logical_rank
            + int(affine["affine_offset"])
        ) % modulus
        candidate = codec.unrank(physical_rank)
        candidates.append((float(score(candidate)), choice, candidate))
    score_keys = deterministic.quantized_score_keys([item[0] for item in candidates])
    # Maximum declared score key; np.argmax returns the smallest choice index on ties.
    expected = candidates[int(np.argmax(score_keys))][2]
    return sequence == expected


def deterministic_payloads(limit: int, count: int, selector_bits: int) -> list[int]:
    if count > limit:
        raise ValueError("requested more payload samples than payload domain")
    fixed = {0, min(1, limit - 1), limit // 3, limit // 2, max(0, limit - 2), limit - 1}
    payloads = {value for value in fixed if 0 <= value < limit}
    index = 0
    while len(payloads) < count:
        digest = hashlib.sha256(
            f"paper2-choice-payload-v1|r={selector_bits}|i={index}".encode("ascii")
        ).digest()
        payloads.add(int.from_bytes(digest, "big") % limit)
        index += 1
    return sorted(payloads)[:count]


def generated_candidates(
    codec: GCHomopolymerCodec,
    fixed_bits: int,
    payload_samples: int,
    namespace_index: int = 0,
    selector_bits_values: Iterable[int] = GENERATED_SELECTOR_BITS,
) -> tuple[pd.DataFrame, dict[str, int]]:
    affine = affine_record(codec, fixed_bits, namespace_index)
    modulus = int(affine["language_modulus"])
    multiplier = int(affine["affine_multiplier"])
    offset = int(affine["affine_offset"])
    rows: list[dict[str, Any]] = []
    for selector_bits in selector_bits_values:
        payload_limit = 1 << (fixed_bits - selector_bits)
        payloads = deterministic_payloads(payload_limit, payload_samples, selector_bits)
        for sample_index, payload in enumerate(payloads):
            for choice in range(1 << selector_bits):
                logical_rank = (payload << selector_bits) | choice
                physical_rank = (multiplier * logical_rank + offset) % modulus
                sequence = codec.unrank(physical_rank)
                rerank = codec.rank(sequence)
                recovered_payload, recovered_choice, recovered_logical, _ = (
                    decode_payload_candidate(
                        codec, sequence, fixed_bits, selector_bits, affine
                    )
                )
                if (
                    rerank != physical_rank
                    or recovered_payload != payload
                    or recovered_choice != choice
                ):
                    raise AssertionError("generated candidate rank/unrank failure")
                rows.append(
                    {
                        "selector_bits": selector_bits,
                        "mapping_namespace_index": namespace_index,
                        "payload_sample_index": sample_index,
                        "payload_rank": payload,
                        "choice_index": choice,
                        "logical_rank": logical_rank,
                        "physical_rank": physical_rank,
                        "sequence": sequence,
                        "sequence_sha256": sha256_text(sequence),
                        "rerank": rerank,
                        "decoded_payload": recovered_payload,
                        "decoded_choice": recovered_choice,
                    }
                )
    frame = pd.DataFrame(rows)
    if frame.duplicated(["selector_bits", "payload_rank", "choice_index"]).any():
        raise ValueError("duplicate logical generated candidate")
    return frame, affine


def generated_feature_row(row: tuple[int, str]) -> dict[str, Any]:
    index, sequence = row
    result = base.pcr_feature_row("generated", str(index), sequence)
    result["generated_row_index"] = index
    return result


def generated_feature_frame(candidates: pd.DataFrame, jobs: int) -> pd.DataFrame:
    unique = candidates[["sequence_sha256", "sequence"]].drop_duplicates().reset_index(drop=True)
    rows = Parallel(n_jobs=jobs, prefer="processes", batch_size=32)(
        delayed(generated_feature_row)((index, row.sequence))
        for index, row in enumerate(unique.itertuples(index=False))
    )
    features = pd.DataFrame(rows).drop(columns=["pool", "sample_id"])
    features = unique.merge(
        features,
        on="sequence_sha256",
        how="inner",
        validate="one_to_one",
    )
    if len(features) != len(unique):
        raise ValueError("generated feature merge changed row count")
    return features


def generated_fiber_analysis(
    candidates: pd.DataFrame,
    features: pd.DataFrame,
    models: dict[tuple[str, str], FrozenRidge],
    pcr: pd.DataFrame,
    bootstrap: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    frame = candidates.merge(
        features.drop(columns=["sequence"]),
        on="sequence_sha256",
        how="inner",
        validate="many_to_one",
    )
    for source_pool in ("GCall", "GCfix"):
        for model_name in MODEL_NAMES:
            frame[f"prediction__{source_pool}__{model_name}"] = models[
                (source_pool, model_name)
            ].predict(frame)

    summary_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    rejection_rows: list[dict[str, Any]] = []
    source_thresholds: dict[tuple[str, str], float] = {}
    for source_pool in ("GCall", "GCfix"):
        source = pcr.loc[pcr["pool"].eq(source_pool)].reset_index(drop=True)
        for model_name in MODEL_NAMES:
            prediction = models[(source_pool, model_name)].predict(source)
            selected_n = int(math.ceil(0.25 * len(prediction)))
            source_thresholds[(source_pool, model_name)] = float(
                np.sort(prediction)[-selected_n]
            )

    for selector_bits in sorted(
        int(value) for value in candidates["selector_bits"].unique()
    ):
        subset = frame.loc[frame["selector_bits"].eq(selector_bits)].copy()
        width = 1 << selector_bits
        payloads = subset["payload_rank"].drop_duplicates().sort_values().tolist()
        subset = subset.sort_values(["payload_rank", "choice_index"], kind="mergesort")
        if len(subset) != len(payloads) * width:
            raise ValueError("generated fiber width mismatch")

        for source_index, source_pool in enumerate(("GCall", "GCfix")):
            other_pool = "GCfix" if source_pool == "GCall" else "GCall"
            for model_index, model_name in enumerate(MODEL_NAMES):
                selector_column = f"prediction__{source_pool}__{model_name}"
                selector_scores = subset[selector_column].to_numpy(float).reshape(
                    len(payloads), width
                )
                choices = deterministic.argmax_smallest(selector_scores, axis=1)
                selector_stability = deterministic.choice_stability_audit(
                    selector_scores, "max"
                )
                selected_offsets = np.arange(len(payloads)) * width + choices
                selected = subset.iloc[selected_offsets].reset_index(drop=True)
                own_selected = selector_scores[np.arange(len(payloads)), choices]
                own_gain = own_selected - selector_scores.mean(axis=1)
                own_ci = bootstrap_mean_interval(
                    own_gain,
                    bootstrap,
                    BASE_SEED
                    + 200_000
                    + selector_bits * 10_000
                    + source_index * 1000
                    + model_index,
                )
                other_column = f"prediction__{other_pool}__{model_name}"
                other_scores = subset[other_column].to_numpy(float).reshape(
                    len(payloads), width
                )
                other_selected = other_scores[np.arange(len(payloads)), choices]
                other_gain = other_selected - other_scores.mean(axis=1)
                other_ci = bootstrap_mean_interval(
                    other_gain,
                    bootstrap,
                    BASE_SEED
                    + 250_000
                    + selector_bits * 10_000
                    + source_index * 1000
                    + model_index,
                )
                selected_rerank_failures = int(
                    (selected["physical_rank"] != selected["rerank"]).sum()
                    + (selected["payload_rank"] != selected["decoded_payload"]).sum()
                )
                summary_rows.append(
                    {
                        "selector_source_pool": source_pool,
                        "cross_evaluator_pool": other_pool,
                        "selector_model": model_name,
                        "selector_bits": selector_bits,
                        "candidates_per_payload": width,
                        "sampled_payloads": len(payloads),
                        "base_fixed_payload_bits": 213,
                        "payload_bits_after_selection": 213 - selector_bits,
                        "payload_bits_per_variable_nt": (213 - selector_bits)
                        / GENERATED_LENGTH,
                        "own_model_mean_gain_vs_fiber": float(own_gain.mean()),
                        "own_model_gain_minimum": float(own_gain.min()),
                        "own_model_gain_q25": float(np.quantile(own_gain, 0.25)),
                        "own_model_gain_median": float(np.median(own_gain)),
                        "own_model_gain_q75": float(np.quantile(own_gain, 0.75)),
                        "own_model_gain_maximum": float(own_gain.max()),
                        "own_gain_ci_2p5": own_ci[0],
                        "own_gain_ci_97p5": own_ci[1],
                        "cross_model_mean_gain_vs_fiber": float(other_gain.mean()),
                        "cross_model_gain_minimum": float(other_gain.min()),
                        "cross_model_gain_q25": float(np.quantile(other_gain, 0.25)),
                        "cross_model_gain_median": float(np.median(other_gain)),
                        "cross_model_gain_q75": float(np.quantile(other_gain, 0.75)),
                        "cross_model_gain_maximum": float(other_gain.max()),
                        "cross_gain_ci_2p5": other_ci[0],
                        "cross_gain_ci_97p5": other_ci[1],
                        "own_model_fiber_score_variance_mean": float(
                            selector_scores.var(axis=1, ddof=0).mean()
                        ),
                        "own_model_fiber_score_variance_median": float(
                            np.median(selector_scores.var(axis=1, ddof=0))
                        ),
                        "cross_model_fiber_score_variance_mean": float(
                            other_scores.var(axis=1, ddof=0).mean()
                        ),
                        "cross_model_fiber_score_variance_median": float(
                            np.median(other_scores.var(axis=1, ddof=0))
                        ),
                        "own_model_negative_gain_fraction": float((own_gain < 0).mean()),
                        "own_model_zero_gain_fraction": float((own_gain == 0).mean()),
                        "fibers_improved_in_own_model_fraction": float((own_gain > 0).mean()),
                        "cross_model_negative_gain_fraction": float((other_gain < 0).mean()),
                        "cross_model_zero_gain_fraction": float((other_gain == 0).mean()),
                        "fibers_improved_in_cross_model_fraction": float(
                            (other_gain > 0).mean()
                        ),
                        "score_decimal_places": deterministic.SCORE_DECIMAL_PLACES,
                        "choices_changed_by_quantization": int(
                            selector_stability["choices_changed_by_quantization"]
                        ),
                        "quantized_extremum_tie_fibers": int(
                            selector_stability["quantized_extremum_tie_fibers"]
                        ),
                        "minimum_raw_top_gap": float(
                            selector_stability["minimum_raw_top_gap"]
                        ),
                        "selected_roundtrip_failures": selected_rerank_failures,
                        "inference_boundary": "bootstrap describes deterministic payload-sample variability, not biological replication",
                    }
                )

                if model_name == "P5_combined_context" and selector_bits in PRIMARY_SELECTOR_BITS:
                    for payload_index in sorted(
                        set([0, 1, len(payloads) // 3, len(payloads) // 2, len(payloads) - 2, len(payloads) - 1])
                    ):
                        row = selected.iloc[payload_index]
                        sample_rows.append(
                            {
                                "selector_source_pool": source_pool,
                                "selector_model": model_name,
                                "selector_bits": selector_bits,
                                "payload_rank": int(row["payload_rank"]),
                                "choice_index": int(row["choice_index"]),
                                "physical_rank": int(row["physical_rank"]),
                                "sequence": row["sequence"],
                                "sequence_sha256": row["sequence_sha256"],
                                "gc_count": int(row["sequence"].count("G") + row["sequence"].count("C")),
                                "max_homopolymer": int(row["max_homopolymer"]),
                                "variable_candidate_pairs": int(row["variable_candidate_pairs"]),
                                "variable_weighted_pairs": int(row["variable_weighted_pairs"]),
                                "source_model_score": float(row[selector_column]),
                                "cross_model_score": float(row[other_column]),
                                "rerank": int(row["rerank"]),
                                "decoded_payload": int(row["decoded_payload"]),
                            }
                        )

                threshold = source_thresholds[(source_pool, model_name)]
                accepted = selector_scores >= threshold
                first_choice = np.argmax(accepted, axis=1)
                success = accepted.any(axis=1)
                attempts = np.where(success, first_choice + 1, width)
                rejection_rows.append(
                    {
                        "selector_source_pool": source_pool,
                        "selector_model": model_name,
                        "selector_bits": selector_bits,
                        "max_attempts": width,
                        "source_public_top25_threshold": threshold,
                        "sampled_payloads": len(payloads),
                        "successful_payloads": int(success.sum()),
                        "failed_payloads": int((~success).sum()),
                        "failure_fraction": float((~success).mean()),
                        "mean_attempts_including_failures_at_cap": float(attempts.mean()),
                        "choice_codec_failure_fraction": 0.0,
                        "boundary": "threshold transfer to generated sequences is descriptive and may be out of distribution",
                    }
                )

    return summary_rows, sample_rows, rejection_rows


def generated_namespace_sensitivity(
    codec: GCHomopolymerCodec,
    fixed_bits: int,
    models: dict[tuple[str, str], FrozenRidge],
    pcr: pd.DataFrame,
    jobs: int,
    bootstrap: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Audit generated-score stability over outcome-blind full-language mappings."""

    frames: list[pd.DataFrame] = []
    for namespace_index in GENERATED_NAMESPACE_INDICES:
        frame, _affine = generated_candidates(
            codec,
            fixed_bits,
            GENERATED_NAMESPACE_PAYLOADS,
            namespace_index=namespace_index,
            selector_bits_values=PRIMARY_SELECTOR_BITS,
        )
        frames.append(frame)
    all_candidates = pd.concat(frames, ignore_index=True)
    all_features = generated_feature_frame(all_candidates, jobs)

    rows: list[dict[str, Any]] = []
    for frame in frames:
        namespace_index = int(frame["mapping_namespace_index"].iloc[0])
        summaries, _samples, _rejections = generated_fiber_analysis(
            frame, all_features, models, pcr, bootstrap
        )
        for row in summaries:
            if (
                row["selector_model"] == "P5_combined_context"
                and int(row["selector_bits"]) in PRIMARY_SELECTOR_BITS
            ):
                rows.append(
                    {
                        "mapping_namespace_index": namespace_index,
                        "selector_source_pool": row["selector_source_pool"],
                        "cross_evaluator_pool": row["cross_evaluator_pool"],
                        "selector_model": row["selector_model"],
                        "selector_bits": int(row["selector_bits"]),
                        "sampled_payloads": int(row["sampled_payloads"]),
                        "own_model_mean_gain_vs_fiber": float(
                            row["own_model_mean_gain_vs_fiber"]
                        ),
                        "cross_model_mean_gain_vs_fiber": float(
                            row["cross_model_mean_gain_vs_fiber"]
                        ),
                        "own_model_negative_gain_fraction": float(
                            row["own_model_negative_gain_fraction"]
                        ),
                        "cross_model_negative_gain_fraction": float(
                            row["cross_model_negative_gain_fraction"]
                        ),
                        "own_model_fiber_score_variance_mean": float(
                            row["own_model_fiber_score_variance_mean"]
                        ),
                        "cross_model_fiber_score_variance_mean": float(
                            row["cross_model_fiber_score_variance_mean"]
                        ),
                        "choices_changed_by_quantization": int(
                            row["choices_changed_by_quantization"]
                        ),
                        "quantized_extremum_tie_fibers": int(
                            row["quantized_extremum_tie_fibers"]
                        ),
                        "minimum_raw_top_gap": float(row["minimum_raw_top_gap"]),
                        "selected_roundtrip_failures": int(
                            row["selected_roundtrip_failures"]
                        ),
                        "interpretation": "outcome-blind full-language affine mapping sensitivity; not biological replication",
                    }
                )

    summary_rows: list[dict[str, Any]] = []
    frame = pd.DataFrame(rows)
    for keys, group in frame.groupby(
        ["selector_source_pool", "cross_evaluator_pool", "selector_bits"], sort=True
    ):
        source_pool, cross_pool, selector_bits = keys
        own = group["own_model_mean_gain_vs_fiber"].to_numpy(float)
        cross = group["cross_model_mean_gain_vs_fiber"].to_numpy(float)
        own_negative = group["own_model_negative_gain_fraction"].to_numpy(float)
        cross_negative = group["cross_model_negative_gain_fraction"].to_numpy(float)
        own_variance = group["own_model_fiber_score_variance_mean"].to_numpy(float)
        cross_variance = group["cross_model_fiber_score_variance_mean"].to_numpy(float)
        summary_rows.append(
            {
                "selector_source_pool": source_pool,
                "cross_evaluator_pool": cross_pool,
                "selector_bits": int(selector_bits),
                "mapping_namespaces": len(group),
                "payloads_per_namespace": GENERATED_NAMESPACE_PAYLOADS,
                "own_gain_minimum": float(own.min()),
                "own_gain_median": float(np.median(own)),
                "own_gain_maximum": float(own.max()),
                "own_gain_positive_namespace_fraction": float((own > 0).mean()),
                "cross_gain_minimum": float(cross.min()),
                "cross_gain_median": float(np.median(cross)),
                "cross_gain_maximum": float(cross.max()),
                "cross_gain_positive_namespace_fraction": float((cross > 0).mean()),
                "own_negative_fiber_fraction_minimum": float(own_negative.min()),
                "own_negative_fiber_fraction_median": float(np.median(own_negative)),
                "own_negative_fiber_fraction_maximum": float(own_negative.max()),
                "cross_negative_fiber_fraction_minimum": float(cross_negative.min()),
                "cross_negative_fiber_fraction_median": float(np.median(cross_negative)),
                "cross_negative_fiber_fraction_maximum": float(cross_negative.max()),
                "own_fiber_score_variance_mean_minimum": float(own_variance.min()),
                "own_fiber_score_variance_mean_median": float(np.median(own_variance)),
                "own_fiber_score_variance_mean_maximum": float(own_variance.max()),
                "cross_fiber_score_variance_mean_minimum": float(cross_variance.min()),
                "cross_fiber_score_variance_mean_median": float(np.median(cross_variance)),
                "cross_fiber_score_variance_mean_maximum": float(cross_variance.max()),
                "choices_changed_by_quantization_total": int(
                    group["choices_changed_by_quantization"].sum()
                ),
                "quantized_extremum_tie_fibers_total": int(
                    group["quantized_extremum_tie_fibers"].sum()
                ),
                "minimum_raw_top_gap_over_namespaces": float(
                    group["minimum_raw_top_gap"].min()
                ),
                "roundtrip_failures_total": int(
                    group["selected_roundtrip_failures"].sum()
                ),
                "interpretation": "algorithmic mapping sensitivity only; not a confidence interval or biological replication",
            }
        )
    return rows, summary_rows


def feature_overlap_rows(
    generated_features: pd.DataFrame,
    pcr: pd.DataFrame,
    models: dict[tuple[str, str], FrozenRidge],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    all_features = sorted(set().union(*MODEL_FEATURES.values()))
    for pool in ("GCall", "GCfix"):
        source_full = pcr.loc[pcr["pool"].eq(pool)]
        p5 = models[(pool, "P5_combined_context")]
        p5_coefficients = dict(
            zip(
                p5.feature_names,
                np.asarray(p5.pipeline.named_steps["model"].coef_).ravel(),
                strict=True,
            )
        )
        references = [
            ("full_training_pool", source_full),
            (
                "codec_eligible_training_subset",
                source_full.loc[source_full["sequence"].map(exact_gc_hp_eligible)],
            ),
        ]
        for source_reference, source in references:
            for feature in all_features:
                source_values = source[feature].to_numpy(float)
                generated_values = generated_features[feature].to_numpy(float)
                source_q01, source_q05, source_q50, source_q95, source_q99 = np.quantile(
                    source_values, [0.01, 0.05, 0.50, 0.95, 0.99]
                )
                generated_q01, generated_q05, generated_q50, generated_q95, generated_q99 = np.quantile(
                    generated_values, [0.01, 0.05, 0.50, 0.95, 0.99]
                )
                pooled_sd = math.sqrt(
                    (
                        float(source_values.var(ddof=0))
                        + float(generated_values.var(ddof=0))
                    )
                    / 2
                )
                standardized_mean_difference = (
                    0.0
                    if pooled_sd == 0
                    else float(generated_values.mean() - source_values.mean()) / pooled_sd
                )
                ks = ks_2samp(source_values, generated_values, method="asymp")
                rows.append(
                    {
                        "source_pool": pool,
                        "source_reference": source_reference,
                        "source_sequences": len(source),
                        "generated_sequences": len(generated_values),
                        "feature": feature,
                        "source_min": float(source_values.min()),
                        "source_mean": float(source_values.mean()),
                        "source_sd_population": float(source_values.std(ddof=0)),
                        "source_q01": float(source_q01),
                        "source_q05": float(source_q05),
                        "source_median": float(source_q50),
                        "source_q95": float(source_q95),
                        "source_q99": float(source_q99),
                        "source_max": float(source_values.max()),
                        "generated_min": float(generated_values.min()),
                        "generated_mean": float(generated_values.mean()),
                        "generated_sd_population": float(generated_values.std(ddof=0)),
                        "generated_q01": float(generated_q01),
                        "generated_q05": float(generated_q05),
                        "generated_median": float(generated_q50),
                        "generated_q95": float(generated_q95),
                        "generated_q99": float(generated_q99),
                        "generated_max": float(generated_values.max()),
                        "generated_within_source_minmax_fraction": float(
                            (
                                (generated_values >= source_values.min())
                                & (generated_values <= source_values.max())
                            ).mean()
                        ),
                        "generated_within_source_q01_q99_fraction": float(
                            (
                                (generated_values >= source_q01)
                                & (generated_values <= source_q99)
                            ).mean()
                        ),
                        "generated_minus_source_standardized_mean_difference": standardized_mean_difference,
                        "two_sample_ks_statistic": float(ks.statistic),
                        "two_sample_ks_p_value_descriptive_only": float(ks.pvalue),
                        "wasserstein_distance_source_units": float(
                            wasserstein_distance(source_values, generated_values)
                        ),
                        "p5_standardized_coefficient": float(
                            p5_coefficients.get(feature, 0.0)
                        ),
                        "p5_active_feature": bool(
                            abs(float(p5_coefficients.get(feature, 0.0))) > 1e-12
                        ),
                    }
                )
    return rows


def positional_nucleotide_frequency_rows(
    candidates: pd.DataFrame, pcr: pd.DataFrame
) -> list[dict[str, Any]]:
    generated = candidates[["sequence_sha256", "sequence"]].drop_duplicates()
    populations = [
        ("generated_unique_candidates", generated["sequence"].astype(str).tolist()),
        (
            "source_GCall",
            pcr.loc[pcr["pool"].eq("GCall"), "sequence"].astype(str).tolist(),
        ),
        (
            "source_GCfix",
            pcr.loc[pcr["pool"].eq("GCfix"), "sequence"].astype(str).tolist(),
        ),
    ]
    rows: list[dict[str, Any]] = []
    for population, sequences in populations:
        array = np.asarray([list(sequence) for sequence in sequences], dtype="U1")
        if array.shape != (len(sequences), GENERATED_LENGTH):
            raise ValueError(f"{population}: unexpected positional array shape")
        for position in range(GENERATED_LENGTH):
            rows.append(
                {
                    "population": population,
                    "sequences": len(sequences),
                    "position_1_based": position + 1,
                    **{
                        f"frequency_{base_name}": float((array[:, position] == base_name).mean())
                        for base_name in "ACGT"
                    },
                }
            )
    return rows


def terminal_window_frequency_rows(
    candidates: pd.DataFrame, pcr: pd.DataFrame
) -> list[dict[str, Any]]:
    generated = candidates[["sequence_sha256", "sequence"]].drop_duplicates()
    populations = [
        ("generated_unique_candidates", generated["sequence"].astype(str).tolist()),
        (
            "source_GCall",
            pcr.loc[pcr["pool"].eq("GCall"), "sequence"].astype(str).tolist(),
        ),
        (
            "source_GCfix",
            pcr.loc[pcr["pool"].eq("GCfix"), "sequence"].astype(str).tolist(),
        ),
    ]
    windows = {
        "first_base": np.arange(0, 1),
        "first_24nt": np.arange(0, 24),
        "last_24nt": np.arange(GENERATED_LENGTH - 24, GENERATED_LENGTH),
        "full_108nt": np.arange(0, GENERATED_LENGTH),
    }
    rows: list[dict[str, Any]] = []
    for population, sequences in populations:
        array = np.asarray([list(sequence) for sequence in sequences], dtype="U1")
        for window_name, positions in windows.items():
            values = array[:, positions].reshape(-1)
            rows.append(
                {
                    "population": population,
                    "sequences": len(sequences),
                    "window": window_name,
                    "positions_1_based": (
                        str(int(positions[0] + 1))
                        if len(positions) == 1
                        else f"{int(positions[0] + 1)}-{int(positions[-1] + 1)}"
                    ),
                    **{
                        f"frequency_{base_name}": float((values == base_name).mean())
                        for base_name in "ACGT"
                    },
                }
            )
    return rows


def independent_bottom_up_total(
    length: int, gc_min: int, gc_max: int, max_hp: int
) -> int:
    layer: dict[tuple[str, int, int], int] = {("", 0, 0): 1}
    for _pos in range(length):
        next_layer: dict[tuple[str, int, int], int] = {}
        for (last, run, gc_count), count in layer.items():
            for nucleotide in "ACGT":
                new_run = run + 1 if nucleotide == last else 1
                if new_run > max_hp:
                    continue
                new_gc = gc_count + int(nucleotide in "GC")
                if new_gc > gc_max:
                    continue
                key = (nucleotide, new_run, new_gc)
                next_layer[key] = next_layer.get(key, 0) + count
        layer = next_layer
    return sum(
        count
        for (_last, _run, gc_count), count in layer.items()
        if gc_min <= gc_count <= gc_max
    )


def reduced_exhaustive_audit() -> dict[str, Any]:
    length, gc_min, gc_max, max_hp = 8, 3, 5, 3
    codec = GCHomopolymerCodec(length, gc_min, gc_max, max_hp)
    accepted: list[str] = []
    for letters in product("ACGT", repeat=length):
        sequence = "".join(letters)
        gc_count = sequence.count("G") + sequence.count("C")
        if gc_min <= gc_count <= gc_max and max_homopolymer(sequence) <= max_hp:
            accepted.append(sequence)
    failures = 0
    for rank, sequence in enumerate(accepted):
        failures += int(codec.unrank(rank) != sequence)
        failures += int(codec.rank(sequence) != rank)
    return {
        "reduced_length": length,
        "reduced_gc_min": gc_min,
        "reduced_gc_max": gc_max,
        "reduced_max_homopolymer": max_hp,
        "reduced_explicit_total": len(accepted),
        "reduced_dp_total": codec.count(),
        "reduced_all_rank_checks": len(accepted),
        "reduced_rank_unrank_failures": failures,
    }


def base_language_rows(codec: GCHomopolymerCodec) -> list[dict[str, Any]]:
    total = codec.count()
    fixed_bits = total.bit_length() - 1
    independent_total = independent_bottom_up_total(
        codec.length, codec.gc_min, codec.gc_max, codec.max_hp
    )
    reduced = reduced_exhaustive_audit()
    branch_masses = {
        nucleotide: codec.count(
            1, nucleotide, 1, int(nucleotide in "GC")
        )
        for nucleotide in "ACGT"
    }
    payload_domain = 1 << fixed_bits
    legacy_covered: dict[str, int] = {}
    branch_start = 0
    for nucleotide in "ACGT":
        branch_end = branch_start + branch_masses[nucleotide]
        legacy_covered[nucleotide] = max(
            0, min(payload_domain, branch_end) - branch_start
        )
        branch_start = branch_end
    return [
        {
            "language": "exact_108nt_gc49_59_hp3",
            "length_nt": codec.length,
            "gc_min": codec.gc_min,
            "gc_max": codec.gc_max,
            "max_homopolymer": codec.max_hp,
            "exact_total": total,
            "log2_exact_total": math.log2(total),
            "fixed_payload_bits": fixed_bits,
            "fixed_payload_rate_bits_per_nt": fixed_bits / codec.length,
            "unused_exact_mass_fraction": 1.0 - (1 << fixed_bits) / total,
            "payload_domain_fraction_of_exact_language": payload_domain / total,
            "rank_dispersion_modulus": total,
            "rank_dispersion_contract": "affine permutation on complete exact language; gcd(multiplier,N)=1; inverse rejects q>=2^K",
            "first_base_branch_masses_equal": len(set(branch_masses.values())) == 1,
            "first_base_branch_mass_A": branch_masses["A"],
            "first_base_branch_mass_C": branch_masses["C"],
            "first_base_branch_mass_G": branch_masses["G"],
            "first_base_branch_mass_T": branch_masses["T"],
            "legacy_lexicographic_prefix_first_base_A_fraction": legacy_covered["A"] / payload_domain,
            "legacy_lexicographic_prefix_first_base_C_fraction": legacy_covered["C"] / payload_domain,
            "legacy_lexicographic_prefix_first_base_G_fraction": legacy_covered["G"] / payload_domain,
            "legacy_lexicographic_prefix_first_base_T_fraction": legacy_covered["T"] / payload_domain,
            "completion_cache_entries": codec._count_cached.cache_info().currsize,
            "independent_bottom_up_total": independent_total,
            "top_down_bottom_up_match": total == independent_total,
            **reduced,
            "claim_boundary": "noiseless constrained-language payload rate; not achieved physical storage density",
        }
    ]


def generated_roundtrip_summary(
    codec: GCHomopolymerCodec,
    candidates: pd.DataFrame,
    fixed_bits: int,
    affine: dict[str, int],
) -> list[dict[str, Any]]:
    candidate_failures = int(
        (candidates["physical_rank"] != candidates["rerank"]).sum()
        + (candidates["payload_rank"] != candidates["decoded_payload"]).sum()
        + (candidates["choice_index"] != candidates["decoded_choice"]).sum()
    )
    try:
        codec.unrank(codec.count())
        exact_out_of_domain_rejected = False
    except ValueError:
        exact_out_of_domain_rejected = True

    fixed_domain = 1 << fixed_bits
    outside_physical_rank = (
        int(affine["affine_multiplier"]) * fixed_domain
        + int(affine["affine_offset"])
    ) % int(affine["language_modulus"])
    outside_sequence = codec.unrank(outside_physical_rank)
    try:
        decode_payload_candidate(codec, outside_sequence, fixed_bits, 0, affine)
        payload_domain_out_of_domain_sequence_rejected = False
    except ValueError:
        payload_domain_out_of_domain_sequence_rejected = True
    try:
        decode_payload_candidate(codec, "A" * codec.length, fixed_bits, 0, affine)
        invalid_constraint_sequence_rejected = False
    except ValueError:
        invalid_constraint_sequence_rejected = True
    try:
        decode_payload_candidate(codec, "A" * (codec.length - 1), fixed_bits, 0, affine)
        wrong_length_sequence_rejected = False
    except ValueError:
        wrong_length_sequence_rejected = True
    try:
        decode_payload_candidate(
            codec, "N" + "A" * (codec.length - 1), fixed_bits, 0, affine
        )
        ambiguous_base_sequence_rejected = False
    except ValueError:
        ambiguous_base_sequence_rejected = True

    canonical_selector_bits = 2
    canonical_payload = 0
    fiber: list[tuple[int, str]] = []
    for choice in range(1 << canonical_selector_bits):
        logical_rank = (canonical_payload << canonical_selector_bits) | choice
        physical_rank = (
            int(affine["affine_multiplier"]) * logical_rank
            + int(affine["affine_offset"])
        ) % int(affine["language_modulus"])
        fiber.append((choice, codec.unrank(physical_rank)))
    score = lambda sequence: -float(codec.rank(sequence) / codec.count())
    score_values = [score(sequence) for _choice, sequence in fiber]
    canonical_choice = int(deterministic.argmax_smallest(score_values))
    canonical_sequence = dict(fiber)[canonical_choice]
    noncanonical_choice, noncanonical_sequence = next(
        item for item in fiber if item[0] != canonical_choice
    )
    decoded_noncanonical = decode_payload_candidate(
        codec,
        noncanonical_sequence,
        fixed_bits,
        canonical_selector_bits,
        affine,
    )
    noncanonical_payload_decoder_accepts = decoded_noncanonical[:2] == (
        canonical_payload,
        noncanonical_choice,
    )
    canonical_verifier_accepts_emitted = verify_canonical_candidate(
        codec,
        canonical_sequence,
        fixed_bits,
        canonical_selector_bits,
        affine,
        score,
    )
    canonical_verifier_rejects_noncanonical = not verify_canonical_candidate(
        codec,
        noncanonical_sequence,
        fixed_bits,
        canonical_selector_bits,
        affine,
        score,
    )
    tie_score = lambda _sequence: 0.0
    tie_verifier_selects_smallest_choice = verify_canonical_candidate(
        codec,
        dict(fiber)[0],
        fixed_bits,
        canonical_selector_bits,
        affine,
        tie_score,
    ) and all(
        not verify_canonical_candidate(
            codec,
            sequence,
            fixed_bits,
            canonical_selector_bits,
            affine,
            tie_score,
        )
        for choice, sequence in fiber
        if choice != 0
    )
    try:
        verify_canonical_candidate(
            codec,
            dict(fiber)[0],
            fixed_bits,
            canonical_selector_bits,
            affine,
            lambda _sequence: float("nan"),
        )
        nonfinite_score_rejected = False
    except ValueError:
        nonfinite_score_rejected = True
    rte_half_integer_vectors_pass = bool(
        np.array_equal(
            deterministic.round_to_even_int64(
                [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
            ),
            np.asarray([-4, -2, -2, 0, 0, 2, 2, 4], dtype=np.int64),
        )
    )
    negative_score_vectors_pass = bool(
        np.array_equal(
            deterministic.quantized_score_keys([-2.0e-12, -1.0e-12, 0.0]),
            np.asarray([-2, -1, 0], dtype=np.int64),
        )
    )
    try:
        deterministic.round_to_even_int64([float(1 << 63)])
        positive_int64_overflow_rejected = False
    except OverflowError:
        positive_int64_overflow_rejected = True
    replacement_score = lambda sequence: float(codec.rank(sequence) / codec.count())
    replacement_values = [replacement_score(sequence) for _choice, sequence in fiber]
    replacement_choice = int(deterministic.argmax_smallest(replacement_values))
    replacement_sequence = dict(fiber)[replacement_choice]
    score_replacement_preserves_payload = (
        decode_payload_candidate(
            codec,
            canonical_sequence,
            fixed_bits,
            canonical_selector_bits,
            affine,
        )[0]
        == decode_payload_candidate(
            codec,
            replacement_sequence,
            fixed_bits,
            canonical_selector_bits,
            affine,
        )[0]
        == canonical_payload
    )

    return [
        {
            "language": "exact_108nt_gc49_59_hp3",
            "fixed_payload_bits": fixed_bits,
            "candidate_rows": len(candidates),
            "candidate_exact_roundtrips": len(candidates) if candidate_failures == 0 else 0,
            "candidate_roundtrip_failures": candidate_failures,
            "duplicate_logical_candidates_within_setting": int(
                candidates.duplicated(
                    ["selector_bits", "payload_rank", "choice_index"]
                ).sum()
            ),
            "exact_language_out_of_domain_rank_rejected": exact_out_of_domain_rejected,
            "payload_domain_out_of_domain_sequence_rejected": payload_domain_out_of_domain_sequence_rejected,
            "invalid_constraint_sequence_rejected": invalid_constraint_sequence_rejected,
            "wrong_length_sequence_rejected": wrong_length_sequence_rejected,
            "ambiguous_base_sequence_rejected": ambiguous_base_sequence_rejected,
            "noncanonical_valid_candidate_accepted_by_payload_decoder": noncanonical_payload_decoder_accepts,
            "canonical_emitted_candidate_accepted_by_optional_verifier": canonical_verifier_accepts_emitted,
            "noncanonical_valid_candidate_rejected_by_optional_verifier": canonical_verifier_rejects_noncanonical,
            "exact_tie_selects_smallest_choice_index": tie_verifier_selects_smallest_choice,
            "nonfinite_score_rejected_by_encoder_verifier": nonfinite_score_rejected,
            "rte_half_integer_vectors_pass": rte_half_integer_vectors_pass,
            "negative_score_vectors_pass": negative_score_vectors_pass,
            "positive_int64_overflow_rejected": positive_int64_overflow_rejected,
            "score_replacement_preserves_payload_decoding": score_replacement_preserves_payload,
            "score_decimal_places": deterministic.SCORE_DECIMAL_PLACES,
            "score_rounding_contract": deterministic.ROUNDING_CONTRACT,
            "payload_decoder_contract": "accept any valid declared fiber member; invert rank modulo N; reject inverse q>=2^K",
            "canonical_verifier_contract": "optional quantized-score-and-smallest-choice tie-rule check; not part of payload inversion",
        }
    ]


def summary_text(
    base_rows: list[dict[str, Any]],
    public_rows: list[dict[str, Any]],
    seed_sensitivity_summary: list[dict[str, Any]],
    generated_rows: list[dict[str, Any]],
    rejection_rows: list[dict[str, Any]],
    overlap_rows: list[dict[str, Any]],
) -> str:
    base_row = base_rows[0]
    lines = [
        "# Reversible Choice Codec Validation",
        "",
        "Status: `COMPLETE_REPRODUCIBLE_ASSAY_CALIBRATED_CHOICE_CODEC`",
        "",
        "## Exact formal result",
        "",
        (
            f"The 108-nt GC {base_row['gc_min']}--{base_row['gc_max']} / homopolymer <= "
            f"{base_row['max_homopolymer']} language contains exactly "
            f"`{base_row['exact_total']}` sequences and supports a total "
            f"{base_row['fixed_payload_bits']}-bit domain "
            f"({base_row['fixed_payload_rate_bits_per_nt']:.6f} bits/nt)."
        ),
        (
            f"An independent bottom-up recurrence returned the same total "
            f"(`{base_row['independent_bottom_up_total']}`), and an explicit reduced-length "
            f"enumeration completed {base_row['reduced_all_rank_checks']}/{base_row['reduced_all_rank_checks']} "
            f"all-rank checks with {base_row['reduced_rank_unrank_failures']} failures."
        ),
        "Reserving `r` low-order selector bits gives every payload `2^r` exact candidates; "
        "the selected sequence reranks to the same physical rank and decoding discards only "
        "those `r` bits.  The rate cost is therefore exactly `r` bits per sequence.",
        "",
        "## Held-out public experimental libraries",
        "",
        "Each direction uses a source-only frozen model and an outcome-blind 2,048-sequence "
        "target codebook satisfying the same GC/homopolymer constraints.  One candidate fiber, "
        "not one candidate sequence or hash seed, is the paired analysis unit.",
        "",
    ]
    for direction in ("GCall_to_GCfix", "GCfix_to_GCall"):
        for selector_bits in PRIMARY_SELECTOR_BITS:
            row = next(
                item
                for item in public_rows
                if item["direction"] == direction
                and item["selector_model"] == "P5_combined_context"
                and int(item["selector_bits"]) == selector_bits
            )
            lines.append(
                f"- {direction}, r={selector_bits}: paired experimental gain "
                f"{row['standardized_paired_gain']:.4f} outcome SD; raw gain "
                f"{row['paired_mean_outcome_gain']:.8f}, 95% fiber bootstrap "
                f"[{row['paired_gain_ci_2p5']:.8f}, {row['paired_gain_ci_97p5']:.8f}]; "
                f"{row['roundtrip_payloads']}/{row['roundtrip_payloads']} payloads round-tripped."
            )
    lines.extend(
        [
            "",
            "These are retrospective exact codebooks assembled from previously assayed sequences. "
            "They test the selector on held-out measured outcomes but are not a prospective synthesis "
            "or PCR experiment on newly generated codewords.",
            "",
            "Across 32 prespecified outcome-blind library/hash mappings:",
        ]
    )
    for direction in ("GCall_to_GCfix", "GCfix_to_GCall"):
        for selector_bits in PRIMARY_SELECTOR_BITS:
            p5 = next(
                row
                for row in seed_sensitivity_summary
                if row["direction"] == direction
                and int(row["selector_bits"]) == selector_bits
                and row["selector_model"] == "P5_combined_context"
            )
            increment = next(
                row
                for row in seed_sensitivity_summary
                if row["direction"] == direction
                and int(row["selector_bits"]) == selector_bits
                and row["selector_model"] == "P5_minus_P2_selected_outcome"
            )
            lines.append(
                f"- {direction}, r={selector_bits}: P5 standardized gain range "
                f"[{p5['minimum_standardized_gain']:.4f}, {p5['maximum_standardized_gain']:.4f}], "
                f"positive in {p5['positive_seed_fraction']:.3f} of mappings; P5-minus-P2 "
                f"range [{increment['minimum_standardized_gain']:.4f}, "
                f"{increment['maximum_standardized_gain']:.4f}], positive in "
                f"{increment['positive_seed_fraction']:.3f}."
            )
    lines.extend(
        [
            "These seed ranges are algorithmic sensitivity analyses, not additional biological replicates.",
            "",
            "## Generated-language rate--utility audit",
            "",
        ]
    )
    for source_pool in ("GCall", "GCfix"):
        for selector_bits in PRIMARY_SELECTOR_BITS:
            row = next(
                item
                for item in generated_rows
                if item["selector_source_pool"] == source_pool
                and item["selector_model"] == "P5_combined_context"
                and int(item["selector_bits"]) == selector_bits
            )
            lines.append(
                f"- {source_pool} selector, r={selector_bits}: retained payload "
                f"{row['payload_bits_after_selection']} bits ({row['payload_bits_per_variable_nt']:.6f} bits/nt); "
                f"own-model gain {row['own_model_mean_gain_vs_fiber']:.8f}; opposite-pool-model "
                f"gain {row['cross_model_mean_gain_vs_fiber']:.8f}; selected round-trip failures "
                f"{row['selected_roundtrip_failures']}."
            )
    max_failure = max(float(row["failure_fraction"]) for row in rejection_rows)
    min_active_overlap = min(
        float(row["generated_within_source_minmax_fraction"])
        for row in overlap_rows
        if bool(row["p5_active_feature"])
    )
    invariant_rows = [
        row
        for row in overlap_rows
        if not bool(row["p5_active_feature"])
        and float(row["generated_within_source_minmax_fraction"]) < 0.5
    ]
    lines.extend(
        [
            "",
            f"Threshold rejection had failure fractions up to {max_failure:.6f} at the tested caps, "
            "whereas fixed-choice encoding remained total on every sampled payload.",
            f"The minimum generated-within-source-min/max fraction across active P5 features was {min_active_overlap:.6f}. "
            f"{len(invariant_rows)} low-overlap invariant feature row(s) had a zero P5 coefficient; "
            "feature-level overlap is nevertheless reported because generated-sequence scores remain computational predictions.",
            "",
            "## Statistical and claim contract",
            "",
            "- Public confidence intervals use 2,000 paired bootstrap resamples of candidate fibers.",
            "- Generated-score intervals use payload-fiber resampling and do not represent biological replication.",
            "- Source-only fitting, regularization and thresholds use no target-pool outcomes.",
            "- P0 composition, P2 assay-context and P5 complete-context selectors use identical fibers.",
            "- No result claims mechanism, causal PCR improvement, emitted-sequence wet-lab validation, NUPACK, recovery, full weighted-language capacity or storage-system superiority.",
            "",
            "## Reproduction",
            "",
            "```bash",
            "$PYTHON analysis_tools/validate_paper2_reversible_choice_codec.py --jobs 12 --bootstrap 2000 --generated-payloads 512",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def environment_record(args: argparse.Namespace, affine: dict[str, int]) -> dict[str, Any]:
    packages = {}
    for name in ("joblib", "numpy", "pandas", "scikit-learn", "scipy"):
        packages[name] = importlib.metadata.version(name)
    return {
        "analysis_date": "2026-07-17",
        "base_seed": BASE_SEED,
        "bootstrap_replicates": args.bootstrap,
        "generated_payload_samples_per_selector_setting": args.generated_payloads,
        "public_library_size": PUBLIC_LIBRARY_SIZE,
        "public_selector_bits": PUBLIC_SELECTOR_BITS,
        "generated_selector_bits": GENERATED_SELECTOR_BITS,
        "generated_length": GENERATED_LENGTH,
        "generated_gc_min": GENERATED_GC_MIN,
        "generated_gc_max": GENERATED_GC_MAX,
        "generated_max_homopolymer": GENERATED_MAX_HOMOPOLYMER,
        "score_key_contract": {
            "input_dtype": "float64",
            "output_dtype": "int64",
            "decimal_places": deterministic.SCORE_DECIMAL_PLACES,
            "integer_scale": deterministic.SCORE_SCALE,
            "rounding": deterministic.ROUNDING_CONTRACT,
            "signed_int64_min": deterministic.INT64_MIN,
            "signed_int64_max": deterministic.INT64_MAX,
            "out_of_range_behavior": "reject after RTE before int64 array construction",
            "candidate_order": "lexicographically minimize (-quantized_key, choice_index)",
            "tie_break": "smallest declared choice index",
            "nonfinite_behavior": "reject",
        },
        "parallel_jobs": args.jobs,
        "logical_cpus": os.cpu_count(),
        "platform": platform.platform(),
        "python": sys.version,
        "packages": packages,
        "affine_permutation": affine,
        "command": (
            "$PYTHON analysis_tools/validate_paper2_reversible_choice_codec.py "
            f"--jobs {args.jobs} --bootstrap {args.bootstrap} "
            f"--generated-payloads {args.generated_payloads}"
        ),
        "claim_boundary": CLAIM_BOUNDARY,
    }


def output_manifest() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(OUT_DIR.iterdir()):
        if path.name == "sha256_manifest.tsv" or not path.is_file():
            continue
        rows.append(
            {
                "path": relative(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    rows.append(
        {
            "path": relative(Path(__file__)),
            "bytes": Path(__file__).stat().st_size,
            "sha256": sha256(Path(__file__)),
        }
    )
    return rows


def main() -> int:
    args = parse_args()
    if args.bootstrap != BOOTSTRAP_REPLICATES:
        raise ValueError(f"frozen run requires --bootstrap {BOOTSTRAP_REPLICATES}")
    if args.generated_payloads != 512:
        raise ValueError("frozen run requires --generated-payloads 512")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    provenance = selection.validate_public_sources()
    _dt, pcr, pickle_audit = selection.load_frames(args.jobs)
    models, coefficient_rows = fit_frozen_models(pcr, args.jobs)

    public_rows, comparison_rows, public_samples, public_roundtrips = public_fiber_analysis(
        pcr, models, args.bootstrap
    )
    seed_sensitivity_rows, seed_sensitivity_summary = public_seed_sensitivity(
        pcr, models
    )

    codec = GCHomopolymerCodec(
        GENERATED_LENGTH,
        GENERATED_GC_MIN,
        GENERATED_GC_MAX,
        GENERATED_MAX_HOMOPOLYMER,
    )
    base_rows = base_language_rows(codec)
    if int(base_rows[0]["fixed_payload_bits"]) != 213:
        raise ValueError("unexpected fixed payload width for exact 108-nt base language")
    candidates, affine = generated_candidates(codec, 213, args.generated_payloads)
    generated_roundtrips = generated_roundtrip_summary(
        codec, candidates, 213, affine
    )
    generated_features = generated_feature_frame(candidates, args.jobs)
    generated_rows, generated_samples, rejection_rows = generated_fiber_analysis(
        candidates, generated_features, models, pcr, args.bootstrap
    )
    overlap_rows = feature_overlap_rows(generated_features, pcr, models)
    positional_rows = positional_nucleotide_frequency_rows(candidates, pcr)
    terminal_rows = terminal_window_frequency_rows(candidates, pcr)
    namespace_rows, namespace_summary = generated_namespace_sensitivity(
        codec, 213, models, pcr, args.jobs, args.bootstrap
    )

    public_rows.sort(
        key=lambda row: (row["direction"], int(row["selector_bits"]), row["selector_model"])
    )
    comparison_rows.sort(
        key=lambda row: (row["direction"], int(row["selector_bits"]), row["comparison"])
    )
    coefficient_rows.sort(
        key=lambda row: (row["source_pool"], row["model"], row["feature"])
    )
    generated_rows.sort(
        key=lambda row: (
            row["selector_source_pool"],
            int(row["selector_bits"]),
            row["selector_model"],
        )
    )
    rejection_rows.sort(
        key=lambda row: (
            row["selector_source_pool"],
            int(row["selector_bits"]),
            row["selector_model"],
        )
    )

    write_table(OUT_DIR / "base_language_summary.tsv", base_rows)
    write_table(OUT_DIR / "frozen_model_coefficients.tsv", coefficient_rows)
    write_table(OUT_DIR / "public_fiber_results.tsv", public_rows)
    write_table(OUT_DIR / "public_fiber_model_comparisons.tsv", comparison_rows)
    write_table(OUT_DIR / "public_fiber_selected_samples.tsv", public_samples)
    write_table(OUT_DIR / "public_fiber_roundtrip.tsv", public_roundtrips)
    write_table(OUT_DIR / "public_seed_sensitivity.tsv", seed_sensitivity_rows)
    write_table(
        OUT_DIR / "public_seed_sensitivity_summary.tsv", seed_sensitivity_summary
    )
    write_table(OUT_DIR / "generated_fiber_results.tsv", generated_rows)
    write_table(OUT_DIR / "generated_selected_samples.tsv", generated_samples)
    write_table(OUT_DIR / "generated_roundtrip_summary.tsv", generated_roundtrips)
    write_table(OUT_DIR / "threshold_rejection_baseline.tsv", rejection_rows)
    write_table(OUT_DIR / "generated_feature_overlap.tsv", overlap_rows)
    write_table(
        OUT_DIR / "generated_positional_nucleotide_frequencies.tsv", positional_rows
    )
    write_table(
        OUT_DIR / "generated_terminal_nucleotide_frequencies.tsv", terminal_rows
    )
    write_table(
        OUT_DIR / "generated_mapping_namespace_sensitivity.tsv", namespace_rows
    )
    write_table(
        OUT_DIR / "generated_mapping_namespace_sensitivity_summary.tsv",
        namespace_summary,
    )
    (OUT_DIR / "analysis_summary.md").write_text(
        summary_text(
            base_rows,
            public_rows,
            seed_sensitivity_summary,
            generated_rows,
            rejection_rows,
            overlap_rows,
        ),
        encoding="utf-8",
    )
    (OUT_DIR / "environment_and_seeds.json").write_text(
        json.dumps(environment_record(args, affine), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (OUT_DIR / "source_provenance.json").write_text(
        json.dumps(
            {
                **provenance,
                "pickle_audit": pickle_audit,
                "analysis_script": {
                    "path": relative(Path(__file__)),
                    "sha256": sha256(Path(__file__)),
                },
                "frozen_prediction_source": {
                    "path": relative(FROZEN_PREDICTIONS),
                    "sha256": sha256(FROZEN_PREDICTIONS),
                },
                "claim_boundary": CLAIM_BOUNDARY,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    write_table(OUT_DIR / "sha256_manifest.tsv", output_manifest())
    print(f"Wrote reversible choice-codec outputs to {OUT_DIR}", flush=True)
    print(f"Manifest SHA-256: {sha256(OUT_DIR / 'sha256_manifest.tsv')}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
