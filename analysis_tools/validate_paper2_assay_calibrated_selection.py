#!/usr/bin/env python3
"""Assay-calibrated sequence-selection audit for Paper 2.

This analysis asks a deployment-facing question that is distinct from molecular
mechanism validation: at an identical retained fraction of a public sequence
pool, does a score calibrated only on training folds select sequences with more
favourable experimental outcomes than composition, structural-rule, or frozen
weighted-pair baselines?

The experiment uses repeated nested cross-validation for within-pool estimates,
source-only fitting for cross-pool transfer, and paired sequence-level bootstrap
intervals conditional on the out-of-fold predictions.  It does not show that
Paper 2 codewords were synthesized, amplified, sequenced, or recovered, and it
does not identify a causal molecular mechanism.
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
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import validate_paper2_public_experimental_data as base  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
OUT_DIR = PAPER_DIR / "bioinformatics_reframe" / "assay_calibrated_selection"

BASE_SEED = 20260718
BOOTSTRAP_REPLICATES = 2000
RETENTIONS = (0.10, 0.25, 0.50)
PRIMARY_RETENTION = 0.25

DT_MODELS = OrderedDict(
    {
        "D0_composition": ["gc_deviation_from_0p5", "max_homopolymer"],
        "D1_structural_rules": [
            "gc_deviation_from_0p5",
            "max_homopolymer",
            "candidate_pairs",
            "longest_stem",
        ],
        "D2_assay_calibrated": [
            "gc_deviation_from_0p5",
            "max_homopolymer",
            "candidate_pairs",
            "longest_stem",
            "weighted_pairs",
        ],
    }
)

PCR_MODELS = OrderedDict(
    {
        "P0_composition": base.PCR_COMPOSITION,
        "P1_variable_structure": base.PCR_COMPOSITION
        + ["variable_candidate_pairs", "variable_longest_stem"],
        "P2_assay_context": base.PCR_COMPOSITION
        + base.PCR_ADAPTER_FEATURES
        + ["variable_candidate_pairs", "variable_longest_stem"],
        "P3_assay_calibrated": base.PCR_COMPOSITION
        + base.PCR_ADAPTER_FEATURES
        + [
            "variable_candidate_pairs",
            "variable_longest_stem",
            "variable_weighted_pairs",
        ],
        "P4_full_structural_context": base.PCR_COMPOSITION
        + base.PCR_ADAPTER_FEATURES
        + [
            "variable_candidate_pairs",
            "variable_longest_stem",
            "variable_weighted_pairs",
            "full_amplicon_candidate_pairs",
            "full_amplicon_longest_stem",
        ],
        "P5_combined_context": base.PCR_COMPOSITION
        + base.PCR_ADAPTER_FEATURES
        + [
            "variable_candidate_pairs",
            "variable_longest_stem",
            "variable_weighted_pairs",
            "full_amplicon_candidate_pairs",
            "full_amplicon_longest_stem",
            "full_amplicon_weighted_pairs",
        ],
        "P6_full_amplicon_replacement": base.PCR_COMPOSITION
        + base.PCR_ADAPTER_FEATURES
        + [
            "full_amplicon_candidate_pairs",
            "full_amplicon_longest_stem",
            "full_amplicon_weighted_pairs",
        ],
    }
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(12, (os.cpu_count() or 2) - 2)),
        help="Parallel outer-fold and bootstrap jobs.",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=BOOTSTRAP_REPLICATES,
        help=f"Paired bootstrap replicates (frozen: {BOOTSTRAP_REPLICATES}).",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def validate_public_sources() -> dict[str, Any]:
    """Validate public inputs without re-freezing the subsequently edited manuscript."""

    public_paths = [
        base.DT_CSV,
        base.DT_SUMMARY,
        base.PCR_FILES["GCall"],
        base.PCR_FILES["GCfix"],
        base.PCR_REPO / "LICENSE.txt",
        base.PCR_SUPP_PDF,
        base.PCR_SOURCE_ZIP,
    ]
    hashes: dict[str, str] = {}
    for path in public_paths:
        expected = base.EXPECTED_SHA256[path]
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"SHA-256 mismatch for {path}: {observed} != {expected}")
        hashes[relative(path)] = observed

    commits: dict[str, str] = {}
    for repo, expected in base.EXPECTED_COMMITS.items():
        observed = base.git_value(repo, "rev-parse", "HEAD")
        if observed != expected:
            raise ValueError(f"commit mismatch for {repo}: {observed} != {expected}")
        commits[relative(repo)] = observed

    base.validate_fast_evaluator()
    return {"hashes": hashes, "commits": commits}


def load_frames(jobs: int) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    pcr_raw, pickle_audit = base.load_pcr_data()
    pcr_features = base.get_pcr_features(pcr_raw, jobs=jobs, force=False)
    pcr = pcr_raw.merge(
        pcr_features,
        on=["pool", "sample_id"],
        how="inner",
        validate="one_to_one",
    )
    if len(pcr) != 23_992:
        raise ValueError("PCR merge changed the frozen row count")
    pcr["sample_id"] = pcr["sample_id"].astype(str).str.zfill(6)
    pcr = pcr.sort_values(["pool", "sample_id"], kind="mergesort").reset_index(drop=True)

    dt = pd.read_csv(base.DT_CSV)
    if len(dt) != 24_472:
        raise ValueError("DT4DDS row count mismatch")
    dt["reference_id"] = dt["reference_id"].astype(str)
    dt = dt.sort_values(["pool", "reference_id"], kind="mergesort").reset_index(drop=True)
    return dt, pcr, pickle_audit


def selected_mask(prediction: np.ndarray, retention: float) -> np.ndarray:
    prediction = np.asarray(prediction, dtype=float)
    if np.isnan(prediction).any():
        raise ValueError("selection prediction contains NaN")
    selected_n = max(1, int(math.ceil(retention * len(prediction))))
    order = np.argsort(-prediction, kind="mergesort")
    mask = np.zeros(len(prediction), dtype=bool)
    mask[order[:selected_n]] = True
    return mask


def finite_pool_rate_penalty(retention: float, length_nt: int) -> tuple[float, float]:
    bits_per_sequence = -math.log2(retention)
    return bits_per_sequence, bits_per_sequence / length_nt


def oof_group(
    dataset: str,
    pool: str,
    endpoint: str,
    frame: pd.DataFrame,
    outcome: np.ndarray,
    event: np.ndarray,
    models: OrderedDict[str, list[str]],
    frozen_proxy_column: str,
    id_column: str,
    seed: int,
    jobs: int,
) -> dict[str, Any]:
    metric_rows, coefficient_rows, predictions = base.run_nested_models(
        dataset=dataset,
        analysis="selection_continuous_oof",
        pool=pool,
        endpoint=endpoint,
        frame=frame,
        y=outcome,
        models=models,
        classification=False,
        jobs=jobs,
        seed=seed,
    )
    # The frozen hypothesis treats lower weighted-pair burden as favourable.
    predictions = {"F0_frozen_low_proxy": -frame[frozen_proxy_column].to_numpy(float), **predictions}
    prediction_rows = pd.DataFrame(
        {
            "dataset": dataset,
            "evaluation": "repeated_nested_oof",
            "pool_or_direction": pool,
            "endpoint": endpoint,
            "sequence_id": frame[id_column].astype(str),
            "outcome": outcome,
            "adverse_event": event,
            **{f"prediction__{name}": value for name, value in predictions.items()},
        }
    )
    return {
        "dataset": dataset,
        "evaluation": "repeated_nested_oof",
        "pool": pool,
        "endpoint": endpoint,
        "length_nt": int(frame["length_nt"].iloc[0]) if "length_nt" in frame else 108,
        "outcome": outcome,
        "event": event,
        "predictions": predictions,
        "prediction_rows": prediction_rows,
        "metric_rows": metric_rows,
        "coefficient_rows": coefficient_rows,
    }


def fit_transfer_ridge(
    source: pd.DataFrame,
    target: pd.DataFrame,
    models: OrderedDict[str, list[str]],
    seed: int,
    jobs: int,
) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
    y_source = source["eff"].to_numpy(float)

    def fit_one(index: int, model_name: str, features: list[str]) -> tuple[str, np.ndarray, list[dict[str, Any]]]:
        inner_seed = seed + index
        estimator = Pipeline([("scale", StandardScaler()), ("model", Ridge())])
        inner = KFold(n_splits=5, shuffle=True, random_state=inner_seed)
        search = GridSearchCV(
            estimator,
            {"model__alpha": base.RIDGE_ALPHA_GRID},
            scoring="neg_mean_squared_error",
            cv=inner,
            n_jobs=1,
            refit=True,
        )
        search.fit(source[features].to_numpy(float), y_source)
        prediction = search.predict(target[features].to_numpy(float))
        coefficients = np.asarray(search.best_estimator_.named_steps["model"].coef_).ravel()
        rows = [
            {
                "dataset": "Gimpel2025_PCR",
                "analysis": "selection_continuous_transfer",
                "pool_or_direction": f"{source['pool'].iloc[0]}_to_{target['pool'].iloc[0]}",
                "endpoint": "relative_efficiency",
                "model": model_name,
                "repeat": "transfer",
                "fold": "all_source",
                "feature": feature,
                "standardized_coefficient": float(coefficient),
                "selected_regularization": float(search.best_params_["model__alpha"]),
            }
            for feature, coefficient in zip(features, coefficients)
        ]
        return model_name, prediction, rows

    fitted = Parallel(n_jobs=min(jobs, len(models)), prefer="processes")(
        delayed(fit_one)(index, name, features)
        for index, (name, features) in enumerate(models.items())
    )
    predictions: dict[str, np.ndarray] = {
        "F0_frozen_low_proxy": -target["variable_weighted_pairs"].to_numpy(float)
    }
    coefficients: list[dict[str, Any]] = []
    for name, prediction, rows in fitted:
        predictions[name] = prediction
        coefficients.extend(rows)
    return predictions, coefficients


def transfer_group(
    source_pool: str,
    target_pool: str,
    pcr: pd.DataFrame,
    seed: int,
    jobs: int,
) -> dict[str, Any]:
    source = pcr.loc[pcr["pool"].eq(source_pool)].reset_index(drop=True)
    target = pcr.loc[pcr["pool"].eq(target_pool)].reset_index(drop=True)
    predictions, coefficient_rows = fit_transfer_ridge(
        source, target, PCR_MODELS, seed=seed, jobs=jobs
    )
    outcome = target["eff"].to_numpy(float)
    event = target["label"].to_numpy(int)
    direction = f"{source_pool}_to_{target_pool}"
    prediction_rows = pd.DataFrame(
        {
            "dataset": "Gimpel2025_PCR",
            "evaluation": "source_only_transfer",
            "pool_or_direction": direction,
            "endpoint": "relative_efficiency",
            "sequence_id": target["sample_id"].astype(str),
            "outcome": outcome,
            "adverse_event": event,
            **{f"prediction__{name}": value for name, value in predictions.items()},
        }
    )
    return {
        "dataset": "Gimpel2025_PCR",
        "evaluation": "source_only_transfer",
        "pool": direction,
        "endpoint": "relative_efficiency",
        "length_nt": 108,
        "outcome": outcome,
        "event": event,
        "predictions": predictions,
        "prediction_rows": prediction_rows,
        "metric_rows": [],
        "coefficient_rows": coefficient_rows,
    }


def comparison_pairs(dataset: str) -> list[tuple[str, str, str]]:
    if dataset == "DT4DDS":
        return [
            ("calibration_vs_frozen_proxy", "F0_frozen_low_proxy", "D2_assay_calibrated"),
            ("structural_rules_vs_composition", "D0_composition", "D1_structural_rules"),
            ("weighted_increment", "D1_structural_rules", "D2_assay_calibrated"),
        ]
    return [
        ("calibration_vs_frozen_proxy", "F0_frozen_low_proxy", "P5_combined_context"),
        ("variable_structure_vs_composition", "P0_composition", "P1_variable_structure"),
        ("assay_context_increment", "P1_variable_structure", "P2_assay_context"),
        ("weighted_increment", "P2_assay_context", "P3_assay_calibrated"),
        ("full_structural_context_increment", "P3_assay_calibrated", "P4_full_structural_context"),
        ("full_weighted_increment", "P4_full_structural_context", "P5_combined_context"),
        ("full_amplicon_replacement_sensitivity", "P3_assay_calibrated", "P6_full_amplicon_replacement"),
    ]


def bootstrap_selection_group(
    group: dict[str, Any], n_bootstrap: int, seed: int, batch_size: int = 20
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    outcome = np.asarray(group["outcome"], dtype=float)
    event = np.asarray(group["event"], dtype=float)
    predictions = group["predictions"]
    model_names = list(predictions)
    n = len(outcome)
    outcome_sd = float(np.std(outcome, ddof=1))
    pool_mean = float(np.mean(outcome))
    pool_event_rate = float(np.mean(event))

    masks: dict[tuple[str, float], np.ndarray] = {}
    selection_rows: list[dict[str, Any]] = []
    for retention in RETENTIONS:
        for model in model_names:
            mask = selected_mask(predictions[model], retention)
            achieved_retention = float(mask.sum() / n)
            bits, bits_per_nt = finite_pool_rate_penalty(
                achieved_retention, group["length_nt"]
            )
            masks[(model, retention)] = mask
            selected_mean = float(np.mean(outcome[mask]))
            selected_event = float(np.mean(event[mask]))
            selection_rows.append(
                {
                    "dataset": group["dataset"],
                    "evaluation": group["evaluation"],
                    "pool_or_direction": group["pool"],
                    "endpoint": group["endpoint"],
                    "model": model,
                    "retention_fraction": retention,
                    "achieved_retention_fraction": achieved_retention,
                    "selected_n": int(mask.sum()),
                    "n_sequences": n,
                    "pool_mean_outcome": pool_mean,
                    "selected_mean_outcome": selected_mean,
                    "outcome_gain_vs_pool": selected_mean - pool_mean,
                    "standardized_gain_vs_pool": (selected_mean - pool_mean) / outcome_sd,
                    "outcome_gain_ci_2p5": math.nan,
                    "outcome_gain_ci_97p5": math.nan,
                    "pool_adverse_event_rate": pool_event_rate,
                    "pool_adverse_events": int(event.sum()),
                    "selected_adverse_event_rate": selected_event,
                    "selected_adverse_events": int(event[mask].sum()),
                    "selected_to_pool_event_rate_ratio": (
                        selected_event / pool_event_rate if pool_event_rate > 0 else math.nan
                    ),
                    "event_rate_change_vs_pool": selected_event - pool_event_rate,
                    "event_change_ci_2p5": math.nan,
                    "event_change_ci_97p5": math.nan,
                    "finite_pool_penalty_bits_per_sequence": bits,
                    "finite_pool_penalty_bits_per_nt": bits_per_nt,
                    "rate_boundary": "finite observed-pool retention cost; not full-language capacity",
                }
            )

    pairs = comparison_pairs(group["dataset"])
    comparison_rows: list[dict[str, Any]] = []
    for retention in RETENTIONS:
        for comparison, baseline, extended in pairs:
            if baseline not in predictions or extended not in predictions:
                continue
            baseline_mask = masks[(baseline, retention)]
            extended_mask = masks[(extended, retention)]
            comparison_rows.append(
                {
                    "dataset": group["dataset"],
                    "evaluation": group["evaluation"],
                    "pool_or_direction": group["pool"],
                    "endpoint": group["endpoint"],
                    "comparison": comparison,
                    "baseline_model": baseline,
                    "extended_model": extended,
                    "retention_fraction": retention,
                    "delta_selected_mean": float(
                        np.mean(outcome[extended_mask]) - np.mean(outcome[baseline_mask])
                    ),
                    "delta_selected_mean_ci_2p5": math.nan,
                    "delta_selected_mean_ci_97p5": math.nan,
                    "delta_selected_event_rate": float(
                        np.mean(event[extended_mask]) - np.mean(event[baseline_mask])
                    ),
                    "delta_event_rate_ci_2p5": math.nan,
                    "delta_event_rate_ci_97p5": math.nan,
                    "n_sequences": n,
                    "bootstrap_unit": "sequence; conditional on frozen OOF or transfer predictions",
                }
            )

    all_masks = np.column_stack(
        [masks[(model, retention)] for retention in RETENTIONS for model in model_names]
    ).astype(float)
    mask_keys = [(model, retention) for retention in RETENTIONS for model in model_names]
    mask_index = {key: index for index, key in enumerate(mask_keys)}
    rng = np.random.default_rng(seed)
    probabilities = np.full(n, 1.0 / n)

    selection_boot = {
        key: {"outcome": np.empty(n_bootstrap), "event": np.empty(n_bootstrap)}
        for key in mask_keys
    }
    comparison_boot = {
        (comparison, retention): {
            "outcome": np.empty(n_bootstrap),
            "event": np.empty(n_bootstrap),
        }
        for retention in RETENTIONS
        for comparison, baseline, extended in pairs
        if baseline in predictions and extended in predictions
    }

    for start in range(0, n_bootstrap, batch_size):
        stop = min(start + batch_size, n_bootstrap)
        weights = rng.multinomial(n, probabilities, size=stop - start).astype(float)
        denominator = weights @ all_masks
        selected_outcome = (weights * outcome) @ all_masks / denominator
        selected_event = (weights * event) @ all_masks / denominator
        sampled_pool_outcome = (weights @ outcome) / n
        sampled_pool_event = (weights @ event) / n
        for key, index in mask_index.items():
            selection_boot[key]["outcome"][start:stop] = (
                selected_outcome[:, index] - sampled_pool_outcome
            )
            selection_boot[key]["event"][start:stop] = (
                selected_event[:, index] - sampled_pool_event
            )
        for retention in RETENTIONS:
            for comparison, baseline, extended in pairs:
                if baseline not in predictions or extended not in predictions:
                    continue
                baseline_index = mask_index[(baseline, retention)]
                extended_index = mask_index[(extended, retention)]
                key = (comparison, retention)
                comparison_boot[key]["outcome"][start:stop] = (
                    selected_outcome[:, extended_index] - selected_outcome[:, baseline_index]
                )
                comparison_boot[key]["event"][start:stop] = (
                    selected_event[:, extended_index] - selected_event[:, baseline_index]
                )

    for row in selection_rows:
        key = (row["model"], float(row["retention_fraction"]))
        row["outcome_gain_ci_2p5"], row["outcome_gain_ci_97p5"] = np.quantile(
            selection_boot[key]["outcome"], [0.025, 0.975]
        )
        row["event_change_ci_2p5"], row["event_change_ci_97p5"] = np.quantile(
            selection_boot[key]["event"], [0.025, 0.975]
        )
    for row in comparison_rows:
        key = (row["comparison"], float(row["retention_fraction"]))
        row["delta_selected_mean_ci_2p5"], row["delta_selected_mean_ci_97p5"] = np.quantile(
            comparison_boot[key]["outcome"], [0.025, 0.975]
        )
        row["delta_event_rate_ci_2p5"], row["delta_event_rate_ci_97p5"] = np.quantile(
            comparison_boot[key]["event"], [0.025, 0.975]
        )
    return selection_rows, comparison_rows


def prediction_metrics(groups: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for group in groups:
        outcome = group["outcome"]
        for model, prediction in group["predictions"].items():
            rows.append(
                {
                    "dataset": group["dataset"],
                    "evaluation": group["evaluation"],
                    "pool_or_direction": group["pool"],
                    "endpoint": group["endpoint"],
                    "model": model,
                    "spearman": float(spearmanr(outcome, prediction).statistic),
                    "n_sequences": len(outcome),
                    "metric_boundary": "sequence-level predictive association; not mechanism or recovery",
                }
            )
    return rows


def coefficient_summary_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(rows)
    keys = [
        "dataset",
        "analysis",
        "pool_or_direction",
        "endpoint",
        "model",
        "feature",
    ]
    output: list[dict[str, Any]] = []
    for values, group in frame.groupby(keys, sort=True, dropna=False):
        coefficients = group["standardized_coefficient"].astype(float).to_numpy()
        output.append(
            {
                **dict(zip(keys, values)),
                "median_standardized_coefficient": float(np.median(coefficients)),
                "coefficient_q25": float(np.quantile(coefficients, 0.25)),
                "coefficient_q75": float(np.quantile(coefficients, 0.75)),
                "negative_fits": int(np.sum(coefficients < 0)),
                "positive_fits": int(np.sum(coefficients > 0)),
                "n_fits": len(coefficients),
                "sign_boundary": "association direction within fitted model; not causal effect",
            }
        )
    return output


def summary_text(
    selection_rows: list[dict[str, Any]], comparison_rows: list[dict[str, Any]]
) -> str:
    primary = [
        row
        for row in comparison_rows
        if float(row["retention_fraction"]) == PRIMARY_RETENTION
        and row["comparison"]
        in {
            "calibration_vs_frozen_proxy",
            "weighted_increment",
            "full_structural_context_increment",
            "full_weighted_increment",
        }
    ]
    lines = [
        "# Assay-Calibrated Selection Utility",
        "",
        "Status: `COMPLETE_REPRODUCIBLE_SELECTION_AUDIT`",
        "",
        "## Question",
        "",
        "At an identical retained fraction of each public sequence pool, does training-fold-only calibration select sequences with more favourable measured outcomes than a frozen low-weighted-pair rule or a matched structural baseline?",
        "",
        "## Primary 25% retention comparisons",
        "",
    ]
    for row in primary:
        lines.append(
            f"- {row['dataset']} / {row['evaluation']} / {row['pool_or_direction']} / "
            f"{row['endpoint']} / {row['comparison']}: "
            f"delta selected mean = {float(row['delta_selected_mean']):.6g}, "
            f"95% paired bootstrap [{float(row['delta_selected_mean_ci_2p5']):.6g}, "
            f"{float(row['delta_selected_mean_ci_97p5']):.6g}]; "
            f"delta adverse-event rate = {float(row['delta_selected_event_rate']):.6g}."
        )
    lines.extend(
        [
            "",
            "## Interpretation contract",
            "",
            "- Positive within-pool selection utility supports assay-specific retrospective calibration, not a universal molecular-risk direction.",
            "- A null weighted increment means the weighted-pair scalar adds no utility beyond candidate count, longest stem and declared assay context in that comparison.",
            "- Source-only transfer is reported separately and determines whether a calibration can move between sequence pools without target tuning.",
            "- The finite-pool retention penalty is computed as -log2(selected_n/N) from the achieved retained fraction and is reported per nucleotide only as finite-pool accounting. It is not full-language capacity or achieved storage density.",
            "- Public sequences were not emitted by the Paper 2 codec. No synthesis, PCR, sequencing or recovery experiment on Paper 2 oligos is claimed.",
            "",
            "## Statistical contract",
            "",
            "- The independent analysis unit is one public reference sequence.",
            "- Within-pool predictions are means over five repeated five-fold outer passes; scaling and ridge tuning occur inside training folds.",
            "- Transfer models are fitted and tuned on the complete source pool without target-pool outcomes.",
            "- Intervals use 2,000 paired sequence-level bootstrap resamples conditional on the frozen predictions and fixed selected sets.",
            "- Retention fractions of 10%, 25% and 50% are all retained in the output; 25% is the primary descriptive operating point.",
            "",
        ]
    )
    return "\n".join(lines)


def environment_record(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "analysis_date": "2026-07-17",
        "python": sys.version,
        "platform": platform.platform(),
        "logical_cpus": os.cpu_count(),
        "parallel_jobs": args.jobs,
        "base_seed": BASE_SEED,
        "bootstrap_replicates": args.bootstrap,
        "retention_fractions": RETENTIONS,
        "primary_retention": PRIMARY_RETENTION,
        "outer_repeats": base.OUTER_REPEATS,
        "outer_folds": base.OUTER_FOLDS,
        "inner_folds": base.INNER_FOLDS,
        "ridge_alpha_grid": base.RIDGE_ALPHA_GRID.tolist(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "pandas", "scipy", "scikit-learn", "joblib"]
        },
        "model_features": {"DT4DDS": DT_MODELS, "Gimpel2025_PCR": PCR_MODELS},
        "command": f"$PYTHON {relative(Path(__file__))} --jobs {args.jobs} --bootstrap {args.bootstrap}",
    }


def output_manifest() -> list[dict[str, Any]]:
    paths = sorted(
        path for path in OUT_DIR.iterdir() if path.is_file() and path.name != "sha256_manifest.tsv"
    )
    return [
        {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in paths
    ]


def main() -> int:
    args = parse_args()
    if args.bootstrap != BOOTSTRAP_REPLICATES:
        raise ValueError(f"frozen analysis requires --bootstrap {BOOTSTRAP_REPLICATES}")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    provenance = validate_public_sources()
    dt, pcr, pickle_audit = load_frames(args.jobs)

    groups: list[dict[str, Any]] = []
    dt_specs = [
        ("mean_log2_cpm_plus1", "missing_any"),
        ("day7_vs_day0_log2fc_plus1", "missing_day7"),
    ]
    dt_pools = ["Genscript_GCall", "Twist_GCall"]
    for pool_index, pool in enumerate(dt_pools):
        frame = dt.loc[dt["pool"].eq(pool)].reset_index(drop=True)
        for endpoint_index, (endpoint, event_column) in enumerate(dt_specs):
            groups.append(
                oof_group(
                    dataset="DT4DDS",
                    pool=pool,
                    endpoint=endpoint,
                    frame=frame,
                    outcome=frame[endpoint].to_numpy(float),
                    event=frame[event_column].astype(int).to_numpy(),
                    models=DT_MODELS,
                    frozen_proxy_column="weighted_pairs",
                    id_column="reference_id",
                    seed=BASE_SEED + pool_index * 1000 + endpoint_index * 100,
                    jobs=args.jobs,
                )
            )

    for pool_index, pool in enumerate(["GCall", "GCfix"]):
        frame = pcr.loc[pcr["pool"].eq(pool)].reset_index(drop=True)
        groups.append(
            oof_group(
                dataset="Gimpel2025_PCR",
                pool=pool,
                endpoint="relative_efficiency",
                frame=frame,
                outcome=frame["eff"].to_numpy(float),
                event=frame["label"].astype(int).to_numpy(),
                models=PCR_MODELS,
                frozen_proxy_column="variable_weighted_pairs",
                id_column="sample_id",
                seed=BASE_SEED + 10_000 + pool_index * 1000,
                jobs=args.jobs,
            )
        )

    for direction_index, (source_pool, target_pool) in enumerate(
        [("GCall", "GCfix"), ("GCfix", "GCall")]
    ):
        groups.append(
            transfer_group(
                source_pool,
                target_pool,
                pcr,
                seed=BASE_SEED + 20_000 + direction_index * 1000,
                jobs=args.jobs,
            )
        )

    bootstrap_results = Parallel(n_jobs=min(args.jobs, len(groups)), prefer="processes")(
        delayed(bootstrap_selection_group)(
            group, args.bootstrap, BASE_SEED + 30_000 + index * 1000
        )
        for index, group in enumerate(groups)
    )
    selection_rows = [row for result in bootstrap_results for row in result[0]]
    comparison_rows = [row for result in bootstrap_results for row in result[1]]
    metric_rows = prediction_metrics(groups)
    coefficient_rows = [row for group in groups for row in group["coefficient_rows"]]
    coefficient_summary = coefficient_summary_rows(coefficient_rows)
    prediction_frame = pd.concat([group["prediction_rows"] for group in groups], ignore_index=True)

    selection_rows.sort(
        key=lambda row: (
            row["dataset"],
            row["evaluation"],
            row["pool_or_direction"],
            row["endpoint"],
            float(row["retention_fraction"]),
            row["model"],
        )
    )
    comparison_rows.sort(
        key=lambda row: (
            row["dataset"],
            row["evaluation"],
            row["pool_or_direction"],
            row["endpoint"],
            float(row["retention_fraction"]),
            row["comparison"],
        )
    )
    metric_rows.sort(
        key=lambda row: (
            row["dataset"], row["evaluation"], row["pool_or_direction"], row["endpoint"], row["model"]
        )
    )
    coefficient_rows.sort(
        key=lambda row: (
            row["dataset"],
            row["analysis"],
            row["pool_or_direction"],
            row["endpoint"],
            row["model"],
            str(row["repeat"]),
            str(row["fold"]),
            row["feature"],
        )
    )

    write_table(OUT_DIR / "selection_metrics.tsv", selection_rows)
    write_table(OUT_DIR / "selection_model_comparisons.tsv", comparison_rows)
    write_table(OUT_DIR / "prediction_metrics.tsv", metric_rows)
    write_table(OUT_DIR / "calibration_coefficients.tsv", coefficient_rows)
    write_table(OUT_DIR / "calibration_coefficient_summary.tsv", coefficient_summary)
    prediction_frame.to_csv(
        OUT_DIR / "sequence_level_predictions.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    (OUT_DIR / "analysis_summary.md").write_text(
        summary_text(selection_rows, comparison_rows), encoding="utf-8"
    )
    (OUT_DIR / "environment_and_seeds.json").write_text(
        json.dumps(environment_record(args), indent=2, sort_keys=True) + "\n",
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
                "claim_boundary": "public retrospective sequence selection; not Paper 2 oligo wet-lab validation, mechanism, recovery, storage density, or full-language capacity",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    write_table(OUT_DIR / "sha256_manifest.tsv", output_manifest())
    print(f"Wrote assay-calibrated selection outputs to {OUT_DIR}", flush=True)
    print(f"Manifest SHA-256: {sha256(OUT_DIR / 'sha256_manifest.tsv')}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
