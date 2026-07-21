#!/usr/bin/env python3
"""SOTA, external-laboratory and two-stage uncertainty audit for Paper 2.

The program adds three evidence layers without changing the exact codec:

1. fixed-commit Gimpel et al. 1D-CNN predictions and the publication's
   GC/homopolymer/mfold rule are evaluated on the identical public fibers;
2. frozen GCall/GCfix selectors are tested on a separately locked external-
   laboratory pool under the same assay conditions; and
3. source-model and target-fiber uncertainty are propagated together for P2
   and P5, with fixed-prediction intervals and within-fiber randomization tests
   retained as explicitly different estimands.

Public measurements remain retrospective with respect to Paper 2 codebooks.
Generated Paper 2 codec outputs are not synthesized or measured here.  The
analysis provides no NUPACK, material, recovery, full 110-nt weighted-language
capacity, causal-mechanism or end-to-end storage-superiority evidence.
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
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import norm
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV, KFold, PredefinedSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import validate_paper2_assay_calibrated_selection as selection  # noqa: E402
import validate_paper2_public_experimental_data as base  # noqa: E402
import validate_paper2_reversible_choice_codec as codec  # noqa: E402
import paper2_deterministic_selection as deterministic  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
INPUT_DIR = ROOT / "external_data" / "paper2_gimpel2025_public_release_af62c57"
OUT_DIR = PAPER_DIR / "bioinformatics_reframe" / "sota_and_external_validation"

UPSTREAM_COMMIT = "af62c57f9a90ecdfdd0f1623441e82bdb7e082c1"
BASE_SEED = 20260718
BOOTSTRAP_REPLICATES = 2000
TWO_STAGE_REPLICATES = 2000
RANDOMIZATION_REPLICATES = 100_000
PUBLIC_BITS = (0, 2, 4)
INFERENTIAL_BITS = (2, 4)
MODEL_NAMES = codec.MODEL_NAMES
TWO_STAGE_MODELS = ("P2_assay_context", "P5_combined_context")

CLAIM_BOUNDARY = (
    "retrospective measured public-sequence selection and an exact computational "
    "choice-codec interface; not prospective codec-output wet lab, mechanism, "
    "NUPACK, material, sequencing, recovery, full 110-nt weighted-language "
    "capacity, achieved storage density, end-to-end superiority, or unpublished-"
    "sister-manuscript evidence"
)

INPUT_HASHES = {
    "External_GCall2GCfix_2perc_regression_plus_probs.csv":
        "b694f853322114eb7bf0dbd249c027520e91d5b1d9db93c5b198b011892ac9de",
    "External_GCfix2GCall_2perc_regression_plus_probs.csv":
        "7bfc2440f3963505a8885cd4cd18b44fe0c99fffa854956793f083f5f6acbcff",
    "external-validation-predictions.csv":
        "2195ede4ce1a7172840c3091553e73a394c777a1d61dbad292b0ba556b984fb2",
    "sequence_data_anonymized_no_duplicates.csv":
        "53e9e4bdd4b740061b1650d3b4c8adcb3012312af20bbfbe5b75748865017864",
    "GCall_seqprops.csv":
        "154ea0a1cf48b368a03037f1681fccf434bb3d3bd66ff4f3197eef55fceff864",
    "GCfix_seqprops.csv":
        "6463d9d91a96cf33e67362df5e9f07985df61e08f6b49a09f2992da0fa9d48dd",
    "GCall_params.csv":
        "376eb0db5c0c391cc9c5941795516ca26c0dbe6ff289b1e8e4b0981828db7918",
    "GCfix_params.csv":
        "bd8b345ca06d21af630656b725951c93bd166afb1f23f67b8f04f36334cd3565",
    "45_external_validation_filtering_run_analysis.ipynb":
        "455fd443ba799224d6cf9185644ce2fbfdc395a4c6cc20aa9c87e296cdeefa16",
}

UPSTREAM_PATHS = {
    "External_GCall2GCfix_2perc_regression_plus_probs.csv":
        "CNN/results_revision/External_GCall2GCfix_2perc_regression_plus_probs.csv",
    "External_GCfix2GCall_2perc_regression_plus_probs.csv":
        "CNN/results_revision/External_GCfix2GCall_2perc_regression_plus_probs.csv",
    "external-validation-predictions.csv":
        "analysis/data/machine_learning_results/external-validation-predictions.csv",
    "sequence_data_anonymized_no_duplicates.csv":
        "analysis/41_external_validation_generation/sequence_data_anonymized_no_duplicates.csv",
    "GCall_seqprops.csv": "analysis/data/internal_datasets/GCall/seqprops.csv",
    "GCfix_seqprops.csv": "analysis/data/internal_datasets/GCfix/seqprops.csv",
    "GCall_params.csv": "analysis/data/internal_datasets/GCall/params.csv",
    "GCfix_params.csv": "analysis/data/internal_datasets/GCfix/params.csv",
    "45_external_validation_filtering_run_analysis.ipynb":
        "analysis/45_external_validation_filtering/run_analysis.ipynb",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(12, (os.cpu_count() or 2) - 2)),
        help="Parallel feature and two-stage-bootstrap jobs.",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=BOOTSTRAP_REPLICATES,
        help="Target-fiber bootstrap replicates.",
    )
    parser.add_argument(
        "--two-stage",
        type=int,
        default=TWO_STAGE_REPLICATES,
        help="Source-plus-target two-stage bootstrap replicates.",
    )
    parser.add_argument(
        "--randomizations",
        type=int,
        default=RANDOMIZATION_REPLICATES,
        help="Within-fiber random-choice replicates.",
    )
    parser.add_argument(
        "--force-features",
        action="store_true",
        help="Recompute the locked external-library feature cache.",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Use small replicate counts and write below a smoke subdirectory.",
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
    fieldnames: list[str] = []
    for row in rows:
        for field in row:
            if field not in fieldnames:
                fieldnames.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_inputs() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, expected in INPUT_HASHES.items():
        path = INPUT_DIR / name
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"SHA-256 mismatch for {path}: {observed} != {expected}")
        rows.append(
            {
                "local_path": relative(path),
                "upstream_repository": "https://github.com/BorgwardtLab/PCR-bias",
                "upstream_commit": UPSTREAM_COMMIT,
                "upstream_path": UPSTREAM_PATHS[name],
                "bytes": path.stat().st_size,
                "sha256": observed,
            }
        )
    readme = INPUT_DIR / "README.md"
    if readme.is_file():
        rows.append(
            {
                "local_path": relative(readme),
                "upstream_repository": "local acquisition note",
                "upstream_commit": UPSTREAM_COMMIT,
                "upstream_path": "not_applicable",
                "bytes": readme.stat().st_size,
                "sha256": sha256(readme),
            }
        )
    write_table(INPUT_DIR / "sha256_manifest.tsv", rows)
    return rows


def zfill_id(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(6)


def load_cross_pool_prediction(
    path: Path, target: pd.DataFrame, direction: str
) -> tuple[pd.DataFrame, dict[str, Any]]:
    frame = pd.read_csv(path)
    required = {
        "seq_id",
        "true_efficiency",
        "binary_label",
        "pred_efficiency",
        "pred_probability",
    }
    if set(frame.columns) != required:
        raise ValueError(f"{path.name}: unexpected schema {frame.columns.tolist()}")
    missing = frame["seq_id"].isna()
    if int(missing.sum()) != 1:
        raise ValueError(f"{path.name}: expected exactly one blank seq_id")
    frame.loc[missing, "seq_id"] = 0
    frame["sample_id"] = frame["seq_id"].astype(int).map(lambda value: f"{value:06d}")
    if frame["sample_id"].duplicated().any() or len(frame) != len(target):
        raise ValueError(f"{path.name}: target identifier cardinality mismatch")
    joined = target[["sample_id", "eff", "label"]].merge(
        frame, on="sample_id", how="left", validate="one_to_one"
    )
    if joined["pred_efficiency"].isna().any():
        raise ValueError(f"{path.name}: missing target prediction after join")
    eff_delta = float(np.max(np.abs(joined["eff"] - joined["true_efficiency"])))
    label_mismatch = int((joined["label"].astype(int) != joined["binary_label"].astype(int)).sum())
    if eff_delta > 2e-7 or label_mismatch:
        raise ValueError(
            f"{path.name}: outcome audit failed delta={eff_delta}, labels={label_mismatch}"
        )
    audit = {
        "artifact": path.name,
        "direction": direction,
        "rows": len(frame),
        "unique_ids": int(frame["sample_id"].nunique()),
        "blank_seq_id_repaired_as": "000000",
        "true_efficiency_max_abs_difference": eff_delta,
        "binary_label_mismatches": label_mismatch,
        "prediction_efficiency_min": float(frame["pred_efficiency"].min()),
        "prediction_efficiency_max": float(frame["pred_efficiency"].max()),
    }
    return frame.set_index("sample_id"), audit


def load_sequence_properties(pool: str, target: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    props = pd.read_csv(INPUT_DIR / f"{pool}_seqprops.csv", dtype={"id": str})
    props["sample_id"] = props["id"].str.zfill(6)
    if len(props) != 12_000 or props["sample_id"].duplicated().any():
        raise ValueError(f"{pool}: sequence-properties cardinality mismatch")
    params = pd.read_csv(INPUT_DIR / f"{pool}_params.csv", dtype={"seq_id": str})
    params["sample_id"] = params["seq_id"].str.zfill(6)
    joined = target[["sample_id", "sequence", "eff"]].merge(
        props[["sample_id", "GC", "hp", "dg"]],
        on="sample_id",
        how="left",
        validate="one_to_one",
    ).merge(
        params[["sample_id", "eff"]].rename(columns={"eff": "params_eff"}),
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    if joined[["GC", "hp", "dg", "params_eff"]].isna().any().any():
        raise ValueError(f"{pool}: missing published sequence property or parameter")
    computed_gc = joined["sequence"].map(lambda s: (s.count("G") + s.count("C")) / len(s))
    computed_hp = joined["sequence"].map(codec.max_homopolymer)
    gc_delta = float(np.max(np.abs(computed_gc - joined["GC"])))
    hp_mismatch = int((computed_hp.astype(int) != joined["hp"].astype(int)).sum())
    eff_delta = float(np.max(np.abs(joined["eff"] - joined["params_eff"])))
    if gc_delta > 1e-14 or hp_mismatch or eff_delta > 2e-14:
        raise ValueError(
            f"{pool}: property audit failed gc={gc_delta}, hp={hp_mismatch}, eff={eff_delta}"
        )
    audit = {
        "artifact": f"{pool}_seqprops_and_params",
        "pool": pool,
        "target_rows": len(target),
        "published_property_rows": len(props),
        "computed_gc_max_abs_difference": gc_delta,
        "homopolymer_mismatches": hp_mismatch,
        "parameter_efficiency_max_abs_difference": eff_delta,
        "published_hard_filter": "GC>0.4; GC<0.6; hp<5; dg>-15 kcal/mol",
    }
    return joined.set_index("sample_id")[["GC", "hp", "dg"]], audit


def load_external_pool() -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    predictions = pd.read_csv(INPUT_DIR / "external-validation-predictions.csv")
    sequences = pd.read_csv(INPUT_DIR / "sequence_data_anonymized_no_duplicates.csv")
    if len(predictions) != 11_904 or len(sequences) != 11_904:
        raise ValueError("external pool row count mismatch")
    if predictions["seq_id_anonymized"].duplicated().any() or sequences[
        "seq_id_anonymized"
    ].duplicated().any():
        raise ValueError("external pool duplicate anonymized identifier")
    merged = predictions.merge(
        sequences[["seq_id_anonymized", "seq"]],
        on="seq_id_anonymized",
        how="inner",
        validate="one_to_one",
    )
    sequence_mismatches = int((merged["sequence"] != merged["seq"]).sum())
    if sequence_mismatches:
        raise ValueError("external sequence tables disagree")
    if int((~merged["has_insertedmotif"]).sum()) != 10_000:
        raise ValueError("external random-sequence subset is not 10,000")
    random = merged.loc[~merged["has_insertedmotif"]].copy()
    random = random.loc[random["eff_Taq"].notna()].copy()
    if len(random) != 9_995:
        raise ValueError(f"expected 9,995 random sequences with Taq outcome, observed {len(random)}")
    if not random["sequence"].map(lambda value: len(value) == 108 and set(value) <= set("ACGT")).all():
        raise ValueError("external random pool contains invalid sequence")
    x, counts = np.unique(random["eff_Taq"].to_numpy(float), return_counts=True)
    cumulative = np.cumsum(counts) / counts.sum()
    threshold = float(x[int(np.argmax(cumulative > 0.02))])
    random["external_low_efficiency"] = (random["eff_Taq"] < threshold).astype(int)
    random["sequence_sha256"] = random["sequence"].map(
        lambda value: hashlib.sha256(value.encode("ascii")).hexdigest()
    )
    if random["sequence_sha256"].duplicated().any():
        raise ValueError("duplicate external random-pool sequence")
    audits = [
        {
            "artifact": "external_pool_join",
            "rows_after_deduplication": len(merged),
            "unique_ids": int(merged["seq_id_anonymized"].nunique()),
            "sequence_mismatches": sequence_mismatches,
            "random_rows": int((~merged["has_insertedmotif"]).sum()),
            "motif_inserted_rows": int(merged["has_insertedmotif"].sum()),
            "random_Taq_nonmissing_rows": len(random),
            "random_Q5_nonmissing_rows": int(
                merged.loc[~merged["has_insertedmotif"], "eff_Q5"].notna().sum()
            ),
        },
        {
            "artifact": "external_Taq_event_definition",
            "eligible_random_rows": len(random),
            "bottom_2pct_threshold": threshold,
            "strictly_below_threshold_events": int(random["external_low_efficiency"].sum()),
            "event_fraction": float(random["external_low_efficiency"].mean()),
        },
    ]
    return random.reset_index(drop=True), audits


def load_external_q5_pool() -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Load the secondary exploratory workflow-shift endpoint without Taq outcomes."""

    predictions = pd.read_csv(INPUT_DIR / "external-validation-predictions.csv")
    sequences = pd.read_csv(INPUT_DIR / "sequence_data_anonymized_no_duplicates.csv")
    merged = predictions.merge(
        sequences[["seq_id_anonymized", "seq"]],
        on="seq_id_anonymized",
        how="inner",
        validate="one_to_one",
    )
    if int((merged["sequence"] != merged["seq"]).sum()):
        raise ValueError("external Q5 sequence tables disagree")
    random = merged.loc[(~merged["has_insertedmotif"]) & merged["eff_Q5"].notna()].copy()
    if len(random) != 9_955:
        raise ValueError(f"expected 9,955 random sequences with Q5 outcome, observed {len(random)}")
    x, counts = np.unique(random["eff_Q5"].to_numpy(float), return_counts=True)
    cumulative = np.cumsum(counts) / counts.sum()
    threshold = float(x[int(np.argmax(cumulative > 0.02))])
    random["external_low_efficiency_Q5"] = (random["eff_Q5"] < threshold).astype(int)
    random["sequence_sha256"] = random["sequence"].map(
        lambda value: hashlib.sha256(value.encode("ascii")).hexdigest()
    )
    eligible = int(random["sequence"].map(codec.exact_gc_hp_eligible).sum())
    if eligible != 2_040:
        raise ValueError(f"expected 2,040 exact-codec-eligible Q5 rows, observed {eligible}")
    return random.reset_index(drop=True), [
        {
            "artifact": "external_Q5_sensitivity_pool",
            "random_Q5_nonmissing_rows": len(random),
            "exact_codec_eligible_rows": eligible,
            "bottom_2pct_threshold": threshold,
            "strictly_below_threshold_events": int(
                random["external_low_efficiency_Q5"].sum()
            ),
            "event_fraction": float(random["external_low_efficiency_Q5"].mean()),
            "endpoint_role": "workflow-shift sensitivity; not the primary external endpoint",
        }
    ]


def locked_external_library(
    random: pd.DataFrame,
    *,
    endpoint: str = "external_Taq",
    library_size: int = codec.PUBLIC_LIBRARY_SIZE,
) -> pd.DataFrame:
    eligible = random.loc[random["sequence"].map(codec.exact_gc_hp_eligible)].copy()
    if len(eligible) < library_size:
        raise ValueError("insufficient exact-GC/HP external sequences")
    if library_size & (library_size - 1):
        raise ValueError("external codebook size must be a power of two")
    namespace = (
        "paper2-choice-external-v1|external_Taq"
        if endpoint == "external_Taq"
        else f"paper2-choice-external-q5-sensitivity-v1|{endpoint}"
    )
    eligible["library_key"] = eligible["sequence_sha256"].map(
        lambda value: hashlib.sha256(
            f"{namespace}|{value}".encode("ascii")
        ).hexdigest()
    )
    library = (
        eligible.sort_values(["library_key", "sequence_sha256"], kind="mergesort")
        .iloc[:library_size]
        .reset_index(drop=True)
    )
    library["base_rank"] = np.arange(len(library), dtype=int)
    if library["sequence"].duplicated().any():
        raise ValueError("locked external library contains duplicate sequence")
    return library


def external_features(
    library: pd.DataFrame,
    output_dir: Path,
    jobs: int,
    force: bool,
    *,
    cache_name: str = "external_locked_library_features.tsv",
    pool_label: str = "External_Taq",
) -> pd.DataFrame:
    path = output_dir / cache_name
    expected = library[["seq_id_anonymized", "sequence_sha256"]].copy()
    expected = expected.rename(columns={"seq_id_anonymized": "sample_id"})
    if path.is_file() and not force:
        cached = pd.read_csv(path, sep="\t", dtype={"sample_id": str})
        joined = expected.merge(
            cached[["sample_id", "sequence_sha256"]],
            on="sample_id",
            suffixes=("_expected", "_cached"),
            validate="one_to_one",
        )
        if len(cached) != len(library) or not joined[
            "sequence_sha256_expected"
        ].eq(joined["sequence_sha256_cached"]).all():
            raise ValueError("external feature cache does not match locked library")
        return cached
    tasks = [
        (pool_label, str(row.seq_id_anonymized), str(row.sequence))
        for row in library[["seq_id_anonymized", "sequence"]].itertuples(index=False)
    ]
    rows = Parallel(n_jobs=jobs, prefer="processes", batch_size=32)(
        delayed(base.pcr_feature_row)(pool, sample_id, sequence)
        for pool, sample_id, sequence in tasks
    )
    features = pd.DataFrame(rows).rename(columns={"sample_id": "sample_id"})
    features.to_csv(path, sep="\t", index=False, float_format="%.12g")
    return features


def bootstrap_distribution(values: np.ndarray, replicates: int, seed: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("bootstrap requires at least two finite values")
    rng = np.random.default_rng(seed)
    result = np.empty(replicates, dtype=float)
    batch = 128
    for start in range(0, replicates, batch):
        stop = min(replicates, start + batch)
        indices = rng.integers(0, len(values), size=(stop - start, len(values)))
        result[start:stop] = values[indices].mean(axis=1)
    return result


def percentile_interval(distribution: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(np.asarray(distribution, dtype=float), [0.025, 0.975])
    return float(low), float(high)


def bca_interval(values: np.ndarray, distribution: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    distribution = np.asarray(distribution, dtype=float)
    observed = float(values.mean())
    less = (np.sum(distribution < observed) + 0.5 * np.sum(distribution == observed)) / len(
        distribution
    )
    less = float(np.clip(less, 1 / (2 * len(distribution)), 1 - 1 / (2 * len(distribution))))
    z0 = float(norm.ppf(less))
    total = float(values.sum())
    jackknife = (total - values) / (len(values) - 1)
    center = float(jackknife.mean())
    delta = center - jackknife
    denominator = 6.0 * float(np.sum(delta**2) ** 1.5)
    acceleration = 0.0 if denominator == 0 else float(np.sum(delta**3) / denominator)
    adjusted: list[float] = []
    for alpha in (0.025, 0.975):
        z_alpha = float(norm.ppf(alpha))
        value = norm.cdf(z0 + (z0 + z_alpha) / (1 - acceleration * (z0 + z_alpha)))
        adjusted.append(float(np.clip(value, 0, 1)))
    low, high = np.quantile(distribution, adjusted)
    return float(low), float(high)


def random_choice_test(
    outcomes: np.ndarray, observed_gain: float, replicates: int, seed: int
) -> dict[str, float]:
    outcomes = np.asarray(outcomes, dtype=float)
    fibers, width = outcomes.shape
    baseline = float(outcomes.mean(axis=1).mean())
    rng = np.random.default_rng(seed)
    estimates = np.empty(replicates, dtype=float)
    batch = 256
    fiber_index = np.arange(fibers)[None, :]
    for start in range(0, replicates, batch):
        stop = min(replicates, start + batch)
        choices = rng.integers(0, width, size=(stop - start, fibers))
        selected = outcomes[fiber_index, choices]
        estimates[start:stop] = selected.mean(axis=1) - baseline
    exceedances_one_sided = int(np.sum(estimates >= observed_gain))
    exceedances_two_sided = int(np.sum(np.abs(estimates) >= abs(observed_gain)))
    one_sided = (1 + exceedances_one_sided) / (replicates + 1)
    two_sided = (1 + exceedances_two_sided) / (replicates + 1)
    return {
        "randomization_null_mean": float(estimates.mean()),
        "randomization_null_sd": float(estimates.std(ddof=1)),
        "randomization_p_one_sided_positive": float(one_sided),
        "randomization_p_two_sided": float(two_sided),
        "one_sided_exceedances_b": exceedances_one_sided,
        "two_sided_exceedances_b": exceedances_two_sided,
        "monte_carlo_p_formula": "(b+1)/(B+1)",
        "minimum_attainable_monte_carlo_p": float(1 / (replicates + 1)),
    }


def summarize_gains(
    gains: np.ndarray,
    outcome_sd: float,
    bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    gains = np.asarray(gains, dtype=float)
    distribution = bootstrap_distribution(gains, bootstrap, seed)
    low, high = percentile_interval(distribution)
    bca_low, bca_high = bca_interval(gains, distribution)
    q25, median, q75 = np.quantile(gains, [0.25, 0.5, 0.75])
    mean = float(gains.mean())
    return {
        "paired_mean_outcome_gain": mean,
        "paired_gain_ci_2p5": low,
        "paired_gain_ci_97p5": high,
        "paired_gain_bca_ci_2p5": bca_low,
        "paired_gain_bca_ci_97p5": bca_high,
        "standardized_paired_gain": mean / outcome_sd,
        "fiber_gain_q25": float(q25),
        "fiber_gain_median": float(median),
        "fiber_gain_q75": float(q75),
        "negative_gain_fraction": float((gains < 0).mean()),
        "zero_gain_fraction": float((gains == 0).mean()),
        "positive_gain_fraction": float((gains > 0).mean()),
    }


def evaluate_library(
    dataset: str,
    direction: str,
    library: pd.DataFrame,
    outcome_column: str,
    event_column: str,
    selectors: OrderedDict[str, dict[str, Any]],
    bootstrap: int,
    randomizations: int,
    seed_offset: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    summary_rows: list[dict[str, Any]] = []
    distribution_rows: list[dict[str, Any]] = []
    randomization_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    outcome_all = library[outcome_column].to_numpy(float)
    event_all = library[event_column].to_numpy(int)
    outcome_sd = float(outcome_all.std(ddof=0))
    if outcome_sd <= 0:
        raise ValueError(f"{dataset}/{direction}: zero outcome SD")
    base_payload_bits = int(math.log2(len(library)))
    if (1 << base_payload_bits) != len(library):
        raise ValueError(f"{dataset}/{direction}: library size is not a power of two")

    for bits_index, bits in enumerate(PUBLIC_BITS):
        width = 1 << bits
        fibers = len(library) // width
        outcomes = outcome_all.reshape(fibers, width)
        events = event_all.reshape(fibers, width)
        fiber_mean = outcomes.mean(axis=1)
        fiber_event_mean = events.mean(axis=1)
        selected_by_model: dict[str, np.ndarray] = {}
        selected_event_by_model: dict[str, np.ndarray] = {}

        for model_index, (name, contract) in enumerate(selectors.items()):
            mode = contract["mode"]
            accepted = np.ones(fibers, dtype=bool)
            if mode in {"max", "min"}:
                scores = np.asarray(contract["values"], dtype=float).reshape(fibers, width)
                choices = (
                    deterministic.argmax_smallest(scores, axis=1)
                    if mode == "max"
                    else deterministic.argmin_smallest(scores, axis=1)
                )
                score_stability = deterministic.choice_stability_audit(scores, mode)
                fiber_score_variance = scores.var(axis=1, ddof=0)
                fiber_score_variance_mean: float | str = float(fiber_score_variance.mean())
                selected_score = scores[np.arange(fibers), choices]
                if width > 1:
                    ordered_scores = np.sort(scores, axis=1)
                    if mode == "max":
                        runner_up_score = ordered_scores[:, -2]
                        score_margin = selected_score - runner_up_score
                    else:
                        runner_up_score = ordered_scores[:, 1]
                        score_margin = runner_up_score - selected_score
                else:
                    runner_up_score = np.full(fibers, np.nan)
                    score_margin = np.full(fibers, np.nan)
            elif mode == "first_pass":
                passes = np.asarray(contract["values"], dtype=bool).reshape(fibers, width)
                accepted = passes.any(axis=1)
                choices = passes.argmax(axis=1)
                score_stability = {
                    "fibers": fibers,
                    "choices_changed_by_quantization": "",
                    "quantized_extremum_tie_fibers": "",
                    "minimum_raw_top_gap": "",
                    "median_raw_top_gap": "",
                }
                fiber_score_variance_mean = ""
                fiber_score_variance = np.full(fibers, np.nan)
                selected_score = np.full(fibers, np.nan)
                runner_up_score = np.full(fibers, np.nan)
                score_margin = np.full(fibers, np.nan)
            else:
                raise ValueError(f"unsupported selector mode {mode}")
            selected_outcome = outcomes[np.arange(fibers), choices]
            selected_event = events[np.arange(fibers), choices]
            gains = selected_outcome - fiber_mean
            event_changes = selected_event - fiber_event_mean
            accepted_gains = gains[accepted]
            accepted_events = selected_event[accepted]
            accepted_event_changes = event_changes[accepted]
            if len(accepted_gains) < 2:
                raise ValueError(f"{dataset}/{direction}/{name}/r={bits}: too few emitted fibers")
            stats = summarize_gains(
                accepted_gains,
                outcome_sd,
                bootstrap,
                BASE_SEED + seed_offset + bits_index * 10_000 + model_index * 100,
            )
            selected_by_model[name] = np.where(accepted, selected_outcome, np.nan)
            selected_event_by_model[name] = np.where(accepted, selected_event, np.nan)
            event_distribution = bootstrap_distribution(
                accepted_event_changes,
                bootstrap,
                BASE_SEED
                + 3_000_000
                + seed_offset
                + bits_index * 10_000
                + model_index * 100,
            )
            event_low, event_high = percentile_interval(event_distribution)
            event_bca_low, event_bca_high = bca_interval(
                accepted_event_changes, event_distribution
            )
            summary_rows.append(
                {
                    "dataset": dataset,
                    "direction": direction,
                    "selector_model": name,
                    "selector_class": contract["selector_class"],
                    "selector_score_definition": contract["score_definition"],
                    "selection_mode": mode,
                    "selector_bits": bits,
                    "candidates_per_payload": width,
                    "base_library_sequences": len(library),
                    "payloads": fibers,
                    "emitted_payloads": int(accepted.sum()),
                    "failed_payloads": int((~accepted).sum()),
                    "payload_retention_fraction": float(accepted.mean()),
                    "base_payload_bits": base_payload_bits,
                    "payload_bits_after_selection": base_payload_bits - bits,
                    "payload_bits_per_variable_nt": (base_payload_bits - bits)
                    / codec.GENERATED_LENGTH,
                    "target_library_mean_outcome": float(outcome_all.mean()),
                    "selected_mean_outcome": float(accepted_gains.mean() + fiber_mean[accepted].mean()),
                    **stats,
                    "target_library_event_rate": float(event_all.mean()),
                    "selected_event_rate": float(accepted_events.mean()),
                    "paired_event_rate_change": float(accepted_event_changes.mean()),
                    "paired_event_rate_change_ci_2p5": event_low,
                    "paired_event_rate_change_ci_97p5": event_high,
                    "paired_event_rate_change_bca_ci_2p5": event_bca_low,
                    "paired_event_rate_change_bca_ci_97p5": event_bca_high,
                    "fiber_score_variance_mean": fiber_score_variance_mean,
                    "fiber_outcome_variance_mean": float(
                        outcomes.var(axis=1, ddof=0)[accepted].mean()
                    ),
                    "score_decimal_places": deterministic.SCORE_DECIMAL_PLACES
                    if mode in {"max", "min"}
                    else "",
                    **score_stability,
                    "gain_estimand": "conditional_on_emitted_payloads"
                    if not accepted.all()
                    else "all_nonoverlapping_payload_fibers",
                    "independent_unit": "one deterministic non-overlapping candidate fiber",
                    "evidence_boundary": contract["evidence_boundary"],
                }
            )
            for fiber in range(fibers):
                rank = int(fiber * width + choices[fiber])
                distribution_rows.append(
                    {
                        "dataset": dataset,
                        "direction": direction,
                        "selector_model": name,
                        "selector_bits": bits,
                        "fiber_id": fiber,
                        "accepted": bool(accepted[fiber]),
                        "choice_index": int(choices[fiber]),
                        "selected_base_rank": rank,
                        "sequence_sha256": library.iloc[rank]["sequence_sha256"],
                        "fiber_mean_outcome": float(fiber_mean[fiber]),
                        "selected_outcome": float(selected_outcome[fiber])
                        if accepted[fiber]
                        else math.nan,
                        "outcome_gain": float(gains[fiber]) if accepted[fiber] else math.nan,
                        "fiber_mean_event": float(fiber_event_mean[fiber]),
                        "selected_event": int(selected_event[fiber]) if accepted[fiber] else "",
                        "event_change": float(event_changes[fiber]) if accepted[fiber] else math.nan,
                        "fiber_score_variance": float(fiber_score_variance[fiber])
                        if mode in {"max", "min"}
                        else math.nan,
                        "selected_score": float(selected_score[fiber])
                        if mode in {"max", "min"}
                        else math.nan,
                        "runner_up_score": float(runner_up_score[fiber])
                        if mode in {"max", "min"} and width > 1
                        else math.nan,
                        "score_margin": float(score_margin[fiber])
                        if mode in {"max", "min"} and width > 1
                        else math.nan,
                    }
                )
            if accepted.all() and mode in {"max", "min"} and bits in INFERENTIAL_BITS:
                test = random_choice_test(
                    outcomes,
                    float(gains.mean()),
                    randomizations,
                    BASE_SEED
                    + 5_000_000
                    + seed_offset
                    + bits_index * 10_000
                    + model_index * 100,
                )
                randomization_rows.append(
                    {
                        "dataset": dataset,
                        "direction": direction,
                        "selector_model": name,
                        "selector_bits": bits,
                        "fibers": fibers,
                        "observed_mean_gain": float(gains.mean()),
                        "randomization_replicates": randomizations,
                        **test,
                        "null": "one uniformly random candidate per fixed fiber",
                    }
                )

        for comparison_index, (left, right) in enumerate(
            [
                ("P5_combined_context", "P2_assay_context"),
                ("P5_combined_context", "published_1dcnn"),
            ]
        ):
            if left not in selected_by_model or right not in selected_by_model:
                continue
            difference = selected_by_model[left] - selected_by_model[right]
            difference = difference[np.isfinite(difference)]
            event_difference = selected_event_by_model[left] - selected_event_by_model[right]
            event_difference = event_difference[np.isfinite(event_difference)]
            distribution = bootstrap_distribution(
                difference,
                bootstrap,
                BASE_SEED
                + 7_000_000
                + seed_offset
                + bits_index * 10_000
                + comparison_index,
            )
            low, high = percentile_interval(distribution)
            event_distribution = bootstrap_distribution(
                event_difference,
                bootstrap,
                BASE_SEED
                + 7_500_000
                + seed_offset
                + bits_index * 10_000
                + comparison_index,
            )
            event_low, event_high = percentile_interval(event_distribution)
            comparison_rows.append(
                {
                    "dataset": dataset,
                    "direction": direction,
                    "selector_bits": bits,
                    "comparison": f"{left}_minus_{right}",
                    "paired_selected_outcome_difference": float(difference.mean()),
                    "difference_ci_2p5": low,
                    "difference_ci_97p5": high,
                    "standardized_difference": float(difference.mean()) / outcome_sd,
                    "compared_fibers": len(difference),
                    "paired_selected_event_risk_difference": float(
                        event_difference.mean()
                    ),
                    "event_risk_difference_ci_2p5": event_low,
                    "event_risk_difference_ci_97p5": event_high,
                    "event_difference_interpretation": "negative favors P5 for the low-efficiency endpoint",
                }
            )
    return summary_rows, distribution_rows, randomization_rows, comparison_rows


def fit_bootstrap_model(
    source_features: np.ndarray,
    source_outcome: np.ndarray,
    alpha_seed: int,
    source_group_ids: np.ndarray,
    original_group_count: int,
) -> GridSearchCV:
    test_fold = fixed_group_test_fold(
        source_group_ids,
        original_group_count=original_group_count,
        alpha_seed=alpha_seed,
    )
    estimator = Pipeline([("scale", StandardScaler()), ("model", Ridge())])
    search = GridSearchCV(
        estimator,
        {"model__alpha": base.RIDGE_ALPHA_GRID},
        scoring="neg_mean_squared_error",
        cv=PredefinedSplit(test_fold),
        n_jobs=1,
        refit=True,
    )
    search.fit(source_features, source_outcome)
    return search


@dataclass(frozen=True)
class FastRidgeCV:
    alpha: float
    scaler_mean: np.ndarray
    scaler_scale: np.ndarray
    ridge_mean: np.ndarray
    coefficient: np.ndarray
    outcome_mean: float

    def predict(self, features: np.ndarray) -> np.ndarray:
        scaled = (np.asarray(features, dtype=float) - self.scaler_mean) / self.scaler_scale
        return (scaled - self.ridge_mean) @ self.coefficient + self.outcome_mean


def standardized_ridge_components(
    features: np.ndarray, outcome: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    features = np.asarray(features, dtype=float)
    outcome = np.asarray(outcome, dtype=float)
    scaler_mean = features.mean(axis=0)
    scaler_scale = features.std(axis=0, ddof=0)
    scaler_scale = np.where(scaler_scale == 0, 1.0, scaler_scale)
    scaled = (features - scaler_mean) / scaler_scale
    ridge_mean = scaled.mean(axis=0)
    centered = scaled - ridge_mean
    outcome_mean = float(outcome.mean())
    return scaler_mean, scaler_scale, ridge_mean, centered, outcome_mean


def fixed_group_test_fold(
    source_group_ids: np.ndarray,
    original_group_count: int,
    alpha_seed: int,
) -> np.ndarray:
    """Assign bootstrap rows to fixed original-sequence CV folds.

    ``source_group_ids`` contains the original row identifier for every row in
    a source bootstrap sample.  All copies of the same original sequence are
    assigned to one fold generated before resampling.  Explicit duplicate rows
    therefore implement integer bootstrap weights without allowing an original
    sequence to occur in both training and validation within a fold.
    """

    groups = np.asarray(source_group_ids)
    if groups.ndim != 1:
        raise ValueError("source_group_ids must be one-dimensional")
    if not np.issubdtype(groups.dtype, np.integer):
        if not np.all(np.isfinite(groups)) or not np.all(groups == np.floor(groups)):
            raise ValueError("source_group_ids must contain integer identifiers")
        groups = groups.astype(np.int64)
    else:
        groups = groups.astype(np.int64, copy=False)
    if original_group_count < 5:
        raise ValueError("at least five original source groups are required")
    if len(groups) == 0 or int(groups.min()) < 0 or int(groups.max()) >= original_group_count:
        raise ValueError("source_group_ids fall outside the original source domain")

    group_to_fold = np.full(original_group_count, -1, dtype=np.int8)
    splitter = KFold(n_splits=5, shuffle=True, random_state=alpha_seed)
    original_ids = np.arange(original_group_count, dtype=np.int64)
    for fold, (_, validation_groups) in enumerate(splitter.split(original_ids)):
        group_to_fold[validation_groups] = fold
    if np.any(group_to_fold < 0):
        raise RuntimeError("fixed source-group fold assignment is incomplete")

    test_fold = group_to_fold[groups].astype(int, copy=False)
    if len(np.unique(test_fold)) != 5:
        raise RuntimeError("a source bootstrap replicate omitted an entire CV fold")
    for fold in range(5):
        train_groups = np.unique(groups[test_fold != fold])
        validation_groups = np.unique(groups[test_fold == fold])
        if np.intersect1d(train_groups, validation_groups).size:
            raise RuntimeError("an original source sequence crossed CV folds")
    return test_fold


def fit_fast_ridge_cv(
    features: np.ndarray,
    outcome: np.ndarray,
    alpha_seed: int,
    source_group_ids: np.ndarray | None = None,
    original_group_count: int | None = None,
) -> FastRidgeCV:
    """Exact contract-equivalent acceleration of StandardScaler + Ridge grid CV.

    Standardization and centering are reused across alpha values within each
    fixed fold.  This removes repeated preprocessing/model-dispatch overhead but
    preserves the nine-alpha, five-fold mean-squared-error selection rule.
    Equivalence to the original scikit-learn pipeline is audited before use.
    """

    features = np.asarray(features, dtype=float)
    outcome = np.asarray(outcome, dtype=float)
    alphas = np.asarray(base.RIDGE_ALPHA_GRID, dtype=float)
    fold_mse = np.zeros((len(alphas), 5), dtype=float)
    if source_group_ids is None:
        splitter = KFold(n_splits=5, shuffle=True, random_state=alpha_seed)
        splits = list(splitter.split(features))
    else:
        if original_group_count is None:
            raise ValueError("original_group_count is required with source_group_ids")
        test_fold = fixed_group_test_fold(
            source_group_ids,
            original_group_count=original_group_count,
            alpha_seed=alpha_seed,
        )
        splits = [
            (np.flatnonzero(test_fold != fold), np.flatnonzero(test_fold == fold))
            for fold in range(5)
        ]
    identity = np.eye(features.shape[1], dtype=float)
    for fold, (train_index, validation_index) in enumerate(splits):
        train_x = features[train_index]
        train_y = outcome[train_index]
        mean, scale, ridge_mean, centered, outcome_mean = standardized_ridge_components(
            train_x, train_y
        )
        centered_y = train_y - outcome_mean
        xtx = centered.T @ centered
        xty = centered.T @ centered_y
        validation_scaled = (features[validation_index] - mean) / scale - ridge_mean
        validation_y = outcome[validation_index]
        for alpha_index, alpha in enumerate(alphas):
            coefficient = np.linalg.solve(xtx + alpha * identity, xty)
            prediction = validation_scaled @ coefficient + outcome_mean
            fold_mse[alpha_index, fold] = float(np.mean((validation_y - prediction) ** 2))
    selected_index = int(np.argmin(fold_mse.mean(axis=1)))
    selected_alpha = float(alphas[selected_index])
    mean, scale, ridge_mean, centered, outcome_mean = standardized_ridge_components(
        features, outcome
    )
    centered_y = outcome - outcome_mean
    coefficient = np.linalg.solve(
        centered.T @ centered + selected_alpha * identity,
        centered.T @ centered_y,
    )
    return FastRidgeCV(
        alpha=selected_alpha,
        scaler_mean=mean,
        scaler_scale=scale,
        ridge_mean=ridge_mean,
        coefficient=coefficient,
        outcome_mean=outcome_mean,
    )


def verify_fast_ridge_cv(
    pcr: pd.DataFrame,
    models: dict[tuple[str, str], codec.FrozenRidge],
    public_libraries: dict[str, pd.DataFrame],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for source_index, (source_pool, direction) in enumerate(
        [("GCall", "GCall_to_GCfix"), ("GCfix", "GCfix_to_GCall")]
    ):
        source = pcr.loc[pcr["pool"].eq(source_pool)].reset_index(drop=True)
        target = public_libraries[direction]
        source_outcome = source["eff"].to_numpy(float)
        bootstrap_rng = np.random.default_rng(
            BASE_SEED + 10_000_000 + source_index * 1_000_000
        )
        sampled = bootstrap_rng.integers(0, len(source), size=len(source))
        for model_name in TWO_STAGE_MODELS:
            names = list(codec.MODEL_FEATURES[model_name])
            alpha_seed = (
                selection.BASE_SEED
                + 20_000
                + source_index * 1000
                + codec.MODEL_ORIGINAL_INDEX[model_name]
            )
            source_features = source[names].to_numpy(float)
            target_features = target[names].to_numpy(float)
            fast_full = fit_fast_ridge_cv(source_features, source_outcome, alpha_seed)
            frozen = models[(source_pool, model_name)]
            full_delta = float(
                np.max(np.abs(fast_full.predict(target_features) - frozen.predict(target)))
            )
            if fast_full.alpha != frozen.alpha or full_delta > 2e-10:
                raise ValueError(
                    f"fast ridge full-data mismatch {source_pool}/{model_name}: "
                    f"alpha {fast_full.alpha}/{frozen.alpha}, delta {full_delta}"
                )
            bootstrap_x = source_features[sampled]
            bootstrap_y = source_outcome[sampled]
            fast_bootstrap = fit_fast_ridge_cv(
                bootstrap_x,
                bootstrap_y,
                alpha_seed,
                source_group_ids=sampled,
                original_group_count=len(source),
            )
            sklearn_bootstrap = fit_bootstrap_model(
                bootstrap_x,
                bootstrap_y,
                alpha_seed,
                source_group_ids=sampled,
                original_group_count=len(source),
            )
            bootstrap_delta = float(
                np.max(
                    np.abs(
                        fast_bootstrap.predict(target_features)
                        - sklearn_bootstrap.predict(target_features)
                    )
                )
            )
            sklearn_alpha = float(sklearn_bootstrap.best_params_["model__alpha"])
            if fast_bootstrap.alpha != sklearn_alpha or bootstrap_delta > 2e-10:
                raise ValueError(
                    f"fast ridge bootstrap mismatch {source_pool}/{model_name}: "
                    f"alpha {fast_bootstrap.alpha}/{sklearn_alpha}, delta {bootstrap_delta}"
                )
            rows.append(
                {
                    "artifact": "fast_ridge_cv_equivalence",
                    "source_pool": source_pool,
                    "model": model_name,
                    "full_selected_alpha": fast_full.alpha,
                    "full_prediction_max_abs_difference": full_delta,
                    "bootstrap_replicate": 0,
                    "bootstrap_selected_alpha": fast_bootstrap.alpha,
                    "bootstrap_prediction_max_abs_difference": bootstrap_delta,
                    "reference": "scikit-learn StandardScaler + GridSearchCV + Ridge",
                    "bootstrap_cv_contract": "fixed original-sequence folds; duplicate copies never cross folds",
                    "bootstrap_unique_source_sequences": int(np.unique(sampled).size),
                    "bootstrap_duplicate_rows": int(len(sampled) - np.unique(sampled).size),
                    "bootstrap_cross_fold_duplicate_sequences": 0,
                }
            )
    return rows


def two_stage_one(
    source_pool: str,
    source_index: int,
    replicate: int,
    source: pd.DataFrame,
    targets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(BASE_SEED + 10_000_000 + source_index * 1_000_000 + replicate)
    sampled = rng.integers(0, len(source), size=len(source))
    source_outcome = source["eff"].to_numpy(float)[sampled]
    fitted: dict[str, FastRidgeCV] = {}
    for model_name in TWO_STAGE_MODELS:
        feature_names = list(codec.MODEL_FEATURES[model_name])
        source_features = source[feature_names].to_numpy(float)[sampled]
        alpha_seed = (
            selection.BASE_SEED
            + 20_000
            + source_index * 1000
            + codec.MODEL_ORIGINAL_INDEX[model_name]
        )
        fitted[model_name] = fit_fast_ridge_cv(
            source_features,
            source_outcome,
            alpha_seed,
            source_group_ids=sampled,
            original_group_count=len(source),
        )

    rows: list[dict[str, Any]] = []
    for target_index, target in enumerate(targets):
        predictions: dict[str, np.ndarray] = {}
        for model_name, model in fitted.items():
            features = target["frame"][list(codec.MODEL_FEATURES[model_name])].to_numpy(float)
            predictions[model_name] = model.predict(features)
        outcome = target["outcome"]
        for bits in INFERENTIAL_BITS:
            width = 1 << bits
            fibers = len(outcome) // width
            outcomes = outcome.reshape(fibers, width)
            fiber_mean = outcomes.mean(axis=1)
            selected: dict[str, np.ndarray] = {}
            gains: dict[str, np.ndarray] = {}
            for model_name in TWO_STAGE_MODELS:
                score = predictions[model_name].reshape(fibers, width)
                choices = deterministic.argmax_smallest(score, axis=1)
                selected[model_name] = outcomes[np.arange(fibers), choices]
                gains[model_name] = selected[model_name] - fiber_mean
            cnn_score = np.asarray(target["cnn_values"], dtype=float).reshape(
                fibers, width
            )
            if target["cnn_mode"] == "max":
                cnn_choices = deterministic.argmax_smallest(cnn_score, axis=1)
            elif target["cnn_mode"] == "min":
                cnn_choices = deterministic.argmin_smallest(cnn_score, axis=1)
            else:
                raise ValueError(f"unsupported fixed CNN selector mode {target['cnn_mode']}")
            cnn_selected = outcomes[np.arange(fibers), cnn_choices]
            target_rng = np.random.default_rng(
                BASE_SEED
                + 20_000_000
                + source_index * 1_000_000
                + target_index * 100_000
                + bits * 10_000
                + replicate
            )
            sampled_fibers = target_rng.integers(0, fibers, size=fibers)
            for model_name in TWO_STAGE_MODELS:
                rows.append(
                    {
                        "source_pool": source_pool,
                        "target_dataset": target["dataset"],
                        "direction": target["direction"],
                        "replicate": replicate,
                        "selector_bits": bits,
                        "estimand": f"{model_name}_mean_gain",
                        "estimate": float(gains[model_name][sampled_fibers].mean()),
                        "selected_alpha": fitted[model_name].alpha,
                        "source_unique_sequences": int(np.unique(sampled).size),
                        "source_cv_contract": "fixed original-sequence folds; bootstrap copies remain grouped",
                    }
                )
            difference = (
                selected["P5_combined_context"] - selected["P2_assay_context"]
            )
            rows.append(
                {
                    "source_pool": source_pool,
                    "target_dataset": target["dataset"],
                    "direction": target["direction"],
                    "replicate": replicate,
                    "selector_bits": bits,
                    "estimand": "P5_minus_P2_selected_outcome",
                    "estimate": float(difference[sampled_fibers].mean()),
                    "selected_alpha": "joint_models",
                    "source_unique_sequences": int(np.unique(sampled).size),
                    "source_cv_contract": "fixed original-sequence folds; bootstrap copies remain grouped",
                }
            )
            p5_minus_cnn = selected["P5_combined_context"] - cnn_selected
            rows.append(
                {
                    "source_pool": source_pool,
                    "target_dataset": target["dataset"],
                    "direction": target["direction"],
                    "replicate": replicate,
                    "selector_bits": bits,
                    "estimand": "P5_minus_released_1dcnn_selected_outcome",
                    "estimate": float(p5_minus_cnn[sampled_fibers].mean()),
                    "selected_alpha": fitted["P5_combined_context"].alpha,
                    "source_unique_sequences": int(np.unique(sampled).size),
                    "source_cv_contract": "fixed original-sequence folds; bootstrap copies remain grouped",
                }
            )
    return rows


def two_stage_chunk(
    source_pool: str,
    source_index: int,
    replicates: list[int],
    source: pd.DataFrame,
    targets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for replicate in replicates:
        rows.extend(two_stage_one(source_pool, source_index, replicate, source, targets))
    return rows


def two_stage_bootstrap(
    pcr: pd.DataFrame,
    public_libraries: dict[str, pd.DataFrame],
    external_library: pd.DataFrame,
    external_q5_library: pd.DataFrame,
    replicates: int,
    jobs: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    tasks: list[tuple[str, int, list[int], pd.DataFrame, list[dict[str, Any]]]] = []
    for source_index, (source_pool, public_direction) in enumerate(
        [("GCall", "GCall_to_GCfix"), ("GCfix", "GCfix_to_GCall")]
    ):
        source = pcr.loc[pcr["pool"].eq(source_pool)].reset_index(drop=True)
        public_library = public_libraries[public_direction]
        targets = [
            {
                "dataset": "Gimpel2025_cross_pool_public_codebook",
                "direction": public_direction,
                "frame": public_library,
                "outcome": public_library["eff"].to_numpy(float),
                "cnn_values": public_library["published_1dcnn"].to_numpy(float),
                "cnn_mode": "max",
            },
            {
                "dataset": "Gimpel2025_external_laboratory_Taq",
                "direction": f"{source_pool}_to_external_Taq",
                "frame": external_library,
                "outcome": external_library["eff_Taq"].to_numpy(float),
                "cnn_values": external_library[f"{source_pool} 2perc"].to_numpy(float),
                "cnn_mode": "min",
            },
            {
                "dataset": "Gimpel2025_external_laboratory_Q5_sensitivity",
                "direction": f"{source_pool}_to_external_Q5",
                "frame": external_q5_library,
                "outcome": external_q5_library["eff_Q5"].to_numpy(float),
                "cnn_values": external_q5_library[f"{source_pool} 2perc"].to_numpy(float),
                "cnn_mode": "min",
            },
        ]
        chunks = [
            chunk.astype(int).tolist()
            for chunk in np.array_split(np.arange(replicates, dtype=int), min(jobs, replicates))
            if len(chunk)
        ]
        tasks.extend((source_pool, source_index, chunk, source, targets) for chunk in chunks)
    nested = Parallel(n_jobs=jobs, prefer="processes", batch_size=1, verbose=10)(
        delayed(two_stage_chunk)(*task) for task in tasks
    )
    rows = [row for block in nested for row in block]
    frame = pd.DataFrame(rows)
    summary: list[dict[str, Any]] = []
    group_columns = [
        "source_pool",
        "target_dataset",
        "direction",
        "selector_bits",
        "estimand",
    ]
    for keys, group in frame.groupby(group_columns, sort=True):
        estimates = group["estimate"].to_numpy(float)
        low, high = np.quantile(estimates, [0.025, 0.975])
        alpha_values = pd.to_numeric(group["selected_alpha"], errors="coerce").dropna()
        alpha_counts = (
            alpha_values.value_counts().sort_index().to_dict() if len(alpha_values) else {}
        )
        summary.append(
            {
                **dict(zip(group_columns, keys, strict=True)),
                "replicates": len(group),
                "two_stage_mean_estimate": float(estimates.mean()),
                "two_stage_median_estimate": float(np.median(estimates)),
                "two_stage_ci_2p5": float(low),
                "two_stage_ci_97p5": float(high),
                "fraction_above_zero": float((estimates > 0).mean()),
                "selected_alpha_counts": json.dumps(alpha_counts, sort_keys=True),
                "resampling_unit_source": "one source sequence",
                "resampling_unit_target": "one fixed candidate fiber",
            }
        )
    return rows, summary


def build_public_libraries(
    pcr: pd.DataFrame,
    models: dict[tuple[str, str], codec.FrozenRidge],
) -> tuple[dict[str, pd.DataFrame], list[dict[str, Any]]]:
    libraries: dict[str, pd.DataFrame] = {}
    audits: list[dict[str, Any]] = []
    prediction_files = {
        "GCall_to_GCfix": "External_GCall2GCfix_2perc_regression_plus_probs.csv",
        "GCfix_to_GCall": "External_GCfix2GCall_2perc_regression_plus_probs.csv",
    }
    for source_pool, target_pool in [("GCall", "GCfix"), ("GCfix", "GCall")]:
        direction = f"{source_pool}_to_{target_pool}"
        target = pcr.loc[pcr["pool"].eq(target_pool)].reset_index(drop=True)
        published, prediction_audit = load_cross_pool_prediction(
            INPUT_DIR / prediction_files[direction], target, direction
        )
        properties, property_audit = load_sequence_properties(target_pool, target)
        library = codec.deterministic_library(target, direction)
        library["published_1dcnn"] = library["sample_id"].map(
            published["pred_efficiency"]
        )
        library["published_1dcnn_low_probability"] = library["sample_id"].map(
            published["pred_probability"]
        )
        library[["published_GC", "published_hp", "published_dg"]] = library[
            "sample_id"
        ].map(properties.to_dict("index")).apply(pd.Series)[["GC", "hp", "dg"]].to_numpy()
        library["published_sota_hard_pass"] = (
            (library["published_GC"] > 0.4)
            & (library["published_GC"] < 0.6)
            & (library["published_hp"] < 5)
            & (library["published_dg"] > -15)
        )
        for model_name in MODEL_NAMES:
            library[model_name] = models[(source_pool, model_name)].predict(library)
        if library[
            ["published_1dcnn", "published_GC", "published_hp", "published_dg"]
        ].isna().any().any():
            raise ValueError(f"{direction}: missing published comparator value")
        audits.extend([prediction_audit, property_audit])
        audits.append(
            {
                "artifact": "public_identical_fiber_library",
                "direction": direction,
                "target_pool": target_pool,
                "eligible_before_hash_selection": int(
                    target["sequence"].map(codec.exact_gc_hp_eligible).sum()
                ),
                "library_sequences": len(library),
                "published_hard_filter_pass_sequences": int(
                    library["published_sota_hard_pass"].sum()
                ),
                "published_hard_filter_pass_fraction": float(
                    library["published_sota_hard_pass"].mean()
                ),
            }
        )
        libraries[direction] = library
    return libraries, audits


def output_manifest(output_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(output_dir.iterdir()):
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
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    bootstrap = 50 if args.smoke else args.bootstrap
    two_stage = 20 if args.smoke else args.two_stage
    randomizations = 2000 if args.smoke else args.randomizations
    output_dir = OUT_DIR / "smoke" if args.smoke else OUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    input_manifest = validate_inputs()
    public_provenance = selection.validate_public_sources()
    _dt, pcr, pickle_audit = selection.load_frames(args.jobs)
    models, coefficient_rows = codec.fit_frozen_models(pcr, args.jobs)
    public_libraries, public_audits = build_public_libraries(pcr, models)
    public_audits.extend(verify_fast_ridge_cv(pcr, models, public_libraries))

    external_random, external_audits = load_external_pool()
    external_library = locked_external_library(external_random)
    feature_frame = external_features(
        external_library, output_dir, args.jobs, args.force_features
    )
    feature_frame = feature_frame.rename(columns={"sample_id": "seq_id_anonymized"})
    external_library = external_library.merge(
        feature_frame.drop(columns=["pool"]),
        on=["seq_id_anonymized", "sequence_sha256"],
        how="left",
        validate="one_to_one",
    )
    if external_library[list(codec.MODEL_FEATURES["P5_combined_context"])].isna().any().any():
        raise ValueError("external library feature join is incomplete")

    external_q5_random, external_q5_audits = load_external_q5_pool()
    external_q5_library = locked_external_library(
        external_q5_random,
        endpoint="external_Q5",
        library_size=1024,
    )
    q5_feature_frame = external_features(
        external_q5_library,
        output_dir,
        args.jobs,
        args.force_features,
        cache_name="external_q5_locked_library_features.tsv",
        pool_label="External_Q5_sensitivity",
    ).rename(columns={"sample_id": "seq_id_anonymized"})
    external_q5_library = external_q5_library.merge(
        q5_feature_frame.drop(columns=["pool"]),
        on=["seq_id_anonymized", "sequence_sha256"],
        how="left",
        validate="one_to_one",
    )
    if external_q5_library[list(codec.MODEL_FEATURES["P5_combined_context"])].isna().any().any():
        raise ValueError("external Q5 library feature join is incomplete")

    all_summary: list[dict[str, Any]] = []
    all_distributions: list[dict[str, Any]] = []
    all_randomizations: list[dict[str, Any]] = []
    all_comparisons: list[dict[str, Any]] = []
    for direction_index, (source_pool, target_pool) in enumerate(
        [("GCall", "GCfix"), ("GCfix", "GCall")]
    ):
        direction = f"{source_pool}_to_{target_pool}"
        library = public_libraries[direction]
        selectors: OrderedDict[str, dict[str, Any]] = OrderedDict()
        for model_name in MODEL_NAMES:
            selectors[model_name] = {
                "mode": "max",
                "values": library[model_name].to_numpy(float),
                "selector_class": "Paper2_internal_matched_control"
                if model_name != "P5_combined_context"
                else "Paper2_primary_selector",
                "evidence_boundary": "source-only fitted Paper 2 ridge score on retrospective measured public fibers",
                "score_definition": "Paper 2 continuous relative-PCR-efficiency ridge score",
            }
        selectors["published_1dcnn"] = {
            "mode": "max",
            "values": library["published_1dcnn"].to_numpy(float),
            "selector_class": "external_published_SOTA_1D_CNN",
            "evidence_boundary": "fixed upstream cross-pool 1D-CNN regression prediction; target outcomes not refitted by Paper 2",
            "score_definition": "released 1D-CNN regression prediction",
        }
        selectors["published_sota_hard_filter"] = {
            "mode": "first_pass",
            "values": library["published_sota_hard_pass"].to_numpy(bool),
            "selector_class": "external_published_SOTA_hard_filter",
            "evidence_boundary": "published GC/homopolymer/mfold rule with deterministic smallest-rank tie-break; conditional if a fiber has no passing candidate",
            "score_definition": "released hard pass rule; smallest-rank tie-break among passers",
        }
        result = evaluate_library(
            "Gimpel2025_cross_pool_public_codebook",
            direction,
            library,
            "eff",
            "label",
            selectors,
            bootstrap,
            randomizations,
            direction_index * 1_000_000,
        )
        all_summary.extend(result[0])
        all_distributions.extend(result[1])
        all_randomizations.extend(result[2])
        all_comparisons.extend(result[3])

    for source_index, source_pool in enumerate(["GCall", "GCfix"]):
        direction = f"{source_pool}_to_external_Taq"
        selectors = OrderedDict()
        for model_name in MODEL_NAMES:
            values = models[(source_pool, model_name)].predict(external_library)
            external_library[f"{source_pool}__{model_name}"] = values
            selectors[model_name] = {
                "mode": "max",
                "values": values,
                "selector_class": "Paper2_internal_matched_control"
                if model_name != "P5_combined_context"
                else "Paper2_primary_selector",
                "evidence_boundary": "fully frozen source-pool ridge score transferred to an external-laboratory measured pool",
                "score_definition": "Paper 2 continuous relative-PCR-efficiency ridge score",
            }
        cnn_column = f"{source_pool} 2perc"
        selectors["published_1dcnn"] = {
            "mode": "min",
            "values": external_library[cnn_column].to_numpy(float),
            "selector_class": "external_published_SOTA_1D_CNN",
            "evidence_boundary": "fixed upstream probability of the positive low-efficiency class on an external-laboratory measured pool",
            "score_definition": "released 1D-CNN low-efficiency-class probability",
        }
        result = evaluate_library(
            "Gimpel2025_external_laboratory_Taq",
            direction,
            external_library,
            "eff_Taq",
            "external_low_efficiency",
            selectors,
            bootstrap,
            randomizations,
            2_000_000 + source_index * 1_000_000,
        )
        all_summary.extend(result[0])
        all_distributions.extend(result[1])
        all_randomizations.extend(result[2])
        all_comparisons.extend(result[3])

    for source_index, source_pool in enumerate(["GCall", "GCfix"]):
        direction = f"{source_pool}_to_external_Q5"
        selectors = OrderedDict()
        for model_name in MODEL_NAMES:
            values = models[(source_pool, model_name)].predict(external_q5_library)
            external_q5_library[f"{source_pool}__{model_name}"] = values
            selectors[model_name] = {
                "mode": "max",
                "values": values,
                "selector_class": "Paper2_internal_matched_control"
                if model_name != "P5_combined_context"
                else "Paper2_primary_selector",
                "evidence_boundary": "fully frozen Taq-trained source-pool ridge score evaluated on the external Q5 workflow-shift sensitivity codebook",
                "score_definition": "Paper 2 continuous Taq relative-PCR-efficiency ridge score",
            }
        cnn_column = f"{source_pool} 2perc"
        selectors["published_1dcnn"] = {
            "mode": "min",
            "values": external_q5_library[cnn_column].to_numpy(float),
            "selector_class": "external_published_SOTA_1D_CNN",
            "evidence_boundary": "fixed upstream low-efficiency-class probability evaluated against the Q5 workflow-shift endpoint",
            "score_definition": "released 1D-CNN low-efficiency-class probability",
        }
        result = evaluate_library(
            "Gimpel2025_external_laboratory_Q5_sensitivity",
            direction,
            external_q5_library,
            "eff_Q5",
            "external_low_efficiency_Q5",
            selectors,
            bootstrap,
            randomizations,
            4_000_000 + source_index * 1_000_000,
        )
        all_summary.extend(result[0])
        all_distributions.extend(result[1])
        all_randomizations.extend(result[2])
        all_comparisons.extend(result[3])

    two_stage_rows, two_stage_summary = two_stage_bootstrap(
        pcr,
        public_libraries,
        external_library,
        external_q5_library,
        two_stage,
        args.jobs,
    )

    join_audits = public_audits + external_audits + external_q5_audits
    join_audits.append(
        {
            "artifact": "locked_external_library",
            "eligible_before_hash_selection": int(
                external_random["sequence"].map(codec.exact_gc_hp_eligible).sum()
            ),
            "library_sequences": len(external_library),
            "library_key": "paper2-choice-external-v1|external_Taq|sequence_sha256",
            "library_outcome_mean": float(external_library["eff_Taq"].mean()),
            "library_outcome_sd_population": float(external_library["eff_Taq"].std(ddof=0)),
            "library_low_efficiency_events": int(
                external_library["external_low_efficiency"].sum()
            ),
        }
    )

    event_rows: list[dict[str, Any]] = []
    summary_frame_for_events = pd.DataFrame(all_summary)
    for keys, group in summary_frame_for_events.loc[
        summary_frame_for_events["selector_bits"].isin(INFERENTIAL_BITS)
        & summary_frame_for_events["selector_model"].isin(
            ["P5_combined_context", "P2_assay_context", "published_1dcnn"]
        )
    ].groupby(["dataset", "direction", "selector_bits"], sort=True):
        dataset, direction, selector_bits = keys
        first = group.iloc[0]
        event_rows.append(
            {
                "dataset": dataset,
                "direction": direction,
                "selector_bits": int(selector_bits),
                "selector_model": "uniform_random_choice_expectation",
                "selector_score_definition": "one uniformly random candidate per fixed fiber",
                "fibers": int(first["payloads"]),
                "random_choice_expected_low_efficiency_rate": float(
                    first["target_library_event_rate"]
                ),
                "selected_low_efficiency_rate": float(
                    first["target_library_event_rate"]
                ),
                "paired_low_efficiency_risk_difference": 0.0,
                "risk_difference_ci_2p5": 0.0,
                "risk_difference_ci_97p5": 0.0,
                "interpretation": "reference expectation, not an estimated selector",
            }
        )
        for row in group.itertuples(index=False):
            event_rows.append(
                {
                    "dataset": dataset,
                    "direction": direction,
                    "selector_bits": int(selector_bits),
                    "selector_model": row.selector_model,
                    "selector_score_definition": row.selector_score_definition,
                    "fibers": int(row.payloads),
                    "random_choice_expected_low_efficiency_rate": float(
                        row.target_library_event_rate
                    ),
                    "selected_low_efficiency_rate": float(row.selected_event_rate),
                    "paired_low_efficiency_risk_difference": float(
                        row.paired_event_rate_change
                    ),
                    "risk_difference_ci_2p5": float(
                        row.paired_event_rate_change_bca_ci_2p5
                    ),
                    "risk_difference_ci_97p5": float(
                        row.paired_event_rate_change_bca_ci_97p5
                    ),
                    "interpretation": "negative favors the selector for the endpoint-matched low-efficiency event",
                }
            )

    event_pairwise_rows: list[dict[str, Any]] = []
    for row in all_comparisons:
        if (
            row["comparison"] == "P5_combined_context_minus_published_1dcnn"
            and int(row["selector_bits"]) in INFERENTIAL_BITS
        ):
            event_pairwise_rows.append(
                {
                    "dataset": row["dataset"],
                    "direction": row["direction"],
                    "selector_bits": int(row["selector_bits"]),
                    "comparison": "FullContext_ridge_minus_released_1dcnn_on_low_efficiency_event",
                    "paired_selected_event_risk_difference": row[
                        "paired_selected_event_risk_difference"
                    ],
                    "event_risk_difference_ci_2p5": row[
                        "event_risk_difference_ci_2p5"
                    ],
                    "event_risk_difference_ci_97p5": row[
                        "event_risk_difference_ci_97p5"
                    ],
                    "compared_fibers": row["compared_fibers"],
                    "interpretation": "negative favors FullContext ridge; endpoint-matched paired comparison",
                }
            )
    join_audits.append(
        {
            "artifact": "locked_external_Q5_sensitivity_library",
            "eligible_before_hash_selection": int(
                external_q5_random["sequence"].map(codec.exact_gc_hp_eligible).sum()
            ),
            "library_sequences": len(external_q5_library),
            "library_key": "paper2-choice-external-q5-sensitivity-v1|external_Q5|sequence_sha256",
            "library_outcome_mean": float(external_q5_library["eff_Q5"].mean()),
            "library_outcome_sd_population": float(
                external_q5_library["eff_Q5"].std(ddof=0)
            ),
            "library_low_efficiency_events": int(
                external_q5_library["external_low_efficiency_Q5"].sum()
            ),
            "endpoint_role": "workflow-shift sensitivity; separate outcome-blind 1,024-sequence codebook",
        }
    )

    write_table(output_dir / "input_manifest.tsv", input_manifest)
    write_table(output_dir / "input_join_audit.tsv", join_audits)
    write_table(output_dir / "source_model_coefficients.tsv", coefficient_rows)
    write_table(output_dir / "fiber_benchmark_summary.tsv", all_summary)
    write_table(output_dir / "fiber_gain_distributions.tsv", all_distributions)
    write_table(output_dir / "fixed_model_randomization_tests.tsv", all_randomizations)
    write_table(output_dir / "fixed_model_pairwise_comparisons.tsv", all_comparisons)
    write_table(output_dir / "endpoint_matched_low_efficiency_summary.tsv", event_rows)
    write_table(
        output_dir / "endpoint_matched_low_efficiency_pairwise.tsv",
        event_pairwise_rows,
    )
    write_table(output_dir / "two_stage_bootstrap_replicates.tsv", two_stage_rows)
    write_table(output_dir / "two_stage_bootstrap_summary.tsv", two_stage_summary)

    environment = {
        "analysis": "Paper 2 SOTA, locked external validation and two-stage uncertainty",
        "command": (
            "$PYTHON analysis_tools/validate_paper2_sota_external_uncertainty.py "
            f"--jobs {args.jobs} --bootstrap {bootstrap} --two-stage {two_stage} "
            f"--randomizations {randomizations}"
            + (" --smoke" if args.smoke else "")
        ),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "pandas", "scipy", "scikit-learn", "joblib"]
        },
        "upstream_commit": UPSTREAM_COMMIT,
        "base_seed": BASE_SEED,
        "target_fiber_bootstrap_replicates": bootstrap,
        "two_stage_bootstrap_replicates": two_stage,
        "two_stage_source_resampling": {
            "unit": "one original source sequence",
            "scheme": "nonparametric bootstrap represented by duplicate-row integer weights",
            "cross_validation": "five folds fixed on original source-sequence IDs before resampling",
            "duplicate_handling": "all bootstrap copies of one original sequence remain in the same fold",
            "cross_fold_duplicate_sequences": 0,
        },
        "measurement_uncertainty_boundary": (
            "intervals condition on the published sequence-level outcome estimates and do not "
            "propagate technical-replicate, batch, normalization or source-assay measurement "
            "uncertainty"
        ),
        "within_fiber_randomization_replicates": randomizations,
        "score_key_contract": {
            "input_dtype": "float64",
            "output_dtype": "int64",
            "decimal_places": deterministic.SCORE_DECIMAL_PLACES,
            "integer_scale": deterministic.SCORE_SCALE,
            "rounding": deterministic.ROUNDING_CONTRACT,
            "signed_int64_min": deterministic.INT64_MIN,
            "signed_int64_max": deterministic.INT64_MAX,
            "candidate_order": "lexicographically minimize (-quantized_key, choice_index)",
            "tie_break": "smallest declared choice index",
            "nonfinite_behavior": "reject",
            "out_of_range_behavior": "reject",
        },
        "jobs": args.jobs,
        "smoke": args.smoke,
        "frozen_full_run": bool(
            not args.smoke
            and bootstrap == BOOTSTRAP_REPLICATES
            and two_stage == TWO_STAGE_REPLICATES
            and randomizations == RANDOMIZATION_REPLICATES
        ),
        "primary_selector_bits": 2,
        "secondary_selector_bits": 4,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "environment_and_seeds.json", environment)
    write_json(
        output_dir / "source_provenance.json",
        {
            "input_manifest": input_manifest,
            "public_source_validation": public_provenance,
            "pickle_opcode_audit": pickle_audit,
            "join_audits": join_audits,
            "external_protocol": relative(
                PAPER_DIR
                / "bioinformatics_reframe"
                / "LOCKED_EXTERNAL_VALIDATION_PROTOCOL_20260718.md"
            ),
            "published_hard_filter": {
                "GC": "strictly greater than 0.4 and strictly less than 0.6",
                "maximum_homopolymer": "less than 5",
                "mfold_free_energy_kcal_per_mol": "greater than -15",
                "source": "fixed upstream filtering notebook and publication Methods",
            },
            "claim_boundary": CLAIM_BOUNDARY,
        },
    )

    summary_frame = pd.DataFrame(all_summary)
    two_frame = pd.DataFrame(two_stage_summary)
    lines = [
        "# Paper 2 SOTA and external-validation analysis",
        "",
        f"Run class: {'SMOKE' if args.smoke else 'FROZEN FULL'}",
        "",
        "## Fixed-model identical-fiber results",
        "",
    ]
    for dataset in [
        "Gimpel2025_cross_pool_public_codebook",
        "Gimpel2025_external_laboratory_Taq",
        "Gimpel2025_external_laboratory_Q5_sensitivity",
    ]:
        lines.append(f"### {dataset}")
        lines.append("")
        subset = summary_frame.loc[
            summary_frame["dataset"].eq(dataset)
            & summary_frame["selector_bits"].isin(INFERENTIAL_BITS)
            & summary_frame["selector_model"].isin(
                ["P5_combined_context", "published_1dcnn", "published_sota_hard_filter"]
            )
        ]
        for row in subset.itertuples(index=False):
            lines.append(
                f"- {row.direction}, r={row.selector_bits}, {row.selector_model}: "
                f"gain={row.paired_mean_outcome_gain:.8f} "
                f"({row.standardized_paired_gain:.4f} SD), 95% fixed-model fiber CI "
                f"[{row.paired_gain_ci_2p5:.8f}, {row.paired_gain_ci_97p5:.8f}], "
                f"payload retention={row.payload_retention_fraction:.4f}."
            )
        lines.append("")
    lines.extend(["## Two-stage source-plus-target intervals", ""])
    for row in two_frame.loc[
        two_frame["estimand"].eq("P5_combined_context_mean_gain")
    ].itertuples(index=False):
        lines.append(
            f"- {row.direction}, r={row.selector_bits}: mean={row.two_stage_mean_estimate:.8f}, "
            f"95% two-stage interval [{row.two_stage_ci_2p5:.8f}, "
            f"{row.two_stage_ci_97p5:.8f}], fraction above zero={row.fraction_above_zero:.4f}."
        )
    lines.extend(
        [
            "",
            "## Statistical boundary",
            "",
            "Source bootstrap copies retained fixed original-sequence CV folds, so no original sequence occurred in both training and validation within a replicate. Intervals condition on the published sequence-level outcomes and do not propagate technical-replicate, batch, normalization or source-assay measurement uncertainty.",
            "",
            "## Evidence boundary",
            "",
            CLAIM_BOUNDARY + ".",
            "",
        ]
    )
    (output_dir / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")

    manifest = output_manifest(output_dir)
    write_table(output_dir / "sha256_manifest.tsv", manifest)
    print(f"SOTA/external audit complete: {output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
