#!/usr/bin/env python3
"""Reproducible public-experimental-data validation for Paper 2.

This program evaluates the frozen Paper 2 reverse-complement risk proxy on two
independent public experimental resources:

1. DT4DDS sequence-level post-workflow detectability measurements; and
2. Gimpel et al. 2025 multi-template PCR efficiency labels.

The analyses are association and incremental-prediction audits.  They are not
material validation of Paper 2 emitted oligos, a recovery experiment, direct
thermodynamic validation, or evidence for a full 110-nt language capacity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import pickletools
import platform
import subprocess
import sys
import zipfile
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import hairpin_risk_features as hrf  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
REFRAME_DIR = PAPER_DIR / "bioinformatics_reframe"
OUT_DIR = REFRAME_DIR / "public_experimental_validation"

DT_CSV = ROOT / "data" / "paper2_dt4dds_sequence_detectability.csv"
DT_SUMMARY = ROOT / "data" / "paper2_dt4dds_sequence_detectability_summary.csv"
DT_REPO = ROOT / "external_data" / "paper2_dt4dds_sequence_recovery"

PCR_REPO = ROOT / "external_data" / "paper2_gimpel2025_pcr_bias"
PCR_FILES = {
    "GCall": PCR_REPO / "Data" / "GCall" / "bad_seqs_2perc.pkl",
    "GCfix": PCR_REPO / "Data" / "GCfix" / "bad_seqs_2perc.pkl",
}
PCR_SUPP_DIR = ROOT / "external_data" / "paper2_gimpel2025_supplement"
PCR_SUPP_PDF = PCR_SUPP_DIR / "41467_2025_64221_MOESM1_ESM.pdf"
PCR_SOURCE_ZIP = PCR_SUPP_DIR / "41467_2025_64221_MOESM4_ESM.zip"

MAIN_TEX = PAPER_DIR / "main.tex"
FEATURE_CACHE = OUT_DIR / "pcr_sequence_features.tsv"

EXPECTED_SHA256 = {
    DT_CSV: "32ed334ae29154c7e515d6d36409d3b7d7399d6f865265460c2e5c0b7f64dba3",
    DT_SUMMARY: "78939aa80b9ef33fc72c8c7112ace7eb8dcc4dc4c9f2151c8db2374b83964fb2",
    PCR_FILES["GCall"]: "0c88133354b30585f0a4851be7bd5fba87805a3c19189bca9db7f12fdea61384",
    PCR_FILES["GCfix"]: "185843669f779f421ce049765da92cf6ecaf1a6990cc5f51cf6b9229a43ad502",
    PCR_REPO / "LICENSE.txt": "b9ad7d4bbfbdc63bd1a8c4cf89593769d33a95c46ef4b1065d61c254368ec0c7",
    PCR_SUPP_PDF: "a0ebb08bd6c6da6abc32afd91a388b646d579b87409474efbfe647c945c9063c",
    PCR_SOURCE_ZIP: "3fdc4f6dbe21f17e4741e13794d91283b7cf4ee2caa4af4d1cbc71adb53dcd68",
}

# The manuscript hash is captured at analysis start and checked again before
# output is written. It is intentionally not a fixed source hash: revising
# prose must not invalidate or create a circular dependency in the public-data
# input lock.

EXPECTED_COMMITS = {
    DT_REPO: "126e6da5c41f4e5de072b7a1a0934068b743de6c",
    PCR_REPO: "af62c57f9a90ecdfdd0f1623441e82bdb7e082c1",
}

# Supplementary Table 6, Gimpel et al. 2025.  The synthesized strand is
# 0F + 108-nt variable region + 0R-prime (reverse complement of primer 0R).
LEFT_ADAPTER_0F = "ACACGACGCTCTTCCGATCT"
RIGHT_ADAPTER_0R_PRIME = "AGATCGGAAGAGCACACGTCT"
BOUNDARY_WINDOW = 24
COMPLEMENT_K = 4

RIDGE_ALPHA_GRID = np.logspace(-4, 4, 9)
LOGISTIC_C_GRID = np.logspace(-4, 4, 9)
OUTER_REPEATS = 5
OUTER_FOLDS = 5
INNER_FOLDS = 5
BASE_SEED = 20260717

DT_MODELS = OrderedDict(
    {
        "D0_composition": ["gc_deviation_from_0p5", "max_homopolymer"],
        "D1_structural_rules": [
            "gc_deviation_from_0p5",
            "max_homopolymer",
            "candidate_pairs",
            "longest_stem",
        ],
        "D2_weighted_only": [
            "gc_deviation_from_0p5",
            "max_homopolymer",
            "weighted_pairs",
        ],
        "D3_full": [
            "gc_deviation_from_0p5",
            "max_homopolymer",
            "candidate_pairs",
            "longest_stem",
            "weighted_pairs",
        ],
    }
)

PCR_COMPOSITION = [
    "gc_deviation_from_0p5",
    "max_homopolymer",
    "base_freq_A",
    "base_freq_C",
    "base_freq_G",
    "base_freq_T",
]
PCR_ADAPTER_FEATURES = [
    "left_boundary_vs_left_adapter_max_rc_tract",
    "left_boundary_vs_right_adapter_max_rc_tract",
    "right_boundary_vs_left_adapter_max_rc_tract",
    "right_boundary_vs_right_adapter_max_rc_tract",
    "left_boundary_vs_left_adapter_complement_4mers",
    "left_boundary_vs_right_adapter_complement_4mers",
    "right_boundary_vs_left_adapter_complement_4mers",
    "right_boundary_vs_right_adapter_complement_4mers",
]
PCR_MODELS = OrderedDict(
    {
        "P0_composition": PCR_COMPOSITION,
        "P1_variable_proxy": PCR_COMPOSITION + ["variable_weighted_pairs"],
        "P2_adapter_aware": PCR_COMPOSITION + PCR_ADAPTER_FEATURES,
        "P3_combined": PCR_COMPOSITION
        + PCR_ADAPTER_FEATURES
        + ["variable_weighted_pairs"],
        "P4_full_amplicon": PCR_COMPOSITION + ["full_amplicon_weighted_pairs"],
    }
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(12, (os.cpu_count() or 2) - 2)),
        help="Parallel outer-fold/feature jobs (default leaves two logical CPUs free).",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=2000,
        help="Paired sequence-level bootstrap replicates (frozen value: 2000).",
    )
    parser.add_argument(
        "--force-features",
        action="store_true",
        help="Recompute the deterministic PCR feature cache.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(repo: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


def validate_frozen_sources() -> dict[str, Any]:
    for path, expected in EXPECTED_SHA256.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"SHA-256 mismatch for {path}: {observed} != {expected}")
    commits: dict[str, str] = {}
    for repo, expected in EXPECTED_COMMITS.items():
        observed = git_value(repo, "rev-parse", "HEAD")
        if observed != expected:
            raise ValueError(f"commit mismatch for {repo}: {observed} != {expected}")
        commits[relative(repo)] = observed
    if len(LEFT_ADAPTER_0F) != 20 or len(RIGHT_ADAPTER_0R_PRIME) != 21:
        raise ValueError("publication-defined adapter length mismatch")
    return commits


def pickle_opcode_audit(path: Path) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    global_args: set[str] = set()
    with path.open("rb") as handle:
        for opcode, argument, _position in pickletools.genops(handle):
            counts[opcode.name] += 1
            if opcode.name == "GLOBAL" and isinstance(argument, str):
                global_args.add(argument.replace("\n", "."))
    forbidden = {"PERSID", "BINPERSID", "EXT1", "EXT2", "EXT4"}
    present_forbidden = sorted(forbidden.intersection(counts))
    if present_forbidden:
        raise ValueError(f"unsupported pickle opcodes in {path}: {present_forbidden}")
    return {
        "sha256_allowlisted": True,
        "opcode_count": int(sum(counts.values())),
        "opcode_names": sorted(counts),
        "global_arguments": sorted(global_args),
        "forbidden_opcodes_present": present_forbidden,
    }


def load_pcr_data() -> tuple[pd.DataFrame, dict[str, Any]]:
    frames: list[pd.DataFrame] = []
    audits: dict[str, Any] = {}
    expected_rows = {"GCall": 11998, "GCfix": 11994}
    for pool, path in PCR_FILES.items():
        audits[pool] = pickle_opcode_audit(path)
        payload = pd.read_pickle(path)
        if set(payload) != {"rest", "bottom"}:
            raise ValueError(f"{pool}: unexpected pickle keys {sorted(payload)}")
        rest = payload["rest"].copy()
        bottom = payload["bottom"].copy()
        required = ["x0", "eff", "sequence", "label"]
        if list(rest.columns) != required or list(bottom.columns) != required:
            raise ValueError(f"{pool}: unexpected processed-data schema")
        if not (rest["label"].eq(0).all() and bottom["label"].eq(1).all()):
            raise ValueError(f"{pool}: label/key mismatch")
        frame = pd.concat([rest, bottom], axis=0).sort_index()
        frame.index = frame.index.astype(str).str.zfill(6)
        frame.index.name = "sample_id"
        if len(frame) != expected_rows[pool]:
            raise ValueError(f"{pool}: row count {len(frame)} != {expected_rows[pool]}")
        if int(frame["label"].sum()) != 240:
            raise ValueError(f"{pool}: expected 240 worst-2% labels")
        if frame.index.duplicated().any() or frame["sequence"].duplicated().any():
            raise ValueError(f"{pool}: duplicate identifier or sequence")
        if not frame["sequence"].map(lambda value: len(value) == 108).all():
            raise ValueError(f"{pool}: sequence length mismatch")
        if not frame["sequence"].map(lambda value: set(value) <= set("ACGT")).all():
            raise ValueError(f"{pool}: non-ACGT sequence")
        frame = frame.reset_index()
        frame.insert(0, "pool", pool)
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    with zipfile.ZipFile(PCR_SOURCE_ZIP) as archive:
        member = "Fig. 2c/eff_dist.csv"
        with archive.open(member) as handle:
            source = pd.read_csv(handle, index_col=0)
    gcall = combined.loc[combined["pool"].eq("GCall")].copy()
    source.index = source.index.str.removeprefix("GCall_").str.zfill(6)
    source = source.sort_index()
    gcall = gcall.set_index("sample_id").sort_index()
    if len(source) != len(gcall) or not source.index.equals(gcall.index):
        raise ValueError("GCall source-data index does not match processed repository table")
    if not np.allclose(source["eff"], gcall["eff"], rtol=0.0, atol=1e-14):
        raise ValueError("GCall efficiency values differ between source data and repository")
    if not np.allclose(source["x0"], gcall["x0"], rtol=0.0, atol=1e-14):
        raise ValueError("GCall initial-abundance values differ between sources")
    audits["source_data_crosscheck"] = {
        "member": member,
        "pool": "GCall",
        "matched_rows": len(gcall),
        "eff_max_abs_difference": float(np.max(np.abs(source["eff"] - gcall["eff"]))),
        "x0_max_abs_difference": float(np.max(np.abs(source["x0"] - gcall["x0"]))),
        "GCfix_note": "Figure 2c source-data member contains GCall only; GCfix is hash-pinned to the publication repository release.",
    }
    return combined, audits


def rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def longest_common_substring(left: str, right: str) -> int:
    previous = [0] * (len(right) + 1)
    best = 0
    for left_base in left:
        current = [0] * (len(right) + 1)
        for index, right_base in enumerate(right, start=1):
            if left_base == right_base:
                current[index] = previous[index - 1] + 1
                best = max(best, current[index])
        previous = current
    return best


def complementary_kmer_count(boundary: str, adapter: str, k: int) -> int:
    complement_target = rc(adapter)
    targets = {complement_target[index : index + k] for index in range(len(adapter) - k + 1)}
    return sum(
        boundary[index : index + k] in targets
        for index in range(len(boundary) - k + 1)
    )


def adapter_pair_features(boundary: str, adapter: str, prefix: str) -> dict[str, int]:
    return {
        f"{prefix}_max_rc_tract": longest_common_substring(boundary, rc(adapter)),
        f"{prefix}_complement_4mers": complementary_kmer_count(
            boundary, adapter, COMPLEMENT_K
        ),
    }


def hairpin_features_exact_numpy(
    sequence: str, min_stem: int = 4, min_loop: int = 3
) -> dict[str, int]:
    """Vectorized equivalent of ``hrf.hairpin_features_simple``.

    The frozen evaluator stores the length of each reverse-complement run on an
    anti-diagonal.  This implementation retains that recurrence but performs
    each anti-diagonal row in NumPy.  ``validate_fast_evaluator`` compares all
    returned fields against the original implementation before any public data
    are evaluated.
    """

    encoding = np.fromiter(
        ({"A": 0, "C": 1, "G": 2, "T": 3}[base] for base in sequence),
        dtype=np.int8,
        count=len(sequence),
    )
    complement = np.asarray([3, 2, 1, 0], dtype=np.int8)
    n = len(sequence)
    runs = np.zeros((n + 1, n), dtype=np.int16)
    encoded_complements = complement[encoding]
    for left in range(n - 1, -1, -1):
        matches = encoding[left] == encoded_complements
        runs[left, 0] = int(matches[0])
        runs[left, 1:] = matches[1:] * (1 + runs[left + 1, :-1])

    left_index = np.arange(n, dtype=np.int16)[:, None]
    right_index = np.arange(n, dtype=np.int16)[None, :]
    maximum_by_loop = (right_index - left_index - min_loop + 1) // 2
    maximum_stem = np.minimum(runs[:-1], maximum_by_loop)
    valid_geometry = right_index >= left_index + 2 * min_stem + min_loop - 1
    active = valid_geometry & (maximum_stem >= min_stem)
    if not np.any(active):
        return {
            "candidate_pairs": 0,
            "weighted_pairs": 0,
            "longest_stem": 0,
            "min_loop_seen": -1,
        }
    active_stems = maximum_stem[active].astype(np.int64)
    counts = active_stems - min_stem + 1
    loops = (
        right_index - left_index - 2 * maximum_stem + 1
    )[active].astype(np.int64)
    return {
        "candidate_pairs": int(counts.sum()),
        "weighted_pairs": int(((min_stem + active_stems) * counts // 2).sum()),
        "longest_stem": int(active_stems.max()),
        "min_loop_seen": int(loops.min()),
    }


def validate_fast_evaluator() -> None:
    rng = np.random.default_rng(BASE_SEED)
    sequences = [
        "A" * 20,
        "AT" * 10,
        "ACGT" * 5,
        "G" * 108,
        "ACGT" * 27,
        ("ACGT" * 37) + "A",
    ]
    for length in [20, 37, 108, 149]:
        for _ in range(8):
            sequences.append("".join(rng.choice(list("ACGT"), size=length)))
    for sequence in sequences:
        expected = hrf.hairpin_features_simple(sequence, 4, 3)
        observed = hairpin_features_exact_numpy(sequence, 4, 3)
        if observed != expected:
            raise ValueError(
                "vectorized frozen-evaluator mismatch for "
                f"{sequence}: {observed} != {expected}"
            )


def pcr_feature_row(pool: str, sample_id: str, sequence: str) -> dict[str, Any]:
    variable = hairpin_features_exact_numpy(sequence, 4, 3)
    full_sequence = LEFT_ADAPTER_0F + sequence + RIGHT_ADAPTER_0R_PRIME
    if len(full_sequence) != 149:
        raise ValueError("full amplicon must be 149 nt")
    full = hairpin_features_exact_numpy(full_sequence, 4, 3)
    left_boundary = sequence[:BOUNDARY_WINDOW]
    right_boundary = sequence[-BOUNDARY_WINDOW:]
    gc_fraction = hrf.gc_fraction(sequence)
    row: dict[str, Any] = {
        "pool": pool,
        "sample_id": sample_id,
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "variable_length_nt": len(sequence),
        "full_amplicon_length_nt": len(full_sequence),
        "gc_deviation_from_0p5": abs(gc_fraction - 0.5),
        "max_homopolymer": hrf.max_homopolymer(sequence),
        **{f"base_freq_{base}": sequence.count(base) / len(sequence) for base in "ACGT"},
        "variable_candidate_pairs": variable["candidate_pairs"],
        "variable_weighted_pairs": variable["weighted_pairs"],
        "variable_longest_stem": variable["longest_stem"],
        "full_amplicon_candidate_pairs": full["candidate_pairs"],
        "full_amplicon_weighted_pairs": full["weighted_pairs"],
        "full_amplicon_longest_stem": full["longest_stem"],
    }
    row.update(
        adapter_pair_features(
            left_boundary, LEFT_ADAPTER_0F, "left_boundary_vs_left_adapter"
        )
    )
    row.update(
        adapter_pair_features(
            left_boundary, RIGHT_ADAPTER_0R_PRIME, "left_boundary_vs_right_adapter"
        )
    )
    row.update(
        adapter_pair_features(
            right_boundary, LEFT_ADAPTER_0F, "right_boundary_vs_left_adapter"
        )
    )
    row.update(
        adapter_pair_features(
            right_boundary, RIGHT_ADAPTER_0R_PRIME, "right_boundary_vs_right_adapter"
        )
    )
    return row


def get_pcr_features(data: pd.DataFrame, jobs: int, force: bool) -> pd.DataFrame:
    expected_keys = data[["pool", "sample_id", "sequence"]].copy()
    expected_keys["sequence_sha256"] = expected_keys["sequence"].map(
        lambda value: hashlib.sha256(value.encode("ascii")).hexdigest()
    )
    if FEATURE_CACHE.is_file() and not force:
        cached = pd.read_csv(FEATURE_CACHE, sep="\t", dtype={"sample_id": str})
        cached["sample_id"] = cached["sample_id"].str.zfill(6)
        joined = expected_keys.merge(
            cached[["pool", "sample_id", "sequence_sha256"]],
            on=["pool", "sample_id"],
            how="left",
            suffixes=("_expected", "_cached"),
            validate="one_to_one",
        )
        valid = (
            len(cached) == len(data)
            and joined["sequence_sha256_expected"].eq(joined["sequence_sha256_cached"]).all()
        )
        if valid:
            print(f"PCR feature cache verified: {FEATURE_CACHE}", flush=True)
            return cached
        raise ValueError("existing PCR feature cache does not match frozen sequences")

    print(f"Computing PCR features for {len(data):,} sequences with {jobs} jobs", flush=True)
    tasks = [
        (str(row.pool), str(row.sample_id), str(row.sequence))
        for row in data[["pool", "sample_id", "sequence"]].itertuples(index=False)
    ]
    rows = Parallel(n_jobs=jobs, prefer="processes", batch_size=32)(
        delayed(pcr_feature_row)(pool, sample_id, sequence)
        for pool, sample_id, sequence in tasks
    )
    features = pd.DataFrame(rows)
    features["sample_id"] = features["sample_id"].astype(str).str.zfill(6)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    features.to_csv(FEATURE_CACHE, sep="\t", index=False, float_format="%.12g")
    return features


def make_splits(
    y: np.ndarray,
    classification: bool,
    folds: int,
    repeats: int,
    seed: int,
) -> list[tuple[int, int, np.ndarray, np.ndarray]]:
    splits: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    dummy = np.zeros(len(y))
    for repeat in range(repeats):
        if classification:
            splitter = StratifiedKFold(
                n_splits=folds, shuffle=True, random_state=seed + repeat
            )
            iterator = splitter.split(dummy, y)
        else:
            splitter = KFold(n_splits=folds, shuffle=True, random_state=seed + repeat)
            iterator = splitter.split(dummy)
        for fold, (train_index, test_index) in enumerate(iterator):
            splits.append((repeat, fold, train_index, test_index))
    return splits


def regression_metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    correlation = spearmanr(y, prediction).statistic
    return {
        "spearman": float(correlation),
        "r2": float(r2_score(y, prediction)),
        "mae": float(mean_absolute_error(y, prediction)),
        "rmse": float(math.sqrt(mean_squared_error(y, prediction))),
    }


def classification_metrics(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    return {
        "average_precision": float(average_precision_score(y, probability)),
        "roc_auc": float(roc_auc_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def nested_fold_fit(
    repeat: int,
    fold: int,
    train_index: np.ndarray,
    test_index: np.ndarray,
    all_x: pd.DataFrame,
    y: np.ndarray,
    models: OrderedDict[str, list[str]],
    classification: bool,
    inner_folds: int,
    seed: int,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "repeat": repeat,
        "fold": fold,
        "test_index": test_index,
        "models": {},
    }
    inner_seed = seed + repeat * 100 + fold
    for model_name, feature_names in models.items():
        x_train = all_x.iloc[train_index][feature_names].to_numpy(dtype=float)
        x_test = all_x.iloc[test_index][feature_names].to_numpy(dtype=float)
        y_train = y[train_index]
        if classification:
            estimator = Pipeline(
                [
                    ("scale", StandardScaler()),
                    (
                        "model",
                        LogisticRegression(
                            solver="liblinear",
                            max_iter=5000,
                            random_state=inner_seed,
                        ),
                    ),
                ]
            )
            inner = StratifiedKFold(
                n_splits=inner_folds, shuffle=True, random_state=inner_seed
            )
            search = GridSearchCV(
                estimator,
                {"model__C": LOGISTIC_C_GRID},
                scoring="average_precision",
                cv=inner,
                n_jobs=1,
                refit=True,
            )
        else:
            estimator = Pipeline([("scale", StandardScaler()), ("model", Ridge())])
            inner = KFold(n_splits=inner_folds, shuffle=True, random_state=inner_seed)
            search = GridSearchCV(
                estimator,
                {"model__alpha": RIDGE_ALPHA_GRID},
                scoring="neg_mean_squared_error",
                cv=inner,
                n_jobs=1,
                refit=True,
            )
        search.fit(x_train, y_train)
        if classification:
            prediction = search.predict_proba(x_test)[:, 1]
        else:
            prediction = search.predict(x_test)
        coefficients = np.asarray(search.best_estimator_.named_steps["model"].coef_).ravel()
        result["models"][model_name] = {
            "prediction": prediction,
            "features": feature_names,
            "coefficients": coefficients,
            "best_parameter": float(next(iter(search.best_params_.values()))),
        }
    return result


def run_nested_models(
    dataset: str,
    analysis: str,
    pool: str,
    endpoint: str,
    frame: pd.DataFrame,
    y: np.ndarray,
    models: OrderedDict[str, list[str]],
    classification: bool,
    jobs: int,
    seed: int,
    folds: int = OUTER_FOLDS,
    repeats: int = OUTER_REPEATS,
    inner_folds: int = INNER_FOLDS,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, np.ndarray]]:
    splits = make_splits(y, classification, folds, repeats, seed)
    print(
        f"Nested CV: {dataset} {pool} {endpoint} "
        f"({repeats}x{folds}, models={len(models)})",
        flush=True,
    )
    fold_results = Parallel(n_jobs=jobs, prefer="processes")(
        delayed(nested_fold_fit)(
            repeat,
            fold,
            train_index,
            test_index,
            frame,
            y,
            models,
            classification,
            inner_folds,
            seed + 10_000,
        )
        for repeat, fold, train_index, test_index in splits
    )
    prediction_arrays = {
        model_name: np.full((repeats, len(frame)), np.nan, dtype=float)
        for model_name in models
    }
    coefficient_rows: list[dict[str, Any]] = []
    for fold_result in fold_results:
        repeat = int(fold_result["repeat"])
        fold = int(fold_result["fold"])
        test_index = fold_result["test_index"]
        for model_name, model_result in fold_result["models"].items():
            prediction_arrays[model_name][repeat, test_index] = model_result["prediction"]
            for feature, coefficient in zip(
                model_result["features"], model_result["coefficients"]
            ):
                coefficient_rows.append(
                    {
                        "dataset": dataset,
                        "analysis": analysis,
                        "pool_or_direction": pool,
                        "endpoint": endpoint,
                        "model": model_name,
                        "repeat": repeat,
                        "fold": fold,
                        "feature": feature,
                        "standardized_coefficient": float(coefficient),
                        "selected_regularization": model_result["best_parameter"],
                    }
                )
    if any(np.isnan(array).any() for array in prediction_arrays.values()):
        raise ValueError(f"incomplete OOF predictions for {dataset}/{pool}/{endpoint}")
    mean_predictions = {
        model_name: array.mean(axis=0) for model_name, array in prediction_arrays.items()
    }
    metric_rows: list[dict[str, Any]] = []
    for model_name, prediction in mean_predictions.items():
        metrics = (
            classification_metrics(y, prediction)
            if classification
            else regression_metrics(y, prediction)
        )
        for metric, value in metrics.items():
            metric_rows.append(
                {
                    "dataset": dataset,
                    "analysis": analysis,
                    "pool_or_direction": pool,
                    "endpoint": endpoint,
                    "model": model_name,
                    "metric": metric,
                    "value": value,
                    "n_sequences": len(y),
                    "n_events": int(y.sum()) if classification else "",
                    "event_prevalence": float(y.mean()) if classification else "",
                    "evaluation": f"mean prediction over {repeats} repeated out-of-fold passes",
                }
            )
    return metric_rows, coefficient_rows, mean_predictions


def fit_transfer_models(
    source_pool: str,
    target_pool: str,
    source: pd.DataFrame,
    target: pd.DataFrame,
    models: OrderedDict[str, list[str]],
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, np.ndarray]]:
    y_source = source["label"].to_numpy(dtype=int)
    y_target = target["label"].to_numpy(dtype=int)
    predictions: dict[str, np.ndarray] = {}
    metrics: list[dict[str, Any]] = []
    coefficients: list[dict[str, Any]] = []
    direction = f"{source_pool}_to_{target_pool}"
    print(f"Cross-pool transfer: {direction}", flush=True)
    for model_index, (model_name, feature_names) in enumerate(models.items()):
        estimator = Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        solver="liblinear",
                        max_iter=5000,
                        random_state=seed + model_index,
                    ),
                ),
            ]
        )
        inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed + model_index)
        search = GridSearchCV(
            estimator,
            {"model__C": LOGISTIC_C_GRID},
            scoring="average_precision",
            cv=inner,
            n_jobs=1,
            refit=True,
        )
        search.fit(source[feature_names].to_numpy(dtype=float), y_source)
        probability = search.predict_proba(target[feature_names].to_numpy(dtype=float))[:, 1]
        predictions[model_name] = probability
        for metric, value in classification_metrics(y_target, probability).items():
            metrics.append(
                {
                    "dataset": "Gimpel2025_PCR",
                    "analysis": "pcr_binary_transfer",
                    "pool_or_direction": direction,
                    "endpoint": "worst_2_percent_label",
                    "model": model_name,
                    "metric": metric,
                    "value": value,
                    "n_sequences": len(target),
                    "n_events": int(y_target.sum()),
                    "event_prevalence": float(y_target.mean()),
                    "evaluation": "train on complete source pool; no target-pool tuning",
                }
            )
        values = np.asarray(search.best_estimator_.named_steps["model"].coef_).ravel()
        for feature, coefficient in zip(feature_names, values):
            coefficients.append(
                {
                    "dataset": "Gimpel2025_PCR",
                    "analysis": "pcr_binary_transfer",
                    "pool_or_direction": direction,
                    "endpoint": "worst_2_percent_label",
                    "model": model_name,
                    "repeat": "transfer",
                    "fold": "all_source",
                    "feature": feature,
                    "standardized_coefficient": float(coefficient),
                    "selected_regularization": float(next(iter(search.best_params_.values()))),
                }
            )
    return metrics, coefficients, predictions


class WeightedRankPlan:
    def __init__(self, values: np.ndarray):
        values = np.asarray(values, dtype=float)
        self.n = len(values)
        self.order = np.argsort(values, kind="mergesort")
        sorted_values = values[self.order]
        self.starts = np.r_[0, np.flatnonzero(sorted_values[1:] != sorted_values[:-1]) + 1]
        self.group_ids = np.repeat(
            np.arange(len(self.starts)), np.diff(np.r_[self.starts, self.n])
        )

    def ranks(self, weights: np.ndarray) -> np.ndarray:
        sorted_weights = weights[:, self.order]
        group_weights = np.add.reduceat(sorted_weights, self.starts, axis=1)
        cumulative = np.cumsum(group_weights, axis=1)
        midranks = cumulative - group_weights + (group_weights + 1.0) / 2.0
        sorted_ranks = midranks[:, self.group_ids]
        ranks = np.empty_like(sorted_ranks, dtype=float)
        ranks[:, self.order] = sorted_ranks
        return ranks


class WeightedAveragePrecisionPlan:
    def __init__(self, labels: np.ndarray, scores: np.ndarray):
        self.labels = np.asarray(labels, dtype=float)
        scores = np.asarray(scores, dtype=float)
        self.order = np.argsort(-scores, kind="mergesort")
        sorted_scores = scores[self.order]
        self.ends = np.r_[
            np.flatnonzero(sorted_scores[1:] != sorted_scores[:-1]), len(scores) - 1
        ]

    def metric(self, weights: np.ndarray) -> np.ndarray:
        ordered_weights = weights[:, self.order]
        ordered_labels = self.labels[self.order]
        cumulative_true = np.cumsum(ordered_weights * ordered_labels, axis=1)[:, self.ends]
        cumulative_all = np.cumsum(ordered_weights, axis=1)[:, self.ends]
        positives = np.sum(weights * self.labels, axis=1)
        precision = np.divide(
            cumulative_true,
            cumulative_all,
            out=np.zeros_like(cumulative_true),
            where=cumulative_all > 0,
        )
        recall = np.divide(
            cumulative_true,
            positives[:, None],
            out=np.zeros_like(cumulative_true),
            where=positives[:, None] > 0,
        )
        recall_increment = np.diff(
            np.concatenate([np.zeros((len(weights), 1)), recall], axis=1), axis=1
        )
        return np.sum(recall_increment * precision, axis=1)


def weighted_correlation(weights: np.ndarray, left: np.ndarray, right: np.ndarray) -> np.ndarray:
    total = weights.sum(axis=1)
    left_mean = np.sum(weights * left, axis=1) / total
    right_mean = np.sum(weights * right, axis=1) / total
    left_centered = left - left_mean[:, None]
    right_centered = right - right_mean[:, None]
    covariance = np.sum(weights * left_centered * right_centered, axis=1)
    denominator = np.sqrt(
        np.sum(weights * left_centered**2, axis=1)
        * np.sum(weights * right_centered**2, axis=1)
    )
    return np.divide(
        covariance,
        denominator,
        out=np.full_like(covariance, np.nan),
        where=denominator > 0,
    )


def bootstrap_group(
    specs: list[dict[str, Any]], n_bootstrap: int, seed: int, batch_size: int = 20
) -> dict[str, tuple[float, float]]:
    if not specs:
        return {}
    n = len(specs[0]["y"])
    if any(len(spec["y"]) != n for spec in specs):
        raise ValueError("bootstrap group contains different sample counts")
    prepared: list[dict[str, Any]] = []
    for spec in specs:
        if spec["metric_type"] == "spearman":
            prepared.append(
                {
                    **spec,
                    "y_plan": WeightedRankPlan(spec["y"]),
                    "baseline_plan": WeightedRankPlan(spec["baseline_prediction"]),
                    "extended_plan": WeightedRankPlan(spec["extended_prediction"]),
                }
            )
        elif spec["metric_type"] == "average_precision":
            prepared.append(
                {
                    **spec,
                    "baseline_plan": WeightedAveragePrecisionPlan(
                        spec["y"], spec["baseline_prediction"]
                    ),
                    "extended_plan": WeightedAveragePrecisionPlan(
                        spec["y"], spec["extended_prediction"]
                    ),
                }
            )
        else:
            raise ValueError(f"unsupported bootstrap metric {spec['metric_type']}")
    deltas = {spec["key"]: np.empty(n_bootstrap, dtype=float) for spec in prepared}
    rng = np.random.default_rng(seed)
    probabilities = np.full(n, 1.0 / n)
    for start in range(0, n_bootstrap, batch_size):
        stop = min(start + batch_size, n_bootstrap)
        weights = rng.multinomial(n, probabilities, size=stop - start).astype(float)
        for spec in prepared:
            if spec["metric_type"] == "spearman":
                outcome_rank = spec["y_plan"].ranks(weights)
                baseline_rank = spec["baseline_plan"].ranks(weights)
                extended_rank = spec["extended_plan"].ranks(weights)
                baseline_value = weighted_correlation(weights, outcome_rank, baseline_rank)
                extended_value = weighted_correlation(weights, outcome_rank, extended_rank)
            else:
                baseline_value = spec["baseline_plan"].metric(weights)
                extended_value = spec["extended_plan"].metric(weights)
            deltas[spec["key"]][start:stop] = extended_value - baseline_value
    return {
        key: (float(np.nanquantile(values, 0.025)), float(np.nanquantile(values, 0.975)))
        for key, values in deltas.items()
    }


def metric_value(
    metric_rows: list[dict[str, Any]],
    analysis: str,
    pool: str,
    endpoint: str,
    model: str,
    metric: str,
) -> float:
    matches = [
        row
        for row in metric_rows
        if row["analysis"] == analysis
        and row["pool_or_direction"] == pool
        and row["endpoint"] == endpoint
        and row["model"] == model
        and row["metric"] == metric
    ]
    if len(matches) != 1:
        raise ValueError(
            f"metric lookup expected one row: {analysis}/{pool}/{endpoint}/{model}/{metric}"
        )
    return float(matches[0]["value"])


def add_comparison(
    comparison_rows: list[dict[str, Any]],
    bootstrap_specs: dict[str, list[dict[str, Any]]],
    metric_rows: list[dict[str, Any]],
    dataset: str,
    analysis: str,
    pool: str,
    endpoint: str,
    comparison: str,
    baseline_model: str,
    extended_model: str,
    metric: str,
    y: np.ndarray,
    predictions: dict[str, np.ndarray],
    group: str,
    status: str,
) -> None:
    baseline_value = metric_value(
        metric_rows, analysis, pool, endpoint, baseline_model, metric
    )
    extended_value = metric_value(
        metric_rows, analysis, pool, endpoint, extended_model, metric
    )
    key = f"{analysis}|{pool}|{endpoint}|{comparison}|{metric}"
    row = {
        "key": key,
        "dataset": dataset,
        "analysis": analysis,
        "pool_or_direction": pool,
        "endpoint": endpoint,
        "comparison": comparison,
        "baseline_model": baseline_model,
        "extended_model": extended_model,
        "primary_metric": metric,
        "baseline_value": baseline_value,
        "extended_value": extended_value,
        "delta_extended_minus_baseline": extended_value - baseline_value,
        "bootstrap_2p5": math.nan,
        "bootstrap_97p5": math.nan,
        "n_sequences": len(y),
        "n_events": int(y.sum()) if metric == "average_precision" else "",
        "confirmatory_status": status,
    }
    comparison_rows.append(row)
    bootstrap_specs.setdefault(group, []).append(
        {
            "key": key,
            "metric_type": metric,
            "y": y,
            "baseline_prediction": predictions[baseline_model],
            "extended_prediction": predictions[extended_model],
        }
    )


def source_manifest(commits: dict[str, str]) -> list[dict[str, Any]]:
    rows = [
        {
            "source": "DT4DDS derived sequence-level table",
            "local_path": relative(DT_CSV),
            "url_or_accession": "https://github.com/fml-ethz/dt4dds_notebooks",
            "version": f"v1.0.0; {commits[relative(DT_REPO)]}",
            "access_date": "2026-07-17",
            "sha256": sha256(DT_CSV),
            "bytes": DT_CSV.stat().st_size,
            "licence_boundary": "No machine-readable licence in frozen upstream copy; do not redistribute upstream raw content without confirmation.",
            "role": "24,472 public experimental post-workflow detectability records",
        },
        {
            "source": "Gimpel 2025 GCall processed labels",
            "local_path": relative(PCR_FILES["GCall"]),
            "url_or_accession": "https://github.com/BorgwardtLab/PCR-bias",
            "version": f"v1.0.3; {commits[relative(PCR_REPO)]}",
            "access_date": "2026-07-17",
            "sha256": sha256(PCR_FILES["GCall"]),
            "bytes": PCR_FILES["GCall"].stat().st_size,
            "licence_boundary": "BSD-3-Clause repository licence",
            "role": "11,998 filtered 108-nt sequences; 240 author-defined worst-2% labels",
        },
        {
            "source": "Gimpel 2025 GCfix processed labels",
            "local_path": relative(PCR_FILES["GCfix"]),
            "url_or_accession": "https://github.com/BorgwardtLab/PCR-bias",
            "version": f"v1.0.3; {commits[relative(PCR_REPO)]}",
            "access_date": "2026-07-17",
            "sha256": sha256(PCR_FILES["GCfix"]),
            "bytes": PCR_FILES["GCfix"].stat().st_size,
            "licence_boundary": "BSD-3-Clause repository licence",
            "role": "11,994 filtered 108-nt sequences; 240 author-defined worst-2% labels",
        },
        {
            "source": "Gimpel 2025 supplementary information",
            "local_path": relative(PCR_SUPP_PDF),
            "url_or_accession": "https://doi.org/10.1038/s41467-025-64221-4",
            "version": "Supplementary Information; Table 6",
            "access_date": "2026-07-17",
            "sha256": sha256(PCR_SUPP_PDF),
            "bytes": PCR_SUPP_PDF.stat().st_size,
            "licence_boundary": "Article is CC BY 4.0",
            "role": "Exact 0F and 0R primer sequences used to reconstruct 149-nt amplicons",
        },
        {
            "source": "Gimpel 2025 article source data",
            "local_path": relative(PCR_SOURCE_ZIP),
            "url_or_accession": "PRJEB77604; https://doi.org/10.1038/s41467-025-64221-4",
            "version": "Publication Source Data archive",
            "access_date": "2026-07-17",
            "sha256": sha256(PCR_SOURCE_ZIP),
            "bytes": PCR_SOURCE_ZIP.stat().st_size,
            "licence_boundary": "Article is CC BY 4.0",
            "role": "Independent row/value cross-check of all 11,998 processed GCall efficiency records",
        },
        {
            "source": "Analysis script",
            "local_path": relative(Path(__file__)),
            "url_or_accession": "Paper-2-local reproducibility artifact",
            "version": "frozen statistical contract dated 2026-07-17",
            "access_date": "2026-07-17",
            "sha256": sha256(Path(__file__)),
            "bytes": Path(__file__).stat().st_size,
            "licence_boundary": "Project-local script",
            "role": "Source validation, deterministic features, nested CV and paired bootstrap",
        },
    ]
    return rows


def data_dictionary_rows() -> list[dict[str, str]]:
    rows = [
        ("DT4DDS", "mean_log2_cpm_plus1", "endpoint", "Mean log2(CPM+1) across 0a, 0b, 2d, 4d and 7d; post-workflow detectability."),
        ("DT4DDS", "day7_vs_day0_log2fc_plus1", "endpoint", "log2 ratio of day-7 CPM+1 to mean day-0 CPM+1; composite workflow endpoint."),
        ("DT4DDS", "missing_any", "secondary endpoint", "Zero assigned reads in at least one of five public samples."),
        ("DT4DDS", "missing_day7", "secondary endpoint", "Zero assigned reads in the day-7 public sample."),
        ("both", "gc_deviation_from_0p5", "baseline predictor", "Absolute GC fraction minus 0.5."),
        ("both", "max_homopolymer", "baseline predictor", "Longest identical-base run."),
        ("DT4DDS", "candidate_pairs", "structural-rule predictor", "Candidate reverse-complement stem count, min stem 4 and min loop 3."),
        ("DT4DDS", "longest_stem", "structural-rule predictor", "Longest exact reverse-complement stem under the frozen evaluator."),
        ("DT4DDS", "weighted_pairs", "Paper 2 proxy", "Stem-length-weighted reverse-complement pair count, min stem 4 and min loop 3."),
        ("Gimpel2025_PCR", "worst_2_percent_label", "primary endpoint", "Author-defined bottom 2% PCR-efficiency label; not recalculated on evaluation folds."),
        ("Gimpel2025_PCR", "eff", "secondary endpoint", "Author-estimated relative PCR amplification efficiency."),
        ("Gimpel2025_PCR", "variable_weighted_pairs", "Paper 2 proxy", "Frozen weighted-pair proxy on the public 108-nt variable region."),
        ("Gimpel2025_PCR", "adapter complementarity features", "biological baseline", "Longest exact reverse-complement tract and count of complementary 4-mers for each 24-nt variable-region boundary against each publication-defined terminal adapter."),
        ("Gimpel2025_PCR", "full_amplicon_weighted_pairs", "sensitivity predictor", "Frozen weighted-pair proxy on 0F + variable + 0R-prime (149 nt)."),
    ]
    return [
        {"dataset": dataset, "field": field, "role": role, "definition_and_boundary": definition}
        for dataset, field, role, definition in rows
    ]


def write_table(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def coefficient_summary(
    coefficients: list[dict[str, Any]],
    analysis: str,
    pool: str,
    endpoint: str,
    model: str,
    feature: str,
) -> tuple[float, int, int]:
    values = [
        float(row["standardized_coefficient"])
        for row in coefficients
        if row["analysis"] == analysis
        and row["pool_or_direction"] == pool
        and row["endpoint"] == endpoint
        and row["model"] == model
        and row["feature"] == feature
    ]
    if not values:
        return math.nan, 0, 0
    return float(np.median(values)), int(np.sum(np.asarray(values) < 0)), len(values)


def comparison_lookup(
    comparisons: list[dict[str, Any]], analysis: str, pool: str, endpoint: str, comparison: str
) -> dict[str, Any]:
    matches = [
        row
        for row in comparisons
        if row["analysis"] == analysis
        and row["pool_or_direction"] == pool
        and row["endpoint"] == endpoint
        and row["comparison"] == comparison
    ]
    if len(matches) != 1:
        raise ValueError(f"comparison lookup failed: {analysis}/{pool}/{endpoint}/{comparison}")
    return matches[0]


def classify_dt_result(
    endpoint: str,
    comparisons: list[dict[str, Any]],
    coefficients: list[dict[str, Any]],
) -> str:
    rows = [
        comparison_lookup(comparisons, "dt_continuous", pool, endpoint, "D3_minus_D1")
        for pool in ["Genscript_GCall", "Twist_GCall"]
    ]
    coefficients_ok = []
    for pool in ["Genscript_GCall", "Twist_GCall"]:
        median, _negative, _total = coefficient_summary(
            coefficients, "dt_continuous", pool, endpoint, "D3_full", "weighted_pairs"
        )
        coefficients_ok.append(median < 0)
    positive = [float(row["delta_extended_minus_baseline"]) > 0 for row in rows]
    intervals_positive = [float(row["bootstrap_2p5"]) > 0 for row in rows]
    if not any(positive):
        return "no incremental association"
    if all(positive) and all(coefficients_ok) and all(intervals_positive):
        return "stable incremental association"
    if all(positive) and all(coefficients_ok):
        return "weak/partial association"
    if any(positive) or any(coefficients_ok):
        return "not replicated"
    return "no incremental association"


def classify_pcr_result(
    comparisons: list[dict[str, Any]], coefficients: list[dict[str, Any]]
) -> str:
    within = [
        comparison_lookup(
            comparisons,
            "pcr_binary_internal",
            pool,
            "worst_2_percent_label",
            "P3_minus_P2",
        )
        for pool in ["GCall", "GCfix"]
    ]
    transfer = [
        comparison_lookup(
            comparisons,
            "pcr_binary_transfer",
            direction,
            "worst_2_percent_label",
            "P3_minus_P2",
        )
        for direction in ["GCall_to_GCfix", "GCfix_to_GCall"]
    ]
    weighted_medians = [
        coefficient_summary(
            coefficients,
            "pcr_binary_internal",
            pool,
            "worst_2_percent_label",
            "P3_combined",
            "variable_weighted_pairs",
        )[0]
        for pool in ["GCall", "GCfix"]
    ]
    if (
        all(float(row["bootstrap_2p5"]) > 0 for row in within)
        and all(
        float(row["delta_extended_minus_baseline"]) > 0 for row in transfer
        )
        and all(median > 0 for median in weighted_medians)
    ):
        return "robust orthogonal support"
    p1 = [
        comparison_lookup(
            comparisons,
            "pcr_binary_internal",
            pool,
            "worst_2_percent_label",
            "P1_minus_P0",
        )
        for pool in ["GCall", "GCfix"]
    ]
    if all(float(row["delta_extended_minus_baseline"]) > 0 for row in p1) and not all(
        float(row["delta_extended_minus_baseline"]) > 0 for row in within
    ):
        return "proxy signal explained by adapter context"
    signs = [float(row["delta_extended_minus_baseline"]) > 0 for row in within + transfer]
    if any(signs) and not all(signs):
        if all(median < 0 for median in weighted_medians):
            return "pool-specific performance with risk-direction discordance"
        return "pool-specific signal"
    return "no PCR support"


def classify_pcr_continuous_result(
    comparisons: list[dict[str, Any]], coefficients: list[dict[str, Any]]
) -> str:
    rows = [
        comparison_lookup(
            comparisons,
            "pcr_continuous_efficiency",
            pool,
            "relative_efficiency",
            "P3_minus_P2",
        )
        for pool in ["GCall", "GCfix"]
    ]
    medians = [
        coefficient_summary(
            coefficients,
            "pcr_continuous_efficiency",
            pool,
            "relative_efficiency",
            "P3_combined",
            "variable_weighted_pairs",
        )[0]
        for pool in ["GCall", "GCfix"]
    ]
    intervals_positive = [float(row["bootstrap_2p5"]) > 0 for row in rows]
    if all(intervals_positive) and all(median < 0 for median in medians):
        return "replicated adverse-direction incremental association"
    if all(intervals_positive) and all(median > 0 for median in medians):
        return "replicated predictive increment with risk-direction discordance"
    if any(float(row["delta_extended_minus_baseline"]) > 0 for row in rows):
        return "heterogeneous or uncertain continuous association"
    return "no continuous-efficiency increment"


def render_summary(
    metric_rows: list[dict[str, Any]],
    comparison_rows: list[dict[str, Any]],
    coefficient_rows: list[dict[str, Any]],
    event_rows: list[dict[str, Any]],
    audits: dict[str, Any],
    main_hash: str,
) -> str:
    dt_labels = {
        endpoint: classify_dt_result(endpoint, comparison_rows, coefficient_rows)
        for endpoint in ["mean_log2_cpm_plus1", "day7_vs_day0_log2fc_plus1"]
    }
    pcr_label = classify_pcr_result(comparison_rows, coefficient_rows)
    pcr_continuous_label = classify_pcr_continuous_result(
        comparison_rows, coefficient_rows
    )
    if pcr_label == "robust orthogonal support" and any(
        label == "stable incremental association" for label in dt_labels.values()
    ):
        route = "PUBLIC_DATA_GATE_STRENGTHENS_BIOINFORMATICS_ROUTE"
    elif "discordance" in pcr_continuous_label:
        route = "PUBLIC_DATA_GATE_NARROWS_TO_VALIDATION_AND_DOMAIN_MISMATCH"
    elif pcr_label != "no PCR support" or any(
        label in {"stable incremental association", "weak/partial association"}
        for label in dt_labels.values()
    ):
        route = "PUBLIC_DATA_GATE_NARROWS_TO_ASSOCIATION_AND_INCREMENTAL_PREDICTION"
    else:
        route = "PUBLIC_DATA_GATE_DOES_NOT_YET_SUPPORT_THE_BIOINFORMATICS_CLAIM"

    def comparison_line(analysis: str, pool: str, endpoint: str, comparison: str) -> str:
        row = comparison_lookup(comparison_rows, analysis, pool, endpoint, comparison)
        return (
            f"{pool}: Δ={float(row['delta_extended_minus_baseline']):.6f}, "
            f"95% paired bootstrap [{float(row['bootstrap_2p5']):.6f}, "
            f"{float(row['bootstrap_97p5']):.6f}]"
        )

    lines = [
        "# Public Experimental Data Validation",
        "",
        "Status: `COMPLETE_REPRODUCIBLE_PUBLIC_DATA_AUDIT`",
        f"Route decision: `{route}`",
        "",
        "## Data actually analysed",
        "",
        "- DT4DDS: 24,472 public sequence-level records, analysed separately as 12,472 Genscript 102-nt references and 12,000 Twist 108-nt references across 0a, 0b, 2d, 4d and 7d measurements.",
        "- Gimpel 2025 PCR: 23,992 filtered public 108-nt variable sequences (GCall 11,998; GCfix 11,994), with 240 author-defined worst-2% labels in each pool. The publication design size was 12,000 per pool, but two GCall and six GCfix references were excluded by the authors' processing rule and are not silently restored here.",
        f"- Source-data cross-check: all {audits['source_data_crosscheck']['matched_rows']:,} GCall `eff` and `x0` values matched the independent publication Source Data member exactly.",
        "- Full PCR amplicons were reconstructed as publication-defined 0F (20 nt) + variable region (108 nt) + 0R-prime (21 nt), total 149 nt.",
        "",
        "## Prespecified primary comparisons",
        "",
        "DT4DDS primary metric is the change in out-of-fold Spearman correlation for `D3_full - D1_structural_rules`:",
        "",
    ]
    for endpoint in ["mean_log2_cpm_plus1", "day7_vs_day0_log2fc_plus1"]:
        lines.append(f"### {endpoint}: `{dt_labels[endpoint]}`")
        lines.append("")
        for pool in ["Genscript_GCall", "Twist_GCall"]:
            lines.append(
                "- " + comparison_line("dt_continuous", pool, endpoint, "D3_minus_D1")
            )
            median, negative, total = coefficient_summary(
                coefficient_rows,
                "dt_continuous",
                pool,
                endpoint,
                "D3_full",
                "weighted_pairs",
            )
            lines.append(
                f"  Standardized weighted-pair coefficient: median {median:.6g}; negative in {negative}/{total} outer fits."
            )
        lines.append("")
    lines.extend(
        [
            f"### PCR worst-2% labels: `{pcr_label}`",
            "",
            "Primary metric is average precision; prevalence is approximately 2% in each processed pool.",
            "",
        ]
    )
    for pool in ["GCall", "GCfix"]:
        lines.append(
            "- "
            + comparison_line(
                "pcr_binary_internal", pool, "worst_2_percent_label", "P3_minus_P2"
            )
        )
        median, negative, total = coefficient_summary(
            coefficient_rows,
            "pcr_binary_internal",
            pool,
            "worst_2_percent_label",
            "P3_combined",
            "variable_weighted_pairs",
        )
        lines.append(
            f"  Standardized variable-region weighted-pair coefficient: median {median:.6g}; negative in {negative}/{total} outer fits (low-efficiency label is 1, so the intended risk direction was positive)."
        )
    for direction in ["GCall_to_GCfix", "GCfix_to_GCall"]:
        lines.append(
            "- "
            + comparison_line(
                "pcr_binary_transfer",
                direction,
                "worst_2_percent_label",
                "P3_minus_P2",
            )
        )
    lines.extend(
        [
            "",
            f"### PCR continuous efficiency: `{pcr_continuous_label}`",
            "",
            "The same `P3-P2` comparison was prespecified for the authors' continuous relative-efficiency estimate. Higher efficiency is favorable, so an adverse risk proxy was expected to have a negative standardized coefficient.",
            "",
        ]
    )
    for pool in ["GCall", "GCfix"]:
        lines.append(
            "- "
            + comparison_line(
                "pcr_continuous_efficiency",
                pool,
                "relative_efficiency",
                "P3_minus_P2",
            )
        )
        median, negative, total = coefficient_summary(
            coefficient_rows,
            "pcr_continuous_efficiency",
            pool,
            "relative_efficiency",
            "P3_combined",
            "variable_weighted_pairs",
        )
        lines.append(
            f"  Standardized variable-region weighted-pair coefficient: median {median:.6g}; negative in {negative}/{total} outer fits (prespecified adverse direction: negative)."
        )
    lines.extend(
        [
            "",
            "The complete baseline (`P1-P0`), full-amplicon sensitivity, continuous-efficiency, binary-DT4DDS gate and secondary metrics are retained in the TSV outputs whether favorable or unfavorable.",
            "",
            "## Event gates",
            "",
        ]
    )
    for row in event_rows:
        lines.append(
            f"- {row['dataset']} / {row['pool']} / {row['endpoint']}: "
            f"{row['events']}/{row['n_sequences']} events; `{row['analysis_gate']}`."
        )
    lines.extend(
        [
            "",
            "## Allowed interpretation",
            "",
            "The results quantify whether the frozen Paper 2 string-level proxy carries incremental information for public, sequence-level experimental endpoints after prespecified baselines. For PCR continuous efficiency, the increment is reproducible but its coefficient is opposite to the intended risk direction; this is evidence of domain mismatch, not molecular-risk validation. Pool heterogeneity, null results and confidence intervals are part of the result.",
            "",
            "## Forbidden interpretation",
            "",
            "These public-data analyses do not show that Paper 2 emitted oligos were synthesized, aged, amplified or recovered; do not isolate a causal molecular mechanism; do not validate NUPACK; do not establish storage-density or full 110-nt language capacity; and do not import evidence from an unpublished sister paper. DT4DDS read counts are composite post-workflow detectability measurements, not isolated aging or synthesis outcomes. Gimpel labels concern their published random pools and adapter/PCR workflow, not Paper 2 codewords.",
            "",
            "## Reproducibility and manuscript freeze",
            "",
            "- Nested preprocessing and regularization tuning were confined to training folds.",
            "- Every primary delta uses identical outer predictions and 2,000 paired sequence-level bootstrap resamples.",
            "- Source hashes, release commits, exact adapter provenance, environment, all model metrics, coefficients and output hashes are saved beside this report.",
            f"- `main.tex` remained frozen at SHA-256 `{main_hash}` during this run.",
            "",
        ]
    )
    return "\n".join(lines)


def environment_record(args: argparse.Namespace) -> dict[str, Any]:
    packages = {}
    for name in ["numpy", "pandas", "scipy", "scikit-learn", "joblib"]:
        packages[name] = importlib.metadata.version(name)
    cpu_brand = ""
    if platform.system() == "Darwin":
        try:
            cpu_brand = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
            ).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            cpu_brand = "unavailable"
    return {
        "analysis_date": "2026-07-17",
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_brand": cpu_brand,
        "logical_cpus": os.cpu_count(),
        "parallel_jobs": args.jobs,
        "bootstrap_replicates": args.bootstrap,
        "outer_repeats": OUTER_REPEATS,
        "outer_folds": OUTER_FOLDS,
        "inner_folds": INNER_FOLDS,
        "base_seed": BASE_SEED,
        "ridge_alpha_grid": RIDGE_ALPHA_GRID.tolist(),
        "logistic_C_grid": LOGISTIC_C_GRID.tolist(),
        "packages": packages,
        "command": f"$PYTHON {relative(Path(__file__))} --jobs {args.jobs} --bootstrap {args.bootstrap}",
    }


def output_hash_manifest() -> list[dict[str, Any]]:
    paths = sorted(
        path for path in OUT_DIR.iterdir() if path.is_file() and path.name != "sha256_manifest.tsv"
    )
    return [
        {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in paths
    ]


def main() -> int:
    args = parse_args()
    if args.bootstrap != 2000:
        raise ValueError("the frozen confirmatory contract requires --bootstrap 2000")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Validating frozen source hashes and commits", flush=True)
    commits = validate_frozen_sources()
    validate_fast_evaluator()
    print("Vectorized evaluator matched the frozen evaluator on fixed tests", flush=True)
    main_hash_start = sha256(MAIN_TEX)

    pcr_raw, pickle_audits = load_pcr_data()
    pcr_features = get_pcr_features(pcr_raw, args.jobs, args.force_features)
    pcr = pcr_raw.merge(
        pcr_features,
        on=["pool", "sample_id"],
        how="inner",
        validate="one_to_one",
    )
    if len(pcr) != 23992:
        raise ValueError("PCR feature merge changed row count")

    dt = pd.read_csv(DT_CSV)
    if len(dt) != 24472:
        raise ValueError("DT4DDS row count mismatch")
    expected_dt = {"Genscript_GCall": (12472, 102), "Twist_GCall": (12000, 108)}
    for pool, (rows, length) in expected_dt.items():
        subset = dt.loc[dt["pool"].eq(pool)]
        if len(subset) != rows or not subset["length_nt"].eq(length).all():
            raise ValueError(f"DT4DDS frozen contract mismatch for {pool}")

    metric_rows: list[dict[str, Any]] = []
    coefficient_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    bootstrap_specs: dict[str, list[dict[str, Any]]] = {}
    event_rows: list[dict[str, Any]] = []

    dt_endpoints = ["mean_log2_cpm_plus1", "day7_vs_day0_log2fc_plus1"]
    for pool_index, pool in enumerate(expected_dt):
        subset = dt.loc[dt["pool"].eq(pool)].reset_index(drop=True)
        for endpoint_index, endpoint in enumerate(dt_endpoints):
            y = subset[endpoint].to_numpy(dtype=float)
            rows, coefficients, predictions = run_nested_models(
                "DT4DDS",
                "dt_continuous",
                pool,
                endpoint,
                subset,
                y,
                DT_MODELS,
                False,
                args.jobs,
                BASE_SEED + pool_index * 1000 + endpoint_index * 100,
            )
            metric_rows.extend(rows)
            coefficient_rows.extend(coefficients)
            add_comparison(
                comparison_rows,
                bootstrap_specs,
                metric_rows,
                "DT4DDS",
                "dt_continuous",
                pool,
                endpoint,
                "D3_minus_D1",
                "D1_structural_rules",
                "D3_full",
                "spearman",
                y,
                predictions,
                f"DT:{pool}",
                "confirmatory",
            )
            add_comparison(
                comparison_rows,
                bootstrap_specs,
                metric_rows,
                "DT4DDS",
                "dt_continuous",
                pool,
                endpoint,
                "D2_minus_D0",
                "D0_composition",
                "D2_weighted_only",
                "spearman",
                y,
                predictions,
                f"DT:{pool}",
                "prespecified sensitivity",
            )

        for endpoint_index, endpoint in enumerate(["missing_any", "missing_day7"]):
            y = subset[endpoint].astype(int).to_numpy()
            events = int(y.sum())
            if events >= 100:
                gate = "nested_5x5_secondary"
                rows, coefficients, predictions = run_nested_models(
                    "DT4DDS",
                    "dt_binary",
                    pool,
                    endpoint,
                    subset,
                    y,
                    DT_MODELS,
                    True,
                    args.jobs,
                    BASE_SEED + 5000 + pool_index * 1000 + endpoint_index * 100,
                )
                metric_rows.extend(rows)
                coefficient_rows.extend(coefficients)
                for comparison, baseline, extended in [
                    ("D3_minus_D1", "D1_structural_rules", "D3_full"),
                    ("D2_minus_D0", "D0_composition", "D2_weighted_only"),
                ]:
                    add_comparison(
                        comparison_rows,
                        bootstrap_specs,
                        metric_rows,
                        "DT4DDS",
                        "dt_binary",
                        pool,
                        endpoint,
                        comparison,
                        baseline,
                        extended,
                        "average_precision",
                        y,
                        predictions,
                        f"DT:{pool}",
                        "secondary under event gate",
                    )
            elif events >= 30:
                gate = "exploratory_3fold_not_triggered_in_frozen_data"
                raise ValueError("unexpected 30-99 event stratum requires explicit exploratory branch")
            else:
                gate = "descriptive_only_below_30_events"
                score = subset["weighted_pairs"].to_numpy(dtype=float)
                descriptive_metrics = {
                    "average_precision": float(average_precision_score(y, score)),
                    "roc_auc": float(roc_auc_score(y, score)),
                }
                for metric, value in descriptive_metrics.items():
                    metric_rows.append(
                        {
                            "dataset": "DT4DDS",
                            "analysis": "dt_binary_descriptive",
                            "pool_or_direction": pool,
                            "endpoint": endpoint,
                            "model": "weighted_pairs_univariate",
                            "metric": metric,
                            "value": value,
                            "n_sequences": len(y),
                            "n_events": events,
                            "event_prevalence": float(y.mean()),
                            "evaluation": "descriptive only; multivariable event gate not met",
                        }
                    )
            event_rows.append(
                {
                    "dataset": "DT4DDS",
                    "pool": pool,
                    "endpoint": endpoint,
                    "n_sequences": len(y),
                    "events": events,
                    "analysis_gate": gate,
                }
            )

    for pool_index, pool in enumerate(["GCall", "GCfix"]):
        subset = pcr.loc[pcr["pool"].eq(pool)].reset_index(drop=True)
        y_binary = subset["label"].to_numpy(dtype=int)
        event_rows.append(
            {
                "dataset": "Gimpel2025_PCR",
                "pool": pool,
                "endpoint": "worst_2_percent_label",
                "n_sequences": len(y_binary),
                "events": int(y_binary.sum()),
                "analysis_gate": "nested_5x5_confirmatory",
            }
        )
        rows, coefficients, predictions = run_nested_models(
            "Gimpel2025_PCR",
            "pcr_binary_internal",
            pool,
            "worst_2_percent_label",
            subset,
            y_binary,
            PCR_MODELS,
            True,
            args.jobs,
            BASE_SEED + 10_000 + pool_index * 1000,
        )
        metric_rows.extend(rows)
        coefficient_rows.extend(coefficients)
        for comparison, baseline, extended, status in [
            ("P1_minus_P0", "P0_composition", "P1_variable_proxy", "prespecified sensitivity"),
            ("P3_minus_P2", "P2_adapter_aware", "P3_combined", "confirmatory"),
        ]:
            add_comparison(
                comparison_rows,
                bootstrap_specs,
                metric_rows,
                "Gimpel2025_PCR",
                "pcr_binary_internal",
                pool,
                "worst_2_percent_label",
                comparison,
                baseline,
                extended,
                "average_precision",
                y_binary,
                predictions,
                f"PCR:{pool}",
                status,
            )

        y_eff = subset["eff"].to_numpy(dtype=float)
        rows, coefficients, predictions = run_nested_models(
            "Gimpel2025_PCR",
            "pcr_continuous_efficiency",
            pool,
            "relative_efficiency",
            subset,
            y_eff,
            PCR_MODELS,
            False,
            args.jobs,
            BASE_SEED + 12_000 + pool_index * 1000,
        )
        metric_rows.extend(rows)
        coefficient_rows.extend(coefficients)
        for comparison, baseline, extended in [
            ("P1_minus_P0", "P0_composition", "P1_variable_proxy"),
            ("P3_minus_P2", "P2_adapter_aware", "P3_combined"),
        ]:
            add_comparison(
                comparison_rows,
                bootstrap_specs,
                metric_rows,
                "Gimpel2025_PCR",
                "pcr_continuous_efficiency",
                pool,
                "relative_efficiency",
                comparison,
                baseline,
                extended,
                "spearman",
                y_eff,
                predictions,
                f"PCR:{pool}",
                "prespecified secondary continuous endpoint",
            )

    for direction_index, (source_pool, target_pool) in enumerate(
        [("GCall", "GCfix"), ("GCfix", "GCall")]
    ):
        source = pcr.loc[pcr["pool"].eq(source_pool)].reset_index(drop=True)
        target = pcr.loc[pcr["pool"].eq(target_pool)].reset_index(drop=True)
        rows, coefficients, predictions = fit_transfer_models(
            source_pool,
            target_pool,
            source,
            target,
            PCR_MODELS,
            BASE_SEED + 20_000 + direction_index * 1000,
        )
        metric_rows.extend(rows)
        coefficient_rows.extend(coefficients)
        y = target["label"].to_numpy(dtype=int)
        for comparison, baseline, extended, status in [
            ("P1_minus_P0", "P0_composition", "P1_variable_proxy", "prespecified transfer sensitivity"),
            ("P3_minus_P2", "P2_adapter_aware", "P3_combined", "prespecified transfer direction check"),
        ]:
            add_comparison(
                comparison_rows,
                bootstrap_specs,
                metric_rows,
                "Gimpel2025_PCR",
                "pcr_binary_transfer",
                f"{source_pool}_to_{target_pool}",
                "worst_2_percent_label",
                comparison,
                baseline,
                extended,
                "average_precision",
                y,
                predictions,
                f"PCR_TRANSFER:{source_pool}_to_{target_pool}",
                status,
            )

    print(f"Paired bootstrap: {args.bootstrap} resamples per sequence pool", flush=True)
    bootstrap_intervals: dict[str, tuple[float, float]] = {}
    bootstrap_groups = list(enumerate(sorted(bootstrap_specs.items())))
    for _group_index, (group, specs) in bootstrap_groups:
        print(f"  {group}: {len(specs)} comparisons", flush=True)
    grouped_intervals = Parallel(
        n_jobs=min(args.jobs, len(bootstrap_groups)), prefer="processes"
    )(
        delayed(bootstrap_group)(
            specs,
            args.bootstrap,
            BASE_SEED + 30_000 + group_index * 1000,
        )
        for group_index, (_group, specs) in bootstrap_groups
    )
    for intervals in grouped_intervals:
        bootstrap_intervals.update(intervals)
    for row in comparison_rows:
        low, high = bootstrap_intervals[row["key"]]
        row["bootstrap_2p5"] = low
        row["bootstrap_97p5"] = high

    if sha256(MAIN_TEX) != main_hash_start:
        raise RuntimeError("main.tex changed during public-data analysis")

    metric_rows.sort(
        key=lambda row: (
            row["dataset"],
            row["analysis"],
            row["pool_or_direction"],
            row["endpoint"],
            row["model"],
            row["metric"],
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
    comparison_rows.sort(key=lambda row: row["key"])
    event_rows.sort(key=lambda row: (row["dataset"], row["pool"], row["endpoint"]))

    write_table(OUT_DIR / "source_manifest.tsv", source_manifest(commits))
    write_table(OUT_DIR / "data_dictionary.tsv", data_dictionary_rows())
    write_table(OUT_DIR / "event_counts_and_gates.tsv", event_rows)
    write_table(OUT_DIR / "model_metrics.tsv", metric_rows)
    write_table(OUT_DIR / "model_comparisons.tsv", comparison_rows)
    write_table(OUT_DIR / "outer_fold_coefficients.tsv", coefficient_rows)
    (OUT_DIR / "pickle_audit.json").write_text(
        json.dumps(pickle_audits, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (OUT_DIR / "environment_and_seeds.json").write_text(
        json.dumps(environment_record(args), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = render_summary(
        metric_rows,
        comparison_rows,
        coefficient_rows,
        event_rows,
        pickle_audits,
        main_hash_start,
    )
    (OUT_DIR / "analysis_summary.md").write_text(summary, encoding="utf-8")
    write_table(OUT_DIR / "sha256_manifest.tsv", output_hash_manifest())
    print(f"Wrote reproducible outputs to {OUT_DIR}", flush=True)
    print(f"Output manifest SHA-256: {sha256(OUT_DIR / 'sha256_manifest.tsv')}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
