#!/usr/bin/env python3
"""Post hoc external-KAPA fiber-mapping sensitivity for Paper 2.

The analysis-plan-locked external mapping remains the primary analysis.  This
script applies 32 additional outcome-blind SHA-256 namespaces to the 2,053
eligible external KAPA sequences, retains the first 2,048 sequences under
each namespace and changes only their ordering/fiber grouping.  Frozen
FullContext, AssayContext and released-CNN scores are then evaluated on the
resulting fibers.

Mappings are deterministic algorithmic sensitivity settings, not biological
or experimental replicates.  Their min--median--max ranges are descriptive;
they are not confidence intervals and are not used to narrow uncertainty.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import paper2_deterministic_selection as deterministic
from joblib import Parallel, delayed


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import validate_paper2_public_experimental_data as base  # noqa: E402
import validate_paper2_reversible_choice_codec as codec  # noqa: E402
import validate_paper2_sota_external_uncertainty as sota  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
SOTA_DIR = PAPER_DIR / "bioinformatics_reframe" / "sota_and_external_validation"
OUT_DIR = PAPER_DIR / "bioinformatics_reframe" / "external_mapping_sensitivity"
MAPPING_INPUTS = OUT_DIR / "external_kapa_mapping_inputs.tsv"
PRIMARY_NAMESPACE = "paper2-choice-external-v1|external_Taq"
SENSITIVITY_NAMESPACES = tuple(
    f"paper2-choice-external-kapa-mapping-sensitivity-v1|mapping={index:02d}"
    for index in range(32)
)
SOURCE_MODELS = ("GCall", "GCfix")
SELECTOR_BITS = (2, 4)
LIBRARY_SIZE = 2048
CLAIM_BOUNDARY = (
    "post hoc deterministic mapping sensitivity on public external-KAPA "
    "measurements; mappings are not biological repeats, do not define a "
    "confidence interval and do not replace the analysis-plan-locked primary mapping"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(12, (os.cpu_count() or 2) - 2)),
        help="Parallel jobs used only when computing features for uncached eligible sequences.",
    )
    parser.add_argument(
        "--rebuild-inputs",
        action="store_true",
        help="Rebuild the frozen 2,053-row mapping input table from public inputs.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def hash_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def write_table(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def predict_from_coefficients(
    frame: pd.DataFrame, coefficients: pd.DataFrame, source: str, model: str
) -> np.ndarray:
    rows = coefficients.loc[
        coefficients["source_pool"].eq(source) & coefficients["model"].eq(model)
    ].copy()
    if rows.empty:
        raise ValueError(f"missing coefficients for {source}/{model}")
    intercepts = rows["raw_intercept"].astype(float).unique()
    if len(intercepts) != 1:
        raise ValueError(f"non-unique raw intercept for {source}/{model}")
    prediction = np.full(len(frame), float(intercepts[0]), dtype=float)
    for row in rows.itertuples(index=False):
        prediction += frame[row.feature].to_numpy(float) * float(row.raw_coefficient)
    return prediction


def build_mapping_inputs(jobs: int) -> pd.DataFrame:
    random, _audits = sota.load_external_pool()
    eligible = random.loc[random["sequence"].map(codec.exact_gc_hp_eligible)].copy()
    if len(eligible) != 2053:
        raise ValueError(f"expected 2,053 eligible KAPA rows, observed {len(eligible)}")

    cached = pd.read_csv(
        SOTA_DIR / "external_locked_library_features.tsv",
        sep="\t",
        dtype={"sample_id": str},
    )
    cached_ids = set(cached["sample_id"].astype(str))
    missing = eligible.loc[
        ~eligible["seq_id_anonymized"].astype(str).isin(cached_ids)
    ].copy()
    if len(missing) != 5:
        raise ValueError(f"expected five uncached eligible sequences, observed {len(missing)}")
    tasks = [
        ("External_KAPA", str(row.seq_id_anonymized), str(row.sequence))
        for row in missing[["seq_id_anonymized", "sequence"]].itertuples(index=False)
    ]
    new_rows = Parallel(n_jobs=min(jobs, len(tasks)), prefer="processes")(
        delayed(base.pcr_feature_row)(*task) for task in tasks
    )
    features = pd.concat([cached, pd.DataFrame(new_rows)], ignore_index=True)
    features["sample_id"] = features["sample_id"].astype(str)
    if len(features) != 2053 or features["sample_id"].duplicated().any():
        raise ValueError("eligible external feature table is not one-to-one")

    merged = eligible.rename(columns={"seq_id_anonymized": "sample_id"}).merge(
        features.drop(columns=["pool"]),
        on=["sample_id", "sequence_sha256"],
        how="left",
        validate="one_to_one",
    )
    required_features = sorted(
        set(codec.MODEL_FEATURES["P2_assay_context"])
        | set(codec.MODEL_FEATURES["P5_combined_context"])
    )
    if merged[required_features].isna().any().any():
        raise ValueError("eligible external feature join is incomplete")

    coefficients = pd.read_csv(SOTA_DIR / "source_model_coefficients.tsv", sep="\t")
    for source in SOURCE_MODELS:
        for model in ("P2_assay_context", "P5_combined_context"):
            merged[f"{source}__{model}"] = predict_from_coefficients(
                merged, coefficients, source, model
            )

    columns = [
        "sample_id",
        "sequence_sha256",
        "eff_Taq",
        "external_low_efficiency",
        "GCall 2perc",
        "GCfix 2perc",
        "GCall__P2_assay_context",
        "GCall__P5_combined_context",
        "GCfix__P2_assay_context",
        "GCfix__P5_combined_context",
    ]
    result = merged[columns].sort_values("sequence_sha256", kind="mergesort").reset_index(drop=True)
    result.to_csv(MAPPING_INPUTS, sep="\t", index=False, float_format="%.15g")
    return result


def ordered_library(frame: pd.DataFrame, namespace: str) -> pd.DataFrame:
    library = frame.copy()
    library["mapping_key"] = library["sequence_sha256"].map(
        lambda value: hash_text(f"{namespace}|{value}")
    )
    library = (
        library.sort_values(["mapping_key", "sequence_sha256"], kind="mergesort")
        .iloc[:LIBRARY_SIZE]
        .reset_index(drop=True)
    )
    if len(library) != LIBRARY_SIZE or library["sequence_sha256"].duplicated().any():
        raise ValueError("mapping did not produce a unique 2,048-sequence library")
    return library


def mapping_metrics(
    library: pd.DataFrame, source: str, bits: int
) -> dict[str, float | int]:
    width = 1 << bits
    fibers = len(library) // width
    outcomes = library["eff_Taq"].to_numpy(float).reshape(fibers, width)
    fiber_mean = outcomes.mean(axis=1)
    p5_score = library[f"{source}__P5_combined_context"].to_numpy(float).reshape(fibers, width)
    p2_score = library[f"{source}__P2_assay_context"].to_numpy(float).reshape(fibers, width)
    cnn_score = library[f"{source} 2perc"].to_numpy(float).reshape(fibers, width)
    rows = np.arange(fibers)
    p5_choice = deterministic.argmax_smallest(p5_score, axis=1)
    p2_choice = deterministic.argmax_smallest(p2_score, axis=1)
    cnn_choice = deterministic.argmin_smallest(cnn_score, axis=1)
    p5_selected = outcomes[rows, p5_choice]
    p2_selected = outcomes[rows, p2_choice]
    cnn_selected = outcomes[rows, cnn_choice]
    p5_gain = p5_selected - fiber_mean
    p5_stability = deterministic.choice_stability_audit(p5_score, "max")
    return {
        "fibers": fibers,
        "fullcontext_gain": float(np.mean(p5_gain)),
        "assaycontext_gain": float(np.mean(p2_selected - fiber_mean)),
        "released_cnn_gain": float(np.mean(cnn_selected - fiber_mean)),
        "fullcontext_minus_assaycontext": float(np.mean(p5_selected - p2_selected)),
        "fullcontext_minus_released_cnn": float(np.mean(p5_selected - cnn_selected)),
        "fullcontext_negative_gain_fraction": float(np.mean(p5_gain < 0)),
        "fullcontext_zero_gain_fraction": float(np.mean(p5_gain == 0)),
        "fullcontext_positive_gain_fraction": float(np.mean(p5_gain > 0)),
        "fullcontext_fiber_score_variance_mean": float(
            p5_score.var(axis=1, ddof=0).mean()
        ),
        "fullcontext_fiber_outcome_variance_mean": float(
            outcomes.var(axis=1, ddof=0).mean()
        ),
        "fullcontext_choices_changed_by_quantization": int(
            p5_stability["choices_changed_by_quantization"]
        ),
        "fullcontext_quantized_extremum_tie_fibers": int(
            p5_stability["quantized_extremum_tie_fibers"]
        ),
        "fullcontext_minimum_raw_top_gap": float(
            p5_stability["minimum_raw_top_gap"]
        ),
    }


def verify_primary_mapping(frame: pd.DataFrame) -> list[dict[str, Any]]:
    library = ordered_library(frame, PRIMARY_NAMESPACE)
    frozen = pd.read_csv(SOTA_DIR / "fiber_benchmark_summary.tsv", sep="\t")
    checks: list[dict[str, Any]] = []
    selector_to_metric = {
        "P5_combined_context": "fullcontext_gain",
        "P2_assay_context": "assaycontext_gain",
        "published_1dcnn": "released_cnn_gain",
    }
    for source in SOURCE_MODELS:
        direction = f"{source}_to_external_Taq"
        for bits in SELECTOR_BITS:
            observed = mapping_metrics(library, source, bits)
            for selector, metric in selector_to_metric.items():
                row = frozen.loc[
                    frozen["dataset"].eq("Gimpel2025_external_laboratory_Taq")
                    & frozen["direction"].eq(direction)
                    & frozen["selector_bits"].eq(bits)
                    & frozen["selector_model"].eq(selector)
                ]
                if len(row) != 1:
                    raise ValueError(f"missing frozen primary row for {source}/{bits}/{selector}")
                expected = float(row.iloc[0]["paired_mean_outcome_gain"])
                difference = float(observed[metric]) - expected
                if abs(difference) > 2e-9:
                    raise ValueError(
                        f"primary mapping mismatch for {source}/r={bits}/{selector}: {difference}"
                    )
                checks.append(
                    {
                        "source_model": source,
                        "selector_bits": bits,
                        "selector": selector,
                        "recomputed_gain": observed[metric],
                        "frozen_gain": expected,
                        "difference": difference,
                        "tolerance": 2e-9,
                        "status": "PASS",
                    }
                )
    return checks


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(rows)
    directional_metrics = (
        "fullcontext_gain",
        "fullcontext_minus_assaycontext",
        "fullcontext_minus_released_cnn",
    )
    descriptive_metrics = (
        "fullcontext_negative_gain_fraction",
        "fullcontext_positive_gain_fraction",
        "fullcontext_fiber_score_variance_mean",
        "fullcontext_fiber_outcome_variance_mean",
        "fullcontext_minimum_raw_top_gap",
    )
    summaries: list[dict[str, Any]] = []
    for (source, bits), group in frame.groupby(["source_model", "selector_bits"], sort=True):
        row: dict[str, Any] = {
            "source_model": source,
            "selector_bits": int(bits),
            "mappings": len(group),
            "mapping_role": "post_hoc_exploratory_algorithmic_sensitivity_not_replicates",
        }
        for metric in directional_metrics:
            values = group[metric].to_numpy(float)
            row[f"{metric}_min"] = float(np.min(values))
            row[f"{metric}_median"] = float(np.median(values))
            row[f"{metric}_max"] = float(np.max(values))
            row[f"{metric}_positive_mappings"] = int(np.sum(values > 0))
        for metric in descriptive_metrics:
            values = group[metric].to_numpy(float)
            row[f"{metric}_min"] = float(np.min(values))
            row[f"{metric}_median"] = float(np.median(values))
            row[f"{metric}_max"] = float(np.max(values))
        row["fullcontext_choices_changed_by_quantization_total"] = int(
            group["fullcontext_choices_changed_by_quantization"].sum()
        )
        row["fullcontext_quantized_extremum_tie_fibers_total"] = int(
            group["fullcontext_quantized_extremum_tie_fibers"].sum()
        )
        summaries.append(row)
    return summaries


def output_manifest() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(OUT_DIR.iterdir()):
        if not path.is_file() or path.name == "sha256_manifest.tsv":
            continue
        rows.append(
            {
                "path": str(path.resolve().relative_to(ROOT.resolve())),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    rows.append(
        {
            "path": str(Path(__file__).resolve().relative_to(ROOT.resolve())),
            "bytes": Path(__file__).stat().st_size,
            "sha256": sha256(Path(__file__)),
        }
    )
    return rows


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.rebuild_inputs or not MAPPING_INPUTS.is_file():
        inputs = build_mapping_inputs(args.jobs)
        input_mode = "rebuilt_from_hash_verified_public_inputs"
    else:
        inputs = pd.read_csv(MAPPING_INPUTS, sep="\t", dtype={"sample_id": str})
        input_mode = "frozen_derived_mapping_inputs"
    if len(inputs) != 2053 or inputs["sequence_sha256"].duplicated().any():
        raise ValueError("mapping input table must contain 2,053 unique eligible sequences")

    primary_checks = verify_primary_mapping(inputs)
    mapping_rows: list[dict[str, Any]] = []
    for mapping_index, namespace in enumerate(SENSITIVITY_NAMESPACES):
        library = ordered_library(inputs, namespace)
        membership_hash = hash_text("\n".join(sorted(library["sequence_sha256"])))
        order_hash = hash_text("\n".join(library["sequence_sha256"]))
        for source in SOURCE_MODELS:
            for bits in SELECTOR_BITS:
                mapping_rows.append(
                    {
                        "mapping_index": mapping_index,
                        "namespace_sha256": hash_text(namespace),
                        "membership_sha256": membership_hash,
                        "ordered_library_sha256": order_hash,
                        "source_model": source,
                        "selector_bits": bits,
                        "library_sequences": len(library),
                        "eligible_sequences": len(inputs),
                        **mapping_metrics(library, source, bits),
                        "mapping_role": "post_hoc_exploratory_algorithmic_sensitivity_not_replicates",
                    }
                )
    summaries = summarize(mapping_rows)
    write_table(OUT_DIR / "primary_mapping_reproduction_check.tsv", primary_checks)
    write_table(OUT_DIR / "external_kapa_mapping_sensitivity.tsv", mapping_rows)
    write_table(OUT_DIR / "external_kapa_mapping_sensitivity_summary.tsv", summaries)

    environment = {
        "analysis": "Paper 2 post hoc external KAPA mapping sensitivity",
        "command": f"$PYTHON analysis_tools/analyze_paper2_external_mapping_sensitivity.py --jobs {args.jobs}",
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "pandas", "joblib"]
        },
        "input_mode": input_mode,
        "eligible_sequences": len(inputs),
        "library_sequences_per_mapping": LIBRARY_SIZE,
        "mapping_namespaces": len(SENSITIVITY_NAMESPACES),
        "source_models": list(SOURCE_MODELS),
        "selector_bits": list(SELECTOR_BITS),
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
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(OUT_DIR / "environment_and_seeds.json", environment)

    summary_frame = pd.DataFrame(summaries)
    lines = [
        "# External KAPA mapping sensitivity",
        "",
        "Status: `POST_HOC_EXPLORATORY_ALGORITHMIC_SENSITIVITY`",
        "",
        "The analysis-plan-locked mapping remains primary. Thirty-two additional namespaces depend only on a fixed index and sequence SHA-256. Ranges below are deterministic min--median--max summaries, not confidence intervals or biological repeats.",
        "",
    ]
    for row in summary_frame.itertuples(index=False):
        lines.append(
            f"- {row.source_model}, r={row.selector_bits}: FullContext gain "
            f"{row.fullcontext_gain_min:.8f}--{row.fullcontext_gain_median:.8f}--{row.fullcontext_gain_max:.8f} "
            f"({row.fullcontext_gain_positive_mappings}/32 positive); FullContext-minus-AssayContext "
            f"{row.fullcontext_minus_assaycontext_min:.8f}--{row.fullcontext_minus_assaycontext_median:.8f}--{row.fullcontext_minus_assaycontext_max:.8f} "
            f"({row.fullcontext_minus_assaycontext_positive_mappings}/32 positive); FullContext-minus-CNN "
            f"{row.fullcontext_minus_released_cnn_min:.8f}--{row.fullcontext_minus_released_cnn_median:.8f}--{row.fullcontext_minus_released_cnn_max:.8f} "
            f"({row.fullcontext_minus_released_cnn_positive_mappings}/32 positive)."
        )
    lines.extend(["", "Evidence boundary: " + CLAIM_BOUNDARY + ".", ""])
    (OUT_DIR / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")
    write_table(OUT_DIR / "sha256_manifest.tsv", output_manifest())
    print(f"External KAPA mapping sensitivity complete: {OUT_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
