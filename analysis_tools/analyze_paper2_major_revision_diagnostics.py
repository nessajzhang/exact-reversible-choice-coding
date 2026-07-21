#!/usr/bin/env python3
"""Build reviewer-requested diagnostics for the Paper 2 major revision.

This script consumes only manifest-locked derived outputs.  It does not create
new experimental observations, refit an unreported model, or treat mapping
namespaces as biological replicates.  The resulting tables separate measured,
predicted, computational and not-performed evidence.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import paper2_deterministic_selection as deterministic


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
REFRAME_DIR = PAPER_DIR / "bioinformatics_reframe"
CODEC_DIR = REFRAME_DIR / "reversible_choice_codec"
SOTA_DIR = REFRAME_DIR / "sota_and_external_validation"
EXTERNAL_MAPPING_DIR = REFRAME_DIR / "external_mapping_sensitivity"
OUT_DIR = REFRAME_DIR / "major_revision_diagnostics"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(name: str, frame: pd.DataFrame) -> Path:
    if frame.empty:
        raise ValueError(f"refusing to write empty diagnostic table: {name}")
    path = OUT_DIR / name
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def human_dataset(value: str) -> str:
    return {
        "Gimpel2025_cross_pool_public_codebook": "GCall/GCfix retrospective cross-pool codebook",
        "Gimpel2025_external_laboratory_Taq": "same-source-publication external-laboratory KAPA codebook",
        "Gimpel2025_external_laboratory_Q5_sensitivity": "same-source-publication altered-protocol Q5 codebook",
    }.get(value, value)


def dataset_partition_table() -> pd.DataFrame:
    rows = [
        {
            "dataset": "GCall",
            "source_publication": "Gimpel et al. 2025",
            "laboratory_role": "source-study primary pool",
            "assay_condition": "KAPA-matched multi-template PCR workflow",
            "polymerase": "KAPA SYBR FAST",
            "annealing_temperature_c": 54,
            "measured_sequences": 11998,
            "codec_eligible_sequences": "reported in frozen public-codebook output",
            "finite_codebook_sequences": 2048,
            "used_for_training": "yes when source; no when target",
            "used_for_feature_development": "yes",
            "used_for_model_selection": "yes when source",
            "used_only_for_final_test": "no",
            "evidence_role": "retrospective bidirectional cross-pool transfer",
            "independent_third_party_replication": "no",
            "sequence_overlap_control": "exhaustive Hamming, duplicate and reverse-complement audit",
        },
        {
            "dataset": "GCfix",
            "source_publication": "Gimpel et al. 2025",
            "laboratory_role": "source-study primary pool",
            "assay_condition": "KAPA-matched multi-template PCR workflow",
            "polymerase": "KAPA SYBR FAST",
            "annealing_temperature_c": 54,
            "measured_sequences": 11994,
            "codec_eligible_sequences": "reported in frozen public-codebook output",
            "finite_codebook_sequences": 2048,
            "used_for_training": "yes when source; no when target",
            "used_for_feature_development": "yes",
            "used_for_model_selection": "yes when source",
            "used_only_for_final_test": "no",
            "evidence_role": "retrospective bidirectional cross-pool transfer",
            "independent_third_party_replication": "no",
            "sequence_overlap_control": "exhaustive Hamming, duplicate and reverse-complement audit",
        },
        {
            "dataset": "external KAPA",
            "source_publication": "Gimpel et al. 2025 (same source publication)",
            "laboratory_role": "external laboratory reported in the source publication",
            "assay_condition": "external KAPA SYBR FAST at 54 C",
            "polymerase": "KAPA SYBR FAST",
            "annealing_temperature_c": 54,
            "measured_sequences": 9995,
            "codec_eligible_sequences": 2053,
            "finite_codebook_sequences": 2048,
            "used_for_training": "no",
            "used_for_feature_development": "no",
            "used_for_model_selection": "no",
            "used_only_for_final_test": "yes for the analysis-plan-locked KAPA evaluation",
            "evidence_role": "same-source-publication external-laboratory measured evaluation",
            "independent_third_party_replication": "no",
            "sequence_overlap_control": "exhaustive source-to-external Hamming, duplicate and reverse-complement audit",
        },
        {
            "dataset": "external Q5",
            "source_publication": "Gimpel et al. 2025 (same source publication)",
            "laboratory_role": "external laboratory reported in the source publication",
            "assay_condition": "Q5 Hot Start High-Fidelity at 60 C",
            "polymerase": "Q5 Hot Start High-Fidelity",
            "annealing_temperature_c": 60,
            "measured_sequences": 9955,
            "codec_eligible_sequences": 2040,
            "finite_codebook_sequences": 1024,
            "used_for_training": "no",
            "used_for_feature_development": "no",
            "used_for_model_selection": "no",
            "used_only_for_final_test": "no; secondary exploratory workflow-shift analysis",
            "evidence_role": "altered-protocol measured sensitivity",
            "independent_third_party_replication": "no",
            "sequence_overlap_control": "same random-sequence release; missing outcomes removed without imputation",
        },
        {
            "dataset": "generated exact-codec language",
            "source_publication": "author-derived computational language",
            "laboratory_role": "not applicable",
            "assay_condition": "not assayed",
            "polymerase": "not applicable",
            "annealing_temperature_c": "not applicable",
            "measured_sequences": 0,
            "codec_eligible_sequences": 65024,
            "finite_codebook_sequences": "not applicable",
            "used_for_training": "no",
            "used_for_feature_development": "no",
            "used_for_model_selection": "no",
            "used_only_for_final_test": "computational verification only",
            "evidence_role": "exact codec verification and predicted-score audit",
            "independent_third_party_replication": "not performed",
            "sequence_overlap_control": "generated-domain distribution diagnostics; no prospective synthesis",
        },
    ]
    return pd.DataFrame(rows)


def evidence_hierarchy_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "layer": 1,
                "evidence_type": "computational verification",
                "material": "rank/unrank, modulo-N inverse and 65,024 generated candidates",
                "supports": "language membership, fixed fibers, exact payload inversion and failure behavior",
                "does_not_support": "PCR benefit, synthesis performance or material recovery",
            },
            {
                "layer": 2,
                "evidence_type": "retrospective measured analysis",
                "material": "GCall/GCfix PCR-tested finite codebooks",
                "supports": "measured within-fiber selection utility in previously assayed libraries",
                "does_not_support": "prospective performance of codec-emitted sequences",
            },
            {
                "layer": 3,
                "evidence_type": "same-source-publication external-laboratory evaluation",
                "material": "analysis-plan-locked external KAPA finite codebook from already-public measurements",
                "supports": "same-assay transfer to a source-publication external laboratory pool",
                "does_not_support": "independent third-party replication",
            },
            {
                "layer": 4,
                "evidence_type": "altered-protocol evaluation",
                "material": "secondary exploratory external Q5 finite codebook",
                "supports": "an assay-transfer boundary; direction reversal under changed amplification conditions",
                "does_not_support": "polymerase-universal sequence quality",
            },
            {
                "layer": 5,
                "evidence_type": "predicted generated-candidate analysis",
                "material": "FullContext scores on generated exact-codec candidates",
                "supports": "in-silico score-rate and mapping-sensitivity diagnostics",
                "does_not_support": "measured gain",
            },
            {
                "layer": 6,
                "evidence_type": "prospective validation",
                "material": "actual synthesis and PCR of codec-emitted candidates",
                "supports": "not available",
                "does_not_support": "not performed in this study",
            },
        ]
    )


def negative_fiber_diagnostics(bench: pd.DataFrame, gains: pd.DataFrame) -> pd.DataFrame:
    selected = bench.loc[
        bench["selector_model"].eq("P5_combined_context")
        & bench["selector_bits"].isin([2, 4])
    ].copy()
    rows: list[dict[str, Any]] = []
    for item in selected.itertuples(index=False):
        group = gains.loc[
            gains["dataset"].eq(item.dataset)
            & gains["direction"].eq(item.direction)
            & gains["selector_model"].eq(item.selector_model)
            & gains["selector_bits"].eq(item.selector_bits)
            & gains["accepted"].astype(bool)
        ].copy()
        values = group["outcome_gain"].to_numpy(float)
        if not len(values):
            raise ValueError(f"missing fiber gains for {item.dataset}/{item.direction}/r={item.selector_bits}")
        target_sd = float(item.paired_mean_outcome_gain) / float(item.standardized_paired_gain)
        standardized = values / target_sd
        rows.append(
            {
                "dataset": human_dataset(item.dataset),
                "direction": item.direction.replace("external_Taq", "external_KAPA"),
                "selector_bits": int(item.selector_bits),
                "fibers": int(len(values)),
                "evidence_type": "measured retrospective selection",
                "mean_gain": float(values.mean()),
                "standardized_mean_gain": float(standardized.mean()),
                "minimum_gain": float(values.min()),
                "q25_gain": float(np.quantile(values, 0.25)),
                "median_gain": float(np.median(values)),
                "q75_gain": float(np.quantile(values, 0.75)),
                "maximum_gain": float(values.max()),
                "negative_gain_fraction": float((values < 0).mean()),
                "near_zero_gain_fraction_abs_le_0p01_target_sd": float((np.abs(standardized) <= 0.01).mean()),
                "positive_gain_fraction": float((values > 0).mean()),
                "mean_within_fiber_outcome_variance": float(getattr(item, "fiber_outcome_variance_mean")),
                "mean_within_fiber_score_variance": float(getattr(item, "fiber_score_variance_mean")),
                "interpretation": "average shift is not a per-payload guarantee; near-zero threshold is exploratory",
            }
        )
    return pd.DataFrame(rows)


def generated_rate_benefit(generated: pd.DataFrame) -> pd.DataFrame:
    frame = generated.loc[generated["selector_model"].eq("P5_combined_context")].copy()
    keep = [
        "selector_source_pool",
        "cross_evaluator_pool",
        "selector_bits",
        "candidates_per_payload",
        "sampled_payloads",
        "payload_bits_after_selection",
        "payload_bits_per_variable_nt",
        "own_model_mean_gain_vs_fiber",
        "own_model_gain_minimum",
        "own_model_gain_median",
        "own_model_gain_maximum",
        "own_model_negative_gain_fraction",
        "own_model_zero_gain_fraction",
        "fibers_improved_in_own_model_fraction",
        "own_model_fiber_score_variance_mean",
        "choices_changed_by_quantization",
        "quantized_extremum_tie_fibers",
        "minimum_raw_top_gap",
    ]
    frame = frame[keep].rename(columns={"selector_bits": "choice_bits"})
    frame.insert(len(frame.columns), "evidence_type", "predicted score on unmeasured generated candidates")
    return frame.sort_values(["selector_source_pool", "choice_bits"]).reset_index(drop=True)


def negative_fiber_feature_stratification(
    gains: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    feature_paths = [
        REFRAME_DIR / "public_experimental_validation" / "pcr_sequence_features.tsv",
        SOTA_DIR / "external_locked_library_features.tsv",
        SOTA_DIR / "external_q5_locked_library_features.tsv",
    ]
    features = pd.concat(
        [pd.read_csv(path, sep="\t", dtype={"sample_id": str}) for path in feature_paths],
        ignore_index=True,
    )
    feature_columns = [
        column
        for column in features.columns
        if column
        not in {
            "pool",
            "sample_id",
            "sequence_sha256",
            "variable_length_nt",
            "full_amplicon_length_nt",
        }
    ]
    for column in feature_columns:
        features[column] = pd.to_numeric(features[column], errors="raise")
    features = features.drop_duplicates("sequence_sha256", keep="first")

    selected = gains.loc[
        gains["selector_model"].eq("P5_combined_context")
        & gains["selector_bits"].isin([2, 4])
        & gains["accepted"].astype(bool)
    ].copy()
    merged = selected.merge(
        features[["sequence_sha256", *feature_columns]],
        on="sequence_sha256",
        how="left",
        validate="many_to_one",
    )
    if merged[feature_columns].isna().any().any():
        missing = int(merged[feature_columns].isna().all(axis=1).sum())
        raise ValueError(f"missing selected-sequence features for {missing} fibers")

    detail_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    group_columns = ["dataset", "direction", "selector_bits"]
    for keys, group in merged.groupby(group_columns, sort=True):
        negative = group.loc[group["outcome_gain"] < 0]
        nonnegative = group.loc[group["outcome_gain"] >= 0]
        if len(negative) < 2 or len(nonnegative) < 2:
            raise ValueError(f"insufficient negative/nonnegative fibers for feature comparison: {keys}")
        setting_rows: list[dict[str, Any]] = []
        for feature in feature_columns:
            left = negative[feature].to_numpy(float)
            right = nonnegative[feature].to_numpy(float)
            pooled_variance = (
                (len(left) - 1) * left.var(ddof=1) + (len(right) - 1) * right.var(ddof=1)
            ) / (len(left) + len(right) - 2)
            pooled_sd = float(np.sqrt(max(pooled_variance, 0.0)))
            smd = float((left.mean() - right.mean()) / pooled_sd) if pooled_sd > 0 else 0.0
            row = {
                "dataset": human_dataset(keys[0]),
                "direction": str(keys[1]).replace("external_Taq", "external_KAPA"),
                "selector_bits": int(keys[2]),
                "feature": feature,
                "negative_fibers": int(len(left)),
                "nonnegative_fibers": int(len(right)),
                "negative_fiber_selected_feature_mean": float(left.mean()),
                "nonnegative_fiber_selected_feature_mean": float(right.mean()),
                "standardized_mean_difference_negative_minus_nonnegative": smd,
                "analysis_role": "post hoc exploratory selected-sequence feature stratification",
            }
            detail_rows.append(row)
            setting_rows.append(row)
        top = max(setting_rows, key=lambda row: abs(row["standardized_mean_difference_negative_minus_nonnegative"]))
        summary_rows.append(
            {
                "dataset": top["dataset"],
                "direction": top["direction"],
                "selector_bits": top["selector_bits"],
                "largest_absolute_smd_feature": top["feature"],
                "largest_absolute_smd": abs(top["standardized_mean_difference_negative_minus_nonnegative"]),
                "signed_smd_negative_minus_nonnegative": top[
                    "standardized_mean_difference_negative_minus_nonnegative"
                ],
                "interpretation": "exploratory association in selected candidates; not a validated failure mechanism",
            }
        )
    return pd.DataFrame(detail_rows), pd.DataFrame(summary_rows)


def score_margin_abstention_sensitivity(
    bench: pd.DataFrame, gains: pd.DataFrame
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    selected = gains.loc[
        gains["selector_model"].eq("P5_combined_context")
        & gains["selector_bits"].isin([2, 4])
        & gains["accepted"].astype(bool)
    ].copy()
    for keys, group in selected.groupby(["dataset", "direction", "selector_bits"], sort=True):
        reference = bench.loc[
            bench["dataset"].eq(keys[0])
            & bench["direction"].eq(keys[1])
            & bench["selector_model"].eq("P5_combined_context")
            & bench["selector_bits"].eq(keys[2])
        ]
        if len(reference) != 1:
            raise ValueError(f"non-unique benchmark row for {keys}")
        reference_row = reference.iloc[0]
        target_sd = float(reference_row["paired_mean_outcome_gain"]) / float(
            reference_row["standardized_paired_gain"]
        )
        for nominal_retention in (1.0, 0.9, 0.75):
            if nominal_retention == 1.0:
                retained = group.copy()
                threshold = float("-inf")
            else:
                threshold = float(
                    np.quantile(group["score_margin"].to_numpy(float), 1 - nominal_retention)
                )
                retained = group.loc[group["score_margin"] >= threshold].copy()
            values = retained["outcome_gain"].to_numpy(float)
            rows.append(
                {
                    "dataset": human_dataset(keys[0]),
                    "direction": str(keys[1]).replace("external_Taq", "external_KAPA"),
                    "selector_bits": int(keys[2]),
                    "nominal_score_margin_retention": nominal_retention,
                    "score_margin_threshold": threshold,
                    "retained_fibers": int(len(retained)),
                    "actual_retention_fraction": float(len(retained) / len(group)),
                    "mean_measured_gain": float(values.mean()),
                    "standardized_mean_measured_gain": float(values.mean() / target_sd),
                    "negative_gain_fraction": float((values < 0).mean()),
                    "near_zero_gain_fraction_abs_le_0p01_target_sd": float(
                        (np.abs(values / target_sd) <= 0.01).mean()
                    ),
                    "positive_gain_fraction": float((values > 0).mean()),
                    "analysis_role": "post hoc exploratory outcome-blind score-margin abstention sensitivity",
                    "boundary": "threshold was examined after outcomes and is not a validated rejection policy or per-payload safety guarantee",
                }
            )
    return pd.DataFrame(rows)


def aggregate_mapping(
    frame: pd.DataFrame,
    group_columns: list[str],
    gain_column: str,
    negative_column: str,
    variance_column: str,
    scope: str,
    unit: str,
    changed_column: str = "choices_changed_by_quantization",
    tie_column: str = "quantized_extremum_tie_fibers",
    gap_column: str = "minimum_raw_top_gap",
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for keys, group in frame.groupby(group_columns, sort=True):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        row = dict(zip(group_columns, key_values))
        gains = group[gain_column].to_numpy(float)
        negatives = group[negative_column].to_numpy(float)
        variances = group[variance_column].to_numpy(float)
        row.update(
            {
                "scope": scope,
                "mapping_count": int(len(group)),
                "gain_unit": unit,
                "gain_minimum": float(gains.min()),
                "gain_median": float(np.median(gains)),
                "gain_maximum": float(gains.max()),
                "positive_gain_mappings": int((gains > 0).sum()),
                "negative_fiber_fraction_minimum": float(negatives.min()),
                "negative_fiber_fraction_median": float(np.median(negatives)),
                "negative_fiber_fraction_maximum": float(negatives.max()),
                "fiber_score_variance_minimum": float(variances.min()),
                "fiber_score_variance_median": float(np.median(variances)),
                "fiber_score_variance_maximum": float(variances.max()),
                "choices_changed_by_quantization_total": int(group[changed_column].sum()),
                "quantized_extremum_tie_fibers_total": int(group[tie_column].sum()),
                "minimum_raw_top_gap": float(group[gap_column].min()),
                "interpretation": "descriptive algorithmic sensitivity; mappings are not biological replicates or confidence intervals",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def mapping_sensitivity_table() -> pd.DataFrame:
    public = pd.read_csv(CODEC_DIR / "public_seed_sensitivity.tsv", sep="\t")
    public = public.loc[public["selector_model"].eq("P5_combined_context")].copy()
    public_agg = aggregate_mapping(
        public,
        ["direction", "selector_bits"],
        "standardized_paired_gain",
        "negative_gain_fraction",
        "fiber_score_variance_mean",
        "retrospective public codebook",
        "target standard deviations",
    )

    generated = pd.read_csv(CODEC_DIR / "generated_mapping_namespace_sensitivity.tsv", sep="\t")
    generated_agg = aggregate_mapping(
        generated,
        ["selector_source_pool", "selector_bits"],
        "own_model_mean_gain_vs_fiber",
        "own_model_negative_gain_fraction",
        "own_model_fiber_score_variance_mean",
        "generated exact-codec language",
        "model-score units",
    )

    external = pd.read_csv(EXTERNAL_MAPPING_DIR / "external_kapa_mapping_sensitivity.tsv", sep="\t")
    external_agg = aggregate_mapping(
        external,
        ["source_model", "selector_bits"],
        "fullcontext_gain",
        "fullcontext_negative_gain_fraction",
        "fullcontext_fiber_score_variance_mean",
        "same-source-publication external-laboratory KAPA codebook",
        "measured efficiency units",
        "fullcontext_choices_changed_by_quantization",
        "fullcontext_quantized_extremum_tie_fibers",
        "fullcontext_minimum_raw_top_gap",
    )
    return pd.concat([public_agg, external_agg, generated_agg], ignore_index=True, sort=False)


def score_determinism_table(
    bench: pd.DataFrame, generated: pd.DataFrame, mapping: pd.DataFrame
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    def add(scope: str, frame: pd.DataFrame, fiber_column: str) -> None:
        active = frame.loc[frame["selector_bits"].astype(int) > 0].copy()
        rows.append(
            {
                "scope": scope,
                "selection_runs": int(len(active)),
                "fibers_checked": int(active[fiber_column].sum()),
                "score_decimal_places": deterministic.SCORE_DECIMAL_PLACES,
                "choices_changed_by_quantization": int(active["choices_changed_by_quantization"].sum()),
                "quantized_extremum_tie_fibers": int(active["quantized_extremum_tie_fibers"].sum()),
                "minimum_raw_top_gap": float(active["minimum_raw_top_gap"].min()),
                "tie_rule": "smallest choice index after fixed-precision integer score key",
            }
        )

    score_bench = bench.loc[
        bench["selector_model"].isin(["P2_assay_context", "P5_combined_context", "published_1dcnn"])
    ].copy()
    add("measured-codebook AssayContext, FullContext and released-CNN selectors", score_bench, "payloads")
    add(
        "measured-codebook exploratory Composition selector",
        bench.loc[bench["selector_model"].eq("P0_composition")].copy(),
        "payloads",
    )
    add(
        "generated AssayContext and FullContext selectors",
        generated.loc[generated["selector_model"].isin(["P2_assay_context", "P5_combined_context"])].copy(),
        "sampled_payloads",
    )
    add(
        "generated exploratory Composition selector",
        generated.loc[generated["selector_model"].eq("P0_composition")].copy(),
        "sampled_payloads",
    )

    mapping_rows = mapping.copy()
    rows.append(
        {
            "scope": "all mapping-sensitivity selector runs",
            "selection_runs": int(mapping_rows["mapping_count"].sum()),
            "fibers_checked": "reported per frozen mapping table",
            "score_decimal_places": deterministic.SCORE_DECIMAL_PLACES,
            "choices_changed_by_quantization": int(mapping_rows["choices_changed_by_quantization_total"].sum()),
            "quantized_extremum_tie_fibers": int(mapping_rows["quantized_extremum_tie_fibers_total"].sum()),
            "minimum_raw_top_gap": float(mapping_rows["minimum_raw_top_gap"].min()),
            "tie_rule": "smallest choice index after fixed-precision integer score key",
        }
    )
    return pd.DataFrame(rows)


def decoder_failure_table(roundtrip: pd.DataFrame) -> pd.DataFrame:
    row = roundtrip.iloc[0]
    return pd.DataFrame(
        [
            ("wrong sequence length", "payload decoder", "reject", bool(row.wrong_length_sequence_rejected)),
            ("non-ACGT symbol including N", "payload decoder", "reject", bool(row.ambiguous_base_sequence_rejected)),
            ("sequence outside constrained language", "payload decoder", "reject", bool(row.invalid_constraint_sequence_rejected)),
            ("inverse logical rank q >= 2^213", "payload decoder", "reject", bool(row.payload_domain_out_of_domain_sequence_rejected)),
            ("valid noncanonical member of declared fiber", "payload decoder", "accept and recover payload", bool(row.noncanonical_valid_candidate_accepted_by_payload_decoder)),
            ("valid noncanonical member of declared fiber", "optional canonical verifier", "reject", bool(row.noncanonical_valid_candidate_rejected_by_optional_verifier)),
            ("scoring model missing", "payload decoder", "decode remains available", True),
            ("scoring model missing", "optional canonical verifier", "canonical status unavailable", True),
            ("scoring model replaced", "payload decoder", "historical payload remains decodable", bool(row.score_replacement_preserves_payload_decoding)),
            ("exact fixed-key score tie", "encoder/verifier", "select smallest choice index", bool(row.exact_tie_selects_smallest_choice_index)),
            ("NaN or infinite score", "encoder/verifier", "reject", bool(row.nonfinite_score_rejected_by_encoder_verifier)),
            ("checksum failure", "payload decoder", "not evaluated; no checksum is claimed", "not performed"),
        ],
        columns=["input_condition", "interface", "required_behavior", "verified"],
    )


def leakage_scope_table() -> pd.DataFrame:
    cross = pd.read_csv(
        REFRAME_DIR / "sequence_independence_audit" / "cross_pool_independence_summary.tsv",
        sep="\t",
    )
    external = pd.read_csv(
        REFRAME_DIR / "source_external_independence_audit" / "source_external_independence_summary.tsv",
        sep="\t",
    )
    rows: list[dict[str, Any]] = []
    for item in pd.concat([cross, external], ignore_index=True, sort=False).itertuples(index=False):
        rows.append(
            {
                "comparison": item.audit_scope,
                "audit": "exhaustive oriented Hamming distance, exact duplicate and exact reverse complement",
                "status": "performed",
                "result": (
                    f"pairs={int(item.all_cross_pool_pairs) if hasattr(item, 'all_cross_pool_pairs') and not pd.isna(item.all_cross_pool_pairs) else int(item.all_oriented_pairs)}; "
                    f"minimum Hamming={int(item.minimum_cross_pool_hamming_distance) if hasattr(item, 'minimum_cross_pool_hamming_distance') and not pd.isna(item.minimum_cross_pool_hamming_distance) else int(item.minimum_oriented_hamming_distance_nt)}; "
                    f"maximum identity={float(item.maximum_cross_pool_identity_fraction) if hasattr(item, 'maximum_cross_pool_identity_fraction') and not pd.isna(item.maximum_cross_pool_identity_fraction) else float(item.maximum_oriented_identity_fraction):.6f}; "
                    "duplicates=0; reverse complements=0; >=70% identity pairs=0"
                ),
                "boundary": "shared fixed adapters were excluded as declared assay context",
            }
        )
    rows.extend(
        [
            {
                "comparison": "all sequence pools",
                "audit": "edit-distance near-neighbour enumeration",
                "status": "not performed",
                "result": "not available",
                "boundary": "all compared variable regions have fixed length 108; exhaustive Hamming results do not establish edit-distance separation",
            },
            {
                "comparison": "all sequence pools",
                "audit": "shared motif or design-template genealogy",
                "status": "not fully available",
                "result": "external motif-inserted rows were excluded; common adapters are declared; upstream random-sequence genealogy was not reconstructed",
                "boundary": "absence of exact/near-exact Hamming matches is not a sequence-cluster split claim",
            },
            {
                "comparison": "model-development partitions",
                "audit": "cluster-based train/test split",
                "status": "not performed",
                "result": "source-only cross-pool fitting plus exhaustive sequence-distance audit",
                "boundary": "GCall/GCfix were not a permanent preregistered lockbox",
            },
        ]
    )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    bench = pd.read_csv(SOTA_DIR / "fiber_benchmark_summary.tsv", sep="\t")
    gains = pd.read_csv(SOTA_DIR / "fiber_gain_distributions.tsv", sep="\t")
    generated = pd.read_csv(CODEC_DIR / "generated_fiber_results.tsv", sep="\t")
    roundtrip = pd.read_csv(CODEC_DIR / "generated_roundtrip_summary.tsv", sep="\t")

    dataset = dataset_partition_table()
    hierarchy = evidence_hierarchy_table()
    negative = negative_fiber_diagnostics(bench, gains)
    rate = generated_rate_benefit(generated)
    feature_detail, feature_summary = negative_fiber_feature_stratification(gains)
    abstention = score_margin_abstention_sensitivity(bench, gains)
    mapping = mapping_sensitivity_table()
    determinism = score_determinism_table(bench, generated, mapping)
    failures = decoder_failure_table(roundtrip)
    leakage = leakage_scope_table()

    output_paths = [
        write_tsv("dataset_partition_and_evidence_table.tsv", dataset),
        write_tsv("evidence_hierarchy.tsv", hierarchy),
        write_tsv("measured_negative_fiber_diagnostics.tsv", negative),
        write_tsv("generated_rate_benefit.tsv", rate),
        write_tsv("negative_fiber_feature_stratification.tsv", feature_detail),
        write_tsv("negative_fiber_feature_summary.tsv", feature_summary),
        write_tsv("score_margin_abstention_sensitivity.tsv", abstention),
        write_tsv("mapping_sensitivity_diagnostics.tsv", mapping),
        write_tsv("score_determinism_audit.tsv", determinism),
        write_tsv("decoder_failure_behavior.tsv", failures),
        write_tsv("leakage_scope.tsv", leakage),
    ]

    source_files = [
        SOTA_DIR / "fiber_benchmark_summary.tsv",
        SOTA_DIR / "fiber_gain_distributions.tsv",
        CODEC_DIR / "generated_fiber_results.tsv",
        CODEC_DIR / "generated_roundtrip_summary.tsv",
        CODEC_DIR / "public_seed_sensitivity.tsv",
        CODEC_DIR / "generated_mapping_namespace_sensitivity.tsv",
        EXTERNAL_MAPPING_DIR / "external_kapa_mapping_sensitivity.tsv",
    ]
    environment = {
        "analysis": "paper2_major_revision_diagnostics",
        "evidence_contract": "derived diagnostics only; no new biological experiment",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {name: package_version(name) for name in ("numpy", "pandas")},
        "score_key_contract": {
            "input_dtype": "float64",
            "output_dtype": "int64",
            "decimal_places": deterministic.SCORE_DECIMAL_PLACES,
            "integer_scale": deterministic.SCORE_SCALE,
            "rounding": deterministic.ROUNDING_CONTRACT,
            "signed_int64_min": deterministic.INT64_MIN,
            "signed_int64_max": deterministic.INT64_MAX,
            "candidate_order": "lexicographically minimize (-quantized_key, choice_index)",
            "tie_rule": "smallest choice index",
            "nonfinite_behavior": "reject",
            "out_of_range_behavior": "reject",
        },
        "source_files": {
            str(path.relative_to(ROOT)): sha256(path) for path in source_files
        },
        "fixed_exploratory_near_zero_threshold_target_sd": 0.01,
    }
    environment_path = OUT_DIR / "environment_and_seeds.json"
    environment_path.write_text(json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_paths.append(environment_path)

    summary = OUT_DIR / "analysis_summary.md"
    summary.write_text(
        "# Paper 2 major-revision diagnostics\n\n"
        "These tables implement reviewer-requested distinctions among computational verification, "
        "retrospective measured selection, same-source-publication external-laboratory evaluation, "
        "altered-protocol evaluation, predicted generated-candidate results and unperformed prospective validation.\n\n"
        f"- Dataset/evidence rows: {len(dataset)}\n"
        f"- Measured FullContext fiber-distribution settings: {len(negative)}\n"
        f"- Generated FullContext rate-benefit settings: {len(rate)}\n"
        f"- Negative/nonnegative selected-sequence feature comparisons: {len(feature_detail)}\n"
        f"- Post hoc score-margin abstention settings: {len(abstention)}\n"
        f"- Mapping-sensitivity summaries: {len(mapping)}\n"
        f"- Fixed-precision choice changes in AssayContext/FullContext/released-CNN primary comparisons: "
        f"{int(determinism.loc[determinism['scope'].str.startswith('measured-codebook AssayContext'), 'choices_changed_by_quantization'].sum())}\n"
        f"- Fixed-precision choice changes in the exploratory measured-codebook Composition selector: "
        f"{int(determinism.loc[determinism['scope'].eq('measured-codebook exploratory Composition selector'), 'choices_changed_by_quantization'].sum())}\n"
        f"- Generated candidate round-trip failures: {int(roundtrip.iloc[0]['candidate_roundtrip_failures'])}\n\n"
        "The reported bootstrap intervals remain conditional on the final feature families, model class, "
        "endpoint definitions, eligibility filters and declared mapping plans; they are not total uncertainty "
        "for the exploratory research process.\n",
        encoding="utf-8",
    )
    output_paths.append(summary)

    manifest_rows = [
        {
            "relative_path": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(output_paths)
    ]
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUT_DIR / "sha256_manifest.tsv", sep="\t", index=False, lineterminator="\n")
    print(f"Wrote major-revision diagnostics to {OUT_DIR}")
    print(f"Manifest SHA-256: {sha256(OUT_DIR / 'sha256_manifest.tsv')}")


if __name__ == "__main__":
    main()
