#!/usr/bin/env python3
"""Verify the current Paper 2 Bioinformatics choice-code manuscript bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
REFRAME = PAPER / "bioinformatics_reframe"
MAIN = PAPER / "main.tex"
SUPP = PAPER / "supplementary_codec_evidence.tex"
ROOT_BIB = PAPER / "references.bib"
AUDIT = REFRAME / "MANUSCRIPT_CONSISTENCY_AUDIT.md"

PUBLIC = REFRAME / "public_experimental_validation"
ASSAY = REFRAME / "assay_calibrated_selection"
CHOICE = REFRAME / "reversible_choice_codec"
SEQUENCE = REFRAME / "sequence_independence_audit"
RUNTIME = REFRAME / "runtime_benchmark"
SOTA = REFRAME / "sota_and_external_validation"
SOURCE_EXTERNAL = REFRAME / "source_external_independence_audit"
EXTERNAL_MAPPING = REFRAME / "external_mapping_sensitivity"
DIAGNOSTICS = REFRAME / "major_revision_diagnostics"
CHANNEL = REFRAME / "channel_error_boundary"

LATEX_LOG_NAMES = ("main.log", "supplementary_codec_evidence.log")
LATEX_ERROR_PATTERNS = (
    r"LaTeX Error",
    r"undefined references",
    r"Citation `[^']+' on page .* undefined",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def one(items: list[dict[str, str]], **wanted: object) -> dict[str, str]:
    matches = [
        row
        for row in items
        if all(str(row.get(key)) == str(value) for key, value in wanted.items())
    ]
    if len(matches) != 1:
        raise AssertionError(f"expected one row for {wanted}, found {len(matches)}")
    return matches[0]


def close(actual: str | float, expected: float, tolerance: float = 5e-12) -> None:
    if abs(float(actual) - expected) > tolerance:
        raise AssertionError(f"{actual} != {expected}")


def visible_abstract_word_count(tex: str) -> int:
    abstract = tex.split("\\begin{abstract}", 1)[1].split("\\end{abstract}", 1)[0]
    abstract = re.sub(
        r"\\(?:noindent|textbf|texttt)\{([^}]*)\}", r" \1 ", abstract
    )
    abstract = re.sub(r"\\[A-Za-z]+|\\.|[{}$^_]", " ", abstract)
    return len(re.findall(r"[A-Za-z0-9]+(?:[.-][A-Za-z0-9]+)*", abstract))


def verify_manifest(directory: Path) -> tuple[int, str]:
    manifest_path = directory / "sha256_manifest.tsv"
    manifest = tsv(manifest_path)
    for row in manifest:
        if "file" in row:
            path = directory / row["file"]
        elif "relative_path" in row:
            path = directory / row["relative_path"]
        else:
            raw = Path(row["path"])
            path = raw if raw.is_absolute() else ROOT / raw
        if not path.is_file():
            raise AssertionError(f"missing manifested output: {path}")
        if int(row["bytes"]) != path.stat().st_size:
            raise AssertionError(f"size mismatch: {path}")
        if row["sha256"] != sha256(path):
            raise AssertionError(f"hash mismatch: {path}")
    return len(manifest), sha256(manifest_path)


def assert_text_contains(text: str, phrases: list[str]) -> None:
    for phrase in phrases:
        if phrase not in text:
            raise AssertionError(f"required text missing: {phrase}")


def verify_references(tex: str) -> int:
    cited: set[str] = set()
    for group in re.findall(r"\\cite\{([^}]+)\}", tex):
        cited.update(key.strip() for key in group.split(","))
    bib_keys = set(re.findall(r"^@\w+\{([^,]+),", ROOT_BIB.read_text(encoding="utf-8"), re.M))
    missing = sorted(cited - bib_keys)
    if missing:
        raise AssertionError(f"missing BibTeX entries: {missing}")
    expected = {
        "cover1973enumerative",
        "adler1983slidingblock",
        "wang2019bioconstrained",
        "immink2019balancedhomopolymer",
        "erlich2017dnafountain",
        "gimpel2025pcrbias",
        "gimpel2023digitaltwin",
        "welzel2023dnaaeon",
        "zulkower2020dnachisel",
        "hammer2017rnablueprint",
        "ruan2026scone",
        "forney1992trellisshaping",
        "wu2022listccdm",
        "dou2024explorer",
        "press2020hedges",
        "zhao2024compositehedges",
        "volkel2023framed",
        "gimpel2026codecbenchmark",
        "zhang2026gungnir",
    }
    if cited != expected:
        raise AssertionError(f"unexpected cited-key set: {sorted(cited)}")
    return len(cited)


def audit_latex_logs(build_dir: Path) -> list[dict[str, str]]:
    """Classify optional compilation logs without weakening PDF requirements.

    Release archives intentionally omit transient LaTeX logs.  An absent log is
    therefore a documented SKIP, whereas a present clean log is PASS and a
    present log containing a declared error pattern is FAIL.
    """

    checks: list[dict[str, str]] = []
    for log_name in LATEX_LOG_NAMES:
        path = build_dir / log_name
        if not path.is_file():
            checks.append(
                {
                    "file": log_name,
                    "status": "SKIP",
                    "detail": "log absent from release archive; compiled-log checks skipped",
                }
            )
            continue
        log = path.read_text(encoding="utf-8", errors="replace")
        matched = [
            pattern
            for pattern in LATEX_ERROR_PATTERNS
            if re.search(pattern, log, flags=re.I)
        ]
        if matched:
            checks.append(
                {
                    "file": log_name,
                    "status": "FAIL",
                    "detail": "matched: " + ", ".join(matched),
                }
            )
        else:
            checks.append(
                {
                    "file": log_name,
                    "status": "PASS",
                    "detail": "present; no declared unresolved-reference or LaTeX-error pattern",
                }
            )
    return checks


def main() -> None:
    manuscript = MAIN.read_text(encoding="utf-8")
    supplement = SUPP.read_text(encoding="utf-8")
    combined = manuscript + "\n" + supplement

    words = visible_abstract_word_count(manuscript)
    if words > 150:
        raise AssertionError(f"visible structured abstract has {words} words")

    prediction = tsv(ASSAY / "prediction_metrics.tsv")
    expected_prediction = {
        ("repeated_nested_oof", "GCall"): 0.35405289911688576,
        ("repeated_nested_oof", "GCfix"): 0.3223807546903752,
        ("source_only_transfer", "GCall_to_GCfix"): 0.3185305405490628,
        ("source_only_transfer", "GCfix_to_GCall"): 0.2955429392937819,
    }
    for (evaluation, pool), value in expected_prediction.items():
        row = one(
            prediction,
            dataset="Gimpel2025_PCR",
            evaluation=evaluation,
            pool_or_direction=pool,
            model="P5_combined_context",
        )
        close(row["spearman"], value)

    fibers = tsv(CHOICE / "public_fiber_results.tsv")
    expected_fibers = {
        ("GCall_to_GCfix", "2"): (512, 0.0023490308967440937, 0.0019002820308026654, 0.002818400437795054, 0.25025491333659844),
        ("GCall_to_GCfix", "4"): (128, 0.003189670592436909, 0.0025525755796258426, 0.003851612187506525, 0.339812787813474),
        ("GCfix_to_GCall", "2"): (512, 0.0021354591246280225, 0.0016621447473999472, 0.002585022697895241, 0.27442918566916336),
        ("GCfix_to_GCall", "4"): (128, 0.002966662283132861, 0.002308550162257849, 0.0036090761360569137, 0.38124762264292344),
    }
    for (direction, selector_bits), expected in expected_fibers.items():
        row = one(
            fibers,
            direction=direction,
            selector_model="P5_combined_context",
            selector_bits=selector_bits,
        )
        if int(row["payloads"]) != expected[0] or int(row["roundtrip_failures"]) != 0:
            raise AssertionError(f"public fiber count/roundtrip mismatch: {direction}/r={selector_bits}")
        for key, value in zip(
            ["paired_mean_outcome_gain", "paired_gain_ci_2p5", "paired_gain_ci_97p5", "standardized_paired_gain"],
            expected[1:],
        ):
            close(row[key], value)

    comparisons = tsv(CHOICE / "public_fiber_model_comparisons.tsv")
    for direction in ["GCall_to_GCfix", "GCfix_to_GCall"]:
        for selector_bits in ["2", "4"]:
            row = one(
                comparisons,
                direction=direction,
                selector_bits=selector_bits,
                comparison="P5_combined_context_minus_P2_assay_context",
            )
            if float(row["difference_ci_2p5"]) <= 0:
                raise AssertionError(f"P5-minus-P2 interval is not positive: {direction}/r={selector_bits}")

    benchmark = tsv(SOTA / "fiber_benchmark_summary.tsv")
    expected_external = {
        ("GCall_to_external_Taq", "2"): (512, 0.002458436092035822, 0.23445735449266403),
        ("GCall_to_external_Taq", "4"): (128, 0.0040853863860975444, 0.38961715835028454),
        ("GCfix_to_external_Taq", "2"): (512, 0.002397661842549399, 0.22866140567706378),
        ("GCfix_to_external_Taq", "4"): (128, 0.0034325725731437894, 0.3273592863408233),
    }
    for (direction, selector_bits), expected in expected_external.items():
        row = one(
            benchmark,
            dataset="Gimpel2025_external_laboratory_Taq",
            direction=direction,
            selector_model="P5_combined_context",
            selector_bits=selector_bits,
        )
        if int(row["payloads"]) != expected[0] or int(row["failed_payloads"]) != 0:
            raise AssertionError(f"external fiber count/retention mismatch: {direction}/r={selector_bits}")
        close(row["paired_mean_outcome_gain"], expected[1])
        close(row["standardized_paired_gain"], expected[2])

    pairwise = tsv(SOTA / "fixed_model_pairwise_comparisons.tsv")
    p5_cnn = [
        row
        for row in pairwise
        if row["comparison"] == "P5_combined_context_minus_published_1dcnn"
        and row["selector_bits"] in {"2", "4"}
        and row["dataset"]
        in {
            "Gimpel2025_cross_pool_public_codebook",
            "Gimpel2025_external_laboratory_Taq",
        }
    ]
    if len(p5_cnn) != 8 or any(float(row["difference_ci_2p5"]) <= 0 for row in p5_cnn):
        raise AssertionError("P5-minus-published-1D-CNN identical-fiber intervals are not all positive")

    two_stage = tsv(SOTA / "two_stage_bootstrap_summary.tsv")
    p5_two_stage = [
        row
        for row in two_stage
        if row["estimand"] == "P5_combined_context_mean_gain"
        and row["selector_bits"] in {"2", "4"}
    ]
    p5_p2_two_stage = [
        row
        for row in two_stage
        if row["estimand"] == "P5_minus_P2_selected_outcome"
        and row["selector_bits"] in {"2", "4"}
    ]
    p5_cnn_two_stage = [
        row
        for row in two_stage
        if row["estimand"] == "P5_minus_released_1dcnn_selected_outcome"
        and row["selector_bits"] in {"2", "4"}
    ]
    positive_datasets = {
        "Gimpel2025_cross_pool_public_codebook",
        "Gimpel2025_external_laboratory_Taq",
    }
    q5_dataset = "Gimpel2025_external_laboratory_Q5_sensitivity"
    p5_positive = [row for row in p5_two_stage if row["target_dataset"] in positive_datasets]
    p5_q5 = [row for row in p5_two_stage if row["target_dataset"] == q5_dataset]
    p5_p2_positive = [
        row for row in p5_p2_two_stage if row["target_dataset"] in positive_datasets
    ]
    p5_p2_q5 = [row for row in p5_p2_two_stage if row["target_dataset"] == q5_dataset]
    p5_cnn_positive = [
        row for row in p5_cnn_two_stage if row["target_dataset"] in positive_datasets
    ]
    p5_cnn_q5 = [
        row for row in p5_cnn_two_stage if row["target_dataset"] == q5_dataset
    ]
    if len(p5_two_stage) != 12 or len(p5_positive) != 8 or any(
        int(row["replicates"]) != 2000 or float(row["two_stage_ci_2p5"]) <= 0
        for row in p5_positive
    ):
        raise AssertionError("FullContext cross-pool/KAPA two-stage intervals changed")
    if len(p5_q5) != 4 or any(float(row["two_stage_ci_97p5"]) >= 0 for row in p5_q5):
        raise AssertionError("Q5 FullContext workflow-shift boundary is not uniformly negative")
    if len(p5_p2_two_stage) != 12 or len(p5_p2_positive) != 8 or any(
        int(row["replicates"]) != 2000 or float(row["two_stage_ci_2p5"]) <= 0
        for row in p5_p2_positive
    ):
        raise AssertionError("FullContext-minus-AssayContext cross-pool/KAPA intervals changed")
    if len(p5_p2_q5) != 4 or any(
        float(row["two_stage_ci_97p5"]) >= 0 for row in p5_p2_q5
    ):
        raise AssertionError("Q5 FullContext-minus-AssayContext boundary changed")
    if len(p5_cnn_two_stage) != 12 or len(p5_cnn_positive) != 8 or any(
        int(row["replicates"]) != 2000 or float(row["two_stage_ci_2p5"]) <= 0
        for row in p5_cnn_positive
    ):
        raise AssertionError("FullContext-minus-CNN cross-pool/KAPA intervals changed")
    if len(p5_cnn_q5) != 4 or any(
        float(row["two_stage_ci_97p5"]) >= 0 for row in p5_cnn_q5
    ):
        raise AssertionError("Q5 FullContext-minus-CNN boundary changed")

    randomization = tsv(SOTA / "fixed_model_randomization_tests.tsv")
    p5_randomization = [
        row
        for row in randomization
        if row["selector_model"] == "P5_combined_context"
        and row["selector_bits"] in {"2", "4"}
    ]
    p5_randomization_positive = [
        row for row in p5_randomization if row["dataset"] in positive_datasets
    ]
    if len(p5_randomization) != 12 or len(p5_randomization_positive) != 8 or any(
        int(row["randomization_replicates"]) != 100000
        or int(row["one_sided_exceedances_b"]) != 0
        or row["monte_carlo_p_formula"] != "(b+1)/(B+1)"
        or abs(
            float(row["randomization_p_one_sided_positive"])
            - float(row["minimum_attainable_monte_carlo_p"])
        )
        > 1e-12
        for row in p5_randomization_positive
    ):
        raise AssertionError("FullContext cross-pool/KAPA random-choice tests changed")
    if sum(1 for _ in (SOTA / "external_locked_library_features.tsv").open(encoding="utf-8")) - 1 != 2048:
        raise AssertionError("locked external library is not 2,048 sequences")
    if sum(1 for _ in (SOTA / "external_q5_locked_library_features.tsv").open(encoding="utf-8")) - 1 != 1024:
        raise AssertionError("locked external Q5 sensitivity library is not 1,024 sequences")

    event_pairwise = tsv(SOTA / "endpoint_matched_low_efficiency_pairwise.tsv")
    event_primary = [row for row in event_pairwise if row["dataset"] in positive_datasets]
    if len(event_primary) != 8 or any(
        float(row["event_risk_difference_ci_2p5"]) > 0
        or float(row["event_risk_difference_ci_97p5"]) < 0
        for row in event_primary
    ):
        raise AssertionError("endpoint-matched FullContext-minus-CNN intervals no longer include zero")

    sensitivity = tsv(CHOICE / "public_seed_sensitivity_summary.tsv")
    for direction in ["GCall_to_GCfix", "GCfix_to_GCall"]:
        for selector_bits in ["2", "4"]:
            for model in ["P5_combined_context", "P5_minus_P2_selected_outcome"]:
                row = one(
                    sensitivity,
                    direction=direction,
                    selector_bits=selector_bits,
                    selector_model=model,
                )
                close(row["positive_seed_fraction"], 1.0)
                if int(row["seeds"]) != 32:
                    raise AssertionError("mapping-sensitivity seed count changed")

    external_mapping = tsv(
        EXTERNAL_MAPPING / "external_kapa_mapping_sensitivity_summary.tsv"
    )
    if len(external_mapping) != 4:
        raise AssertionError("external KAPA mapping summary must contain four source-width rows")
    external_metrics = (
        "fullcontext_gain",
        "fullcontext_minus_assaycontext",
        "fullcontext_minus_released_cnn",
    )
    for row in external_mapping:
        if int(row["mappings"]) != 32:
            raise AssertionError("external KAPA mapping count changed")
        if row["mapping_role"] != "post_hoc_exploratory_algorithmic_sensitivity_not_replicates":
            raise AssertionError("external KAPA mapping evidence role changed")
        for metric in external_metrics:
            if int(row[f"{metric}_positive_mappings"]) != 32:
                raise AssertionError(f"external KAPA {metric} is not positive in all mappings")
            if float(row[f"{metric}_min"]) <= 0:
                raise AssertionError(f"external KAPA {metric} minimum is not positive")
    primary_mapping_checks = tsv(
        EXTERNAL_MAPPING / "primary_mapping_reproduction_check.tsv"
    )
    if len(primary_mapping_checks) != 12 or any(
        row["status"] != "PASS"
        or abs(float(row["difference"])) > float(row["tolerance"])
        for row in primary_mapping_checks
    ):
        raise AssertionError("analysis-plan-locked external mapping reproduction check failed")
    if sum(
        1
        for _ in (EXTERNAL_MAPPING / "external_kapa_mapping_inputs.tsv").open(
            encoding="utf-8"
        )
    ) - 1 != 2053:
        raise AssertionError("external KAPA mapping input table is not 2,053 sequences")

    base = one(tsv(CHOICE / "base_language_summary.tsv"), language="exact_108nt_gc49_59_hp3")
    if (
        int(base["length_nt"]) != 108
        or int(base["fixed_payload_bits"]) != 213
        or int(base["rank_dispersion_modulus"]) != int(base["exact_total"])
        or base["rank_dispersion_contract"]
        != "affine permutation on complete exact language; gcd(multiplier,N)=1; inverse rejects q>=2^K"
        or int(base["completion_cache_entries"]) != 41101
        or base["top_down_bottom_up_match"] != "True"
        or int(base["reduced_all_rank_checks"]) != 45128
        or int(base["reduced_rank_unrank_failures"]) != 0
    ):
        raise AssertionError("exact base-language invariant changed")
    close(base["fixed_payload_rate_bits_per_nt"], 213 / 108)
    close(base["legacy_lexicographic_prefix_first_base_A_fraction"], 0.40870571610114614)
    close(base["legacy_lexicographic_prefix_first_base_C_fraction"], 0.40870571610114614)
    close(base["legacy_lexicographic_prefix_first_base_G_fraction"], 0.18258856779770768)
    close(base["legacy_lexicographic_prefix_first_base_T_fraction"], 0.0)

    roundtrip = one(tsv(CHOICE / "generated_roundtrip_summary.tsv"), language="exact_108nt_gc49_59_hp3")
    if (
        int(roundtrip["candidate_rows"]) != 65024
        or int(roundtrip["candidate_exact_roundtrips"]) != 65024
        or int(roundtrip["candidate_roundtrip_failures"]) != 0
        or roundtrip["exact_language_out_of_domain_rank_rejected"] != "True"
        or roundtrip["payload_domain_out_of_domain_sequence_rejected"] != "True"
        or roundtrip["invalid_constraint_sequence_rejected"] != "True"
        or roundtrip["noncanonical_valid_candidate_accepted_by_payload_decoder"] != "True"
        or roundtrip["canonical_emitted_candidate_accepted_by_optional_verifier"] != "True"
        or roundtrip["noncanonical_valid_candidate_rejected_by_optional_verifier"] != "True"
        or roundtrip["wrong_length_sequence_rejected"] != "True"
        or roundtrip["ambiguous_base_sequence_rejected"] != "True"
        or roundtrip["exact_tie_selects_smallest_choice_index"] != "True"
        or roundtrip["nonfinite_score_rejected_by_encoder_verifier"] != "True"
        or roundtrip.get("rte_half_integer_vectors_pass") != "True"
        or roundtrip.get("negative_score_vectors_pass") != "True"
        or roundtrip.get("positive_int64_overflow_rejected") != "True"
        or roundtrip["score_replacement_preserves_payload_decoding"] != "True"
        or int(roundtrip["score_decimal_places"]) != 12
    ):
        raise AssertionError("generated roundtrip invariant changed")

    positional = [
        row
        for row in tsv(CHOICE / "generated_positional_nucleotide_frequencies.tsv")
        if row["population"] == "generated_unique_candidates"
    ]
    if len(positional) != 108:
        raise AssertionError("generated positional-frequency audit does not cover 108 positions")
    first = one(positional, position_1_based="1")
    for nucleotide in "ACGT":
        value = float(first[f"frequency_{nucleotide}"])
        if not 0.245 <= value <= 0.255:
            raise AssertionError(f"generated first-base {nucleotide} remains biased: {value}")
    all_position_frequencies = [
        float(row[f"frequency_{nucleotide}"])
        for row in positional
        for nucleotide in "ACGT"
    ]
    if min(all_position_frequencies) < 0.23 or max(all_position_frequencies) > 0.27:
        raise AssertionError("generated positional frequencies exceed the declared dispersion audit range")

    generated_mapping = tsv(CHOICE / "generated_mapping_namespace_sensitivity_summary.tsv")
    if len(generated_mapping) != 4 or any(
        int(row["mapping_namespaces"]) != 8
        or float(row["own_gain_positive_namespace_fraction"]) != 1.0
        or float(row["cross_gain_positive_namespace_fraction"]) != 1.0
        or int(row["roundtrip_failures_total"]) != 0
        for row in generated_mapping
    ):
        raise AssertionError("generated full-language namespace sensitivity changed")

    diagnostic_dataset = tsv(DIAGNOSTICS / "dataset_partition_and_evidence_table.tsv")
    diagnostic_hierarchy = tsv(DIAGNOSTICS / "evidence_hierarchy.tsv")
    diagnostic_negative = tsv(DIAGNOSTICS / "measured_negative_fiber_diagnostics.tsv")
    diagnostic_rate = tsv(DIAGNOSTICS / "generated_rate_benefit.tsv")
    diagnostic_feature = tsv(DIAGNOSTICS / "negative_fiber_feature_stratification.tsv")
    diagnostic_abstention = tsv(DIAGNOSTICS / "score_margin_abstention_sensitivity.tsv")
    diagnostic_mapping = tsv(DIAGNOSTICS / "mapping_sensitivity_diagnostics.tsv")
    diagnostic_score = tsv(DIAGNOSTICS / "score_determinism_audit.tsv")
    diagnostic_failures = tsv(DIAGNOSTICS / "decoder_failure_behavior.tsv")
    diagnostic_leakage = tsv(DIAGNOSTICS / "leakage_scope.tsv")
    if len(diagnostic_dataset) != 5 or len(diagnostic_hierarchy) != 6:
        raise AssertionError("dataset/evidence hierarchy diagnostics are incomplete")
    if len(diagnostic_negative) != 12:
        raise AssertionError("negative-fiber diagnostics must contain 12 P5 settings")
    matched_negative = [row for row in diagnostic_negative if "Q5" not in row["dataset"]]
    q5_negative = [row for row in diagnostic_negative if "Q5" in row["dataset"]]
    if len(matched_negative) != 8 or not all(
        0.17 <= float(row["negative_gain_fraction"]) <= 0.31
        for row in matched_negative
    ):
        raise AssertionError("matched-condition negative-fiber ranges changed")
    if len(q5_negative) != 4 or not all(
        float(row["mean_gain"]) < 0
        and 0.47 <= float(row["negative_gain_fraction"]) <= 0.72
        for row in q5_negative
    ):
        raise AssertionError("Q5 negative-fiber boundary changed")
    if len(diagnostic_rate) != 14 or {
        (row["selector_source_pool"], int(row["choice_bits"]))
        for row in diagnostic_rate
    } != {(pool, bits) for pool in ("GCall", "GCfix") for bits in range(7)}:
        raise AssertionError("generated r=0..6 rate-benefit table is incomplete")
    if len(diagnostic_feature) != 240 or len(diagnostic_abstention) != 36:
        raise AssertionError("negative-fiber feature/abstention diagnostics are incomplete")
    q5_abstention = [
        row
        for row in diagnostic_abstention
        if "Q5" in row["dataset"] and float(row["nominal_score_margin_retention"]) == 0.75
    ]
    if len(q5_abstention) != 4 or any(
        float(row["mean_measured_gain"]) >= 0 for row in q5_abstention
    ):
        raise AssertionError("score-margin abstention unexpectedly closes the Q5 boundary")
    if len(diagnostic_mapping) != 12:
        raise AssertionError("mapping diagnostics must contain 12 scope/source/width rows")
    primary_score_rows = [
        row
        for row in diagnostic_score
        if row["scope"].startswith(
            "measured-codebook AssayContext, FullContext and released-CNN"
        )
    ]
    if len(primary_score_rows) != 1 or int(
        primary_score_rows[0]["choices_changed_by_quantization"]
    ) != 0:
        raise AssertionError("fixed-precision keys changed a primary measured-codebook choice")
    if not any(
        row["input_condition"] == "NaN or infinite score"
        and row["verified"] == "True"
        for row in diagnostic_failures
    ):
        raise AssertionError("non-finite score failure path is not recorded")
    if not any(row["status"] == "not performed" for row in diagnostic_leakage):
        raise AssertionError("unperformed leakage controls are not disclosed")

    overlap = tsv(CHOICE / "generated_feature_overlap.tsv")
    eligible_active = [
        row
        for row in overlap
        if row["source_reference"] == "codec_eligible_training_subset"
        and row["p5_active_feature"] == "True"
    ]
    for source_pool in ["GCall", "GCfix"]:
        rows = [row for row in eligible_active if row["source_pool"] == source_pool]
        if not rows or max(
            abs(float(row["generated_minus_source_standardized_mean_difference"]))
            for row in rows
        ) > 0.04:
            raise AssertionError(f"eligible-subset feature SMD changed for {source_pool}")

    thresholds = tsv(CHOICE / "threshold_rejection_baseline.tsv")
    expected_thresholds = {
        ("GCall", "2"): 0.341796875,
        ("GCall", "4"): 0.009765625,
        ("GCfix", "2"): 0.314453125,
        ("GCfix", "4"): 0.005859375,
    }
    for (pool, bits), value in expected_thresholds.items():
        row = one(
            thresholds,
            selector_source_pool=pool,
            selector_model="P5_combined_context",
            selector_bits=bits,
        )
        close(row["failure_fraction"], value)
        close(row["choice_codec_failure_fraction"], 0.0)

    sequence = one(
        tsv(SEQUENCE / "cross_pool_independence_summary.tsv"),
        audit_scope="GCall_variable_regions_vs_GCfix_variable_regions",
    )
    if (
        int(sequence["all_cross_pool_pairs"]) != 143904012
        or int(sequence["minimum_cross_pool_hamming_distance"]) != 54
        or float(sequence["maximum_cross_pool_identity_fraction"]) != 0.5
        or any(
            int(sequence[key]) != 0
            for key in [
                "cross_pool_exact_duplicates",
                "cross_pool_exact_reverse_complements",
                "pairs_identity_at_least_90_percent",
                "pairs_identity_at_least_80_percent",
                "pairs_identity_at_least_70_percent",
            ]
        )
    ):
        raise AssertionError("cross-pool sequence-independence invariant changed")

    source_external_rows = tsv(
        SOURCE_EXTERNAL / "source_external_independence_summary.tsv"
    )
    expected_source_external = {
        "GCall_variable_regions_vs_external_locked_codebook_2048": (
            24571904,
            55,
            0.49074074074074076,
        ),
        "GCfix_variable_regions_vs_external_locked_codebook_2048": (
            24563712,
            56,
            0.48148148148148145,
        ),
        "GCall_variable_regions_vs_external_codec_eligible_2053": (
            24631894,
            55,
            0.49074074074074076,
        ),
        "GCfix_variable_regions_vs_external_codec_eligible_2053": (
            24623682,
            56,
            0.48148148148148145,
        ),
    }
    for scope, expected in expected_source_external.items():
        row = one(source_external_rows, audit_scope=scope)
        if (
            int(row["all_oriented_pairs"]) != expected[0]
            or int(row["minimum_oriented_hamming_distance_nt"]) != expected[1]
            or abs(float(row["maximum_oriented_identity_fraction"]) - expected[2])
            > 1e-15
            or any(
                int(row[key]) != 0
                for key in [
                    "exact_duplicates",
                    "exact_reverse_complements",
                    "pairs_identity_at_least_90_percent",
                    "pairs_identity_at_least_80_percent",
                    "pairs_identity_at_least_70_percent",
                ]
            )
        ):
            raise AssertionError(f"source-external independence changed: {scope}")

    runtime_rows = tsv(RUNTIME / "runtime_summary.tsv")
    runtime_env = json.loads((RUNTIME / "environment.json").read_text(encoding="utf-8"))
    if int(runtime_env["base_count_cache_entries"]) != 41101:
        raise AssertionError("runtime cache-entry count changed")
    if not all(row["all_roundtrips_passed"] == "True" for row in runtime_rows):
        raise AssertionError("a timed roundtrip failed")
    runtime_by_key = {
        (row["operation"], int(row["selector_bits"])): row for row in runtime_rows
    }
    required_runtime_keys = {
        ("base_unrank", 0),
        ("base_rank", 0),
        ("P5_choice_encode", 0),
        ("P5_choice_encode", 2),
        ("P5_choice_encode", 4),
        ("choice_decode", 0),
        ("choice_decode", 2),
        ("choice_decode", 4),
    }
    if set(runtime_by_key) != required_runtime_keys:
        raise AssertionError("runtime operation set changed")

    def runtime_ms(operation: str, selector_bits: int) -> float:
        return float(
            runtime_by_key[(operation, selector_bits)][
                "median_milliseconds_per_payload_or_item"
            ]
        )

    rank_ms = runtime_ms("base_rank", 0)
    unrank_ms = runtime_ms("base_unrank", 0)
    encode_ms = [runtime_ms("P5_choice_encode", r) for r in (0, 2, 4)]
    decode_ms = [runtime_ms("choice_decode", r) for r in (0, 2, 4)]
    row_peak_rss_mib = max(float(row["peak_process_rss_mib"]) for row in runtime_rows)
    peak_rss_mib = float(runtime_env["peak_process_rss_mib"])
    if peak_rss_mib < row_peak_rss_mib:
        raise AssertionError("environment peak RSS is below a timed-row peak")
    assert_text_contains(
        supplement,
        [
            f"Base unrank & 0 & 1 & 4,096 & {unrank_ms:.3f}",
            f"Base rank & 0 & 1 & 4,096 & {rank_ms:.3f}",
            f"FullContext choice encode & 4 & 16 & 128 & {encode_ms[2]:.3f}",
            f"Payload decode & 4 & 16 & 128 & {decode_ms[2]:.3f}",
        ],
    )

    channel_contract = one(
        tsv(CHANNEL / "channel_error_audit_contract.tsv"),
        audit="single_edit_channel_boundary",
    )
    expected_channel_contract = {
        "language": "exact_108nt_gc49_59_hp3",
        "selector_bits": "2",
        "candidate_width": "4",
        "sources": "2",
        "messages_per_source": "32",
        "emitted_codewords": "64",
        "total_edit_operations": "55552",
        "crc_name": "CRC-16/CCITT-FALSE",
        "crc_check_vector_123456789": "0x29b1",
        "message_bits": "195",
        "crc_bits": "16",
        "choice_bits": "2",
        "wrong_payload_acceptances_without_checksum": "12093",
        "wrong_payload_acceptances_after_crc16": "0",
        "wrong_payload_acceptances_after_crc16_and_optional_canonical_verifier": "0",
    }
    for key, expected in expected_channel_contract.items():
        if channel_contract[key] != expected:
            raise AssertionError(f"channel-error contract changed for {key}")
    if channel_contract["scope"] != (
        "error-detection composition audit only; no correction, dropout, "
        "clustering, consensus or file recovery"
    ):
        raise AssertionError("channel-error inference boundary changed")

    edit_summary = tsv(CHANNEL / "single_edit_summary.tsv")
    if len(edit_summary) != 6:
        raise AssertionError("single-edit summary must contain six source/edit rows")
    substitution_wrong = {
        row["selector_source_pool"]: int(row["accepted_wrong_payload"])
        for row in edit_summary
        if row["edit_type"] == "substitution"
    }
    if substitution_wrong != {"GCall": 6046, "GCfix": 6047}:
        raise AssertionError("single-substitution silent-misdecode counts changed")
    for row in edit_summary:
        if row["edit_type"] in {"insertion", "deletion"} and int(
            row["rejected_length"]
        ) != int(row["edit_operations"]):
            raise AssertionError("fixed-length insertion/deletion rejection changed")
        if int(row["crc16_residual_wrong_payload"]) != 0 or int(
            row["combined_crc16_and_canonical_residual_wrong_payload"]
        ) != 0:
            raise AssertionError("CRC-composition finite-sample residual changed")
    crc_vector = one(
        tsv(CHANNEL / "crc16_reference_vectors.tsv"),
        algorithm="CRC-16/CCITT-FALSE",
        input_ascii="123456789",
    )
    if crc_vector["expected_hex"] != "0x29b1" or crc_vector["passed"] != "True":
        raise AssertionError("CRC-16/CCITT-FALSE check vector failed")

    required_boundaries = [
        "gains were predicted",
        "none was synthesized or assayed",
        "were not independent laboratories",
        "not independent third-party replication",
        "not physical storage rates",
        "no NUPACK calculation",
        "but not synthesis errors",
        "negative domain boundary",
        "already-public measurements",
        "exclude technical-replicate, batch, normalization and assay-measurement uncertainty",
        "single-edit/CRC audit",
        "secondary exploratory Q5 workflow-shift sensitivity",
        "Across 32 additional post hoc outcome-blind mappings",
        "optional canonical verifier",
        "modulo-$N$",
        "conditional on final feature families",
        "score-margin abstention",
        "round-to-nearest with ties to even",
        "Hamming-separated PCR pools",
    ]
    assert_text_contains(combined, required_boundaries)
    assert_text_contains(
        supplement,
        [
            "9,955 Q5-complete rows",
            "2,034 Q5-complete sequences and 14 Q5-missing sequences",
            "yielding 2,040 eligible rows",
            "256 complete fibers at $r=2$ and 64 at $r=4$",
            "threshold was 0.9314214756",
            "199/9,955",
        ],
    )
    forbidden_positive_claims = [
        r"NUPACK (?:confirmed|validated)",
        r"recovered (?:the|a) 4,096-byte",
        r"full 110-nt weighted-language capacity is",
        r"public sequences were emitted by",
        r"generated codewords (?:showed|demonstrated) measured PCR",
        r"achieved physical storage density",
    ]
    for pattern in forbidden_positive_claims:
        if re.search(pattern, combined, flags=re.I):
            raise AssertionError(f"forbidden positive claim found: {pattern}")
    forbidden_internal_text = [
        "AUTHOR_INPUT_NEEDED",
        "unpublished sister",
        "NUPACK bridge",
        "material-return",
        "requested by the target journal",
    ]
    assert not any(phrase in combined for phrase in forbidden_internal_text), forbidden_internal_text
    if "\\author{}" not in manuscript or "\\author{}" not in supplement:
        raise AssertionError("anonymous author blocks are not empty")
    for heading in (
        "\\section{System and methods}",
        "\\section{Algorithm}",
        "\\section{Implementation}",
        "\\section{Results}",
        "\\section{Discussion}",
    ):
        if heading not in manuscript:
            raise AssertionError(f"required Bioinformatics section missing: {heading}")
    if "sequence-distinct" in manuscript:
        raise AssertionError("ambiguous sequence-distinct wording remains in main manuscript")

    citation_count = verify_references(combined)

    required_files = [
        PAPER / "figures" / "assay_calibrated_selection.pdf",
        PAPER / "figures" / "reversible_choice_codec.pdf",
        PAPER / "figures" / "paired_two_stage_sensitivity.pdf",
        PAPER / "figures" / "paired_two_stage_sensitivity_source_data.tsv",
        PAPER / "tests" / "test_paper2_deterministic_selection.py",
        PAPER / "tests" / "test_paper2_grouped_bootstrap.py",
        PAPER / "tests" / "test_paper2_channel_error_boundary.py",
        PAPER / "tests" / "data" / "paper2_score_quantization_vectors.tsv",
        CHANNEL / "channel_error_audit_contract.tsv",
        CHANNEL / "single_edit_summary.tsv",
        CHANNEL / "crc16_reference_vectors.tsv",
        PAPER / "build" / "main.pdf",
        PAPER / "build" / "supplementary_codec_evidence.pdf",
        PAPER / "oup_preflight" / "build" / "main_oup_preflight.pdf",
        REFRAME / "PDF_VISUAL_QC_20260722.md",
    ]
    missing = [str(path) for path in required_files if not path.is_file()]
    if missing:
        raise AssertionError(f"missing manuscript artifacts: {missing}")

    figure_two_source = tsv(PAPER / "figures" / "reversible_choice_codec_source_data.tsv")
    panel_d = [row for row in figure_two_source if row["panel"] == "d"]
    matched_panel_d = [
        row for row in panel_d if row["series"] in {"Cross-pool", "External KAPA"}
    ]
    q5_panel_d = [row for row in panel_d if row["series"] == "Q5 altered protocol"]
    if len(matched_panel_d) != 8 or any(float(row["low"]) <= 0 for row in matched_panel_d):
        raise AssertionError("Figure 2d matched-assay intervals are incomplete or not positive")
    if len(q5_panel_d) != 4 or any(float(row["high"]) >= 0 for row in q5_panel_d):
        raise AssertionError("Figure 2d Q5 intervals are incomplete or not negative")
    assert_text_contains(
        manuscript,
        [
            "recovers the payload from any candidate",
            "measured selected-minus-fiber-mean gains",
            "absolute FullContext gains were modest",
            "does not guarantee detection under arbitrary corruptions",
            "negative secondary exploratory Q5 gains on a common zero axis",
        ],
    )
    if "recover any fiber member" in manuscript:
        raise AssertionError("abstract still describes recovery of a fiber member rather than its payload")

    latex_log_checks = audit_latex_logs(PAPER / "build")
    for check in latex_log_checks:
        print(
            f"{check['status']}_LATEX_LOG {check['file']}: {check['detail']}"
        )
    failed_latex_logs = [
        check for check in latex_log_checks if check["status"] == "FAIL"
    ]
    if failed_latex_logs:
        details = "; ".join(
            f"{check['file']}: {check['detail']}" for check in failed_latex_logs
        )
        raise AssertionError(f"LaTeX compilation-log audit failed: {details}")

    manifest_results = {
        name: verify_manifest(path)
        for name, path in {
            "public_experimental_validation": PUBLIC,
            "assay_calibrated_selection": ASSAY,
            "reversible_choice_codec": CHOICE,
            "sequence_independence_audit": SEQUENCE,
            "runtime_benchmark": RUNTIME,
            "sota_and_external_validation": SOTA,
            "source_external_independence_audit": SOURCE_EXTERNAL,
            "external_mapping_sensitivity": EXTERNAL_MAPPING,
            "major_revision_diagnostics": DIAGNOSTICS,
            "channel_error_boundary": CHANNEL,
        }.items()
    }
    local_runtime = RUNTIME / "latest_local"
    if (local_runtime / "sha256_manifest.tsv").is_file():
        manifest_results["runtime_benchmark_latest_local_nonfrozen"] = verify_manifest(
            local_runtime
        )

    report_lines = [
        "# Paper 2 Manuscript Consistency Audit",
        "",
        "Date: 2026-07-22 (Asia/Shanghai)",
        "",
        "Status: `PASS_SUBMISSION_TECHNICAL_OUTPUT_AND_BOUNDARY_CONSISTENCY`",
        "",
        f"- Visible structured abstract count: {words} words, including headings and placeholders (maximum recommended: 150).",
        "- FullContext within-pool and source-only transfer Spearman values match the frozen assay table.",
        "- Public FullContext choice-fiber gains and exact round trips match both transfer directions at `r=2` and `r=4`.",
        "- Analysis-plan-locked external KAPA gains match both frozen source models and widths on 2,048 outcome-blind sequences drawn from already-public measurements.",
        "- Corrected grouped-source plus target-fiber bootstrap intervals are positive in all eight cross-pool/KAPA FullContext settings; resampled duplicate source rows retain their original fold, equivalent to integer weighting.",
        "- FullContext-minus-released-1D-CNN continuous-efficiency intervals are positive in all eight cross-pool/KAPA comparisons after the corrected two-stage propagation.",
        "- Endpoint-matched FullContext-minus-CNN low-efficiency intervals include zero in all eight cross-pool/KAPA comparisons; no cross-endpoint dominance is asserted.",
        "- The Q5 workflow-shift FullContext, FullContext-minus-AssayContext and FullContext-minus-CNN two-stage intervals are negative in all four source-width settings.",
        "- Figure 2d places eight positive matched-assay and four negative exploratory-Q5 FullContext intervals on one zero-effect axis; its source-data rows are checked directly.",
        "- Cross-pool/KAPA FullContext random-choice tests retain 100,000 draws, zero exceedances and the declared minimum attainable Monte Carlo P formula.",
        "- FullContext and FullContext-minus-AssayContext gains are positive in all 32 outcome-blind public-codebook mappings.",
        "- FullContext gain, FullContext-minus-AssayContext and FullContext-minus-CNN are positive in all 32 post hoc external KAPA mappings for both source models and widths; mappings are not treated as replicates.",
        "- The data-partition/evidence hierarchy, r=0..6 rate-benefit, negative-fiber distributions, feature stratification, score-margin abstention and fixed-precision selector audits are complete and manifest-locked.",
        "- Fixed-precision score keys changed no AssayContext, FullContext or released-CNN primary choice; exploratory Composition differences are disclosed.",
        "- Exact 108-nt total, modulo-N rank dispersion, 213-bit domain, 45,128 reduced all-rank checks and 65,024 generated candidate round trips are internally consistent.",
        "- Generated first-base and all-position composition audits exclude the former lexicographic-prefix bias; eight generated mapping namespaces remain positive with zero round-trip failures.",
        "- Exhaustive 143,904,012-pair fixed-length Hamming audit retains zero exact duplicates or reverse complements, minimum distance 54/108 and no pair at 70% identity; no edit-distance, genealogy or cluster-independence claim is made.",
        "- Exhaustive source-external audits retain zero duplicates, reverse complements or pairs at 70% identity across both the 2,048 lock and all 2,053 eligible sequences.",
        "- Threshold-rejection failure fractions and zero choice-codec failures match the frozen table.",
        "- Runtime prose and SI table match the current fixed-input, machine-specific benchmark; all timed round trips passed.",
        "- The deterministic single-edit audit covers 55,552 substitutions, insertions and deletions over 64 emitted codewords; 12,093 wrong-payload acceptances without a checksum fell to zero after CRC-16/CCITT-FALSE in this finite composition audit, which is not a channel or recovery claim.",
        "- Reported intervals condition on already-public sequence-level measurements and do not propagate technical-replicate, batch, normalization or assay-measurement uncertainty.",
        f"- All {citation_count} active citation keys exist in the paper-local BibTeX file and match the separately verified reference set.",
        "- Required public/generated, computational/material and observational/causal boundaries are present; forbidden positive-claim patterns are absent.",
        "- LaTeX compilation-log audit: "
        + "; ".join(
            f"{check['file']}={check['status']} ({check['detail']})"
            for check in latex_log_checks
        )
        + ". Compiled PDFs remain required even when release-excluded logs are skipped.",
        "- The compiled manuscript, seven-page OUP preflight and Supplementary pages are covered by the rendered-page audit in `PDF_VISUAL_QC_20260722.md`.",
        "",
        "## Verified manifests",
        "",
    ]
    for name, (count, digest) in manifest_results.items():
        report_lines.append(f"- `{name}`: {count} entries; manifest SHA-256 `{digest}`.")
    report_lines.extend(
        [
            "",
            "## Current artifact hashes",
            "",
            f"- `main.tex`: `{sha256(MAIN)}`",
            f"- `supplementary_codec_evidence.tex`: `{sha256(SUPP)}`",
            f"- `build/main.pdf`: `{sha256(PAPER / 'build' / 'main.pdf')}`",
            f"- `build/supplementary_codec_evidence.pdf`: `{sha256(PAPER / 'build' / 'supplementary_codec_evidence.pdf')}`",
            f"- `figures/assay_calibrated_selection.pdf`: `{sha256(PAPER / 'figures' / 'assay_calibrated_selection.pdf')}`",
            f"- `figures/reversible_choice_codec.pdf`: `{sha256(PAPER / 'figures' / 'reversible_choice_codec.pdf')}`",
            f"- `figures/paired_two_stage_sensitivity.pdf`: `{sha256(PAPER / 'figures' / 'paired_two_stage_sensitivity.pdf')}`",
            "",
            "## Submission boundary",
            "",
            "Scientific and technical consistency is a narrow PASS, not an acceptance guarantee. The public repository and BSD-3-Clause licence are machine-verifiable. A persistent archive URL or identifier for the manuscript-matched software and test-data snapshot remains outstanding because the author deferred that external publication action; a DOI is one possible identifier, not the only route. Journal upload also requires author-controlled identity/contact/funding/conflict/CRediT fields and author verification of any policy-required declarations and final wording. The generated codec still has no prospective emitted-codeword wet-lab validation.",
            "",
        ]
    )
    AUDIT.write_text("\n".join(report_lines), encoding="utf-8")
    print("PASS_SUBMISSION_TECHNICAL_OUTPUT_AND_BOUNDARY_CONSISTENCY")
    print(AUDIT)


if __name__ == "__main__":
    main()
