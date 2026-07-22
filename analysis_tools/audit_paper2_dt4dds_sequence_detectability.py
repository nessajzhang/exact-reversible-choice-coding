#!/usr/bin/env python3
"""External DT4DDS sequence-risk versus detectability calibration for Paper 2."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import random
import statistics
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import hairpin_risk_features as hrf  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "external_data" / "paper2_dt4dds_sequence_recovery"
DETAIL = ROOT / "data" / "paper2_dt4dds_sequence_detectability.csv"
SUMMARY = ROOT / "data" / "paper2_dt4dds_sequence_detectability_summary.csv"
NOTE = (
    ROOT
    / "papers"
    / "paper2_thermodynamic_risk_coding"
    / "dt4dds_sequence_detectability_calibration_20260713.md"
)

SAMPLES = ["0a", "0b", "2d", "4d", "7d"]
POOLS = {
    "Genscript_GCall": {
        "design_blob": "cbe3467327585f2dc7483d869d8608e65093c041",
        "scaf_blobs": [
            "4723f905c9232a556f040ba76f4a5065a1b9d25d",
            "69239e1a6d7c5236f121d6be18082130da3bed67",
            "45bb6bf5f97525f9eb367dd1d134bdadfac31666",
            "db59ca0c47ec58f022f17fa97705d24170bb733e",
            "595f11341ba1f9528126f881c07e6be48be6ba58",
        ],
        "reference_count": 12_472,
        "missing_counts": [110, 136, 113, 100, 100],
    },
    "Twist_GCall": {
        "design_blob": "a87f05318c6bbab78d255984dbc1ee03f6f783a8",
        "scaf_blobs": [
            "ce538ea62b6275d9ac068c95f48f9f8ac5e5dc5e",
            "f7733b5fd030863bdd6c1faf178a97ef8fbff8c7",
            "c1665a638d2b9f0de09098672e8986c251e0b5a6",
            "ce31a9206e756a049239b377564fa99a6beb2ddb",
            "2808731cc1bf2c03f0264cdeae001377e43b7647",
        ],
        "reference_count": 12_000,
        "missing_counts": [9, 12, 16, 17, 11],
    },
}

BOUNDARY = (
    "external observational calibration of post-workflow reference detectability; "
    "not isolated molecular aging, not a randomized sequence intervention, not direct "
    "hairpin thermodynamics or NUPACK validation, not Paper 2 emitted-codec validation, "
    "not sequence or file recovery, and not storage-density evidence"
)


def git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    if len(records) != len(set(records)) or any(set(seq) - set("ACGT") for seq in records.values()):
        raise ValueError("invalid FASTA records")
    return records


def parse_scafstats(path: Path) -> dict[str, int]:
    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().lstrip("#").strip().split("\t")
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
        return {row["name"]: int(row["assignedReads"]) for row in reader}


def rank_average(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop]] = (start + stop - 1) / 2.0 + 1.0
        start = stop
    return ranks


def spearman_from_ranks(left_rank: np.ndarray, right_rank: np.ndarray) -> float:
    return float(np.corrcoef(left_rank, right_rank)[0, 1])


def auc_from_ranks(feature_rank: np.ndarray, labels: np.ndarray) -> float:
    positives = int(labels.sum())
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return math.nan
    rank_sum = float(feature_rank[labels.astype(bool)].sum())
    return (rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def permutation_continuous(
    feature_rank: np.ndarray, outcome_rank: np.ndarray, permutations: int, seed: int
) -> tuple[float, float]:
    left = feature_rank - feature_rank.mean()
    right = outcome_rank - outcome_rank.mean()
    denominator = math.sqrt(float(np.dot(left, left)) * float(np.dot(right, right)))
    rng = np.random.default_rng(seed)
    null = np.asarray(
        [float(np.dot(left, rng.permutation(right))) / denominator for _ in range(permutations)]
    )
    return float(np.quantile(null, 0.025)), float(np.quantile(null, 0.975))


def permutation_auc(
    feature_rank: np.ndarray, labels: np.ndarray, permutations: int, seed: int
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    null = np.asarray(
        [auc_from_ranks(feature_rank, rng.permutation(labels)) for _ in range(permutations)]
    )
    return float(np.quantile(null, 0.025)), float(np.quantile(null, 0.975))


def load_pool(pool: str, config: dict[str, object]) -> list[dict[str, object]]:
    sample_dirs = [SOURCE / "data" / "Aging" / f"{sample}_{pool}" for sample in SAMPLES]
    design_paths = [sample_dir / "design_files.fasta" for sample_dir in sample_dirs]
    scaf_paths = [sample_dir / "scafstats.txt" for sample_dir in sample_dirs]
    if any(git_blob_sha1(path) != config["design_blob"] for path in design_paths):
        raise ValueError(f"{pool}: design blob mismatch")
    for path, expected in zip(scaf_paths, config["scaf_blobs"]):
        if git_blob_sha1(path) != expected:
            raise ValueError(f"{pool}: scafstats blob mismatch for {path}")
    references = parse_fasta(design_paths[0])
    if len(references) != config["reference_count"]:
        raise ValueError(f"{pool}: reference count mismatch")
    counts = [parse_scafstats(path) for path in scaf_paths]
    names = set(references)
    if any(set(sample_counts) - names for sample_counts in counts):
        raise ValueError(f"{pool}: unexpected scafstats names")
    missing = [len(names - set(sample_counts)) for sample_counts in counts]
    if missing != config["missing_counts"]:
        raise ValueError(f"{pool}: missing-count vector mismatch: {missing}")
    totals = [sum(sample_counts.values()) for sample_counts in counts]

    rows: list[dict[str, object]] = []
    for name, sequence in references.items():
        raw_counts = [sample_counts.get(name, 0) for sample_counts in counts]
        cpm = [count / total * 1_000_000.0 for count, total in zip(raw_counts, totals)]
        day0_cpm = statistics.fmean(cpm[:2])
        day7_cpm = cpm[-1]
        features = hrf.hairpin_features_simple(sequence, 4, 3)
        gc_fraction = hrf.gc_fraction(sequence)
        rows.append(
            {
                "pool": pool,
                "reference_id": name,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "length_nt": len(sequence),
                **{f"assigned_reads_{sample}": raw_counts[index] for index, sample in enumerate(SAMPLES)},
                **{f"cpm_{sample}": f"{cpm[index]:.9f}" for index, sample in enumerate(SAMPLES)},
                "mean_log2_cpm_plus1": f"{statistics.fmean(math.log2(value + 1.0) for value in cpm):.9f}",
                "day7_vs_day0_log2fc_plus1": f"{math.log2((day7_cpm + 1.0) / (day0_cpm + 1.0)):.9f}",
                "missing_any": any(count == 0 for count in raw_counts),
                "missing_day7": raw_counts[-1] == 0,
                "gc_fraction": f"{gc_fraction:.9f}",
                "gc_deviation_from_0p5": f"{abs(gc_fraction - 0.5):.9f}",
                "max_homopolymer": hrf.max_homopolymer(sequence),
                "weighted_pairs": features["weighted_pairs"],
                "candidate_pairs": features["candidate_pairs"],
                "longest_stem": features["longest_stem"],
                "claim_boundary": BOUNDARY,
            }
        )
    return rows


def analyze(pool_rows: list[dict[str, object]], permutations: int, seed_offset: int) -> list[dict[str, str]]:
    pool = str(pool_rows[0]["pool"])
    metrics = ["weighted_pairs", "candidate_pairs", "max_homopolymer", "gc_deviation_from_0p5"]
    outcomes = [
        ("mean_log2_cpm_plus1", "spearman"),
        ("day7_vs_day0_log2fc_plus1", "spearman"),
        ("missing_any", "auc"),
        ("missing_day7", "auc"),
    ]
    output: list[dict[str, str]] = []
    for metric_index, metric in enumerate(metrics):
        feature = np.asarray([float(row[metric]) for row in pool_rows])
        feature_rank = rank_average(feature)
        for outcome_index, (outcome, statistic_name) in enumerate(outcomes):
            seed = 20260713 + seed_offset + metric_index * 10 + outcome_index
            if statistic_name == "spearman":
                values = np.asarray([float(row[outcome]) for row in pool_rows])
                outcome_rank = rank_average(values)
                statistic = spearman_from_ranks(feature_rank, outcome_rank)
                null_low, null_high = permutation_continuous(
                    feature_rank, outcome_rank, permutations, seed
                )
                event_count = ""
            else:
                labels = np.asarray([row[outcome] is True for row in pool_rows], dtype=int)
                statistic = auc_from_ranks(feature_rank, labels)
                null_low, null_high = permutation_auc(feature_rank, labels, permutations, seed)
                event_count = str(int(labels.sum()))
            output.append(
                {
                    "status": "PASS_EXTERNAL_DETECTABILITY_CALIBRATION",
                    "pool": pool,
                    "metric": metric,
                    "outcome": outcome,
                    "statistic": statistic_name,
                    "n_sequences": str(len(pool_rows)),
                    "event_count": event_count,
                    "observed": f"{statistic:.9f}",
                    "permutation_null_2p5": f"{null_low:.9f}",
                    "permutation_null_97p5": f"{null_high:.9f}",
                    "outside_null_envelope": str(statistic < null_low or statistic > null_high),
                    "claim_boundary": BOUNDARY,
                }
            )
    return output


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--detail", type=Path, default=DETAIL)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    parser.add_argument("--note", type=Path, default=NOTE)
    args = parser.parse_args()

    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, str]] = []
    for pool_index, (pool, config) in enumerate(POOLS.items()):
        pool_rows = load_pool(pool, config)
        detail_rows.extend(pool_rows)
        summary_rows.extend(analyze(pool_rows, args.permutations, pool_index * 100))
    write_csv(args.detail, detail_rows)
    write_csv(args.summary, summary_rows)

    weighted_rows = [row for row in summary_rows if row["metric"] == "weighted_pairs"]
    args.note.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Paper 2 DT4DDS sequence-detectability calibration",
        "",
        "Overall status: `PASS_EXTERNAL_DETECTABILITY_CALIBRATION`",
        "",
        "## Frozen source",
        "",
        "- Upstream: `fml-ethz/dt4dds_notebooks`, tag `v1.0.0`, commit `126e6da5c41f4e5de072b7a1a0934068b743de6c`.",
        "- Pools: 12,472 Genscript GCall references of length 102 nt and 12,000 Twist GCall references of length 108 nt, each observed at `0a`, `0b`, `2d`, `4d`, and `7d`.",
        "- Every design FASTA and scafstats Git blob SHA-1 is checked before analysis.",
        "- The upstream repository does not expose a machine-readable license; redistribution terms must be confirmed before public bundling.",
        "",
        "## Analysis",
        "",
        f"- Sequence rows: {len(detail_rows):,}; summary rows: {len(summary_rows)}; deterministic permutations per row: {args.permutations}.",
        "- Coverage is normalized to counts per million within each sample.",
        "- Continuous outcomes are mean log2(CPM+1) and day-7 versus mean-day-0 log2 fold change.",
        "- Binary outcomes are any-timepoint missing detection and day-7 missing detection.",
        "- All associations are computed separately within each supplier pool.",
        "",
        "## Weighted-pair results",
        "",
    ]
    for row in weighted_rows:
        lines.append(
            f"- `{row['pool']}` / `{row['outcome']}`: {row['statistic']} `{row['observed']}`, "
            f"permutation 95% envelope `[{row['permutation_null_2p5']}, {row['permutation_null_97p5']}]`, "
            f"outside envelope `{row['outside_null_envelope']}`."
        )
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            BOUNDARY + ".",
            "The source experiment combines material aging with library preparation, PCR, alignment, and sequencing-depth effects. Missing detection and CPM therefore cannot be interpreted as isolated molecular dropout, empirical file recovery, or a causal benefit from imposing the Paper 2 proxy.",
            "",
        ]
    )
    args.note.write_text("\n".join(lines), encoding="utf-8")
    print("PASS_EXTERNAL_DETECTABILITY_CALIBRATION")
    print(f"detail_rows={len(detail_rows)} summary_rows={len(summary_rows)}")
    for row in weighted_rows:
        print(
            f"{row['pool']} {row['outcome']} {row['statistic']}={row['observed']} "
            f"outside={row['outside_null_envelope']}"
        )


if __name__ == "__main__":
    main()
