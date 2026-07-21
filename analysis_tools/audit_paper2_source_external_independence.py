#!/usr/bin/env python3
"""Exhaustive GCall/GCfix versus locked-external sequence audit for Paper 2.

The 108-nt variable regions are compared exactly.  Shared fixed adapters are
reported separately as assay context and are not counted as sequence-set
independence.  The audit covers both all 2,053 exact-codec-eligible external
sequences and the outcome-blind 2,048-sequence locked codebook.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
from joblib import Parallel, delayed


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import audit_paper2_pcr_sequence_independence as distance  # noqa: E402
import validate_paper2_public_experimental_data as public_data  # noqa: E402
import validate_paper2_reversible_choice_codec as codec  # noqa: E402
import validate_paper2_sota_external_uncertainty as external  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
OUT_DIR = PAPER_DIR / "bioinformatics_reframe" / "source_external_independence_audit"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(10, (os.cpu_count() or 2) - 2)),
    )
    parser.add_argument("--block-size", type=int, default=64)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def audit_pair(
    source_pool: str,
    source,
    target_name: str,
    target,
    jobs: int,
    block_size: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    source_sequences = source["sequence"].astype(str).tolist()
    target_sequences = target["sequence"].astype(str).tolist()
    source_ids = source["sample_id"].astype(str).tolist()
    target_ids = target["seq_id_anonymized"].astype(str).tolist()

    exact = set(source_sequences).intersection(target_sequences)
    reverse_complements = set(map(distance.reverse_complement, source_sequences)).intersection(
        target_sequences
    )
    source_packed = distance.pack_sequences(source_sequences)
    target_packed = distance.pack_sequences(target_sequences)
    tasks = [
        (
            start,
            source_packed[start : min(start + block_size, len(source_packed))],
            target_packed,
        )
        for start in range(0, len(source_packed), block_size)
    ]
    results = Parallel(n_jobs=jobs, prefer="threads", batch_size=1)(
        delayed(distance.compare_block)(*task) for task in tasks
    )

    source_nearest = np.empty(len(source_sequences), dtype=np.uint8)
    target_nearest = np.full(len(target_sequences), distance.LENGTH + 1, dtype=np.uint8)
    counts = {threshold: 0 for threshold in distance.THRESHOLDS}
    closest = (distance.LENGTH + 1, -1, -1)
    for start, source_minimum, target_minimum, block_counts, block_closest in results:
        source_nearest[start : start + len(source_minimum)] = source_minimum
        target_nearest = np.minimum(target_nearest, target_minimum)
        for threshold, count in block_counts.items():
            counts[threshold] += count
        closest = min(closest, block_closest)

    if counts[0] != len(exact):
        raise ValueError("packed exact-match count disagrees with set intersection")
    minimum, source_index, target_index = closest
    summary = {
        "audit_scope": f"{source_pool}_variable_regions_vs_{target_name}",
        "sequence_length_nt": distance.LENGTH,
        "source_pool": source_pool,
        "source_sequences": len(source_sequences),
        "external_cohort": target_name,
        "external_sequences": len(target_sequences),
        "all_oriented_pairs": len(source_sequences) * len(target_sequences),
        "source_within_duplicates": len(source_sequences) - len(set(source_sequences)),
        "external_within_duplicates": len(target_sequences) - len(set(target_sequences)),
        "exact_duplicates": len(exact),
        "exact_reverse_complements": len(reverse_complements),
        "minimum_oriented_hamming_distance_nt": int(minimum),
        "maximum_oriented_identity_fraction": float(
            (distance.LENGTH - minimum) / distance.LENGTH
        ),
        "closest_source_sample_id": source_ids[source_index],
        "closest_external_sample_id": target_ids[target_index],
        "pairs_identity_at_least_90_percent": counts[10],
        "pairs_identity_at_least_80_percent": counts[21],
        "pairs_identity_at_least_70_percent": counts[32],
        "fixed_adapter_handling": "excluded_shared_declared_assay_context",
        "interpretation": "exhaustive identity audit over oriented 108-nt variable regions",
    }
    distributions: list[dict[str, Any]] = []
    for side, values in ((source_pool, source_nearest), (target_name, target_nearest)):
        for quantile in (0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0):
            hamming = float(np.quantile(values, quantile, method="nearest"))
            distributions.append(
                {
                    "audit_scope": summary["audit_scope"],
                    "population": side,
                    "quantile": quantile,
                    "nearest_hamming_distance_nt": hamming,
                    "nearest_identity_fraction": 1 - hamming / distance.LENGTH,
                }
            )
    return summary, distributions


def main() -> None:
    args = parse_args()
    if args.jobs < 1 or args.block_size < 1:
        raise ValueError("jobs and block size must be positive")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pcr, _audits = public_data.load_pcr_data()
    external_random, _external_audits = external.load_external_pool()
    eligible = external_random.loc[
        external_random["sequence"].map(codec.exact_gc_hp_eligible)
    ].copy()
    locked = external.locked_external_library(external_random)
    if len(eligible) != 2053 or len(locked) != 2048:
        raise ValueError("external eligible/locked cohort cardinality changed")

    summaries: list[dict[str, Any]] = []
    distributions: list[dict[str, Any]] = []
    for source_pool in ("GCall", "GCfix"):
        source = pcr.loc[pcr["pool"].eq(source_pool)].reset_index(drop=True)
        for target_name, target in (
            ("external_codec_eligible_2053", eligible),
            ("external_locked_codebook_2048", locked),
        ):
            summary, rows = audit_pair(
                source_pool, source, target_name, target, args.jobs, args.block_size
            )
            summaries.append(summary)
            distributions.extend(rows)

    adapter_rows = [
        {
            "construct_scope": "GCall_GCfix_and_external_feature_reconstruction",
            "left_adapter_name": "0F",
            "left_adapter_5prime_to_3prime": public_data.LEFT_ADAPTER_0F,
            "left_adapter_length_nt": len(public_data.LEFT_ADAPTER_0F),
            "right_adapter_name": "0R-prime (synthesized reverse-complemented strand)",
            "right_adapter_5prime_to_3prime": public_data.RIGHT_ADAPTER_0R_PRIME,
            "right_adapter_length_nt": len(public_data.RIGHT_ADAPTER_0R_PRIME),
            "full_construct_length_nt": (
                len(public_data.LEFT_ADAPTER_0F)
                + codec.GENERATED_LENGTH
                + len(public_data.RIGHT_ADAPTER_0R_PRIME)
            ),
            "feature_reconstruction_contract": "the same publication-defined terminal adapters are prepended/appended to every 108-nt variable region before P2/FullContext feature calculation",
            "sequence_independence_contract": "the shared 41 adapter nucleotides are declared assay context and excluded from identity thresholds",
        }
    ]

    write_tsv(OUT_DIR / "source_external_independence_summary.tsv", summaries)
    write_tsv(OUT_DIR / "nearest_source_external_distance_distribution.tsv", distributions)
    write_tsv(OUT_DIR / "adapter_reconstruction_contract.tsv", adapter_rows)

    environment = {
        "command": (
            f"$PYTHON analysis_tools/{Path(__file__).name} --jobs {args.jobs} "
            f"--block-size {args.block_size}"
        ),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "jobs": args.jobs,
        "block_size": args.block_size,
        "distance_definition": "exact oriented Hamming distance over 108-nt variable regions",
        "threshold_distance_cutoffs": {"90pct": 10, "80pct": 21, "70pct": 32},
    }
    (OUT_DIR / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    lines = [
        "# Source--external fixed-length Hamming-separation audit",
        "",
        "Status: `COMPLETE_EXHAUSTIVE_SOURCE_EXTERNAL_AUDIT`",
        "",
    ]
    for row in summaries:
        lines.append(
            f"- {row['source_pool']} vs {row['external_cohort']}: "
            f"{row['all_oriented_pairs']:,} exact pair comparisons; minimum Hamming "
            f"distance {row['minimum_oriented_hamming_distance_nt']}/108 "
            f"(maximum identity {row['maximum_oriented_identity_fraction']:.3%}); "
            f"duplicates={row['exact_duplicates']}, reverse complements="
            f"{row['exact_reverse_complements']}, pairs at >=70% identity="
            f"{row['pairs_identity_at_least_70_percent']}."
        )
    lines.extend(
        [
            "",
            "The fixed 0F and 0R-prime sequences are assay context, not variable content. This audit establishes only the reported fixed-length Hamming separation; it does not test edit-distance neighbourhoods, design genealogy or cluster-level independence and does not create biological replication or prospective codec-output evidence.",
            "",
        ]
    )
    (OUT_DIR / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")

    manifest_rows: list[dict[str, Any]] = []
    for path in sorted(OUT_DIR.iterdir()):
        if path.is_file() and path.name != "sha256_manifest.tsv":
            manifest_rows.append(
                {
                    "path": str(path.relative_to(ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    manifest_rows.append(
        {
            "path": str(Path(__file__).relative_to(ROOT)),
            "bytes": Path(__file__).stat().st_size,
            "sha256": sha256(Path(__file__)),
        }
    )
    write_tsv(OUT_DIR / "sha256_manifest.tsv", manifest_rows)
    print("\n".join(lines).strip())
    print(f"manifest_sha256={sha256(OUT_DIR / 'sha256_manifest.tsv')}")


if __name__ == "__main__":
    main()
