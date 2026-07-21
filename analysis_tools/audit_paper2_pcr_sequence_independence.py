#!/usr/bin/env python3
"""Audit fixed-length Hamming separation for the Paper 2 PCR analysis.

Bioinformatics requires machine-learning studies on biological sequences to
describe identical and near-identical examples across data roles. The GCall and GCfix libraries are synthetic random
108-nt variable regions rather than homologous biological sequences, but that
fact does not replace a direct sequence audit.  This program therefore checks
all cross-pool pairs exactly using two-bit packed Hamming distance.

The fixed PCR adapters are shared assay context and are intentionally excluded
from the distance calculation.  The audit concerns only the variable sequence
used as the model's sequence-bearing analysis unit. Edit distance, design
genealogy and cluster-level independence are outside this audit.
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

import validate_paper2_public_experimental_data as public_data  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
OUT_DIR = PAPER_DIR / "bioinformatics_reframe" / "sequence_independence_audit"
LENGTH = 108
THRESHOLDS = (0, 10, 21, 32)
BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(10, (os.cpu_count() or 2) - 2)),
        help="Thread workers for exact cross-pool distance blocks.",
    )
    parser.add_argument(
        "--block-size",
        type=int,
        default=64,
        help="Number of GCall sequences per all-GCfix comparison block.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty output: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def pack_sequences(sequences: list[str]) -> np.ndarray:
    words = (LENGTH + 63) // 64
    packed = np.zeros((len(sequences), 2, words), dtype=np.uint64)
    for row, sequence in enumerate(sequences):
        if len(sequence) != LENGTH or set(sequence) - set(BASE_CODE):
            raise ValueError(f"invalid variable sequence at row {row}")
        for position, nucleotide in enumerate(sequence):
            code = BASE_CODE[nucleotide]
            word = position // 64
            bit = np.uint64(position % 64)
            mask = np.uint64(1) << bit
            if code & 1:
                packed[row, 0, word] |= mask
            if code & 2:
                packed[row, 1, word] |= mask
    return packed


def compare_block(
    start: int,
    left: np.ndarray,
    right: np.ndarray,
) -> tuple[int, np.ndarray, np.ndarray, dict[int, int], tuple[int, int, int]]:
    xor_low = np.bitwise_xor(left[:, None, 0, :], right[None, :, 0, :])
    xor_high = np.bitwise_xor(left[:, None, 1, :], right[None, :, 1, :])
    mismatch = np.bitwise_or(xor_low, xor_high)
    distances = np.bitwise_count(mismatch).sum(axis=2, dtype=np.uint16)

    left_minimum = distances.min(axis=1).astype(np.uint8)
    right_minimum = distances.min(axis=0).astype(np.uint8)
    threshold_counts = {
        threshold: int(np.count_nonzero(distances <= threshold))
        for threshold in THRESHOLDS
    }
    flat_index = int(np.argmin(distances))
    local_left, right_index = np.unravel_index(flat_index, distances.shape)
    closest = (
        int(distances[local_left, right_index]),
        start + int(local_left),
        int(right_index),
    )
    return start, left_minimum, right_minimum, threshold_counts, closest


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main() -> None:
    args = parse_args()
    if args.jobs < 1 or args.block_size < 1:
        raise ValueError("jobs and block size must be positive")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frame, _audits = public_data.load_pcr_data()
    gcall = frame.loc[frame["pool"].eq("GCall")].reset_index(drop=True)
    gcfix = frame.loc[frame["pool"].eq("GCfix")].reset_index(drop=True)
    left_sequences = gcall["sequence"].astype(str).tolist()
    right_sequences = gcfix["sequence"].astype(str).tolist()

    if len(left_sequences) != 11_998 or len(right_sequences) != 11_994:
        raise ValueError("unexpected public PCR pool sizes")
    if len(set(left_sequences)) != len(left_sequences):
        raise ValueError("within-GCall duplicate sequence")
    if len(set(right_sequences)) != len(right_sequences):
        raise ValueError("within-GCfix duplicate sequence")

    exact_intersection = set(left_sequences).intersection(right_sequences)
    reverse_complement_intersection = set(map(reverse_complement, left_sequences)).intersection(
        right_sequences
    )

    left_packed = pack_sequences(left_sequences)
    right_packed = pack_sequences(right_sequences)
    tasks = [
        (
            start,
            left_packed[start : min(start + args.block_size, len(left_packed))],
            right_packed,
        )
        for start in range(0, len(left_packed), args.block_size)
    ]
    results = Parallel(n_jobs=args.jobs, prefer="threads", batch_size=1)(
        delayed(compare_block)(*task) for task in tasks
    )

    left_nearest = np.empty(len(left_sequences), dtype=np.uint8)
    right_nearest = np.full(len(right_sequences), LENGTH + 1, dtype=np.uint8)
    threshold_counts = {threshold: 0 for threshold in THRESHOLDS}
    closest = (LENGTH + 1, -1, -1)
    for start, block_minimum, target_minimum, counts, block_closest in results:
        left_nearest[start : start + len(block_minimum)] = block_minimum
        right_nearest = np.minimum(right_nearest, target_minimum)
        for threshold, count in counts.items():
            threshold_counts[threshold] += count
        closest = min(closest, block_closest)

    if threshold_counts[0] != len(exact_intersection):
        raise ValueError("packed exact-match count disagrees with set intersection")
    minimum_distance, left_index, right_index = closest
    maximum_identity = (LENGTH - minimum_distance) / LENGTH

    summary_rows: list[dict[str, Any]] = [
        {
            "audit_scope": "GCall_variable_regions_vs_GCfix_variable_regions",
            "sequence_length_nt": LENGTH,
            "gcall_sequences": len(left_sequences),
            "gcfix_sequences": len(right_sequences),
            "all_cross_pool_pairs": len(left_sequences) * len(right_sequences),
            "within_gcall_duplicates": len(left_sequences) - len(set(left_sequences)),
            "within_gcfix_duplicates": len(right_sequences) - len(set(right_sequences)),
            "cross_pool_exact_duplicates": len(exact_intersection),
            "cross_pool_exact_reverse_complements": len(reverse_complement_intersection),
            "minimum_cross_pool_hamming_distance": minimum_distance,
            "maximum_cross_pool_identity_fraction": maximum_identity,
            "closest_gcall_sample_id": gcall.iloc[left_index]["sample_id"],
            "closest_gcfix_sample_id": gcfix.iloc[right_index]["sample_id"],
            "pairs_identity_at_least_90_percent": threshold_counts[10],
            "pairs_identity_at_least_80_percent": threshold_counts[21],
            "pairs_identity_at_least_70_percent": threshold_counts[32],
            "fixed_adapter_handling": "excluded_shared_assay_context",
            "interpretation": (
                "exhaustive sequence-independence audit; synthetic variable regions have no "
                "cross-pool identical or near-identical examples at the reported thresholds"
            ),
        }
    ]
    write_tsv(OUT_DIR / "cross_pool_independence_summary.tsv", summary_rows)

    distribution_rows: list[dict[str, Any]] = []
    for pool, distances in (("GCall", left_nearest), ("GCfix", right_nearest)):
        for quantile in (0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0):
            distribution_rows.append(
                {
                    "pool": pool,
                    "nearest_sequence_in_other_pool": "GCfix" if pool == "GCall" else "GCall",
                    "quantile": quantile,
                    "nearest_hamming_distance_nt": float(
                        np.quantile(distances, quantile, method="nearest")
                    ),
                    "nearest_identity_fraction": float(
                        1 - np.quantile(distances, quantile, method="nearest") / LENGTH
                    ),
                }
            )
    write_tsv(OUT_DIR / "nearest_cross_pool_distance_distribution.tsv", distribution_rows)

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
        "input_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in public_data.PCR_FILES.values()
        },
        "distance_definition": (
            "exact Hamming distance over oriented 108-nt variable regions; common fixed adapters excluded"
        ),
    }
    (OUT_DIR / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    summary = f"""# PCR fixed-length Hamming-separation audit

Status: `COMPLETE_EXHAUSTIVE_CROSS_POOL_AUDIT`

All `{len(left_sequences) * len(right_sequences):,}` oriented GCall--GCfix variable-region pairs were compared exactly. The closest pair differed at `{minimum_distance}` of `{LENGTH}` positions (maximum identity `{maximum_identity:.3%}`). There were `{len(exact_intersection)}` exact duplicates, `{len(reverse_complement_intersection)}` exact cross-pool reverse complements, and `{threshold_counts[10]}` pairs with at least 90% identity.

The sequence-level analysis unit is the 108-nt synthetic variable region. The 41 shared adapter nucleotides are declared assay context and were excluded. This audit establishes the reported fixed-length Hamming separation only; it does not test edit-distance neighbourhoods, design genealogy or cluster-level independence and does not create biological replication or prospective validation.
"""
    (OUT_DIR / "analysis_summary.md").write_text(summary, encoding="utf-8")

    manifest_rows = []
    for path in sorted(OUT_DIR.iterdir()):
        if path.is_file() and path.name != "sha256_manifest.tsv":
            manifest_rows.append(
                {
                    "path": str(path.relative_to(ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    write_tsv(OUT_DIR / "sha256_manifest.tsv", manifest_rows)

    print(summary.strip())
    print(f"manifest_sha256={sha256(OUT_DIR / 'sha256_manifest.tsv')}")


if __name__ == "__main__":
    main()
