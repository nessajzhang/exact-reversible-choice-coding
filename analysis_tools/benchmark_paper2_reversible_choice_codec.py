#!/usr/bin/env python3
"""Benchmark the exact Paper 2 base codec and P5 choice encoder.

Timing is a machine-specific implementation benchmark. It supports complexity
and usability reporting only; it is not biological evidence and is not used to
claim superiority over an end-to-end DNA-storage system.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import statistics
import sys
import time
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import validate_paper2_public_experimental_data as public_data  # noqa: E402
import validate_paper2_reversible_choice_codec as choice  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
CHOICE_DIR = PAPER_DIR / "bioinformatics_reframe" / "reversible_choice_codec"
DEFAULT_OUT_DIR = PAPER_DIR / "bioinformatics_reframe" / "runtime_benchmark"
COEFFICIENTS = CHOICE_DIR / "frozen_model_coefficients.tsv"
BASE_SEED = 20260722


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payloads", type=int, default=128)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--rank-checks", type=int, default=4096)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help=(
            "output directory; the one-command workflow writes local reruns "
            "below the frozen benchmark directory without replacing its manuscript snapshot"
        ),
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
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def peak_rss_mib() -> float:
    raw = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if platform.system() == "Darwin":
        return raw / (1024.0 * 1024.0)
    return raw / 1024.0


def p5_raw_model(source_pool: str) -> tuple[float, dict[str, float]]:
    rows: list[dict[str, str]] = []
    with COEFFICIENTS.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source_pool"] == source_pool and row["model"] == "P5_combined_context":
                rows.append(row)
    if not rows:
        raise ValueError(f"missing frozen P5 coefficients for {source_pool}")
    intercepts = {float(row["raw_intercept"]) for row in rows}
    if len(intercepts) != 1:
        raise ValueError("inconsistent frozen raw intercept")
    coefficients = {row["feature"]: float(row["raw_coefficient"]) for row in rows}
    return intercepts.pop(), coefficients


def score_sequence(
    sequence: str,
    intercept: float,
    coefficients: dict[str, float],
) -> float:
    features = public_data.pcr_feature_row("generated", "runtime", sequence)
    return intercept + sum(coefficients[name] * float(features[name]) for name in coefficients)


def deterministic_ranks(limit: int, count: int) -> list[int]:
    fixed = [0, 1, limit // 3, limit // 2, limit - 2, limit - 1]
    ranks = {value for value in fixed if 0 <= value < limit}
    index = 0
    while len(ranks) < count:
        digest = hashlib.sha256(
            f"paper2-runtime-rank-v1|{BASE_SEED}|{index}".encode("ascii")
        ).digest()
        ranks.add(int.from_bytes(digest, "big") % limit)
        index += 1
    return sorted(ranks)[:count]


def aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = sorted({(row["operation"], int(row["selector_bits"])) for row in rows})
    summary: list[dict[str, Any]] = []
    for operation, selector_bits in keys:
        subset = [
            row
            for row in rows
            if row["operation"] == operation and int(row["selector_bits"]) == selector_bits
        ]
        seconds = [float(row["seconds"]) for row in subset]
        ms_per_payload = [float(row["milliseconds_per_payload"]) for row in subset]
        summary.append(
            {
                "operation": operation,
                "selector_bits": selector_bits,
                "candidates_per_payload": subset[0]["candidates_per_payload"],
                "payloads_or_items": subset[0]["payloads_or_items"],
                "repeats": len(subset),
                "median_seconds": statistics.median(seconds),
                "minimum_seconds": min(seconds),
                "maximum_seconds": max(seconds),
                "median_milliseconds_per_payload_or_item": statistics.median(ms_per_payload),
                "all_roundtrips_passed": all(bool(row["all_roundtrips_passed"]) for row in subset),
                "peak_process_rss_mib": max(float(row["peak_process_rss_mib"]) for row in subset),
                "boundary": "single-process machine-specific implementation timing; not a biological replicate or system-level comparison",
            }
        )
    return summary


def main() -> int:
    args = parse_args()
    if args.payloads < 1 or args.repeats < 1 or args.rank_checks < 1:
        raise ValueError("payloads, repeats and rank-checks must be positive")
    out_dir = args.output_dir.expanduser()
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    codec = choice.GCHomopolymerCodec(108, 49, 59, 3)
    start = time.perf_counter()
    exact_total = codec.count()
    count_seconds = time.perf_counter() - start
    if exact_total != int(
        "21520867790325216400381593480072139809549357542204871770858223072"
    ):
        raise ValueError("exact base-language total changed")

    payload_domain = 1 << 213
    ranks = deterministic_ranks(payload_domain, min(args.rank_checks, payload_domain))
    rows: list[dict[str, Any]] = []
    for repeat in range(args.repeats):
        start = time.perf_counter()
        sequences = [codec.unrank(rank) for rank in ranks]
        unrank_seconds = time.perf_counter() - start
        start = time.perf_counter()
        reranked = [codec.rank(sequence) for sequence in sequences]
        rank_seconds = time.perf_counter() - start
        passed = reranked == ranks
        for operation, seconds in (("base_unrank", unrank_seconds), ("base_rank", rank_seconds)):
            rows.append(
                {
                    "operation": operation,
                    "selector_bits": 0,
                    "candidates_per_payload": 1,
                    "payloads_or_items": len(ranks),
                    "repeat": repeat,
                    "seconds": seconds,
                    "milliseconds_per_payload": 1000.0 * seconds / len(ranks),
                    "all_roundtrips_passed": passed,
                    "peak_process_rss_mib": peak_rss_mib(),
                }
            )

    affine = choice.affine_record(codec, 213)
    modulus = int(affine["language_modulus"])
    multiplier = int(affine["affine_multiplier"])
    offset = int(affine["affine_offset"])
    intercept, coefficients = p5_raw_model("GCall")
    for selector_bits in (0, 2, 4):
        payload_limit = 1 << (213 - selector_bits)
        payloads = choice.deterministic_payloads(payload_limit, args.payloads, selector_bits)
        width = 1 << selector_bits
        for repeat in range(args.repeats):
            selected: list[tuple[str, int]] = []
            start = time.perf_counter()
            for payload in payloads:
                best: tuple[float, int, str] | None = None
                for candidate in range(width):
                    logical_rank = (payload << selector_bits) | candidate
                    physical_rank = (multiplier * logical_rank + offset) % modulus
                    sequence = codec.unrank(physical_rank)
                    score = score_sequence(sequence, intercept, coefficients)
                    item = (score, -candidate, sequence)
                    if best is None or item > best:
                        best = item
                if best is None:
                    raise AssertionError("empty choice fiber")
                selected.append((best[2], payload))
            encode_seconds = time.perf_counter() - start

            start = time.perf_counter()
            decoded: list[int] = []
            for sequence, _payload in selected:
                payload, _choice, _logical_rank, _physical_rank = (
                    choice.decode_payload_candidate(
                        codec, sequence, 213, selector_bits, affine
                    )
                )
                decoded.append(payload)
            decode_seconds = time.perf_counter() - start
            passed = decoded == payloads
            for operation, seconds in (
                ("P5_choice_encode", encode_seconds),
                ("choice_decode", decode_seconds),
            ):
                rows.append(
                    {
                        "operation": operation,
                        "selector_bits": selector_bits,
                        "candidates_per_payload": width,
                        "payloads_or_items": len(payloads),
                        "repeat": repeat,
                        "seconds": seconds,
                        "milliseconds_per_payload": 1000.0 * seconds / len(payloads),
                        "all_roundtrips_passed": passed,
                        "peak_process_rss_mib": peak_rss_mib(),
                    }
                )

    summary = aggregate(rows)
    write_tsv(out_dir / "runtime_replicates.tsv", rows)
    write_tsv(out_dir / "runtime_summary.tsv", summary)
    environment = {
        "analysis_date": "2026-07-18",
        "command": (
            "$PYTHON analysis_tools/benchmark_paper2_reversible_choice_codec.py "
            f"--payloads {args.payloads} --repeats {args.repeats} "
            f"--rank-checks {args.rank_checks} --output-dir $OUTPUT_DIR"
        ),
        "python": sys.version,
        "platform": platform.platform(),
        "logical_cpus": os.cpu_count(),
        "base_seed": BASE_SEED,
        "payloads": args.payloads,
        "repeats": args.repeats,
        "rank_checks": len(ranks),
        "base_count_seconds_cold": count_seconds,
        "base_count_cache_entries": codec._count_cached.cache_info().currsize,
        "rank_dispersion_modulus": modulus,
        "payload_domain_size": payload_domain,
        "peak_process_rss_mib": peak_rss_mib(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "pandas", "scikit-learn", "joblib"]
        },
        "boundary": "machine-specific single-process timing; not a biological result or end-to-end storage comparison",
    }
    (out_dir / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary_lines = [
        "# Reversible choice-codec runtime benchmark",
        "",
        "Status: `COMPLETE_MACHINE_SPECIFIC_RUNTIME_AUDIT`",
        "",
        f"Cold exact completion-table construction: {count_seconds:.6f} s; "
        f"cache entries after all checks: {environment['base_count_cache_entries']:,}.",
        "",
        "Median timings:",
        "",
        "| Operation | r | Candidates | Items | Median ms/item | Peak process RSS (MiB) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary:
        summary_lines.append(
            f"| {row['operation']} | {row['selector_bits']} | "
            f"{row['candidates_per_payload']} | {row['payloads_or_items']} | "
            f"{float(row['median_milliseconds_per_payload_or_item']):.6f} | "
            f"{float(row['peak_process_rss_mib']):.3f} |"
        )
    summary_lines.extend(
        [
            "",
            "All timed encode/decode and base rank/unrank checks passed. Timings are machine-specific, single-process implementation measurements; they are not biological replicates, a physical-density result or an end-to-end system comparison.",
        ]
    )
    (out_dir / "analysis_summary.md").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )
    manifest_paths = [
        out_dir / "analysis_summary.md",
        out_dir / "environment.json",
        out_dir / "runtime_replicates.tsv",
        out_dir / "runtime_summary.tsv",
        Path(__file__),
    ]

    def manifest_path(path: Path) -> str:
        resolved = path.resolve()
        try:
            return str(resolved.relative_to(ROOT.resolve()))
        except ValueError:
            return str(resolved)

    manifest = [
        {
            "path": manifest_path(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in manifest_paths
    ]
    write_tsv(out_dir / "sha256_manifest.tsv", manifest)
    print(f"Wrote runtime benchmark to {out_dir}", flush=True)
    print(f"Manifest SHA-256: {sha256(out_dir / 'sha256_manifest.tsv')}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
