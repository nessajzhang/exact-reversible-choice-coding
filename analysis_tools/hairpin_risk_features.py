#!/usr/bin/env python3
"""Legacy implementation of deterministic reverse-complement burden features.

The filename and ``hairpin_features*`` symbols are retained for replay of
historical analyses.  New public code should import
``self_complementarity_features.weighted_pair_features``.  These sequence-only
features are not thermodynamic folding predictions or observed hairpins.
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path


BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(1 for x in seq if x in "GC") / len(seq)


def max_homopolymer(seq: str) -> int:
    best = 0
    cur = 0
    prev = None
    for ch in seq:
        if ch == prev:
            cur += 1
        else:
            prev = ch
            cur = 1
        best = max(best, cur)
    return best


def hairpin_features(seq: str, min_stem: int, min_loop: int) -> dict[str, int | float]:
    n = len(seq)
    candidate_pairs = 0
    weighted_pairs = 0
    longest_stem = 0
    min_loop_seen = None

    for i in range(n):
        for j in range(i + min_stem + min_loop, n - min_stem + 1):
            stem = 0
            while i + stem < j and j + min_stem - 1 - stem < n:
                if seq[i + stem] != rc(seq[j:j + min_stem + stem])[-1 - stem]:
                    break
                stem += 1
                if stem >= min_stem:
                    candidate_pairs += 1
                    weighted_pairs += stem
                    longest_stem = max(longest_stem, stem)
                    loop = j - i - stem
                    if loop >= min_loop:
                        min_loop_seen = loop if min_loop_seen is None else min(min_loop_seen, loop)
            # The while condition above extends from the right window; stop if impossible.

    return {
        "candidate_pairs": candidate_pairs,
        "weighted_pairs": weighted_pairs,
        "longest_stem": longest_stem,
        "min_loop_seen": -1 if min_loop_seen is None else min_loop_seen,
    }


def hairpin_features_simple(seq: str, min_stem: int, min_loop: int) -> dict[str, int | float]:
    """Exact reverse-complement stem count using anti-diagonal DP.

    For a left stem starting at i and a right stem ending at q, the paired bases
    follow the anti-diagonal (i, q), (i+1, q-1), ... .  The DP entry
    run[i][q] stores how many consecutive Watson-Crick matches continue along
    that anti-diagonal.  This avoids repeated substring reverse-complement
    comparisons in batch experiments.
    """
    n = len(seq)
    run = [[0] * n for _ in range(n + 1)]
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    for i in range(n - 1, -1, -1):
        for q in range(n):
            if seq[i] == comp.get(seq[q], ""):
                run[i][q] = 1 + (run[i + 1][q - 1] if q > 0 else 0)

    candidate_pairs = 0
    weighted_pairs = 0
    longest_stem = 0
    min_loop_seen = None
    for i in range(n):
        for q in range(i + 2 * min_stem + min_loop - 1, n):
            max_by_loop = (q - i - min_loop + 1) // 2
            max_stem = min(run[i][q], max_by_loop)
            if max_stem < min_stem:
                continue
            count = max_stem - min_stem + 1
            candidate_pairs += count
            weighted_pairs += (min_stem + max_stem) * count // 2
            longest_stem = max(longest_stem, max_stem)
            loop = q - i - 2 * max_stem + 1
            min_loop_seen = loop if min_loop_seen is None else min(min_loop_seen, loop)
    return {
        "candidate_pairs": candidate_pairs,
        "weighted_pairs": weighted_pairs,
        "longest_stem": longest_stem,
        "min_loop_seen": -1 if min_loop_seen is None else min_loop_seen,
    }


def read_sequences(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    if path.suffix.lower() in {".fa", ".fasta"}:
        name = None
        chunks: list[str] = []
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name = line[1:].strip() or f"seq{len(records)}"
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            records.append((name, "".join(chunks).upper()))
    else:
        for idx, line in enumerate(path.read_text(encoding="utf-8").splitlines()):
            line = line.strip()
            if not line:
                continue
            parts = line.replace(",", "\t").split()
            if len(parts) == 1:
                records.append((f"seq{idx}", parts[0].upper()))
            else:
                records.append((parts[0], parts[1].upper()))
    return records


def random_sequences(count: int, length: int, seed: int) -> list[tuple[str, str]]:
    rng = random.Random(seed)
    return [(f"random_{i}", "".join(rng.choice(BASES) for _ in range(length))) for i in range(count)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path)
    parser.add_argument("--random-count", type=int)
    parser.add_argument("--length", type=int, default=110)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--min-stem", type=int, default=4)
    parser.add_argument("--min-loop", type=int, default=3)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    if args.input:
        records = read_sequences(args.input)
    elif args.random_count:
        records = random_sequences(args.random_count, args.length, args.seed)
    else:
        raise SystemExit("provide --input or --random-count")

    rows = []
    for name, seq in records:
        features = hairpin_features_simple(seq, args.min_stem, args.min_loop)
        rows.append(
            {
                "id": name,
                "sequence": seq,
                "length": len(seq),
                "gc_fraction": f"{gc_fraction(seq):.6f}",
                "max_homopolymer": max_homopolymer(seq),
                "min_stem": args.min_stem,
                "min_loop": args.min_loop,
                **features,
            }
        )

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
