#!/usr/bin/env python3
"""Minimal deterministic encode/decode example for the Paper 2 shaping layer.

The toy score below is deliberately not the manuscript's fitted PCR selector.
It demonstrates the scorer interface and scorer-independent payload inversion
without requiring public assay inputs.
"""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analysis_tools"))

import paper2_deterministic_selection as deterministic  # noqa: E402
import validate_paper2_reversible_choice_codec as choice  # noqa: E402


def toy_score(sequence: str) -> float:
    """Prefer a GC count near 54; this is not an assay model."""

    gc = sequence.count("G") + sequence.count("C")
    return -abs(gc - 54) / 108.0


def main() -> int:
    selector_bits = 2
    payload = 42
    codec = choice.GCHomopolymerCodec(
        choice.GENERATED_LENGTH,
        choice.GENERATED_GC_MIN,
        choice.GENERATED_GC_MAX,
        choice.GENERATED_MAX_HOMOPOLYMER,
    )
    affine = choice.affine_record(codec, fixed_bits=213, namespace_index=0)
    modulus = int(affine["language_modulus"])

    candidates: list[str] = []
    for index in range(1 << selector_bits):
        logical_rank = (payload << selector_bits) | index
        physical_rank = (
            int(affine["affine_multiplier"]) * logical_rank
            + int(affine["affine_offset"])
        ) % modulus
        candidates.append(codec.unrank(physical_rank))

    keys = deterministic.quantized_score_keys([toy_score(x) for x in candidates])
    choice_index = int(keys.argmax())  # np.argmax returns the smallest index on a tie.
    emitted = candidates[choice_index]

    physical_rank = codec.rank(emitted)
    logical_rank = (
        int(affine["affine_inverse"])
        * (physical_rank - int(affine["affine_offset"]))
    ) % modulus
    if logical_rank >= 1 << 213:
        raise ValueError("received sequence lies outside the declared logical domain")
    decoded_payload = logical_rank >> selector_bits
    decoded_choice_index = logical_rank & ((1 << selector_bits) - 1)

    assert decoded_payload == payload
    assert decoded_choice_index == choice_index
    print(f"payload={payload}")
    print(f"choice_index={choice_index}")
    print(f"sequence={emitted}")
    print(f"decoded_payload={decoded_payload}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
