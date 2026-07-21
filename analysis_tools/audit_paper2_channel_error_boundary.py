#!/usr/bin/env python3
"""Audit the single-edit and checksum boundary of the Paper 2 shaping layer.

This analysis does not turn the reversible choice construction into a channel
code.  It quantifies how deterministic, FullContext-selected generated
codewords behave under every single substitution, insertion and deletion in a
fixed sample, and illustrates composition with CRC-16/CCITT-FALSE.  Insertions
and deletions are rejected by the declared fixed-length interface; substitutions
may remain inside the constrained language and silently decode to another
payload unless an external integrity layer detects them.
"""

from __future__ import annotations

import argparse
import binascii
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


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import paper2_deterministic_selection as deterministic  # noqa: E402
import validate_paper2_assay_calibrated_selection as selection  # noqa: E402
import validate_paper2_reversible_choice_codec as choice  # noqa: E402


PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
OUT_DIR = PAPER_DIR / "bioinformatics_reframe" / "channel_error_boundary"
BASE_SEED = 20260722
FIXED_BITS = 213
SELECTOR_BITS = 2
CRC_BITS = 16
PAYLOAD_BITS = FIXED_BITS - SELECTOR_BITS
MESSAGE_BITS = PAYLOAD_BITS - CRC_BITS
FROZEN_MESSAGE_SAMPLES = 32
SOURCES = ("GCall", "GCfix")
MODEL_NAME = "P5_combined_context"
DNA = "ACGT"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(12, (os.cpu_count() or 2) - 2)),
        help="Parallel feature-extraction jobs.",
    )
    parser.add_argument(
        "--message-samples",
        type=int,
        default=FROZEN_MESSAGE_SAMPLES,
        help=f"Deterministic CRC-valid messages per source (frozen: {FROZEN_MESSAGE_SAMPLES}).",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


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


def crc16_ccitt_false(data: bytes) -> int:
    """CRC-16/CCITT-FALSE: poly 0x1021, init 0xffff, no reflection/xorout."""

    return int(binascii.crc_hqx(data, 0xFFFF))


def message_bytes(message: int) -> bytes:
    if message < 0 or message >= 1 << MESSAGE_BITS:
        raise ValueError("message outside fixed CRC payload domain")
    return message.to_bytes((MESSAGE_BITS + 7) // 8, "big")


def payload_with_crc16(message: int) -> int:
    return (message << CRC_BITS) | crc16_ccitt_false(message_bytes(message))


def crc16_payload_valid(payload: int) -> bool:
    if payload < 0 or payload >= 1 << PAYLOAD_BITS:
        return False
    message = payload >> CRC_BITS
    observed = payload & ((1 << CRC_BITS) - 1)
    return observed == crc16_ccitt_false(message_bytes(message))


def deterministic_messages(count: int) -> list[int]:
    limit = 1 << MESSAGE_BITS
    fixed = {0, 1, limit // 3, limit // 2, limit - 2, limit - 1}
    messages = {value for value in fixed if 0 <= value < limit}
    index = 0
    while len(messages) < count:
        digest = hashlib.sha256(
            f"paper2-channel-error-message-v1|i={index}".encode("ascii")
        ).digest()
        messages.add(int.from_bytes(digest, "big") % limit)
        index += 1
    return sorted(messages)[:count]


def sequence_for_logical_rank(
    codec: choice.GCHomopolymerCodec,
    affine: dict[str, int],
    logical_rank: int,
) -> str:
    physical_rank = (
        int(affine["affine_multiplier"]) * logical_rank
        + int(affine["affine_offset"])
    ) % int(affine["language_modulus"])
    return codec.unrank(physical_rank)


def feature_scores(
    sequences: list[str],
    models: dict[tuple[str, str], choice.FrozenRidge],
    jobs: int,
) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    unique = list(dict.fromkeys(sequences))
    candidates = pd.DataFrame(
        {
            "sequence": unique,
            "sequence_sha256": [sha256_text(sequence) for sequence in unique],
        }
    )
    features = choice.generated_feature_frame(candidates, jobs)
    score_maps: dict[str, dict[str, float]] = {}
    for source in SOURCES:
        scores = models[(source, MODEL_NAME)].predict(features)
        score_maps[source] = dict(zip(features["sequence"], scores, strict=True))
    return features, score_maps


def decode_status(
    codec: choice.GCHomopolymerCodec,
    affine: dict[str, int],
    sequence: str,
    original_payload: int,
) -> tuple[str, int | None, int | None]:
    if len(sequence) != codec.length:
        return "rejected_length", None, None
    try:
        rank = codec.rank(sequence)
    except ValueError:
        return "rejected_constraint_or_alphabet", None, None
    logical_rank = (
        int(affine["affine_inverse"])
        * (rank - int(affine["affine_offset"]))
    ) % int(affine["language_modulus"])
    if logical_rank >= 1 << FIXED_BITS:
        return "rejected_outside_power_of_two_domain", None, None
    payload = logical_rank >> SELECTOR_BITS
    index = logical_rank & ((1 << SELECTOR_BITS) - 1)
    status = "accepted_same_payload" if payload == original_payload else "accepted_wrong_payload"
    return status, payload, index


def edit_records(
    source: str,
    sample_index: int,
    original_payload: int,
    sequence: str,
    codec: choice.GCHomopolymerCodec,
    affine: dict[str, int],
) -> tuple[list[dict[str, Any]], list[tuple[int, str, int]]]:
    rows: list[dict[str, Any]] = []
    accepted: list[tuple[int, str, int]] = []

    def add_record(
        edit_type: str,
        position: int,
        original_base: str,
        edited_base: str,
        mutated: str,
    ) -> None:
        status, decoded_payload, decoded_choice = decode_status(
            codec, affine, mutated, original_payload
        )
        row = {
            "selector_source_pool": source,
            "selector_model": MODEL_NAME,
            "selector_bits": SELECTOR_BITS,
            "message_sample_index": sample_index,
            "original_payload_hex": hex(original_payload),
            "original_sequence_sha256": sha256_text(sequence),
            "edit_type": edit_type,
            "edit_position_0_based": position,
            "original_base": original_base,
            "edited_base": edited_base,
            "mutated_length": len(mutated),
            "mutated_sequence_sha256": sha256_text(mutated),
            "decoder_status": status,
            "decoded_payload_hex": "" if decoded_payload is None else hex(decoded_payload),
            "decoded_choice_index": "" if decoded_choice is None else decoded_choice,
            "crc16_payload_valid": ""
            if decoded_payload is None
            else crc16_payload_valid(decoded_payload),
            "canonical_verifier_accepts": "not_evaluated"
            if decoded_payload is not None
            else "not_applicable",
            "combined_crc16_and_canonical_accepts_wrong_payload": "not_evaluated"
            if status == "accepted_wrong_payload"
            else False,
        }
        row_index = len(rows)
        rows.append(row)
        if decoded_payload is not None:
            accepted.append((row_index, mutated, decoded_payload))

    for position, original_base in enumerate(sequence):
        for base in DNA:
            if base == original_base:
                continue
            add_record(
                "substitution",
                position,
                original_base,
                base,
                sequence[:position] + base + sequence[position + 1 :],
            )
    for position in range(len(sequence) + 1):
        for base in DNA:
            add_record(
                "insertion",
                position,
                "",
                base,
                sequence[:position] + base + sequence[position:],
            )
    for position, original_base in enumerate(sequence):
        add_record(
            "deletion",
            position,
            original_base,
            "",
            sequence[:position] + sequence[position + 1 :],
        )
    return rows, accepted


def canonical_sequences_for_payloads(
    payloads: set[int],
    source: str,
    codec: choice.GCHomopolymerCodec,
    affine: dict[str, int],
    score_map: dict[str, float],
) -> dict[int, str]:
    result: dict[int, str] = {}
    for payload in sorted(payloads):
        candidates = [
            sequence_for_logical_rank(
                codec, affine, (payload << SELECTOR_BITS) | candidate_index
            )
            for candidate_index in range(1 << SELECTOR_BITS)
        ]
        keys = deterministic.quantized_score_keys([score_map[value] for value in candidates])
        result[payload] = candidates[int(np.argmax(keys))]
    return result


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(rows)
    output: list[dict[str, Any]] = []
    for (source, edit_type), group in frame.groupby(
        ["selector_source_pool", "edit_type"], sort=True
    ):
        status = group["decoder_status"]
        wrong = status.eq("accepted_wrong_payload")
        same = status.eq("accepted_same_payload")
        crc_pass = group["crc16_payload_valid"].eq(True)
        verifier_pass = group["canonical_verifier_accepts"].eq(True)
        combined = group["combined_crc16_and_canonical_accepts_wrong_payload"].eq(True)
        output.append(
            {
                "selector_source_pool": source,
                "selector_model": MODEL_NAME,
                "selector_bits": SELECTOR_BITS,
                "edit_type": edit_type,
                "emitted_codewords": int(group["message_sample_index"].nunique()),
                "edit_operations": len(group),
                "rejected_length": int(status.eq("rejected_length").sum()),
                "rejected_constraint_or_alphabet": int(
                    status.eq("rejected_constraint_or_alphabet").sum()
                ),
                "rejected_outside_power_of_two_domain": int(
                    status.eq("rejected_outside_power_of_two_domain").sum()
                ),
                "accepted_same_payload": int(same.sum()),
                "accepted_wrong_payload": int(wrong.sum()),
                "wrong_payload_acceptance_fraction_all_edits": float(wrong.mean()),
                "crc16_detected_wrong_payload": int((wrong & ~crc_pass).sum()),
                "crc16_residual_wrong_payload": int((wrong & crc_pass).sum()),
                "canonical_verifier_detected_accepted_corruption": int(
                    ((wrong | same) & ~verifier_pass).sum()
                ),
                "canonical_verifier_residual_wrong_payload": int(
                    (wrong & verifier_pass).sum()
                ),
                "combined_crc16_and_canonical_residual_wrong_payload": int(combined.sum()),
                "inference_boundary": "deterministic single-edit audit; not a stochastic DNA-channel or recovery benchmark",
            }
        )
    return output


def output_manifest(output_dir: Path) -> list[dict[str, Any]]:
    files = sorted(
        path
        for path in output_dir.iterdir()
        if path.is_file() and path.name != "sha256_manifest.tsv"
    )
    files.append(Path(__file__).resolve())
    return [
        {"path": relative(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in files
    ]


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    if args.message_samples != FROZEN_MESSAGE_SAMPLES:
        raise ValueError(f"frozen run requires --message-samples {FROZEN_MESSAGE_SAMPLES}")
    if crc16_ccitt_false(b"123456789") != 0x29B1:
        raise RuntimeError("CRC-16/CCITT-FALSE reference vector failed")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    selection.validate_public_sources()
    _dt, pcr, _pickle_audit = selection.load_frames(args.jobs)
    models, _coefficients = choice.fit_frozen_models(pcr, args.jobs)

    codec = choice.GCHomopolymerCodec(
        choice.GENERATED_LENGTH,
        choice.GENERATED_GC_MIN,
        choice.GENERATED_GC_MAX,
        choice.GENERATED_MAX_HOMOPOLYMER,
    )
    affine = choice.affine_record(codec, FIXED_BITS, namespace_index=0)
    messages = deterministic_messages(args.message_samples)
    payloads = [payload_with_crc16(message) for message in messages]
    if not all(crc16_payload_valid(payload) for payload in payloads):
        raise RuntimeError("a constructed payload failed its CRC-16 contract")

    original_fiber_rows: list[dict[str, Any]] = []
    original_sequences: list[str] = []
    for sample_index, (message, payload) in enumerate(zip(messages, payloads, strict=True)):
        for candidate_index in range(1 << SELECTOR_BITS):
            sequence = sequence_for_logical_rank(
                codec, affine, (payload << SELECTOR_BITS) | candidate_index
            )
            original_sequences.append(sequence)
            original_fiber_rows.append(
                {
                    "message_sample_index": sample_index,
                    "message_hex": hex(message),
                    "payload_hex": hex(payload),
                    "choice_index": candidate_index,
                    "sequence": sequence,
                    "sequence_sha256": sha256_text(sequence),
                }
            )
    original_features, original_scores = feature_scores(original_sequences, models, args.jobs)
    del original_features

    emitted_rows: list[dict[str, Any]] = []
    all_edit_rows: list[dict[str, Any]] = []
    accepted_by_source: dict[str, list[tuple[int, str, int]]] = {source: [] for source in SOURCES}
    for source in SOURCES:
        for sample_index, payload in enumerate(payloads):
            fiber = [
                row
                for row in original_fiber_rows
                if row["message_sample_index"] == sample_index
            ]
            keys = deterministic.quantized_score_keys(
                [original_scores[source][str(row["sequence"])] for row in fiber]
            )
            selected = fiber[int(np.argmax(keys))]
            sequence = str(selected["sequence"])
            emitted_rows.append(
                {
                    "selector_source_pool": source,
                    "selector_model": MODEL_NAME,
                    "selector_bits": SELECTOR_BITS,
                    "message_sample_index": sample_index,
                    "message_bits": MESSAGE_BITS,
                    "crc_bits": CRC_BITS,
                    "payload_bits_after_choice": PAYLOAD_BITS,
                    "message_hex": selected["message_hex"],
                    "payload_hex": selected["payload_hex"],
                    "selected_choice_index": selected["choice_index"],
                    "sequence": sequence,
                    "sequence_sha256": selected["sequence_sha256"],
                    "fullcontext_score": original_scores[source][sequence],
                }
            )
            rows, accepted = edit_records(
                source,
                sample_index,
                payload,
                sequence,
                codec,
                affine,
            )
            offset = len(all_edit_rows)
            all_edit_rows.extend(rows)
            accepted_by_source[source].extend(
                (offset + row_index, mutated, decoded_payload)
                for row_index, mutated, decoded_payload in accepted
            )

    decoded_payload_union = {
        payload
        for accepted in accepted_by_source.values()
        for _row_index, _mutated, payload in accepted
    }
    canonical_candidate_sequences = [
        sequence_for_logical_rank(
            codec, affine, (payload << SELECTOR_BITS) | candidate_index
        )
        for payload in sorted(decoded_payload_union)
        for candidate_index in range(1 << SELECTOR_BITS)
    ]
    _canonical_features, canonical_scores = feature_scores(
        canonical_candidate_sequences, models, args.jobs
    )
    for source in SOURCES:
        canonical = canonical_sequences_for_payloads(
            {payload for _row, _sequence, payload in accepted_by_source[source]},
            source,
            codec,
            affine,
            canonical_scores[source],
        )
        for row_index, mutated, decoded_payload in accepted_by_source[source]:
            verifier_accepts = mutated == canonical[decoded_payload]
            all_edit_rows[row_index]["canonical_verifier_accepts"] = verifier_accepts
            is_wrong = all_edit_rows[row_index]["decoder_status"] == "accepted_wrong_payload"
            crc_pass = all_edit_rows[row_index]["crc16_payload_valid"] is True
            all_edit_rows[row_index][
                "combined_crc16_and_canonical_accepts_wrong_payload"
            ] = bool(is_wrong and crc_pass and verifier_accepts)

    summary_rows = summarize(all_edit_rows)
    total_edits = len(all_edit_rows)
    total_wrong = sum(row["decoder_status"] == "accepted_wrong_payload" for row in all_edit_rows)
    total_crc_residual = sum(
        row["decoder_status"] == "accepted_wrong_payload"
        and row["crc16_payload_valid"] is True
        for row in all_edit_rows
    )
    total_combined_residual = sum(
        row["combined_crc16_and_canonical_accepts_wrong_payload"] is True
        for row in all_edit_rows
    )
    contract_rows = [
        {
            "audit": "single_edit_channel_boundary",
            "language": "exact_108nt_gc49_59_hp3",
            "selector_bits": SELECTOR_BITS,
            "candidate_width": 1 << SELECTOR_BITS,
            "sources": len(SOURCES),
            "messages_per_source": args.message_samples,
            "emitted_codewords": len(emitted_rows),
            "all_substitutions_per_codeword": choice.GENERATED_LENGTH * 3,
            "all_insertions_per_codeword": (choice.GENERATED_LENGTH + 1) * 4,
            "all_deletions_per_codeword": choice.GENERATED_LENGTH,
            "total_edit_operations": total_edits,
            "crc_name": "CRC-16/CCITT-FALSE",
            "crc_polynomial": "0x1021",
            "crc_initial_value": "0xffff",
            "crc_check_vector_123456789": "0x29b1",
            "message_bits": MESSAGE_BITS,
            "crc_bits": CRC_BITS,
            "choice_bits": SELECTOR_BITS,
            "retained_message_rate_bits_per_variable_nt": MESSAGE_BITS
            / choice.GENERATED_LENGTH,
            "wrong_payload_acceptances_without_checksum": total_wrong,
            "wrong_payload_acceptances_after_crc16": total_crc_residual,
            "wrong_payload_acceptances_after_crc16_and_optional_canonical_verifier": total_combined_residual,
            "scope": "error-detection composition audit only; no correction, dropout, clustering, consensus or file recovery",
        }
    ]

    write_table(OUT_DIR / "channel_error_audit_contract.tsv", contract_rows)
    write_table(OUT_DIR / "emitted_crc16_codewords.tsv", emitted_rows)
    write_table(OUT_DIR / "single_edit_outcomes.tsv", all_edit_rows)
    write_table(OUT_DIR / "single_edit_summary.tsv", summary_rows)
    write_table(
        OUT_DIR / "crc16_reference_vectors.tsv",
        [
            {
                "algorithm": "CRC-16/CCITT-FALSE",
                "input_ascii": "123456789",
                "expected_hex": "0x29b1",
                "observed_hex": hex(crc16_ccitt_false(b"123456789")),
                "passed": crc16_ccitt_false(b"123456789") == 0x29B1,
            }
        ],
    )
    write_json(
        OUT_DIR / "environment_and_seeds.json",
        {
            "analysis": "Paper 2 deterministic channel-error boundary audit",
            "command": (
                "$PYTHON analysis_tools/audit_paper2_channel_error_boundary.py "
                f"--jobs {args.jobs} --message-samples {args.message_samples}"
            ),
            "python": sys.version,
            "platform": platform.platform(),
            "packages": {
                name: importlib.metadata.version(name)
                for name in ["numpy", "pandas", "scikit-learn", "joblib"]
            },
            "base_seed": BASE_SEED,
            "message_samples_per_source": args.message_samples,
            "source_models": list(SOURCES),
            "model": MODEL_NAME,
            "fixed_bits": FIXED_BITS,
            "selector_bits": SELECTOR_BITS,
            "crc_bits": CRC_BITS,
            "scope": contract_rows[0]["scope"],
        },
    )
    lines = [
        "# Paper 2 channel-error boundary audit",
        "",
        f"The deterministic audit covered {len(emitted_rows)} FullContext-selected generated codewords and all {total_edits:,} single-base edit operations.",
        f"The noiseless payload decoder accepted {total_wrong:,} edits with a wrong payload; CRC-16/CCITT-FALSE reduced this to {total_crc_residual:,}, and CRC plus the optional source-specific canonical verifier left {total_combined_residual:,} in this sample.",
        "",
        "Insertions and deletions are rejected by the fixed-length interface. This analysis is an error-detection composition audit, not a DNA-channel model, correction experiment, dropout benchmark, consensus analysis or file-recovery result.",
        "",
    ]
    (OUT_DIR / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")
    write_table(OUT_DIR / "sha256_manifest.tsv", output_manifest(OUT_DIR))
    print(f"Channel-error boundary audit complete: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
