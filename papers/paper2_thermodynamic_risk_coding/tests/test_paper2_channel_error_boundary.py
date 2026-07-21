#!/usr/bin/env python3
"""Unit tests for the deterministic single-edit/CRC boundary audit."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analysis_tools"))

import audit_paper2_channel_error_boundary as audit  # noqa: E402
import validate_paper2_reversible_choice_codec as choice  # noqa: E402


class ChannelErrorBoundaryTests(unittest.TestCase):
    def test_crc16_reference_vector(self) -> None:
        self.assertEqual(audit.crc16_ccitt_false(b"123456789"), 0x29B1)

    def test_crc_payload_construction(self) -> None:
        for message in [0, 1, (1 << audit.MESSAGE_BITS) // 2, (1 << audit.MESSAGE_BITS) - 1]:
            payload = audit.payload_with_crc16(message)
            self.assertTrue(audit.crc16_payload_valid(payload))
            self.assertEqual(payload >> audit.CRC_BITS, message)

    def test_fixed_length_and_noiseless_decode_contract(self) -> None:
        codec = choice.GCHomopolymerCodec(
            choice.GENERATED_LENGTH,
            choice.GENERATED_GC_MIN,
            choice.GENERATED_GC_MAX,
            choice.GENERATED_MAX_HOMOPOLYMER,
        )
        affine = choice.affine_record(codec, audit.FIXED_BITS, namespace_index=0)
        payload = audit.payload_with_crc16(0)
        sequence = audit.sequence_for_logical_rank(
            codec, affine, payload << audit.SELECTOR_BITS
        )
        status, decoded_payload, decoded_choice = audit.decode_status(
            codec, affine, sequence, payload
        )
        self.assertEqual(status, "accepted_same_payload")
        self.assertEqual(decoded_payload, payload)
        self.assertEqual(decoded_choice, 0)
        deleted = sequence[:-1]
        self.assertEqual(
            audit.decode_status(codec, affine, deleted, payload)[0],
            "rejected_length",
        )


if __name__ == "__main__":
    unittest.main()
