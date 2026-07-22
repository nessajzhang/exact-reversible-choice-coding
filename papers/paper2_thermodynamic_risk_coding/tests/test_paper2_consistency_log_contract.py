#!/usr/bin/env python3
"""Tests for PASS/SKIP/FAIL handling of optional manuscript build logs."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "analysis_tools"))

import verify_paper2_bioinformatics_choice_consistency as consistency  # noqa: E402


class CompilationLogContractTests(unittest.TestCase):
    def test_release_excluded_logs_are_explicit_skips(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            checks = consistency.audit_latex_logs(Path(directory))
        self.assertEqual([row["status"] for row in checks], ["SKIP", "SKIP"])
        self.assertTrue(all("release archive" in row["detail"] for row in checks))

    def test_present_clean_logs_pass(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            build = Path(directory)
            for name in consistency.LATEX_LOG_NAMES:
                (build / name).write_text(
                    "Output written on manuscript.pdf.\n", encoding="utf-8"
                )
            checks = consistency.audit_latex_logs(build)
        self.assertEqual([row["status"] for row in checks], ["PASS", "PASS"])

    def test_present_error_log_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            build = Path(directory)
            (build / "main.log").write_text(
                "! LaTeX Error: broken input.\n", encoding="utf-8"
            )
            (build / "supplementary_codec_evidence.log").write_text(
                "Output written on supplement.pdf.\n", encoding="utf-8"
            )
            checks = consistency.audit_latex_logs(build)
        self.assertEqual([row["status"] for row in checks], ["FAIL", "PASS"])


if __name__ == "__main__":
    unittest.main()
