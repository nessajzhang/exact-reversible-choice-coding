#!/usr/bin/env python3

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analysis_tools"))

import build_paper2_public_release as release  # noqa: E402


class ReleaseHygieneTests(unittest.TestCase):
    def test_clean_release_tree_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "analysis.py").write_text("value = 1\n", encoding="utf-8")
            release.check_no_transient_artifacts(root)

    def test_python_cache_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "analysis_tools" / "__pycache__"
            cache.mkdir(parents=True)
            (cache / "module.cpython-312.pyc").write_bytes(b"cache")
            with self.assertRaisesRegex(RuntimeError, "__pycache__"):
                release.check_no_transient_artifacts(root)

    def test_public_runner_does_not_create_python_cache(self) -> None:
        common = (
            ROOT
            / "papers"
            / "paper2_thermodynamic_risk_coding"
            / "paper2_repro_common.sh"
        )
        with tempfile.TemporaryDirectory() as directory:
            probe_dir = Path(directory)
            (probe_dir / "paper2_cache_probe.py").write_text(
                "VALUE = 1\n", encoding="utf-8"
            )
            env = os.environ.copy()
            env["PYTHONPATH"] = str(probe_dir)
            subprocess.run(
                [
                    "bash",
                    "-c",
                    'source "$1"; paper2_initialize_context; '
                    '"$2" -c "import paper2_cache_probe"',
                    "paper2-release-hygiene-test",
                    str(common),
                    sys.executable,
                ],
                cwd=ROOT,
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertFalse((probe_dir / "__pycache__").exists())


if __name__ == "__main__":
    unittest.main()
