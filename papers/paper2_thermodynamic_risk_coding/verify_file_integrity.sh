#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/paper2_repro_common.sh"
paper2_initialize_context
cd "${PAPER2_ROOT}"

VERIFY_PYTHON="${PAPER2_VERIFY_PYTHON:-${PYTHON:-}}"
if [[ -n "${VERIFY_PYTHON}" ]] && ! command -v "${VERIFY_PYTHON}" >/dev/null 2>&1; then
  printf 'Ignoring unavailable manifest interpreter: %s\n' "${VERIFY_PYTHON}" >&2
  VERIFY_PYTHON=""
fi
if [[ -z "${VERIFY_PYTHON}" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    VERIFY_PYTHON="python3"
  elif command -v python >/dev/null 2>&1; then
    VERIFY_PYTHON="python"
  else
    printf '%s\n' 'A basic Python interpreter is required for the standard-library manifest parser.' >&2
    exit 2
  fi
fi

# This verifier uses only the Python standard library; it does not require the
# canonical scientific environment or any third-party package.
"${VERIFY_PYTHON}" analysis_tools/verify_paper2_output_manifests.py

if [[ -f SHA256SUMS.txt ]]; then
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum -c SHA256SUMS.txt
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 -c SHA256SUMS.txt
  else
    printf '%s\n' 'Neither sha256sum nor shasum is available.' >&2
    exit 2
  fi
  if find . \( -type d -name '__pycache__' -o -type f \( \
      -name '*.pyc' -o -name '*.pyo' -o -name '.DS_Store' -o -name 'Thumbs.db' \
    \) \) -print -quit | grep -q .; then
    printf '%s\n' 'Transient interpreter/editor/OS artifact found in release tree.' >&2
    exit 3
  fi
  printf '%s\n' 'Paper 2 root SHA-256 and release-hygiene checks passed.'
else
  printf '%s\n' 'SKIP_ROOT_SHA256: source workspace has no release-level SHA256SUMS.txt.'
fi
