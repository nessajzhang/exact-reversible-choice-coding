#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/paper2_repro_common.sh"
paper2_prepare_environment
cd "${PAPER2_ROOT}"

"${PAPER2_PYTHON}" analysis_tools/download_paper2_public_inputs.py
"${SCRIPT_DIR}/run_full_audit.sh" "$@"
