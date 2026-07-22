#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/paper2_repro_common.sh"
paper2_prepare_environment
cd "${PAPER2_ROOT}"

"${PAPER2_PYTHON}" analysis_tools/verify_paper2_output_manifests.py
"${PAPER2_PYTHON}" -m unittest \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_deterministic_selection.py \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_grouped_bootstrap.py \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_channel_error_boundary.py \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_consistency_log_contract.py \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_environment_contract.py \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_release_hygiene.py
printf '%s\n' 'Paper 2 integrity and deterministic-contract checks passed.'
