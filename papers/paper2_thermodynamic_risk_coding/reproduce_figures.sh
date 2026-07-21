#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/paper2_repro_common.sh"
paper2_prepare_environment
cd "${PAPER2_ROOT}"

"${PAPER2_PYTHON}" analysis_tools/make_paper2_assay_calibrated_selection_figure.py
"${PAPER2_PYTHON}" analysis_tools/plot_paper2_reversible_choice_codec.py
printf '%s\n' 'Paper 2 figure reproduction passed.'
