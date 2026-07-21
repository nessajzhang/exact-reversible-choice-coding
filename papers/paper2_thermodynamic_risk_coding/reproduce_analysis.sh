#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/paper2_repro_common.sh"
paper2_prepare_environment
cd "${PAPER2_ROOT}"

MODE="public-inputs"
if [[ "${1:-}" == "--from-frozen" ]]; then
  MODE="frozen"
elif [[ $# -gt 0 ]]; then
  printf 'Unknown option: %s\n' "$1" >&2
  exit 2
fi

JOBS="${JOBS:-12}"
SEQ_JOBS="${SEQ_JOBS:-10}"
"${PAPER2_PYTHON}" -m unittest \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_deterministic_selection.py \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_grouped_bootstrap.py \
  papers/paper2_thermodynamic_risk_coding/tests/test_paper2_channel_error_boundary.py -v

if [[ "${MODE}" == "public-inputs" ]]; then
  "${PAPER2_PYTHON}" analysis_tools/validate_paper2_public_experimental_data.py --jobs "${JOBS}" --bootstrap 2000
  "${PAPER2_PYTHON}" analysis_tools/validate_paper2_assay_calibrated_selection.py --jobs "${JOBS}" --bootstrap 2000
  "${PAPER2_PYTHON}" analysis_tools/audit_paper2_pcr_sequence_independence.py --jobs "${SEQ_JOBS}" --block-size 64
  "${PAPER2_PYTHON}" analysis_tools/audit_paper2_source_external_independence.py --jobs "${SEQ_JOBS}" --block-size 64
  "${PAPER2_PYTHON}" analysis_tools/validate_paper2_reversible_choice_codec.py --jobs "${JOBS}" --bootstrap 2000 --generated-payloads 512
  "${PAPER2_PYTHON}" analysis_tools/validate_paper2_sota_external_uncertainty.py --jobs "${JOBS}" --bootstrap 2000 --two-stage 2000 --randomizations 100000
  "${PAPER2_PYTHON}" analysis_tools/audit_paper2_channel_error_boundary.py --jobs "${JOBS}" --message-samples 32
fi

"${PAPER2_PYTHON}" analysis_tools/analyze_paper2_external_mapping_sensitivity.py --jobs "${JOBS}"
"${PAPER2_PYTHON}" analysis_tools/analyze_paper2_major_revision_diagnostics.py

if [[ "${MODE}" == "public-inputs" ]]; then
  LOCAL_RUNTIME_DIR="${PAPER2_LOCAL_RUNTIME_DIR:-${SCRIPT_DIR}/bioinformatics_reframe/runtime_benchmark/latest_local}"
  "${PAPER2_PYTHON}" analysis_tools/benchmark_paper2_reversible_choice_codec.py \
    --payloads 128 --repeats 3 --rank-checks 4096 --output-dir "${LOCAL_RUNTIME_DIR}"
fi

printf 'Paper 2 analysis reproduction passed (%s mode).\n' "${MODE}"
