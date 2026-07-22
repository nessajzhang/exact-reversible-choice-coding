#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKIP_TEX=0
MODE="full"
ANALYSIS_MODE="public-inputs"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-tex) SKIP_TEX=1 ;;
    --analysis-only) MODE="analysis" ;;
    --figures-only) MODE="figures" ;;
    --from-frozen) ANALYSIS_MODE="frozen" ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
  shift
done

# Verify the immutable release before any reproduction command updates
# environment records, figures or compiled artifacts.
"${SCRIPT_DIR}/verify_integrity.sh"

if [[ "${MODE}" == "analysis" ]]; then
  if [[ "${ANALYSIS_MODE}" == "frozen" ]]; then
    "${SCRIPT_DIR}/reproduce_analysis.sh" --from-frozen
  else
    "${SCRIPT_DIR}/reproduce_analysis.sh"
  fi
  source "${SCRIPT_DIR}/paper2_repro_common.sh"
  paper2_prepare_environment
  cd "${PAPER2_ROOT}"
  "${PAPER2_PYTHON}" analysis_tools/verify_paper2_output_manifests.py
  exit 0
fi
if [[ "${MODE}" == "figures" ]]; then
  "${SCRIPT_DIR}/reproduce_figures.sh"
  exit 0
fi

if [[ "${ANALYSIS_MODE}" == "frozen" ]]; then
  "${SCRIPT_DIR}/reproduce_analysis.sh" --from-frozen
else
  "${SCRIPT_DIR}/reproduce_analysis.sh"
fi
"${SCRIPT_DIR}/reproduce_figures.sh"

TEX_COMPILED=0
if [[ "${SKIP_TEX}" == "1" ]]; then
  printf '%s\n' 'Analysis and figure reproduction passed; TeX compilation skipped by --skip-tex.'
elif command -v latexmk >/dev/null 2>&1 && command -v bibtex >/dev/null 2>&1; then
  "${SCRIPT_DIR}/compile_manuscript.sh" --strict
  TEX_COMPILED=1
else
  "${SCRIPT_DIR}/compile_manuscript.sh"
fi

source "${SCRIPT_DIR}/paper2_repro_common.sh"
paper2_prepare_environment
cd "${PAPER2_ROOT}"
"${PAPER2_PYTHON}" analysis_tools/verify_paper2_output_manifests.py
if [[ "${TEX_COMPILED}" == "1" ]]; then
  "${PAPER2_PYTHON}" analysis_tools/verify_paper2_bioinformatics_choice_consistency.py
elif [[ "${SKIP_TEX}" != "1" ]]; then
  printf '%s\n' 'Analysis reproduction passed; TeX compilation was unavailable, so PDF consistency checks were skipped.'
fi
printf '%s\n' 'Paper 2 full audit completed.'
