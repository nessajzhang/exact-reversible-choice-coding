#!/usr/bin/env bash

paper2_initialize_context() {
  PAPER2_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  PAPER2_ROOT="$(cd "${PAPER2_SCRIPT_DIR}/../.." && pwd)"
  PAPER2_CANONICAL_PYTHON_VERSION="3.12.10"
  export PAPER2_SCRIPT_DIR PAPER2_ROOT PAPER2_CANONICAL_PYTHON_VERSION
  # Verification and reproduction must not mutate a clean release tree with
  # interpreter caches. This applies to every Python process launched below.
  export PYTHONDONTWRITEBYTECODE=1
  export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1784332800}"
  export FORCE_SOURCE_DATE="${FORCE_SOURCE_DATE:-1}"
}

paper2_python_matches_frozen_environment() {
  local candidate_python="$1"
  "${candidate_python}" \
    "${PAPER2_ROOT}/analysis_tools/verify_paper2_environment.py" --quiet \
    >/dev/null 2>&1
}

paper2_python_has_contract_dependencies() {
  local candidate_python="$1"
  "${candidate_python}" -c \
    'import joblib, matplotlib, numpy, pandas, scipy, sklearn' \
    >/dev/null 2>&1
}

paper2_prepare_environment() {
  paper2_initialize_context

  local requested_python="${PYTHON:-python3}"
  if paper2_python_matches_frozen_environment "${requested_python}"; then
    PAPER2_PYTHON="${requested_python}"
  elif command -v uv >/dev/null 2>&1; then
    local env_dir="${PAPER2_UV_ENV:-${HOME}/.cache/paper2-thermodynamic-risk-coding/py31210}"
    if [[ ! -x "${env_dir}/bin/python" ]] || \
       ! "${env_dir}/bin/python" -c 'import sys; raise SystemExit(sys.version_info[:3] != (3, 12, 10))' >/dev/null 2>&1; then
      mkdir -p "$(dirname "${env_dir}")"
      # Python 3.12.10 is the last 3.12 maintenance release with regular
      # binary installers and is available through uv's managed-Python path.
      # Install it explicitly so the bootstrap does not depend on implicit
      # interpreter-download behavior changing across uv releases.
      uv python install "${PAPER2_CANONICAL_PYTHON_VERSION}"
      uv venv --clear --python "${PAPER2_CANONICAL_PYTHON_VERSION}" "${env_dir}"
    fi
    if ! paper2_python_matches_frozen_environment "${env_dir}/bin/python"; then
      uv pip install --python "${env_dir}/bin/python" \
        -r "${PAPER2_SCRIPT_DIR}/requirements.txt"
    fi
    if ! paper2_python_matches_frozen_environment "${env_dir}/bin/python"; then
      "${env_dir}/bin/python" \
        "${PAPER2_ROOT}/analysis_tools/verify_paper2_environment.py" || true
      printf '%s\n' 'Unable to establish the frozen Paper 2 Python environment.' >&2
      return 2
    fi
    PAPER2_PYTHON="${env_dir}/bin/python"
  else
    if command -v "${requested_python}" >/dev/null 2>&1; then
      "${requested_python}" \
        "${PAPER2_ROOT}/analysis_tools/verify_paper2_environment.py" || true
    fi
    printf '%s\n' "${requested_python} does not match the frozen environment and uv is unavailable." >&2
    return 2
  fi
  PAPER2_ENVIRONMENT_MODE="canonical-frozen"
  export PAPER2_PYTHON PAPER2_ENVIRONMENT_MODE
}

paper2_prepare_contract_environment() {
  paper2_initialize_context
  local requested_python="${PAPER2_TEST_PYTHON:-${PYTHON:-python3}}"
  if paper2_python_has_contract_dependencies "${requested_python}"; then
    PAPER2_CONTRACT_PYTHON="${requested_python}"
    PAPER2_CONTRACT_ENVIRONMENT_MODE="provided-compatible"
  else
    paper2_prepare_environment
    PAPER2_CONTRACT_PYTHON="${PAPER2_PYTHON}"
    PAPER2_CONTRACT_ENVIRONMENT_MODE="canonical-fallback"
  fi
  export PAPER2_CONTRACT_PYTHON PAPER2_CONTRACT_ENVIRONMENT_MODE
}
