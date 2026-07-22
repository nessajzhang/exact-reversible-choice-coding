#!/usr/bin/env bash

paper2_python_matches_frozen_environment() {
  local candidate_python="$1"
  "${candidate_python}" \
    "${PAPER2_ROOT}/analysis_tools/verify_paper2_environment.py" --quiet \
    >/dev/null 2>&1
}

paper2_prepare_environment() {
  PAPER2_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  PAPER2_ROOT="$(cd "${PAPER2_SCRIPT_DIR}/../.." && pwd)"
  export PAPER2_SCRIPT_DIR PAPER2_ROOT
  # Verification and reproduction must not mutate a clean release tree with
  # interpreter caches. This applies to every Python process launched below.
  export PYTHONDONTWRITEBYTECODE=1
  export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1784332800}"
  export FORCE_SOURCE_DATE="${FORCE_SOURCE_DATE:-1}"

  local requested_python="${PYTHON:-python3}"
  if paper2_python_matches_frozen_environment "${requested_python}"; then
    PAPER2_PYTHON="${requested_python}"
  elif command -v uv >/dev/null 2>&1; then
    local env_dir="${PAPER2_UV_ENV:-${HOME}/.cache/paper2-thermodynamic-risk-coding/py312}"
    if [[ ! -x "${env_dir}/bin/python" ]] || \
       ! "${env_dir}/bin/python" -c 'import sys; raise SystemExit(sys.version_info[:3] != (3, 12, 13))' >/dev/null 2>&1; then
      mkdir -p "$(dirname "${env_dir}")"
      uv venv --clear --python 3.12.13 "${env_dir}"
    fi
    if ! paper2_python_matches_frozen_environment "${env_dir}/bin/python"; then
      uv pip install --python "${env_dir}/bin/python" \
        joblib==1.5.3 matplotlib==3.11.0 numpy==2.5.0 pandas==3.0.3 \
        scikit-learn==1.9.0 scipy==1.18.0
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
