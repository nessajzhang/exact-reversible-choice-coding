#!/usr/bin/env bash

paper2_prepare_environment() {
  PAPER2_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  PAPER2_ROOT="$(cd "${PAPER2_SCRIPT_DIR}/../.." && pwd)"
  export PAPER2_SCRIPT_DIR PAPER2_ROOT
  export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1784332800}"
  export FORCE_SOURCE_DATE="${FORCE_SOURCE_DATE:-1}"

  local requested_python="${PYTHON:-python3}"
  if "${requested_python}" -c 'import joblib, matplotlib, numpy, pandas, scipy, sklearn' >/dev/null 2>&1; then
    PAPER2_PYTHON="${requested_python}"
  elif command -v uv >/dev/null 2>&1; then
    local env_dir="${PAPER2_UV_ENV:-${HOME}/.cache/paper2-thermodynamic-risk-coding/py312}"
    if [[ ! -x "${env_dir}/bin/python" ]]; then
      mkdir -p "$(dirname "${env_dir}")"
      uv venv --python 3.12 "${env_dir}"
    fi
    if ! "${env_dir}/bin/python" -c 'import joblib, matplotlib, numpy, pandas, scipy, sklearn' >/dev/null 2>&1; then
      uv pip install --python "${env_dir}/bin/python" \
        joblib==1.5.3 matplotlib==3.11.0 numpy==2.5.0 pandas==3.0.3 \
        scikit-learn==1.9.0 scipy==1.18.0
    fi
    PAPER2_PYTHON="${env_dir}/bin/python"
  else
    printf '%s\n' "${requested_python} lacks the required packages and uv is unavailable." >&2
    return 2
  fi
  export PAPER2_PYTHON
}
