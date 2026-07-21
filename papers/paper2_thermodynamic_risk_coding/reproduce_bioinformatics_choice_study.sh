#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARGS=("$@")
if [[ "${PAPER2_SKIP_UPSTREAM:-0}" == "1" ]]; then
  ARGS+=("--from-frozen")
fi
exec "${SCRIPT_DIR}/run_full_audit.sh" "${ARGS[@]}"
