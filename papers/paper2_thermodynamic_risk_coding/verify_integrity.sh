#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="all"
if [[ "${1:-}" == "--files-only" ]]; then
  MODE="files"
elif [[ "${1:-}" == "--contracts-only" ]]; then
  MODE="contracts"
elif [[ $# -gt 0 ]]; then
  printf 'Unknown option: %s\n' "$1" >&2
  exit 2
fi

if [[ "${MODE}" != "contracts" ]]; then
  "${SCRIPT_DIR}/verify_file_integrity.sh"
fi
if [[ "${MODE}" != "files" ]]; then
  "${SCRIPT_DIR}/verify_contracts.sh"
fi
printf 'Paper 2 integrity verification passed (%s).\n' "${MODE}"
