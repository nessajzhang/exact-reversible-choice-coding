#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STRICT=0
if [[ "${1:-}" == "--strict" ]]; then
  STRICT=1
elif [[ $# -gt 0 ]]; then
  printf 'Unknown option: %s\n' "$1" >&2
  exit 2
fi

if ! command -v latexmk >/dev/null 2>&1 || ! command -v bibtex >/dev/null 2>&1; then
  printf '%s\n' 'TeX compilation skipped: latexmk and BibTeX are both required.' >&2
  if [[ "${STRICT}" == "1" ]]; then exit 2; else exit 0; fi
fi

mkdir -p "${SCRIPT_DIR}/build"
(
  cd "${SCRIPT_DIR}"
  latexmk -norc -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
  latexmk -norc -pdf -interaction=nonstopmode -halt-on-error -outdir=build supplementary_codec_evidence.tex
)
printf '%s\n' 'Paper 2 manuscript and Supplement compilation passed.'
