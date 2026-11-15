#!/usr/bin/env bash
# A minimal Codex Cloud preflight; does not run expensive mathematical verifiers.
set -euo pipefail

case "${1:-}" in
  ""|--tex) ;;
  *) echo "Usage: bash scripts/cloud_preflight.sh [--tex]" >&2; exit 2 ;;
esac

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

command -v python3 >/dev/null || { echo "Missing python3" >&2; exit 1; }
python3 - <<'PY'
import sys
import sympy
print("Python:", sys.version.split()[0])
print("SymPy:", sympy.__version__)
if sympy.__version__ != "1.14.0":
    raise SystemExit("Expected sympy==1.14.0 from requirements.txt")
PY

python3 scripts/check_project.py

if [[ "${1:-}" == "--tex" ]]; then
  for program in latexmk lualatex kpsewhich; do
    command -v "$program" >/dev/null || {
      echo "Missing $program; install a LuaLaTeX toolchain in the Cloud environment" >&2
      exit 1
    }
  done
  for texfile in ltjsarticle.cls luatexja-preset.sty; do
    [[ -n "$(kpsewhich "$texfile")" ]] || {
      echo "Missing $texfile (LuaTeX-ja / Japanese TeX packages)" >&2
      exit 1
    }
  done
  latexmk research/plan_ja.tex
fi

echo "PASS: cloud preflight (scope: environment and structural checks only)"
