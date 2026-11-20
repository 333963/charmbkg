#!/bin/bash
# PEAK_lxy010 v2 — odtworzenie (18.11.2026 wieczór). Bieg lokalny. ~10 min na 4 rdzeniach (bootstrap 200).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); PY=/home/k/Documents/charm/.venv_h5/bin/python
mkdir -p "$HERE/out" "$HERE/logs"; cd "$HERE"; export PYTHONDONTWRITEBYTECODE=1
"$PY" check_inputs.py > out/check_inputs.txt 2>&1                 # kontrola wejść (wiersze, pliki, fold, klucz fd_nd)
"$PY" fit_v2.py fd_nd "${NBOOT:-200}" > logs/fit_fd_nd.log 2>&1   # kształty MC per cięcie L_xy (musi być pierwsze)
for s in mb_nd data_427929 data_428770; do "$PY" fit_v2.py $s "${NBOOT:-200}" > logs/fit_$s.log 2>&1 & done
"$PY" refl_v2.py "${NBOOT:-200}" > logs/refl.log 2>&1 &
wait
"$PY" mb_truth.py > logs/mb_truth.log 2>&1
"$PY" summarize_v2.py > logs/summarize.log 2>&1                   # → out/tables_v2.md, out/summary_v2.json
echo OK
