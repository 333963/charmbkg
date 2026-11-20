#!/bin/bash
# PEAK_lxy010 — odtworzenie obliczeń lokalnych (18.11.2026). Nic nie wysyła, nic nie pobiera. ~12 min (z bootstrapem 200).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); W=$(cd "$HERE/.." && pwd); PY=$W/../.venv_h5/bin/python
mkdir -p "$HERE/out" "$HERE/logs"; cd "$HERE"
# 1) D0 z tagiem D* z lokalnych h5 v2.0 runu 427929 (226/458 plików; L_xy ≥ 0,10), score D0_H z modeli foldów — ten sam moduł co etap hyp_dstar
if [ ! -s out/dstar_data_427929_v20.parquet ]; then
  PYTHONPATH=$W/Q_lxplus/hyp_extract/src PYTHONDONTWRITEBYTECODE=1 "$PY" -u -m charmtag.hyp_extract --sample data_427929 \
    --h5 "$W/DLZ_dplus/DATA_2026-11-17/h5/data_427929/*_none.h5" --out out/dstar_data_427929_v20.parquet \
    --models "$W/../runs/lxplus/models/D0_H" --score-mode model --dstar-only > logs/extract_427929.log 2>&1
fi
export PYTHONDONTWRITEBYTECODE=1
"$PY" fit_lxy.py "${NBOOT:-200}" > logs/fit_lxy.log 2>&1      # → out/lxy_rows.csv, out/lxy_boot.csv
"$PY" lxy_trend.py > logs/lxy_trend.log 2>&1                  # → out/lxy_trend.csv
"$PY" refl_lxy.py > logs/refl_lxy.log 2>&1                    # → out/refl_rows.csv, out/refl_shape.csv
"$PY" summarize.py > logs/summarize.log 2>&1                  # → out/tables.md, out/summary.json
echo OK
