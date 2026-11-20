#!/usr/bin/env python3
"""MC_templates/analyze.py — wrapper: uruchamia Q_lxplus/mctmpl/mctmpl_analyze.py (jedno źródło kodu, to samo co w węźle merge na condorze).
Użycie: .venv_h5/bin/python charm_work/MC_templates/analyze.py [--in MC_templates/runs/lxplus] [--out MC_templates/out] [--min-n X]
Lokalnie 1 rdzeń, kilka minut; gotowy wynik z condora leży też w runs/lxplus/ana/ (WYNIKI.md, results.json) po fetch."""
import os
import runpy
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if "--in" not in sys.argv:
    sys.argv += ["--in", os.path.join(HERE, "runs", "lxplus")]
if "--out" not in sys.argv:
    sys.argv += ["--out", os.path.join(HERE, "out")]
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
runpy.run_path(os.path.join(HERE, "..", "Q_lxplus", "mctmpl", "mctmpl_analyze.py"), run_name="__main__")
