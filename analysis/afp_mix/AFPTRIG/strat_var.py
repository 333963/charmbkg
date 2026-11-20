#!/usr/bin/env python3
"""AFPTRIG/strat_var.py — wrażliwość na ziarnistość klas nSelBL (nominał: lib.NSEL_EDGES). → strat_var.json"""
import json, numpy as np, lib
from analysis import Run, nominal_pool
import analysis
ED = {"coarse": [0, 10, 20, 30, 50, 10000], "nominal": list(lib.NSEL_EDGES), "fine": list(range(0, 60, 3)) + [60, 70, 80, 100, 10000],
      "vfine": list(range(0, 80, 2)) + [80, 90, 100, 120, 10000]}
OUT = {}
for run in lib.RUNS:
    R = Run(run); O = OUT[str(run)] = {}
    for stream in ("afp_any", "pt4_L1AFP", "noalg_L1AFP", "nonAFP"):
        pool = nominal_pool(run, stream); O[stream] = {}
        for nm, e in ED.items():
            e = np.array(e); analysis.NCLS = len(e) - 1
            R.cls = np.clip(np.searchsorted(e, R.d["nsel"], side="right") - 1, 0, analysis.NCLS - 1)
            r = [R.est(s, stream, pool) for s in lib.SIDES]
            O[stream][nm] = {"n_acc": sum(x["n_acc"] for x in r), "err": float(np.hypot(r[0]["err"], r[1]["err"])), "n_cr_nopool": sum(x["n_cr_nopool"] for x in r)}
        print(run, stream, {k: f"{v['n_acc']:.1f}±{v['err']:.1f} (bez puli {v['n_cr_nopool']:.0f})" for k, v in O[stream].items()}, flush=True)
        R.cache.clear()
json.dump(OUT, open("strat_var.json", "w"), indent=1)
