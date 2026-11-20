#!/usr/bin/env python3
"""AFPTRIG/leak.py — przeciek własnych protonów nad cięcie: nadwyżka danych nad szablonem (normalizacja w Δ > 0,10)
w 0,02 < Δ ≤ 0,045 (przerwa) i 0,045 < Δ ≤ 0,10 (dolna część CR), jako ułamek N_tag(SR), w przedziałach μ. → leak.json"""
import json, numpy as np, lib
from analysis import Run, nominal_pool, INF
OUT = {}
for run in lib.RUNS:
    R = Run(run); O = OUT[str(run)] = {}
    for stream in ("afp_any", "pt4_L1AFP", "nonAFP"):
        pool = nominal_pool(run, stream); rows = []
        mu = R.d["mu"]; ref = R.S[stream] & R.has & (R.d["tag_A"] | R.d["tag_C"])
        e = np.quantile(mu[ref], [0, 1 / 3, 2 / 3, 1]); e[-1] += 1e-6
        for lo, hi in [(e[0], e[-1])] + list(zip(e[:-1], e[1:])):
            w = (R.has & (mu >= lo) & (mu < hi)).astype(float); r = {"mu_lo": float(lo), "mu_hi": float(hi), "mu_mean": float(mu[ref & (mu >= lo) & (mu < hi)].mean())}
            for nm, sr in (("gap", (0.02, 0.045)), ("crlow", (0.045, 0.10))):
                o = p = v = 0.0
                for s in lib.SIDES:
                    x = R.est(s, stream, pool, sr=sr, cr=(0.10, INF), w=w); o += x["n_sr"]; p += x["n_acc"]; v += x["err"] ** 2 + x["n_sr"]
                r[nm] = {"obs": o, "pred": p, "excess": o - p, "excess_err": float(np.sqrt(v))}
            nsr = sum(R.est(s, stream, pool, w=w)["n_sr"] for s in lib.SIDES); nacc = sum(R.est(s, stream, pool, cr=(0.10, INF), w=w)["n_acc"] for s in lib.SIDES)
            r["n_sr"] = nsr; r["n_acc_far"] = nacc
            for nm in ("gap", "crlow"):
                r[nm]["frac_of_true"] = r[nm]["excess"] / (nsr - nacc); r[nm]["frac_err"] = r[nm]["excess_err"] / (nsr - nacc)
            rows.append(r)
            print(run, stream, f"μ {lo:.4f}-{hi:.4f} <{r['mu_mean']:.4f}> N_SR {nsr:.0f} gap {100*r['gap']['frac_of_true']:.2f}±{100*r['gap']['frac_err']:.2f} %  CRlow {100*r['crlow']['frac_of_true']:.2f}±{100*r['crlow']['frac_err']:.2f} %", flush=True)
        O[stream] = rows; R.cache.clear()
json.dump(OUT, open("leak.json", "w"), indent=1)
