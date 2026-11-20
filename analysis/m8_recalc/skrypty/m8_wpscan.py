#!/usr/bin/env python3
"""M8_przeliczenie pkt 1 — skan punktu pracy drzewa H dla D0 modelem M8 (FITBIAS, bez zmian) obok M0. POST HOC.
Wejście: WPSCAN wps_hist.npz (tylko odczyt). Moduły fitu: FITBIAS/skrypty/fb_fit.py, fb_var.py (import, bez zmian).
Na (run, bin, ε): fit M0 i M8 (κ_MC × 1,364) repliki 0 + 30 replik bootstrapu po zdarzeniach (każda z HESSE → błąd stosunku błędów).
Pass/fail (A = ε ≥ ε_t, B = WP70 ∧ ¬A): fit łączny, E = N_A + N_B, repliki.
  m8_wpscan.py [--nproc 7]  → out/wpscan_m8.json"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FB = os.path.join(HERE, "..", "..", "FITBIAS_2026-11-19", "skrypty")
sys.path.insert(0, FB)
sys.dont_write_bytecode = True
import fb_fit as F  # noqa: E402
import fb_var as V  # noqa: E402

KSC = 1.364
MODELS = ("M0", "M8")
SP = "D0"


def cat_s(*a, **k):
    c = V.wcat(*a, **k)
    c["kap_mc"] = c["kap"]
    c["kap"] = c["kap"] * KSC
    return c


def task(t):
    bn, run = t
    ed = V.INFO["eps_desc"]
    res = {"bin": bn, "run": run, "S": {}, "J": {}}
    for e in ed:
        k = cat_s(SP, bn, e)
        D = V.wdata(k, V.SRC[run])
        o = {"eps": k["e"], "kap": k["kap"], "frefl": k["frefl"], "n_rows": float(D[0].sum())}
        for m in MODELS:
            f = V.fit_set(D[0][None], [k], SP, m)
            r = V.sres(f)
            r["p"] = F.pval(f["chi2"], f["ndf"])
            bN, bE = [], []
            for i in range(1, D.shape[0]):
                fb = V.fit_set(D[i][None], [k], SP, m, start=f["x"], hesse=True)
                ok = fb["valid"] and fb.get("err") is not None
                bN.append(fb["N"][0] if ok else None)
                bE.append(fb["err"][0] if ok else None)
            r["boot_N"], r["boot_err"] = bN, bE
            o[m] = r
        res["S"][str(e)] = o
    for et in (0.5, 0.4, 0.35):
        a, b = cat_s(SP, bn, et), cat_s(SP, bn, 0.7, et)
        DA, DB = V.wdata(a, V.SRC[run]), V.wdata(b, V.SRC[run])
        o = {"f_mc": a["e"] / (a["e"] + b["e"])}
        for m in MODELS:
            f = V.fit_set(np.stack([DA[0], DB[0]]), [a, b], SP, m)
            r = V.jres(f)
            r["p"] = F.pval(f["chi2"], f["ndf"])
            bf, bE = [], []
            for i in range(1, DA.shape[0]):
                fb = V.jres(V.fit_set(np.stack([DA[i], DB[i]]), [a, b], SP, m, start=f["x"], hesse=False))
                ok = fb["valid"] and fb["f"] is not None
                bf.append(fb["f"] if ok else None)
                bE.append(fb["E"] if ok else None)
            r["boot_f"], r["boot_E"] = bf, bE
            o[m] = r
        res["J"][str(et)] = o
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()
    tasks = [(bn, run) for bn in F.BINS for run in V.SRC]
    out = []
    with Pool(a.nproc, initializer=V._init) as p:
        for i, r in enumerate(p.imap_unordered(task, tasks)):
            out.append(r)
            print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(tasks)}", flush=True)
    out.sort(key=lambda r: (F.BINS.index(r["bin"]), r["run"]))
    fn = os.path.join(HERE, "..", "out", "wpscan_m8.json")
    json.dump({"res": out, "kscale": KSC, "sec": time.time() - t0}, open(fn, "w"))
    print("zapisano", fn)


if __name__ == "__main__":
    main()
