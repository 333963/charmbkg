#!/usr/bin/env python3
"""FITBIAS_Dplus — zabawki przy czystości DANYCH (KRYTERIA §3.4): tło = fit exp(Czeb. 3) do PASM danych WP70 (okno ±45 MeV zaślepione),
wstrzyknięcie N_inj (= znany uzysk M0) × [sygnał fd(dmu_t, sc_t) + κ_MC (DMM + CM)]; fit każdym ustawieniem A0–A8. Okno sygnału danych nieużywane.
  fd_toys.py [--ntoy 200] [--nproc 7] → out/toys.json"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fd_lib as D  # noqa: E402
F, V = D.F, D.V
SC_T = (1.05, 0.10)


def task(t):
    bn, run, ntoy, chunk = t
    V._init()
    sigT, k = D.sig_tmpl(bn)
    T = D.tmpls(D.comps("full", bn))
    blind = np.abs(D.CTR - D.C["m0"]) < D.C["win"]
    Hd = V.wdata(k, V.SRC[run])[0]
    fb = F.fit(Hd[None], [sigT], [[]], D.FE, D.C["m0"], D.C["nbkg"], shape="fix", win=D.C["win"], blind=blind, hesse=False)
    bkg = fb["mu"][0]
    ninj = D.M0_DATA[run][D.BINS.index(bn)]
    s0 = sigT.bins(D.FE, D.C["m0"])
    info = {"N_inj": ninj, "bkg_win": float(bkg[blind].sum()), "SB_win": ninj * float(s0[blind].sum()) / float(bkg[blind].sum()), "chi2_sb": fb["chi2"], "ndf_sb": fb["ndf"], "sb_valid": fb["valid"]}
    rng = np.random.default_rng(20261120 + 1000 * D.BINS.index(bn) + (70000 if run == "428770" else 0) + 13 * chunk)
    rows = {m: [] for m in D.SET}
    nan = float("nan")
    for _ in range(ntoy):
        dm = rng.normal(0, D.GC[0])
        s = min(max(rng.normal(*SC_T), 0.80), 1.40)
        H = rng.poisson(bkg + ninj * (sigT.bins(D.FE, D.C["m0"], dm, s) + T["pk"][1] * T["pk"][0].bins(D.FE, D.C["m0"], dm, 1.0))).astype(float)
        for m in D.SET:
            r = D.res(D.fitA(m, H, sigT, T))
            rows[m].append([r["N"], r["err"] if r["err"] is not None else nan, float(r["valid"]), float(r["at_limit"]), r["sc"], r["dmu"], r["p"], dm, s])
    return {"bin": bn, "run": run, "chunk": chunk, "info": info, "rows": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ntoy", type=int, default=200); ap.add_argument("--nproc", type=int, default=7); ap.add_argument("--nchunk", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    tasks = [(bn, run, a.ntoy // a.nchunk, ch) for bn in D.BINS for run in V.SRC for ch in range(a.nchunk)]
    with Pool(a.nproc) as p:
        R = []
        for i, r in enumerate(p.imap_unordered(task, tasks)):
            R.append(r)
            print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(tasks)}", flush=True)
    M = {}
    for r in sorted(R, key=lambda r: (D.BINS.index(r["bin"]), r["run"], r["chunk"])):
        o = M.setdefault((r["bin"], r["run"]), {"bin": r["bin"], "run": r["run"], "info": r["info"], "rows": {m: [] for m in D.SET}})
        for m in D.SET:
            o["rows"][m] += r["rows"][m]
    json.dump({"res": list(M.values()), "cols": ["N", "err", "valid", "at_limit", "sc", "dmu", "p", "dmu_t", "sc_t"], "sc_t": SC_T}, open(os.path.join(D.BASE, "out", "toys.json"), "w"))
    print(f"zapisano out/toys.json, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
