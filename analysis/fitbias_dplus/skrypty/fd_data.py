#!/usr/bin/env python3
"""FITBIAS_Dplus — DANE (po wyborze; KRYTERIA §5): uzyski D+ WP70 w 427929 i 428770 wszystkimi ustawieniami (κ × 1), bootstrap 30 replik po zdarzeniach,
p(χ²), rozkład χ² po obszarach masy, reszty 25 MeV; diagnostyka 3–5 GeV (opisowa): sc wolne vs stałe, dodatkowy Gauss w 2,010 GeV (skan normalizacji).
  fd_data.py [--nproc 6] → out/data.json"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fd_lib as D  # noqa: E402
F, V = D.F, D.V
REG = [("1,70–1,78", 1.70, 1.78), ("1,78–1,825", 1.78, 1.825), ("okno 1,825–1,915", 1.825, 1.915), ("1,915–2,00", 1.915, 2.00), ("2,00–2,05", 2.00, 2.05), ("2,05–2,20", 2.05, 2.20)]


def gauss_tmpl(mu, sig):
    from scipy.stats import norm
    return F.TmplE(D.FE, np.diff(norm.cdf(D.FE, mu, sig)) + 1e-12)


def task(t):
    bn, run = t
    V._init()
    sigT, k = D.sig_tmpl(bn)
    T = D.tmpls(D.comps("full", bn))
    Dt = V.wdata(k, V.SRC[run])                                   # [31, 100]
    H = Dt[0]
    out = {"bin": bn, "run": run, "n_cand": float(H.sum()), "kap": {x: T[x][1] for x in T}, "models": {}, "h25": H.reshape(-1, 5).sum(1).tolist()}
    for m in list(D.SET):
        S = D.SET[m]
        f = D.fitA(m, H, sigT, T, keep_mu=True)
        r = D.res(f)
        mu = f["mu"][0]
        h = D.cutr(H, S["rng"])
        e = np.arange(S["rng"][0], S["rng"][1] + F.BW / 2, F.BW)
        c = 0.5 * (e[1:] + e[:-1])
        r["reg"] = [{"name": nm, "chi2": float(np.sum(((h - mu) ** 2 / mu)[(c > lo) & (c < hi)])), "nb": int(((c > lo) & (c < hi)).sum()), "sum_res": float(np.sum((h - mu)[(c > lo) & (c < hi)])),
                     "sum_mu": float(np.sum(mu[(c > lo) & (c < hi)]))} for nm, lo, hi in REG]
        nb = (len(h) // 5) * 5
        r["pull25"] = ((h[:nb].reshape(-1, 5).sum(1) - mu[:nb].reshape(-1, 5).sum(1)) / np.sqrt(mu[:nb].reshape(-1, 5).sum(1))).tolist()
        r["lo25"] = float(e[0])
        bN = []
        for i in range(1, Dt.shape[0]):
            fb = D.fitA(m, Dt[i], sigT, T, start=f["x"], hesse=False)
            bN.append(fb["N"][0] if fb["valid"] else None)
        r["boot"] = bN
        g = [x for x in bN if x is not None]
        r["boot_rms"] = float(np.std(g, ddof=1)) if len(g) >= 20 else None
        # diagnostyka opisowa: kształt stały (sc = 1, dmu = 0)
        ff = D.fitA(m, H, sigT, T, shape="fix")
        r["fix"] = {"N": ff["N"][0], "chi2": ff["chi2"], "ndf": ff["ndf"], "p": F.pval(ff["chi2"], ff["ndf"]), "fval": ff["fval"]}
        r["fval"] = f["fval"]
        # diagnostyka opisowa: dodatkowy Gauss w 2,010 GeV (σ = 8 MeV), normalizacja = κ_g × N_S, skan κ_g
        if S["rng"][1] >= 2.05:
            G = gauss_tmpl(2.0103, 0.008)
            sc_ = []
            for kg in np.linspace(-0.06, 0.12, 19):
                T2 = dict(T); T2["pk2"] = (G, float(kg))
                S2 = dict(S); S2["tied"] = tuple(S["tied"]) + ("pk2",)
                D.SET["_tmp"] = S2
                f2 = D.fitA("_tmp", H, sigT, T2, start=f["x"], hesse=False)
                sc_.append([float(kg), f2["fval"], f2["N"][0], f2["chi2"]])
            D.SET.pop("_tmp", None)
            r["dstar_scan"] = sc_
        out["models"][m] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=6)
    a = ap.parse_args()
    with Pool(a.nproc) as p:
        R = p.map(task, [(bn, run) for bn in D.BINS for run in V.SRC])
    json.dump(R, open(os.path.join(D.BASE, "out", "data.json"), "w"))
    print("zapisano out/data.json")


if __name__ == "__main__":
    main()
