#!/usr/bin/env python3
"""M8_przeliczenie pkt 3b — G_fit = [N/σ_N]_sieć / [N/σ_N]_drzewo przy WP70 z fitów jednej kategorii (wps_hist.npz), M0 vs M8; D0. POST HOC.
Sieć: κ_MC z prawdy mb_nd przy WP70 sieci (out/R_m8.json, kap_cum.net) — ROZSZERZENIE (przepis M8 ma κ tylko dla drzewa); dane: × 1,364. Drzewo: M8 bez zmian.
Dodatkowo domknięcie N_fit/N_prawda na mb_nd dla sieci (κ_MC).   → out/gfit_m8.json, tabele/gfit_m8.md"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, "..", "..", "FITBIAS_2026-11-19", "skrypty"))
import fb_fit as F  # noqa: E402
import fb_var as V  # noqa: E402

NET, TREE, KSC = "D0.tr.H+e04.mc.d5.s1", "D0.xgb.H.sb.n.s1", 1.364
RJ = {(r["bin"], r["wp"]): r for r in json.load(open(os.path.join(HERE, "..", "out", "R_m8.json")))["res"] if r["tag"] == NET}


def one(D, k, m):
    f = V.fit_set(D[0][None], [k], "D0", m)
    r = V.sres(f)
    r["p"] = F.pval(f["chi2"], f["ndf"])
    b = []
    for i in range(1, D.shape[0]):
        fb = V.fit_set(D[i][None], [k], "D0", m, start=f["x"], hesse=False)
        if fb["valid"]:
            b.append(fb["N"][0])
    r["boot_rms"] = float(np.std(b)) if len(b) > 5 else None
    r["sig"] = max(r["err"] or 0.0, r["boot_rms"] or 0.0)
    return r


def task(t):
    bn, cls = t
    V.TAG["D0"] = NET if cls == "siec" else TREE
    k = V.wcat("D0", bn, 0.7)
    kmc = RJ[(bn, 0.7)]["kap_cum"]["net" if cls == "siec" else "tree"]
    res = {"bin": bn, "cls": cls, "eps": k["e"], "kap_T_interp": k["kap"], "kap_truth": kmc}
    if cls == "siec":
        k["kap"] = kmc                                               # drzewo: κ wg przepisu M8 (interpolacja T)
    kd = dict(k); kd["kap"] = k["kap"] * KSC
    for run, src in V.SRC.items():
        D = V.wdata(kd, src)
        res[run] = {m: one(D, kd, m) for m in ("M0", "M8")}
    H = V.wdata(k, "mb|all")
    nt = float(V.wdata(k, "mb|sig").sum())
    res["mc"] = {"N_true": nt}
    for m in ("M0", "M8"):
        f = V.fit_set(H[None], [k], "D0", m)
        res["mc"][m] = {"N": f["N"][0], "err": (f.get("err") or [None])[0], "b": f["N"][0] / nt - 1, "p": F.pval(f["chi2"], f["ndf"]), "valid": f["valid"], "at_limit": f["at_limit"]}
    return res


def main():
    with Pool(6, initializer=V._init) as p:
        out = p.map(task, [(bn, c) for bn in F.BINS for c in ("drzewo", "siec")])
    json.dump(out, open(os.path.join(HERE, "..", "out", "gfit_m8.json"), "w"), indent=1)
    R = {(o["bin"], o["cls"]): o for o in out}
    L = ["# M8_przeliczenie — pkt 3b: G_fit(D⁰, WP70) z fitów jednej kategorii, sieć `tr.H+e04.mc.d5.s1` vs drzewo H; post hoc", "",
         "σ = max(HESSE, bootstrap 30 replik). Sieć z M8: κ z prawdy mb_nd przy WP70 sieci (rozszerzenie przepisu), dane × 1,364.", "",
         "| próbka | pT | model | drzewo: N ± σ (p) | sieć: N ± σ (p) | G_fit = [N/σ]_sieć/[N/σ]_drzewo |", "|---|---|---|---|---|---|"]
    for run in V.SRC:
        for bn in F.BINS:
            for m in ("M0", "M8"):
                t, n = R[(bn, "drzewo")][run][m], R[(bn, "siec")][run][m]
                g = (n["N"] / n["sig"]) / (t["N"] / t["sig"])
                L.append(f"| {run} | {bn} | {m} | {t['N']:.0f} ± {t['sig']:.0f} ({t['p']:.2g}) | {n['N']:.0f} ± {n['sig']:.0f} ({n['p']:.2g}){'' if n['valid'] and not n['at_limit'] else ' [nieważny/granica]'} | {g:.2f} |")
    L += ["", "Domknięcie fitu jednej kategorii na mb_nd przy WP70 (b = N_fit/N_prawda − 1 ± HESSE; κ_MC):", "", "| pT | klasyfikator | κ_MC (prawda; dla drzewa obok κ z interpolacji T użyte w M8) | b M0 | b M8 |", "|---|---|---|---|---|"]
    for bn in F.BINS:
        for c in ("drzewo", "siec"):
            o = R[(bn, c)]
            L.append(f"| {bn} | {c} | {o['kap_truth']:.2f}" + (f" ({o['kap_T_interp']:.2f})" if c == "drzewo" else "") + " | " +
                     " | ".join(f"{o['mc'][m]['b']:+.3f} ± {(o['mc'][m]['err'] or float('nan')) / o['mc']['N_true']:.3f}" for m in ("M0", "M8")) + " |")
    open(os.path.join(HERE, "..", "tabele", "gfit_m8.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
