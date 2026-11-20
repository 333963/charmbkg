#!/usr/bin/env python3
"""FITBIAS_Dplus — składowe objęte szablonem (DMM, CM_Ds, CM_oth) poza próbką, ustawienie A2: rzeczywiste zdarzenia jednej połowy mb_nd wobec szablonu z drugiej.
Pseudo-dane = N_prawda × [sygnał fd + κ szablon] − przewidywany udział składowej + rzeczywiste zdarzenia składowej + WYGŁADZONE tło (PartComb + Comb) tej połowy
(cokół jest konieczny: bez tła N_B stoi na granicy 0 i fit ma sztuczne obciążenie ujemne −2…−5 % — pierwsza wersja w mc.json['dec']['AtoB'|'BtoA'] dla składowych modelowanych jest przez to NIEWAŻNA).
→ out/mc_oos.json"""
import json
import os
import sys
import zlib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fd_lib as D  # noqa: E402
F, V = D.F, D.V
V._init()
NREP = 40
out = []
for bn in D.BINS:
    sigT, _ = D.sig_tmpl(bn)
    cm = {w: D.comps(w, bn) for w in ("A", "B")}
    T = {w: D.tmpls(cm[w]) for w in cm}
    for m in ("A2",):
        for cn in ("DMM", "CM_Ds", "CM_oth", "pk"):
            r = {"bin": bn, "model": m, "comp": cn}
            for tw, dw in (("A", "B"), ("B", "A")):
                c = cm[dw]
                nt = float(c["Sig"].sum())
                ped = D.smooth(c["PC"] + c["Comb"], sigT, nb=3)                      # ta sama rodzina co model → kontrola = 0
                base = nt * sigT.bins(D.FE, D.C["m0"]) + nt * T[tw]["pk"][1] * T[tw]["pk"][0].bins(D.FE, D.C["m0"]) + ped - nt * T[tw]["pk"][1] * cm[tw][cn] / cm[tw]["pk"].sum()
                f0 = D.fitA(m, base + nt * T[tw]["pk"][1] * cm[tw][cn] / cm[tw]["pk"].sum(), sigT, T[tw], shape="fix", hesse=False)      # kontrola: Asimov → 0
                f = D.fitA(m, base + c[cn], sigT, T[tw], shape="fix", hesse=False)
                rng = np.random.default_rng(zlib.crc32(f"{bn}|{cn}|{tw}".encode()))
                bN = [D.fitA(m, base + rng.poisson(c[cn]), sigT, T[tw], shape="fix", start=f["x"], hesse=False)["N"][0] for _ in range(NREP)]
                r[f"{tw}to{dw}"] = {"spur": (f["N"][0] - nt) / nt, "err": float(np.std(bN)) / nt, "ctrl": (f0["N"][0] - nt) / nt, "n": float(c[cn].sum()), "n_pred": float(nt * T[tw]["pk"][1] * cm[tw][cn].sum() / cm[tw]["pk"].sum())}
            a, b = r["AtoB"], r["BtoA"]
            r["mean"], r["mean_err"] = 0.5 * (a["spur"] + b["spur"]), 0.5 * float(np.hypot(a["err"], b["err"]))
            out.append(r)
            print(bn, cn, f"A→B {100 * a['spur']:+.1f} ± {100 * a['err']:.1f} (kontrola {100 * a['ctrl']:+.2f}; n {a['n']:.0f} vs {a['n_pred']:.0f}) | B→A {100 * b['spur']:+.1f} ± {100 * b['err']:.1f} | średnia {100 * r['mean']:+.1f} ± {100 * r['mean_err']:.1f}")
json.dump(out, open(os.path.join(D.BASE, "out", "mc_oos.json"), "w"), indent=1)
