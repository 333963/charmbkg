#!/usr/bin/env python3
"""FITBIAS diagnoza T2 — udział „sygnałopodobny” składowych tła: N_spur/N_składowej na fd_nd (duża statystyka) i mb_nd; κ·frakcja = obciążenie.
Szablon sygnału stały, tło exp nominalne. Wyjście: out/diag_T2.json + tabele/diag_T2.md"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import fb_fit as F
import fb_diag_T as D
D._init()
NAMES = ["SigC", "SigB", "ReflSwap", "ReflOther", "CM_411", "CM_421", "CM_431", "CM_413", "CM_423", "CM_4122", "CM_other", "DMM_C", "DMM_B",
         "PC_411", "PC_421", "PC_423", "PC_413", "PC_431", "PC_4122", "PC_other", "MatchedBL", "OtherNoB", "OtherB", "Unknown"]
out = []
for sp in ("D0", "Dplus"):
    c = F.SP[sp]; lo, hi = c["mr"]
    fe = np.arange(lo, hi + F.BW / 2, F.BW); edges = np.arange(c["fit"][0], c["fit"][1] + F.BW / 2, F.BW)
    for ib in range(3):
        for cat in ("WP80", "L80-50", "WP50"):
            mb, fd = D.comp(D.M, sp, ib, D.CATS[cat]), D.comp(D.FD, sp, ib, D.CATS[cat])
            sT = F.TmplE(fe, fd["Sig"]); nt = mb["Sig"].sum()
            for k in range(4, 20):
                row = {"sp": sp, "bin": F.BINS[ib], "cat": cat, "comp": NAMES[k], "kappa_mb": float(mb["fine"][k].sum() / nt), "kappa_fd": float(fd["fine"][k].sum() / fd["Sig"].sum())}
                for nm, src in (("fd", fd), ("mb", mb)):
                    h = D.cut(src["fine"][k], sp)
                    if h.sum() < 200:
                        continue
                    f = F.fit(h[None], [sT], [[]], edges, c["m0"], c["nbkg"], shape="fix", win=c["win"])
                    row[nm] = {"n": float(h.sum()), "frac": f["N"][0] / src["fine"][k].sum(), "frac_err": (f["err"][0] / src["fine"][k].sum()) if "err" in f else None, "chi2": f["chi2"], "ndf": f["ndf"]}
                if "fd" in row or "mb" in row:
                    out.append(row)
json.dump(out, open(os.path.join(HERE, "..", "out", "diag_T2.json"), "w"))
f = lambda x, n=3: "—" if x is None else f"{x:+.{n}f}".replace(".", ",")
L = ["# FITBIAS — diagnoza T2: frakcja sygnałopodobna składowych tła (N_spur/N_składowej; szablon sygnału stały, exp nominalne) na fd_nd i mb_nd; obciążenie = κ_mb × frakcja", "",
     "| gat. | pT | kat. | składowa | κ_mb | κ_fd | frakcja fd | frakcja mb | obciążenie b = κ_mb·frakcja_fd | b z mb |", "|---|---|---|---|---|---|---|---|---|---|"]
for r in out:
    a, b = r.get("fd"), r.get("mb")
    if r["kappa_mb"] < 0.02: continue
    L.append(f"| {r['sp']} | {r['bin']} | {r['cat']} | {r['comp']} | {r['kappa_mb']:.3f} | {r['kappa_fd']:.3f} | " + (f"{f(a['frac'])} ± {f(a['frac_err'])[1:]}" if a else "—") + " | " + (f"{f(b['frac'])} ± {f(b['frac_err'])[1:]}" if b and b['frac_err'] else "—")
             + " | " + (f"{f(r['kappa_mb'] * a['frac'])} ± {f(r['kappa_mb'] * a['frac_err'])[1:]}" if a and a['frac_err'] else "—") + " | " + (f"{f(r['kappa_mb'] * b['frac'])} ± {f(r['kappa_mb'] * b['frac_err'])[1:]}" if b and b['frac_err'] else "—") + " |")
open(os.path.join(HERE, "..", "tabele", "diag_T2.md"), "w").write("\n".join(L) + "\n")
