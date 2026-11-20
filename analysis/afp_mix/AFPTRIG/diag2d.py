#!/usr/bin/env python3
"""AFPTRIG/diag2d.py — rozróżnienie „przeciek sygnału” vs „niedopasowanie puli”: stosunek dane/szablon w (ξ_p, ξ_cal).
Szablon: pula centralna × widmo protonów niemożliwych, normalizacja per komórka (blok LB × nSelBL) w ξ_cal > 0,20 ∧ Δ > 0,045.
Protony przypadkowe mają ξ_p niezależne od ξ_cal, więc niedopasowanie puli daje ten sam stosunek r(ξ_cal) dla każdego ξ_p,
a przeciek własnych protonów daje nadwyżkę wzdłuż ξ_cal ≈ ξ_p + stała. → diag2d.json"""
import json, sys, numpy as np, lib
from analysis import Run, NCLS, nominal_pool, POOL_THR
XP = np.array([0.03, 0.045, 0.06, 0.08, 0.10001]); XC = np.array([0.05, 0.065, 0.08, 0.10, 0.12, 0.14, 0.17, 0.20, 0.30, 10.0])
OUT = {}
for run in lib.RUNS:
    R = Run(run); d = R.d; O = OUT[str(run)] = {}
    for stream in ("pt4_L1AFP", "noalg_L1AFP", "afp_any", "nonAFP"):
        pool = nominal_pool(run, stream); D2 = np.zeros((len(XP) - 1, len(XC) - 1)); P2 = np.zeros_like(D2); V2 = np.zeros_like(D2)
        for s in lib.SIDES:
            g = R.groups(s, stream, 8, "lb"); cell = g * NCLS + R.cls; nc = int(cell.max()) + 1
            dm = R.S[stream] & R.has & d[f"tag_{s}"]; pm = R.S[pool] & R.has & (d[f"np_{s}"] == 0); qm = R.protons(s, stream, "impossible")
            xcd, xpd = d[f"xical_{s}"][dm], d[f"xi1_{s}"][dm]
            bc = np.digitize(xcd, XC) - 1; bp = np.digitize(xpd, XP) - 1; ok = (bc >= 0) & (bc < len(XC) - 1) & (bp >= 0) & (bp < len(XP) - 1)
            np.add.at(D2, (bp[ok], bc[ok]), 1)
            nb = len(XC) - 1
            Hd = np.zeros((nc, nb)); np.add.at(Hd, (cell[dm][ok], bc[ok]), 1)       # dane per komórka × bin ξ_cal (wszystkie ξ_p)
            bcp = np.digitize(d[f"xical_{s}"][pm], XC) - 1; okp = bcp >= 0
            Hp = np.zeros((nc, nb)); np.add.at(Hp, (cell[pm][okp], bcp[okp]), 1)
            far = np.arange(nb) >= 7                                                 # ξ_cal > 0,20 → Δ > 0,10 dla każdego ξ_p
            k = np.where(Hp[:, far].sum(1) > 0, Hd[:, far].sum(1) / np.maximum(Hp[:, far].sum(1), 1), 0.0)
            # widmo protonów per blok LB
            ng = int(g.max()) + 1
            hp = np.zeros((ng, len(XP) - 1)); bq = np.digitize(d[f"xi1_{s}"][qm], XP) - 1; np.add.at(hp, (g[qm], np.clip(bq, 0, len(XP) - 2)), 1)
            hp = hp / np.maximum(hp.sum(1, keepdims=True), 1)
            cg = np.arange(nc) // NCLS
            for c in range(nc):
                if k[c] > 0:
                    P2 += k[c] * np.outer(hp[cg[c]], Hp[c]); V2 += (k[c] ** 2 / max(Hd[c, far].sum(), 1)) * np.outer(hp[cg[c]], Hp[c]) ** 2
        # maska fizyczna: bin w całości w Δ > 0,02
        r = np.where(P2 > 0, D2 / np.where(P2 > 0, P2, 1), np.nan); e = r * np.sqrt(1 / np.maximum(D2, 1) + V2 / np.maximum(P2, 1e-9) ** 2)
        O[stream] = {"xp_edges": XP.tolist(), "xc_edges": XC.tolist(), "data": D2.tolist(), "pred": P2.tolist(), "ratio": r.tolist(), "ratio_err": e.tolist()}
        print("==", run, stream, "pool", pool)
        print("   ξ_cal bins:", " ".join(f"{XC[i]:.3f}-{XC[i+1]:.2f}" for i in range(len(XC) - 1)))
        for i in range(len(XP) - 1):
            print(f"   ξ_p {XP[i]:.3f}-{XP[i+1]:.3f}:", " ".join((f"{r[i, j]:5.2f}±{e[i, j]:.2f}" if XC[j] - XP[i + 1] >= 0.02 - 1e-9 else "   (SR/gap)") for j in range(len(XC) - 1)))
        R.cache.clear()
json.dump(OUT, open("diag2d.json", "w"), indent=1)
