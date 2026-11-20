#!/usr/bin/env python3
"""KONSTR A — zamknięcie pomiaru odbicia na mb_nd z prawdą (KRYTERIA §A a–c). Wyjście out/ka_mc.json, tabele/A_mc.md"""
import json, os, sys
import numpy as np
from scipy.stats import chi2 as C2
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import ka_lib as K
T = os.path.dirname(H)
MB, FD = np.load(f"{T}/in/kx_mb_nd.npz"), np.load(f"{T}/in/kx_fd_nd.npz")
rng = np.random.default_rng(20261120)
out, L = {}, ["| bin | WP | test | A_R | prawda (t1) | A/prawda | χ²/ndf (p) wobec prawdy per przedział | s_R (wolne) | status |", "|---|---|---|---|---|---|---|---|---|"]
for bn in K.BINS:
    for e in (0.7, 0.5, 0.9):
        g = {c: K.cum(MB[f"G|mb_nd|{bn}|{c}"], e) for c in ("t0", "t1", "t2", "t3")}
        fd1 = K.cum(FD[f"G|fd_nd|{bn}|t1"], e)
        tr = g["t1"].sum()
        trs = g["t1"].reshape(14, 5, 110).sum((1, 2))
        thin = lambda h: rng.binomial(h.astype(int), 1 / 1.364).astype(float)
        tests = {"a: wszystko": (sum(g.values()), tr, trs), "b: samo tło (t2+t3)": (g["t2"] + g["t3"], 0.0, 0 * trs), "b': tło + sygnał poprawny": (g["t0"] + g["t2"] + g["t3"], 0.0, 0 * trs)}
        t0, t1 = thin(g["t0"]), thin(g["t1"])
        tests["c: czystość danych (S×1/1,364)"] = (t0 + t1 + g["t2"] + g["t3"], t1.sum(), t1.reshape(14, 5, 110).sum((1, 2)))
        for tn, (G2, tru, trus) in tests.items():
            r = K.measure(G2, fd1)
            rf = K.measure(G2, fd1, free_s=True)
            u = K.USE
            c2 = float(np.sum(((r["y"][u] - trus[u]) / r["ey"][u]) ** 2)); p = float(C2.sf(c2, int(u.sum())))
            if tru > 0:
                rat, z = r["A"] / tru, (r["A"] - tru) / r["eA"]
                ok = p >= 0.01 and abs(z) <= 2 and abs(rat - 1) <= 0.10
                rs = f"{rat:.3f} ± {r['eA'] / tru:.3f}"
            else:
                z = r["A"] / r["eA"]; ok = abs(z) <= 2; rs = f"Z = {z:+.1f}"
            out[f"{bn}|{e}|{tn}"] = {"A": r["A"], "eA": r["eA"], "truth": float(tru), "chi2": c2, "ndf": int(u.sum()), "p": p, "s_free": rf["s"], "A_free": rf["A"], "ok": bool(ok), "eff": r["eff"],
                                     "y": r["y"].tolist(), "ey": r["ey"].tolist(), "truth_sl": np.asarray(trus, float).tolist()}
            L.append(f"| {bn} | {int(e * 100)} | {tn} | {r['A']:.0f} ± {r['eA']:.0f} | {tru:.0f} | {rs} | {c2:.1f}/{int(u.sum())} ({p:.2g}) | {rf['s']:.2f} (A {rf['A']:.0f}) | {'OK' if ok else 'NIE'} |")
json.dump(out, open(f"{T}/out/ka_mc.json", "w"), indent=1)
open(f"{T}/tabele/A_mc.md", "w").write("# A — zamknięcie pomiaru odbicia z bliźniaków na mb_nd (prawda)\n\n" + "\n".join(L) + "\n")
print("\n".join(L))
