#!/usr/bin/env python3
"""KONSTR A — zamknięcie estymatora E2 na mb_nd (testy a, b, c; 50 realizacji przerzedzenia dla c). → out/ka_mc2.json, tabele/A_mc2.md"""
import json, os, sys
import numpy as np
from scipy.stats import chi2 as C2
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import ka_lib as K
T = os.path.dirname(H)
MB, FD = np.load(f"{T}/in/kx_mb_nd.npz"), np.load(f"{T}/in/kx_fd_nd.npz")
rng = np.random.default_rng(11)
out, L = {}, ["| bin | WP | a: A_R / prawda (χ² wobec prawdy, p) | b: samo tło, A_R (Z; χ², p) | c: czystość danych, ⟨A/prawda⟩ ± rms (50) | χ² kształtu fd (p) | status |", "|---|---|---|---|---|---|---|"]
for bn in K.BINS:
    for e in (0.7, 0.5, 0.9):
        g = {c: K.cum(MB[f"G|mb_nd|{bn}|{c}"], e) for c in ("t0", "t1", "t2", "t3")}
        fd1 = K.cum(FD[f"G|fd_nd|{bn}|t1"], e)
        u = K.USE
        trs = g["t1"].reshape(14, 5, 110).sum((1, 2)); tr = trs.sum()
        a = K.measure2(sum(g.values()), fd1)
        c2a = float(np.sum(((a["y"][u] - trs[u]) / a["ey"][u]) ** 2)); pa = float(C2.sf(c2a, 9))
        b = K.measure2(g["t2"] + g["t3"], fd1)
        c2b = float(np.sum((b["y"][u] / b["ey"][u]) ** 2)); pb = float(C2.sf(c2b, 9))
        rc = []
        for _ in range(50):
            t0, t1 = rng.binomial(g["t0"].astype(int), 1 / 1.364), rng.binomial(g["t1"].astype(int), 1 / 1.364)
            c = K.measure2(t0 + t1 + g["t2"] + g["t3"], fd1)
            rc.append([c["A"] / t1.sum(), (c["A"] - t1.sum()) / c["eA"]])
        rc = np.array(rc)
        za = (a["A"] - tr) / a["eA"]
        ok = pa >= 0.01 and abs(za) <= 2 and abs(a["A"] / tr - 1) <= 0.10 and abs(b["A"] / b["eA"]) <= 2 and abs(rc[:, 0].mean() - 1) <= 0.10
        out[f"{bn}|{e}"] = {"A": a["A"], "eA": a["eA"], "truth": float(tr), "chi2_a": c2a, "p_a": pa, "A_b": b["A"], "eA_b": b["eA"], "chi2_b": c2b, "p_b": pb, "c_mean": float(rc[:, 0].mean()), "c_rms": float(rc[:, 0].std()),
                             "c_pull_rms": float(rc[:, 1].std()), "chi2_shape": a["chi2"], "ok": bool(ok)}
        L.append(f"| {bn} | {int(e * 100)} | {a['A']:.0f} ± {a['eA']:.0f} / {tr:.0f} = {a['A'] / tr:.3f} ± {a['eA'] / tr:.3f} ({c2a:.1f}/9, {pa:.2g}) | {b['A']:.0f} ± {b['eA']:.0f} ({b['A'] / b['eA']:+.1f}; {c2b:.1f}/9, {pb:.2g}) | "
                 f"{rc[:, 0].mean():.3f} ± {rc[:, 0].std():.3f} (rms pull {rc[:, 1].std():.2f}) | {a['chi2']:.1f}/8 ({C2.sf(a['chi2'], 8):.2g}) | {'OK' if ok else 'NIE'} |")
json.dump(out, open(f"{T}/out/ka_mc2.json", "w"), indent=1)
open(f"{T}/tabele/A_mc2.md", "w").write("# A — zamknięcie estymatora E2 (fit per przedział: sygnał w m_swap + pol2) na mb_nd\n\n" + "\n".join(L) + "\n")
print("\n".join(L))
