#!/usr/bin/env python3
"""FITBIAS diagnoza sygnału: kształt sygnału prawdy mb_nd wobec szablonu fd_nd tej samej kategorii, BEZ tła w modelu (N stałe = N_prawda w zakresie fitu):
χ² przy (dmu, sc) = (0, 1) oraz najlepsze (dmu, sc). Wyjście: tabele/diag_sig.md"""
import os, sys, json
import numpy as np
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import fb_fit as F, fb_diag_W as W
W._init()
L = ["# FITBIAS — kształt sygnału prawdy mb_nd vs szablon fd_nd (bez tła; normalizacja do liczby w zakresie fitu)", "",
     "| gat. | pT | kat. | N mb (zakres fitu) | χ²/ndf (0; 1) | p | dmu [MeV] | sc | χ²/ndf (dopasowane) | udział sygnału poza oknem ±win: mb / fd | udział poza zakresem fitu: mb / fd |", "|---|---|---|---|---|---|---|---|---|---|---|"]
out = []
for sp in ("D0", "Dplus"):
    c = F.SP[sp]; lo, hi = c["mr"]; edges = np.arange(c["fit"][0], c["fit"][1] + F.BW / 2, F.BW); ctr = 0.5 * (edges[1:] + edges[:-1]); wm = np.abs(ctr - c["m0"]) < c["win"]
    for bn in F.BINS:
        for n in ("WP90", "WP70", "A50", "B50", "L90-70"):
            b = W.build(sp, bn, n)
            h = F.rebin1(b["comp"]["Sig"], lo, *c["fit"]); N = h.sum()
            def chi(p):
                s = b["sigT"].bins(edges, c["m0"], p[0] / 1000, p[1]); mu = np.maximum(N * s / s.sum(), 1e-6)
                return 2 * np.sum(mu - h + np.where(h > 0, h * np.log(np.where(h > 0, h, 1) / mu), 0))
            c0 = chi([0, 1]); r = minimize(chi, [0, 1], method="Nelder-Mead")
            t0 = b["sigT"].bins(edges, c["m0"]); 
            fo_mb, fo_fd = 1 - h[wm].sum() / N, 1 - t0[wm].sum() / t0.sum()
            L.append(f"| {sp} | {bn} | {n} | {N:.0f} | {c0:.0f}/{len(h)-1} | {F.pval(c0, len(h)-1):.2g} | {r.x[0]:+.2f} | {r.x[1]:.3f} | {r.fun:.0f}/{len(h)-3} | {fo_mb:.3f} / {fo_fd:.3f} | {1 - N / b['N_true']:.3f} / {1 - t0.sum():.3f} |")
open(os.path.join(HERE, "..", "tabele", "diag_sig.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
