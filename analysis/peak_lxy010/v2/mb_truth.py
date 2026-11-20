"""PEAK_lxy010 v2 — wariant MC: sygnał z prawdy mb_nd, m(Kπ) dopasowane swobodnie (fit_mc: 2 Gaussy bez tła) i porównane z fd_nd przy tym samym
cięciu L_xy. (W fit_v2.py sygnał mb był dopasowany jak dane, z tłem Czebyszewa — dla czystego sygnału m(Kπ) tło jest zdegenerowane z szerokim Gaussem,
dlatego dla m(Kπ) używana jest ta tabela; dla Δm zostaje rows_mb_nd_truth.csv.)  Wymóg n ≥ 300.  → out/mb_truth_free.csv"""
import sys, os
import numpy as np, pandas as pd
sys.argv = [sys.argv[0], "mb_nd_skip"]
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"
src = open(f"{HERE}/fit_v2.py").read().split("t0 = time.time()")[0]      # tylko definicje (load, cut, binsel, xm, ...)
exec(compile(src, "fit_v2_defs", "exec"))
MC = pd.read_csv(f"{OUT}/mc_fd.csv").set_index(["sel", "bin", "cut", "obs"])
S = load("mb_nd", truth=True); rows = []
for sel in SELS:
    for lo, hi in PTB:
        for c in CUTS:
            x = xm(cut(binsel(S, sel, lo, hi), c)); f = MC.loc[(sel, f"pT{lo}-{hi}", c, "mKpi")]
            if len(x) < 300:
                rows.append(dict(sel=sel, bin=f"pT{lo}-{hi}", cut=c, n=len(x), ok=False)); continue
            r = fit_mc(x, 1.76, 1.98)
            rows.append(dict(sel=sel, bin=f"pT{lo}-{hi}", cut=c, n=len(x), ok=r["ok"], dmu=(r["mu"] - f.mu) * 1e3, e_dmu=np.hypot(r["e_mu"], f.e_mu) * 1e3,
                             rs=r["seff"] / f.seff, mean_diff=(x.mean() - 1.86484) * 1e3, rms=x.std() * 1e3))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/mb_truth_free.csv", index=False); print(T.round(3).to_string())
