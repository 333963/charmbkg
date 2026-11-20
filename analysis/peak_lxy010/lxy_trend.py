"""PEAK_lxy010: trend pozycji/szerokości piku w plasterkach L_xy POWYŻEJ 0,25 mm (jedyny zakres, gdzie lokalnie jest MC).
Dane 428770 pełna lumi (v1.2, WP70 D0_H) vs MC fd_nd; kształt MC dopasowany osobno w każdym plasterku.
Służy do oceny, czy ekstrapolacja do 0,10–0,25 mm jest bezpieczna (czy μ_MC, σ_MC i dane−MC zależą od L_xy).
Wynik: out/lxy_trend.csv"""
import sys, os
import numpy as np, pandas as pd
W = "/home/k/Documents/charm/charm_work"
sys.path.insert(0, f"{W}/PEAK_pos"); sys.argv = sys.argv[:1]
from fitlib import fit_mc, fit_data          # noqa: E402
from common import THR, M0, pass_wp          # noqa: E402
from dmfit import fit_dm, DM0               # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"
COLS = ["m", "pt", "dm", "score", "hypo_from_tag", "dstar_conflict", "Lxy", "truth", "hypo_swap"]


def sel(L):
    return L[(L.hypo_from_tag == 0) & (L.dstar_conflict == 0) & pass_wp(L.score.to_numpy(), L.pt.to_numpy(), THR["D0"], 0.7)]


S = pd.read_parquet(f"{W}/HYP_d0/runs/lxplus_dstar/fd_nd.parquet", columns=COLS)
S = sel(S[np.isin(S.truth, [0, 1]) & (S.hypo_swap == 0)])
D = sel(pd.read_parquet(f"{W}/HYP_d0/runs/lxplus_dstar/data_428770.parquet", columns=COLS[:7]))
SL = [(0.25, 0.40), (0.40, 0.70), (0.70, 1e9), (0.25, 1e9)]
rows = []
for lo, hi in [(3, 5), (5, 8), (8, 30), (3, 30)]:
    for a, b in SL:
        s = S[(S.pt >= lo) & (S.pt < hi) & (S.Lxy >= a) & (S.Lxy < b)]; d = D[(D.pt >= lo) & (D.pt < hi) & (D.Lxy >= a) & (D.Lxy < b)]
        shd = fit_dm(s.dm.to_numpy()[(np.abs(s.m - M0["D0"]) < 0.04).to_numpy()])
        shm = fit_mc(s.m.to_numpy()[(np.abs(s.dm - DM0) < 0.0012).to_numpy()], 1.76, 1.98)
        r1 = fit_dm(d.dm.to_numpy()[(np.abs(d.m - M0["D0"]) < 0.04).to_numpy()], shape=shd)
        r2 = fit_data(d.m.to_numpy()[(np.abs(d.dm - DM0) < 0.0012).to_numpy()], 1.76, 1.98, shm, bkg="cheb2", mu0=M0["D0"])
        for obs, r, sh, sc in (("dm", r1, shd, 1e6), ("mKpi", r2, shm, 1e3)):
            rows.append(dict(bin=f"pT{lo}-{hi}", lxy=f"{a:.2f}-{b:.2f}" if b < 1e8 else f">={a:.2f}", obs=obs, n_mc=len(s),
                             mu_mc=sh["mu"] * sc, e_mu_mc=sh["e_mu"] * sc, seff_mc=sh["seff"] * sc,
                             Ns=r["Ns"], Z=r["Ns"] / r["e_Ns"], dmu=(r["mu"] - sh["mu"]) * sc, e_dmu=np.hypot(r["e_mu"], sh["e_mu"]) * sc,
                             k=r["k"], e_k=r["e_k"], ok=r["ok"]))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/lxy_trend.csv", index=False)
pd.set_option("display.width", 250)
print(T.round(3).to_string())
