"""PEAK_pos: kompaktowe ramki (dane dall + MC fd sygnał S) z wynikiem XGB H (model test=mc z FEAT_charmnet, OOF/fold-model).
Kolejność wierszy odtwarza FEAT_charmnet/scripts/study.load (bez jej importu, żeby nie ładować mb).
Użycie: python prep.py D0|Dplus
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/home/k/Projects/charm/charm-bkg")
from charmbkg import config as C  # noqa
FC = "/home/k/Documents/charm/charm_work/FEAT_charmnet"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
sp = sys.argv[1]
cfg = C.load(); m0 = cfg["species"][sp]["m0"]
KEEP = ["m", "pt", "eta", "phi", "charge", "truth", "evkey", "weight", "m_swap", "hypo_swap", "dm",
        "d_pt0", "d_eta0", "d_phi0", "d_pt1", "d_eta1", "d_phi1", "d_pt2", "d_eta2", "d_phi2", "d_trk0"]
z = np.load(os.path.join(FC, "out", f"{sp}_mc", "scores_H.npz"))
zd = np.load(os.path.join(FC, "out", f"{sp}_data", "scores_H.npz")) if os.path.exists(os.path.join(FC, "out", f"{sp}_data", "scores_H.npz")) else None
for s in ["fd_nd", "data_428770"]:
    L = pd.read_pickle(os.path.join(FC, "frames", f"{s}_{sp}.pkl"))
    ok = np.ones(len(L), bool)
    if sp == "Dplus":
        ok &= L["dstar_ok"].to_numpy(bool)
    if s == "fd_nd":
        t = L["truth"].to_numpy(float)
        sel = np.isin(t, [0, 1]) & ok
        if sp == "D0":
            sel &= L["hypo_swap"].to_numpy() != 1
        key, name = "S", "mc"
    else:
        sel = ok.copy()
        if sp == "D0":
            sel &= np.abs(L["m_swap"].to_numpy() - m0) > 0.025
        key, name = "dall", "data"
    cols = [c for c in KEEP if c in L.columns]
    A = L.loc[sel, cols].reset_index(drop=True)
    assert len(A) == len(z[key]), (len(A), len(z[key]))
    A["score"] = z[key]
    if zd is not None:
        A["score_dtrain"] = zd[key]
    if s == "fd_nd" and sp == "D0":  # odbicia (hypo_swap==1) bez wyniku — do oceny resztkowej refleksji po wecie
        R = L.loc[np.isin(L["truth"].to_numpy(float), [0, 1]) & (L["hypo_swap"].to_numpy() == 1), cols].reset_index(drop=True)
        R.to_parquet(os.path.join(OUT, f"{sp}_mc_refl.parquet"))
    A.to_parquet(os.path.join(OUT, f"{sp}_{name}.parquet"))
    print(sp, s, len(L), "->", len(A), cols)
    del L
