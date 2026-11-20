"""FEAT_charmnet: modelowanie cech dane vs MC w dalekich sidebandach (te same weta co tło treningowe testu data).
MC = mb_nd, wszystkie kandydaty (dowolny truth — SB danych też zawiera wszystko), przeważone do (pT,|η|) danych.
Per cecha: KS (ważony), różnica median / IQR, C2ST 1D (AUC HistGB OOF); wielowymiarowo C2ST dla H vs H+grupa.
Użycie: python modelling.py SPECIES
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from study import BASE, GROUPS, H, cfg, load  # noqa: E402
from charmbkg.samples import reweight  # noqa: E402


def wks(x, wx, y, wy):
    v = np.r_[x, y]
    o = np.argsort(v, kind="stable")
    a = np.r_[wx / wx.sum(), np.zeros(len(y))][o]
    b = np.r_[np.zeros(len(x)), wy / wy.sum()][o]
    return float(np.max(np.abs(np.cumsum(a) - np.cumsum(b))))


def c2st(Xd, Xm, wm, seed=0):
    X = np.vstack([Xd, Xm])
    y = np.r_[np.ones(len(Xd)), np.zeros(len(Xm))]
    w = np.r_[np.ones(len(Xd)) / len(Xd), wm / wm.sum()] * len(y) / 2
    f = np.random.default_rng(seed).integers(0, 3, len(y))
    s = np.zeros(len(y))
    for k in range(3):
        c = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.1, max_leaf_nodes=31, random_state=seed)
        c.fit(X[f != k], y[f != k], sample_weight=w[f != k])
        s[f == k] = c.predict_proba(X[f == k])[:, 1]
    return float(roc_auc_score(y, s, sample_weight=w))


def main(sp):
    m0 = cfg["species"][sp]["m0"]
    S, Bmc, Bda, dall, mb = load(sp)
    far = {"D0": 0.100, "Dplus": 0.110}[sp]
    sel = mb["ok"].to_numpy() & (mb["dm0"].to_numpy() > far)
    if sp == "D0":
        sel &= np.abs(mb["m_swap"].to_numpy() - m0) > 0.025
    M = mb[sel].reset_index(drop=True)
    rng = np.random.default_rng(3)
    n = 300000
    D = Bda.iloc[np.sort(rng.choice(len(Bda), min(n, len(Bda)), replace=False))].reset_index(drop=True)
    M = M.iloc[np.sort(rng.choice(len(M), min(n, len(M)), replace=False))].reset_index(drop=True)
    D["weight"] = 1.0
    wm = reweight(M, D, cfg["training"]["reweight_bins"])
    allnew = list(dict.fromkeys(f for g in GROUPS[sp].values() for f in g))
    res = {"n_data": len(D), "n_mc": len(M), "features": {}, "c2st": {}, "per_bin_iso": {}}
    for f in H + allnew:
        x, y = D[f].to_numpy(float), M[f].to_numpy(float)
        okx, oky = np.isfinite(x), np.isfinite(y)
        q = lambda v, w: np.interp([0.25, 0.5, 0.75], np.cumsum(w[np.argsort(v)]) / w.sum(), np.sort(v))  # noqa: E731
        qd, qm = q(x[okx], np.ones(okx.sum())), q(y[oky], wm[oky])
        res["features"][f] = dict(new=f not in H, ks=wks(x[okx], np.ones(okx.sum()), y[oky], wm[oky]),
                                  med_shift_iqr=float((qd[1] - qm[1]) / max(qd[2] - qd[0], 1e-9)),
                                  med_data=float(qd[1]), med_mc=float(qm[1]),
                                  c2st1d=c2st(np.nan_to_num(x[:, None], nan=-9), np.nan_to_num(y[:, None], nan=-9), wm))
        print(sp, f, res["features"][f], flush=True)
    for gname, feats in [("H", H)] + [("H+" + k, H + v) for k, v in GROUPS[sp].items()]:
        res["c2st"][gname] = c2st(np.nan_to_num(D[feats].to_numpy(float), nan=-9), np.nan_to_num(M[feats].to_numpy(float), nan=-9), wm)
        print(sp, gname, res["c2st"][gname], flush=True)
    for lo, hi in ((3, 5), (5, 8), (8, 30)):
        a, b = (D.pt >= lo) & (D.pt < hi), (M.pt >= lo) & (M.pt < hi)
        res["per_bin_iso"][f"{lo}-{hi}"] = {f: dict(ks=wks(D.loc[a, f].to_numpy(float), np.ones(a.sum()), M.loc[b, f].to_numpy(float), wm[b.to_numpy()]),
                                                    mean_data=float(D.loc[a, f].mean()), mean_mc=float(np.average(M.loc[b, f], weights=wm[b.to_numpy()])))
                                            for f in ("ptcone40_cn", "ptcone100_cn")}
    # izolacja sygnału: fd (sygnał) vs tło — informacyjnie, średnie
    res["signal_iso_mean"] = {f: float(S[f].mean()) for f in ("ptcone40_cn", "ptcone100_cn")}
    os.makedirs(os.path.join(BASE, "out"), exist_ok=True)
    json.dump(res, open(os.path.join(BASE, "out", f"modelling_{sp}.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
