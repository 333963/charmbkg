"""DLZ_d0: testy rzeźbienia masy tła D0 po cięciu WP (modele testu mc), ponad R_sculpt z FEAT/sculpt.py:
 (1) fałszywy sygnał: fit masy tła MC truth≥8 (bez sygnału z konstrukcji) po cięciu WP70/WP50 modelem Gauss(σ z sygnału MC) + exp(pol2)
     w oknie 1,75–2,00 GeV → N_fake ± σ, pull = N_fake/σ, oraz N_fake/N_tła w ±55 MeV; to samo dla dalekich+bliskich SB danych nie jest
     możliwe (sygnał w oknie) → dla danych tylko R pasma (sculpt.py).
 (2) kształt: ε_B(m) w 25 binach, χ²/ndf dopasowania liniowego i kwadratowego w całym oknie (tło MC), nachylenie.
 (3) regresja masy: R² (OOF, HistGB) przewidywania m z cech na tle MC truth≥8 (całe okno) — H vs H+grupa: ile informacji o masie niosą cechy.
Wyjście: out/d0_sculpt.json.
"""
import json
import os
import sys

import numpy as np
from iminuit import Minuit
from scipy.stats import norm
from sklearn.ensemble import HistGradientBoostingRegressor

sys.argv = sys.argv[:1]
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.argv = ["d0_run.py", "none"]
import d0_run  # noqa: E402  (ustawia BASE, grupy)
study = d0_run.study
from sculpt import thr_for, passed  # noqa: E402

D0 = d0_run.BASE
sp = "D0"
s = study.cfg["species"][sp]
M0, LO, HI = s["m0"], s["mass_range"][0], s["mass_range"][1]
E = np.linspace(LO, HI, 51)  # 5 MeV


def fake_fit(h, sig_sigma):
    c = 0.5 * (E[1:] + E[:-1])
    x = np.linspace(LO, HI, 1001)
    t = (x - M0) / 0.1

    def lam(N, nb, b1, b2):
        fb = np.exp(b1 * t + b2 * t * t)
        cb = np.concatenate([[0], np.cumsum(0.5 * (fb[1:] + fb[:-1]) * np.diff(x))])
        return np.clip(N * np.diff(norm.cdf(E, M0, sig_sigma)) + nb * np.diff(np.interp(E, x, cb / cb[-1])), 1e-9, None)

    def nll(N, nb, b1, b2):
        l_ = lam(N, nb, b1, b2)
        return float(np.sum(l_ - h * np.log(l_)))
    m = Minuit(nll, N=0.0, nb=h.sum(), b1=-0.3, b2=0.0)
    m.errordef = 0.5
    m.limits["nb"] = (0, None)
    m.migrad()
    m.hesse()
    win = np.abs(c - M0) < s["signal_window"]
    l_ = lam(*m.values)
    return dict(N=float(m.values["N"]), err=float(m.errors["N"]), pull=float(m.values["N"] / m.errors["N"]), n_win=float(h[win].sum()),
                frac=float(m.values["N"] / max(h[win].sum(), 1)), chi2_ndf=float(np.sum((h - l_) ** 2 / np.maximum(l_, 1)) / (len(h) - 4)), valid=bool(m.valid))


def main():
    S, Bmc, Bda, dall, mb = study.load(sp)
    odir = os.path.join(D0, "out", "D0_mc")
    mS = S["m"].to_numpy(float)
    wS = S["weight"].to_numpy(float)
    sig_sigma = float(np.sqrt(np.average((mS[np.abs(mS - M0) < 0.06] - M0) ** 2, weights=wS[np.abs(mS - M0) < 0.06])))
    mB, ptB = Bmc["m"].to_numpy(float), Bmc["pt"].to_numpy(float)
    res = {"sig_sigma": sig_sigma, "fake": {}, "shape": {}, "massreg": {}}
    e25 = np.linspace(LO, HI, 26)
    c25 = 0.5 * (e25[1:] + e25[:-1])
    for f in sorted(os.listdir(odir)):
        if not (f.startswith("scores_") and f.endswith(".npz")):
            continue
        n = f[7:-4]
        z = np.load(os.path.join(odir, f))
        res["fake"][n], res["shape"][n] = {}, {}
        for wp in (0.7, 0.5):
            t = thr_for(z["S"], wS, S["pt"].to_numpy(float), wp)
            ok = passed(z["B"], ptB, t)
            for lo, hi in study.BINS + [(3, 30)]:
                a = (ptB >= lo) & (ptB < hi)
                h, _ = np.histogram(mB[a & ok], E)
                res["fake"][n][f"WP{int(wp*100)}|{lo}-{hi}"] = fake_fit(h.astype(float), sig_sigma)
            n0, _ = np.histogram(mB, e25)
            n1, _ = np.histogram(mB[ok], e25)
            e = n1 / np.maximum(n0, 1)
            err = np.sqrt(np.maximum(e * (1 - e), 1e-12) / np.maximum(n0, 1))
            out = {"eps": e.tolist(), "err": err.tolist()}
            for deg in (1, 2):
                A = np.vander(c25 - M0, deg + 1, increasing=True)
                W = 1 / err**2
                beta = np.linalg.solve(A.T @ (A * W[:, None]), A.T @ (W * e))
                out[f"chi2_pol{deg}"] = float(np.sum(W * (e - A @ beta) ** 2))
                out[f"ndf_pol{deg}"] = 25 - deg - 1
                if deg == 1:
                    out["slope_rel_per_100MeV"] = float(beta[1] * 0.1 / beta[0])
            res["shape"][n][f"WP{int(wp*100)}"] = out
        fk = res["fake"][n]["WP70|3-30"]
        print(n, "WP70 fake N=%.0f±%.0f pull=%.1f frac=%.4f" % (fk["N"], fk["err"], fk["pull"], fk["frac"]),
              {k: round(v["pull"], 1) for k, v in res["fake"][n].items()}, "chi2 pol1/pol2 WP70:", round(res["shape"][n]["WP70"]["chi2_pol1"], 1), round(res["shape"][n]["WP70"]["chi2_pol2"], 1), flush=True)
    # regresja masy na tle MC
    rng = np.random.default_rng(5)
    idx = np.sort(rng.choice(len(Bmc), min(200000, len(Bmc)), replace=False))
    y = mB[idx]
    fold = Bmc["fold"].to_numpy()[idx]
    sets = {"H": study.H}
    for k, v in d0_run.NEW.items():
        sets["H+" + k] = study.H + v
    for k, v in d0_run.NEW.items():
        if k in ("kin2", "kinall", "ip", "time", "pidbb"):
            sets["only_" + k] = v
    for nm, feats in sets.items():
        X = np.nan_to_num(Bmc[feats].to_numpy(np.float32)[idx], nan=-9)
        pred = np.zeros(len(y))
        for k in range(5):
            r = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.1, max_leaf_nodes=31, random_state=0)
            r.fit(X[fold != k], y[fold != k])
            pred[fold == k] = r.predict(X[fold == k])
        r2 = 1 - np.mean((y - pred) ** 2) / np.var(y)
        res["massreg"][nm] = dict(r2=float(r2), rmse_MeV=float(1e3 * np.sqrt(np.mean((y - pred) ** 2))), rms_m_MeV=float(1e3 * np.std(y)))
        print("massreg", nm, res["massreg"][nm], flush=True)
    json.dump(res, open(os.path.join(D0, "out", "d0_sculpt.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
