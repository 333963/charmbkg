"""FEAT_charmnet: rzeźbienie masy tła po cięciu WP (progi per bin pT z sygnału OOF).
R_sculpt = ε_B(okno) obserwowane / ε_B(okno) z dopasowania liniowego ε_B(m) w sidebandach (|Δm| > sideband_min).
Tło: mb truth≥8 (prawda MC, oba testy) oraz dane (test data: OOF; test mc: model MC na danych) — dla danych okno
zawiera sygnał, więc dla danych: fit liniowy w dalekich SB (>100/110 MeV) → pasmo [okno, 100/110 MeV) (bez sygnału poza ogonami).
Użycie: python sculpt.py SPECIES TEST
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from study import BASE, BINS, cfg, load  # noqa: E402


def thr_for(sS, wS, ptS, wp):
    t = {}
    for lo, hi in BINS:
        a = (ptS >= lo) & (ptS < hi)
        s, w = sS[a], wS[a]
        o = np.argsort(-s)
        c = np.cumsum(w[o]) / w.sum()
        t[(lo, hi)] = s[o][min(np.searchsorted(c, wp), len(o) - 1)]
    return t


def passed(s, pt, t):
    ok = np.zeros(len(s), bool)
    for (lo, hi), v in t.items():
        ok |= (pt >= lo) & (pt < hi) & (s >= v)
    return ok


def eps_m(m, ok, edges):
    n0, _ = np.histogram(m, edges)
    n1, _ = np.histogram(m[ok], edges)
    e = n1 / np.maximum(n0, 1)
    err = np.sqrt(np.maximum(e * (1 - e), 1e-12) / np.maximum(n0, 1))
    return e, err, n0, n1


def closure(m, ok, sp, region_lo, region_hi, fitmin=None):
    """linia w SB (|Δm|>sideband_min) → przewidywanie dla pasma |Δm| in [region_lo, region_hi); zwraca obs/pred ± błąd."""
    s = cfg["species"][sp]
    m0, sbm = s["m0"], (fitmin or s["sideband_min"])
    lo, hi = s["mass_range"]
    edges = np.linspace(lo, hi, 26)
    e, err, n0, n1 = eps_m(m, ok, edges)
    c = 0.5 * (edges[1:] + edges[:-1])
    sb = (np.abs(c - m0) > sbm) & (n0 > 20)
    A = np.c_[np.ones(sb.sum()), c[sb] - m0]
    W = 1 / err[sb] ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    beta = cov @ (A.T @ (W * e[sb]))
    d = np.abs(m - m0)
    reg = (d >= region_lo) & (d < region_hi)
    obs = ok[reg].mean()
    x = np.c_[np.ones(reg.sum()), m[reg] - m0]
    pred = float((x @ beta).mean())
    vp = float(np.mean(x, 0) @ cov @ np.mean(x, 0))
    eo = np.sqrt(obs * (1 - obs) / reg.sum())
    r = obs / pred
    return dict(obs=float(obs), pred=pred, R=float(r), R_err=float(r * np.hypot(eo / obs, np.sqrt(vp) / pred)),
                slope_rel_per_100MeV=float(beta[1] * 0.1 / beta[0]), eps_m=e.tolist(), edges=edges.tolist())


def main(sp, test):
    S, Bmc, Bda, dall, mb = load(sp)
    s = cfg["species"][sp]
    odir = os.path.join(BASE, "out", f"{sp}_{test}")
    res = {}
    for f in sorted(os.listdir(odir)):
        if not (f.startswith("scores_") and f.endswith(".npz")):
            continue
        n = f[7:-4]
        z = np.load(os.path.join(odir, f))
        r = {}
        for wp in (0.8, 0.7, 0.5):
            t = thr_for(z["S"], S["weight"].to_numpy(float), S["pt"].to_numpy(float), wp)
            sB = z["B"] if test == "mc" else z["mb"]
            okm = passed(sB, Bmc["pt"].to_numpy(float), t)
            rr = {"mc_truth8_SR": closure(Bmc["m"].to_numpy(float), okm, sp, 0, s["signal_window"])}
            if "dall" in z.files:
                okd = passed(z["dall"], dall["pt"].to_numpy(float), t)
                far = {"D0": 0.100, "Dplus": 0.110}[sp]
                rr["data_band"] = closure(dall["m"].to_numpy(float), okd, sp, s["signal_window"], far, fitmin=far)
            for k in rr:
                rr[k] = {q: v for q, v in rr[k].items() if (wp == 0.7 or q not in ("eps_m", "edges"))}
            r[f"WP{int(wp*100)}"] = rr
        res[n] = r
        print(sp, test, n, {w: {k: (round(v["R"], 3), round(v["R_err"], 3), round(v["slope_rel_per_100MeV"], 3)) for k, v in q.items()} for w, q in r.items()}, flush=True)
    json.dump(res, open(os.path.join(odir, "sculpt.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
