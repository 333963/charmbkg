#!/usr/bin/env python3
"""Krok 4: rozkłady aktywności (nSelBL, lead_pt, ntrk500, ptcone40/100, pT) per klasa trygera vs 427929 vs MC:
(a) zdarzenia z kandydatem (ev_*), (b) kandydaci tła = pasma boczne masy inkluzywnego D0 (all_*), (c) kandydaci sygnału D* (wagi SBS).
C2ST (HGB, 5-fold GroupKFold po zdarzeniu) klasa vs 427929, vs MC (mb_nd), vs J12only — na pasmach bocznych. Wynik: out/s4.json."""
import numpy as np
from scipy import stats
from sklearn.ensemble import HistGradientBoostingClassifier as HGB
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
import afp_common as A
import bkg_common as K
from iso_common import load_pairs, sbs_w, C
import s2_universal_rw as U
from s1_iso import samples

NMAX = 60000
VARS = ["pt", "ptcone40", "ptcone100", "e_nSelBL", "e_lead_pt", "e_ntrk500"]
FS = {"zdarzenie(3)": ["e_nSelBL", "e_lead_pt", "e_ntrk500"], "kandydat(3)": ["pt", "ptcone40", "ptcone100"], "wszystkie(6)": VARS, "wszystkie+score(7)": VARS + ["score"]}
CL = ["AFPnoalg", "AFPnoalg_any", "AFPsptrk4", "AFPboth", "AFP&J12", "AFP&MBTS", "AFPonly", "AFPstream", "none", "MBTS!A", "J12only!A"]


def c2st(Xa, Xb, ga, gb, seed=0):
    r = np.random.default_rng(seed)
    n = min(len(Xa), len(Xb), NMAX)
    ia, ib = r.choice(len(Xa), n, replace=False), r.choice(len(Xb), n, replace=False)
    X = np.vstack([Xa[ia], Xb[ib]]).astype(np.float32)
    y = np.r_[np.ones(n), np.zeros(n)]
    _, ga2 = np.unique(ga[ia], return_inverse=True); _, gb2 = np.unique(gb[ib], return_inverse=True)
    g = np.r_[ga2, gb2 + ga2.max() + 1]
    oof = np.zeros(2 * n); aucs = []
    for tr, te in GroupKFold(5).split(X, y, g):
        clf = HGB(max_iter=100, learning_rate=0.1, max_leaf_nodes=15, min_samples_leaf=max(20, min(50, n // 40)), random_state=seed)
        clf.fit(X[tr], y[tr]); oof[te] = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], oof[te]))
    return {"n": int(n), "auc": float(np.mean(aucs)), "err": float(np.std(aucs, ddof=1) / np.sqrt(5))}


def load_all(z, extra_mask=None):
    D = {v: z["all_" + v] for v in VARS + ["score", "m", "event", "run"]}
    D["sb"] = np.abs(D["m"] - K.M0) >= K.SB_MIN
    return D


def stat(x, w=None):
    x = np.asarray(x, float)
    if w is None:
        return {"n": int(len(x)), "mean": float(x.mean()), "med": float(np.median(x)), "q25": float(np.quantile(x, .25)), "q75": float(np.quantile(x, .75)),
                "mean_err": float(x.std() / np.sqrt(max(len(x), 1)))}
    T = w.sum(); mu = (w * x).sum() / T
    return {"n": float(T), "mean": float(mu), "mean_err": float(np.sqrt((w ** 2 * (x - mu) ** 2).sum()) / T)}


def main():
    z8 = np.load(f"{A.R}/bt_data_428770.npz", allow_pickle=True)
    z9 = np.load(f"{A.R}/bt_data_427929.npz", allow_pickle=True)
    zm = np.load(f"{A.R}/bt_mb_nd.npz", allow_pickle=True)
    res = {"ev": {}, "bkg": {}, "sig": {}, "c2st": {}, "ks": {}, "hist": {}}
    # (a) zdarzenia z kandydatem
    ce = A.classes8(A.bits(z8, "ev_")); ce["AFPnoalg_any"] = ce["_an"]
    ce9 = A.classes9(A.bits(z9, "ev_"), z9["ev_grl"] == 1)
    for k in CL + ["J12only", "MBTSbased", "other", "all"]:
        res["ev"]["8:" + k] = {v: stat(z8["ev_" + v][ce[k]]) for v in ("nSelBL", "lead_pt", "ntrk500")}
        res["ev"]["8:" + k]["frac_lead_lt4"] = float((z8["ev_lead_pt"][ce[k]] < 4).mean())
        res["ev"]["8:" + k]["frac_npu0"] = float((z8["ev_npu"][ce[k]] == 0).mean())
    for k in ("MB427929", "MB9&AFPnoalg", "MB9&AFPsptrk4"):
        res["ev"]["9:" + k] = {v: stat(z9["ev_" + v][ce9[k]]) for v in ("nSelBL", "lead_pt", "ntrk500")}
        res["ev"]["9:" + k]["frac_lead_lt4"] = float((z9["ev_lead_pt"][ce9[k]] < 4).mean())
        res["ev"]["9:" + k]["frac_npu0"] = float((z9["ev_npu"][ce9[k]] == 0).mean())
    res["ev"]["MC mb_nd"] = {v: stat(zm["ev_" + v]) for v in ("nSelBL", "lead_pt", "ntrk500")} if "ev_nSelBL" in zm.files else {}
    # (b) tło: SB inkluzywnego D0
    D8, D9a, DM = load_all(z8), load_all(z9), load_all(zm)
    c8 = A.classes8(A.bits(z8, "all_e_")); c8["AFPnoalg_any"] = c8["_an"]
    c9 = A.classes9(A.bits(z9, "all_e_"), z9["all_e_grl"] == 1)
    REF = {"427929": (D9a, c9["MB427929"] & D9a["sb"]), "MC": (DM, DM["sb"]), "J12only": (D8, c8["J12only!A"] & D8["sb"]),
           "427929 j20&!J12": (D9a, c9["MB9&j20&!J12"] & D9a["sb"])}
    HED = {"e_nSelBL": np.arange(0, 121, 5.0), "e_lead_pt": np.r_[np.arange(0, 12, 1.0), 14, 17, 20, 25, 30, 40], "ptcone40": np.r_[np.arange(0, 10, 0.5), 12, 15, 20],
           "ptcone100": np.r_[np.arange(0, 30, 2.0), 35, 40, 50], "pt": np.r_[np.arange(3, 8, 0.5), 9, 10, 12, 15, 20, 30], "e_ntrk500": np.arange(0, 161, 8.0)}
    def hh(D, s):
        return {v: np.histogram(D[v][s], ed)[0].tolist() for v, ed in HED.items()}
    res["hist"]["edges"] = {v: ed.tolist() for v, ed in HED.items()}
    for rn, (D, s) in REF.items():
        res["bkg"][rn] = {v: stat(D[v][s]) for v in VARS + ["score"]}
        res["hist"]["bkg|" + rn] = hh(D, s)
    for k in CL:
        s = c8[k] & D8["sb"]
        res["bkg"]["8:" + k] = {v: stat(D8[v][s]) for v in VARS + ["score"]}
        res["hist"]["bkg|8:" + k] = hh(D8, s)
        for rn, (D, sr) in REF.items():
            if rn == "J12only" and k == "J12only!A":
                continue
            res["ks"][f"{k}|{rn}"] = {v: float(stats.ks_2samp(D8[v][s], D[v][sr][:: max(1, sr.sum() // 300000)]).statistic) for v in VARS + ["score"]}
            if s.sum() < 400:
                continue
            for fn, fs in FS.items():
                r = c2st(np.column_stack([D8[f][s] for f in fs]), np.column_stack([D[f][sr] for f in fs]), D8["event"][s], D["event"][sr] * 1000 + D["run"][sr] % 1000)
                res["c2st"][f"{k}|{rn}|{fn}"] = r
                print(k, rn, fn, r, flush=True)
    # odniesienia między próbkami
    for a, b in (("427929", "MC"), ("J12only", "427929"), ("J12only", "MC"), ("427929 j20&!J12", "427929"), ("427929 j20&!J12", "MC")):
        Da, sa = REF[a]; Db, sb = REF[b]
        res["ks"][f"{a}|{b}"] = {v: float(stats.ks_2samp(Da[v][sa][:: max(1, sa.sum() // 300000)], Db[v][sb][:: max(1, sb.sum() // 300000)]).statistic) for v in VARS + ["score"]}
        for fn, fs in FS.items():
            r = c2st(np.column_stack([Da[f][sa] for f in fs]), np.column_stack([Db[f][sb] for f in fs]), Da["event"][sa], Db["event"][sb] * 1000 + Db["run"][sb] % 1000)
            res["c2st"][f"{a}|{b}|{fn}"] = r
            print(a, b, fn, r, flush=True)
    # test zerowy C2ST: dwie losowe połowy 427929 o liczności klasy AFPsptrk4
    Da, sa = REF["427929"]; idx = np.flatnonzero(sa); rr = np.random.default_rng(1); rr.shuffle(idx)
    n0 = int((c8["AFPsptrk4"] & D8["sb"]).sum())
    h1, h2 = idx[:n0], idx[n0:2 * n0]
    res["c2st"]["NULL 427929 polowy|wszystkie(6)"] = c2st(np.column_stack([Da[f][h1] for f in VARS]), np.column_stack([Da[f][h2] for f in VARS]), Da["event"][h1], Da["event"][h2] + 10 ** 12)
    print("NULL", res["c2st"]["NULL 427929 polowy|wszystkie(6)"])
    # (c) sygnał D* (wagi SBS) — średnie
    P8, P9, S = samples()
    M = load_pairs("bt_fd_nd"); mm = U.mc_sel(M, 0.040, 0.6)
    SV = ["pt", "ptcone40", "ptcone100", "e_nSelBL", "e_lead_pt", "e_ntrk500", "rel40"]
    res["sig"]["MC fd_nd"] = {v: stat(M[v][mm]) for v in SV}
    for name in ["8:AFPstream", "8:AFPonly", "8:AFPsptrk4", "8:AFPnoalg_any", "8:AFP&J12", "8:none", "8:MBTS!A", "8:J12only!A", "9:MB427929", "9:MB9&j20&!J12"]:
        D, cm = S[name]; m = U.base(D, 0.040, 0.6) & cm
        w = sbs_w(D["dm"][m], D["pt"][m], np.ones(int(m.sum()), bool))
        res["sig"][name] = {v: stat(D[v][m], w) for v in SV}
        res["hist"]["sig|" + name] = {v: {"h": np.histogram(D[v][m], ed, weights=w)[0].tolist(), "h2": np.histogram(D[v][m], ed, weights=w ** 2)[0].tolist()} for v, ed in HED.items()}
    res["hist"]["sig|MC"] = {v: np.histogram(M[v][mm], ed)[0].tolist() for v, ed in HED.items()}
    A.jdump(res, f"{A.OUT}/s4.json")
    print("zapisano out/s4.json")


if __name__ == "__main__":
    main()
