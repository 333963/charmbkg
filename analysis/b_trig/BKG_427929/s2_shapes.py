"""Krok 2: kształty tła w pasmach bocznych (|m−m0| ≥ 75 MeV, okno 1,75–2,00) dane vs MC (inkluzywne mb_nd, bez wag) w obu runach
+ dane–dane; KS, χ² w binach, stosunki w binach score; C2ST (HGB, 5-fold GroupKFold po zdarzeniu).
Część b: pary D0←D* (18 cech H + η) w paśmie bliskim 55–75 MeV, Δm > 150 MeV — jedyne miejsce z pełnymi cechami.
Wyjście: out/s2.json, out/s2.md."""
import json, os, sys
import numpy as np
from scipy import stats
from sklearn.ensemble import HistGradientBoostingClassifier as HGB
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
import bkg_common as C

rng = np.random.default_rng(2)
NMAX = 100_000
S = {t: C.load_all(t) for t in ("d427", "d428", "mb")}
for D in S.values():
    D["sb"] = np.abs(D["m"] - C.M0) >= C.SB_MIN
    D["gid"] = D["run"].astype(np.int64) * 100_000_000_003 + D["event"].astype(np.int64)
PAIRS = {"427929 vs MC": ("d427", "mb"), "428770 vs MC": ("d428", "mb"), "427929 vs 428770": ("d427", "d428")}
REG = {"3-30": (3.0, 30.0), **C.BINS}
SHAPEV = ["score", "pt", "ptcone40", "ptcone100", "e_nSelBL", "e_lead_pt", "e_ntrk500", "m"]


def chi2_shape(a, b, edges):
    ha, hb = np.histogram(a, edges)[0].astype(float), np.histogram(b, edges)[0].astype(float)
    na, nb = ha.sum(), hb.sum()
    ok = (ha + hb) > 0
    c2 = np.sum((ha[ok] / na - hb[ok] / nb) ** 2 / (ha[ok] / na ** 2 + hb[ok] / nb ** 2))
    return float(c2), int(ok.sum() - 1), ha, hb


def c2st(Xa, Xb, ga, gb, nmax=NMAX, seed=0):
    r = np.random.default_rng(seed)
    n = min(len(Xa), len(Xb), nmax)
    ia, ib = r.choice(len(Xa), n, replace=False), r.choice(len(Xb), n, replace=False)
    X = np.vstack([Xa[ia], Xb[ib]]).astype(np.float32)
    y = np.r_[np.ones(n), np.zeros(n)]
    g = np.r_[ga[ia], gb[ib] + 7]                      # różne próbki → różne grupy
    oof = np.zeros(2 * n)
    aucs = []
    for tr, te in GroupKFold(5).split(X, y, g):
        clf = HGB(max_iter=100, learning_rate=0.1, max_leaf_nodes=31, min_samples_leaf=50, random_state=seed)
        clf.fit(X[tr], y[tr])
        oof[te] = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], oof[te]))
    return {"n_per_class": int(n), "auc": float(np.mean(aucs)), "err_folds": float(np.std(aucs, ddof=1) / np.sqrt(5)),
            "auc_pooled": float(roc_auc_score(y, oof))}


FS = {"score": ["score"], "kand (score+ptcone)": ["score", "ptcone40", "ptcone100"], "zdarzenie (nSelBL, lead_pt, ntrk500)": ["e_nSelBL", "e_lead_pt", "e_ntrk500"],
      "wszystkie bez score": C.VARS, "wszystkie + score": C.VARS + ["score"]}
res = {"shape": {}, "c2st": {}, "ratio_score": {}}
sedges = np.linspace(0, 1, 21)
for pn, (a, b) in PAIRS.items():
    A, B = S[a], S[b]
    for rn, (lo, hi) in REG.items():
        sa = A["sb"] & (A["pt"] >= lo) & (A["pt"] < hi)
        sb = B["sb"] & (B["pt"] >= lo) & (B["pt"] < hi)
        for v in SHAPEV:
            xa, xb = A[v][sa].astype(float), B[v][sb].astype(float)
            ks = stats.ks_2samp(xa, xb)
            edges = sedges if v == "score" else np.unique(np.quantile(np.r_[xa[:: max(1, len(xa) // 200000)], xb[:: max(1, len(xb) // 200000)]], np.linspace(0, 1, 21)))
            edges = edges.copy(); edges[0], edges[-1] = -np.inf, np.inf
            c2, ndf, ha, hb = chi2_shape(xa, xb, edges)
            res["shape"][f"{pn}|{rn}|{v}"] = {"ks": float(ks.statistic), "ks_p": float(ks.pvalue), "chi2": c2, "ndf": ndf, "na": int(len(xa)), "nb": int(len(xb)),
                                              "mean_a": float(np.mean(xa)), "mean_b": float(np.mean(xb)), "med_a": float(np.median(xa)), "med_b": float(np.median(xb))}
            if v == "score":
                ra = (ha / ha.sum()) / (hb / hb.sum())
                er = ra * np.sqrt(1 / np.maximum(ha, 1) + 1 / np.maximum(hb, 1))
                res["ratio_score"][f"{pn}|{rn}"] = {"ratio": ra.tolist(), "err": er.tolist()}
        for fn, fs in FS.items():
            Xa, Xb = np.column_stack([A[f][sa] for f in fs]), np.column_stack([B[f][sb] for f in fs])
            r = c2st(Xa, Xb, A["gid"][sa], B["gid"][sb])
            res["c2st"][f"{pn}|{rn}|{fn}"] = r
            print(pn, rn, fn, r, flush=True)
# wariant: MC truth ≥ 8 zamiast inkluzywnego (score, bin 3–5 i 5–8) — wpływ składu
mbD = S["mb"]
for rn, (lo, hi) in REG.items():
    for a in ("d427", "d428"):
        A = S[a]
        sa = A["sb"] & (A["pt"] >= lo) & (A["pt"] < hi)
        sb = mbD["sb"] & (mbD["pt"] >= lo) & (mbD["pt"] < hi) & (mbD["truth"] >= 8)
        ks = stats.ks_2samp(A["score"][sa], mbD["score"][sb])
        c2, ndf, ha, hb = chi2_shape(A["score"][sa], mbD["score"][sb], np.r_[-np.inf, sedges[1:-1], np.inf])
        res["shape"][f"{a} vs MC t8|{rn}|score"] = {"ks": float(ks.statistic), "chi2": c2, "ndf": ndf}

# ---------------------------------------------------------------- część b: pary, 18 cech
F18 = ["log_pt", "abs_eta", "chi2", "Lxy", "impact_sig", "abs_impact_z0", "d_max_abs_d0", "pointing", "abs_cos_theta_star", "abs_impact_sig",
       "d_min_abs_d0", "d_max_abs_z0", "d_max_abs_eta", "dedx_K", "log_p_K", "npts_K", "dedx_pi_max", "log_p_pi_min"]


def load_pairs(tag):
    z = np.load(os.path.join(C.RUNS, C.FILES[tag]), allow_pickle=True)
    cols = ["run", "event", "file_idx", "row", "m", "pt", "eta", "score", "dm", "ptcone40", "ptcone100", "e_nSelBL", "e_lead_pt", "e_ntrk500", "e_grl", "truth"] + F18
    P = {c: z[c] for c in cols}
    keep = (np.abs(P["m"] - C.M0) >= 0.055) & (P["dm"] > 0.150)
    if tag == "d427":
        keep &= (P["e_grl"] == 1) & (z["e_tp_HLT_mb_mbts_L1MBTS_1"] == 1)
    P = {k: v[keep] for k, v in P.items()}
    key = np.stack([P["run"], P["event"], P["file_idx"], P["row"]], 1)
    _, i = np.unique(key, axis=0, return_index=True)                 # jeden wiersz na kandydata D0
    P = {k: v[i] for k, v in P.items()}
    P["gid"] = P["run"].astype(np.int64) * 100_000_000_003 + P["event"].astype(np.int64) + P["file_idx"].astype(np.int64) * 1_000_003
    return P


PP = {t: load_pairs(t) for t in ("d427", "d428", "mb")}
res["pairs_n"] = {t: int(len(P["m"])) for t, P in PP.items()}
FSP = {"score": ["score"], "18 cech H": F18, "18 cech + ptcone + zdarzenie": F18 + ["ptcone40", "ptcone100", "e_nSelBL", "e_lead_pt", "e_ntrk500"],
       "18 cech bez log_pt": [f for f in F18 if f != "log_pt"]}
for pn, (a, b) in PAIRS.items():
    A, B = PP[a], PP[b]
    for rn, (lo, hi) in {"3-30": (3.0, 30.0), "3-5": (3.0, 5.0), "5-30": (5.0, 30.0)}.items():
        sa, sb = (A["pt"] >= lo) & (A["pt"] < hi), (B["pt"] >= lo) & (B["pt"] < hi)
        for fn, fs in FSP.items():
            r = c2st(np.column_stack([A[f][sa] for f in fs]), np.column_stack([B[f][sb] for f in fs]), A["gid"][sa], B["gid"][sb])
            res["c2st"][f"pary|{pn}|{rn}|{fn}"] = r
            print("pary", pn, rn, fn, r, flush=True)
        if rn == "3-30":
            for v in ["eta", "score"] + F18:
                ks = stats.ks_2samp(A[v][sa].astype(float), B[v][sb].astype(float))
                res["shape"][f"pary|{pn}|{v}"] = {"ks": float(ks.statistic), "ks_p": float(ks.pvalue), "na": int(sa.sum()), "nb": int(sb.sum())}
json.dump(res, open(os.path.join(C.OUT, "s2.json"), "w"), indent=1)
