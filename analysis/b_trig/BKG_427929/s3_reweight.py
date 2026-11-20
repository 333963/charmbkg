"""Krok 3: przeważenie MC „SALAD-podobne”: klasyfikator dane-SB vs MC-SB (HGB, 5-fold GroupKFold po zdarzeniu), w = p/(1−p),
trenowany WYŁĄCZNIE na pasmach bocznych (|m−m0| ≥ 75 MeV); wagi dla MC w SR = średnia (logit) z 5 modeli (MC-SR nie bierze udziału
w treningu), dla MC w SB = out-of-fold. Test: ε_B(SR) przewidziane (MC ważone) vs zmierzone (interpolacja pol2 danych, s1).
Wyjście: out/s3.json, out/s3.md, out/s3_w_<run>.npz (wagi)."""
import json, os
import numpy as np
from scipy import stats
from sklearn.ensemble import HistGradientBoostingClassifier as HGB, HistGradientBoostingRegressor as HGR
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, r2_score
import bkg_common as C

NREP = 400
NDATA = 400_000
CLIP = (0.05, 20.0)
thr = C.thresholds()
rep1 = np.load(os.path.join(C.OUT, "s1_rep.npz"))
s1 = json.load(open(os.path.join(C.OUT, "s1.json")))
mb = C.load_all("mb")
mb["gid"] = mb["event"].astype(np.int64) + mb["run"].astype(np.int64) * 100_000_000_003
mb_sb = np.abs(mb["m"] - C.M0) >= C.SB_MIN
mb_sr = np.abs(mb["m"] - C.M0) < C.SR_HALF
FSETS = {
    "V0 pT+zdarzenie": ["pt", "e_nSelBL", "e_lead_pt", "e_ntrk500"],
    "V1 bez score": C.VARS,
    "V2 + score": C.VARS + ["score"],
    "V3 + score + m": C.VARS + ["score", "m"],
}


def wks(x1, x2, w2):
    """KS: próbka 1 nieważona vs próbka 2 ważona."""
    o1, o2 = np.sort(x1), np.argsort(x2)
    x2s, c2 = x2[o2], np.cumsum(w2[o2]) / w2.sum()
    grid = np.r_[o1[:: max(1, len(o1) // 5000)], x2s[:: max(1, len(x2s) // 5000)]]
    f1 = np.searchsorted(o1, grid, side="right") / len(o1)
    f2 = np.r_[0, c2][np.searchsorted(x2s, grid, side="right")]
    return float(np.abs(f1 - f2).max())


def X_of(D, fs, sel):
    return np.column_stack([D[f][sel] for f in fs]).astype(np.float32)


res = {}
rng = np.random.default_rng(3)
# wyciek masy: czy cechy przewidują m na tle MC (truth ≥ 8, całe okno)?
leak = {}
bk = mb["truth"] >= 8
for fn, fs in FSETS.items():
    if "m" in fs:
        continue
    X, y, g = X_of(mb, fs, bk), mb["m"][bk], mb["gid"][bk]
    oof = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, y, g):
        oof[te] = HGR(max_iter=200, learning_rate=0.1, min_samples_leaf=100, random_state=0).fit(X[tr], y[tr]).predict(X[te])
    leak[fn] = {"r2_oof": float(r2_score(y, oof)), "corr": float(np.corrcoef(y, oof)[0, 1])}
    print("leak", fn, leak[fn], flush=True)
res["mass_leak"] = leak

for run in ("d427", "d428"):
    D = C.load_all(run)
    D["gid"] = D["event"].astype(np.int64) + D["run"].astype(np.int64) * 100_000_000_003 + 11
    d_sb = np.flatnonzero(np.abs(D["m"] - C.M0) >= C.SB_MIN)
    d_tr = rng.choice(d_sb, min(NDATA, len(d_sb)), replace=False)
    m_tr = np.flatnonzero(mb_sb)
    m_sr = np.flatnonzero(mb_sr)
    wsave = {}
    for fn, fs in FSETS.items():
        X = np.vstack([X_of(D, fs, d_tr), X_of(mb, fs, m_tr)])
        y = np.r_[np.ones(len(d_tr)), np.zeros(len(m_tr))]
        g = np.r_[D["gid"][d_tr], mb["gid"][m_tr]]
        Xsr = X_of(mb, fs, m_sr)
        if "m" in fs:      # SR: prawdziwa masa kandydata — drzewa interpolują wagę przez przerwę 1,79–1,94 (stała kawałkami)
            pass
        oof = np.zeros(len(y))
        lsr = np.zeros(len(m_sr))
        for tr, te in GroupKFold(5).split(X, y, g):
            clf = HGB(max_iter=200, learning_rate=0.1, max_leaf_nodes=31, min_samples_leaf=100, l2_regularization=1.0, random_state=1).fit(X[tr], y[tr])
            oof[te] = clf.predict_proba(X[te])[:, 1]
            p = np.clip(clf.predict_proba(Xsr)[:, 1], 1e-6, 1 - 1e-6)
            lsr += np.log(p / (1 - p)) / 5
        auc = float(roc_auc_score(y, oof))
        prior = len(m_tr) / len(d_tr)
        p_sb = np.clip(oof[len(d_tr):], 1e-6, 1 - 1e-6)
        w_sb = p_sb / (1 - p_sb) * prior
        w_sr = np.exp(lsr) * prior
        norm = w_sb.mean()
        w_sb, w_sr = np.clip(w_sb / norm, *CLIP), np.clip(w_sr / norm, *CLIP)
        wsave[fn + "|sb"], wsave[fn + "|sr"] = w_sb.astype(np.float32), w_sr.astype(np.float32)
        R = {"auc_sb": auc, "n_data_train": int(len(d_tr)), "n_mc_train": int(len(m_tr)), "clip_frac": float(np.mean((w_sb <= CLIP[0]) | (w_sb >= CLIP[1]))),
             "neff_sb": float(w_sb.sum() ** 2 / (w_sb ** 2).sum() / len(w_sb)), "neff_sr": float(w_sr.sum() ** 2 / (w_sr ** 2).sum() / len(w_sr)),
             "w_mean_sr_over_sb": float(w_sr.mean() / w_sb.mean()), "corr_w_m_sb": float(np.corrcoef(w_sb, mb["m"][m_tr])[0, 1]), "bins": {}}
        c2 = {"t8": [0.0, 0.0], "t2": [0.0, 0.0]}
        for b, (lo, hi) in C.BINS.items():
            T = np.array(thr[b])
            pt_sr, pt_sb = mb["pt"][m_sr], mb["pt"][m_tr]
            e_data = np.array(s1["eps"][f"{run}|{b}"]["pol2"])
            e_data_sym = np.array(s1["eps"][f"{run}|{b}"]["sym"])
            rd = rep1[f"{run}|{b}|pol2"]
            B = {"eps_data_pol2": e_data.tolist(), "stat_data_rel": (rd.std(0) / e_data).tolist()}
            # SB closure (OOF): MC inkl. ważone vs dane, zliczenie w SB (ten sam zakres co `sym`)
            sb_b = (pt_sb >= lo) & (pt_sb < hi) & (np.abs(mb["m"][m_tr] - C.M0) <= C.SB_MAX)
            sc_sb = mb["score"][m_tr][sb_b]
            ww = w_sb[sb_b]
            B["sb_mc_w"] = [(ww * (sc_sb >= t)).sum() / ww.sum() for t in T]
            B["sb_mc_0"] = [float((sc_sb >= t).mean()) for t in T]
            B["sb_data"] = e_data_sym.tolist()
            dsel = (D["pt"][d_sb] >= lo) & (D["pt"][d_sb] < hi)
            B["ks_score_sb_before"] = wks(D["score"][d_sb][dsel], sc_sb, np.ones(len(sc_sb)))
            B["ks_score_sb_after"] = wks(D["score"][d_sb][dsel], sc_sb, ww)
            for tn, tmin in (("t8", 8), ("t2", 2)):
                s = (pt_sr >= lo) & (pt_sr < hi) & (mb["truth"][m_sr] >= tmin)
                sc, w = mb["score"][m_sr][s], w_sr[s]
                P = (sc[:, None] >= T[None, :]).astype(float)
                e0, ew = P.mean(0), (w[:, None] * P).sum(0) / w.sum()
                k = rng.poisson(1.0, (NREP, len(w)))
                kw = k * w[None, :]
                rw = (kw @ P) / kw.sum(1)[:, None]
                r0 = (k @ P) / k.sum(1)[:, None]
                Rb, Ra = e_data / e0, e_data / ew
                cov_a = np.cov((rd / rw[rng.permutation(NREP)]).T)
                cov_b = np.cov((rd / r0[rng.permutation(NREP)]).T)
                c2[tn][0] += float((Rb - 1) @ np.linalg.solve(cov_b, Rb - 1))
                c2[tn][1] += float((Ra - 1) @ np.linalg.solve(cov_a, Ra - 1))
                B[tn] = {"n_sr": int(s.sum()), "neff_over_n": float(w.sum() ** 2 / (w ** 2).sum() / len(w)), "eps_mc_0": e0.tolist(), "eps_mc_w": ew.tolist(),
                         "R_before": Rb.tolist(), "R_after": Ra.tolist(), "err_after": np.sqrt(np.diag(cov_a)).tolist(), "err_before": np.sqrt(np.diag(cov_b)).tolist(),
                         "stat_mc_after_rel": (rw.std(0) / ew).tolist()}
            R["bins"][b] = B
        R["chi2"] = {tn: {"before": c2[tn][0], "after": c2[tn][1], "ndf": 15, "p_before": float(stats.chi2.sf(c2[tn][0], 15)), "p_after": float(stats.chi2.sf(c2[tn][1], 15))} for tn in c2}
        res[f"{run}|{fn}"] = R
        print(run, fn, "AUC", round(auc, 4), "neff_sr", round(R["neff_sr"], 3), R["chi2"], flush=True)
        for b in C.BINS:
            print("   ", b, "R_before", np.round(R["bins"][b]["t2"]["R_before"], 3), "R_after", np.round(R["bins"][b]["t2"]["R_after"], 3), "±", np.round(R["bins"][b]["t2"]["err_after"], 3),
                  "| SB mc_w/data", np.round(np.array(R["bins"][b]["sb_mc_w"]) / np.array(R["bins"][b]["sb_data"]), 3), flush=True)
    np.savez_compressed(os.path.join(C.OUT, f"s3_w_{run}.npz"), **wsave)
json.dump(res, open(os.path.join(C.OUT, "s3.json"), "w"), indent=1, default=float)
