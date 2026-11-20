"""D23/ANALIZA krok 5: sufit dotrenowania klasyfikatora w 2–3 GeV (szybki XGB, 18 cech H, CV 5-krotne po zdarzeniach).
Wariant A: sygnał fd_nd vs pasma boczne danych 427929 (może uczyć się różnic dane/MC -> górne ograniczenie).
Wariant B: sygnał fd_nd vs tło fd_nd (truth >= 8) -> bez różnic dane/MC, oceniany na pasmach danych.
Odniesienie: istniejący score H (trening pT >= 3). Miara: odrzucanie tła (pasma danych, ważone w_ps) przy stałym ε_S i oczekiwane Z.
Okno sygnału danych nieużywane. Zapis: out/retrain.json.  Użycie: python s05_retrain.py"""
import json
import os
import time

import numpy as np
import pandas as pd
import xgboost as xgb

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
IN = os.path.join(os.path.dirname(TOP), "runs", "lxplus")
M0 = 1.86966
FEATS = ['log_pt', 'abs_eta', 'chi2', 'Lxy', 'impact_sig', 'abs_impact_z0', 'd_max_abs_d0', 'pointing', 'abs_cos_theta_star', 'abs_impact_sig', 'd_min_abs_d0',
         'd_max_abs_z0', 'd_max_abs_eta', 'dedx_K', 'log_p_K', 'npts_K', 'dedx_pi_max', 'log_p_pi_min']
PAR = dict(objective="binary:logistic", eta=0.1, max_depth=6, subsample=0.8, colsample_bytree=0.8, min_child_weight=5, tree_method="hist", nthread=8, eval_metric="logloss")
NR = 400
EPS = [0.5, 0.4, 0.3, 0.2, 0.1]


def wq(s, w, x):
    o = np.argsort(-s)
    c = np.cumsum(w[o]) / w.sum()
    return float(s[o][min(np.searchsorted(c, x), len(o) - 1)])


def main():
    t0 = time.time()
    rng = np.random.default_rng(5)
    out = {}
    mc = pd.read_parquet(os.path.join(IN, "dp_fd_nd.parquet"))
    mc = mc[mc.dstar_veto == 0]
    sig = mc[mc.is_sig == 1].copy()
    mcb = mc[mc.truth >= 8].copy()
    sig["g"] = (sig.file_idx.to_numpy(np.int64) * 10_000_000 + sig.entry.to_numpy(np.int64)) % 5
    mcb["g"] = (mcb.file_idx.to_numpy(np.int64) * 10_000_000 + mcb.entry.to_numpy(np.int64)) % 5
    d = pd.read_parquet(os.path.join(IN, "dp_data_427929.parquet"), columns=FEATS + ["run", "event", "m", "score", "dstar_veto", "w_ps", "ps_tier", "grl", "trig_HLT_mb_mbts_L1MBTS_1"])
    d = d[(d.dstar_veto == 0) & (d.grl == 1) & (d["trig_HLT_mb_mbts_L1MBTS_1"] == 1) & (np.abs(d.m - M0) > 0.060)]
    # przerzedzenie tier 0 dla szybkości (waga skorygowana)
    t0m = (d.ps_tier == 0).to_numpy()
    keep = ~t0m | (rng.random(len(d)) < 0.6)
    d = d[keep].copy()
    d["w"] = d.w_ps.to_numpy(float) * np.where(d.ps_tier == 0, 1 / 0.6, 1.0)
    d["g"] = (d.event.to_numpy(np.int64) * 2654435761 % 1000003) % 5
    out["n"] = dict(sig=len(sig), mc_bkg_rows=len(mcb), mc_bkg_sumw=float(mcb.w_ps.sum()), data_side_rows=len(d), data_side_sumw=float(d.w.sum()))
    print(out["n"], flush=True)
    Xs, Xd, Xb = sig[FEATS].to_numpy(np.float32), d[FEATS].to_numpy(np.float32), mcb[FEATS].to_numpy(np.float32)
    wd, wb = d.w.to_numpy(float), mcb.w_ps.to_numpy(float)
    oofS = {"A": np.zeros(len(sig)), "B": np.zeros(len(sig))}
    oofD = {"A": np.zeros(len(d)), "B": np.zeros(len(d))}
    for k in range(5):
        trs, tes = (sig.g != k).to_numpy(), (sig.g == k).to_numpy()
        trd, ted = (d.g != k).to_numpy(), (d.g == k).to_numpy()
        trb = (mcb.g != k).to_numpy()
        for var, Xbk, wbk in (("A", Xd[trd], wd[trd]), ("B", Xb[trb], wb[trb])):
            X = np.vstack([Xs[trs], Xbk])
            y = np.r_[np.ones(trs.sum()), np.zeros(len(Xbk))]
            w = np.r_[np.full(trs.sum(), wbk.sum() / trs.sum()), wbk]
            w = w / w.mean()
            bst = xgb.train(PAR, xgb.DMatrix(X, y, weight=w), NR)
            oofS[var][tes] = bst.predict(xgb.DMatrix(Xs[tes]))
            oofD[var][ted] = bst.predict(xgb.DMatrix(Xd[ted]))
        print(f"fold {k} ({time.time()-t0:.0f} s)", flush=True)
    ws = np.ones(len(sig))
    scS = {"H": sig.score.to_numpy(float), **oofS}
    scD = {"H": d.score.to_numpy(float), **oofD}
    md = d.m.to_numpy(float)
    lo, hi = md < M0, md > M0
    # normalizacja: B(±45) = k * N_pasm; k z T1 (WP50, 427929); S = 0,6 * 4331 * 0,93 * eps
    WP = json.load(open(os.path.join(TOP, "out", "wp.json")))
    r50 = [r for r in WP["wp"]["427929"] if r["wp"] == 0.5][0]
    kB = r50["B"] / r50["n_side"]
    out["kB"] = kB
    res = {}
    for name in ("H", "A", "B"):
        res[name] = {}
        for e in EPS:
            t = wq(scS[name], ws, e)
            ps = scD[name] >= t
            nb = float(wd[ps].sum())
            nrow = int(ps.sum())
            rel = float(np.sqrt((wd[ps] ** 2).sum()) / max(nb, 1e-9))
            S = 0.6 * 4331 * 0.93 * e
            B = kB * nb
            res[name][str(e)] = dict(thr=t, eff_b=nb / wd.sum(), n_side=nb, n_rows=nrow, rel_err=rel, B=B, S=S, Z=S / np.sqrt(S + B), SB=S / B,
                                     lo_over_hi=float(wd[ps & lo].sum() / max(wd[ps & hi].sum(), 1e-9)))
        print(name, {e: f"eff_b {v['eff_b']:.2e} ±{100*v['rel_err']:.0f}% Z {v['Z']:.1f} S/B {v['SB']:.3f} lo/hi {v['lo_over_hi']:.2f}" for e, v in res[name].items()})
    out["res"] = res
    out["lo_over_hi_all"] = float(wd[lo].sum() / wd[hi].sum())
    # AUC ważone (sygnał vs pasma danych)
    for name in ("H", "A", "B"):
        s = np.r_[scS[name], scD[name]]
        y = np.r_[np.ones(len(sig)), np.zeros(len(d))]
        w = np.r_[ws, wd]
        o = np.argsort(s)
        cw = np.cumsum(w[o] * (1 - y[o]))
        auc = float((w[o] * y[o] * cw).sum() / (ws.sum() * wd.sum()))
        out.setdefault("auc", {})[name] = auc
    print("AUC", out["auc"], "lo/hi wszystkie", out["lo_over_hi_all"])
    json.dump(out, open(os.path.join(TOP, "out", "retrain.json"), "w"), indent=1)
    print(f"gotowe {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
