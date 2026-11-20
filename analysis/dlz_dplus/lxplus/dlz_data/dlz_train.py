#!/usr/bin/env python3
"""DLZ_DATA: trening H, H+DLZ, H+dalitz wg DLZ_dplus/PROCEDURA.md §2 na pełnych próbkach (h5 v2.0, L_xy ≥ 0,25 mm).
Sygnał fd truth 0/1 po weto D*; tło = SB danych 428770 |m−m0| > 70 MeV; przeważenie tła do (pT,|η|) sygnału; 5 foldów po zdarzeniach
(MC: splitmix64(sid, evkey), dane: splitmix64(run, event)); hiperparametry XGB H; max_bkg_train = DLZ_CAP (3·10⁶), max_sig 2·10⁶;
early stopping na 10 % (jak FEAT/DLZ). Wyjście (układ = lokalny DLZ_dplus/DATA_*): models/<tag>/fold{k}.json + meta.json, thr.json,
oof/base.npz, oof/<tag>.npz (S, mb), oof/S_feat.npz.
Użycie: dlz_train.py --parts DIR --out DIR   (parts: frm_fd_*.npz, frm_mb_*.npz, frm_sb_*.npz)
Env: DLZ_SETS="H:12345 HDLZ:12345 Hdal:12345 H:23456 HDLZ:23456", DLZ_CAP, DLZ_NJ, DLZ_FAST=1 (test: 40 drzew)."""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dlzlib as D  # noqa: E402
from charmbkg.samples import reweight  # noqa: E402

FEATS = D.H + D.DLZ8


def cat(parts, tag):
    fs = sorted(glob.glob(os.path.join(parts, f"frm_{tag}_*.npz")))
    if not fs:
        raise SystemExit(f"brak frm_{tag}_*.npz w {parts}")
    zs = [np.load(f) for f in fs]
    return {k: np.concatenate([z[k] for z in zs]) for k in ("X", "m", "pt", "abs_eta", "weight", "truth", "run", "event", "evkey")}, len(fs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sets = [x.split(":") for x in os.environ.get("DLZ_SETS", "H:12345 HDLZ:12345 Hdal:12345 H:23456 HDLZ:23456").split()]
    cap = int(os.environ.get("DLZ_CAP", "3000000"))
    nj = int(os.environ.get("DLZ_NJ", "8"))
    fast = os.environ.get("DLZ_FAST", "0") == "1"
    S, nfd = cat(a.parts, "fd")
    B, nsb = cat(a.parts, "sb")
    MB, nmb = cat(a.parts, "mb")
    S["fold"] = D.fold_of(np.full(len(S["m"]), 1), S["evkey"], 5)
    MB["fold"] = D.fold_of(np.full(len(MB["m"]), 2), MB["evkey"], 5)
    B["fold"] = D.fold_of(B["run"], B["event"], 5)
    viol = 0
    for Z, key in ((S, S["evkey"]), (MB, MB["evkey"]), (B, B["event"])):
        g = pd.DataFrame({"k": key, "r": Z["run"], "f": Z["fold"]}).groupby(["r", "k"])["f"].nunique()
        viol += int((g > 1).sum())
    print(f"S={len(S['m'])} ({nfd} cz.) B_sb70={len(B['m'])} ({nsb} cz.) mb={len(MB['m'])} ({nmb} cz.); naruszenia foldów: {viol}; cap={cap} nj={nj}", flush=True)
    if viol:
        raise SystemExit("G1 FAIL: zdarzenie w > 1 foldzie")
    if len(S["m"]) > 2_000_000:
        raise SystemExit("max_sig_train przekroczone — nieoczekiwane")
    os.makedirs(os.path.join(a.out, "oof"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "models"), exist_ok=True)
    np.savez_compressed(os.path.join(a.out, "oof", "base.npz"), S_m=S["m"], S_pt=S["pt"], S_w=S["weight"], S_key=S["evkey"], S_fold=S["fold"],
                        mb_m=MB["m"], mb_pt=MB["pt"], mb_truth=MB["truth"], mb_key=MB["evkey"], mb_w=MB["weight"], fold_viol=viol)
    np.savez_compressed(os.path.join(a.out, "oof", "S_feat.npz"), feat=S["X"], names=np.array(FEATS))
    df = lambda Z: pd.DataFrame({"pt": Z["pt"], "abs_eta": Z["abs_eta"], "weight": Z["weight"].astype(float)})  # noqa: E731
    dS, dB = df(S), df(B)
    dB["weight"] = 1.0
    wb_all = reweight(dB, dS, D.cfg["training"]["reweight_bins"])
    p = dict(D.cfg["training"]["xgb"])
    es = p.pop("early_stopping_rounds")
    p.pop("n_jobs")
    if fast:
        p["n_estimators"] = 40
    T = {}
    for name, seed in sets:
        seed = int(seed)
        cols = [FEATS.index(f) for f in D.VARIANTS[name]]
        tag = f"{name}_s{seed}"
        mdir = os.path.join(a.out, "models", tag)
        os.makedirs(mdir, exist_ok=True)
        rng = np.random.default_rng(seed)
        XS, XB, XM = S["X"][:, cols], B["X"][:, cols], MB["X"][:, cols]
        oS, oM = np.full(len(XS), np.nan, np.float32), np.full(len(XM), np.nan, np.float32)
        best, atr, aoof, imp = [], [], [], {}
        ws_all = S["weight"].astype(float)
        for k in range(5):
            t0 = time.time()
            is_ = np.flatnonzero(S["fold"] != k)
            ib = np.flatnonzero(B["fold"] != k)
            if len(ib) > cap:
                ib = np.sort(rng.choice(ib, cap, replace=False))
            X = np.vstack([XS[is_], XB[ib]])
            y = np.r_[np.ones(len(is_)), np.zeros(len(ib))]
            ws, wb = ws_all[is_], wb_all[ib]
            w = np.r_[ws / ws.sum(), wb / wb.sum()] * len(y) / 2
            va = rng.random(len(y)) < 0.1
            m = xgb.XGBClassifier(objective="binary:logistic", eval_metric="logloss", early_stopping_rounds=es, n_jobs=nj, random_state=seed, **p)
            m.fit(X[~va], y[~va], sample_weight=w[~va], eval_set=[(X[va], y[va])], sample_weight_eval_set=[w[va]], verbose=False)
            bo = m.get_booster()
            bo.save_model(os.path.join(mdir, f"fold{k}.json"))
            best.append(int(m.best_iteration))
            sub = rng.random((~va).sum()) < min(1.0, 300000 / max((~va).sum(), 1))
            atr.append(float(roc_auc_score(y[~va][sub], m.predict_proba(X[~va][sub])[:, 1], sample_weight=w[~va][sub])))
            so = S["fold"] == k
            bo_ = np.flatnonzero(B["fold"] == k)[:200000]
            ps, pb = m.predict_proba(XS[so])[:, 1], m.predict_proba(XB[bo_])[:, 1]
            aoof.append(float(roc_auc_score(np.r_[np.ones(len(ps)), np.zeros(len(pb))], np.r_[ps, pb],
                                            sample_weight=np.r_[ws_all[so] / ws_all[so].sum(), wb_all[bo_] / wb_all[bo_].sum()])))
            oS[so] = ps
            sm = MB["fold"] == k
            if sm.any():
                oM[sm] = m.predict_proba(XM[sm])[:, 1]
            if k == 0:
                g = bo.get_score(importance_type="gain")
                tot = sum(g.values()) or 1
                imp = {D.VARIANTS[name][int(f[1:])]: round(v / tot, 4) for f, v in sorted(g.items(), key=lambda x: -x[1])}
            print(f"  {tag} fold {k}: n_s={len(is_)} n_b={len(ib)} best_it={best[-1]} AUC tr={atr[-1]:.4f} oof={aoof[-1]:.4f} {time.time()-t0:.0f}s", flush=True)
        np.savez_compressed(os.path.join(a.out, "oof", tag + ".npz"), S=oS, mb=oM)
        json.dump({"features": D.VARIANTS[name], "seed": seed, "best_it": best, "auc_train": atr, "auc_oof": aoof, "importance_fold0": imp,
                   "n_S": len(XS), "n_B_sb70": len(XB), "cap": cap, "fold_violations": viol, "fast": fast}, open(os.path.join(mdir, "meta.json"), "w"), indent=1)
        th = D.thresholds(oS, ws_all, S["pt"].astype(float))
        T[tag] = {f"{lo}-{hi}": {f"{wp:.2f}": th[(lo, hi, wp)] for wp in D.WPS} for lo, hi in D.BINS}
        json.dump(T, open(os.path.join(a.out, "thr.json"), "w"), indent=1)
    print("OK trening", list(T), flush=True)


if __name__ == "__main__":
    main()
