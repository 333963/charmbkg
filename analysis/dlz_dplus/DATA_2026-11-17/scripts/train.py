"""DLZ_dplus/DATA: trening H, H+DLZ, H+dalitz wg PROCEDURA §2 (odstępstwa: KRYTERIA.md §4) z zapisem modeli.

Sygnał = fd truth 0/1 po weto D*; tło = SB danych 428770 |m−m0| > 70 MeV po weto D* (pasmo 45–70 wyłączone);
przeważenie tła do (pT,|η|) sygnału; 5 foldów po zdarzeniach (MC: splitmix64(sid, evkey); dane: splitmix64(run, event));
hiperparametry XGB H; cap tła 300k/fold. Wyjście: models/<wariant>_s<seed>/fold{k}.json + meta.json, oof/<wariant>_s<seed>.npz
(S, mb, dloc = lokalna ramka danych po weto), oof/base.npz (m, pt, waga, truth, klucze).
Użycie: python train.py [wariant:seed ...]
"""
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

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FR = "/home/k/Documents/charm/charm_work/FEAT_charmnet/frames"
NJ = int(os.environ.get("DLZ_NJ", "3"))
CAP = 300000


def load():
    fr = {}
    for s in ("fd_nd", "mb_nd", "data_428770"):
        L = pd.read_pickle(os.path.join(FR, f"{s}_Dplus.pkl"))
        L = L[L["dstar_ok"].to_numpy(bool)].reset_index(drop=True)
        if s.startswith("data"):
            L["fold"] = D.data_fold(L)
            L["key"] = L["event"].to_numpy().astype(np.int64)
        else:
            L["fold"] = D.fold_of(np.full(len(L), D.SID[s]), L["evkey"].to_numpy(np.int64), 5)
            L["key"] = L["evkey"].to_numpy(np.int64)
        L["dm0"] = np.abs(L["m"].to_numpy() - D.M0)
        fr[s] = L
    return fr


def main():
    jobs = [a.split(":") for a in sys.argv[1:]] or [["H", "12345"], ["HDLZ", "12345"], ["Hdal", "12345"], ["H", "23456"], ["HDLZ", "23456"]]
    fr = load()
    fd, mb, da = fr["fd_nd"], fr["mb_nd"], fr["data_428770"]
    sig = np.isin(fd["truth"].to_numpy(), [0, 1])
    S = fd[sig].reset_index(drop=True)
    bsel = da["dm0"].to_numpy() > D.SP["sideband_min"]
    B = da[bsel].reset_index(drop=True)
    # G1: żaden klucz zdarzenia w > 1 foldzie
    viol = 0
    for L in (S, B, mb, da):
        g = pd.DataFrame({"k": L["key"], "f": L["fold"]}).groupby("k")["f"].nunique()
        viol += int((g > 1).sum())
    print(f"S={len(S)} B(SB>70MeV)={len(B)} mb={len(mb)} data={len(da)}; naruszenia foldów: {viol}", flush=True)
    os.makedirs(os.path.join(BASE, "oof"), exist_ok=True)
    np.savez_compressed(os.path.join(BASE, "oof", "base.npz"),
                        S_m=S["m"], S_pt=S["pt"], S_w=S["weight"], S_key=S["key"], S_fold=S["fold"],
                        mb_m=mb["m"], mb_pt=mb["pt"], mb_truth=mb["truth"], mb_key=mb["key"], mb_w=mb["weight"],
                        d_m=da["m"], d_pt=da["pt"], d_key=da["key"], d_fold=da["fold"],
                        fd_m=fd["m"], fd_pt=fd["pt"], fd_truth=fd["truth"], fd_key=fd["key"], fd_w=fd["weight"], fold_viol=viol)
    p = dict(D.cfg["training"]["xgb"])
    es = p.pop("early_stopping_rounds")
    p.pop("n_jobs")
    wb_all = reweight(B, S, D.cfg["training"]["reweight_bins"])
    for name, seed in jobs:
        seed = int(seed)
        feats = D.VARIANTS[name]
        tag = f"{name}_s{seed}"
        mdir = os.path.join(BASE, "models", tag)
        of = os.path.join(BASE, "oof", tag + ".npz")
        if os.path.exists(of):
            print(tag, "cache", flush=True)
            continue
        os.makedirs(mdir, exist_ok=True)
        rng = np.random.default_rng(seed)
        XS, XB = S[feats].to_numpy(np.float32), B[feats].to_numpy(np.float32)
        out = {"S": np.full(len(S), np.nan, np.float32), "mb": np.full(len(mb), np.nan, np.float32),
               "dloc": np.full(len(da), np.nan, np.float32), "fdall": np.full(len(fd), np.nan, np.float32)}
        best, auc_tr, auc_oof, imp = [], [], [], {}
        for k in range(5):
            t0 = time.time()
            is_ = np.flatnonzero(S["fold"].to_numpy() != k)
            ib = np.flatnonzero(B["fold"].to_numpy() != k)
            if len(ib) > CAP:
                ib = np.sort(rng.choice(ib, CAP, replace=False))
            X = np.vstack([XS[is_], XB[ib]])
            y = np.r_[np.ones(len(is_)), np.zeros(len(ib))]
            ws, wb = S["weight"].to_numpy(float)[is_], wb_all[ib]
            w = np.r_[ws / ws.sum(), wb / wb.sum()] * len(y) / 2
            va = rng.random(len(y)) < 0.1
            m = xgb.XGBClassifier(objective="binary:logistic", eval_metric="logloss", early_stopping_rounds=es, n_jobs=NJ, random_state=seed, **p)
            m.fit(X[~va], y[~va], sample_weight=w[~va], eval_set=[(X[va], y[va])], sample_weight_eval_set=[w[va]], verbose=False)
            bo = m.get_booster()
            bo.save_model(os.path.join(mdir, f"fold{k}.json"))
            best.append(int(m.best_iteration))
            ptr = m.predict_proba(X[~va])[:, 1]
            auc_tr.append(float(roc_auc_score(y[~va], ptr, sample_weight=w[~va])))
            so = S["fold"].to_numpy() == k
            bo_ = np.flatnonzero(B["fold"].to_numpy() == k)
            if len(bo_) > 100000:
                bo_ = bo_[:100000]
            ps, pb = m.predict_proba(XS[so])[:, 1], m.predict_proba(XB[bo_])[:, 1]
            wso, wbo = S["weight"].to_numpy(float)[so], wb_all[bo_]
            auc_oof.append(float(roc_auc_score(np.r_[np.ones(len(ps)), np.zeros(len(pb))], np.r_[ps, pb],
                                               sample_weight=np.r_[wso / wso.sum(), wbo / wbo.sum()])))
            out["S"][so] = ps
            for nm, A in (("mb", mb), ("dloc", da), ("fdall", fd)):
                sel = A["fold"].to_numpy() == k
                out[nm][sel] = m.predict_proba(A.loc[sel, feats].to_numpy(np.float32))[:, 1]
            if k == 0:
                g = bo.get_score(importance_type="gain")
                tot = sum(g.values()) or 1
                imp = {feats[int(f[1:])]: round(v / tot, 4) for f, v in sorted(g.items(), key=lambda x: -x[1])}
            print(f"  {tag} fold {k}: n_s={len(is_)} n_b={len(ib)} best_it={best[-1]} AUC tr={auc_tr[-1]:.4f} oof={auc_oof[-1]:.4f} {time.time()-t0:.0f}s", flush=True)
        np.savez_compressed(of, **out)
        json.dump({"features": feats, "seed": seed, "best_it": best, "auc_train": auc_tr, "auc_oof": auc_oof, "importance_fold0": imp,
                   "n_S": len(S), "n_B_sb70": len(B), "cap": CAP, "fold_violations": viol},
                  open(os.path.join(mdir, "meta.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
