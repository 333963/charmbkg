"""FEAT_charmnet: XGB H (18 cech, hiperparametry H z charmbkg default.yaml) vs H + grupy cech CharmNet.

Test mc  : sygnał fd truth 0/1 (D0: hypo_swap≠1; D+: weto D*), tło trening = mb truth≥8 (cały zakres mas),
           ε_B = mb truth≥8 w oknie |m−m0| < 55/45 MeV (OOF).
Test data: tło trening = dalekie SB danych (D0 >100 MeV + weto refl. |m_swap−m0|>25 MeV, D+ >110 MeV + weto D*),
           ε_B = te same SB danych (OOF), ε_S = fd (OOF). Dodatkowo wyniki na mb truth≥8 (rzeźbienie, ε_B MC).
Foldy: splitmix64(sample_id, evkey) (charmbkg.samples.fold_of), 5 foldów. Bootstrap Poissona po zdarzeniach (te same repliki
dla wszystkich wariantów → błąd ΔZ sparowany).
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/k/Projects/charm/charm-bkg")
from charmbkg import config as C, features as F  # noqa: E402
from charmbkg.samples import fold_of, reweight  # noqa: E402
import xgboost as xgb  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
cfg = C.load()
H = F.feature_list(cfg, "H")
SID = {"fd_nd": 1, "mb_nd": 2, "data_428770": 3}
GROUPS = {
    "D0": {
        "iso2": ["ptcone40_cn", "ptcone100_cn"],
        "iso4": ["ptcone40_cn", "ptcone60_cn", "ptcone80_cn", "ptcone100_cn"],
        "absimp": ["absImpact_cn"],
        "trkz0": ["d0AbsZ0_cn", "d1AbsZ0_cn"],
        "ptratio": ["dMinPtRatio_cn", "dMaxPtRatio_cn"],
        "apalpha": ["apAlpha_cn"],
    },
    "Dplus": {
        "iso2": ["ptcone40_cn", "ptcone100_cn"],
        "iso4": ["ptcone40_cn", "ptcone60_cn", "ptcone80_cn", "ptcone100_cn"],
        "absimp": ["absImpact_cn"],
        "maxdR": ["maxdR_cn"],
        "ptratio": ["dMinPt_cn", "dMinPtRatio_cn", "dMedPtRatio_cn", "dMaxPtRatio_cn"],
        "dalitz": ["mKpi1_over_mD_cn", "mKpi2_over_mD_cn"],
        "pid": [f"d{i}_dEdx{q}_cn" for i in range(3) for q in ("DistPerp", "Valid")],
    },
}
for sp, g in GROUPS.items():  # ALL = wszystkie cechy CharmNet „true” nieobecne w H; ALL_noiso = bez izolacji
    allf = [f for k, v in g.items() if k != "iso4" for f in v]
    g["ALL"] = list(dict.fromkeys(allf))
    g["ALL_noiso"] = [f for f in g["ALL"] if "ptcone" not in f]


def variants(sp):
    v = {"H": H, "H_seed2": H}
    for k, f in GROUPS[sp].items():
        v["H+" + k] = H + f
    return v


def load(sp):
    sp_cfg = cfg["species"][sp]
    fr = {s: pd.read_pickle(os.path.join(BASE, "frames", f"{s}_{sp}.pkl")) for s in SID}
    for s, L in fr.items():
        L["sid"] = SID[s]
        L["fold"] = fold_of(np.full(len(L), SID[s]), L["evkey"].to_numpy(np.int64), 5)
        L["dm0"] = np.abs(L["m"].to_numpy() - sp_cfg["m0"])
        ok = np.ones(len(L), bool)
        if sp == "Dplus":
            ok &= L["dstar_ok"].to_numpy(bool)
        L["ok"] = ok
    fd, mb, da = fr["fd_nd"], fr["mb_nd"], fr["data_428770"]
    t = fd["truth"].to_numpy(float)
    sig = np.isin(t, [0, 1]) & fd["ok"].to_numpy()
    if sp == "D0":
        sig &= fd["hypo_swap"].to_numpy() != 1
    S = fd[sig].reset_index(drop=True)
    tb = mb["truth"].to_numpy(float)
    Bmc = mb[(tb >= 8) & (tb < 100) & mb["ok"].to_numpy()].reset_index(drop=True)
    far = {"D0": 0.100, "Dplus": 0.110}[sp]
    dsel = da["ok"].to_numpy() & (da["dm0"].to_numpy() > far)
    if sp == "D0":
        dsel &= np.abs(da["m_swap"].to_numpy() - sp_cfg["m0"]) > 0.025
    Bda = da[dsel].reset_index(drop=True)
    # (dla rzeźbienia i modelowania) pełne dane po wetach, mb wszystkie kand. po wetach
    dall = da[da["ok"].to_numpy() & ((np.abs(da["m_swap"].to_numpy() - sp_cfg["m0"]) > 0.025) if sp == "D0" else True)].reset_index(drop=True)
    return S, Bmc, Bda, dall, mb


def fit_predict(S, Btr, feats, seed, cap, apply_sets, tag, log):
    """5-fold: trening na foldach ≠k (tło podpróbkowane do `cap`), predykcja OOF dla S, Btr i każdego zbioru w apply_sets
    (zbiory spoza treningu dostają wynik modelu swojego foldu — ten sam schemat co charmbkg.Ensemble)."""
    p = dict(cfg["training"]["xgb"])
    es = p.pop("early_stopping_rounds")
    p.pop("n_jobs")
    rng = np.random.default_rng(1000 + seed)
    out = {"S": np.full(len(S), np.nan, np.float32), "B": np.full(len(Btr), np.nan, np.float32)}
    for k in apply_sets:
        out[k] = np.full(len(apply_sets[k]), np.nan, np.float32)
    wb_all = reweight(Btr, S, cfg["training"]["reweight_bins"])
    XS, XB = S[feats].to_numpy(np.float32), Btr[feats].to_numpy(np.float32)
    best = []
    for k in range(5):
        t0 = time.time()
        ts, tb = S["fold"].to_numpy() != k, Btr["fold"].to_numpy() != k
        ib = np.flatnonzero(tb)
        if cap and len(ib) > cap:
            ib = np.sort(rng.choice(ib, cap, replace=False))
        is_ = np.flatnonzero(ts)
        X = np.vstack([XS[is_], XB[ib]])
        y = np.r_[np.ones(len(is_)), np.zeros(len(ib))]
        ws = S["weight"].to_numpy(float)[is_]
        wb = wb_all[ib]
        w = np.r_[ws / ws.sum(), wb / wb.sum()] * len(y) / 2
        va = rng.random(len(y)) < 0.1
        m = xgb.XGBClassifier(objective="binary:logistic", eval_metric="logloss", early_stopping_rounds=es,
                              n_jobs=8, random_state=seed, **p)
        m.fit(X[~va], y[~va], sample_weight=w[~va], eval_set=[(X[va], y[va])], sample_weight_eval_set=[w[va]], verbose=False)
        best.append(int(m.best_iteration))
        sel = S["fold"].to_numpy() == k
        out["S"][sel] = m.predict_proba(XS[sel])[:, 1]
        sel = Btr["fold"].to_numpy() == k
        out["B"][sel] = m.predict_proba(XB[sel])[:, 1]
        for nm, A in apply_sets.items():
            sel = A["fold"].to_numpy() == k
            if sel.any():
                out[nm][sel] = m.predict_proba(A.loc[sel, feats].to_numpy(np.float32))[:, 1]
        if k == 0:
            g = m.get_booster().get_score(importance_type="gain")
            tot = sum(g.values()) or 1
            imp = {feats[int(f[1:])] if f.startswith("f") and f[1:].isdigit() else f: round(v / tot, 4) for f, v in g.items()}
        log(f"  {tag} fold {k}: n_s={len(is_)} n_b={len(ib)} best_it={best[-1]} {time.time()-t0:.0f}s")
    return out, best, dict(sorted(imp.items(), key=lambda x: -x[1]))


BINS = [(3, 5), (5, 8), (8, 30)]
WPS = [0.80, 0.70, 0.50]


def zproxy(sS, wS, ptS, sB, wB, ptB):
    """dict bin→wp→(eps_S, eps_B, Z); progi z sygnału per bin (ważone), ε_B ważone."""
    r = {}
    for lo, hi in BINS:
        a = (ptS >= lo) & (ptS < hi)
        b = (ptB >= lo) & (ptB < hi)
        s, w = sS[a], wS[a]
        o = np.argsort(-s)
        c = np.cumsum(w[o]) / w.sum()
        sb, wb = sB[b], wB[b]
        for wp in WPS:
            thr = s[o][min(np.searchsorted(c, wp), len(o) - 1)]
            es = w[s >= thr].sum() / w.sum()
            eb = wb[sb >= thr].sum() / max(wb.sum(), 1e-12)
            r[(lo, hi, wp)] = (es, eb, es / np.sqrt(max(eb, 1e-12)))
    return r


def boot_weights(keys, nb, seed):
    """Poisson(1) na zdarzenie (klucz), te same repliki dla wszystkich wariantów."""
    u, inv = np.unique(keys, return_inverse=True)
    rng = np.random.default_rng(seed)
    return rng.poisson(1.0, (nb, len(u))).astype(np.float32)[:, inv]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--test", choices=["mc", "data"], required=True)
    ap.add_argument("--variants", default="all")
    ap.add_argument("--cap", type=int, default=600000)
    ap.add_argument("--nboot", type=int, default=200)
    a = ap.parse_args()
    sp = a.species
    m0 = cfg["species"][sp]["m0"]
    win = cfg["species"][sp]["signal_window"]
    odir = os.path.join(BASE, "out", f"{sp}_{a.test}")
    os.makedirs(odir, exist_ok=True)
    lf = open(os.path.join(odir, "run.log"), "a")

    def log(s):
        print(s, flush=True)
        lf.write(s + "\n")
        lf.flush()
    S, Bmc, Bda, dall, mb = load(sp)
    log(f"{sp} {a.test}: S={len(S)} Bmc={len(Bmc)} (SR {int((Bmc.dm0 < win).sum())}) Bdata_farSB={len(Bda)} data_all={len(dall)}")
    V = variants(sp)
    if a.variants != "all":
        V = {k: V[k] for k in a.variants.split(",")}
    Btr = Bmc if a.test == "mc" else Bda
    apply_sets = {"mb": Bmc, "dall": dall} if a.test == "data" else {"dall": dall}
    scores = {}
    meta = {}
    for name, feats in V.items():
        f = os.path.join(odir, f"scores_{name}.npz")
        if os.path.exists(f):
            z = np.load(f)
            scores[name] = {k: z[k] for k in z.files}
            meta[name] = json.load(open(f.replace(".npz", ".json")))
            log(f"{name}: z cache")
            continue
        t0 = time.time()
        out, best, imp = fit_predict(S, Btr, feats, 2 if name == "H_seed2" else 1, a.cap, apply_sets, name, log)
        np.savez_compressed(f, **out)
        meta[name] = {"features": feats, "best_it": best, "importance_fold0": imp, "sec": round(time.time() - t0)}
        json.dump(meta[name], open(f.replace(".npz", ".json"), "w"), indent=1)
        scores[name] = out
        log(f"{name}: {time.time()-t0:.0f}s best_it={best}")
    # --- Zproxy + bootstrap
    if a.test == "mc":
        bsel = Bmc["dm0"].to_numpy() < win
        evalB = Bmc[bsel]
        sB = {n: scores[n]["B"][bsel] for n in scores}
    else:
        evalB = Bda
        sB = {n: scores[n]["B"] for n in scores}
    wS, ptS = S["weight"].to_numpy(float), S["pt"].to_numpy(float)
    wB, ptB = evalB["weight"].to_numpy(float), evalB["pt"].to_numpy(float)
    bwS = boot_weights(S["evkey"].to_numpy(np.int64), a.nboot, 11)
    bwB = boot_weights(evalB["evkey"].to_numpy(np.int64) + 10**9 * evalB["sid"].to_numpy(np.int64), a.nboot, 12)
    nom = {n: zproxy(scores[n]["S"], wS, ptS, sB[n], wB, ptB) for n in scores}
    reps = {n: [zproxy(scores[n]["S"], wS * bwS[r], ptS, sB[n], wB * bwB[r], ptB) for r in range(a.nboot)] for n in scores}
    rows = []
    for n in scores:
        for key, (es, eb, z) in nom[n].items():
            zr = np.array([reps[n][r][key][2] for r in range(a.nboot)])
            dz = np.array([reps[n][r][key][2] / reps["H"][r][key][2] for r in range(a.nboot)]) if "H" in reps else zr * np.nan
            rows.append(dict(variant=n, bin=f"{key[0]}-{key[1]}", wp=key[2], eps_S=es, eps_B=eb, Z=z, Z_err=float(np.std(zr)),
                             ratio_vs_H=z / nom["H"][key][2] if "H" in nom else np.nan, ratio_err=float(np.std(dz)),
                             nB_pass=float(wB[(ptB >= key[0]) & (ptB < key[1]) & (sB[n] >= 0)].sum() * eb)))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(odir, "zproxy.csv"), index=False)
    log(T[T.wp == 0.7].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    json.dump(meta, open(os.path.join(odir, "meta.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
