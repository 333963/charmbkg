#!/usr/bin/env python3
"""SALAD427 — studium: SALAD v1.0 na runie 427929 (MBTS) i w klasach trygera 428770; transfer wag; closure adwersarialny (B1).

Wejście: części s427_<próbka>_<g>.npz (s427_extract.py). Tylko pasma boczne i tło — w SR danych liczone są wyłącznie
interpolacje z pasm bocznych (exp(pol2) ML), żadnego fitu piku.
Zadania (--tasks main,adv):
  main: (a) tagger H427: HGB fd (truth 0/1) vs SB danych MB427929, K foldów po zdarzeniach, score OOF dla wszystkiego; progi WP
            z fd OOF per bin pT (procedura B1 `b1_dataeff_d0.py`, tu z k-fold zamiast parzystości i z jednym taggerem dla klas);
        (b) per klasa {MB427929, J12only, MBTSbased, AFPonly}: SALAD v1.0 (`fit_salad_v1`, hgb) K-fold po zdarzeniach (wagi OOF na MC,
            2 ziarna) + bieg „prod” (jeden fit na całym SB, bramki v1.0: n_eff/N(SR), N_SB, C2ST pasma);
            eps_B: SB (OOF), pasmo walidacyjne, SR (MC ważone) vs interpolacja SB danych; C2ST przed/po; n_eff/N; rzeźbienie masy;
        (c) transfer: wagi klasy a -> dane klasy b (wszystkie pary);
        (d) D0: wariant MC truth >= 8 dla MB427929 (nominał D0 = truth >= 2).
  adv:  B1 na MC v2.0: pseudo-dane = MC (parzystość 1) zniekształcone (smooth, cat_global, cat_local, cat_step, cuts_local),
        zmniejszone do statystyki SB klasy MB427929; „MC” = parzystość 0; tagger H427.
Wyjście: s427_<gatunek>.json, s427_<gatunek>_adv.json, s427_<gatunek>_w.npz (wagi OOF; zostaje na EOS).
Użycie: s427_study.py --species D0|Dplus --parts DIR --out DIR [--tasks main,adv] [--threads 8]   (S427_FAST=1: test)"""
import argparse
import glob
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s427_lib as SL  # noqa: E402

SL.add_src()
warnings.filterwarnings("ignore")
from charmbkg.bkggen import salad_v1 as SV  # noqa: E402
from charmbkg.bkggen.feta import c2st as _c2st  # noqa: E402

FAST = os.environ.get("S427_FAST", "0") == "1"
K = 2 if FAST else SL.K_FOLDS
NBOOT = 20 if FAST else 400
SEEDS = (0,) if FAST else (0, 1)
FIT_CAP = 1_500_000                       # wierszy na stronę przekazywanych do fit_salad_v1 w foldach (fit bierze z nich max_n)
MIN_SB = 200 if FAST else 2000            # poniżej tej liczby kandydatów SB klasa jest pomijana (zapis `skipped`)
PFX = {"D0": "D0", "Dplus": "Dp"}
T0 = time.time()


def say(*a):
    print(f"[{time.time() - T0:6.0f}s]", *a, flush=True)


# ----------------------------------------------------------------------------------------------------- wczytanie
def load(parts, sample, sp):
    fl = sorted(glob.glob(os.path.join(parts, f"s427_{sample}_*.npz")))
    if not fl:
        raise SystemExit(f"brak części s427_{sample}_*.npz w {parts}")
    acc = {}
    for f in fl:
        z = np.load(f, allow_pickle=False)
        for k in z.files:
            if k.startswith(PFX[sp] + "_"):
                acc.setdefault(k[len(PFX[sp]) + 1:], []).append(z[k])
    D = {k: np.concatenate(v) for k, v in acc.items()}
    D["n_parts"] = len(fl)
    return D


def sel(D, mask):
    return {k: (v[mask] if isinstance(v, np.ndarray) and len(v) == len(mask) else v) for k, v in D.items()}


def prep(D, spp, sp):
    """Weto refleksji D0, regiony, fold, bin pT."""
    if sp == "D0":
        D = sel(D, np.abs(D["mswap"].astype(float) - spp["m0"]) > SL.REFL_VETO)
    D["reg"] = SV.regions(D["m"], spp)
    D["fold"] = (SL.mix64(D["uid"], 7) % np.uint64(K)).astype(np.int8)
    D["ib"] = SL.ptbin(D["pt"])
    return sel(D, (D["reg"] != "out") & (D["ib"] >= 0) & (D["ib"] < 3))


def frame(D, feats, mask=None):
    X = D["X"] if mask is None else D["X"][mask]
    df = pd.DataFrame(X.astype(np.float64), columns=feats)
    df["m"] = D["m"] if mask is None else D["m"][mask]
    return df


# ----------------------------------------------------------------------------------------------------- tagger H427
def hgb(X, y, seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(max_iter=30 if FAST else 200, learning_rate=0.08, max_leaf_nodes=31, l2_regularization=1.0,
                                          early_stopping=True, validation_fraction=0.15, random_state=seed).fit(X, y)


def tagger(S, Dmb, others, rng):
    """K modeli fd(fold != k) vs SB danych MB427929 (fold != k); score OOF dla fd, danych i MC (fold wg hasha zdarzenia)."""
    from sklearn.metrics import roc_auc_score
    sig = np.isin(S["truth"], [0, 1])
    for D in [S, Dmb] + others:
        D["score"] = np.full(len(D["m"]), np.nan, np.float32)
    aucs = []
    for k in range(K):
        xs = S["X"][sig & (S["fold"] != k)]
        ib = np.flatnonzero((Dmb["reg"] == "SB") & (Dmb["fold"] != k))
        if len(ib) > 250_000:
            ib = rng.choice(ib, 250_000, replace=False)
        clf = hgb(np.vstack([xs, Dmb["X"][ib]]), np.r_[np.ones(len(xs)), np.zeros(len(ib))], seed=k)
        for D in [S, Dmb] + others:
            te = D["fold"] == k
            if te.any():
                D["score"][te] = clf.predict_proba(D["X"][te])[:, 1]
        a = (Dmb["reg"] == "SB") & (Dmb["fold"] == k)
        b = sig & (S["fold"] == k)
        if a.any() and b.any():
            aucs.append(float(roc_auc_score(np.r_[np.ones(b.sum()), np.zeros(a.sum())], np.r_[S["score"][b], Dmb["score"][a]])))
        say(f"tagger fold {k}: n_sig {len(xs)} n_sb {len(ib)} AUC_OOF {aucs[-1] if aucs else float('nan'):.4f}")
    thr = np.array([[np.quantile(S["score"][sig & (S["ib"] == i)], 1 - wp) if (sig & (S["ib"] == i)).any() else np.nan
                     for wp in SL.WPS] for i in range(3)])
    for D in [S, Dmb] + others:
        D["ncat"] = (D["score"][:, None] >= thr[D["ib"]]).sum(1).astype(np.int8)
    return thr, aucs


# ----------------------------------------------------------------------------------------------------- SALAD
def oof_weights(cfg, sp, D, MC, feats, seed):
    """Wagi OOF dla MC: model foldu k uczony na danych SB i MC z foldów != k, stosowany do MC z foldu k."""
    w = np.zeros(len(MC["m"]))
    aucs, nclip = [], 0
    rs = np.random.default_rng(1000 + seed)
    for k in range(K):
        # fit v1.0 i tak losuje max_n (400 tys.) na klasę; wstępne obcięcie do FIT_CAP wierszy tylko oszczędza pamięć i czas
        d = frame(D, feats, cap((D["fold"] != k) & (D["reg"] != "SR"), rs))
        rw = SV.fit_salad_v1(cfg, sp, d, frame(MC, feats, cap(MC["fold"] != k, rs)), kind="hgb", seed=seed, c2st_band=False)
        te = MC["fold"] == k
        wr = np.exp(rw.logit(MC["X"][te].astype(float), MC["m"][te]))
        nclip += int(((wr < rw.clip[0]) | (wr > rw.clip[1])).sum())
        w[te] = np.clip(wr, *rw.clip)
        aucs.append(rw.diag["val_auc_SB"])
    return w / w.mean(), {"val_auc_SB": aucs, "frac_clipped": nclip / max(len(w), 1)}


def cap(mask, rng, n=None):
    """Maska ograniczona losowo do n wierszy (FIT_CAP)."""
    n = n or FIT_CAP
    i = np.flatnonzero(mask)
    if len(i) <= n:
        return mask
    out = np.zeros(len(mask), bool)
    out[rng.choice(i, n, replace=False)] = True
    return out


def c2(Xa, Xb, wb=None, cfg=None, seed=0):
    c = cfg["c2st"]
    if len(Xa) < 200 or len(Xb) < 200:
        return None
    return _c2st(np.asarray(Xa, float), np.asarray(Xb, float), None, wb, n_folds=int(c["n_folds"]),
                 max_n=5000 if FAST else int(c["max_n"]), seed=seed)


def evaluate(G, D, MC, w, rng):
    """eps_B per bin pT i WP: dane (SB, pasmo zliczone; SR z interpolacji SB) vs MC surowe i ważone (SB, pasmo, SR)."""
    out = {}
    nw = len(SL.WPS)
    for i, bn in enumerate(SL.BIN_NAMES):
        d = np.ones(len(D["m"]), bool) if bn == "all" else D["ib"] == i
        q = np.ones(len(MC["m"]), bool) if bn == "all" else MC["ib"] == i
        Hd = G.cat_hist(D["m"][d & (D["reg"] != "SR")], D["ncat"][d & (D["reg"] != "SR")])     # SR danych nie wchodzi nigdzie
        r = {"n_data_SB": int((d & (D["reg"] == "SB")).sum()), "n_data_band": int((d & (D["reg"] == "band")).sum()),
             "n_mc_SR": int((q & (MC["reg"] == "SR")).sum()), "n_mc_SB": int((q & (MC["reg"] == "SB")).sum())}
        ed = {mth: G.eps_all(Hd, mth) for mth in ("pol2", "pol1", "sym", "lr")}
        bs = np.array([G.eps_all(rng.poisson(Hd).astype(float), "pol2") for _ in range(NBOOT)])
        r["data_SR"] = {**{mth: ed[mth].tolist() for mth in ed}, "pol2_stat": np.nanstd(bs, 0).tolist(),
                        "syst_interp": np.nanmax(np.abs(np.array([ed[mth] for mth in ("pol1", "sym", "lr")]) - ed["pol2"]), 0).tolist(),
                        "n_pred_all": G.n_sr(Hd.sum(0)), "n_sb_bins": float(Hd.sum(0)[G.sb].sum())}
        for reg in ("SB", "band"):
            dd = d & (D["reg"] == reg)
            n = max(int(dd.sum()), 1)
            e = [float((D["ncat"][dd] > j).mean()) if dd.any() else float("nan") for j in range(nw)]
            r[f"data_{reg}"] = {"eps": e, "err": [float(np.sqrt(max(x * (1 - x), 0) / n)) for x in e]}
        for tag, ww in (("raw", np.ones(len(w))), ("salad", w)):
            m = {}
            for reg in ("SB", "band", "SR"):
                qq = q & (MC["reg"] == reg)
                ee = [SL.weff(ww[qq], MC["ncat"][qq] > j) for j in range(nw)]
                m[reg] = {"eps": [x[0] for x in ee], "err": [x[1] for x in ee], "neff_over_N": SL.neff(ww[qq]) / max(int(qq.sum()), 1)}
            Hm = G.cat_hist(MC["m"][q], MC["ncat"][q], ww[q])
            m["pol2_on_mc"] = G.eps_all(Hm, "pol2").tolist()   # estymator na MC (obciążenie interpolacji)
            m["n_SR_over_SBbins"] = float(ww[q & (MC["reg"] == "SR")].sum() / max(Hm.sum(0)[G.sb].sum(), 1e-12))
            r[f"mc_{tag}"] = m
        r["data_nSR_pred_over_SBbins"] = r["data_SR"]["n_pred_all"] / max(r["data_SR"]["n_sb_bins"], 1e-12)
        out[bn] = r
    return out


def sculpt(G, D, MC, w):
    """Rzeźbienie masy przez wagi: korelacja w z m, średnie wagi w regionach, chi2 kształtu masy dane vs MC w binach SB."""
    sb = MC["reg"] == "SB"
    o = {"corr_w_m_SB": float(np.corrcoef(w[sb], MC["m"][sb])[0, 1]) if sb.sum() > 10 else None,
         "mean_w": {r: float(w[MC["reg"] == r].mean()) if (MC["reg"] == r).any() else None for r in ("SR", "band", "SB")}}
    hd = np.histogram(D["m"][D["reg"] == "SB"], G.edges)[0].astype(float)
    for tag, ww in (("raw", np.ones(len(w))), ("salad", w)):
        hm = np.histogram(MC["m"], G.edges, weights=ww)[0]
        h2 = np.histogram(MC["m"], G.edges, weights=ww ** 2)[0]
        s = hd[G.sb].sum() / max(hm[G.sb].sum(), 1e-12)
        chi = ((hd - s * hm) ** 2 / np.maximum(hd + s * s * h2, 1e-12))[G.sb]
        L, R = G.sb & (G.cen < G.m0), G.sb & (G.cen > G.m0)
        o[f"mass_SB_{tag}"] = {"chi2": float(chi.sum()), "ndf": int(G.sb.sum() - 1),
                               "left_over_right_data": float(hd[L].sum() / max(hd[R].sum(), 1)),
                               "left_over_right_mc": float(hm[L].sum() / max(hm[R].sum(), 1e-12))}
    return o


def c2st_block(cfg, D, MC, w, rng):
    o = {}
    for reg in ("SB", "band"):
        a, b = D["X"][D["reg"] == reg], MC["X"][MC["reg"] == reg]
        o[f"raw_{reg}"] = c2(a, b, None, cfg)
        o[f"salad_{reg}"] = c2(a, b, w[MC["reg"] == reg], cfg)
    a = D["X"][D["reg"] == "SB"]
    h = rng.random(len(a)) < 0.5
    o["floor_SB"] = c2(a[h], a[~h], None, cfg)
    return o


def neff_block(MC, w):
    o = {}
    for reg in ("SR", "band", "SB"):
        o[reg] = {bn: SL.neff(w[(MC["reg"] == reg) & ((MC["ib"] == i) | (bn == "all"))]) /
                  max(int(((MC["reg"] == reg) & ((MC["ib"] == i) | (bn == "all"))).sum()), 1) for i, bn in enumerate(SL.BIN_NAMES)}
    return o


# ----------------------------------------------------------------------------------------------------- zadanie main
def task_main(a, cfg, sp, spp, feats, G, S, MCall, DATA, thr, rng):
    tmin = SV.mc_truth_min(cfg, sp)
    bkg = (MCall["truth"] >= tmin) & ~np.isin(MCall["truth"], [0, 1])
    MC = sel(MCall, bkg)
    res = {"classes": {}, "transfer": {}, "mc_truth_min": tmin, "n_mc_bkg": int(bkg.sum()),
           "truth_counts_mc": {int(k): int(v) for k, v in zip(*np.unique(MCall["truth"], return_counts=True))}}
    W = {}
    for c in SL.CLASSES:
        D = DATA[c]
        nsb = int((D["reg"] == "SB").sum())
        r = {"n_data": int(len(D["m"])), "n_data_SB": nsb, "n_data_band": int((D["reg"] == "band").sum())}
        if nsb < MIN_SB:
            r["skipped"] = f"N_SB = {nsb} < {MIN_SB}"
            res["classes"][c] = r
            say(c, "pominięta:", r["skipped"])
            continue
        # bieg „prod” = v1.0 jak zamrożony: jeden fit na całym SB, bramki z yaml
        rw = SV.fit_salad_v1(cfg, sp, frame(D, feats, D["reg"] != "SR"), frame(MC, feats), kind="hgb", seed=0, c2st_band=True)
        r["prod"] = {k: rw.diag.get(k) for k in ("val_auc_SB", "neff_over_N_SR", "w_quantiles_SR", "c2st_band", "ok", "n_fit_mc",
                                               "n_fit_data", "sha256_mc", "sha256_data_sb")}
        wp = SV.salad_weights(rw, frame(MC, feats))
        r["prod"]["eval"] = evaluate(G, D, MC, wp, rng)
        say(c, "prod", {k: r["prod"][k] for k in ("val_auc_SB", "neff_over_N_SR", "c2st_band", "ok")})
        r["kfold"] = {}
        for seed in SEEDS:
            w, dg = oof_weights(cfg, sp, D, MC, feats, seed)
            blk = {"diag": dg, "eval": evaluate(G, D, MC, w, rng), "neff_over_N": neff_block(MC, w)}
            if seed == SEEDS[0]:
                W[c] = w
                blk["c2st"] = c2st_block(cfg, D, MC, w, rng)
                blk["sculpt"] = sculpt(G, D, MC, w)
            r["kfold"][f"seed{seed}"] = blk
            e = blk["eval"]["all"]
            say(c, f"kfold seed {seed}: AUC_SB {np.mean(dg['val_auc_SB']):.3f} n_eff/N(SR) {blk['neff_over_N']['SR']['all']:.2f}",
                "SALAD/pol2 @WP:", [round(x / y, 3) if y else None for x, y in zip(e["mc_salad"]["SR"]["eps"], e["data_SR"]["pol2"])],
                "C2ST", blk.get("c2st"))
        res["classes"][c] = r
        json.dump(res, open(os.path.join(a.out, f"s427_{sp}.json.part"), "w"))
    # transfer wag a -> b
    for ca, w in W.items():
        res["transfer"][ca] = {}
        for cb in W:
            if cb == ca:
                continue
            D = DATA[cb]
            t = {"eval": evaluate(G, D, MC, w, rng)}
            for reg in ("SB", "band"):
                t[f"c2st_salad_{reg}"] = c2(D["X"][D["reg"] == reg], MC["X"][MC["reg"] == reg], w[MC["reg"] == reg], cfg)
            res["transfer"][ca][cb] = t
            e = t["eval"]["all"]
            say(f"transfer {ca} -> {cb}: SALAD/pol2", [round(x / y, 3) if y else None for x, y in zip(e["mc_salad"]["SR"]["eps"], e["data_SR"]["pol2"])],
                "C2ST SB", t["c2st_salad_SB"])
    # dane vs dane (18 cech, SB): jak różne są klasy
    res["c2st_data_data_SB"] = {}
    ref = DATA["MB427929"]
    for c in SL.CLASSES[1:]:
        if c in W and "MB427929" in W:
            res["c2st_data_data_SB"][c] = c2(DATA[c]["X"][DATA[c]["reg"] == "SB"], ref["X"][ref["reg"] == "SB"], None, cfg)
    # D0: wariant kontrolny MC truth >= 8 (bez swapów) dla MB427929
    if sp == "D0" and "MB427929" in W:
        MC8 = sel(MCall, (MCall["truth"] >= 8) & (MCall["truth"] < 100))
        w8, dg = oof_weights(cfg, sp, DATA["MB427929"], MC8, feats, 0)
        res["variant_truth8_MB427929"] = {"diag": dg, "eval": evaluate(G, DATA["MB427929"], MC8, w8, rng), "neff_over_N": neff_block(MC8, w8)}
    np.savez_compressed(os.path.join(a.out, f"s427_{sp}_w.npz"), uid=MC["uid"], m=MC["m"].astype(np.float32),
                        **{f"w_{c}": w.astype(np.float32) for c, w in W.items()})
    return res


# ----------------------------------------------------------------------------------------------------- zadanie adv (B1)
CAT_F = {"MisMatch": 1.5, "OtherNoB": 0.8, "OtherB": 1.5}
ADV_MODES = ("smooth", "cat_global", "cat_local", "cat_step", "cuts_local")


def cat_factor(t):
    return np.where((t >= 8) & (t <= 11), CAT_F["MisMatch"], np.where(t == 12, CAT_F["OtherNoB"], np.where((t >= 13) & (t <= 15), CAT_F["OtherB"], 1.0)))


def w_smooth(B, feats, m0, strength=0.8):
    t = (B["m"] - m0) / 0.1
    co = {"Lxy": (0.3, 0.0), "chi2": (-0.25, 0.0), "log_pt": (0.1, 0.3), "dedx_K": (0.2, 0.0), "abs_impact_sig": (0.0, -0.2),
          "d_max_abs_d0": (0.15, 0.1), "log_p_pi_min": (0.0, 0.15)}
    lw = np.zeros(len(t))
    for f, (c0, c1) in co.items():
        v = B["X"][:, feats.index(f)].astype(float)
        v = np.where(np.isfinite(v), v, np.nanmedian(v))
        z = (v - v.mean()) / (v.std() + 1e-12)
        lw += strength * (c0 + c1 * t) * np.clip(z, -4, 4)
    return np.minimum(np.exp(lw - np.quantile(lw, 0.97)), 1.0)


def w_mode(B, feats, mode, m0, sigma=0.035):
    m, t = B["m"], B["truth"]
    g = np.exp(-0.5 * ((m - m0) / sigma) ** 2)
    if mode == "smooth":
        return w_smooth(B, feats, m0)
    if mode == "cat_global":
        return cat_factor(t) / 1.5
    if mode == "cat_local":
        return (1 + (cat_factor(t) - 1) * g) / 1.5
    if mode == "cat_step":
        return np.where(((t >= 8) & (t <= 11)) & (m < m0), 1.5, 1.0) / 1.5
    if mode == "cuts_local":
        f = np.where(B["X"][:, feats.index("Lxy")] > 0.3, 1.5, 1.0) * np.where(B["X"][:, feats.index("abs_impact_sig")] > 3, 0.7, 1.0)
        return (1 + (f - 1) * g) / 1.5
    raise ValueError(mode)


def task_adv(a, cfg, sp, spp, feats, G, MCall, n_sb_target, rng):
    tmin = SV.mc_truth_min(cfg, sp)
    MC = sel(MCall, (MCall["truth"] >= tmin) & ~np.isin(MCall["truth"], [0, 1]))
    par = (SL.mix64(MC["uid"], 11) % np.uint64(2)).astype(int)
    A, B = sel(MC, par == 0), sel(MC, par == 1)
    mm = lambda X: (X["truth"] >= 8) & (X["truth"] <= 11)  # noqa: E731
    out = {"n_A": int(len(A["m"])), "n_sb_target": int(n_sb_target), "modes": {}, "cat_factors": CAT_F, "sigma_local": 0.035}
    for mode in ADV_MODES:
        P = sel(B, rng.random(len(B["m"])) < w_mode(B, feats, mode, spp["m0"]))
        nsb = int((P["reg"] == "SB").sum())
        if n_sb_target and nsb > n_sb_target:                       # statystyka pseudo-danych jak SB klasy MB427929
            P = sel(P, rng.random(len(P["m"])) < n_sb_target / nsb)
        r = {"n_pseudo": int(len(P["m"])), "n_pseudo_SB": int((P["reg"] == "SB").sum()), "n_pseudo_SR": int((P["reg"] == "SR").sum()),
             "f_MM_SR_pseudo": float(mm(P)[P["reg"] == "SR"].mean()), "f_MM_SB_pseudo": float(mm(P)[P["reg"] == "SB"].mean()),
             "f_MM_SR_A": float(mm(A)[A["reg"] == "SR"].mean())}
        rw = SV.fit_salad_v1(cfg, sp, frame(P, feats, P["reg"] != "SR"), frame(A, feats), kind="hgb", seed=0, c2st_band=False)
        w = SV.salad_weights(rw, frame(A, feats))
        r.update(val_auc_SB=rw.diag["val_auc_SB"], neff_over_N_SR=rw.diag["neff_over_N_SR"])
        for reg in ("SB", "SR", "band"):
            xa, xb = P["X"][P["reg"] == reg], A["X"][A["reg"] == reg]
            r[f"c2st_raw_{reg}"] = c2(xa, xb, None, cfg)
            r[f"c2st_salad_{reg}"] = c2(xa, xb, w[A["reg"] == reg], cfg)
        xs = P["X"][P["reg"] == "SR"]
        h = rng.random(len(xs)) < 0.5
        r["c2st_floor_SR"] = c2(xs[h], xs[~h], None, cfg)
        r["eps"] = {}
        for i, bn in enumerate(SL.BIN_NAMES):
            p = (P["reg"] == "SR") & ((P["ib"] == i) | (bn == "all"))
            q = (A["reg"] == "SR") & ((A["ib"] == i) | (bn == "all"))
            row = {"true": [], "true_err": [], "salad": [], "salad_err": [], "mcraw": []}
            for j in range(len(SL.WPS)):
                et, st = SL.weff(np.ones(p.sum()), P["ncat"][p] > j)
                es, ss = SL.weff(w[q], A["ncat"][q] > j)
                er, _ = SL.weff(np.ones(q.sum()), A["ncat"][q] > j)
                for k_, v in zip(row, (et, st, es, ss, er)):
                    row[k_].append(v)
            r["eps"][bn] = row
        e = r["eps"]["all"]
        say(f"adv {mode}: C2ST SALAD SR {r['c2st_salad_SR']} pred/true @WP", [round(x / y, 3) if y else None for x, y in zip(e["salad"], e["true"])],
            "f_MM", round(r["f_MM_SR_pseudo"], 4), round(r["f_MM_SR_A"], 4))
        out["modes"][mode] = r
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True, choices=SL.SPECIES)
    ap.add_argument("--parts", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tasks", default="main,adv")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    sp = a.species
    cfg = SV.load_cfg()
    if FAST:
        cfg["hgb"]["max_iter"] = 30
        cfg["min_n_data_sb"] = 0
    spp = SV.species_params(cfg, sp)
    feats = list(cfg["features"])
    G = SL.MassGrid(spp)
    rng = np.random.default_rng(427929)
    S = prep(load(a.parts, "fd_nd", sp), spp, sp)
    MCall = prep(load(a.parts, "mb_nd", sp), spp, sp)
    d427 = prep(load(a.parts, "data_427929", sp), spp, sp)
    d428 = prep(load(a.parts, "data_428770", sp), spp, sp)
    DATA = {"MB427929": sel(d427, SL.class_masks(d427["tbits"], d427["grl"], "data_427929")["MB427929"])}
    for c, m in SL.class_masks(d428["tbits"], d428["grl"], "data_428770").items():
        DATA[c] = sel(d428, m)
    del d427, d428
    say(sp, "wczytane: fd", len(S["m"]), "mb", len(MCall["m"]), {c: len(D["m"]) for c, D in DATA.items()})
    thr, aucs = tagger(S, DATA["MB427929"], [MCall] + [DATA[c] for c in SL.CLASSES[1:]], rng)
    meta = {"version": SL.VERSION, "species": sp, "fast": FAST, "k_folds": K, "nboot": NBOOT, "seeds": list(SEEDS), "wps": list(SL.WPS),
            "pt_bins": SL.PT_BINS, "regions": {k: (list(v) if isinstance(v, tuple) else v) for k, v in spp.items()},
            "salad_cfg_version": cfg.get("version"), "feats": feats, "thr": thr.tolist(), "tagger_auc_oof": aucs,
            "n": {"fd": int(len(S["m"])), "mb": int(len(MCall["m"])), **{c: int(len(D["m"])) for c, D in DATA.items()}},
            "gates": {k: cfg.get(k) for k in ("min_neff_over_N_SR", "min_n_data_sb", "max_c2st_band")}, "refl_veto": SL.REFL_VETO if sp == "D0" else None}
    tasks = a.tasks.split(",")
    if "main" in tasks:
        res = task_main(a, cfg, sp, spp, feats, G, S, MCall, DATA, thr, rng)
        res["meta"] = {**meta, "wall_s": round(time.time() - T0)}
        json.dump(res, open(os.path.join(a.out, f"s427_{sp}.json"), "w"), indent=1)
        say("zapisano", f"s427_{sp}.json")
    if "adv" in tasks:
        adv = task_adv(a, cfg, sp, spp, feats, G, MCall, int((DATA["MB427929"]["reg"] == "SB").sum()), rng)
        adv["meta"] = {**meta, "wall_s": round(time.time() - T0)}
        json.dump(adv, open(os.path.join(a.out, f"s427_{sp}_adv.json"), "w"), indent=1)
        say("zapisano", f"s427_{sp}_adv.json")


if __name__ == "__main__":
    sys.exit(main())
