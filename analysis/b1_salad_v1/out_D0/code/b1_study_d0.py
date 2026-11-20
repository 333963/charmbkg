#!/usr/bin/env python
"""B1-D0 — closure SALAD v1.0 dla D0 (podpróbki S_bkggen/sub, po wecie refleksji |m_swap−m0|>25 MeV) gładki vs ADWERSARIALNY.
Ta sama logika i metryki co b1_study.py (D⁺, etap B1); różnice: --species (domyślnie D0: SR |Δm|<55, pasmo 55–75, SB >75 MeV),
okno wycieku 75–95 MeV, tryby dodatkowe specyficzne dla D0:
  swap_local      : pseudo-dane = tło (truth ≥ 8, parzystość 1) + swapy K↔π (truth 100/101; pula: mb par. 1 + fd par. 1) w udziale
                    nominalnym f_nom(region) = N_swap/N_bkg z mb (po wecie) × (1 + (fac−1)·g(m)), --swap_factor 1.5 | 0.5 (σ 35 MeV)
  swap_plus       : jak wyżej, f = f_nom + 0.01·g(m)  (+1 pp swapów pod pikiem → czułość Δε_B / pp)
  real_t27        : pseudo-dane = wszystko, co w danych nie jest sygnałem: truth 2–7 (częściowo dopasowane) + ≥ 8 + 100/101 (parz. 1),
                    bez zniekształcenia; MC SALAD = truth ≥ 8 (parz. 0) — realistyczny skład D0 (truth 2–7 ≈ 1 % tła w SR)
Dodatkowo w każdym biegu: ε_B@WP per kategoria truth w SR (MC A, surowe) → czułość analityczna (ε_c − ε_B)/ε_B na +1 pp.
"""
from __future__ import annotations
import argparse, json, os, sys, time, warnings
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ.get("CHARMTAG_SRC", "/work/b1d0/in/src")
sys.path.insert(0, SRC)
warnings.filterwarnings("ignore")
from charmtag.bkggen import salad_v1 as SV  # noqa: E402
from charmtag.bkggen.feta import c2st as _c2st  # noqa: E402

CFG = SV.load_cfg(); CFG["c2st"]["max_n"] = 60_000
H = list(CFG["features"]); PT_BINS = [3, 5, 8, 30]; WPS = (0.8, 0.7, 0.5)
SUB = os.path.join(HERE, "sub")
SPN, SP, M0 = "D0", None, None
CATS = {"MisMatch": (8, 11), "OtherNoB": (12, 12), "OtherB": (13, 15)}


def set_species(sp):
    global SPN, SP, M0
    SPN = sp; SP = SV.species_params(CFG, sp); M0 = SP["m0"]


def load(name):
    z = np.load(f"{SUB}/{name}.npz")
    return pd.DataFrame({k: z[k] for k in z.files})


def c2st(Xa, Xb, wa=None, wb=None, seed=0):
    c = CFG["c2st"]
    return _c2st(np.asarray(Xa, float), np.asarray(Xb, float), wa, wb, n_folds=c["n_folds"], max_n=c["max_n"], seed=seed)


def hgb(X, y, seed=0, it=200):
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(max_iter=it, learning_rate=0.08, max_leaf_nodes=31, l2_regularization=1.0,
                                          early_stopping=True, validation_fraction=0.15, random_state=seed).fit(X, y)


def wp_thresholds(s, pt):
    return {(a, b): {wp: float(np.quantile(s[(pt >= a) & (pt < b)], 1 - wp)) for wp in WPS} for a, b in zip(PT_BINS[:-1], PT_BINS[1:])}


def eff(score, pt, thr, w=None, wp=0.8):
    w = np.ones(len(score)) if w is None else w
    out, tn, tk = {}, 0.0, 0.0
    for (a, b), t in thr.items():
        s = (pt >= a) & (pt < b); n, k = w[s].sum(), w[s & (score >= t[wp])].sum()
        out[f"{a}-{b}"] = float(k / n) if n > 0 else float("nan"); tn += n; tk += k
    out["all"] = float(tk / tn) if tn > 0 else float("nan")
    return out


# --------------------------------------------------------------------------- zniekształcenia (jak B1)
def w_smooth(B, feats, strength=0.8):
    t = (B["m"].to_numpy(float) - M0) / 0.1
    co = {"Lxy": (0.3, 0.0), "chi2": (-0.25, 0.0), "log_pt": (0.1, 0.3), "dedx_K": (0.2, 0.0),
          "abs_impact_sig": (0.0, -0.2), "d_max_abs_d0": (0.15, 0.1), "log_p_pi_min": (0.0, 0.15)}
    lw = np.zeros(len(B))
    for f, (c0, c1) in co.items():
        if f in feats:
            v = B[f].to_numpy(float); v = np.where(np.isfinite(v), v, np.nanmedian(v)); z = (v - v.mean()) / (v.std() + 1e-12)
            lw += strength * (c0 + c1 * t) * np.clip(z, -4, 4)
    return np.minimum(np.exp(lw - np.quantile(lw, 0.97)), 1.0)


CAT_F = {"MisMatch": 1.5, "OtherNoB": 0.8, "OtherB": 1.5}


def cat_factor(truth, f=CAT_F):
    t = np.asarray(truth)
    return np.where((t >= 8) & (t <= 11), f["MisMatch"], np.where(t == 12, f["OtherNoB"], np.where((t >= 13) & (t <= 15), f["OtherB"], 1.0)))


def g_local(m, sigma=0.035):
    return np.exp(-0.5 * ((np.asarray(m, float) - M0) / sigma) ** 2)


STRENGTH = [0.8]


def w_mode(B, feats, mode, sigma=0.035):
    m, t = B["m"].to_numpy(float), B["truth"].to_numpy()
    if mode in ("smooth",):
        return w_smooth(B, feats, STRENGTH[0])
    if mode == "cat_global":
        return cat_factor(t) / 1.5
    if mode == "cat_local":
        return (1 + (cat_factor(t) - 1) * g_local(m, sigma)) / 1.5
    if mode == "cat_local_MM":      # tylko MisMatch* ×2 pod pikiem (czułość na jedną kategorię)
        f = np.where((t >= 8) & (t <= 11), 2.0, 1.0)
        return (1 + (f - 1) * g_local(m, sigma)) / 2.0
    if mode == "cat_step":
        f = np.where(((t >= 8) & (t <= 11)) & (m < M0), 1.5, 1.0)
        return f / 1.5
    if mode == "cat_local+smooth":
        return w_smooth(B, feats, STRENGTH[0]) * (1 + (cat_factor(t) - 1) * g_local(m, sigma)) / 1.5
    if mode == "cuts_local":
        f = np.where(B["Lxy"].to_numpy(float) > 0.3, 1.5, 1.0) * np.where(B["abs_impact_sig"].to_numpy(float) > 3, 0.7, 1.0)
        return (1 + (f - 1) * g_local(m, sigma)) / 1.5
    if mode in ("swap_local", "swap_plus", "real_t27", "none"):
        return np.ones(len(B))
    raise ValueError(mode)


def make_pseudo(B1, feats, mode, rng, sigma=0.035):
    w = w_mode(B1, feats, mode, sigma)
    keep = rng.random(len(B1)) < w
    return B1[keep].reset_index(drop=True)


def add_swaps(P, pool, f_nom, rng, factor=1.0, plus=0.0, sigma=0.035):
    """Dodaj swapy do pseudo-danych: w binach m (10 MeV) N_sw(bin) = N_P(bin)·f_nom(region)·(1+(factor−1)g) + N_P(bin)·plus·g."""
    edges = np.arange(SP["mass_range"][0], SP["mass_range"][1] + 1e-9, 0.010)
    mP, mW = P["m"].to_numpy(float), pool["m"].to_numpy(float)
    rW = SV.regions(mW, SP)
    take = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mc = 0.5 * (lo + hi); reg = SV.regions([mc], SP)[0]
        nP = int(((mP >= lo) & (mP < hi)).sum()); g = float(g_local(mc, sigma))
        n = int(round(nP * (f_nom[reg] * (1 + (factor - 1) * g) + plus * g)))
        cand = np.flatnonzero((mW >= lo) & (mW < hi))
        if n > 0 and len(cand):
            take.append(rng.choice(cand, n, replace=n > len(cand)))
    idx = np.concatenate(take) if take else np.array([], int)
    I = pool.iloc[idx].reset_index(drop=True)
    return pd.concat([P, I[P.columns]], ignore_index=True), int(len(I)), int(len(np.unique(idx)))


def inject_signal(P, S_inj, frac, rng, win=(0.075, 0.095)):
    m = P["m"].to_numpy(float); edge = (np.abs(m - M0) >= win[0]) & (np.abs(m - M0) < win[1])
    n = int(round(frac * edge.sum() / (1 - frac)))
    idx = rng.choice(len(S_inj), n, replace=n > len(S_inj))
    I = S_inj.iloc[idx].reset_index(drop=True)
    I["m"] = M0 + rng.choice([-1, 1], n) * rng.uniform(win[0], win[1], n)
    I["truth"] = -1
    return pd.concat([P, I[P.columns]], ignore_index=True), n


def frac_cat(P, r, reg, lo, hi):
    t = P["truth"].to_numpy(); s = r == reg
    return float(((t >= lo) & (t <= hi))[s].mean()) if s.any() else float("nan")


# --------------------------------------------------------------------------- closure
def closure(a):
    t0 = time.time(); rng = np.random.default_rng(a.seed)
    S, Ball = load("sig"), load("mcbkg")
    feats = [f for f in H if f in Ball and np.isfinite(Ball[f].to_numpy(float)).any()]
    Bbkg = Ball[(Ball["truth"] >= 8) & (Ball["truth"] < 100)].reset_index(drop=True)
    A = Ball[(Ball["parity"] == 0) & (Ball["truth"] >= a.mc_truth_min)].reset_index(drop=True)   # MC SALAD + „MC raw”: truth ≥ mc_truth_min (8: bez swapów; 2: + częściowo dopasowane 2–7 i swapy 100/101)
    if a.mc_truth_min >= 8:
        A = A[A["truth"] < 100].reset_index(drop=True)
    tag = f"{a.mode}_{a.kind}" + (f"_s{a.strength:g}" if a.strength != 0.8 else "") + (f"_leak{a.leak:g}" if a.leak else "") + (f"_x{a.swap_factor:g}" if a.mode == "swap_local" else "") + (f"_mc{a.mc_truth_min}" if a.mc_truth_min != 8 else "") + (f"_seed{a.seed}" if a.seed else "")
    res = {"species": SPN, "mode": a.mode, "kind": a.kind, "leak": a.leak, "sigma": a.sigma, "swap_factor": a.swap_factor, "strength": a.strength, "mc_truth_min": a.mc_truth_min, "n_A": int(len(A)), "seed": a.seed}
    if a.mode == "real_t27":
        P = Ball[(Ball["parity"] == 1) & (Ball["truth"] >= 2)].reset_index(drop=True)   # 2–7, 8–15, 100/101 (bez sygnału 0/1)
    else:
        P = make_pseudo(Bbkg[Bbkg["parity"] == 1].reset_index(drop=True), feats, a.mode, rng, a.sigma)
    if a.mode in ("swap_local", "swap_plus"):
        rB = SV.regions(Ball["m"], SP); isw = Ball["truth"].isin([100, 101]).to_numpy(); ib = ((Ball["truth"] >= 8) & (Ball["truth"] < 100)).to_numpy()
        f_nom = {reg: float((isw & (rB == reg)).sum() / max((ib & (rB == reg)).sum(), 1)) for reg in ("SR", "band", "SB")}
        pool = pd.concat([Ball[isw & (Ball["parity"] == 1)], S[S["truth"].isin([100, 101]) & (S["parity"] == 1)]], ignore_index=True)
        P, n_sw, n_uniq = add_swaps(P, pool, f_nom, rng, factor=a.swap_factor if a.mode == "swap_local" else 1.0,
                                    plus=0.01 if a.mode == "swap_plus" else 0.0, sigma=a.sigma)
        res.update(f_swap_nom=f_nom, n_swaps_added=n_sw, n_swaps_unique=n_uniq, n_swap_pool=int(len(pool)))
    res["n_pseudo"] = int(len(P))
    S0 = S[(S.parity == 0) & S.truth.isin([0, 1])].reset_index(drop=True)
    S1 = S[(S.parity == 1) & S.truth.isin([0, 1])].reset_index(drop=True)
    n_inj = 0; P_fit = P
    if a.leak:
        P_fit, n_inj = inject_signal(P, S1, a.leak, rng, win=(SP["sb_min"], SP["sb_min"] + 0.020))
    res["n_injected"] = int(n_inj)
    rA, rP, rF = SV.regions(A["m"], SP), SV.regions(P["m"], SP), SV.regions(P_fit["m"], SP)
    XA, XP = A[feats].to_numpy(float), P[feats].to_numpy(float)
    res["n_pseudo_SR"] = int((rP == "SR").sum())
    for reg in ("SR", "SB"):
        res[f"frac_MisMatch_{reg}_pseudo"] = frac_cat(P, rP, reg, 8, 11); res[f"frac_OtherB_{reg}_pseudo"] = frac_cat(P, rP, reg, 13, 15)
        res[f"frac_swap_{reg}_pseudo"] = frac_cat(P, rP, reg, 100, 101); res[f"frac_t27_{reg}_pseudo"] = frac_cat(P, rP, reg, 2, 7)
        res[f"frac_MisMatch_{reg}_A"] = frac_cat(A, rA, reg, 8, 11); res[f"frac_OtherB_{reg}_A"] = frac_cat(A, rA, reg, 13, 15)
    res["c2st_raw_SB"] = c2st(XP[rP == "SB"], XA[rA == "SB"], seed=a.seed)
    res["c2st_raw_SR"] = c2st(XP[rP == "SR"], XA[rA == "SR"], seed=a.seed)
    print(tag, "raw", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items() if "c2st" in k or "frac" in k or "swap" in k}, f"{time.time()-t0:.0f}s", flush=True)
    CFG["mc_truth_min"] = a.mc_truth_min
    rw = SV.fit_salad_v1(CFG, SPN, P_fit[rF != "SR"], A, kind=a.kind, seed=a.seed, c2st_band=False)
    wA = SV.salad_weights(rw, A)
    res.update({k: rw.diag[k] for k in ("val_auc_SB", "neff_over_N_SR", "w_quantiles_SR")})
    res["c2st_salad_SR"] = c2st(XP[rP == "SR"], XA[rA == "SR"], None, wA[rA == "SR"], seed=a.seed)
    res["c2st_salad_band"] = c2st(XP[rP == "band"], XA[rA == "band"], None, wA[rA == "band"], seed=a.seed)
    res["c2st_salad_SB"] = c2st(XP[rP == "SB"], XA[rA == "SB"], None, wA[rA == "SB"], seed=a.seed)
    res["c2st_sbcopy_SR"] = c2st(XP[rP == "SR"], XP[rP == "SB"], seed=a.seed)
    half = rng.random((rP == "SR").sum()) < 0.5
    res["c2st_floor_SR"] = c2st(XP[rP == "SR"][half], XP[rP == "SR"][~half], seed=a.seed)
    print(tag, "c2st", {k: round(v, 3) for k, v in res.items() if isinstance(v, float) and ("c2st" in k or "neff" in k)}, f"{time.time()-t0:.0f}s", flush=True)
    Xs = S0[feats].to_numpy(float); Xb = XP[rP == "SB"]
    if len(Xb) > 250_000: Xb = Xb[rng.choice(len(Xb), 250_000, replace=False)]
    clf = hgb(np.vstack([Xs, Xb]), np.r_[np.ones(len(Xs)), np.zeros(len(Xb))], seed=a.seed)
    sc = lambda X: clf.predict_proba(X)[:, 1]  # noqa: E731
    thr = wp_thresholds(sc(S1[feats].to_numpy(float)), S1["pt"].to_numpy(float))
    sA, sP = sc(XA), sc(XP); ptA, ptP = A["pt"].to_numpy(float), P["pt"].to_numpy(float)
    tA, tP = A["truth"].to_numpy(), P["truth"].to_numpy()
    for wp in WPS:
        tru = eff(sP[rP == "SR"], ptP[rP == "SR"], thr, None, wp)
        res[f"eff_true_SR@{wp}"] = tru
        res[f"eff_salad_SR@{wp}"] = eff(sA[rA == "SR"], ptA[rA == "SR"], thr, wA[rA == "SR"], wp)
        res[f"eff_mcraw_SR@{wp}"] = eff(sA[rA == "SR"], ptA[rA == "SR"], thr, None, wp)
        res[f"eff_sbcopy@{wp}"] = eff(sP[rP == "SB"], ptP[rP == "SB"], thr, None, wp)
        res[f"eff_salad_SB@{wp}"] = eff(sA[rA == "SB"], ptA[rA == "SB"], thr, wA[rA == "SB"], wp)
        res[f"eff_true_SB@{wp}"] = eff(sP[rP == "SB"], ptP[rP == "SB"], thr, None, wp)
        n = (rP == "SR").sum(); res[f"eff_true_relerr@{wp}"] = float(np.sqrt((1 - tru["all"]) / max(tru["all"] * n, 1)))
        ws = wA[rA == "SR"]; p = res[f"eff_salad_SR@{wp}"]["all"]
        neff = ws.sum() ** 2 / (ws ** 2).sum(); res[f"eff_salad_relerr@{wp}"] = float(np.sqrt((1 - p) / max(p * neff, 1)))
        # ε per kategoria w SR (MC A surowe, + swapy/t2-7 z pseudo) → czułość analityczna na +1 pp
        cat_eff = {}
        for name, (lo, hi) in CATS.items():
            s = (rA == "SR") & (tA >= lo) & (tA <= hi)
            cat_eff[name] = eff(sA[s], ptA[s], thr, None, wp)["all"] if s.sum() > 50 else float("nan")
        for name, (lo, hi) in (("swap", (100, 101)), ("t27", (2, 7))):
            s = (rP == "SR") & (tP >= lo) & (tP <= hi)
            cat_eff[name] = eff(sP[s], ptP[s], thr, None, wp)["all"] if s.sum() > 50 else float("nan")
        cat_eff["bkg_pseudo_SR"] = res[f"eff_true_SR@{wp}"]["all"]
        cat_eff["bkg_A"] = res[f"eff_mcraw_SR@{wp}"]["all"]
        res[f"eff_cat_SR@{wp}"] = cat_eff
        res[f"sens_per_pp@{wp}"] = {k: float((v - cat_eff["bkg_A"]) / cat_eff["bkg_A"] * 0.01) for k, v in cat_eff.items() if k != "bkg_A"}
    sl = np.zeros(len(A), bool)
    for (a_, b_), t in thr.items():
        sl |= (ptA >= a_) & (ptA < b_) & (sA >= t[0.8])
    res["mean_w_siglike_SR"] = float(wA[(rA == "SR") & sl].mean()); res["mean_w_siglike_SB"] = float(wA[(rA == "SB") & sl].mean())
    res["mean_w_siglike_band"] = float(wA[(rA == "band") & sl].mean())
    print(tag, {k: round(v["all"], 4) for k, v in res.items() if k.startswith("eff_") and isinstance(v, dict) and "all" in v}, f"{time.time()-t0:.0f}s", flush=True)
    print(tag, "cat_eff@0.8", {k: round(v, 4) for k, v in res["eff_cat_SR@0.8"].items()}, flush=True)
    os.makedirs(a.out, exist_ok=True)
    json.dump(res, open(f"{a.out}/closure_{tag}.json", "w"), indent=1)
    if a.leak or a.mode == "smooth":
        np.savez_compressed(f"{a.out}/wSR_{tag}.npz", w=wA[rA == "SR"].astype(np.float32))


def report(a):
    import glob
    files = sorted(glob.glob(f"{a.out}/closure_*.json"))
    R = [json.load(open(f)) for f in files]
    out = [f"# B1-D0 — closure SALAD v1.0 ({SPN}, podpróbki S_bkggen/sub po wecie refleksji; pseudo-dane = MC parzystość 1 zniekształcone; prawda = pseudo SR)\n"]
    out.append("| tryb | kind | siła | wyciek | swap× | MC truth≥ | seed | n_pseudo SR | f_MM SR pseudo / A | f_swap SR | C2ST raw SB | raw SR | SALAD SB | **SALAD SR** | SALAD pasmo | SB-kopia SR | podłoga | n_eff/N | AUC_SB | ε_B pred/true @80 (SALAD / SB-kopia / MC raw) | @70 | @50 | σ_stat prawdy @80/70/50 | w̄ sygn.-pod. SR/pasmo/SB |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in R:
        cells = []
        for wp in WPS:
            t = r[f"eff_true_SR@{wp}"]["all"]
            cells.append(" / ".join(f"{r[f'eff_{s}@{wp}']['all'] / t:.3f}" for s in ("salad_SR", "sbcopy", "mcraw_SR")))
        se = " / ".join(f"{r[f'eff_true_relerr@{wp}']:.3f}" for wp in WPS)
        out.append(f"| {r['mode']} | {r['kind']} | {r.get('strength', 0.8):g} | {r['leak'] or 0:g} | {r.get('swap_factor', 1):g} | {r.get('mc_truth_min', 8)} | {r.get('seed', 0)} | {r['n_pseudo_SR']} | {r['frac_MisMatch_SR_pseudo']:.3f} / {r['frac_MisMatch_SR_A']:.3f} | {r['frac_swap_SR_pseudo']:.5f} | "
                   f"{r['c2st_raw_SB']:.3f} | {r['c2st_raw_SR']:.3f} | {r['c2st_salad_SB']:.3f} | **{r['c2st_salad_SR']:.3f}** | {r['c2st_salad_band']:.3f} | "
                   f"{r['c2st_sbcopy_SR']:.3f} | {r['c2st_floor_SR']:.3f} | {r['neff_over_N_SR']:.2f} | {r['val_auc_SB']:.3f} | "
                   f"{cells[0]} | {cells[1]} | {cells[2]} | {se} | {r['mean_w_siglike_SR']:.3f} / {r['mean_w_siglike_band']:.3f} / {r['mean_w_siglike_SB']:.3f} |")
    out.append("\n## Per bin pT: ε_B(SR) SALAD/prawda @80 / @70\n")
    out.append("| tryb | kind | wyciek | swap× | 3-5 | 5-8 | 8-30 | all |"); out.append("|---|---|---|---|---|---|---|---|")
    for r in R:
        cells = [" / ".join(f"{r[f'eff_salad_SR@{wp}'][b] / r[f'eff_true_SR@{wp}'][b]:.3f}" for wp in (0.8, 0.7)) for b in ("3-5", "5-8", "8-30", "all")]
        out.append(f"| {r['mode']} | {r['kind']} | {r['leak'] or 0:g} | {r.get('swap_factor', 1):g} | " + " | ".join(cells) + " |")
    out.append("\n## ε_B@WP per kategoria w SR (MC A surowe; swap z biegu swap_plus, t2-7 z real_t27) i czułość analityczna Δε_B/ε_B na +1 pp\n")
    pick = lambda mode, kind="hgb": next((r for r in R if r["mode"] == mode and r["kind"] == kind and not r["leak"] and r.get("seed", 0) == 0 and r.get("strength", 0.8) == 0.8 and r.get("mc_truth_min", 8) == 8), None)  # noqa: E731
    r0, rs, rt = pick("smooth"), pick("swap_plus"), pick("real_t27")
    if r0:
        out.append("| WP | bkg (A) | MisMatch | OtherNoB | OtherB | swap (po wecie) | t2-7 | czułość/pp: MM / OtherB / swap / t2-7 |"); out.append("|---|---|---|---|---|---|---|---|")
        for wp in WPS:
            ce, s = r0[f"eff_cat_SR@{wp}"], r0[f"sens_per_pp@{wp}"]
            esw = rs[f"eff_cat_SR@{wp}"]["swap"] if rs else float("nan"); et = rt[f"eff_cat_SR@{wp}"]["t27"] if rt else float("nan")
            b = ce["bkg_A"]
            out.append(f"| {wp} | {b:.4f} | {ce['MisMatch']:.4f} | {ce['OtherNoB']:.4f} | {ce['OtherB']:.4f} | {esw:.4f} | {et:.4f} | "
                       f"{100*s['MisMatch']:+.1f} % / {100*s['OtherB']:+.1f} % / {(esw-b)/b:+.1f} % / {(et-b)/b:+.1f} % |")
    open(f"{a.out}/b1_closure_{SPN}.md", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="closure"); ap.add_argument("--mode", default="smooth"); ap.add_argument("--kind", default="hgb")
    ap.add_argument("--species", default="D0")
    ap.add_argument("--leak", type=float, default=0.0); ap.add_argument("--sigma", type=float, default=0.035)
    ap.add_argument("--swap_factor", type=float, default=1.0); ap.add_argument("--strength", type=float, default=0.8); ap.add_argument("--mc_truth_min", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0); ap.add_argument("--out", default=os.path.join(HERE, "out"))
    a = ap.parse_args(); set_species(a.species); STRENGTH[0] = a.strength
    {"closure": closure, "report": report}[a.stage](a)


if __name__ == "__main__":
    main()
