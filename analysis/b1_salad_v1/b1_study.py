#!/usr/bin/env python
"""B1 — closure SALAD v1.0 (D⁺, podpróbki S6) gładki vs ADWERSARIALNY (skład niegładki w m) + test wycieku sygnału.

Pseudo-dane = MC tła (truth ≥ 8, parzystość 1) z zniekształceniem; „MC” = parzystość 0. Prawda = pseudo SR.
Tryby (--mode):
  smooth      : jak S8 (w = exp(Σ(c0+c1·t)·z_j), akceptacja-odrzucenie) — gładkie w m
  cat_global  : skład per truth_category, STAŁY w m (MisMatch* ×1,5, OtherNoB ×0,8, Other*B ×1,5) — SALAD powinien skorygować
  cat_local   : jak wyżej, ale tylko pod pikiem: w = 1 + (f_c−1)·g(m), g = exp(−(m−m0)²/2σ²), σ = 35 MeV — niewidoczne w SB
  cat_step    : skład z progiem: MisMatch* ×1,5 dla m < m0 (lewy SB + lewa połowa SR), ×1 dla m > m0 — SALAD widzi tylko lewy SB
  cat_local+smooth : cat_local nałożone na smooth
  cuts_local  : cięcia fizyczne lokalne: Lxy>0,3 ×1,5; |IPsig|>3 ×0,7 (profil g(m))
Metryki jak S4/S8: C2ST raw SB/SR, SALAD SR, SALAD pasmo, SB-kopia, podłoga; ε_B pred/true @WP80/70/50 (HGB fd vs pseudo-SB).
--leak f: do pseudo-danych SB dodaj fd (cechy sygnału, m ∼ U(70–90 MeV od m0, obie strony)) tak, by stanowiły f·N(70–90) —
          mierzy, czy SALAD „uczy się” sygnału (zmiana wag w SR, C2ST, ε_B).
"""
from __future__ import annotations
import argparse, json, os, sys, time, warnings
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ.get("CHARMTAG_SRC", "/work/b1/in/AA_merge/src")
sys.path.insert(0, SRC)
warnings.filterwarnings("ignore")
from charmtag.bkggen import salad_v1 as SV  # noqa: E402
from charmtag.bkggen.feta import c2st as _c2st  # noqa: E402

CFG = SV.load_cfg(); CFG["c2st"]["max_n"] = 60_000   # studium: szybciej (podłoga AUC ±0,003)
SPN = "Dplus"; SP = SV.species_params(CFG, SPN); M0 = SP["m0"]
H = list(CFG["features"]); PT_BINS = [3, 5, 8, 30]; WPS = (0.8, 0.7, 0.5)
SUB = os.path.join(HERE, "sub")


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


# --------------------------------------------------------------------------- zniekształcenia
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
    return np.where((t >= 8) & (t <= 11), f["MisMatch"], np.where(t == 12, f["OtherNoB"], np.where(t >= 13, f["OtherB"], 1.0)))


def g_local(m, sigma=0.035):
    return np.exp(-0.5 * ((np.asarray(m, float) - M0) / sigma) ** 2)


def w_mode(B, feats, mode, sigma=0.035):
    m, t = B["m"].to_numpy(float), B["truth"].to_numpy()
    if mode == "smooth":
        return w_smooth(B, feats)
    if mode == "cat_global":
        return cat_factor(t) / 1.5
    if mode == "cat_local":
        return (1 + (cat_factor(t) - 1) * g_local(m, sigma)) / 1.5
    if mode == "cat_step":
        f = np.where(((t >= 8) & (t <= 11)) & (m < M0), 1.5, 1.0)
        return f / 1.5
    if mode == "cat_local+smooth":
        return w_smooth(B, feats) * (1 + (cat_factor(t) - 1) * g_local(m, sigma)) / 1.5
    if mode == "cuts_local":
        f = np.where(B["Lxy"].to_numpy(float) > 0.3, 1.5, 1.0) * np.where(B["abs_impact_sig"].to_numpy(float) > 3, 0.7, 1.0)
        return (1 + (f - 1) * g_local(m, sigma)) / 1.5
    raise ValueError(mode)


def make_pseudo(B1, feats, mode, rng, sigma=0.035):
    w = w_mode(B1, feats, mode, sigma)
    keep = rng.random(len(B1)) < w
    return B1[keep].reset_index(drop=True)


def inject_signal(P, S_inj, frac, rng):
    """Dodaj fd (cechy) z m ∼ U w 70–90 MeV od m0 (obie strony), tak by stanowiły frac·N_pseudo(70–90)."""
    m = P["m"].to_numpy(float); edge = (np.abs(m - M0) >= 0.070) & (np.abs(m - M0) < 0.090)
    n = int(round(frac * edge.sum() / (1 - frac)))
    idx = rng.choice(len(S_inj), n, replace=n > len(S_inj))
    I = S_inj.iloc[idx].reset_index(drop=True)
    I["m"] = M0 + rng.choice([-1, 1], n) * rng.uniform(0.070, 0.090, n)
    I["truth"] = -1
    return pd.concat([P, I[P.columns]], ignore_index=True), n


# --------------------------------------------------------------------------- closure
def closure(a):
    t0 = time.time(); rng = np.random.default_rng(a.seed)
    S, B = load("sig"), load("mcbkg")
    B = B[B["truth"] >= 8].reset_index(drop=True)
    feats = [f for f in H if f in B and np.isfinite(B[f].to_numpy(float)).any()]
    A = B[B["parity"] == 0].reset_index(drop=True)
    P = make_pseudo(B[B["parity"] == 1].reset_index(drop=True), feats, a.mode, rng, a.sigma)
    tag = f"{a.mode}_{a.kind}" + (f"_leak{a.leak:g}" if a.leak else "")
    res = {"mode": a.mode, "kind": a.kind, "leak": a.leak, "sigma": a.sigma, "n_A": int(len(A)), "n_pseudo": int(len(P))}
    S0, S1 = S[S.parity == 0].reset_index(drop=True), S[S.parity == 1].reset_index(drop=True)
    n_inj = 0
    P_fit = P
    if a.leak:
        P_fit, n_inj = inject_signal(P, S1, a.leak, rng)       # wstrzyknięty sygnał = fd parzystość 1 (nie w klasyfikatorze ε)
    res["n_injected"] = int(n_inj)
    rA, rP = SV.regions(A["m"], SP), SV.regions(P["m"], SP)
    rF = SV.regions(P_fit["m"], SP)
    XA, XP = A[feats].to_numpy(float), P[feats].to_numpy(float)
    res["n_pseudo_SR"] = int((rP == "SR").sum())
    res["frac_MisMatch_SR_pseudo"] = float(((P["truth"] >= 8) & (P["truth"] <= 11))[rP == "SR"].mean())
    res["frac_MisMatch_SB_pseudo"] = float(((P["truth"] >= 8) & (P["truth"] <= 11))[rP == "SB"].mean())
    res["frac_MisMatch_SR_A"] = float(((A["truth"] >= 8) & (A["truth"] <= 11))[rA == "SR"].mean())
    res["c2st_raw_SB"] = c2st(XP[rP == "SB"], XA[rA == "SB"], seed=a.seed)
    res["c2st_raw_SR"] = c2st(XP[rP == "SR"], XA[rA == "SR"], seed=a.seed)
    print(tag, "raw", {k: v for k, v in res.items() if "c2st" in k or "frac" in k}, f"{time.time()-t0:.0f}s", flush=True)
    # SALAD v1: dane = pseudo (SB + pasmo), MC = A (cały zakres); fit tylko na SB
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
    # ε_B: klasyfikator fd(p0) vs pseudo-SB (jak dane SB); WP z fd(p1)
    Xs = S0[feats].to_numpy(float); Xb = XP[rP == "SB"]
    if len(Xb) > 250_000: Xb = Xb[rng.choice(len(Xb), 250_000, replace=False)]
    clf = hgb(np.vstack([Xs, Xb]), np.r_[np.ones(len(Xs)), np.zeros(len(Xb))], seed=a.seed)
    sc = lambda X: clf.predict_proba(X)[:, 1]  # noqa: E731
    thr = wp_thresholds(sc(S1[feats].to_numpy(float)), S1["pt"].to_numpy(float))
    sA, sP = sc(XA), sc(XP); ptA, ptP = A["pt"].to_numpy(float), P["pt"].to_numpy(float)
    for wp in WPS:
        tru = eff(sP[rP == "SR"], ptP[rP == "SR"], thr, None, wp)
        res[f"eff_true_SR@{wp}"] = tru
        res[f"eff_salad_SR@{wp}"] = eff(sA[rA == "SR"], ptA[rA == "SR"], thr, wA[rA == "SR"], wp)
        res[f"eff_mcraw_SR@{wp}"] = eff(sA[rA == "SR"], ptA[rA == "SR"], thr, None, wp)
        res[f"eff_sbcopy@{wp}"] = eff(sP[rP == "SB"], ptP[rP == "SB"], thr, None, wp)
        res[f"eff_salad_SB@{wp}"] = eff(sA[rA == "SB"], ptA[rA == "SB"], thr, wA[rA == "SB"], wp)
        res[f"eff_true_SB@{wp}"] = eff(sP[rP == "SB"], ptP[rP == "SB"], thr, None, wp)
        n = (rP == "SR").sum(); res[f"eff_true_relerr@{wp}"] = float(np.sqrt((1 - tru["all"]) / max(tru["all"] * n, 1)))
        # błąd predykcji SALAD (ważony, n_eff)
        k = sA[rA == "SR"] >= 0; ws = wA[rA == "SR"]; p = res[f"eff_salad_SR@{wp}"]["all"]
        neff = ws.sum() ** 2 / (ws ** 2).sum(); res[f"eff_salad_relerr@{wp}"] = float(np.sqrt((1 - p) / max(p * neff, 1)))
    # wagi na MC „sygnałopodobnym” (score ≥ próg WP80) — test wycieku: czy SALAD podbija takie zdarzenia w SR
    sl = np.zeros(len(A), bool)
    for (a_, b_), t in thr.items():
        sl |= (ptA >= a_) & (ptA < b_) & (sA >= t[0.8])
    res["mean_w_siglike_SR"] = float(wA[(rA == "SR") & sl].mean()); res["mean_w_siglike_SB"] = float(wA[(rA == "SB") & sl].mean())
    res["mean_w_siglike_band"] = float(wA[(rA == "band") & sl].mean())
    print(tag, {k: round(v["all"], 4) for k, v in res.items() if k.startswith("eff_") and isinstance(v, dict)}, f"{time.time()-t0:.0f}s", flush=True)
    os.makedirs(a.out, exist_ok=True)
    json.dump(res, open(f"{a.out}/closure_{tag}.json", "w"), indent=1)
    if a.leak:   # zapisz też wagi SR do porównania z biegiem bez wycieku
        np.savez_compressed(f"{a.out}/wSR_{tag}.npz", w=wA[rA == "SR"].astype(np.float32))
    elif a.mode == "smooth":
        np.savez_compressed(f"{a.out}/wSR_{tag}.npz", w=wA[rA == "SR"].astype(np.float32))


def report(a):
    import glob
    files = sorted(glob.glob(f"{a.out}/closure_*.json"))
    R = [json.load(open(f)) for f in files]
    out = ["# B1 — closure SALAD v1.0 (D⁺, podpróbki S6; pseudo-dane = MC B parzystość 1 zniekształcone; prawda = pseudo SR)\n"]
    out.append("| tryb | kind | wyciek | n_pseudo SR | f_MM SR pseudo / A | C2ST raw SB | raw SR | SALAD SB | **SALAD SR** | SALAD pasmo | SB-kopia SR | podłoga | n_eff/N | AUC_SB | ε_B pred/true @80 (SALAD / SB-kopia / MC raw) | @70 | @50 | σ_stat prawdy @80/70/50 | w̄ sygn.-pod. SR/pasmo/SB |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in R:
        cells = []
        for wp in WPS:
            t = r[f"eff_true_SR@{wp}"]["all"]
            cells.append(" / ".join(f"{r[f'eff_{s}@{wp}']['all'] / t:.3f}" for s in ("salad_SR", "sbcopy", "mcraw_SR")))
        se = " / ".join(f"{r[f'eff_true_relerr@{wp}']:.3f}" for wp in WPS)
        out.append(f"| {r['mode']} | {r['kind']} | {r['leak'] or 0:g} | {r['n_pseudo_SR']} | {r['frac_MisMatch_SR_pseudo']:.3f} / {r['frac_MisMatch_SR_A']:.3f} | "
                   f"{r['c2st_raw_SB']:.3f} | {r['c2st_raw_SR']:.3f} | {r['c2st_salad_SB']:.3f} | **{r['c2st_salad_SR']:.3f}** | {r['c2st_salad_band']:.3f} | "
                   f"{r['c2st_sbcopy_SR']:.3f} | {r['c2st_floor_SR']:.3f} | {r['neff_over_N_SR']:.2f} | {r['val_auc_SB']:.3f} | "
                   f"{cells[0]} | {cells[1]} | {cells[2]} | {se} | {r['mean_w_siglike_SR']:.3f} / {r['mean_w_siglike_band']:.3f} / {r['mean_w_siglike_SB']:.3f} |")
    out.append("\n## Per bin pT: ε_B(SR) SALAD/prawda @80 / @70\n")
    out.append("| tryb | kind | wyciek | 3-5 | 5-8 | 8-30 | all |"); out.append("|---|---|---|---|---|---|---|")
    for r in R:
        cells = [" / ".join(f"{r[f'eff_salad_SR@{wp}'][b] / r[f'eff_true_SR@{wp}'][b]:.3f}" for wp in (0.8, 0.7)) for b in ("3-5", "5-8", "8-30", "all")]
        out.append(f"| {r['mode']} | {r['kind']} | {r['leak'] or 0:g} | " + " | ".join(cells) + " |")
    open(f"{a.out}/b1_closure_Dplus.md", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="closure"); ap.add_argument("--mode", default="smooth"); ap.add_argument("--kind", default="hgb")
    ap.add_argument("--leak", type=float, default=0.0); ap.add_argument("--sigma", type=float, default=0.035)
    ap.add_argument("--seed", type=int, default=0); ap.add_argument("--out", default=os.path.join(HERE, "out"))
    a = ap.parse_args()
    {"closure": closure, "report": report}[a.stage](a)


if __name__ == "__main__":
    main()
