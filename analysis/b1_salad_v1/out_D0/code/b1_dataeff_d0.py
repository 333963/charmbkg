#!/usr/bin/env python
"""B1-D0: kontrola na DANYCH (jak S8 `data_eff` / B1 §2 pkt 4): ε_B(SR) z fitu pol2 SB→okno po cięciu WP vs predykcja SALAD v1.0.
Klasyfikator H: fd truth 0/1 (parz. 0) vs dane SB (data p0 ∪ p1, parz. 0); progi WP z fd parz. 1 (biny pT); ocena na danych parz. 1.
ε_B(dane, SR) = N_pol2(SR | po cięciu) / N_pol2(SR | przed cięciem), pol2 dopasowany do |m−m0| > 75 MeV (SB) — sygnał (S/B ≈ 0,6 %) poza fitem.
SALAD: fit v1.0 (hgb i sk_mlp) na danych parz. 1 SB × MC mb truth ≥ 8 → ε_B(SR) ważone. MC raw i SB-kopia (dane parz. 1 SB) dla porównania."""
import os, sys, json, time, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.environ.get("CHARMTAG_SRC", "/work/b1d0/in/src")); warnings.filterwarnings("ignore")
from charmtag.bkggen import salad_v1 as SV
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b1_study_d0 as ST
ST.set_species("D0"); SP, M0 = ST.SP, ST.M0; CFG = ST.CFG
load = ST.load
t0 = time.time(); rng = np.random.default_rng(0)
S, B, D = load("sig"), load("mcbkg"), load("data")
feats = [f for f in ST.H if f in B and np.isfinite(B[f].to_numpy(float)).any()]
MCMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 8
B = B[(B.truth >= MCMIN) & ((B.truth < 100) if MCMIN >= 8 else True)].reset_index(drop=True); CFG["mc_truth_min"] = MCMIN
S0, S1 = S[(S.parity == 0) & S.truth.isin([0, 1])], S[(S.parity == 1) & S.truth.isin([0, 1])]
D0, D1 = D[D.parity == 0].reset_index(drop=True), D[D.parity == 1].reset_index(drop=True)
r0, r1, rB = SV.regions(D0.m, SP), SV.regions(D1.m, SP), SV.regions(B.m, SP)
Xb = D0.loc[r0 == "SB", feats].to_numpy(float)
if len(Xb) > 250_000: Xb = Xb[rng.choice(len(Xb), 250_000, replace=False)]
clf = ST.hgb(np.vstack([S0[feats].to_numpy(float), Xb]), np.r_[np.ones(len(S0)), np.zeros(len(Xb))])
sc = lambda X: clf.predict_proba(X[feats].to_numpy(float))[:, 1]  # noqa: E731
thr = ST.wp_thresholds(sc(S1), S1.pt.to_numpy(float))
s1, sB = sc(D1), sc(B); pt1, ptB = D1.pt.to_numpy(float), B.pt.to_numpy(float); m1 = D1.m.to_numpy(float)
print("clf ok", f"{time.time()-t0:.0f}s", flush=True)


def passed(s, pt, wp):
    k = np.zeros(len(s), bool)
    for (a, b), t in thr.items():
        k |= (pt >= a) & (pt < b) & (s >= t[wp])
    return k


def pol2_sr(m, nb=50):
    """Fit pol2 (ważony Poissonem) do histogramu w SB (|m−m0|>sb_min) → całka w SR i błąd (z kowariancji)."""
    lo, hi = SP["mass_range"]; h, e = np.histogram(m, bins=nb, range=(lo, hi)); c = 0.5 * (e[1:] + e[:-1]); x = c - M0
    sb = np.abs(x) > SP["sb_min"]; sr = np.abs(x) < SP["sr_half"]
    V = np.vander(x[sb], 3); w = 1 / np.sqrt(np.maximum(h[sb], 1))
    p, res, *_ = np.linalg.lstsq(V * w[:, None], h[sb] * w, rcond=None)
    cov = np.linalg.inv((V * w[:, None]).T @ (V * w[:, None]))
    Vsr = np.vander(x[sr], 3); n = float((Vsr @ p).sum()); err = float(np.sqrt(Vsr.sum(0) @ cov @ Vsr.sum(0)))
    chi2 = float((((V @ p - h[sb]) * w) ** 2).sum()); return n, err, chi2 / max(sb.sum() - 3, 1), int(h[sr].sum())


out = {"mc_truth_min": MCMIN, "n_data_p1": int(len(D1)), "n_data_p1_SR": int((r1 == "SR").sum()), "feats": feats}
N0, e0, chi0, raw0 = pol2_sr(m1); out["pol2_pre"] = dict(N_SR=N0, err=e0, chi2ndf=chi0, N_SR_raw=raw0)
sal = {}
for kind in ("hgb", "sk_mlp"):
    rw = SV.fit_salad_v1(CFG, "D0", D1[r1 != "SR"], B, kind=kind, seed=0, c2st_band=True)
    sal[kind] = (SV.salad_weights(rw, B), {k: rw.diag[k] for k in ("val_auc_SB", "neff_over_N_SR", "c2st_band")})
    print(kind, sal[kind][1], f"{time.time()-t0:.0f}s", flush=True)
for wp in ST.WPS:
    k1 = passed(s1, pt1, wp); kB = passed(sB, ptB, wp)
    N, e, chi, raw = pol2_sr(m1[k1])
    eff_data = N / N0; err = eff_data * np.sqrt((e / N)**2 + (e0 / N0)**2)
    row = dict(eff_data_SR=eff_data, err=err, chi2ndf_post=chi, N_post_SR_fit=N, N_post_SR_raw=raw,
               eff_sbcopy=float(k1[r1 == "SB"].mean()), eff_mcraw_SR=float(kB[rB == "SR"].mean()))
    for kind, (w, _) in sal.items():
        row[f"eff_salad_{kind}_SR"] = float((w[rB == "SR"] * kB[rB == "SR"]).sum() / w[rB == "SR"].sum())
    row["ratio_salad_hgb"] = row["eff_salad_hgb_SR"] / eff_data; row["ratio_salad_sk_mlp"] = row["eff_salad_sk_mlp_SR"] / eff_data
    row["ratio_sbcopy"] = row["eff_sbcopy"] / eff_data; row["ratio_mcraw"] = row["eff_mcraw_SR"] / eff_data
    out[f"wp{wp}"] = row
    print(wp, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
out["salad_diag"] = {k: v[1] for k, v in sal.items()}
json.dump(out, open(f"out/data_eff_D0{'' if MCMIN == 8 else f'_mc{MCMIN}'}.json", "w"), indent=1)
print("done", f"{time.time()-t0:.0f}s")
