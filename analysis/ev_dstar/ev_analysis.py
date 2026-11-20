#!/usr/bin/env python3
"""EV analiza: anomalia 5–8/WP95 (B2) na zdarzeniowych parach D0←D* + SF_rel per bin.
Wejście: katalog z ev_data.npz, ev_fd.npz [, ev_mb.npz] (ev_extract.py / merge z lxplus).
Pytanie B2: w danych ε(WP95) w oknie ±55 MeV jest ~0,78·MC, ale ≈MC w oknie ±36 → co siedzi w pierścieniu 36–55 MeV?
Metoda (per bin pT D0: 3–5, 5–8, 8–30; regiony masy D0: core |m−m0|<36, ring 36–55, sb 55–75 MeV):
  1. ε_dane = N_pass/(N_pass+N_fail) z JEDNOCZESNEGO fitu Δm pass/fail (dstar_refit, kształt sygnału z MC fd truth per bin)
     dla WP 95/90/80/70/50 i każdego regionu; ε_MC z truth (fd, D0 i D* truth 0/1, swapy osobno). SF = ε_dane/ε_MC.
  2. Udział sygnału Δm w pierścieniu vs rdzeniu (N_sig ring / N_sig core) — dane vs MC → czy ogony masy D0 w danych są większe.
  3. Skład pierścienia w MC fd (truth): prawdziwe D0 (ogony), swapy K↔π (truth 100/101), tło — i ich ε(WP95).
  4. Pierścień, dane: rozkłady score, m_swap, pT(π_s), Lxy, pointing, |cosθ*| dla sygnału Δm (SR − α·SB) vs MC truth,
     osobno pass/fail WP95 → gdzie H odrzuca.
Wyjście: ev_analysis.json + ev_analysis.md. Użycie: ev_analysis.py --dir DIR [--label podprobka|full]
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import dstar_refit as DR  # kopia Q_lxplus/v_fix/py/dstar_refit.py (etap V)

BINS = {"3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
REG = {"core": (0.0, 0.036), "ring": (0.036, 0.055), "nom": (0.0, 0.055), "sb": (0.055, 0.075)}
WPS = ["0.95", "0.90", "0.80", "0.70", "0.50"]
E_DM = np.arange(139.6, 170.0001, 0.1)       # MeV
SR_DM, SB_DM = (143.5, 147.5), (150.0, 165.0)
FEATS_SHOW = ["score", "m_swap_dist", "pis_pt", "Lxy", "pointing", "abs_cos_theta_star", "impact_sig", "chi2", "d_min_abs_d0", "dedx_K"]


def load(path):
    z = np.load(path, allow_pickle=False)
    d = {k: z[k] for k in z.files if k != "meta"}
    d["meta"] = json.loads(str(z["meta"]))
    return d


def thresholds(fd, m0):
    """Progi WP per bin z fd truth 0/1 (|m−m0|<55): 0,95/0,90 liczone; 0,80/0,70/0,50 z meta (kontrola)."""
    thr = {}
    sc, pt, tr, m = fd["all_score"], fd["all_pt"], fd["all_truth"], fd["all_m"]
    sig = np.isin(tr, [0, 1]) & (np.abs(m - m0) < 0.055)
    meta_thr = fd["meta"].get("thresholds") or {}
    for b, (lo, hi) in BINS.items():
        s = sc[sig & (pt >= lo) & (pt < hi)]
        thr[b] = {wp: float(np.quantile(s, 1 - float(wp))) for wp in WPS}
        for wp, v in (meta_thr.get(b) or {}).items():
            thr[b][f"meta_{wp}"] = float(v)
    return thr


def fit_eps(dm_pass, dm_fail, shape):
    hp, _ = np.histogram(dm_pass * 1e3, E_DM)
    hf, _ = np.histogram(dm_fail * 1e3, E_DM)
    if hp.sum() + hf.sum() < 50:
        return {"eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False}
    try:
        r = DR.fit_passfail(hp, hf, E_DM, shape)
    except Exception as ex:  # noqa: BLE001
        return {"eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False, "err": str(ex)}
    return {"eps": float(r.get("eps", np.nan)), "eps_err": float(r.get("eps_err", np.nan)), "N": float(r.get("N_tot", np.nan)),
            "N_err": float(r.get("N_tot_err", np.nan)), "valid": bool(r.get("valid", False)), "chi2": float(r.get("chi2", np.nan)),
            "ndf": int(r.get("ndf", 0))}


def fit_n(dm, shape):
    h, _ = np.histogram(dm * 1e3, E_DM)
    if h.sum() < 50:
        return {"N": 0.0, "N_err": 0.0, "valid": False}
    try:
        return DR.fit_single(h, E_DM, shape)
    except Exception as ex:  # noqa: BLE001
        return {"N": 0.0, "N_err": 0.0, "valid": False, "err": str(ex)}


def alpha_sb(dm):
    """α = tło(SR)/tło(SB) z fitu samego tła Δm poza pikiem (do odejmowania SB w rozkładach cech)."""
    h, _ = np.histogram(dm * 1e3, E_DM)
    h, e = DR._rebin(h, E_DM)
    c = 0.5 * (e[1:] + e[:-1])
    side = (c < 142.5) | (c > 148.5)
    from iminuit import Minuit
    def nll(nb, a, b, c1, c2):
        lam = np.clip(nb * np.diff(DR._bkg_cdf(e, a, b, c1, c2)), 1e-9, None)[side]
        return float(np.sum(lam - h[side] * np.log(lam)))
    m = Minuit(nll, nb=h.sum(), a=0.5, b=0.1, c1=0.0, c2=0.0); m.errordef = 0.5
    m.limits["a"] = (0.01, 3); m.limits["b"] = (-1, 3); m.limits["c1"] = (-0.5, 2); m.limits["c2"] = (-0.05, 0.2); m.limits["nb"] = (0, None)
    m.migrad()
    lam = m.values["nb"] * np.diff(DR._bkg_cdf(e, *[m.values[k] for k in ("a", "b", "c1", "c2")]))
    B = lambda lo, hi: float(lam[(c > lo) & (c < hi)].sum())
    return B(*SR_DM) / max(B(*SB_DM), 1e-9), bool(m.valid)


def wmean(x, w):
    s = w.sum()
    return float((x * w).sum() / s) if s > 0 else np.nan


def feat_compare(d, mask_d, f, mask_f, alpha, feats):
    """Dane: SR − α·SB (Δm) vs MC truth sygnał (SR). Zwraca średnie i udziały pass WP95."""
    dm_d = d["dm"] * 1e3
    sr = mask_d & (dm_d > SR_DM[0]) & (dm_d < SR_DM[1])
    sb = mask_d & (dm_d > SB_DM[0]) & (dm_d < SB_DM[1])
    out = {"n_sr": int(sr.sum()), "n_sb": int(sb.sum()), "alpha": alpha, "S_est": float(sr.sum() - alpha * sb.sum())}
    for c in feats:
        xd = d[c]; xf = f[c]
        wd = np.where(sr, 1.0, 0.0) - alpha * np.where(sb, 1.0, 0.0)
        out[c] = {"data_sw": wmean(xd, wd), "mc": float(xf[mask_f].mean()) if mask_f.any() else np.nan,
                  "mc_sd": float(xf[mask_f].std()) if mask_f.sum() > 1 else np.nan}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    d = load(os.path.join(a.dir, "ev_data.npz")); f = load(os.path.join(a.dir, "ev_fd.npz"))
    m0 = float(f["meta"]["m0"])
    for x in (d, f):
        x["m_swap_dist"] = np.abs(x["m_swap"] - m0)
    thr = thresholds(f, m0)
    res = {"label": a.label or a.dir, "m0": m0, "thresholds": thr, "n_pairs": {"data": int(len(d["m"])), "fd": int(len(f["m"]))}, "bins": {}}
    f_sig = np.isin(f["truth"], [0, 1]) & np.isin(f["truth_ds"], [0, 1])
    f_swap = np.isin(f["truth"], [100, 101]) & np.isin(f["truth_ds"], [0, 1])
    for b, (lo, hi) in BINS.items():
        bd = (d["pt"] >= lo) & (d["pt"] < hi); bf = (f["pt"] >= lo) & (f["pt"] < hi)
        shape = DR.shape_from_mc(f["dm"][bf & f_sig & (np.abs(f["m"] - m0) < 0.055)] * 1e3)
        R = {"shape_mc": shape, "regions": {}}
        for reg, (r0, r1) in REG.items():
            md = bd & (np.abs(d["m"] - m0) >= r0) & (np.abs(d["m"] - m0) < r1)
            mf = bf & (np.abs(f["m"] - m0) >= r0) & (np.abs(f["m"] - m0) < r1)
            rr = {"n_data": int(md.sum()), "n_mc_sig": int((mf & f_sig).sum()), "n_mc_swap": int((mf & f_swap).sum()),
                  "n_mc_bkg": int((mf & ~f_sig & ~f_swap).sum()), "N_sig_data": fit_n(d["dm"][md], shape), "wp": {}}
            for wp in WPS:
                t = thr[b][wp]
                pas = d["score"] >= t
                fit = fit_eps(d["dm"][md & pas], d["dm"][md & ~pas], shape)
                eps_mc = float((f["score"][mf & f_sig] >= t).mean()) if (mf & f_sig).any() else np.nan
                eps_mc_swap = float((f["score"][mf & f_swap] >= t).mean()) if (mf & f_swap).any() else np.nan
                rr["wp"][wp] = {"eps_data": fit["eps"], "eps_data_err": fit["eps_err"], "valid": fit["valid"], "eps_mc": eps_mc,
                                "eps_mc_swap": eps_mc_swap,
                                "SF": fit["eps"] / eps_mc if eps_mc and np.isfinite(fit["eps"]) else np.nan,
                                "SF_err": fit["eps_err"] / eps_mc if eps_mc and np.isfinite(fit["eps_err"]) else np.nan}
            R["regions"][reg] = rr
        # udział sygnału Δm w pierścieniu vs rdzeń: dane (fit) vs MC (truth)
        Nc, Nr = R["regions"]["core"]["N_sig_data"], R["regions"]["ring"]["N_sig_data"]
        nc_mc, nr_mc = R["regions"]["core"]["n_mc_sig"], R["regions"]["ring"]["n_mc_sig"]
        R["ring_over_core"] = {"data": Nr["N"] / Nc["N"] if Nc["N"] > 0 else np.nan,
                               "data_err": (Nr["N"] / Nc["N"]) * np.hypot(Nr["N_err"] / max(Nr["N"], 1e-9), Nc["N_err"] / max(Nc["N"], 1e-9)) if Nc["N"] > 0 and Nr["N"] > 0 else np.nan,
                               "mc_true": nr_mc / nc_mc if nc_mc else np.nan,
                               "mc_true_plus_swap": (nr_mc + R["regions"]["ring"]["n_mc_swap"]) / (nc_mc + R["regions"]["core"]["n_mc_swap"]) if nc_mc else np.nan}
        # cechy w pierścieniu: dane SR−α·SB vs MC truth; pass/fail WP95
        mring_d = bd & (np.abs(d["m"] - m0) >= 0.036) & (np.abs(d["m"] - m0) < 0.055)
        mring_f = bf & (np.abs(f["m"] - m0) >= 0.036) & (np.abs(f["m"] - m0) < 0.055) & f_sig
        alpha, okA = alpha_sb(d["dm"][mring_d])
        t95 = thr[b]["0.95"]
        R["ring_features"] = {"alpha_valid": okA,
                              "all": feat_compare(d, mring_d, f, mring_f, alpha, FEATS_SHOW),
                              "pass95": feat_compare(d, mring_d & (d["score"] >= t95), f, mring_f & (f["score"] >= t95), alpha, FEATS_SHOW),
                              "fail95": feat_compare(d, mring_d & (d["score"] < t95), f, mring_f & (f["score"] < t95), alpha, FEATS_SHOW)}
        # skład pierścienia w MC po WP95: czego H nie przyjmuje
        R["ring_mc_fail95"] = {"true_D0": int((mring_f & (f["score"] < t95)).sum()),
                               "swap": int((bf & (np.abs(f["m"] - m0) >= 0.036) & (np.abs(f["m"] - m0) < 0.055) & f_swap & (f["score"] < t95)).sum()),
                               "true_D0_pass": int((mring_f & (f["score"] >= t95)).sum())}
        res["bins"][b] = R
        print(b, "shape", {k: round(v, 3) for k, v in shape.items()}, "ring/core data %.2f±%.2f MC %.2f" % (
            R["ring_over_core"]["data"], R["ring_over_core"]["data_err"], R["ring_over_core"]["mc_true"]))
    out = a.out or os.path.join(a.dir, "ev_analysis")
    json.dump(res, open(out + ".json", "w"), indent=1, default=float)
    # markdown
    L = [f"# EV — anomalia 5–8/WP95 i SF z tagu D* ({res['label']})", "",
         f"Pary D0←D*: dane {res['n_pairs']['data']}, fd {res['n_pairs']['fd']}. Regiony |m−m0|: core <36, ring 36–55, nom <55, sb 55–75 MeV. "
         "ε_dane z jednoczesnego fitu Δm pass/fail (dstar_refit), ε_MC z truth fd.", ""]
    L += ["## ε(WP95) i SF per region", "", "| bin | region | N_sig Δm (dane) | n_MC sig/swap/bkg | ε_dane 95 | ε_MC 95 | SF 95 | SF 80 | SF 70 | SF 50 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for b, R in res["bins"].items():
        for reg, rr in R["regions"].items():
            w95 = rr["wp"]["0.95"]
            def sf(wp):
                x = rr["wp"][wp]; return "—" if not np.isfinite(x["SF"]) else f"{x['SF']:.2f}±{x['SF_err']:.2f}"
            L.append(f"| {b} | {reg} | {rr['N_sig_data']['N']:.0f}±{rr['N_sig_data']['N_err']:.0f} | {rr['n_mc_sig']}/{rr['n_mc_swap']}/{rr['n_mc_bkg']} | "
                     f"{w95['eps_data']:.2f}±{w95['eps_data_err']:.2f} | {w95['eps_mc']:.2f} | {sf('0.95')} | {sf('0.80')} | {sf('0.70')} | {sf('0.50')} |")
    L += ["", "## Udział sygnału Δm w pierścieniu względem rdzenia (ogony masy D0)", "", "| bin | dane N_ring/N_core | MC prawdziwe D0 | MC D0+swapy |", "|---|---|---|---|"]
    for b, R in res["bins"].items():
        r = R["ring_over_core"]
        L.append(f"| {b} | {r['data']:.3f}±{r['data_err']:.3f} | {r['mc_true']:.3f} | {r['mc_true_plus_swap']:.3f} |")
    L += ["", "## Pierścień 36–55 MeV: cechy (dane SR−α·SB vs MC truth), wszystkie / pass95 / fail95", ""]
    for b, R in res["bins"].items():
        rf = R["ring_features"]
        L += [f"### {b}  (α={rf['all']['alpha']:.3f}, S_est all/pass/fail = {rf['all']['S_est']:.0f}/{rf['pass95']['S_est']:.0f}/{rf['fail95']['S_est']:.0f}; "
              f"MC fail95: prawdziwe D0 {R['ring_mc_fail95']['true_D0']}, swapy {R['ring_mc_fail95']['swap']}, pass {R['ring_mc_fail95']['true_D0_pass']})", "",
              "| cecha | dane all | MC all (σ) | dane fail95 | MC fail95 | dane pass95 | MC pass95 |", "|---|---|---|---|---|---|---|"]
        for c in FEATS_SHOW:
            g = lambda k: rf[k][c]  # noqa: E731
            L.append(f"| {c} | {g('all')['data_sw']:.3f} | {g('all')['mc']:.3f} ({g('all')['mc_sd']:.3f}) | {g('fail95')['data_sw']:.3f} | {g('fail95')['mc']:.3f} | {g('pass95')['data_sw']:.3f} | {g('pass95')['mc']:.3f} |")
        L.append("")
    open(out + ".md", "w").write("\n".join(L))
    print("zapisano", out + ".{json,md}")


if __name__ == "__main__":
    main()
