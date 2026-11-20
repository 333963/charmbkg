#!/usr/bin/env python3
"""Krok 2: SF klasyfikatora H (tag D*) per klasa trygera, WP95/80/50, biny 5–8, 8–30 (+3–5 informacyjnie).
Definicje jak btrig_analyze.cmd_sf: |m(Kπ)−m0| < 55 MeV, progi = kwantyle score prawdziwych D0 w fd_nd (all_*), ε_MC bez wag.
SBS zawsze; fit Δm pass/fail (profil) gdy S_SBS ≥ 120. Wynik: out/s2.json."""
from multiprocessing import Pool
import numpy as np
import afp_common as A
from iso_common import load_pairs, sbs_w, eff_w, fit_pf, C, DR
from s1_iso import samples
import bkg_common as K

WPS = ["0.95", "0.80", "0.50"]
BINS = {"3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}


def _fit(job):
    key, dp, df_, shape = job
    return key, fit_pf(dp, df_, shape, interval=True)


def main():
    thr_all = K.thresholds()
    thr = {b: {wp: thr_all[b][K.WPS.index(wp)] for wp in WPS} for b in BINS}
    M = load_pairs("bt_fd_nd")
    f_sig = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < C.NOM_WIN)
    D8, D9, S = samples()
    res = {"thr": thr, "rows": {}}
    jobs = []
    for b, (lo, hi) in BINS.items():
        mf = f_sig & (M["pt"] >= lo) & (M["pt"] < hi)
        shape = DR.shape_from_mc(M["dm"][mf] * 1e3)
        for wp in WPS:
            t = thr[b][wp]
            emc = float((M["score"][mf] >= t).mean()); emc_e = float(np.sqrt(emc * (1 - emc) / mf.sum()))
            for name, (D, cm) in S.items():
                md = cm & (np.abs(D["m"] - C.M0_D0) < C.NOM_WIN) & (D["pt"] >= lo) & (D["pt"] < hi)
                r = {"ncand": int(md.sum()), "eps_mc": emc}
                if md.sum() > 20:
                    pas = D["score"] >= t
                    w = sbs_w(D["dm"][md], D["pt"][md], pas[md])
                    e, er, T, eT = eff_w(w, pas[md].astype(float))
                    r.update({"eps": e, "err": er, "S": T, "S_err": eT, "SF": e / emc if np.isfinite(e) else np.nan,
                              "SF_err": float(np.hypot(er / emc, e / emc * emc_e / emc)) if np.isfinite(e) else np.nan})
                    if T >= 120:
                        jobs.append(((name, b, wp), D["dm"][md & pas], D["dm"][md & ~pas], shape))
                res["rows"][f"{name}|{b}|{wp}"] = r
    print("fity:", len(jobs), flush=True)
    with Pool(4) as P:
        for (name, b, wp), f in P.imap_unordered(_fit, jobs):
            r = res["rows"][f"{name}|{b}|{wp}"]
            r["fit"] = f
            if np.isfinite(f["eps"]):
                r["SF_fit"] = f["eps"] / r["eps_mc"]; r["SF_fit_err"] = f["err"] / r["eps_mc"]
    A.jdump(res, f"{A.OUT}/s2.json")
    for k, r in res["rows"].items():
        if "S" in r and k.split("|")[0] in ("8:AFPstream", "8:AFPonly", "8:AFPsptrk4", "8:AFPnoalg_any", "8:J12only!A", "8:MBTS!A", "8:none", "8:other", "9:MB427929", "8:all"):
            s = f"{k:32s} n {r['ncand']:6d} S {r['S']:6.0f}±{r['S_err']:3.0f} eps {r['eps']:.3f}±{r['err']:.3f} mc {r['eps_mc']:.3f} SF {r['SF']:.2f}±{r['SF_err']:.2f}"
            if "fit" in r:
                s += f" | FIT N {r['fit']['N']:.0f} SF {r.get('SF_fit', float('nan')):.2f}±{r.get('SF_fit_err', float('nan')):.2f} det {r['fit']['det']}"
            print(s)


if __name__ == "__main__":
    main()
