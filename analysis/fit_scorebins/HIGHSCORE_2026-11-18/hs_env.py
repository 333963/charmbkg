"""HIGHSCORE krok 3b/4 — otoczenie kandydata (ptcone40) i skład MC. EKSPLORACYJNE, PO FAKCIE.
(a) MC: ułamek > WP25 / > WP50 i udział nieprompt w funkcji ptcone40 przy ustalonym pT;
(b) test: MC przeważone w pT (drobne biny) i/lub ptcone40 (rozkład sygnału danych z regionu score >= WP70, odjęte pasma) —
    przewidywane ułamki 6 kategorii vs wydajki danych (chi2 z kowariancją bootstrapową wydajków; wagi ustalone z nominału).
python hs_env.py [NBOOT] -> out/env.json"""
import json, os, sys
import numpy as np
from scipy.stats import chi2 as CHI2
import hs_common as C
R = C.R
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 200
out = {"mc_table": [], "tests": []}
M = C.load("fd_nd"); mu = C.MU["8-30"]; W = 3.0
thr = C.thr_for("8-30", R.WP_MAIN)
sm = C.base_mask(M, "8-30", isdata=False) & (np.abs(M["dm"] - mu) < W)
mpt, msc, mcone, mnp = M["pt"][sm], M["score"][sm], M["ptcone40"][sm], M["np"][sm]
mcat = R.cat_index(msc, thr)
print("MC: f(>WP25), f(>WP50), udział nieprompt vs ptcone40 [GeV] przy ustalonym pT (wszystkie | tylko prompt)")
for lo, hi in ((8, 10), (10, 12), (12, 15), (15, 30)):
    for clo, chi in ((-1, 2), (2, 6), (6, 15), (15, 1e9)):
        s = (mpt >= lo) & (mpt < hi) & (mcone >= clo) & (mcone < chi); p = s & ~mnp
        if s.sum() < 30:
            continue
        row = {"pt": [lo, hi], "cone": [clo, min(chi, 999)], "n": int(s.sum()), "f25": float((mcat[s] == 5).mean()), "f50": float((mcat[s] >= 4).mean()),
               "f25_err": float(np.sqrt((mcat[s] == 5).mean() * (1 - (mcat[s] == 5).mean()) / s.sum())),
               "np": float(mnp[s].mean()), "f25_prompt": float((mcat[p] == 5).mean()), "f50_prompt": float((mcat[p] >= 4).mean())}
        out["mc_table"].append(row)
        print(f"  pT {lo:2d}-{hi:2d} ptcone40 {clo:4.0f}-{min(chi,999):4.0f}: n={row['n']:5d} f25={row['f25']:.3f}±{row['f25_err']:.3f} f50={row['f50']:.3f} np={row['np']:.3f} | prompt: f25={row['f25_prompt']:.3f} f50={row['f50_prompt']:.3f}")
print("MC: f(>WP25) prompt vs nieprompt:", round(float((mcat[~mnp] == 5).mean()), 3), round(float((mcat[mnp] == 5).mean()), 3),
      " f(>WP50):", round(float((mcat[~mnp] >= 4).mean()), 3), round(float((mcat[mnp] >= 4).mean()), 3), " udział np:", round(float(mnp.mean()), 3))
fp, fn_ = (mcat[~mnp] == 5).mean(), (mcat[mnp] == 5).mean(); p0 = mnp.mean()
# R_top(k) = [(1-p0) fp + k p0 fn] / [(1-p0) + k p0] / f_top ;  rozwiązanie dla R = 0,72
ft = (mcat == 5).mean()
ks = np.linspace(0, 200, 20001); Rk = ((1 - p0) * fp + ks * p0 * fn_) / ((1 - p0) + ks * p0) / ft
out["np_solve"] = {"f25_prompt": float(fp), "f25_np": float(fn_), "p0": float(p0), "R_min(k->inf)": float(fn_ / ft), "R(k=2)": float(np.interp(2, ks, Rk)),
                   "R(k=0)": float(Rk[0]), "k_for_0.72": (float(ks[np.argmin(np.abs(Rk - 0.72))]) if Rk.min() < 0.72 else None)}
print("nieprompt: ", out["np_solve"])

for run in R.RUNS:
    d = C.load(f"data_{run}")
    wpt_all, _ = C.pt_weights(run); wpt = wpt_all[sm]; wpt = wpt / wpt.mean()
    # rozkład ptcone40 sygnału danych w regionie >= WP70 (okno ±2 − k·pasma), biny = decyle MC
    sd = C.base_mask(d, "8-30", isdata=True, run=run)
    r70 = sd & (d["score"] >= thr[2])
    ib, ok = C.dm_bin(d["dm"][r70]); h = np.bincount(ib[ok], minlength=R.NB).astype(float)
    cen = R.CEN; sbb = ((cen > 141) & (cen < 143)) | ((cen > 149) & (cen < 155)); wb = np.abs(cen - mu) < 2.0
    p = C.fit_bkg(h, np.abs(cen - mu) > 3.0, 2); bs = p[0] * R.bkg_shape(*p[1:]); k = bs[wb].sum() / bs[sbb].sum()
    iw = r70 & (np.abs(d["dm"] - mu) < 2.0); isb = r70 & (((d["dm"] > 141) & (d["dm"] < 143)) | ((d["dm"] > 149) & (d["dm"] < 155)))
    m70 = msc >= thr[2]
    ed = np.quantile(mcone[m70], np.linspace(0, 1, 11)[1:-1]); nbc = len(ed) + 1
    pd_ = np.bincount(np.searchsorted(ed, d["ptcone40"][iw], side="right"), minlength=nbc) - k * np.bincount(np.searchsorted(ed, d["ptcone40"][isb], side="right"), minlength=nbc)
    pd_ = np.clip(pd_, 0, None) / pd_.sum()
    bc = np.searchsorted(ed, mcone, side="right")

    def rake(w0, use_pt, use_cone, sel):
        w = np.where(sel, w0, 0.0).astype(float)
        if use_pt:
            w = w * wpt
        if use_cone:
            for _ in range(8):
                cur = np.bincount(bc[m70], weights=w[m70], minlength=nbc); cur = cur / cur.sum()
                w = w * (pd_ / np.clip(cur, 1e-9, None))[bc]
                if use_pt:           # utrzymanie widma pT (w regionie >= WP70) po przeważeniu ptcone
                    pb = np.clip(np.searchsorted(C.PT_EDGES[run], mpt, side="right") - 1, 0, len(C.PT_EDGES[run]) - 2)
                    tgt = np.bincount(pb[m70], weights=(np.where(sel, w0, 0.0) * wpt)[m70], minlength=len(C.PT_EDGES[run]) - 1)
                    cur = np.bincount(pb[m70], weights=w[m70], minlength=len(C.PT_EDGES[run]) - 1)
                    w = w * (tgt / np.clip(cur, 1e-9, None))[pb]
        return w
    for name, (lo, hi) in (("8-30", (8, 30)), ("8-12", (8, 12)), ("12-30", (12, 30))):
        sdd = sd & (d["pt"] >= lo) & (d["pt"] < hi)
        cd = R.cat_index(d["score"][sdd], thr); ib, ok = C.dm_bin(d["dm"][sdd]); flat = (cd * R.NB + ib)[ok]; inv = d["inv"][sdd][ok]
        S, P, nco = C.count_yields(np.bincount(flat, minlength=6 * R.NB).reshape(6, R.NB).astype(float), mu, W, 3.0)
        rng = np.random.default_rng(4242 + run)
        SB = np.array([C.count_yields(np.bincount(flat, weights=rng.poisson(1.0, d["nev"]).astype(float)[inv], minlength=6 * R.NB).reshape(6, R.NB), mu, W, 3.0, nco, P)[0] for _ in range(NB)])
        cov = np.cov(SB.T)
        sel = (mpt >= lo) & (mpt < hi)
        f0 = np.bincount(mcat[sel], minlength=6) / sel.sum()
        for tag, up, uc in (("MC", 0, 0), ("MC pT", 1, 0), ("MC ptcone40", 0, 1), ("MC pT+ptcone40", 1, 1)):
            w = rake(np.ones(len(mpt)), up, uc, sel)
            f = np.bincount(mcat, weights=w, minlength=6) / w.sum(); neff = w.sum() ** 2 / (w ** 2).sum()
            # chi2 na 5 kategoriach czystych+pośrednich (1..5; kat. 0 = odrzucone ma S/B ~ 0 i wchodzi przez kowariancję) — pełne 6, N profilowane
            Ci = np.linalg.inv(cov); N = (f @ Ci @ S) / (f @ Ci @ f)
            Ct = cov + N * N * (np.diag(f) - np.outer(f, f)) / neff; Ci = np.linalg.inv(Ct); N = (f @ Ci @ S) / (f @ Ci @ f)
            r = S - N * f; c2 = float(r @ Ci @ r)
            # obserwable czyste z błędem bootstrapowym danych ⊕ MC (dwumian, n_eff)
            def obs(Sv, fv):
                return np.array([(Sv[5] / Sv[4:].sum()) / (fv[5] / fv[4:].sum()), (Sv[5] / Sv[3:].sum()) / (fv[5] / fv[3:].sum())])
            o = obs(S, f); ob = np.array([obs(s, f) for s in SB]); 
            q50 = f[5] / f[4:].sum(); q70 = f[5] / f[3:].sum()
            emc = np.array([np.sqrt(q50 * (1 - q50) / (neff * f[4:].sum())) / q50, np.sqrt(q70 * (1 - q70) / (neff * f[3:].sum())) / q70]) * o
            err = np.sqrt(ob.std(0, ddof=1) ** 2 + emc ** 2)
            t = {"run": run, "bin": name, "mc": tag, "chi2": c2, "ndf": 5, "p": float(CHI2.sf(c2, 5)), "neff": float(neff), "f": f.tolist(),
                 "pred_Rtop": float(f[5] / f0[5]), "pred_rho50": float((f[5] / f[4:].sum()) / (f0[5] / f0[4:].sum())),
                 "rho50_data_over_pred": [float(o[0]), float(err[0])], "rho70_data_over_pred": [float(o[1]), float(err[1])]}
            out["tests"].append(t)
            print(f"{run} pT {name:5s} {tag:15s}: chi2/5={c2:5.1f} (p={t['p']:.2g}) n_eff={neff:6.0f} przewid. Rtop={t['pred_Rtop']:.3f} rho50={t['pred_rho50']:.3f} | dane/przewid.: rho50={o[0]:.3f}±{err[0]:.3f} rho70={o[1]:.3f}±{err[1]:.3f}")
json.dump(out, open(os.path.join(C.OUT, "env.json"), "w"), indent=1)
