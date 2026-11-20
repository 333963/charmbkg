"""HIGHSCORE kroki 2 i 4 — kształt rozkładu score sygnału w drobnych kwantylach MC, modele transformacji score (logit),
skład MC (prompt / nieprompt / swap). EKSPLORACYJNE, PO FAKCIE. python hs_shape.py [NBOOT] -> out/shape.json"""
import json, os, sys
import multiprocessing as mp
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as CHI2
import hs_common as C
R = C.R
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 300
WPS = [0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]
WPS6 = list(R.WP_MAIN)
JOBS = [(428770, (8, 30), WPS), (428770, (8, 12), WPS), (428770, (12, 30), WPS), (427929, (8, 30), WPS6), (428770, (5, 8), WPS6)]


def lgt(s):
    s = np.clip(s, 1e-7, 1 - 1e-7)
    return np.log(s / (1 - s))


def job(a):
    run, b, wps = a
    parent = "5-8" if b[0] < 8 else "8-30"
    d = C.load(f"data_{run}"); M = C.load("fd_nd")
    mu = C.MU[parent]; W = 3.0
    smr = C.base_mask(M, b, isdata=False) & (np.abs(M["dm"] - mu) < W)
    thr = np.quantile(M["score"][smr], 1 - np.array(wps))
    nc = len(thr) + 1
    sd = C.base_mask(d, b, isdata=True, run=run)
    cd = R.cat_index(d["score"][sd], thr)
    ib, ok = C.dm_bin(d["dm"][sd]); flat = (cd * R.NB + ib)[ok]; inv = d["inv"][sd][ok]
    H = np.bincount(flat, minlength=nc * R.NB).reshape(nc, R.NB).astype(float)
    S, P, nco = C.count_yields(H, mu, W, 3.0)
    rng = np.random.default_rng(77 + run)
    SB = []
    for _ in range(NB):
        w = rng.poisson(1.0, d["nev"]).astype(float)[inv]
        SB.append(C.count_yields(np.bincount(flat, weights=w, minlength=nc * R.NB).reshape(nc, R.NB), mu, W, 3.0, nco, P)[0])
    cov = np.cov(np.array(SB).T)
    # --- MC: sygnał (+ swap osobno) w tym samym oknie
    lg = lgt(M["score"][smr]); lt = lgt(thr)
    isnp = M["np"][smr].astype(float)
    if parent == "8-30":
        wpt_all, pttab = C.pt_weights(run)
        wpt = wpt_all[smr]; wpt = wpt / wpt.mean()
    else:
        wpt, pttab = np.ones(smr.sum()), None
    msw = (np.abs(M["m"] - R.M0_D0) < R.NOM_WIN) & (M["pt"] >= b[0]) & (M["pt"] < b[1]) & M["swap"] & (np.abs(M["dm"] - mu) < W)
    cat_sw = R.cat_index(M["score"][msw], thr)
    f_sw = np.bincount(cat_sw, minlength=nc) / max(msw.sum(), 1)
    x_sw = msw.sum() / smr.sum()
    z = np.random.default_rng(5).normal(size=len(lg))

    def frac(lgx, w):
        return np.bincount(np.searchsorted(lt, lgx, side="right"), weights=w, minlength=nc) / w.sum()

    def chi2f(f, w):
        neff = w.sum() ** 2 / (w ** 2).sum()
        Ci = np.linalg.inv(cov)
        N = (f @ Ci @ S) / (f @ Ci @ f)
        Ct = cov + N * N * (np.diag(f) - np.outer(f, f)) / neff
        Ci = np.linalg.inv(Ct)
        N = (f @ Ci @ S) / (f @ Ci @ f)
        r = S - N * f
        return float(r @ Ci @ r), float(N)

    models = {
        "identycznosc": (0, lambda p, w: frac(lg, w), []),
        "przesuniecie logit a": (1, lambda p, w: frac(lg + p[0], w), [[0.0], [-0.5]]),
        "a + b*logit": (2, lambda p, w: frac(p[0] + p[1] * lg, w), [[0.0, 1.0], [-0.3, 0.8], [0.3, 0.7]]),
        "kompresja ogona (l0, b)": (2, lambda p, w: frac(np.where(lg > p[0], p[0] + p[1] * (lg - p[0]), lg), w), [[2.0, 0.7], [3.0, 0.5], [4.0, 0.3], [1.0, 0.8]]),
        "rozmycie logit sigma": (1, lambda p, w: frac(lg + abs(p[0]) * z, w), [[0.5], [1.5]]),
        "a + rozmycie sigma": (2, lambda p, w: frac(lg + p[0] + abs(p[1]) * z, w), [[0.0, 0.5], [-0.3, 1.0]]),
        "nieprompt x k": (1, lambda p, w: frac(lg, w * (1 + (abs(p[0]) - 1) * isnp)), [[1.0], [3.0]]),
        "swap dodany (ulamek MC, 0 par.)": (0, lambda p, w: (frac(lg, w) + x_sw * f_sw) / (1 + x_sw), []),
        "swap x (swobodny)": (1, lambda p, w: (frac(lg, w) + abs(p[0]) * f_sw) / (1 + abs(p[0])), [[0.1], [0.5]]),
    }
    res = {"run": run, "bin": list(b), "wps": wps, "thr": thr.tolist(), "S": S.tolist(), "S_err": np.sqrt(np.diag(cov)).tolist(),
           "cov": cov.tolist(), "n_mc": int(smr.sum()), "pt_tab": pttab, "x_swap_mc": float(x_sw), "f_swap": f_sw.tolist(),
           "neff_ptw": float(wpt.sum() ** 2 / (wpt ** 2).sum()), "models": {}}
    f0 = frac(lg, np.ones(len(lg)))
    res["f_mc"] = f0.tolist(); res["f_mc_ptw"] = frac(lg, wpt).tolist()
    for tag, w in (("", np.ones(len(lg))), ("ptw", wpt)):
        npf = np.bincount(np.searchsorted(lt, lg, side="right"), weights=w * isnp, minlength=nc) / np.bincount(np.searchsorted(lt, lg, side="right"), weights=w, minlength=nc)
        res["np_frac" + tag] = npf.tolist(); res["np_tot" + tag] = float((w * isnp).sum() / w.sum())
        for name, (npar, fn, starts) in models.items():
            best = None
            if npar == 0:
                c, N = chi2f(fn([], w), w); best = (c, [], N)
            else:
                for st in starts:
                    o = minimize(lambda p: chi2f(fn(p, w), w)[0], st, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 600})
                    if best is None or o.fun < best[0]:
                        best = (float(o.fun), [float(x) for x in o.x], chi2f(fn(o.x, w), w)[1])
            ndf = nc - 1 - npar
            f = fn(best[1], w)
            res["models"][name + (" | MC pT-waz." if tag else "")] = {"chi2": best[0], "ndf": ndf, "p": float(CHI2.sf(best[0], ndf)), "par": best[1], "N": best[2],
                                                                    "pred_R": (f / f0).tolist()}
    return res


if __name__ == "__main__":
    C.load("fd_nd"); [C.load(f"data_{r}") for r in R.RUNS]
    with mp.get_context("fork").Pool(5) as P:
        out = P.map(job, JOBS, chunksize=1)
    json.dump(out, open(os.path.join(C.OUT, "shape.json"), "w"), indent=1)
    for r in out:
        S = np.array(r["S"]); f = np.array(r["f_mc"]); N = S.sum()
        print(f"\n== {r['run']} pT {r['bin']}  N={N:.0f}  n_mc={r['n_mc']} neff_ptw={r['neff_ptw']:.0f} x_swap={r['x_swap_mc']:.3f} np_tot={r['np_tot']:.3f} np_tot_ptw={r['np_totptw']:.3f}")
        if r["pt_tab"]:
            print("  pT w:", np.round(r["pt_tab"]["w"], 2).tolist(), "frac_data:", np.round(r["pt_tab"]["frac_data"], 3).tolist())
        print("  kat(cum eff)  R=f_d/f_mc ± err | R wzgl. MC pT-waz. | np_frac | f_swap")
        for c in range(len(S)):
            lab = "<%.2f" % r["wps"][0] if c == 0 else (">%.2f" % r["wps"][-1] if c == len(S) - 1 else "%.2f-%.2f" % (r["wps"][c - 1], r["wps"][c]))
            print(f"  {lab:10s} S={S[c]:7.1f}±{r['S_err'][c]:5.1f}  R={S[c]/N/f[c]:.2f}±{r['S_err'][c]/N/f[c]:.2f}  Rptw={S[c]/N/r['f_mc_ptw'][c]:.2f}  np={r['np_frac'][c]:.3f}  sw={r['f_swap'][c]:.3f}")
        for k, m in r["models"].items():
            print(f"  {k:48s} chi2/ndf={m['chi2']:6.1f}/{m['ndf']:2d} p={m['p']:.2g} par={np.round(m['par'],3).tolist()} predR_top={m['pred_R'][-1]:.3f}")
