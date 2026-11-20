"""HIGHSCORE krok 1 — odporność deficytu > WP25 (EKSPLORACYJNE, PO FAKCIE). python hs_robust.py [NBOOT] -> out/robust.json"""
import json, os, sys, time
import multiprocessing as mp
import numpy as np
import hs_common as C
R = C.R

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 200


def V(name, run, parent="8-30", pt=None, dsel=None, msel=None, W=3.0, side=3.0, emax=None, nco=None, mwin=R.NOM_WIN, fit=False):
    return dict(name=name, run=run, parent=parent, pt=pt, dsel=dsel, msel=msel, W=W, side=side, emax=emax, nco=nco, mwin=mwin, fit=fit)


def variants():
    L = []
    for run in R.RUNS:
        L.append(V("nominal (okno ±3, pasma >3)", run, fit=True))
        L.append(V("5-8 GeV", run, parent="5-8", fit=True))
        for W in (1.5, 2.0, 4.0):
            L.append(V(f"okno ±{W}", run, W=W, side=max(W, 3.0)))
        L.append(V("pasma >5 MeV", run, side=5.0))
        L.append(V("pasma tylko <155 MeV", run, emax=155.0))
        L.append(V("tło: c2=0", run, nco=1))
        L.append(V("tło: sam próg u^a", run, nco=0))
        L.append(V("m(Kπ) ±30 MeV", run, mwin=0.030))
        L.append(V("m(Kπ) ±20 MeV", run, mwin=0.020))
        for lo, hi in ((8, 12), (12, 30), (8, 10), (10, 12), (12, 15), (15, 20), (20, 30), (15, 30)):
            L.append(V(f"pT {lo}-{hi}", run, pt=(lo, hi)))
        L.append(V("|eta|<1", run, dsel=("abs_eta", 0, 1.0), msel=("abs_eta", 0, 1.0)))
        L.append(V("|eta|>=1", run, dsel=("abs_eta", 1.0, 9), msel=("abs_eta", 1.0, 9)))
        L.append(V("npu=0 (dane)", run, dsel=("e_npu", -1, 0.5)))
    r = 428770
    for c in ("J12only", "MBTSbased", "other"):
        L.append(V(f"klasa {c}", r, dsel=("cls_" + c,)))
        L.append(V(f"5-8 GeV klasa {c}", r, parent="5-8", dsel=("cls_" + c,)))
    for lo, hi in ((0, 45), (45, 65), (65, 999)):
        L.append(V(f"nSelBL {lo}-{hi} (dane i MC)", r, dsel=("e_nSelBL", lo, hi), msel=("e_nSelBL", lo, hi)))
        L.append(V(f"nSelBL {lo}-{hi} (dane; MC całe)", r, dsel=("e_nSelBL", lo, hi)))
    for lo, hi in ((0, 30), (30, 45), (45, 999)):
        L.append(V(f"nSelBL {lo}-{hi} (dane i MC)", 427929, dsel=("e_nSelBL", lo, hi), msel=("e_nSelBL", lo, hi)))
    for lo, hi in ((-1, 2.0), (2.0, 6.0), (6.0, 1e9)):
        L.append(V(f"ptcone40 {lo}-{hi:g} (dane i MC)", r, dsel=("ptcone40", lo, hi), msel=("ptcone40", lo, hi)))
    for c in ("J12only", "MBTSbased"):
        for lo, hi in ((8, 12), (12, 30)):
            L.append(V(f"klasa {c}, pT {lo}-{hi}", r, pt=(lo, hi), dsel=("cls_" + c,)))
    L.append(V("LB < 160", r, dsel=("e_lb", 0, 160)))
    L.append(V("LB >= 160", r, dsel=("e_lb", 160, 1e9)))
    return L


def sel(d, s):
    if s is None:
        return np.ones(len(d["m"]), bool)
    if len(s) == 1:
        return d[s[0]]
    return (d[s[0]] >= s[1]) & (d[s[0]] < s[2])


def job(v):
    t0 = time.time()
    d = C.load(f"data_{v['run']}"); M = C.load("fd_nd")
    b = v["pt"] or v["parent"]
    mu = C.MU[v["parent"]]
    thr = C.thr_for(v["parent"], R.WP_MAIN)
    sd = C.base_mask(d, b, v["mwin"], True, v["run"]) & sel(d, v["dsel"])
    sm = C.base_mask(M, b, v["mwin"], False) & sel(M, v["msel"]) & (np.abs(M["dm"] - mu) < v["W"])
    cd = R.cat_index(d["score"][sd], thr); cm = R.cat_index(M["score"][sm], thr)
    ib, ok = C.dm_bin(d["dm"][sd])
    flat = (cd * R.NB + ib)[ok]; inv = d["inv"][sd][ok]
    H = np.bincount(flat, minlength=6 * R.NB).reshape(6, R.NB).astype(float)
    ncof = None if v["nco"] is None else np.full(6, v["nco"])
    S, P, nco = C.count_yields(H, mu, v["W"], v["side"], ncof, None, v["emax"])
    nm = np.bincount(cm, minlength=6).astype(float)
    fm = nm / nm.sum()
    nom = C.derived(S, fm)
    res = {"name": v["name"], "run": v["run"], "bin": str(b), "S": S.tolist(), "f_mc": fm.tolist(), "n_mc": int(nm.sum()),
           "n_pairs": int(sd.sum()), "nco": nco.tolist()}
    if v["fit"]:
        try:
            ff = R.CatFit(H, R.official()["shapes"][v["parent"]]).fit(hesse=False)
            Mf = C.base_mask(M, b, v["mwin"], False)
            fmf = np.bincount(R.cat_index(M["score"][Mf], thr), minlength=6) / Mf.sum()
            res["fit_R"] = (ff["S"] / ff["S"].sum() / fmf).tolist(); res["fit_S"] = ff["S"].tolist()
        except Exception as e:  # noqa: BLE001
            res["fit_R"] = str(e)
    rng = np.random.default_rng(1234 + v["run"])
    keys = ["rho50", "rho70", "SF50", "SF70", "N"]
    acc = {k: [] for k in keys}; RR = []
    for _ in range(NB):
        w = rng.poisson(1.0, d["nev"]).astype(float)[inv]
        Hb = np.bincount(flat, weights=w, minlength=6 * R.NB).reshape(6, R.NB)
        Sb, _, _ = C.count_yields(Hb, mu, v["W"], v["side"], nco, P, v["emax"])
        fb = np.bincount(cm, weights=rng.poisson(1.0, len(cm)), minlength=6).astype(float); fb /= fb.sum()
        dd = C.derived(Sb, fb)
        RR.append(dd["R"]); [acc[k].append(dd[k]) for k in keys]
    RR = np.array(RR)
    res["R"] = nom["R"].tolist(); res["R_err"] = RR.std(0, ddof=1).tolist()
    for k in keys:
        res[k] = float(nom[k]); res[k + "_err"] = float(np.std(acc[k], ddof=1))
    res["sec"] = round(time.time() - t0, 1)
    return res


if __name__ == "__main__":
    t0 = time.time()
    C.load("fd_nd"); [C.load(f"data_{r}") for r in R.RUNS]
    with mp.get_context("fork").Pool(8) as P:
        out = P.map(job, variants(), chunksize=1)
    json.dump(out, open(os.path.join(C.OUT, "robust.json"), "w"), indent=1)
    for r in out:
        print(f"{r['run']} {r['name']:36s} N={r['N']:7.0f}±{r['N_err']:4.0f} Rtop={r['R'][5]:.3f}±{r['R_err'][5]:.3f} R4={r['R'][4]:.3f}±{r['R_err'][4]:.3f} "
              f"rho50={r['rho50']:.3f}±{r['rho50_err']:.3f} rho70={r['rho70']:.3f}±{r['rho70_err']:.3f} SF50={r['SF50']:.3f}±{r['SF50_err']:.3f} fmc_top={r['f_mc'][5]:.3f} nmc={r['n_mc']}"
              + (f" fitRtop={r['fit_R'][5]:.3f}" if isinstance(r.get('fit_R'), list) else ""))
    print("czas", round(time.time() - t0), "s")
