#!/usr/bin/env python3
"""CATFIT — proby wg KRYTERIA.md §6–7 (POST HOC wzgledem XSEC_final; tylko piki inkluzywne 427929, bez AFP).
  cf_run.py mc | toys | data [--nproc 8]   ->  out/{mc,toys,data}.json"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cf_lib as C  # noqa: E402
F = C.F
FIN = json.load(open(os.path.join(C.CW, "XSEC_final_2026-11-20", "out", "final.json")))["cells"]
MODELS = {"D0": ["M8"], "Dplus": ["M0", "A2"]}
WPS = [0.7, 0.6, 0.5, 0.4, 0.3]
NM4 = list(C.K4)


def cats_of(sp, bn, K, ks=1.0):
    return [C.cat(sp, bn, *K[n], kscale=ks) for n in K]


def single(f):
    return {"N": f["N"][0], "err": (f.get("err") or [None])[0], "valid": f["valid"], "at_limit": f["at_limit"], "sc": f["sc"], "dmu": f["dmu"], "chi2": f["chi2"], "ndf": f["ndf"],
            "p": F.pval(f["chi2"], f["ndf"]), "fval": f["fval"]}


def cons_of(cats, SF, SFe, which):
    e = np.array([c["e"] for c in cats]); e = e / e.sum()
    out = []
    for j in ((1,), (1, 2), (1, 2, 3))[:{"F70": 1, "F3": 3}[which]]:
        ix = list(range(j[-1], 4))
        E = e[ix].sum()
        out.append((ix, E * SF, E * SFe))
    return out


# ------------------------------------------------------------------------------------------------ MC
def t_mc(t):
    sp, bn, mn = t
    C.init()
    rng = np.random.default_rng(9100 + 10 * C.BINS.index(bn) + (100 if sp == "D0" else 0) + (7 if mn == "A2" else 0))
    res = {"sp": sp, "bin": bn, "model": mn}
    cats = cats_of(sp, bn, C.K4)
    H = np.stack([C.hist(c, "mb|all") for c in cats])
    tr = np.array([C.hist(c, "mb|sig").sum() for c in cats])
    f = C.fitK(H, cats, sp, mn)
    s = C.summ(f, NM4)
    reps = [rng.poisson(H).astype(float) for _ in range(40)]
    rT, rS, rf = [], [], []
    for R in reps:
        g = C.fitK(R, cats, sp, mn, start=f["x"], hesse=False)
        if g["valid"]:
            n = np.array(g["N"]); rT.append(n.sum()); rS.append(n[1:].sum()); rf.append(n / n.sum())
    res["K4"] = dict(s, truth=dict(zip(NM4, tr.tolist())), truth_tot=float(tr.sum()), rep_rms_tot=float(np.std(rT, ddof=1)), rep_rms_sub=float(np.std(rS, ddof=1)),
                     rep_rms_f=dict(zip(NM4, np.std(np.array(rf), axis=0, ddof=1).tolist())), n_rep=len(rT), fd_e=dict(zip(NM4, [c["e"] for c in cats])), ncand=dict(zip(NM4, H.sum(1).tolist())))
    res["WP"] = {}
    for wp in WPS + [1.0]:
        w = C.cat(sp, bn, *C.WPC[wp])
        Hw = C.hist(w, "mb|all")
        g = C.fitK(Hw[None], [w], sp, mn)
        rN = []
        for i in range(40):
            R = rng.poisson(Hw).astype(float)
            gg = C.fitK(R[None], [w], sp, mn, start=g["x"], hesse=False)
            if gg["valid"]:
                rN.append(gg["N"][0])
        res["WP"][str(wp)] = dict(single(g), truth=float(C.hist(w, "mb|sig").sum()), rep_rms=float(np.std(rN, ddof=1)), n_rep=len(rN), e=w["e"])
    return res


# ------------------------------------------------------------------------------------------------ zabawki
def t_toys(t):
    sp, bn, mn, dord, ntoy, seed = t
    C.init()
    c = F.SP[sp]
    S = C.SET[sp][mn]
    fe = C.full_edges(sp)
    ctr = 0.5 * (fe[1:] + fe[:-1])
    blind = np.abs(ctr - c["m0"]) < c["win"]
    cats = cats_of(sp, bn, C.K4)
    w70 = C.cat(sp, bn, 6, 19)
    e = np.array([k["e"] for k in cats]); e = e / e.sum()
    E70 = e[1:].sum()
    fin = FIN[f"{sp}|{bn}"]
    Ninj = fin["N"] / E70 * e
    SFr = fin["SF_err"] / fin["SF"]
    bkg, info = [], {}
    for nm, k in zip(NM4, cats):
        Hd = C.hist(k, C.SRC)[0]                       # TYLKO PASMA (okno zaslepione)
        fb = C.CF.fit(Hd[None], [k["sigT"]], [[]], fe, c["m0"], S["nbkg"] + dord, shape="fix", win=c["win"], blind=blind, hesse=False)
        bkg.append(fb["mu"][0])
        info[nm] = {"bkg_win": float(fb["mu"][0][blind].sum()), "N_inj": float(Ninj[len(bkg) - 1]), "sb_valid": fb["valid"], "chi2_sb": fb["chi2"], "ndf_sb": fb["ndf"]}
    sw = lambda k, dm, s: k["sigT"].bins(fe, c["m0"], dm, s) + (k["frefl"] * k["refT"].bins(fe, c["m0"], dm, 1.0) if k["refT"] is not None else 0.0) + k["kap"] * k["peakT"].bins(fe, c["m0"], dm, 1.0)  # noqa: E731
    rng = np.random.default_rng(seed)
    rows = []
    nan = float("nan")
    for _ in range(ntoy):
        dm, s = rng.normal(0, C.GC[sp][0]), min(max(rng.normal(1, 0.10), 0.75), 1.4)
        Hs = np.stack([rng.poisson(bkg[i] + Ninj[i] * sw(cats[i], dm, s)).astype(float) for i in range(4)])
        d = rng.normal(0, SFr)
        r = []
        f1 = C.fitK(Hs[1:].sum(0)[None], [w70], sp, mn)
        r += [f1["N"][0], (f1.get("err") or [nan])[0], float(f1["valid"] and not f1["at_limit"])]
        f0 = C.fitK(Hs, cats, sp, mn)
        r += [f0["Ntot"], f0.get("Ntot_err", nan), float(f0["valid"] and not f0["at_limit"])] + list(np.array(f0["N"]) / max(f0["Ntot"], 1e-9))
        cons = [(ix, cen * (1 + d), sg) for ix, cen, sg in cons_of(cats, 1.0, SFr, "F70")]
        f7 = C.fitK(Hs, cats, sp, mn, cons=cons, start=f0["x"])
        r += [f7["Ntot"], f7.get("Ntot_err", nan), float(f7["valid"] and not f7["at_limit"])]
        ff = C.fitK(Hs, cats, sp, mn, fracs=e)
        r += [ff["Ntot"], ff.get("Ntot_err", nan), float(ff["valid"] and not ff["at_limit"]), d]
        rows.append(r)
    return {"sp": sp, "bin": bn, "model": mn, "dord": dord, "seed": seed, "info": info, "N70_inj": float(Ninj[1:].sum()), "Ntot_inj": float(Ninj.sum()), "SF_rel": SFr, "rows": rows}


# ------------------------------------------------------------------------------------------------ dane
def t_data(t):
    sp, bn, mn = t
    C.init()
    fin = FIN[f"{sp}|{bn}"]
    SF, SFe = fin["SF"], fin["SF_err"]
    res = {"sp": sp, "bin": bn, "model": mn, "SF": SF, "SF_err": SFe}
    cats = cats_of(sp, bn, C.K4)
    D = np.stack([C.hist(c, C.SRC) for c in cats], axis=1)          # [nrep+1, 4, nfull]
    nrep = D.shape[0] - 1
    e = np.array([c["e"] for c in cats]); e = e / e.sum()
    res["fd_e"] = dict(zip(NM4, e.tolist())); res["ncand"] = dict(zip(NM4, D[0].sum(1).tolist()))
    res["kap"] = dict(zip(NM4, [c["kap"] for c in cats])); res["frefl"] = dict(zip(NM4, [c["frefl"] for c in cats]))
    f0 = C.fitK(D[0], cats, sp, mn)
    o = C.summ(f0, NM4)
    bT, bS, bf = [], [], []
    for i in range(1, nrep + 1):
        g = C.fitK(D[i], cats, sp, mn, start=f0["x"], hesse=False)
        if g["valid"]:
            n = np.array(g["N"]); bT.append(n.sum()); bS.append(n[1:].sum()); bf.append(n / n.sum())
    o.update({"boot_tot": float(np.std(bT, ddof=1)), "boot_sub": float(np.std(bS, ddof=1)), "boot_f": dict(zip(NM4, np.std(np.array(bf), axis=0, ddof=1).tolist())), "n_boot": len(bT)})
    bE = np.array(bf)[:, 1:].sum(1)
    o["boot_E70"] = float(np.std(bE, ddof=1))
    o["var"] = {v: (lambda g: {"Ntot": g["Ntot"], "valid": g["valid"], "at_limit": g["at_limit"], "p": F.pval(g["chi2"], g["ndf"]), "N": g["N"]})(C.fitK(D[0], cats, sp, v)) for v in C.VARS[sp][mn]}
    g = C.fitK(D[0], cats, sp, mn, nbkg_add=[2, 0, 0, 0])
    o["L+2"] = C.summ(g, NM4)
    if C.SET[sp][mn]["peak"]:
        for tag, ks in (("k+5.8", 1.058), ("k-5.8", 0.942)):
            o[tag] = C.fitK(D[0], cats_of(sp, bn, C.K4, ks), sp, mn)["Ntot"]
    if sp == "D0":
        alt, rows = None, []
        for sc0 in (0.80, 0.85, 0.90, 0.945, 1.0, 1.05, 1.10, 1.20):
            g = C.fitK(D[0], cats, sp, mn, start=f0["x"], sc_start=sc0)
            r = {"sc0": sc0, "Ntot": g["Ntot"], "sc": g["sc"], "dfval": g["fval"] - f0["fval"], "valid": g["valid"], "at_limit": g["at_limit"]}
            rows.append(r)
            if g["valid"] and not g["at_limit"] and abs(g["Ntot"] / f0["Ntot"] - 1) > 0.01 and r["dfval"] < 1.0 and (alt is None or abs(r["Ntot"] - f0["Ntot"]) > abs(alt["Ntot"] - f0["Ntot"])):
                alt = r
        o["altmin"], o["alt_starts"] = alt, rows
    res["F0"] = o
    for which in ("F70", "F3"):
        cons = cons_of(cats, SF, SFe, which)
        g = C.fitK(D[0], cats, sp, mn, cons=cons, start=f0["x"])
        oo = C.summ(g, NM4)
        g2 = C.fitK(D[0], cats, sp, mn, fracs=np.array(g["N"]) / g["Ntot"])          # czesc statystyczna: ulamki zamrozone w minimum
        oo["Ntot_err_fixfrac"] = g2.get("Ntot_err")
        bT = []
        for i in range(1, nrep + 1):
            gg = C.fitK(D[i], cats, sp, mn, cons=cons, start=g["x"], hesse=False)
            if gg["valid"]:
                bT.append(gg["Ntot"])
        oo["boot_tot"] = float(np.std(bT, ddof=1)); oo["n_boot"] = len(bT)
        oo["var"] = {v: C.fitK(D[0], cats, sp, v, cons=cons)["Ntot"] for v in C.VARS[sp][mn]}
        res[which] = oo
    g = C.fitK(D[0], cats, sp, mn, fracs=e)
    res["Ffix"] = {"Ntot": g["Ntot"], "Ntot_err": g.get("Ntot_err"), "valid": g["valid"], "p": F.pval(g["chi2"], g["ndf"]), "sc": g["sc"]}
    # kontrole: K3, K5
    for nm, K in (("K3", C.K3), ("K5", C.K5)):
        ck = cats_of(sp, bn, K)
        g = C.fitK(np.stack([C.hist(c, C.SRC)[0] for c in ck]), ck, sp, mn)
        res[nm] = C.summ(g, list(K))
    # pojedyncze punkty pracy (B2) + WP70
    res["WP"] = {}
    for wp in WPS:
        w = C.cat(sp, bn, *C.WPC[wp])
        Dw = C.hist(w, C.SRC)
        g = C.fitK(Dw[0][None], [w], sp, mn)
        oo = single(g)
        bN = []
        for i in range(1, nrep + 1):
            gg = C.fitK(Dw[i][None], [w], sp, mn, start=g["x"], hesse=False)
            if gg["valid"]:
                bN.append(gg["N"][0])
        oo["boot"] = float(np.std(bN, ddof=1)); oo["n_boot"] = len(bN); oo["e"] = w["e"]
        oo["var"] = {v: C.fitK(Dw[0][None], [w], sp, v)["N"][0] for v in C.VARS[sp][mn]}
        res["WP"][str(wp)] = oo
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["mc", "toys", "data"]); ap.add_argument("--nproc", type=int, default=8); ap.add_argument("--ntoy", type=int, default=150)
    a = ap.parse_args()
    t0 = time.time()
    cells = [(sp, bn, mn) for sp in C.SPS for bn in C.BINS for mn in MODELS[sp]]
    if a.cmd == "mc":
        fn, tasks = t_mc, cells
    elif a.cmd == "data":
        fn, tasks = t_data, cells
    else:
        fn, tasks = t_toys, []
        for ic, (sp, bn, mn) in enumerate(cells):
            for ch in range(a.ntoy // 25):
                tasks.append((sp, bn, mn, 0, 25, 31000 + 100 * ic + ch))
            for ch in range(2):
                tasks.append((sp, bn, mn, 1, 30, 51000 + 100 * ic + ch))
    out = []
    with Pool(a.nproc) as p:
        for i, r in enumerate(p.imap_unordered(fn, tasks)):
            out.append(r)
            print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(tasks)} {r['sp']} {r['bin']} {r['model']}", flush=True)
    out.sort(key=lambda r: (r["sp"], C.BINS.index(r["bin"]), r["model"], r.get("dord", 0), r.get("seed", 0)))
    fo = os.path.join(C.BASE, "out", f"{a.cmd}.json")
    json.dump({"res": out, "sec": time.time() - t0}, open(fo, "w"), default=float)
    print(f"zapisano {fo}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
