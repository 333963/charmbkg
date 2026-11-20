#!/usr/bin/env python3
"""SFABS metoda A (drabina): SF70 = SF(L) x SF_rel(70|L) z fitu lacznego kategorii rozlacznych P = [WP70, .) i F_L = [WP_L, WP70).
Kryteria: ../KRYTERIA.md par. A. Modul fitu: FITBIAS fb_fit/fb_var (import, bez zmian). D+: M0, D0: M8 (kappa x 1,364 na danych).
  sfa_ladder.py mc   OUT     b_L (mb_nd, pelne tlo) + wstrzykniecia G_A1
  sfa_ladder.py data OUT     dane, 2 runy, bootstrap 30
  sfa_ladder.py toys OUT DATAJSON   G_A2 (tlo z pasm danych + sygnal MC), L = 0,90
"""
import json
import math
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, "/home/k/Documents/charm/charm_work/FITBIAS_2026-11-19/skrypty")
import fb_fit as F  # noqa: E402
import fb_var as V  # noqa: E402

MODEL = {"Dplus": "M0", "D0": "M8"}
KSC = {"D0": 1.364, "Dplus": 1.068}
LS = (0.9, 0.95, 0.8)
SRC = {"427929": "427929|MB427929", "428770": "428770|all"}
NREP = 30


def kappa(r, e7, eL):
    """stale kappa z r = SF70/SFL."""
    den = eL * (1 - e7) - r * e7 * (1 - eL)
    return (eL - r * e7) / den if den > 0 else float("inf")


def sf70_k(r, e7, eL):
    k = min(max(kappa(r, e7, eL), 0.5), 2.0)
    return (1 - k * (1 - e7)) / e7, k


def cats(sp, bn, L, data=False):
    a, b = V.wcat(sp, bn, 0.7), V.wcat(sp, bn, L, 0.7)
    if data and MODEL[sp] == "M8":
        a["kap"] *= KSC[sp]; b["kap"] *= KSC[sp]
    return a, b


def jf(H2, cs, sp, start=None, hesse=True):
    f = V.fit_set(H2, cs, sp, MODEL[sp], start=start, hesse=hesse)
    r = V.jres(f)
    r["x"] = f["x"]; r["p"] = F.pval(f["chi2"], f["ndf"])
    return r


def sig_weights(sp, bn, mode):
    """wagi sygnalu prawdy per przedzial score (19): 'zero' | 'step' (SF70 = 0,90) | 'kap' (kappa = 1,3 dla eps >= 0,5)."""
    ed = V.INFO["eps_desc"]
    eps = V.Z[f"{sp}|{V.TAG[sp]}|{bn}|eps"]
    j70, j50 = ed.index(0.7), ed.index(0.5)
    w = np.ones(19)
    if mode == "step":
        e7 = float(eps[j70]); w[j70 + 1:] = 0.9; w[:j70 + 1] = (1 - 0.9 * e7) / (1 - e7)
    elif mode == "kap":
        e5 = float(eps[j50]); w[:j50 + 1] = 1.3; w[j50 + 1:] = (1 - 1.3 * (1 - e5)) / e5
    return w


def t_mc(t):
    sp, bn, L = t
    V._init()
    a, b = cats(sp, bn, L)
    key = a["key"]
    e7, eL = a["e"], a["e"] + b["e"]
    f_mc = e7 / eL
    all_, sig = V.Z[key + "|mb|all"].astype(float), V.Z[key + "|mb|sig"].astype(float)
    res = {"sp": sp, "bin": bn, "L": L, "e7": e7, "eL": eL, "f_mc": f_mc}
    rng = np.random.default_rng(4242 + int(L * 100) + 7 * F.BINS.index(bn) + (300 if sp == "D0" else 0))
    seeds = rng.integers(0, 2 ** 31, NREP)
    out = {}
    for mode in ("zero", "step", "kap"):
        w = sig_weights(sp, bn, mode)
        full = all_ + (w[:, None] - 1.0) * sig
        sw = w[:, None] * sig
        HA, HB = V.full5(V.wsel(full, *a["j"])), V.full5(V.wsel(full, *b["j"]))
        nA, nB = float(V.wsel(sw, *a["j"]).sum()), float(V.wsel(sw, *b["j"]).sum())
        f_true = nA / (nA + nB)
        tot = float(sw.sum())
        o = {"f_true": f_true, "SF70_true": (nA / tot) / (float(V.wsel(sig, *a["j"]).sum()) / float(sig.sum())), "r_true": f_true / (float(V.wsel(sig, *a["j"]).sum()) / float(V.wsel(sig, *a["j"]).sum() + V.wsel(sig, *b["j"]).sum())),
             "NA_true": nA, "NB_true": nB}
        r = jf(np.stack([HA, HB]), [a, b], sp)
        o.update({k: r[k] for k in ("f", "f_err", "NA", "NB", "valid", "at_limit", "p", "sc")})
        fr = []
        for s in seeds:
            g = np.random.default_rng(int(s))
            rr = jf(np.stack([g.poisson(HA).astype(float), g.poisson(HB).astype(float)]), [a, b], sp, start=r["x"], hesse=False)
            fr.append(rr["f"] if rr["valid"] and rr["f"] is not None else np.nan)
        o["f_rep"] = fr
        out[mode] = o
    res["modes"] = out
    return res


def t_data(t):
    sp, bn, L, run = t
    V._init()
    a, b = cats(sp, bn, L, data=True)
    HA, HB = V.wdata(a, SRC[run]), V.wdata(b, SRC[run])
    r = jf(np.stack([HA[0], HB[0]]), [a, b], sp)
    res = {"sp": sp, "bin": bn, "L": L, "run": run, "e7": a["e"], "eL": a["e"] + b["e"], "f_mc": a["e"] / (a["e"] + b["e"]), **{k: r.get(k) for k in ("f", "f_err", "NA", "NB", "NA_err", "NB_err", "valid", "at_limit", "p", "sc", "dmu", "chi2", "ndf")}}
    fr = []
    for i in range(1, HA.shape[0]):
        rr = jf(np.stack([HA[i], HB[i]]), [a, b], sp, start=r["x"], hesse=False)
        fr.append(rr["f"] if rr["valid"] and rr["f"] is not None else np.nan)
    res["f_rep"] = fr
    res["n_rows"] = [float(HA[0].sum()), float(HB[0].sum())]
    return res


def t_toys(t):
    sp, bn, run, d, ntoy = t
    V._init()
    L = 0.9
    a, b = cats(sp, bn, L, data=True)
    c, S = F.SP[sp], V.SET[MODEL[sp]]
    rng_ = c[S["rng"]]
    edges = np.arange(rng_[0], rng_[1] + F.BW / 2, F.BW)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    blind = np.abs(ctr - c["m0"]) < c["win"]
    HA, HB = V.cutr(V.wdata(a, SRC[run])[0], sp, rng_), V.cutr(V.wdata(b, SRC[run])[0], sp, rng_)
    bk = []
    for H, k in ((HA, a), (HB, b)):
        f = F.fit(H[None], [k["sigT"]], [[]], edges, c["m0"], c["nbkg"] + S["dn"], blind=blind, win=c["win"], hesse=False)
        bk.append(f["mu"][0])
    shp = []
    for k in (a, b):
        s = k["sigT"].bins(edges, c["m0"])
        if k["refT"] is not None:
            s = s + k["frefl"] * k["refT"].bins(edges, c["m0"])
        if S["peak"]:
            s = s + k["kap"] * k["peakT"].bins(edges, c["m0"])
        shp.append(s)
    e7, eL = a["e"], a["e"] + b["e"]
    ntot = d["NA"] + d["NB"]
    res = {"sp": sp, "bin": bn, "run": run, "N": ntot}
    for lab, kap in (("k1.0", 1.0), ("k1.3", 1.3)):
        fA = (1 - kap * (1 - e7)) / (1 - kap * (1 - eL))
        f_true = fA
        rng = np.random.default_rng(99 + int(kap * 10) + 17 * F.BINS.index(bn) + (1000 if sp == "D0" else 0) + int(run) % 97)
        ff, fe = [], []
        for i in range(ntoy):
            T = np.stack([rng.poisson(bk[0] + ntot * f_true * shp[0]), rng.poisson(bk[1] + ntot * (1 - f_true) * shp[1])]).astype(float)
            fz = F.fit(T, [a["sigT"], b["sigT"]], [([(k["refT"], k["frefl"])] if k["refT"] is not None else []) + ([(k["peakT"], k["kap"])] if S["peak"] else []) for k in (a, b)],
                       edges, c["m0"], c["nbkg"] + S["dn"], shape=S["shape"], gc=V.GC[sp], win=c["win"])
            fz.pop("mu", None)
            r = V.jres(fz)
            if r["valid"] and r.get("f_err"):
                ff.append(r["f"]); fe.append(r["f_err"])
        ff, fe = np.array(ff), np.array(fe)
        pull = (ff - f_true) / fe
        res[lab] = {"f_true": f_true, "n": len(ff), "mean": float(ff.mean()), "rms": float(ff.std(ddof=1)), "hesse": float(fe.mean()), "pull_mean": float(pull.mean()), "pull_rms": float(pull.std(ddof=1)),
                    "bias_rel": float(ff.mean() / f_true - 1), "bias_rel_err": float(ff.std(ddof=1) / f_true / math.sqrt(len(ff)))}
    return res


def main():
    cmd, out = sys.argv[1:3]
    V._init()
    sps = ("Dplus", "D0")
    if cmd == "mc":
        tasks = [(sp, bn, L) for sp in sps for bn in F.BINS for L in LS]; fn = t_mc
    elif cmd == "data":
        tasks = [(sp, bn, L, run) for sp in sps for bn in F.BINS for L in LS for run in SRC]; fn = t_data
    else:
        D = {(r["sp"], r["bin"], r["run"]): r for r in json.load(open(sys.argv[3])) if r["L"] == 0.9}
        tasks = [(sp, bn, run, D[(sp, bn, run)], 40) for sp in sps for bn in F.BINS for run in SRC]; fn = t_toys
    res = []
    with Pool(int(sys.argv[4]) if len(sys.argv) > 4 and cmd != "toys" else 5) as p:
        for i, r in enumerate(p.imap_unordered(fn, tasks)):
            res.append(r)
            print(f"[{cmd}] {i + 1}/{len(tasks)} {r['sp']} {r['bin']} {r.get('L', '')} {r.get('run', '')}", flush=True)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
