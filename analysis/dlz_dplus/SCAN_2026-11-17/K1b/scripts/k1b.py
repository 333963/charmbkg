"""K1b: zysk σ(N)/N H+DLZ vs H dla D+ z tłem exp(pol2/3/4), F-test rzędu, obciążenie i pull z pseudo-danych. KRYTERIA_K1b.md.
Wejście (odczyt): DATA_2026-11-17/runs/lxplus/ana/{res,rep}_RUN.*   Użycie: python k1b.py [--nrep 200] [--ntoy 200] [--nproc 6]
"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.stats import chi2 as CHI2, f as FD

HERE = os.path.dirname(os.path.abspath(__file__))
K1B = os.path.dirname(HERE)
SCAN = os.path.dirname(K1B)
sys.path.insert(0, os.path.join(SCAN, "scripts"))
import scan as SC  # noqa: E402  (fit_bkg_only, ścieżki; nie uruchamia main)
FL = SC.FL
DATA = SC.DATA
RUNS, BINS, TAGS = SC.RUNS, SC.BINS, SC.TAGS
WPS = [0.5, 0.7, 0.8]
ORD = [2, 3, 4]
NBIN = FL.NB
RES = {}


def res_of(run):
    if run not in RES:
        RES[run] = json.load(open(os.path.join(DATA, "runs", "lxplus", "ana", f"res_{run}.json")))
    return RES[run]


def deviance(h, lam):
    lam = np.asarray(lam)
    t = np.where(h > 0, h * np.log(np.maximum(h, 1e-300) / lam), 0.0)
    return float(2 * np.sum(lam - h + t))


def nominal(a):
    """fity danych rzędu 2,3,4 w obu wariantach kształtu + F-test."""
    run, b, t, wp = a
    res = res_of(run)
    k = SC.key(t, b, wp)
    h = np.asarray(res["hists"][k], float)
    sh = res["shapes"][k]
    out = {"lab": f"{run}|{b}|{t}|{wp}"}
    for fv, free in (("free", True), ("fix", False)):
        prev, fits = None, {}
        for o in ORD:
            st = None
            if prev is not None:
                st = dict(N=prev["N"], dmu=prev["dmu"], sc=prev["sc"], nb=prev["nb"], **{f"b{i+1}": x for i, x in enumerate(prev["b"])})
            r = FL.fit(h, sh, order=o, free_shape=free, start=st)
            npar = (4 if free else 2) + o
            lam = r.pop("lam")
            r.update(chi2=float(r["chi2_ndf"] * (NBIN - npar)), ndf=NBIN - npar, dev=deviance(h, lam), rel=r["N_err"] / max(r["N"], 1e-9))
            fits[o] = r
            prev = r
        ft = {}
        for o1, o2 in ((2, 3), (3, 4), (2, 4)):
            d = fits[o1]["chi2"] - fits[o2]["chi2"]
            F = d / (o2 - o1) / (fits[o2]["chi2"] / fits[o2]["ndf"])
            ft[f"{o1}->{o2}"] = dict(F=float(F), p=float(FD.sf(F, o2 - o1, fits[o2]["ndf"])) if F > 0 else 1.0,
                                     p_lr=float(CHI2.sf(max(fits[o1]["dev"] - fits[o2]["dev"], 0), o2 - o1)))
        if ft["2->3"]["p"] < 0.05:
            sel = 4 if ft["3->4"]["p"] < 0.05 else 3
        else:
            sel = 4 if ft["2->4"]["p"] < 0.05 else 2
        out[fv] = dict(fits={str(o): fits[o] for o in ORD}, ftest=ft, sel=sel)
    return out


def heavy(a):
    """bootstrap binowy (rzędy 2,3,4) + pseudo-dane (prawda exp(pol6), fit rzędu 2,3,4) dla jednej komórki, w wariancie kształtu fv."""
    run, b, t, wp, fv, S0, nrep, ntoy, seed, starts = a
    free = fv == "free"
    res = res_of(run)
    k = SC.key(t, b, wp)
    h = np.asarray(res["hists"][k], float)
    sh = res["shapes"][k]
    rng = np.random.default_rng(seed)
    out = {"lab": f"{run}|{b}|{t}|{wp}", "boot": {}, "toy": {}}
    # bootstrap binowy
    hb = rng.poisson(h, (nrep, NBIN)).astype(float)
    for o in ORD:
        st = starts[str(o)]
        out["boot"][str(o)] = np.array([FL.fit_fast(hb[i], sh, st, order=o, free_shape=free) for i in range(nrep)]).tolist()
    # prawda
    ntrue = S0 * sh["eps_mc"]
    sig = ntrue * np.diff(FL.sig_cdf(FL.E, sh, 0.0, 1.0))
    bA, vA, cA = SC.fit_bkg_only(h, SC.SB45, 6, sig)
    bC, vC, cC = SC.fit_bkg_only(h, SC.SB70, 6, sig)
    out["truth"] = dict(N_true=float(ntrue), A_valid=vA, A_chi2sb=cA, C_valid=vC, C_chi2sb=cC, B_win_A=float(bA[~SC.SB45].sum()), B_win_C=float(bC[~SC.SB45].sum()))
    toys = rng.poisson(sig + bA, (ntoy, NBIN)).astype(float)
    for o in ORD:
        asim = {}
        for nm, bk in (("A", bA), ("C", bC)):
            prev = None
            for oo in range(2, o + 1):   # dojście do rzędu o przez niższe (stabilny start)
                st = dict(N=ntrue, dmu=0.0, sc=1.0, nb=bk.sum(), b1=-0.5) if prev is None else dict(N=prev["N"], dmu=prev["dmu"], sc=prev["sc"], nb=prev["nb"], **{f"b{i+1}": x for i, x in enumerate(prev["b"])})
                prev = FL.fit(sig + bk, sh, order=oo, free_shape=free, start=st)
            asim[nm] = prev
        ra = asim["A"]
        start = [ra["N"], ra["dmu"] if free else 0.0, ra["sc"] if free else 1.0, ra["nb"]] + ra["b"]
        r = np.array([FL.fit_fast(toys[i], sh, start, order=o, free_shape=free) for i in range(ntoy)])
        ok = (r[:, 4] > 0) & np.isfinite(r[:, 0]) & (r[:, 1] > 0)
        x = r[ok, 0] / ntrue
        pull = (r[ok, 0] - ntrue) / r[ok, 1]
        n = int(ok.sum())
        out["toy"][str(o)] = dict(asimov_A=float(ra["N"] / ntrue), asimov_C=float(asim["C"]["N"] / ntrue), n_ok=n, f_ok=float(ok.mean()),
                                  b=float(x.mean()) if n else float("nan"), b_err=float(x.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan"),
                                  med=float(np.median(x)) if n else float("nan"), rms=float(x.std(ddof=1)) if n > 1 else float("nan"),
                                  pull_mean=float(pull.mean()) if n else float("nan"), pull_w=float(pull.std(ddof=1)) if n > 1 else float("nan"),
                                  pull_w_err=float(pull.std(ddof=1) / np.sqrt(2 * n)) if n > 1 else float("nan"),
                                  relerr_hesse=float(np.mean(r[ok, 1]) / ntrue) if n else float("nan"),
                                  pull_w_centered=float(np.std((r[ok, 0] - r[ok, 0].mean()) / r[ok, 1], ddof=1)) if n > 1 else float("nan"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrep", type=int, default=200)
    ap.add_argument("--ntoy", type=int, default=200)
    ap.add_argument("--nproc", type=int, default=6)
    a = ap.parse_args()
    cells = [(run, b, t, wp) for run in RUNS for b in BINS for t in TAGS for wp in WPS]
    with Pool(a.nproc) as p:
        nom = {r["lab"]: r for r in p.imap_unordered(nominal, cells)}
    print("fity nominalne gotowe", flush=True)
    out = {"nrep": a.nrep, "ntoy": a.ntoy, "nominal": nom, "cells": {}}
    jobs = []
    fitvar = {}
    for run in RUNS:
        res = res_of(run)
        for b in BINS:
            fv = "fix" if b == "3-5" else "free"
            if fv == "free" and any(nom[f"{run}|{b}|{t}|0.7"]["free"]["fits"][str(nom[f"{run}|{b}|{t}|0.7"]["free"]["sel"])]["at_limit"] for t in TAGS):
                fv = "fix"
            fitvar[f"{run}|{b}"] = fv
            F2 = res["fits"] if fv == "free" else res["fits_fixshape"]
            S0 = F2[SC.key("H", b, 0.7)]["N"] / res["shapes"][SC.key("H", b, 0.7)]["eps_mc"]
            for t in TAGS:
                for wp in WPS:
                    n = nom[f"{run}|{b}|{t}|{wp}"][fv]["fits"]
                    starts = {str(o): [n[str(o)]["N"], n[str(o)]["dmu"] if fv == "free" else 0.0, n[str(o)]["sc"] if fv == "free" else 1.0, n[str(o)]["nb"]] + n[str(o)]["b"] for o in ORD}
                    jobs.append((run, b, t, wp, fv, S0, a.nrep, a.ntoy, 20261118 + len(jobs), starts))
    out["fitvar"] = fitvar
    with Pool(a.nproc) as p:
        for i, r in enumerate(p.imap_unordered(heavy, jobs, chunksize=1)):
            out["cells"][r["lab"]] = r
            print("komórka", i, r["lab"], flush=True)
    json.dump(out, open(os.path.join(K1B, "out", "k1b.json"), "w"))
    print("zapisano", flush=True)


if __name__ == "__main__":
    os.makedirs(os.path.join(K1B, "out"), exist_ok=True)
    main()
