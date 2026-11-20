#!/usr/bin/env python3
"""CATFIT — wspolne: kategorie score (przedzialy 0..18 histogramow WPSCAN), szablony, ustawienia modeli, fit K kategorii. POST HOC. Bez AFP; tylko piki inkluzywne.
Przedzial i: 0 = ponizej WP95 (eps 100-95 %), i = 1..18: [WP(95-5(i-1)), WP(95-5i)), 18 = powyzej WP10. Skumulowany WPe = przedzialy ed.index(e)+1 .. 18."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "..")
CW = os.path.join(BASE, "..")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(CW, "BKGCAT_2026-11-20", "skrypty"))
sys.path.insert(0, os.path.join(CW, "FITBIAS_2026-11-19", "skrypty"))
import bc_lib as L  # noqa: E402
import fb_fit as F  # noqa: E402
import cf_fit as CF  # noqa: E402

SPS, BINS = ["D0", "Dplus"], F.BINS
TAG = L.TAG
SRC = "427929|MB427929"
# kategorie: nazwa -> (i0, i1) polotwarty zakres przedzialow
K4 = {"L": (0, 6), "C": (6, 10), "B": (10, 14), "A": (14, 19)}
K5 = {"L1": (0, 2), "L2": (2, 6), "C": (6, 10), "B": (10, 14), "A": (14, 19)}
K3 = {"C": (6, 10), "B": (10, 14), "A": (14, 19)}
WPC = {0.7: (6, 19), 0.6: (8, 19), 0.5: (10, 19), 0.4: (12, 19), 0.3: (14, 19), 0.8: (4, 19), 0.9: (2, 19), 1.0: (0, 19)}
SET = {"D0": {"M8": dict(shape="gauss", peak=True, nbkg=2, rng=(1.70, 2.05), kind="exp"),
              "M4": dict(shape="gauss", peak=True, nbkg=2, rng=(1.75, 2.00), kind="exp"),
              "M9": dict(shape="gauss", peak=True, nbkg=3, rng=(1.75, 2.00), kind="exp")},
       "Dplus": {"M0": dict(shape="free", peak=False, nbkg=3, rng=(1.70, 2.20), kind="exp"),
                 "M0-1": dict(shape="free", peak=False, nbkg=2, rng=(1.70, 2.20), kind="exp"),
                 "M0+1": dict(shape="free", peak=False, nbkg=4, rng=(1.70, 2.20), kind="exp"),
                 "M0pol": dict(shape="free", peak=False, nbkg=3, rng=(1.70, 2.20), kind="pol"),
                 "M0rng": dict(shape="free", peak=False, nbkg=3, rng=(1.75, 2.10), kind="exp"),
                 "A2": dict(shape="gauss", peak=True, nbkg=3, rng=(1.70, 2.20), kind="exp"),
                 "A2+1": dict(shape="gauss", peak=True, nbkg=4, rng=(1.70, 2.20), kind="exp"),
                 "A2-1": dict(shape="gauss", peak=True, nbkg=2, rng=(1.70, 2.20), kind="exp"),
                 "A2rng": dict(shape="gauss", peak=True, nbkg=3, rng=(1.75, 2.10), kind="exp"),
                 "A4": dict(shape="free", peak=True, nbkg=3, rng=(1.70, 2.20), kind="exp")}}
NOM = {"D0": "M8", "Dplus": "M0"}
VARS = {"D0": {"M8": ["M4", "M9"]}, "Dplus": {"M0": ["M0-1", "M0+1", "M0pol", "M0rng"], "A2": ["A2-1", "A2+1", "A2rng", "A4"]}}
GC = {"D0": (0.0008, 0.10), "Dplus": (0.0005, 0.10)}
HSC = {"D0": 1.3, "Dplus": 1.0}            # skala HESSE w regule bledu stat. (final.json)
Z = INFO = BK = None
_cache = {}


def init():
    global Z, INFO, BK
    if Z is None:
        Z = np.load(L.WPF, allow_pickle=False)
        INFO = json.loads(str(Z["meta"]))
        BK = L.BK(L.IN)


def full_edges(sp):
    return L.full_edges(sp)


def full5(h1):
    return L.h5(h1)


def cat(sp, bn, i0, i1, kscale=1.0):
    """Skladniki modelu kategorii = przedzialy score [i0, i1)."""
    init()
    ck = (sp, bn, i0, i1)
    if ck not in _cache:
        lo = F.SP[sp]["mr"][0]
        key = f"{sp}|{TAG[sp]}|{bn}"
        swa = Z[key + "|fd|sig_w"]
        sw = swa[i0:i1].sum(0)
        o = {"sigT": F.tmpl1(sw, lo), "refT": None, "frefl": 0.0, "key": key, "sl": (i0, i1), "e": float(sw.sum() / swa.sum()), "n_fd": int(Z[key + "|fd|sig_n"][i0:i1].sum())}
        if F.SP[sp]["refl"]:
            rw = Z[key + "|fd|refl_w"][i0:i1].sum(0)
            o["refT"], o["frefl"] = F.tmpl1(rw, lo), float(rw.sum() / sw.sum())
        d = BK.dense(sp, bn)
        pk, ns = 0.0, 0.0
        for k, h in d.items():
            g = L.group(k)
            if g in ("DMM", "CM"):
                pk = pk + h[i0:i1].sum(0)
            elif g == "Sig":
                ns += h[i0:i1].sum()
        o["peakT"] = F.TmplE(full_edges(sp), L.h5(pk))
        o["kap"] = float(pk.sum() / ns)
        o["npk"] = float(pk.sum())
        _cache[ck] = o
    o = dict(_cache[ck])
    o["kap"] *= kscale
    return o


def hist(c, src):
    """hist 5 MeV pelnego zakresu; dane: [nrep+1, nfull]."""
    i0, i1 = c["sl"]
    return full5(Z[f"{c['key']}|{src}"][..., i0:i1, :].sum(-2))


def cutr(h5, sp, rng):
    return L.cutr(h5, sp, rng)


def fitK(H5, cats, sp, mname, start=None, hesse=True, cons=None, fracs=None, nbkg_add=None, sc_start=None, keep_mu=False):
    """H5 [ncat, nfull]; cats lista slownikow cat(). nbkg_add: lista dodatkow rzedu per kategoria."""
    S = SET[sp][mname]
    rng = S["rng"]
    edges = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    tied = []
    for k in cats:
        t = [(k["refT"], k["frefl"])] if k["refT"] is not None else []
        if S["peak"]:
            t.append((k["peakT"], k["kap"]))
        tied.append(t)
    nb = S["nbkg"] if nbkg_add is None else [S["nbkg"] + a for a in nbkg_add]
    f = CF.fit(cutr(np.asarray(H5, float), sp, rng), [k["sigT"] for k in cats], tied, edges, F.SP[sp]["m0"], nb, kind=S["kind"], shape=S["shape"], gc=GC[sp], win=F.SP[sp]["win"],
               start=start, hesse=hesse, cons=cons, fracs=fracs, sc_start=sc_start)
    if not keep_mu:
        f.pop("mu", None)
    return f


def summ(f, names):
    """Skrot wyniku fitu K kategorii."""
    N = np.array(f["N"])
    o = {"N": dict(zip(names, f["N"])), "Ntot": f["Ntot"], "Ntot_err": f.get("Ntot_err"), "valid": f["valid"], "at_limit": f["at_limit"], "sc": f["sc"], "dmu": f["dmu"],
         "sc_err": f.get("sc_err"), "chi2": f["chi2"], "ndf": f["ndf"], "p": F.pval(f["chi2"], f["ndf"]), "fval": f["fval"],
         "p_cat": {n: F.pval(c, d) for n, c, d in zip(names, f["chi2_cat"], f["ndf_cat"])}, "chi2_cat": dict(zip(names, f["chi2_cat"])), "ndf_cat": f["ndf_cat"][0],
         "err": dict(zip(names, f["err"])) if "err" in f else None}
    T = N.sum()
    o["f"] = {n: float(v / T) for n, v in zip(names, N)} if T > 0 else None
    if "cov_N" in f and T > 0:
        C = np.array(f["cov_N"])
        fe = {}
        for i, n in enumerate(names):
            g = -N[i] / T ** 2 * np.ones(len(N)); g[i] += 1 / T
            fe[n] = float(np.sqrt(max(g @ C @ g, 0)))
        o["f_err"] = fe
        # ulamki skumulowane od najciasniejszej (ostatnia = A): E_j = sum(N[j:])/T
        E, Ee = {}, {}
        for j in range(1, len(N)):
            g = -N[j:].sum() / T ** 2 * np.ones(len(N)); g[j:] += 1 / T
            E[names[j]] = float(N[j:].sum() / T); Ee[names[j]] = float(np.sqrt(max(g @ C @ g, 0)))
        o["E"], o["E_err"] = E, Ee
        # suma bez pierwszej kategorii (N70 dla K4)
        o["Nsub"] = float(N[1:].sum()); o["Nsub_err"] = float(np.sqrt(max(C[1:, 1:].sum(), 0)))
    return o


def q(*a):
    return float(np.sqrt(sum(x * x for x in a)))
