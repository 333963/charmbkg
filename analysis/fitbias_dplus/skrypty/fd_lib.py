#!/usr/bin/env python3
"""FITBIAS_Dplus (20.11.2026, projekt B; POST HOC) — wspólne: ustawienia A0–A8 (KRYTERIA §2), szablony związane wprost z drzewa H (bkgcat.npz / połówki), fit.
Importuje bez zmian: FITBIAS_2026-11-19/skrypty/{fb_fit,fb_var}.py, BKGCAT_2026-11-20/skrypty/bc_lib.py. Tylko D+, WP70, piki inkluzywne, bez AFP."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "..")
CW = os.path.join(BASE, "..")
sys.path.insert(0, os.path.join(CW, "BKGCAT_2026-11-20", "skrypty"))
sys.path.insert(0, os.path.join(CW, "FITBIAS_2026-11-19", "skrypty"))
import bc_lib as L  # noqa: E402
import fb_fit as F  # noqa: E402
import fb_var as V  # noqa: E402

SP = "Dplus"
C = F.SP[SP]
BINS = F.BINS
GC = (0.0005, 0.10)
FE = V.full_edges(SP)                       # 1,70–2,20, 5 MeV
CTR = 0.5 * (FE[1:] + FE[:-1])
SET = {"A0": dict(shape="free", tied=(), nbkg=3, rng=(1.70, 2.20)),
       "A1": dict(shape="gauss", tied=(), nbkg=3, rng=(1.70, 2.20)),
       "A2": dict(shape="gauss", tied=("pk",), nbkg=3, rng=(1.70, 2.20)),
       "A3": dict(shape="gauss", tied=("pk",), nbkg=3, rng=(1.75, 2.10)),
       "A4": dict(shape="free", tied=("pk",), nbkg=3, rng=(1.70, 2.20)),
       "A5": dict(shape="gauss", tied=("pk", "pc"), nbkg=3, rng=(1.70, 2.20)),
       "A6": dict(shape="gauss", tied=(), nbkg=3, rng=(1.75, 2.10)),
       "A7": dict(shape="gauss", tied=("pk",), nbkg=2, rng=(1.78, 2.05)),
       "A8": dict(shape="gauss", tied=("pk",), nbkg=4, rng=(1.70, 2.20))}
SIMPLE = ["A0", "A1", "A6", "A4", "A2", "A8", "A3", "A7", "A5"]
FILES = {"full": L.IN, "A": os.path.join(BASE, "in", "bkgcat_halfA.npz"), "B": os.path.join(BASE, "in", "bkgcat_halfB.npz")}
M0_DATA = {"427929": (3547.0, 1651.0, 687.0), "428770": (4118.0, 3451.0, 6458.0)}   # znane uzyski M0 (FITBIAS/tabele/data.md) — tylko do czystości zabawek
_BK = {}


def bk(which):
    if which not in _BK:
        _BK[which] = L.BK(FILES[which])
    return _BK[which]


def comps(which, bn, elo=0.7, ehi=None):
    """→ {nazwa: hist 5 MeV pełnego zakresu}: Sig, DMM, CM_Ds, CM_oth, CM, PC, PC|<rodzic>, NoCharm, CombB, Comb, pk (= DMM + CM), all."""
    c = bk(which).cat(SP, bn, elo, ehi)
    o = {}
    z = np.zeros(len(FE) - 1)

    def add(k, h):
        o[k] = o.get(k, z) + h
    for k, h1 in c.items():
        h = L.h5(h1)
        g = L.group(k)
        p = k.split("|")
        add("all", h)
        if g == "Sig":
            add("Sig", h)
        elif g == "DMM":
            add("DMM", h); add("pk", h)
        elif g == "CM":
            add("CM", h); add("pk", h); add("CM_Ds" if p[1] == "431" else "CM_oth", h)
        elif g == "PC":
            add("PC", h); add("PC|" + p[1], h)
        elif g == "NoCharm":
            add("NoCharm", h); add("Comb", h)
        else:
            add("CombB", h); add("Comb", h)
    for k in ("Sig", "DMM", "CM", "CM_Ds", "CM_oth", "PC", "NoCharm", "CombB", "Comb", "pk", "all"):
        o.setdefault(k, z.copy())
    return o


def tmpls(cm):
    """Szablony związane z N_S z rozkładu `cm` (comps): {'pk': (TmplE, κ), 'pc': (TmplE, κ)}."""
    ns = cm["Sig"].sum()
    return {"pk": (F.TmplE(FE, cm["pk"]), float(cm["pk"].sum() / ns)), "pc": (F.TmplE(FE, cm["PC"]), float(cm["PC"].sum() / ns))}


def cutr(h, rng):
    i0, i1 = int(round((rng[0] - FE[0]) / F.BW)), int(round((rng[1] - FE[0]) / F.BW))
    return h[..., i0:i1]


def fitA(name, H5, sigT, T, start=None, hesse=True, shape=None, keep_mu=False, kscale=1.0):
    """H5: hist 5 MeV pełnego zakresu [100]. T z tmpls(). → dict fb_fit.fit (N = normalizacja szablonu sygnału na pełnym zakresie)."""
    S = SET[name]
    rng = S["rng"]
    edges = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    tied = [[(T[t][0], T[t][1] * kscale) for t in S["tied"]]]
    f = F.fit(cutr(np.asarray(H5, float), rng)[None], [sigT], tied, edges, C["m0"], S["nbkg"], shape=shape or S["shape"], gc=GC, win=C["win"], start=start, hesse=hesse)
    if not keep_mu:
        f.pop("mu", None)
    return f


def res(f):
    return {"N": f["N"][0], "err": (f.get("err") or [None])[0], "valid": f["valid"], "at_limit": f["at_limit"], "sc": f["sc"], "dmu": f["dmu"], "sc_err": f.get("sc_err"), "dmu_err": f.get("dmu_err"),
            "chi2": f["chi2"], "ndf": f["ndf"], "p": F.pval(f["chi2"], f["ndf"]), "bkg_win": f["bkg_win"][0]}


def sig_tmpl(bn):
    """Szablon sygnału fd_nd przy WP70 (jak fb_var.wcat) + słownik kategorii W (do wdata)."""
    k = V.wcat(SP, bn, 0.7)
    return k["sigT"], k


def smooth(h, sigT, nb=5):
    f = F.fit(np.asarray(h, float)[None], [sigT], [[]], FE, C["m0"], nb, shape="fix", win=C["win"], fixN0=True, hesse=False)
    return f["mu"][0]
