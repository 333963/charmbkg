"""K1c: szybka implementacja TEGO SAMEGO fitu co DATA_2026-11-17/scripts/fitlib.py (te same wzory, siatka 2001 pkt, granice, clip),
z ndtr zamiast scipy.stats.norm.cdf i z góry policzonymi potęgami t. Zgodność z fitlib.fit sprawdzana w k1c.py (log)."""
import numpy as np
from iminuit import Minuit
from scipy.special import ndtr

M0 = 1.86966
LO, HI, NB = 1.70, 2.20, 100
E = np.linspace(LO, HI, NB + 1)
C = 0.5 * (E[1:] + E[:-1])
_X = np.linspace(LO, HI, 2001)
_T = (_X - M0) / 0.1
_TP = np.stack([_T ** (i + 1) for i in range(6)])   # potęgi t^1..t^6
_DX = np.diff(_X)
WIN = 0.045


def sig_bins(sh, dmu, sc):
    mu = sh["mu"] + dmu
    s1 = sh["s"] * sc
    c = sh["f"] * ndtr((E - mu) / s1) + (1 - sh["f"]) * ndtr((E - mu) / (s1 * sh["k"]))
    return np.diff(c)


def bkg_bins(b):
    ex = np.dot(np.asarray(b, float), _TP[:len(b)])
    fb = np.exp(ex - ex.max())
    tr = 0.5 * (fb[1:] + fb[:-1]) * _DX
    s = tr.reshape(NB, 20).sum(1)
    return s / s.sum()


def fit(h, sh, order=2, free_shape=True, start=None, ncall=20000, retry=True):
    """jak fitlib.fit (bez hipotezy N=0). start: lista [N, dmu, sc, nb, b1..] lub None (domyślny start fitlib)."""
    h = np.asarray(h, float)
    names = ["N", "dmu", "sc", "nb"] + [f"b{i+1}" for i in range(order)]

    def lam_of(p):
        return np.clip(p[0] * sig_bins(sh, p[1], p[2]) + p[3] * bkg_bins(p[4:]), 1e-9, None)

    def nll(*p):
        lam = lam_of(p)
        return float(np.sum(lam - h * np.log(lam)))
    if start is None:
        sb = np.abs(C - M0) > 0.07
        sr = np.abs(C - M0) < 0.03
        start = [max(h[sr].sum() - h[sb].mean() * sr.sum(), 10.0), 0.0, 1.0, h.sum(), -0.5] + [0.0] * (order - 1)
    start = list(start) + [0.0] * (len(names) - len(start))
    m = Minuit(nll, *start[:len(names)], name=names)
    m.errordef = 0.5
    m.limits["nb"] = (0, None)
    m.limits["dmu"] = (-0.006, 0.006)
    m.limits["sc"] = (0.7, 1.5)
    if not free_shape:
        m.fixed["dmu"] = True
        m.fixed["sc"] = True
    m.migrad(ncall=ncall)
    if retry and not m.valid:
        m.migrad(ncall=ncall)
    m.hesse()
    p = [m.values[n] for n in names]
    lam = lam_of(p)
    pear = (h - lam) ** 2 / np.maximum(lam, 1)
    npar = m.nfit
    lim = bool(free_shape and (abs(p[1]) > 0.0059 or p[2] < 0.705 or p[2] > 1.495))
    win = np.abs(C - M0) < WIN
    cs = p[0] * sig_bins(sh, p[1], p[2])
    cb = p[3] * bkg_bins(p[4:])
    return dict(N=float(p[0]), N_err=float(m.errors["N"]), dmu=float(p[1]), sc=float(p[2]), nb=float(p[3]), b=[float(x) for x in p[4:]],
                chi2=float(pear.sum()), ndf=int(len(h) - npar), valid=bool(m.valid), at_limit=lim, SB=float(cs[win].sum() / max(cb[win].sum(), 1e-9)),
                par=[float(x) for x in p])


def fit_fast(h, sh, start, order=2, free_shape=True):
    """jak fitlib.fit_fast: (N, N_err, dmu, sc, ok)."""
    r = fit(h, sh, order=order, free_shape=free_shape, start=start, ncall=5000, retry=False)
    return r["N"], r["N_err"], r["dmu"], r["sc"], bool(r["valid"] and not r["at_limit"])
