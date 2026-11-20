#!/usr/bin/env python3
"""FITBIAS — moduł fitu. Kopia dn4_ana.fit_cats (Q_lxplus/dnet4/dn4_ana.py, bez zmian w oryginale) rozszerzona o:
  * shape = "free" | "fix" | "gauss" (ograniczenie gaussowskie dmu, sc: gc = (σ_dmu [GeV], σ_sc));
  * dowolną liczbę składowych związanych z sygnałem w kategorii: tied[c] = [(szablon, κ), ...]  → N_c·[sig + Σ κ·szablon] (przesuwane o dmu, bez skali);
  * fixN0 (fit samego tła), rząd tła 0..n.
Model nominalny (shape="free", tied = [(odbicie, f_refl)]) jest identyczny z dn4_ana.fit_cats."""
import math

import numpy as np

BW = 0.005
GRID = 0.002
DMU_MAX, SC_LIM = 0.006, (0.7, 1.5)
SP = {"D0": {"m0": 1.86484, "win": 0.055, "fit": (1.75, 2.00), "alt": (1.70, 2.05), "nbkg": 2, "refl": True, "mr": (1.70, 2.05)},
      "Dplus": {"m0": 1.86966, "win": 0.045, "fit": (1.70, 2.20), "alt": (1.75, 2.10), "nbkg": 3, "refl": False, "mr": (1.70, 2.20)}}
BINS = ["3-5", "5-8", "8-30"]


class Tmpl:
    """Jak dn4_ana.Tmpl: dystrybuanta na siatce 2 MeV z (m, w), PCHIP; całka po [lo, hi] = 1."""

    def __init__(self, m, w, lo, hi):
        from scipy.interpolate import PchipInterpolator
        self.g = np.arange(lo - 0.02, hi + 0.02 + GRID / 2, GRID)
        h, _ = np.histogram(m, self.g, weights=w)
        tot = h.sum()
        self.n = float(tot)
        c = np.concatenate([[0.0], np.cumsum(h) / tot]) if tot > 0 else np.linspace(0, 1, len(self.g))
        self.f = PchipInterpolator(self.g, c)

    def bins(self, edges, m0, dmu=0.0, sc=1.0):
        x = np.clip((edges - m0 - dmu) / sc + m0, self.g[0], self.g[-1])
        return np.diff(self.f(x))


class TmplE(Tmpl):
    """Szablon z histogramu o zadanych krawędziach (np. 5 MeV): PCHIP przez dystrybuantę w krawędziach; całka po całym histogramie = 1."""

    def __init__(self, edges, h):
        from scipy.interpolate import PchipInterpolator
        h = np.asarray(h, float)
        tot = h.sum()
        self.n = float(tot)
        self.g = np.asarray(edges, float)
        c = np.concatenate([[0.0], np.cumsum(h) / tot]) if tot > 0 else np.linspace(0, 1, len(self.g))
        self.f = PchipInterpolator(self.g, c)


def tmpl1(h, lo, mbw=0.001):
    """Szablon z histogramu 1 MeV (WPSCAN) — identycznie jak wps_fit.tmpl."""
    ctr = lo + (np.arange(len(h)) + 0.5) * mbw
    return Tmpl(ctr, h, lo, lo + len(h) * mbw)


def cheb(x, a):
    t0, t1, s = np.ones_like(x), x, np.zeros_like(x)
    for ak in a:
        s = s + ak * t1
        t0, t1 = t1, 2 * x * t1 - t0
    return s


def bshape(x, a, kind="exp"):
    pol = cheb(x, a)
    return np.exp(np.clip(pol, -50, 50)) if kind == "exp" else np.maximum(1.0 + pol, 1e-9)


def fit(H, sigT, tied, edges, m0, nbkg, kind="exp", shape="free", gc=None, start=None, blind=None, win=0.045, fixN0=False, hesse=True, sig_shape=None):
    """H [ncat, nb]; sigT [ncat]; tied [ncat] listy (szablon, κ). → dict jak fit_cats (+ pełne x, mu modelu).
    sig_shape (opcjonalnie): funkcja (c, dmu, sc) → kształt sygnału+związanych w binach (zastępuje sigT/tied; do szybkich zabawek)."""
    from iminuit import Minuit
    H = np.asarray(H, float)
    ncat, nb = H.shape
    ctr = 0.5 * (edges[1:] + edges[:-1])
    x = 2 * (ctr - edges[0]) / (edges[-1] - edges[0]) - 1
    use = np.ones(nb, bool) if blind is None else ~blind
    npc = 2 + nbkg
    free_shape = shape in ("free", "gauss") and blind is None and not fixN0
    # składowe związane przesuwają się tylko o dmu: przy shape="fix" wszystko stałe → policz raz
    cache = {}

    def sshape(c, dmu, sc):
        if sig_shape is not None:
            return sig_shape(c, dmu, sc)
        k = (c, dmu, sc)
        if not free_shape and k in cache:
            return cache[k]
        s = sigT[c].bins(edges, m0, dmu, sc)
        for T, kap in tied[c]:
            if T is not None and kap != 0:
                s = s + kap * T.bins(edges, m0, dmu, 1.0)
        if not free_shape:
            cache[k] = s
        return s

    def model(p):
        dmu, sc = p[0], p[1]
        mu = np.empty((ncat, nb))
        for c in range(ncat):
            q = p[2 + c * npc: 2 + (c + 1) * npc]
            sh = bshape(x, q[2:], kind)
            mu[c] = q[0] * sshape(c, dmu, sc) + q[1] * sh / sh[use].sum()
        return np.maximum(mu, 1e-9)

    def cost(p):
        mu = model(p)[:, use]
        n = H[:, use]
        v = 2.0 * float(np.sum(mu - n) + np.sum(np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / mu), 0.0)))
        if shape == "gauss":
            v += (p[0] / gc[0]) ** 2 + ((p[1] - 1.0) / gc[1]) ** 2
        return v
    names = ["dmu", "sc"]
    p0 = [0.0, 1.0]
    lim = [(-DMU_MAX, DMU_MAX), SC_LIM]
    inwin = np.abs(ctr - m0) < 0.03
    for c in range(ncat):
        ntot = float(H[c][use].sum())
        nwin = float(H[c][inwin & use].sum())
        nsb = float(H[c][~inwin & use].sum()) * inwin.sum() / max((~inwin & use).sum(), 1)
        s0 = 0.0 if (blind is not None or fixN0) else max(nwin - nsb, 0.0) * 1.2
        names += [f"NS{c}", f"NB{c}"] + [f"a{c}_{k + 1}" for k in range(nbkg)]
        p0 += [s0, max(ntot - s0, 1.0)] + [0.0] * nbkg
        lim += [(-5 * math.sqrt(nwin + 1) - 10, 1.5 * ntot + 20), (0.0, 2.0 * ntot + 20)] + [(-30.0, 30.0)] * nbkg
    if start is not None and len(start) == len(p0):
        p0 = list(start)
    m = Minuit(cost, p0, name=names)
    m.errordef = 1.0
    m.strategy = 1
    for nme, l in zip(names, lim):
        m.limits[nme] = l
    m.fixed["dmu"] = True; m.fixed["sc"] = True
    if blind is not None or fixN0:
        for c in range(ncat):
            m.values[f"NS{c}"] = 0.0; m.fixed[f"NS{c}"] = True

    def good():
        fm = m.fmin
        return fm is not None and not fm.has_reached_call_limit and (m.valid or fm.edm < 0.02) and not fm.is_above_max_edm
    m.migrad(ncall=20000)
    if free_shape:
        m.fixed["dmu"] = False; m.fixed["sc"] = False
        m.migrad(ncall=20000)
    for it in range(3):
        if good():
            break
        if it == 2:
            m.strategy = 2
        m.simplex(ncall=20000)
        m.migrad(ncall=40000)
    cov = None
    if hesse:
        try:
            m.hesse()
            fm = m.fmin
            if m.covariance is not None and fm.has_covariance and fm.has_posdef_covar and not fm.has_made_posdef_covar:
                cov = np.array(m.covariance)
        except Exception:  # noqa: BLE001
            cov = None
    v = np.array(m.values)
    iN = [2 + c * npc for c in range(ncat)]
    mu = model(v)
    ndf = int(use.sum() * ncat - m.nfit)
    pear = float(np.sum(((H - mu) ** 2 / mu)[:, use]))
    out = {"N": [float(v[i]) for i in iN], "NB": [float(v[i + 1]) for i in iN], "dmu": float(v[0]), "sc": float(v[1]),
           "valid": bool(good() and (cov is not None or not hesse)), "edm": float(m.fmin.edm), "fval": float(m.fval), "chi2": pear, "ndf": ndf, "x": v.tolist(),
           "at_limit": bool(free_shape and (abs(v[0]) > 0.98 * DMU_MAX or v[1] < SC_LIM[0] * 1.01 or v[1] > SC_LIM[1] * 0.99)), "mu": mu}
    if cov is not None:
        out["cov_N"] = cov[np.ix_(iN, iN)].tolist()
        out["err"] = [float(math.sqrt(max(cov[i, i], 0))) for i in iN]
        out["dmu_err"], out["sc_err"] = float(math.sqrt(max(cov[0, 0], 0))), float(math.sqrt(max(cov[1, 1], 0)))
    wmask = np.abs(ctr - m0) < win
    out["bkg_win"] = []
    out["sig_win_frac"] = []
    for ci, i in enumerate(iN):
        sh = bshape(x, v[i + 2:i + 2 + nbkg], kind)
        out["bkg_win"].append(float(v[i + 1] * sh[wmask].sum() / sh[use].sum()))
    return out


def rebin1(h, lo, flo, fhi, mbw=0.001):
    """hist 1 MeV od lo → biny 5 MeV w [flo, fhi] (jak wps_fit.rebin)."""
    k = int(round(BW / mbw))
    i0, i1 = int(round((flo - lo) / mbw)), int(round((fhi - lo) / mbw))
    x = h[..., i0:i1]
    return x.reshape(x.shape[:-1] + ((i1 - i0) // k, k)).sum(-1).astype(np.float64)


def pval(chi2, ndf):
    from scipy.stats import chi2 as C2
    return float(C2.sf(chi2, ndf))
