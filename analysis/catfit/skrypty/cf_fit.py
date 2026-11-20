#!/usr/bin/env python3
"""CATFIT (20.11.2026, projekt B; POST HOC wzgledem XSEC_final) — wlasna kopia FITBIAS fb_fit.fit (oryginal bez zmian) rozszerzona o:
  * nbkg per kategoria (lista) — do wariantu "rzad + 2 tylko w kategorii L";
  * cons = [(indeksy_kategorii_licznika, srodek, sigma), ...]: kara Gaussa na E = sum(N_licznik)/sum(N) (ulamki skumulowane);
  * fracs = [f_0..f_{K-1}]: ulamki stale, N_c = N_tot * f_c (parametr NS0 = N_tot, pozostale NS zamrozone).
Bez cons/fracs i ze stalym nbkg model i minimalizacja sa identyczne z fb_fit.fit."""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(CW, "FITBIAS_2026-11-19", "skrypty"))
import fb_fit as F  # noqa: E402

BW, DMU_MAX, SC_LIM = F.BW, F.DMU_MAX, F.SC_LIM
bshape = F.bshape


def fit(H, sigT, tied, edges, m0, nbkg, kind="exp", shape="free", gc=None, start=None, blind=None, win=0.045, fixN0=False, hesse=True, cons=None, fracs=None, sc_start=None):
    from iminuit import Minuit
    H = np.asarray(H, float)
    ncat, nb = H.shape
    nbk = [nbkg] * ncat if np.isscalar(nbkg) else list(nbkg)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    x = 2 * (ctr - edges[0]) / (edges[-1] - edges[0]) - 1
    use = np.ones(nb, bool) if blind is None else ~blind
    off = [2]
    for c in range(ncat):
        off.append(off[-1] + 2 + nbk[c])
    iN = off[:-1]
    free_shape = shape in ("free", "gauss") and blind is None and not fixN0
    cache = {}
    fr = None if fracs is None else np.asarray(fracs, float)

    def sshape(c, dmu, sc):
        k = (c, dmu, sc)
        if k in cache:
            return cache[k]
        s = sigT[c].bins(edges, m0, dmu, sc)
        for T, kap in tied[c]:
            if T is not None and kap != 0:
                s = s + kap * T.bins(edges, m0, dmu, 1.0)
        if free_shape and len(cache) >= 4 * ncat:      # pamiec kilku ostatnich punktow (dmu, sc); wartosci funkcji bez zmian
            cache.clear()
        cache[k] = s
        return s

    def nsig(p):
        if fr is not None:
            return p[iN[0]] * fr
        return np.array([p[i] for i in iN])

    def model(p):
        dmu, sc = p[0], p[1]
        mu = np.empty((ncat, nb))
        ns = nsig(p)
        for c in range(ncat):
            q = p[off[c]:off[c + 1]]
            sh = bshape(x, q[2:], kind)
            mu[c] = ns[c] * sshape(c, dmu, sc) + q[1] * sh / sh[use].sum()
        return np.maximum(mu, 1e-9)

    def cost(p):
        mu = model(p)[:, use]
        n = H[:, use]
        v = 2.0 * float(np.sum(mu - n) + np.sum(np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / mu), 0.0)))
        if shape == "gauss":
            v += (p[0] / gc[0]) ** 2 + ((p[1] - 1.0) / gc[1]) ** 2
        if cons:
            ns = nsig(p)
            tot = ns.sum()
            if tot > 1e-6:
                for ix, cen, sig in cons:
                    v += ((ns[ix].sum() / tot - cen) / sig) ** 2
        return v
    names = ["dmu", "sc"]
    p0 = [0.0, 1.0 if sc_start is None else sc_start]
    lim = [(-DMU_MAX, DMU_MAX), SC_LIM]
    inwin = np.abs(ctr - m0) < 0.03
    s0s = []
    for c in range(ncat):
        ntot = float(H[c][use].sum())
        nwin = float(H[c][inwin & use].sum())
        nsb = float(H[c][~inwin & use].sum()) * inwin.sum() / max((~inwin & use).sum(), 1)
        s0 = 0.0 if (blind is not None or fixN0) else max(nwin - nsb, 0.0) * 1.2
        s0s.append(s0)
        names += [f"NS{c}", f"NB{c}"] + [f"a{c}_{k + 1}" for k in range(nbk[c])]
        p0 += [s0, max(ntot - s0, 1.0)] + [0.0] * nbk[c]
        lim += [(-5 * math.sqrt(nwin + 1) - 10, 1.5 * ntot + 20), (0.0, 2.0 * ntot + 20)] + [(-30.0, 30.0)] * nbk[c]
    if fr is not None:
        p0[iN[0]] = sum(s0s)
        lim[iN[0]] = (-1e5, 3.0 * float(H.sum()) + 20)
    if start is not None and len(start) == len(p0):
        p0 = list(start)
        if sc_start is not None:
            p0[1] = sc_start
    m = Minuit(cost, p0, name=names)
    m.errordef = 1.0
    m.strategy = 1
    for nme, l in zip(names, lim):
        m.limits[nme] = l
    m.fixed["dmu"] = True; m.fixed["sc"] = True
    if blind is not None or fixN0:
        for c in range(ncat):
            m.values[f"NS{c}"] = 0.0; m.fixed[f"NS{c}"] = True
    if fr is not None:
        for c in range(1, ncat):
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
    mu = model(v)
    ndf = int(use.sum() * ncat - m.nfit)
    pear_c = ((H - mu) ** 2 / mu)[:, use].sum(1)
    ns = nsig(v)
    out = {"N": [float(a) for a in ns], "NB": [float(v[i + 1]) for i in iN], "dmu": float(v[0]), "sc": float(v[1]),
           "valid": bool(good() and (cov is not None or not hesse)), "edm": float(m.fmin.edm), "fval": float(m.fval), "chi2": float(pear_c.sum()), "ndf": ndf, "x": v.tolist(),
           "chi2_cat": [float(a) for a in pear_c], "ndf_cat": [int(use.sum() - 2 - nbk[c]) for c in range(ncat)],
           "at_limit": bool(free_shape and (abs(v[0]) > 0.98 * DMU_MAX or v[1] < SC_LIM[0] * 1.01 or v[1] > SC_LIM[1] * 0.99)), "mu": mu}
    out["Ntot"] = float(ns.sum())
    if cov is not None:
        if fr is not None:
            e = math.sqrt(max(cov[iN[0], iN[0]], 0))
            out["Ntot_err"] = float(e)
            out["err"] = [float(e * f) for f in fr]
        else:
            C = cov[np.ix_(iN, iN)]
            out["cov_N"] = C.tolist()
            out["err"] = [float(math.sqrt(max(C[i, i], 0))) for i in range(ncat)]
            out["Ntot_err"] = float(math.sqrt(max(C.sum(), 0)))
        out["dmu_err"], out["sc_err"] = float(math.sqrt(max(cov[0, 0], 0))), float(math.sqrt(max(cov[1, 1], 0)))
    return out
