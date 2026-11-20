"""PEAK_pos: binned extended Poisson NLL (iminuit). Sygnał = 2 Gaussy o wspólnym μ, kształt (s1, s2, f) z MC,
w danych wolne μ i skala k (σ_i → k·σ_i). Tło: Czebyszew stopnia n (cheb) albo exp(a·x)·(1+b·x) (exp1)."""
import numpy as np
from iminuit import Minuit
from scipy.special import erf
from numpy.polynomial import chebyshev as Ch


def gcdf(x, mu, s):
    return 0.5 * (1 + erf((x - mu) / (np.sqrt(2) * s)))


def dg_binfrac(e, mu, s1, s2, f):
    return f * np.diff(gcdf(e, mu, s1)) + (1 - f) * np.diff(gcdf(e, mu, s2))


def seff(s1, s2, f):
    return float(np.sqrt(f * s1**2 + (1 - f) * s2**2))


def fit_mc(m, lo, hi, bw=0.002, w=None):
    """MC prawda: 2 Gaussy, wszystko wolne (bez tła). Zwraca mu, s1, s2, f, błędy, σ_eff."""
    e = np.arange(lo, hi + 1e-9, bw)
    n, _ = np.histogram(m, e, weights=w)
    N0 = n.sum()
    mu0, sd = np.average(m[(m > lo) & (m < hi)]), 0.012

    def nll(N, mu, s1, s2, f):
        lam = N * dg_binfrac(e, mu, s1, s2, f) + 1e-12
        return float(np.sum(lam - n * np.log(lam)))
    M = Minuit(nll, N=N0, mu=mu0, s1=sd * 0.7, s2=sd * 2, f=0.6)
    for kk, vv in dict(N=(0, 3 * N0), mu=(lo, hi), s1=(0.002, 0.05), s2=(0.004, 0.15), f=(0.05, 0.98)).items():
        M.limits[kk] = vv
    M.errordef = 0.5; M.migrad(); M.hesse()
    v, er = M.values.to_dict(), M.errors.to_dict()
    if v["s1"] > v["s2"]:  # porządek: s1 = wąski
        v["s1"], v["s2"], v["f"] = v["s2"], v["s1"], 1 - v["f"]
        er["s1"], er["s2"] = er["s2"], er["s1"]
    # błąd σ_eff z propagacji (kowariancja)
    p = np.array([v["s1"], v["s2"], v["f"]])
    return dict(**{k: float(x) for k, x in v.items()}, **{"e_" + k: float(x) for k, x in er.items()},
                seff=seff(*p), ok=bool(M.valid))


def bkg_shape(x, lo, hi, kind, c):
    u = 2 * (x - lo) / (hi - lo) - 1
    if kind == "exp1":
        return np.exp(c[0] * u) * (1 + c[1] * u)
    return Ch.chebval(u, np.r_[1.0, c])


def fit_data(m, lo, hi, shape, bkg="cheb2", bw=0.002, mu0=1.865, k_fixed=None, mu_fixed=None, refl=None):
    """dane: N_s, μ, k (skala σ), N_b, parametry tła. shape = dict(s1, s2, f) z MC.
    refl = (tablica mas refleksji z MC, R = N_refl/N_sig) → składnik o kształcie z MC i normalizacji R·N_s (keep_both)."""
    e = np.arange(lo, hi + 1e-9, bw)
    n, _ = np.histogram(m, e)
    c = 0.5 * (e[1:] + e[:-1])
    kind = "exp1" if bkg == "exp1" else "cheb"
    nc = 2 if bkg == "exp1" else int(bkg[-1])
    s1, s2, f = shape["s1"], shape["s2"], shape["f"]
    if refl is not None:
        hr, _ = np.histogram(refl[0], e)
        hr = hr / max(np.sum((refl[0] >= lo) & (refl[0] < hi)), 1)   # udział refleksji w oknie
        hr = hr * refl[1]
    else:
        hr = 0.0
    names = ["Ns", "mu", "k", "Nb"] + [f"c{i}" for i in range(nc)]

    def nll(*p):
        Ns, mu, k, Nb = p[:4]
        b = bkg_shape(c, lo, hi, kind, p[4:])
        if np.any(b <= 0):
            return 1e12
        lam = Ns * (dg_binfrac(e, mu, k * s1, k * s2, f) + hr) + Nb * b / b.sum() + 1e-12
        if np.any(lam <= 0):
            return 1e12
        return float(np.sum(lam - n * np.log(lam)))
    # start tła: liniowe z krańców
    st = [max(n.sum() * 0.05, 10), mu0, 1.0, n.sum() * 0.95] + [0.0] * nc
    if kind == "cheb":
        st[4] = (n[-5:].mean() - n[:5].mean()) / max(n.mean(), 1) / 2
    else:
        st[4] = np.log(max(n[-5:].mean(), 1) / max(n[:5].mean(), 1)) / 2
    M = Minuit(nll, *st, name=names)
    M.errordef = 0.5
    M.limits["Ns"] = (-0.5 * n.sum(), n.sum()); M.limits["Nb"] = (0, 2 * n.sum())
    M.limits["mu"] = (mu0 - 0.03, mu0 + 0.03); M.limits["k"] = (0.4, 2.5)
    if k_fixed is not None:
        M.values["k"] = k_fixed; M.fixed["k"] = True
    if mu_fixed is not None:
        M.values["mu"] = mu_fixed; M.fixed["mu"] = True
    M.migrad(); M.migrad(); M.hesse()
    v, er = M.values.to_dict(), M.errors.to_dict()
    lam = None
    return dict(**{k: float(x) for k, x in v.items()}, **{"e_" + k: float(x) for k, x in er.items()},
                ok=bool(M.valid), nll=float(M.fval), n=int(n.sum()), seff=float(v["k"] * seff(s1, s2, f)),
                e_seff=float(er["k"] * seff(s1, s2, f)))
