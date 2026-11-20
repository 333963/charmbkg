"""DLZ_dplus/DATA: fit masy D+ (KRYTERIA §3): 2 Gaussy o kształcie z MC (wolne N_S, przesunięcie dμ, skala szerokości) + exp(pol-k).
Binned Poisson ML, biny 5 MeV, okno 1,70–2,20 GeV."""
import numpy as np
from iminuit import Minuit
from scipy.stats import norm

M0 = 1.86966
LO, HI, NB = 1.70, 2.20, 100
E = np.linspace(LO, HI, NB + 1)
C = 0.5 * (E[1:] + E[:-1])
WIN = 0.045
_X = np.linspace(LO, HI, 2001)
_T = (_X - M0) / 0.1


def mc_shape(m, w=None):
    """podwójny Gauss dopasowany do sygnału MC (truth 0/1) w ±80 MeV → mu, s, k, f."""
    e = np.linspace(M0 - 0.08, M0 + 0.08, 81)
    h, _ = np.histogram(m, e, weights=w)

    def nll(n, mu, s, k, f):
        cdf = f * norm.cdf(e, mu, s) + (1 - f) * norm.cdf(e, mu, s * k)
        lam = np.clip(n * np.diff(cdf), 1e-9, None)
        return float(np.sum(lam - h * np.log(lam)))
    mm = Minuit(nll, n=h.sum(), mu=M0, s=0.010, k=2.0, f=0.7)
    mm.errordef = 0.5
    mm.limits["s"] = (0.003, 0.03)
    mm.limits["k"] = (1.2, 6)
    mm.limits["f"] = (0.05, 1)
    mm.migrad()
    v = {p: float(mm.values[p]) for p in ("mu", "s", "k", "f")}
    v["rms_eff"] = float(np.sqrt(v["f"] * v["s"] ** 2 + (1 - v["f"]) * (v["s"] * v["k"]) ** 2))
    return v


def sig_cdf(xe, sh, dmu, sc):
    mu = sh["mu"] + dmu
    return sh["f"] * norm.cdf(xe, mu, sh["s"] * sc) + (1 - sh["f"]) * norm.cdf(xe, mu, sh["s"] * sh["k"] * sc)


def bkg_cdf(xe, b, lo=LO, hi=HI):
    x = _X if (lo == LO and hi == HI) else np.linspace(lo, hi, 2001)
    t = (x - M0) / 0.1
    ex = sum(bi * t ** (i + 1) for i, bi in enumerate(b))
    fb = np.exp(ex - ex.max())
    cb = np.concatenate([[0], np.cumsum(0.5 * (fb[1:] + fb[:-1]) * np.diff(x))])
    return np.interp(xe, x, cb / cb[-1])


def fit(h, sh, order=2, edges=E, free_shape=True, start=None):
    """zwraca dict: N, N_err, dmu, sc (+błędy), NB, S/B w ±45 MeV, Z = √(2ΔlnL), χ²/ndf (całe okno i same SB)."""
    h = np.asarray(h, float)
    lo, hi = edges[0], edges[-1]
    names = ["N", "dmu", "sc", "nb"] + [f"b{i+1}" for i in range(order)]

    def lam_of(p):
        N, dmu, sc, nb = p[:4]
        cs = sig_cdf(edges, sh, dmu, sc)
        return np.clip(N * np.diff(cs) + nb * np.diff(bkg_cdf(edges, p[4:], lo, hi)), 1e-9, None)

    def nll(*p):
        lam = lam_of(p)
        return float(np.sum(lam - h * np.log(lam)))
    c = 0.5 * (edges[1:] + edges[:-1])
    sb = np.abs(c - M0) > 0.07
    n_sr_est = max(h[np.abs(c - M0) < 0.03].sum() - h[sb].mean() * (np.abs(c - M0) < 0.03).sum(), 10.0)
    st = start or dict(N=n_sr_est, dmu=0.0, sc=1.0, nb=h.sum(), b1=-0.5, b2=0.0, b3=0.0)
    m = Minuit(nll, *[st.get(n, 0.0) for n in names], name=names)
    m.errordef = 0.5
    m.limits["nb"] = (0, None)
    m.limits["dmu"] = (-0.006, 0.006)
    m.limits["sc"] = (0.7, 1.5)
    if not free_shape:
        m.fixed["dmu"] = True
        m.fixed["sc"] = True
    m.migrad(ncall=20000)
    if not m.valid:
        m.migrad(ncall=20000)
    m.hesse()
    p = [m.values[n] for n in names]
    lam = lam_of(p)
    nll1 = m.fval
    # hipoteza N = 0
    m0 = Minuit(nll, *[0.0 if n == "N" else (st.get(n, 0.0) if n in ("dmu", "sc") else m.values[n]) for n in names], name=names)
    m0.errordef = 0.5
    m0.fixed["N"] = True
    m0.fixed["dmu"] = True
    m0.fixed["sc"] = True
    m0.limits["nb"] = (0, None)
    m0.migrad(ncall=20000)
    z2 = 2 * (m0.fval - nll1)
    win = np.abs(c - M0) < WIN
    cb = np.diff(bkg_cdf(edges, p[4:], lo, hi)) * p[3]
    cs = np.diff(sig_cdf(edges, sh, p[1], p[2])) * p[0]
    pear = (h - lam) ** 2 / np.maximum(lam, 1)
    sbb = np.abs(c - M0) > WIN
    return dict(N=float(p[0]), N_err=float(m.errors["N"]), dmu=float(p[1]), dmu_err=float(m.errors["dmu"]), sc=float(p[2]), sc_err=float(m.errors["sc"]),
                nb=float(p[3]), b=[float(x) for x in p[4:]], S_win=float(cs[win].sum()), B_win=float(cb[win].sum()),
                SB=float(cs[win].sum() / max(cb[win].sum(), 1e-9)), Z=float(np.sign(p[0]) * np.sqrt(max(z2, 0))),
                chi2_ndf=float(pear.sum() / (len(h) - m.nfit)), chi2_sb=float(pear[sbb].sum()), n_sb=int(sbb.sum()), npar_b=order + 1,
                valid=bool(m.valid), n=float(h.sum()), lam=lam.tolist(),
                at_limit=bool(free_shape and (abs(p[1]) > 0.0059 or p[2] < 0.705 or p[2] > 1.495)))


def fit_fast(h, sh, start, order=2, edges=E, free_shape=True):
    """szybki fit do bootstrapu: start z nominału, bez fitu N=0; zwraca (N, N_err, dmu, sc)."""
    h = np.asarray(h, float)
    lo, hi = edges[0], edges[-1]
    names = ["N", "dmu", "sc", "nb"] + [f"b{i+1}" for i in range(order)]

    def nll(*p):
        lam = np.clip(p[0] * np.diff(sig_cdf(edges, sh, p[1], p[2])) + p[3] * np.diff(bkg_cdf(edges, p[4:], lo, hi)), 1e-9, None)
        return float(np.sum(lam - h * np.log(lam)))
    m = Minuit(nll, *start, name=names)
    m.errordef = 0.5
    m.limits["nb"] = (0, None)
    m.limits["dmu"] = (-0.006, 0.006)
    m.limits["sc"] = (0.7, 1.5)
    if not free_shape:
        m.fixed["dmu"] = True
        m.fixed["sc"] = True
    m.migrad(ncall=5000)
    m.hesse()
    lim = free_shape and (abs(m.values["dmu"]) > 0.0059 or m.values["sc"] < 0.705 or m.values["sc"] > 1.495)
    return float(m.values["N"]), float(m.errors["N"]), float(m.values["dmu"]), float(m.values["sc"]), bool(m.valid and not lim)
