"""PEAK_lxy010: fit Δm(D*) — kopia 1:1 funkcji fit_dm z PEAK_pos/fit_dstar.py (tamten plik wykonuje się przy imporcie).
2 Gaussy (wolny kształt gdy shape=None, inaczej wolne μ i skala k) + tło u^a·(1+c1·u[+c2·u²]), u = Δm − m_π."""
import sys
import numpy as np
from iminuit import Minuit
sys.path.insert(0, "/home/k/Documents/charm/charm_work/PEAK_pos")
from fitlib import dg_binfrac, seff  # noqa: E402
MPI, DM0 = 0.13957039, 0.14543


def fit_dm(x, lo=0.1396, hi=0.160, bw=0.0002, shape=None, nb=1):   # kopia 1:1 z PEAK_pos/fit_dstar.py
    e = np.arange(lo, hi + 1e-12, bw); n, _ = np.histogram(x, e); c = 0.5 * (e[1:] + e[:-1]); u = (c - MPI) * 1e3

    def bk(a, cs):
        return np.power(np.maximum(u, 1e-6), a) * (1 + sum(ci * u ** (i + 1) for i, ci in enumerate(cs)))
    if shape is None:
        names = ["Ns", "mu", "s1", "s2", "f", "Nb", "a"] + [f"c{i}" for i in range(nb)]

        def nll(*p):
            Ns, mu, s1, s2, f, Nb, a = p[:7]; b = bk(a, p[7:])
            if np.any(b <= 0): return 1e12
            lam = Ns * dg_binfrac(e, mu, s1, s2, f) + Nb * b / b.sum() + 1e-12
            return float(np.sum(lam - n * np.log(lam)))
        st = [n.sum() * 0.3, DM0, 0.0004, 0.001, 0.6, n.sum() * 0.7, 0.5] + [0.0] * nb
    else:
        names = ["Ns", "mu", "k", "Nb", "a"] + [f"c{i}" for i in range(nb)]

        def nll(*p):
            Ns, mu, k, Nb, a = p[:5]; b = bk(a, p[5:])
            if np.any(b <= 0): return 1e12
            lam = Ns * dg_binfrac(e, mu, k * shape["s1"], k * shape["s2"], shape["f"]) + Nb * b / b.sum() + 1e-12
            if np.any(lam <= 0): return 1e12
            return float(np.sum(lam - n * np.log(lam)))
        st = [n.sum() * 0.3, DM0, 1.0, n.sum() * 0.7, 0.5] + [0.0] * nb
    M = Minuit(nll, *st, name=names); M.errordef = 0.5
    M.limits["Ns"] = (0, n.sum()); M.limits["Nb"] = (0, 2 * n.sum()); M.limits["mu"] = (DM0 - 0.002, DM0 + 0.002); M.limits["a"] = (0, 3)
    if shape is None:
        M.limits["s1"] = (0.0001, 0.003); M.limits["s2"] = (0.0002, 0.006); M.limits["f"] = (0.05, 0.98)
    else:
        M.limits["k"] = (0.4, 2.5)
    M.migrad(); M.migrad(); M.hesse()
    v, er = M.values.to_dict(), M.errors.to_dict()
    r = {**{k: float(x) for k, x in v.items()}, **{"e_" + k: float(x) for k, x in er.items()}, "ok": bool(M.valid), "n": int(n.sum())}
    if shape is None:
        if r["s1"] > r["s2"]:
            r["s1"], r["s2"], r["f"] = r["s2"], r["s1"], 1 - r["f"]
        r["seff"] = seff(r["s1"], r["s2"], r["f"])
    else:
        r["seff"] = r["k"] * seff(shape["s1"], shape["s2"], shape["f"]); r["e_seff"] = r["e_k"] * seff(shape["s1"], shape["s2"], shape["f"])
    return r
