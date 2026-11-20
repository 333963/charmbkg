"""XSEC_atlas — wspolne: wczytanie tabel INCAL (lokalne), zmienne ciec w stylu ATLAS z kolumn tabeli."""
import glob, json, os
import numpy as np

CW = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INC = os.path.join(CW, "INCAL_2026-11-19", "runs", "lxplus", "parts")
MK, MPI, MPHI = 0.493677, 0.13957039, 1.019461
M0 = {"D0": 1.86484, "Dplus": 1.86966}
BINS = [("3-5", 3, 5), ("5-8", 5, 8), ("8-30", 8, 30)]


def load(sample, pref, mode="d150", cols=None):
    """Skleja tabele <pref>* (np. 'D0__sig_') ze wszystkich czesci probki."""
    out = {}
    for f in sorted(glob.glob(os.path.join(INC, f"inc_{mode}_{sample}_*.npz"))):
        z = np.load(f, allow_pickle=False)
        for k in z.files:
            if k.startswith(pref) and (cols is None or k[len(pref):] in cols):
                out.setdefault(k[len(pref):], []).append(z[k])
    return {k: np.concatenate(v) for k, v in out.items()}


def p4(pt, eta, phi, m):
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return np.stack([np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m ** 2), px, py, pz], -1)


def minv(*ps):
    s = sum(ps)
    return np.sqrt(np.maximum(s[..., 0] ** 2 - (s[..., 1:] ** 2).sum(-1), 0))


def cts(pd, pD):
    """cos kata miedzy corka w ukladzie spoczynkowym D a kierunkiem lotu D."""
    E, p = pD[..., 0], pD[..., 1:]
    M = np.sqrt(np.maximum(E ** 2 - (p ** 2).sum(-1), 1e-12))
    pm = np.sqrt((p ** 2).sum(-1))
    u = p / pm[..., None]
    pl = (pd[..., 1:] * u).sum(-1)
    g, bg = E / M, pm / M
    pls = g * pl - bg * pd[..., 0]
    pt2 = (pd[..., 1:] ** 2).sum(-1) - pl ** 2
    return pls / np.sqrt(np.maximum(pls ** 2 + pt2, 1e-12))


def derived(T, sp):
    """Dodaje: ptmin corek, maxdR, cts_K (z tabeli), D+: cts_pi_min, dphi (min |m(KK) - m_phi|)."""
    nd = 2 if sp == "D0" else 3
    pt = np.stack([T[f"d{i}_pt"] for i in range(nd)], 1).astype(float)
    eta = np.stack([T[f"d{i}_eta"] for i in range(nd)], 1).astype(float)
    phi = np.stack([T[f"d{i}_phi"] for i in range(nd)], 1).astype(float)
    T["ptmin"] = pt.min(1)
    dr = 0
    for i in range(nd):
        for j in range(i + 1, nd):
            dp = np.abs(phi[:, i] - phi[:, j]); dp = np.where(dp > np.pi, 2 * np.pi - dp, dp)
            dr = np.maximum(dr, np.hypot(eta[:, i] - eta[:, j], dp))
    T["maxdR"] = dr
    if sp == "Dplus":
        K = p4(pt[:, 0], eta[:, 0], phi[:, 0], MK)
        P = [p4(pt[:, i], eta[:, i], phi[:, i], MPI) for i in (1, 2)]
        D = K + P[0] + P[1]
        T["cts_pi_min"] = np.minimum(cts(P[0], D), cts(P[1], D))
        T["cts_K_chk"] = cts(K, D)
        KK = [p4(pt[:, i], eta[:, i], phi[:, i], MK) for i in (1, 2)]
        T["dphi"] = np.minimum(np.abs(minv(K, KK[0]) - MPHI), np.abs(minv(K, KK[1]) - MPHI))
    return T
