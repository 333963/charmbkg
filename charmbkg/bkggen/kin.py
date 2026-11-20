"""Kinematyka 2-ciał D0→Kπ: (m, pT, η, cosθ*, φ*) ↔ pędy córek. Pozwala generować minimalny zestaw zmiennych
bez łamania ścisłych więzów kinematycznych (log p_K, log p_π, max|η| córek, |cosθ*| wyprowadzane)."""
import numpy as np
from ..constants import M_K as MK, M_PI as MPI


def _vec(pt, eta, phi):
    return np.stack([pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)], -1)


def _frame(pD):
    z = pD / np.linalg.norm(pD, axis=-1, keepdims=True)
    beam = np.array([0.0, 0.0, 1.0])
    x = beam - (z @ beam)[:, None] * z
    x /= np.linalg.norm(x, axis=-1, keepdims=True)
    y = np.cross(z, x)
    return x, y, z


def angles(m, pt, eta, phi, kpt, keta, kphi):
    """cosθ*, φ* kaonu w układzie spoczynkowym D (oś z' = kierunek D, x' w płaszczyźnie (D, wiązka))."""
    pD = _vec(pt, eta, phi)
    E = np.sqrt((pD ** 2).sum(-1) + m ** 2)
    k = _vec(kpt, keta, kphi)
    Ek = np.sqrt((k ** 2).sum(-1) + MK ** 2)
    b = pD / E[:, None]
    b2 = (b ** 2).sum(-1)
    g = 1 / np.sqrt(1 - b2)
    bk = (b * k).sum(-1)
    kst = k + (((g - 1) * bk / b2) - g * Ek)[:, None] * b
    x, y, z = _frame(pD)
    kn = np.linalg.norm(kst, axis=-1)
    return (kst * z).sum(-1) / kn, np.arctan2((kst * y).sum(-1), (kst * x).sum(-1))


def daughters(m, pt, eta, cts, phis):
    """Zwraca (p_K, η_K, p_π, η_π) w lab (φ_D = 0, wynik niezależny od φ_D)."""
    n = len(m)
    pD = _vec(pt, eta, np.zeros(n))
    E = np.sqrt((pD ** 2).sum(-1) + m ** 2)
    ps = np.sqrt(np.clip((m ** 2 - (MK + MPI) ** 2) * (m ** 2 - (MK - MPI) ** 2), 0, None)) / (2 * m)
    x, y, z = _frame(pD)
    st = np.sqrt(np.clip(1 - cts ** 2, 0, 1))
    u = (st * np.cos(phis))[:, None] * x + (st * np.sin(phis))[:, None] * y + cts[:, None] * z
    b = pD / E[:, None]
    b2 = (b ** 2).sum(-1)
    g = 1 / np.sqrt(1 - b2)
    out = []
    for sgn, mm in ((1, MK), (-1, MPI)):
        q = sgn * ps[:, None] * u
        Eq = np.sqrt(ps ** 2 + mm ** 2)
        bq = (b * q).sum(-1)
        ql = q + (((g - 1) * bq / b2) + g * Eq)[:, None] * b
        p = np.linalg.norm(ql, axis=-1)
        out += [p, np.arctanh(ql[:, 2] / p)]
    return tuple(out)
