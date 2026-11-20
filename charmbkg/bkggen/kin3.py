"""Kinematyka 3-ciał D⁺→Kππ (S6, dla generatorów S2/S4): (m, pT, η, φ, s12, s13, cosθ*, φ*, ψ) ↔ pędy córek.

Konwencja (zgodna z kin.py S2 dla 2 ciał):
  - układ helicity D: z' = kierunek D w lab, x' w płaszczyźnie (D, wiązka), y' = z'×x';
  - s12 = m²(d0 d1), s13 = m²(d0 d2) — zmienne Dalitza (d0 = K, d1/d2 = π w kolejności indeksów h5);
  - (cosθ*, φ*) = kierunek d0 w spoczynku D względem (x', y', z');
  - ψ = obrót płaszczyzny rozpadu wokół kierunku d0 (0 ⇔ d1 w półpłaszczyźnie zawierającej z' po stronie +z').
Rozpad 3-ciał w spoczynku = 2 zmienne Dalitza + 3 kąty Eulera; dla tła bez spinu (D spin 0) kąty są jednorodne
w (cosθ*, φ*, ψ) tylko dla fazy — tło kombinatoryczne NIE (akceptancja). Generator modeluje je jak inne cechy.
Wszystkie wielkości w GeV. Dokładność odwrócenia ~1e-9 (test_candrot.py).
"""
from __future__ import annotations

import numpy as np

from ..constants import M_K as MK, M_PI as MPI
MASSES3 = (MK, MPI, MPI)


def vec(pt, eta, phi):
    return np.stack([pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)], -1)


def frame(pD):
    """Układ helicity: (x', y', z') dla każdego wiersza; pD (n, 3)."""
    z = pD / np.linalg.norm(pD, axis=-1, keepdims=True)
    beam = np.array([0.0, 0.0, 1.0])
    x = beam - (z @ beam)[:, None] * z
    nx = np.linalg.norm(x, axis=-1, keepdims=True)
    x = np.where(nx > 1e-12, x / np.maximum(nx, 1e-300), np.array([1.0, 0.0, 0.0]))
    y = np.cross(z, x)
    return x, y, z


def _boost(q, Eq, b, sign):
    """Boost Lorentza wektora (Eq, q) o prędkość b (n,3); sign=+1 spoczynek→lab, −1 lab→spoczynek."""
    b2 = (b ** 2).sum(-1)
    g = 1 / np.sqrt(1 - b2)
    bq = (b * q).sum(-1)
    qn = q + (((g - 1) * bq / np.maximum(b2, 1e-300)) + sign * g * Eq)[:, None] * b
    En = g * (Eq + sign * bq)
    return qn, En


def lab_to_eta_phi(q):
    pt = np.hypot(q[:, 0], q[:, 1])
    return pt, np.arcsinh(q[:, 2] / pt), np.arctan2(q[:, 1], q[:, 0])


def _perp_ref(e1, z):
    """Wektor jednostkowy ⟂ e1, najbliższy z' (gdy e1 ∥ z' → dowolny ⟂)."""
    r = z - (e1 * z).sum(-1, keepdims=True) * e1
    n = np.linalg.norm(r, axis=-1, keepdims=True)
    alt = np.cross(e1, np.array([1.0, 0.0, 0.0]))
    alt /= np.maximum(np.linalg.norm(alt, axis=-1, keepdims=True), 1e-300)
    return np.where(n > 1e-9, r / np.maximum(n, 1e-300), alt)


def angles3(pts, etas, phis, masses=MASSES3):
    """Córki w lab (listy 3 tablic) → dict: m, pt, eta, phi, s12, s13, cts, phis, psi."""
    q = [vec(pts[i], etas[i], phis[i]) for i in range(3)]
    E = [np.sqrt((q[i] ** 2).sum(-1) + masses[i] ** 2) for i in range(3)]
    P, Et = q[0] + q[1] + q[2], E[0] + E[1] + E[2]
    m = np.sqrt(np.maximum(Et ** 2 - (P ** 2).sum(-1), 0))
    s = lambda i, j: (E[i] + E[j]) ** 2 - ((q[i] + q[j]) ** 2).sum(-1)  # noqa: E731
    b = P / Et[:, None]
    x, y, z = frame(P)
    k1, _ = _boost(q[0], E[0], b, -1)
    k2, _ = _boost(q[1], E[1], b, -1)
    e1 = k1 / np.linalg.norm(k1, axis=-1, keepdims=True)
    cts = (e1 * z).sum(-1)
    ph = np.arctan2((e1 * y).sum(-1), (e1 * x).sum(-1))
    r = _perp_ref(e1, z)
    t = np.cross(e1, r)
    e2 = k2 - (k2 * e1).sum(-1, keepdims=True) * e1
    psi = np.arctan2((e2 * t).sum(-1), (e2 * r).sum(-1))
    pt, eta, phi = lab_to_eta_phi(P)
    return {"m": m, "pt": pt, "eta": eta, "phi": phi, "s12": s(0, 1), "s13": s(0, 2), "cts": cts, "phis": ph,
            "psi": psi}


def daughters3(m, pt, eta, phi, s12, s13, cts, phis, psi, masses=MASSES3):
    """Odwrotność angles3 → (pts, etas, phis) córek w lab (listy po 3 tablice)."""
    m1, m2, m3 = masses
    m, s12, s13 = (np.asarray(a, float) for a in (m, s12, s13))
    s23 = m ** 2 + m1 ** 2 + m2 ** 2 + m3 ** 2 - s12 - s13
    E1 = (m ** 2 + m1 ** 2 - s23) / (2 * m)
    E2 = (m ** 2 + m2 ** 2 - s13) / (2 * m)
    E3 = (m ** 2 + m3 ** 2 - s12) / (2 * m)
    p1, p2, p3 = (np.sqrt(np.clip(e ** 2 - mm ** 2, 0, None)) for e, mm in ((E1, m1), (E2, m2), (E3, m3)))
    c12 = np.clip((p3 ** 2 - p1 ** 2 - p2 ** 2) / np.maximum(2 * p1 * p2, 1e-300), -1, 1)
    s12a = np.sqrt(1 - c12 ** 2)
    pD = vec(pt, eta, phi)
    ED = np.sqrt((pD ** 2).sum(-1) + m ** 2)
    x, y, z = frame(pD)
    st = np.sqrt(np.clip(1 - cts ** 2, 0, 1))
    e1 = (st * np.cos(phis))[:, None] * x + (st * np.sin(phis))[:, None] * y + cts[:, None] * z
    r = _perp_ref(e1, z)
    t = np.cross(e1, r)
    e2 = np.cos(psi)[:, None] * r + np.sin(psi)[:, None] * t
    k1 = p1[:, None] * e1
    k2 = p2[:, None] * (c12[:, None] * e1 + s12a[:, None] * e2)
    k3 = -k1 - k2
    b = pD / ED[:, None]
    out = [], [], []
    for k, Ek in ((k1, E1), (k2, E2), (k3, E3)):
        ql, _ = _boost(k, Ek, b, +1)
        for lst, v in zip(out, lab_to_eta_phi(ql)):
            lst.append(v)
    return out


def dalitz_norm(m, s12, s13, masses=MASSES3):
    """Znormalizowane współrzędne Dalitza na [0,1]² (kwadratowy Dalitz: m'(s12), θ' w układzie 12) — pod flow."""
    m1, m2, m3 = masses
    m12 = np.sqrt(np.maximum(s12, 0))
    lo, hi = m1 + m2, m - m3
    mp = np.clip((m12 - lo) / np.maximum(hi - lo, 1e-12), 0, 1)
    # cosθ w układzie (12): kąt między d2 a d0 w spoczynku pary (12)
    s23 = m ** 2 + m1 ** 2 + m2 ** 2 + m3 ** 2 - s12 - s13
    E2 = (s12 - m1 ** 2 + m2 ** 2) / (2 * np.maximum(m12, 1e-12))
    E3 = (m ** 2 - s12 - m3 ** 2) / (2 * np.maximum(m12, 1e-12))
    p2 = np.sqrt(np.clip(E2 ** 2 - m2 ** 2, 0, None))
    p3 = np.sqrt(np.clip(E3 ** 2 - m3 ** 2, 0, None))
    c = (m2 ** 2 + m3 ** 2 + 2 * E2 * E3 - s23) / np.maximum(2 * p2 * p3, 1e-12)
    return mp, np.arccos(np.clip(c, -1, 1)) / np.pi
