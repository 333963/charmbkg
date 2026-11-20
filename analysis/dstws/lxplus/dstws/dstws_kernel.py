#!/usr/bin/env python3
"""DSTWS — jądro: kombinacje (kandydat D0→Kπ) × (tor z drzewa Tracks) → Δm = m(Kππ_s) − m(Kπ). numba, gdy jest; bez numba ten sam kod
w czystym Pythonie (wolno — tylko testy).

Typy kombinacji (kolumna typ):
  0 RS   tor o ładunku PRZECIWNYM do K (D*+ → D0 π+, D0 → K− π+), to samo zdarzenie, tor nie jest córką tego D0
  1 WS   tor o ładunku TAKIM SAMYM jak K
  2 ROT  tor RS obrócony w φ o 2π·k/(NROT+1), k = 1..NROT (kolumna k)
  3 MIX  tor o ładunku RS z INNEGO zdarzenia (k = 1..NMIX: sąsiedzi w porządku krotności torów w tej samej porcji pliku)
Δm liczone z 4-pędu kandydata D0 (pT, η, φ, m z dopasowania wierzchołka) i toru w PV z masą pionu — NIE z ponownego dopasowania
wierzchołka D* (nominał CharmAnalysis); różnicę mierzy kolumna dm_nom (Δm nominalnej pary, gdy tor = wolny pion nominalnego D*).
Córki D0 znajdowane wśród torów po (pT, η, φ): |pT/pT_c − 1| < MATCH_REL i ΔR < MATCH_DR, najbliższe pT (trackId wskazuje inny kontener).
Flagi (bity): 1 = tor jest wolnym pionem NOMINALNEJ pary D0←D* tego D0 (nom ≥ 0 = indeks pary w tablicy par nominalnych),
              2 = tor jest wolnym pionem jakiegokolwiek nominalnego D* w zdarzeniu z prawdą |pdg| = 211 i rodzicem |413| (MC),
              4 = tor jest córką (K lub π) innego kandydata — nie liczone (rezerwa).
"""
import math

import numpy as np

try:
    from numba import njit
    HAVE_NUMBA = True
except Exception:  # noqa: BLE001
    HAVE_NUMBA = False

    def njit(*a, **k):
        if a and callable(a[0]):
            return a[0]
        return lambda f: f

M_PI = 0.13957
MATCH_DR = 0.10
MATCH_REL = 0.05
MATCH_NOM = 0.02    # dopasowanie toru do wolnego pionu nominalnego D* (pT po dopasowaniu wierzchołka odbiega zwykle < 1,5 %)
NROT = 5
NMIX = 5
T_RS, T_WS, T_ROT, T_MIX = 0, 1, 2, 3


@njit(cache=True)
def _dphi(a, b):
    d = a - b
    while d > math.pi:
        d -= 2 * math.pi
    while d < -math.pi:
        d += 2 * math.pi
    return d


@njit(cache=True)
def _dm(dpx, dpy, dpz, dE, dm0, ppt, peta, pphi):
    px = ppt * math.cos(pphi); py = ppt * math.sin(pphi); pz = ppt * math.sinh(peta)
    e = math.sqrt(px * px + py * py + pz * pz + M_PI * M_PI)
    m2 = (dE + e) ** 2 - (dpx + px) ** 2 - (dpy + py) ** 2 - (dpz + pz) ** 2
    return math.sqrt(m2 if m2 > 0 else 0.0) - dm0


@njit(cache=True)
def _match(tpt, teta, tphi, a, b, pt, eta, phi, max_dr, max_rel):
    """Indeks toru [a, b) najbliższego w pT (|pT/pt − 1| < max_rel, ΔR < max_dr) albo −1."""
    best = -1; br = max_rel
    if pt <= 0:
        return -1
    for j in range(a, b):
        r = abs(tpt[j] / pt - 1.0)
        if r < br:
            de = teta[j] - eta; dp = _dphi(tphi[j], phi)
            if de * de + dp * dp < max_dr * max_dr:
                br = r; best = j
    return best


@njit(cache=True)
def combine(tpt, teta, tphi, tq, off, d_ev, D, noff, NP, toff, TP, pm, dm_max, pt_min, nrot, cap):
    """D[i] = (pt, eta, phi, m, qK, k_pt, k_eta, k_phi, pi_pt, pi_eta, pi_phi). NP[noff[i]:noff[i+1]] = (pt, eta, phi) wolnych pionów par
    nominalnych D0 i; TP[toff[e]:toff[e+1]] = (pt, eta, phi) prawdziwych wolnych pionów D* w zdarzeniu e; pm[e, k] = zdarzenie-partner
    do mieszania (−1 brak). Zwraca (n, c_d0, c_typ, c_k, c_dm, c_ppt, c_peta, c_dr, c_flag, c_nom, nmatch[nD0]); n > cap → powtórzyć z cap ≥ n."""
    nd = D.shape[0]
    c_d0 = np.empty(cap, np.int32); c_typ = np.empty(cap, np.int8); c_k = np.empty(cap, np.int8)
    c_dm = np.empty(cap, np.float32); c_ppt = np.empty(cap, np.float32); c_peta = np.empty(cap, np.float32)
    c_dr = np.empty(cap, np.float32); c_flag = np.empty(cap, np.int8); c_nom = np.empty(cap, np.int32)
    nmatch = np.zeros(nd, np.int8)
    n = 0
    for i in range(nd):
        e = d_ev[i]; a = off[e]; b = off[e + 1]
        pt = D[i, 0]; eta = D[i, 1]; phi = D[i, 2]; m0 = D[i, 3]; qk = D[i, 4]
        dpx = pt * math.cos(phi); dpy = pt * math.sin(phi); dpz = pt * math.sinh(eta)
        dE = math.sqrt(dpx * dpx + dpy * dpy + dpz * dpz + m0 * m0)
        jk = _match(tpt, teta, tphi, a, b, D[i, 5], D[i, 6], D[i, 7], MATCH_DR, MATCH_REL)
        jp = _match(tpt, teta, tphi, a, b, D[i, 8], D[i, 9], D[i, 10], MATCH_DR, MATCH_REL)
        nmatch[i] = (1 if jk >= 0 else 0) + (1 if jp >= 0 else 0)
        for j in range(a, b):
            if j == jk or j == jp or tpt[j] < pt_min:
                continue
            rs = tq[j] * qk < 0
            x = _dm(dpx, dpy, dpz, dE, m0, tpt[j], teta[j], tphi[j])
            if x < dm_max:
                fl = 0; nom = -1
                for q in range(noff[i], noff[i + 1]):
                    if abs(tpt[j] / NP[q, 0] - 1.0) < MATCH_NOM:
                        de = teta[j] - NP[q, 1]; dp = _dphi(tphi[j], NP[q, 2])
                        if de * de + dp * dp < MATCH_DR * MATCH_DR:
                            fl |= 1; nom = q
                            break
                for q in range(toff[e], toff[e + 1]):
                    if abs(tpt[j] / TP[q, 0] - 1.0) < MATCH_NOM:
                        de = teta[j] - TP[q, 1]; dp = _dphi(tphi[j], TP[q, 2])
                        if de * de + dp * dp < MATCH_DR * MATCH_DR:
                            fl |= 2
                            break
                if n < cap:
                    de = teta[j] - eta; dp = _dphi(tphi[j], phi)
                    c_d0[n] = i; c_typ[n] = T_RS if rs else T_WS; c_k[n] = 0; c_dm[n] = x; c_ppt[n] = tpt[j]; c_peta[n] = teta[j]
                    c_dr[n] = math.sqrt(de * de + dp * dp); c_flag[n] = fl; c_nom[n] = nom
                n += 1
            if rs:
                for k in range(1, nrot + 1):
                    ph = tphi[j] + 2 * math.pi * k / (nrot + 1)
                    x = _dm(dpx, dpy, dpz, dE, m0, tpt[j], teta[j], ph)
                    if x < dm_max:
                        if n < cap:
                            de = teta[j] - eta; dp = _dphi(ph, phi)
                            c_d0[n] = i; c_typ[n] = T_ROT; c_k[n] = k; c_dm[n] = x; c_ppt[n] = tpt[j]; c_peta[n] = teta[j]
                            c_dr[n] = math.sqrt(de * de + dp * dp); c_flag[n] = 0; c_nom[n] = -1
                        n += 1
        for k in range(pm.shape[1]):
            e2 = pm[e, k]
            if e2 < 0:
                continue
            for j in range(off[e2], off[e2 + 1]):
                if tpt[j] < pt_min or tq[j] * qk >= 0:
                    continue
                x = _dm(dpx, dpy, dpz, dE, m0, tpt[j], teta[j], tphi[j])
                if x < dm_max:
                    if n < cap:
                        de = teta[j] - eta; dp = _dphi(tphi[j], phi)
                        c_d0[n] = i; c_typ[n] = T_MIX; c_k[n] = k + 1; c_dm[n] = x; c_ppt[n] = tpt[j]; c_peta[n] = teta[j]
                        c_dr[n] = math.sqrt(de * de + dp * dp); c_flag[n] = 0; c_nom[n] = -1
                    n += 1
    return n, c_d0, c_typ, c_k, c_dm, c_ppt, c_peta, c_dr, c_flag, c_nom, nmatch


def mix_partners(ntrk, nmix=NMIX):
    """pm[e, k]: zdarzenia o najbliższej krotności torów (sąsiedzi w porządku posortowanym, cyklicznie); < nmix + 1 zdarzeń → −1."""
    n = len(ntrk)
    pm = np.full((n, nmix), -1, np.int64)
    if n <= nmix:
        return pm
    order = np.argsort(ntrk, kind="stable")
    pos = np.empty(n, np.int64); pos[order] = np.arange(n)
    for k in range(nmix):
        pm[:, k] = order[(pos + k + 1) % n]
    return pm


def run(tpt, teta, tphi, tq, off, d_ev, D, noff, NP, toff, TP, pm, dm_max=0.170, pt_min=0.0, nrot=NROT):
    cap = max(1024, 64 * len(D))
    while True:
        r = combine(tpt, teta, tphi, tq, off, d_ev, D, noff, NP, toff, TP, pm, dm_max, pt_min, nrot, cap)
        if r[0] <= cap:
            n = r[0]
            return tuple(x[:n] for x in r[1:10]) + (r[10],)
        cap = int(r[0]) + 16


def _selftest():
    """Zdarzenie 0: D0 (pT 5, η 0,3, φ 1) z córkami K−(q=−1), π+; tory: K, π, wolny pion RS (+) blisko, tor WS (−) blisko, daleki tor.
    Zdarzenie 1: jeden tor (+) blisko kierunku D0 → MIX."""
    def p4(pt, eta, phi, m):
        return np.array([pt * math.cos(phi), pt * math.sin(phi), pt * math.sinh(eta), math.sqrt((pt * math.cosh(eta)) ** 2 + m * m)])
    k = (3.0, 0.35, 1.05); pi = (2.0, 0.22, 0.93)
    dpt, deta, dphi, dm0 = 5.0, 0.3, 1.0, 1.86484          # 4-pęd D0 zadany wprost (córki służą tylko do wycięcia torów)
    d = p4(dpt, deta, dphi, dm0)
    tr = [k + (-1.0,), pi + (1.0,), (0.38, deta + 0.01, dphi - 0.01, 1.0), (0.36, deta - 0.02, dphi + 0.02, -1.0), (1.0, -1.5, -2.0, 1.0),
          (0.40, deta, dphi + 0.01, 1.0)]
    T = np.array(tr)
    off = np.array([0, 5, 6])
    D = np.array([[dpt, deta, dphi, dm0, -1.0, k[0], k[1], k[2], pi[0], pi[1], pi[2]]])
    NP = np.array([[0.38, deta + 0.01, dphi - 0.01]]); noff = np.array([0, 1])
    TP = np.array([[0.38, deta + 0.011, dphi - 0.012]]); toff = np.array([0, 1, 1])
    pm = np.array([[1], [0]])
    c_d0, typ, kk, dm, ppt, peta, dr, fl, nom, nmatch = run(T[:, 0].copy(), T[:, 1].copy(), T[:, 2].copy(), T[:, 3].copy(), off, np.array([0]), D,
                                                           noff, NP, toff, TP, pm)
    assert nmatch[0] == 2, nmatch
    rs = np.nonzero(typ == T_RS)[0]; ws = np.nonzero(typ == T_WS)[0]; mx = np.nonzero(typ == T_MIX)[0]
    assert len(rs) == 1 and len(ws) == 1 and len(mx) == 1, (typ, dm)
    assert fl[rs[0]] == 3 and nom[rs[0]] == 0 and fl[ws[0]] == 0, fl
    # niezależne przeliczenie Δm
    x = p4(0.38, deta + 0.01, dphi - 0.01, M_PI) + d
    ref = math.sqrt(x[3] ** 2 - x[0] ** 2 - x[1] ** 2 - x[2] ** 2) - dm0
    assert abs(dm[rs[0]] - ref) < 1e-6 and M_PI < ref < 0.17, (dm[rs[0]], ref)
    assert not np.any(typ == T_ROT)       # obrót o ≥ 60° wyrzuca pion poza Δm < 0,17
    pmx = mix_partners(np.array([5, 9, 7, 30, 8, 6, 100]), 2)
    assert pmx.shape == (7, 2) and np.all(pmx >= 0) and np.all(pmx[:, 0] != np.arange(7))
    print("dstws_kernel selftest OK (numba=%s)" % HAVE_NUMBA)


if __name__ == "__main__":
    _selftest()
