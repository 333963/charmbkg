#!/usr/bin/env python3
"""TRKJET — jądro: dżety torowe anty-kT R = 0,4 (WŁASNA implementacja, schemat E, odległości w (y, φ); fastjet nie ma w q_env LCG_110_swan),
stożki i półkule wokół D0. numba, gdy jest (lxplus: 0.65); bez numba działa ten sam kod w czystym Pythonie (wolno — tylko testy).

Na parę D0←D* liczone są dwie kolekcje dżetów z torów pT > 0,5 GeV, |η| < 2,5 (masa pionu):
  A  — D0 jako JEDNA pseudocząstka (pT, η, φ, m kandydata) zamiast dopasowanych córek; reszta torów (wolny pion zostaje, jeśli przechodzi
       cięcia). Dżet zawierający D0 jest więc zawsze określony → z = pT(D0)/pT(dżetu z D0).
  B  — tory BEZ córek D0 i BEZ wolnego pionu („otoczenie”): dżet obok D0 nie zawiera samego D0 → zmienna niezależna od pT(D0).
  O  — jak B, ale tylko tory POZA stożkiem izolacji ΔR(tor, D0) > 0,4 → zmienne dżetowe niezależne od ptcone40 (bez błędnego koła).
Córki/wolny pion znajdowane w drzewie Tracks po (pT, η, φ): najbliższy tor z ΔR < MATCH_DR i |pT/pT_córki − 1| < MATCH_REL
(trackId z DMesons wskazuje inny kontener — V20_FORMAT §3). Kolumny wyjścia: NAMES (float32).
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

R_JET = 0.4
PT_MIN = 0.5
ETA_MAX = 2.5
MATCH_DR = 0.10     # kierunek córki po dopasowaniu wierzchołka odbiega od toru w PV (ogon do ~0,1 dla pT < 0,5 GeV; test 17.11)
MATCH_REL = 0.05    # pT zgodne zwykle do 10⁻⁵
M_PI = 0.13957
NONE_DR = 99.0

NAMES = [
    # kolekcja A (D0 jako pseudocząstka)
    "a_n5", "a_j1_pt", "a_j1_eta", "a_j1_phi", "a_j2_pt", "a_j2_eta", "a_j2_phi", "a_dr_j1",
    "a_jD_pt", "a_jD_eta", "a_jD_phi", "a_jD_n", "a_dr_jD", "a_jD_rank", "z", "a_drn5", "a_ptn5",
    # kolekcja B (bez córek D0 i wolnego pionu)
    "b_n5", "b_j1_pt", "b_j1_eta", "b_j1_phi", "b_j2_pt", "b_j2_eta", "b_j2_phi", "b_dr_j1",
    "b_drn2", "b_ptn2", "b_drn3", "b_ptn3", "b_drn5", "b_ptn5", "b_drn8", "b_ptn8", "b_ptmax04", "b_ptmax08",
    # stożki wokół D0: tory pT > 0,5, |η| < 2,5, bez córek D0 i wolnego pionu
    "c_n02", "c_s02", "c_n04", "c_s04", "c_n10", "c_s10",
    # walidacja ptcone (R = 0,4 / 1,0): p2 = wszystkie tory drzewa (pT > 0,2), p5 = pT > 0,5 (bez cięcia η); bez córek D0;
    # „x” = dodatkowo bez wolnego pionu
    "v_s40_p2", "v_s40_p2x", "v_s40_p5", "v_s40_p5x", "v_s100_p2", "v_s100_p5",
    # półkule względem D0 (tory pT > 0,5, |η| < 2,5, bez córek i wolnego pionu): away |Δφ| > π/2, tow |Δφ| < π/2
    "h_s_away", "h_n_away", "h_s_tow", "h_n_tow",
    # wiodący tor spoza córek (te same cięcia) i ΔR do D0; liczba torów po cięciach (z córkami)
    "t1_pt", "t1_dr", "ntrk05",
    # dopasowanie córek: liczba dopasowanych córek D0 (0–2), wolny pion (0/1), max ΔR dopasowania;
    # nd_lo = liczba córek D0 z pT < 0,2 GeV (nie ma ich w drzewie Tracks → nie wymagają wycięcia; komplet: nmatch + nd_lo == 2)
    "nmatch", "pis_match", "match_dr_max", "nd_lo",
    # kolekcja O („outside”): dżety z torów kolekcji B leżących POZA stożkiem izolacji (ΔR(tor, D0) > 0,4) → zmienne dżetowe
    # NIEZALEŻNE od zawartości stożka ptcone40 (test hipotezy bez błędnego koła): wiodący, najbliższy z pT > 3 / 5, max pT w ΔR < 0,8
    "o_n5", "o_j1_pt", "o_dr_j1", "o_drn3", "o_ptn3", "o_drn5", "o_ptn5", "o_ptmax08",
]
NOUT = len(NAMES)
IX = {n: i for i, n in enumerate(NAMES)}


@njit(cache=False)
def _dphi(a, b):
    d = a - b
    while d > math.pi:
        d -= 2.0 * math.pi
    while d < -math.pi:
        d += 2.0 * math.pi
    return d


@njit(cache=False)
def _nn(i, n, alive, y, ph, R2):
    """Geometryczny najbliższy sąsiad cząstki i wśród żywych (ΔR² < R²; inaczej −1 i R² = odległość do wiązki)."""
    bd = R2; bj = -1
    for j in range(n):
        if j == i or not alive[j]:
            continue
        dy = y[i] - y[j]
        if dy * dy >= bd:
            continue
        dp = _dphi(ph[i], ph[j])
        d = dy * dy + dp * dp
        if d < bd:
            bd = d; bj = j
    return bj, bd


@njit(cache=False)
def antikt(pt, eta, phi, mass, flag, R):
    """Anty-kT, schemat E, O(N²) (jak fastjet N2Plain: per cząstka geometryczny NN, d_iJ = min(kt_i, kt_NN)·ΔR²_NN, wiązka = R²·kt_i;
    lemat Cacciari–Salam: globalne minimum d_ij jest zawsze parą (i, NN(i))). Wejście: (pT, η, φ, m) i flag (1 = D0).
    Zwraca (jpt, jeta[pseudopospieszność], jphi, jn, jflag) malejąco po pT. Zgodność z antikt_naive: _selftest."""
    n = len(pt)
    px = np.empty(n); py = np.empty(n); pz = np.empty(n); en = np.empty(n)
    y = np.empty(n); ph = np.empty(n); kt = np.empty(n); cnt = np.ones(n, np.int64); fl = np.zeros(n, np.int64)
    alive = np.ones(n, np.bool_)
    nn = np.full(n, -1, np.int64); nd = np.empty(n)
    R2 = R * R
    for i in range(n):
        px[i] = pt[i] * math.cos(phi[i]); py[i] = pt[i] * math.sin(phi[i]); pz[i] = pt[i] * math.sinh(eta[i])
        en[i] = math.sqrt(px[i] ** 2 + py[i] ** 2 + pz[i] ** 2 + mass[i] ** 2)
        y[i] = 0.5 * math.log((en[i] + pz[i]) / (en[i] - pz[i]))
        ph[i] = phi[i]; kt[i] = 1.0 / (pt[i] * pt[i]); fl[i] = flag[i]
    for i in range(n):
        nn[i], nd[i] = _nn(i, n, alive, y, ph, R2)
    jpt = np.empty(n); jeta = np.empty(n); jphi = np.empty(n); jn = np.empty(n, np.int64); jf = np.empty(n, np.int64)
    nj = 0
    left = n
    while left > 0:
        best = 1e300; bi = -1
        for i in range(n):
            if not alive[i]:
                continue
            k = kt[i]
            if nn[i] >= 0 and kt[nn[i]] < k:
                k = kt[nn[i]]
            d = nd[i] * k
            if d < best:
                best = d; bi = i
        bj = nn[bi]
        if bj < 0:
            p2 = px[bi] ** 2 + py[bi] ** 2
            ptj = math.sqrt(p2)
            jpt[nj] = ptj
            jeta[nj] = math.asinh(pz[bi] / ptj) if ptj > 0 else 0.0
            jphi[nj] = math.atan2(py[bi], px[bi]); jn[nj] = cnt[bi]; jf[nj] = fl[bi]
            nj += 1
            alive[bi] = False
            for k in range(n):
                if alive[k] and nn[k] == bi:
                    nn[k], nd[k] = _nn(k, n, alive, y, ph, R2)
        else:
            px[bi] += px[bj]; py[bi] += py[bj]; pz[bi] += pz[bj]; en[bi] += en[bj]
            cnt[bi] += cnt[bj]; fl[bi] = max(fl[bi], fl[bj])
            p2 = px[bi] ** 2 + py[bi] ** 2
            kt[bi] = 1.0 / p2
            ph[bi] = math.atan2(py[bi], px[bi])
            y[bi] = 0.5 * math.log((en[bi] + pz[bi]) / (en[bi] - pz[bi]))
            alive[bj] = False
            nn[bi], nd[bi] = _nn(bi, n, alive, y, ph, R2)
            for k in range(n):
                if not alive[k] or k == bi:
                    continue
                if nn[k] == bi or nn[k] == bj:
                    nn[k], nd[k] = _nn(k, n, alive, y, ph, R2)
                else:
                    dy = y[k] - y[bi]
                    dp = _dphi(ph[k], ph[bi])
                    d = dy * dy + dp * dp
                    if d < nd[k]:
                        nd[k] = d; nn[k] = bi
        left -= 1
    o = np.argsort(-jpt[:nj])
    return jpt[:nj][o], jeta[:nj][o], jphi[:nj][o], jn[:nj][o], jf[:nj][o]


@njit(cache=False)
def antikt_naive(pt, eta, phi, mass, flag, R):
    """Anty-kT, schemat E. Wejście: tablice cząstek (pT, η, φ, m) i flag (1 = D0). Zwraca (jpt, jeta[pseudorapidity], jphi, jn, jflag)
    posortowane malejąco po pT. Naiwne O(N³) — tylko odniesienie do testu antikt()."""
    n = len(pt)
    px = np.empty(n); py = np.empty(n); pz = np.empty(n); en = np.empty(n)
    y = np.empty(n); ph = np.empty(n); kt = np.empty(n); cnt = np.ones(n, np.int64); fl = np.zeros(n, np.int64)
    alive = np.ones(n, np.bool_)
    for i in range(n):
        px[i] = pt[i] * math.cos(phi[i]); py[i] = pt[i] * math.sin(phi[i]); pz[i] = pt[i] * math.sinh(eta[i])
        en[i] = math.sqrt(px[i] ** 2 + py[i] ** 2 + pz[i] ** 2 + mass[i] ** 2)
        y[i] = 0.5 * math.log((en[i] + pz[i]) / (en[i] - pz[i]))
        ph[i] = phi[i]; kt[i] = 1.0 / (pt[i] * pt[i]); fl[i] = flag[i]
    jpt = np.empty(n); jeta = np.empty(n); jphi = np.empty(n); jn = np.empty(n, np.int64); jf = np.empty(n, np.int64)
    nj = 0
    R2 = R * R
    left = n
    while left > 0:
        best = 1e300; bi = -1; bj = -1
        for i in range(n):
            if not alive[i]:
                continue
            if kt[i] < best:
                best = kt[i]; bi = i; bj = -1
            for j in range(i + 1, n):
                if not alive[j]:
                    continue
                dy = y[i] - y[j]
                dp = _dphi(ph[i], ph[j])
                d = min(kt[i], kt[j]) * (dy * dy + dp * dp) / R2
                if d < best:
                    best = d; bi = i; bj = j
        if bj < 0:
            p2 = px[bi] ** 2 + py[bi] ** 2
            ptj = math.sqrt(p2)
            jpt[nj] = ptj
            jeta[nj] = math.asinh(pz[bi] / ptj) if ptj > 0 else 0.0
            jphi[nj] = math.atan2(py[bi], px[bi]); jn[nj] = cnt[bi]; jf[nj] = fl[bi]
            nj += 1
            alive[bi] = False
        else:
            px[bi] += px[bj]; py[bi] += py[bj]; pz[bi] += pz[bj]; en[bi] += en[bj]
            cnt[bi] += cnt[bj]; fl[bi] = max(fl[bi], fl[bj])
            p2 = px[bi] ** 2 + py[bi] ** 2
            kt[bi] = 1.0 / p2
            ph[bi] = math.atan2(py[bi], px[bi])
            y[bi] = 0.5 * math.log((en[bi] + pz[bi]) / (en[bi] - pz[bi]))
            alive[bj] = False
        left -= 1
    o = np.argsort(-jpt[:nj])
    return jpt[:nj][o], jeta[:nj][o], jphi[:nj][o], jn[:nj][o], jf[:nj][o]


@njit(cache=False)
def _dr(e1, p1, e2, p2):
    de = e1 - e2
    dp = _dphi(p1, p2)
    return math.sqrt(de * de + dp * dp)


@njit(cache=False)
def _match(tpt, teta, tphi, used, dpt, deta, dphi, max_dr, max_rel):
    """Indeks toru najbliższego córce (ΔR < max_dr, |pT/pT_c − 1| < max_rel, jeszcze nieużytego); −1 brak. Zwraca (indeks, ΔR)."""
    bi = -1; bd = max_dr
    if not (dpt > 0):
        return -1, 0.0
    for i in range(len(tpt)):
        if used[i]:
            continue
        if abs(tpt[i] / dpt - 1.0) > max_rel:
            continue
        d = _dr(teta[i], tphi[i], deta, dphi)
        if d < bd:
            bd = d; bi = i
    return bi, (bd if bi >= 0 else 0.0)


@njit(cache=False)
def one_pair(tpt, teta, tphi, D, out, R, pt_min, eta_max, max_dr, max_rel, ix):
    """tpt [GeV], teta, tphi — WSZYSTKIE tory zdarzenia; D = (pt, eta, phi, m, d0pt, d0eta, d0phi, d1pt, d1eta, d1phi, ppt, peta, pphi).
    Wypełnia out[NOUT] (kolejność NAMES; ix = tablica indeksów NAMES w kolejności użycia poniżej = arange(NOUT))."""
    n = len(tpt)
    Dpt, Deta, Dphi, Dm = D[0], D[1], D[2], D[3]
    used = np.zeros(n, np.bool_)       # córki D0
    isp = np.zeros(n, np.bool_)        # wolny pion
    nm = 0; mdr = 0.0
    for k in range(2):
        i, d = _match(tpt, teta, tphi, used, D[4 + 3 * k], D[5 + 3 * k], D[6 + 3 * k], max_dr, max_rel)
        if i >= 0:
            used[i] = True; nm += 1; mdr = max(mdr, d)
    i, d = _match(tpt, teta, tphi, used, D[10], D[11], D[12], max_dr, max_rel)
    pm = 0
    if i >= 0:
        isp[i] = True; pm = 1; mdr = max(mdr, d)

    sel = np.zeros(n, np.bool_)
    ntrk = 0
    for i in range(n):
        if tpt[i] > pt_min and abs(teta[i]) < eta_max:
            sel[i] = True; ntrk += 1

    # ---- stożki, półkule, tor wiodący
    c_n = np.zeros(3); c_s = np.zeros(3)
    v = np.zeros(6)
    hs_a = 0.0; hn_a = 0.0; hs_t = 0.0; hn_t = 0.0
    t1 = 0.0; t1dr = NONE_DR
    for i in range(n):
        if used[i]:
            continue
        d = _dr(teta[i], tphi[i], Deta, Dphi)
        p = tpt[i]
        if d < 0.4:
            v[0] += p
            if not isp[i]:
                v[1] += p
            if p > 0.5:
                v[2] += p
                if not isp[i]:
                    v[3] += p
        if d < 1.0:
            v[4] += p
            if p > 0.5:
                v[5] += p
        if isp[i] or not sel[i]:
            continue
        if d < 0.2:
            c_n[0] += 1; c_s[0] += p
        if d < 0.4:
            c_n[1] += 1; c_s[1] += p
        if d < 1.0:
            c_n[2] += 1; c_s[2] += p
        if abs(_dphi(tphi[i], Dphi)) > 0.5 * math.pi:
            hs_a += p; hn_a += 1
        else:
            hs_t += p; hn_t += 1
        if p > t1:
            t1 = p; t1dr = d

    # ---- kolekcja A: D0 jako pseudocząstka
    na = 1
    for i in range(n):
        if sel[i] and not used[i]:
            na += 1
    apt = np.empty(na); aeta = np.empty(na); aphi = np.empty(na); am = np.empty(na); afl = np.zeros(na, np.int64)
    apt[0] = Dpt; aeta[0] = Deta; aphi[0] = Dphi; am[0] = Dm; afl[0] = 1
    k = 1
    for i in range(n):
        if sel[i] and not used[i]:
            apt[k] = tpt[i]; aeta[k] = teta[i]; aphi[k] = tphi[i]; am[k] = M_PI; k += 1
    jpt, jeta, jphi, jn, jf = antikt(apt, aeta, aphi, am, afl, R)
    n5 = 0; drn = NONE_DR; ptn = 0.0; iD = -1
    for j in range(len(jpt)):
        if jf[j] == 1:
            iD = j
        if jpt[j] > 5.0:
            n5 += 1
            d = _dr(jeta[j], jphi[j], Deta, Dphi)
            if d < drn:
                drn = d; ptn = jpt[j]
    o = 0
    out[o] = n5; o += 1
    out[o] = jpt[0]; out[o + 1] = jeta[0]; out[o + 2] = jphi[0]; o += 3
    if len(jpt) > 1:
        out[o] = jpt[1]; out[o + 1] = jeta[1]; out[o + 2] = jphi[1]
    else:
        out[o] = 0.0; out[o + 1] = 0.0; out[o + 2] = 0.0
    o += 3
    out[o] = _dr(jeta[0], jphi[0], Deta, Dphi); o += 1
    out[o] = jpt[iD]; out[o + 1] = jeta[iD]; out[o + 2] = jphi[iD]; out[o + 3] = jn[iD]; o += 4
    out[o] = _dr(jeta[iD], jphi[iD], Deta, Dphi); o += 1
    out[o] = iD; o += 1
    out[o] = Dpt / jpt[iD]; o += 1
    out[o] = drn; out[o + 1] = ptn; o += 2

    # ---- kolekcja B: bez córek i wolnego pionu
    nb = 0
    for i in range(n):
        if sel[i] and not used[i] and not isp[i]:
            nb += 1
    bn5 = 0
    b1 = np.zeros(3); b2 = np.zeros(3); bdr1 = NONE_DR
    thr = np.array([2.0, 3.0, 5.0, 8.0])
    bdrn = np.full(4, NONE_DR); bptn = np.zeros(4)
    bmax04 = 0.0; bmax08 = 0.0
    if nb > 0:
        bpt = np.empty(nb); beta = np.empty(nb); bphi = np.empty(nb); bm = np.full(nb, M_PI); bfl = np.zeros(nb, np.int64)
        k = 0
        for i in range(n):
            if sel[i] and not used[i] and not isp[i]:
                bpt[k] = tpt[i]; beta[k] = teta[i]; bphi[k] = tphi[i]; k += 1
        jpt, jeta, jphi, jn, jf = antikt(bpt, beta, bphi, bm, bfl, R)
        b1[0] = jpt[0]; b1[1] = jeta[0]; b1[2] = jphi[0]
        bdr1 = _dr(jeta[0], jphi[0], Deta, Dphi)
        if len(jpt) > 1:
            b2[0] = jpt[1]; b2[1] = jeta[1]; b2[2] = jphi[1]
        for j in range(len(jpt)):
            d = _dr(jeta[j], jphi[j], Deta, Dphi)
            if jpt[j] > 5.0:
                bn5 += 1
            for q in range(4):
                if jpt[j] > thr[q] and d < bdrn[q]:
                    bdrn[q] = d; bptn[q] = jpt[j]
            if d < 0.4 and jpt[j] > bmax04:
                bmax04 = jpt[j]
            if d < 0.8 and jpt[j] > bmax08:
                bmax08 = jpt[j]
    out[o] = bn5; o += 1
    out[o] = b1[0]; out[o + 1] = b1[1]; out[o + 2] = b1[2]; o += 3
    out[o] = b2[0]; out[o + 1] = b2[1]; out[o + 2] = b2[2]; o += 3
    out[o] = bdr1; o += 1
    for q in range(4):
        out[o] = bdrn[q]; out[o + 1] = bptn[q]; o += 2
    out[o] = bmax04; out[o + 1] = bmax08; o += 2
    for q in range(3):
        out[o] = c_n[q]; out[o + 1] = c_s[q]; o += 2
    for q in range(6):
        out[o] = v[q]; o += 1
    out[o] = hs_a; out[o + 1] = hn_a; out[o + 2] = hs_t; out[o + 3] = hn_t; o += 4
    out[o] = t1; out[o + 1] = t1dr; out[o + 2] = ntrk; o += 3
    ndlo = 0
    for k in range(2):
        if D[4 + 3 * k] < 0.2:
            ndlo += 1
    out[o] = nm; out[o + 1] = pm; out[o + 2] = mdr; out[o + 3] = ndlo; o += 4

    # ---- kolekcja O: jak B, ale tylko tory poza stożkiem 0,4 wokół D0
    no = 0
    for i in range(n):
        if sel[i] and not used[i] and not isp[i] and _dr(teta[i], tphi[i], Deta, Dphi) > 0.4:
            no += 1
    on5 = 0; o1 = 0.0; odr1 = NONE_DR; odr3 = NONE_DR; opt3 = 0.0; odr5 = NONE_DR; opt5 = 0.0; omax = 0.0
    if no > 0:
        opt = np.empty(no); oeta = np.empty(no); ophi = np.empty(no); om = np.full(no, M_PI); ofl = np.zeros(no, np.int64)
        k = 0
        for i in range(n):
            if sel[i] and not used[i] and not isp[i] and _dr(teta[i], tphi[i], Deta, Dphi) > 0.4:
                opt[k] = tpt[i]; oeta[k] = teta[i]; ophi[k] = tphi[i]; k += 1
        jpt, jeta, jphi, jn, jf = antikt(opt, oeta, ophi, om, ofl, R)
        o1 = jpt[0]; odr1 = _dr(jeta[0], jphi[0], Deta, Dphi)
        for j in range(len(jpt)):
            d = _dr(jeta[j], jphi[j], Deta, Dphi)
            if jpt[j] > 5.0:
                on5 += 1
                if d < odr5:
                    odr5 = d; opt5 = jpt[j]
            if jpt[j] > 3.0 and d < odr3:
                odr3 = d; opt3 = jpt[j]
            if d < 0.8 and jpt[j] > omax:
                omax = jpt[j]
    out[o] = on5; out[o + 1] = o1; out[o + 2] = odr1; out[o + 3] = odr3; out[o + 4] = opt3
    out[o + 5] = odr5; out[o + 6] = opt5; out[o + 7] = omax; o += 8
    return o


@njit(cache=False)
def process(tpt, teta, tphi, off, pair_ev, DD, out, R, pt_min, eta_max, max_dr, max_rel):
    """tpt/teta/tphi — spłaszczone tory wybranych zdarzeń (off[e]:off[e+1]); pair_ev[i] = indeks zdarzenia pary i; DD[i] = 13 liczb pary."""
    ix = np.arange(out.shape[1])
    for i in range(len(pair_ev)):
        e = pair_ev[i]
        a, b = off[e], off[e + 1]
        one_pair(tpt[a:b], teta[a:b], tphi[a:b], DD[i], out[i], R, pt_min, eta_max, max_dr, max_rel, ix)


def run(tpt, teta, tphi, off, pair_ev, DD):
    out = np.zeros((len(pair_ev), NOUT), np.float64)
    if len(pair_ev):
        process(np.ascontiguousarray(tpt, np.float64), np.ascontiguousarray(teta, np.float64), np.ascontiguousarray(tphi, np.float64),
                np.ascontiguousarray(off, np.int64), np.ascontiguousarray(pair_ev, np.int64), np.ascontiguousarray(DD, np.float64),
                out, R_JET, PT_MIN, ETA_MAX, MATCH_DR, MATCH_REL)
    return out.astype(np.float32)


def _selftest():
    """Test bez danych: (1) dwa wąskie pęki → 2 dżety o znanym pT; (2) D0 z córkami w pęku → z, stożki bez córek."""
    rng = np.random.default_rng(1)
    # pęk 1 wokół (0,0): D0 (5 GeV) + 3 tory po 2 GeV w ΔR < 0,15; pęk 2 wokół (0.5, 3.0): 4 tory po 3 GeV; + miękkie tło daleko
    d0 = (3.0, 0.02, 0.03); d1 = (2.0, -0.03, -0.02)
    tr = [d0, d1, (2.0, 0.10, 0.05), (2.0, -0.08, 0.10), (2.0, 0.05, -0.12), (0.3, 0.2, 0.2), (0.15, 0.01, 0.0)]
    tr += [(3.0, 0.5 + rng.normal(0, 0.03), 3.0 + rng.normal(0, 0.03)) for _ in range(4)]
    tr += [(0.6, 2.0, -2.0), (0.7, -2.2, 1.5), (0.4, 1.0, 1.0)]
    t = np.array(tr)
    px = d0[0] * np.cos(d0[2]) + d1[0] * np.cos(d1[2]); py = d0[0] * np.sin(d0[2]) + d1[0] * np.sin(d1[2])
    Dpt = float(np.hypot(px, py)); Dphi = float(np.arctan2(py, px)); Deta = 0.0
    DD = np.array([[Dpt, Deta, Dphi, 1.865, *d0, *d1, 0.15, 0.01, 0.0]])
    o = run(t[:, 0], t[:, 1], t[:, 2], np.array([0, len(t)]), np.array([0]), DD)[0]
    g = lambda k: float(o[IX[k]])  # noqa: E731
    assert g("nmatch") == 2 and g("pis_match") == 1, (g("nmatch"), g("pis_match"))
    assert abs(g("c_s04") - 6.0) < 1e-5 and g("c_n04") == 3, (g("c_s04"), g("c_n04"))        # 3 tory po 2 GeV; 0,3 GeV i pion odcięte
    assert abs(g("v_s40_p2") - (6.0 + 0.3 + 0.15)) < 1e-5 and abs(g("v_s40_p2x") - 6.3) < 1e-5
    assert abs(g("a_jD_pt") - (Dpt + 6.0)) < 0.15 and g("a_jD_n") == 4, (g("a_jD_pt"), g("a_jD_n"))
    assert abs(g("z") - Dpt / g("a_jD_pt")) < 1e-6
    assert abs(g("a_j1_pt") - 12.0) < 0.1 and g("a_jD_rank") == 1 and g("a_n5") == 2, (g("a_j1_pt"), g("a_jD_rank"))
    assert abs(g("b_j1_pt") - 12.0) < 0.1 and abs(g("b_ptmax04") - 6.0) < 0.1 and g("b_drn5") < 0.1 and g("b_drn8") > 2.5
    assert abs(g("h_s_away") - 12.6) < 1e-4 and g("h_n_away") == 5, (g("h_s_away"), g("h_n_away"))
    assert g("t1_pt") == 3.0 and g("ntrk05") == 11
    assert abs(g("o_j1_pt") - 12.0) < 0.1 and g("o_drn5") > 2.5 and g("o_ptmax08") == 0.0 and g("o_n5") == 1, (g("o_j1_pt"), g("o_drn5"))
    # niezmienniczość: Σ pT dżetów ≈ Σ pT cząstek (schemat E, wąskie pęki), każda cząstka w dokładnie jednym dżecie
    n = 60
    pt = rng.exponential(1.0, n) + 0.5; eta = rng.uniform(-2.5, 2.5, n); phi = rng.uniform(-np.pi, np.pi, n)
    jpt, jeta, jphi, jn, jf = antikt(pt, eta, phi, np.full(n, M_PI), np.zeros(n, np.int64), 0.4)
    assert jn.sum() == n and np.all(np.diff(jpt) <= 0)
    # O(N²) ≡ naiwny O(N³): gęste zdarzenia (dużo łączeń), różne N
    for n2, seed in ((5, 2), (40, 3), (150, 4), (300, 5)):
        r2 = np.random.default_rng(seed)
        p = r2.exponential(1.2, n2) + 0.5; e = r2.normal(0, 1.2, n2); f = r2.uniform(-np.pi, np.pi, n2)
        fl = np.zeros(n2, np.int64); fl[0] = 1
        x = antikt(p, e, f, np.full(n2, M_PI), fl, 0.4); yv = antikt_naive(p, e, f, np.full(n2, M_PI), fl, 0.4)
        assert len(x[0]) == len(yv[0]) and np.allclose(x[0], yv[0], rtol=1e-9) and np.array_equal(x[3], yv[3]) and np.array_equal(x[4], yv[4]), n2
        assert np.allclose(x[1], yv[1]) and np.allclose(x[2], yv[2])
    # żadne dwa dżety bliżej niż R (własność anty-kT dla twardych dżetów nie jest ścisła; sprawdzamy dżet wiodący)
    d = np.hypot(jeta[1:] - jeta[0], (jphi[1:] - jphi[0] + np.pi) % (2 * np.pi) - np.pi)
    assert d.min() > 0.4 * 0.5
    print(f"tj_kernel selftest OK (numba={HAVE_NUMBA}, NOUT={NOUT})")


if __name__ == "__main__":
    _selftest()
