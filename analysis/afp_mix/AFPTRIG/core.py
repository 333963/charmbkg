#!/usr/bin/env python3
"""AFPTRIG/core.py — szybka wersja estymatora kształtu na tablicach (wspólna dla danych i toy).

prep_pairs(): dla każdego zdarzenia centralnego ułamek protonów jego bloku LB dający Δ w zadanych przedziałach.
shape_estimate(): N_acc = Σ_c N_CR,c · R_c po komórkach c (blok LB × klasa aktywności), z błędami.
"""
import numpy as np

DELTA, M_NOM = 0.02, 0.025


def prep_pairs(xc, gc, xp, gp, wp, los, his):
    """a[i, k] = Σ_j∈blok(g_i) wp_j·1[lo_k < xc_i − xp_j ≤ hi_k] / Σ_j wp_j ;  oraz per blok: lista (xp, wp) do błędu protonów."""
    a = np.zeros((len(xc), len(los)))
    for g in np.unique(gc):
        mc = np.flatnonzero(gc == g); mp = gp == g
        if not mp.any():
            a[mc] = np.nan; continue
        o = np.argsort(xp[mp]); xs = xp[mp][o]; cw = np.concatenate([[0.0], np.cumsum(wp[mp][o])])
        for k, (lo, hi) in enumerate(zip(los, his)):
            i1 = np.searchsorted(xs, xc[mc] - hi, side="left"); i2 = np.searchsorted(xs, xc[mc] - lo, side="left")
            a[mc, k] = (cw[i2] - cw[i1]) / cw[-1]
    return a


def proton_relvar(xc, wc, gc, xp, gp, wp, sr, cr):
    """Względna wariancja R = T(SR)/T(CR) od skończonej liczby protonów, per blok LB (estymator ilorazowy po protonach)."""
    out = {}
    for g in np.unique(gc):
        mc = gc == g; mp = gp == g
        if not mp.any() or not mc.any():
            out[int(g)] = np.nan; continue
        o = np.argsort(xc[mc]); xs = xc[mc][o]; cw = np.concatenate([[0.0], np.cumsum(wc[mc][o])])
        x = xp[mp]; w = wp[mp]
        b = []
        for lo, hi in (sr, cr):
            j1 = np.searchsorted(xs, x + lo, side="right"); j2 = np.searchsorted(xs, x + hi, side="right")
            b.append(cw[j2] - cw[j1])
        S0, S1 = (w * b[0]).sum(), (w * b[1]).sum()
        if S0 <= 0 or S1 <= 0:
            out[int(g)] = np.nan; continue
        R = S0 / S1
        out[int(g)] = float(((w * (b[0] - R * b[1])) ** 2).sum() / S1 ** 2 / R ** 2)
    return out


def shape_estimate(d_cell, d_w, d_in_cr, c_cell, c_w, a_sr, a_cr, cell2g=None, prv=None, ncell=None):
    """d_*: zdarzenia danych (komórka, waga, czy w CR); c_*: pula centralna (komórka, waga, ułamki par SR/CR).
    → dict(n_acc, err_cr (Poisson N_CR), err_tc (stat puli centralnej), err_tp (stat protonów), n_cr, n_cr_nopool)."""
    nc = int(max(d_cell.max(initial=0), c_cell.max(initial=0)) + 1) if ncell is None else ncell
    ok = np.isfinite(a_sr) & np.isfinite(a_cr)
    N = np.bincount(d_cell[d_in_cr], d_w[d_in_cr], nc); N2 = np.bincount(d_cell[d_in_cr], d_w[d_in_cr] ** 2, nc)
    S0 = np.bincount(c_cell[ok], (c_w * a_sr)[ok], nc); S1 = np.bincount(c_cell[ok], (c_w * a_cr)[ok], nc)
    good = S1 > 0
    R = np.where(good, S0 / np.where(good, S1, 1), 0.0)
    v = np.bincount(c_cell[ok], ((c_w * (a_sr - R[c_cell] * a_cr)) ** 2)[ok], nc) / np.where(good, S1, 1) ** 2
    n_acc = (N * R)[good].sum()
    out = {"n_acc": float(n_acc), "err_cr": float(np.sqrt((N2 * R * R)[good].sum())), "err_tc": float(np.sqrt((N * N * v)[good].sum())),
           "n_cr": float(N.sum()), "n_cr_nopool": float(N[~good].sum()), "R_eff": float(n_acc / max(N[good].sum(), 1e-300))}
    if prv is not None and cell2g is not None:
        g = cell2g[:nc]; tp = 0.0
        for gg, rv in prv.items():
            s = (N * R)[good & (g == gg)].sum()
            if np.isfinite(rv):
                tp += s * s * rv
        out["err_tp"] = float(np.sqrt(tp))
    else:
        out["err_tp"] = 0.0
    out["err"] = float(np.sqrt(out["err_cr"] ** 2 + out["err_tc"] ** 2 + out["err_tp"] ** 2))
    return out
