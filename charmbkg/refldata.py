"""Szablon refleksji K↔π D0 z DANYCH do fitu masy (etap S5, addytywnie; metoda S: bkggen.control.reflection).

Kandydat-refleksja (zła hipoteza K/π prawdziwego D0) ma masę po zamianie hipotez m_swap = masa poprawnego partnera ≈ m0.
Szablon(m) = hist(m | |m_swap−m0| < win) − (win/szer. SB)·hist(m | m_swap w SB) — odjęcie kombinatoryki w m_swap,
po tym samym cięciu WP co fit (score kandydata-refleksji). Ujemne biny → 0, wygładzenie Gaussem (smooth_bins).
Normalizacja (`fit.refl_data.norm`):
  mc_ratio — N_refl = r_MC · N_sig (jak nominalnie; z danych tylko kształt)
  data     — N_refl = ΣSzablon / f_win, f_win = ułamek sygnału w |m−m0| < win (kształt MC); r = N_refl/N_sig iterowane (3×)
Konfiguracja: fit.reflection_source: mc (domyślnie, bez zmian) | data;  fit.refl_data: {...} (DEFAULTS niżej).
"""
from __future__ import annotations

import logging

import numpy as np

from . import config as C

log = logging.getLogger("charmbkg")
DEFAULTS = {"sample": "data_none", "win": 0.025, "sb": [0.050, 0.100], "smooth_bins": 1.5, "min_peak": 200, "norm": "mc_ratio",
            "n_iter": 3, "half_both": False, "centre": "mc"}


def opts(fc) -> dict:
    o = dict(DEFAULTS)
    o.update(fc.get("refl_data") or {})
    return o


def swap_mass_rows(cfg, key, species, rows: np.ndarray) -> np.ndarray:
    """m_swap dla wierszy (globalnych, jak `row` w scores) próbki `key`; czyta tylko kinematykę córek K, π."""
    from .samples import _rr_slice, colmap
    sp = cfg["species"][species]
    k = int(sp.get("kaon_index", 0))
    p = int((sp.get("pion_indices") or [1])[0])
    cm = colmap(cfg, C.sample_path(cfg, key), int(sp["n_daughters"]))
    cols = [cm.get(f"{q}{i}") for i in (k, p) for q in ("d_pt", "d_eta", "d_phi")]
    if any(c is None for c in cols):
        raise KeyError("refldata: brak kolumn kinematyki córek")
    upt = float((cfg.get("units") or {}).get("d_pt", 1.0))
    src = C.source(cfg, key)
    rows = np.asarray(rows, np.int64)
    order = np.argsort(rows, kind="stable")
    rs = rows[order]
    out = np.full(len(rows), np.nan)
    chunk = int(cfg["io"]["chunk_rows"])
    from .bkgsrc import M_K, M_PI, _fv
    i = 0
    while i < len(rs):
        start = int(rs[i])
        stop = start + chunk
        j = int(np.searchsorted(rs, stop, "left"))
        raw = _rr_slice(src, cols, cfg["io"]["dataset"], start, min(stop, int(rs[j - 1]) + 1))
        loc = rs[i:j] - start
        a = raw.iloc[loc]
        g = lambda c: a[c].to_numpy(dtype=float)  # noqa: E731
        ck, cp = cols[:3], cols[3:]
        s = _fv(g(ck[0]) * upt, g(ck[1]), g(ck[2]), M_PI) + _fv(g(cp[0]) * upt, g(cp[1]), g(cp[2]), M_K)
        out[order[i:j]] = np.sqrt(np.maximum(s[:, 0] ** 2 - (s[:, 1:] ** 2).sum(1), 0))
        i = j
    return out


def template(m, ms, m0, lo, hi, nbins, o):
    """(histogram szablonu ≥0 wygładzony, N_peak surowe, Σ szablonu przed wygładzeniem)."""
    d = np.abs(np.asarray(ms) - m0)
    win, (s0, s1) = float(o["win"]), o["sb"]
    pk = d < win
    sb = (d >= s0) & (d < s1)
    w = np.where(pk, 1.0, -win / (s1 - s0))
    if o.get("half_both", True):  # obie hipotezy pary w oknie (m i m_swap): para daje 2 kandydatów, z których 1 to sygnał → waga ½
        w = np.where(pk & (np.abs(np.asarray(m) - m0) < win), 0.5, w)
    sel = pk | sb
    h, e = np.histogram(np.asarray(m)[sel], nbins, (lo, hi), weights=w[sel])
    h = np.clip(h, 0, None)
    tot = float(h.sum())
    sb_ = float(o.get("smooth_bins") or 0)
    if sb_ > 0 and tot > 0:
        x = np.arange(-int(4 * sb_) - 1, int(4 * sb_) + 2)
        k = np.exp(-0.5 * (x / sb_) ** 2)
        h = np.convolve(h, k / k.sum(), "same")
        h *= tot / max(h.sum(), 1e-12)
    return h, int(pk.sum()), tot


def hybrid(h_data, h_mc, edges, m0, win):
    """Centrum |m−m0| < win z MC (tam metoda z danych miesza się z sygnałem: poprawne D0 z m_swap w oknie/SB),
    ogony z danych; centrum MC skalowane stosunkiem ogonów dane/MC. Zwraca (szablon, ogony_dane, ułamek centrum MC)."""
    c = np.abs(0.5 * (np.asarray(edges[1:]) + np.asarray(edges[:-1])) - m0) < win
    hd = np.asarray(h_data, float).copy()
    hm = np.asarray(h_mc, float)
    od, om = hd[~c].sum(), hm[~c].sum()
    hd[c] = hm[c] * (od / om) if om > 0 else 0.0
    return hd, float(od), float(hm[c].sum() / max(hm.sum(), 1e-12))


def f_window(shape, m0, win):
    from .massfit import _gcdf
    mu, s1 = shape["mu"], shape["s1"]
    if shape.get("model") == "gauss2":
        f, r = shape["f"], shape["r"]
        F = lambda x: f * _gcdf(x, mu, s1) + (1 - f) * _gcdf(x, mu, s1 * r)  # noqa: E731
    else:
        F = lambda x: _gcdf(x, mu, s1)  # noqa: E731
    return float(F(m0 + win) - F(m0 - win))


def fit_with_template(m, lo, hi, m0, shape, fc, tmpl, rr_mc, tsum, o):
    """Fit danych z kształtem refleksji z danych; norm=data → r iterowane z ΣSzablonu/f_win."""
    from .massfit import fit_data
    rr = float(rr_mc)
    info = {"refl_source": "data", "refl_norm": o["norm"], "refl_ratio_mc": float(rr_mc), "template_sum": tsum}
    if o["norm"] == "data" and tsum > 0:
        fw = f_window(shape, m0, float(o["win"]))
        n_refl = tsum / max(fw, 1e-6)
        info.update({"f_win": fw, "N_refl_data": n_refl})
        r = None
        for _ in range(int(o["n_iter"])):
            r = fit_data(m, lo, hi, m0, shape, fc, tmpl, rr)
            rr = float(np.clip(n_refl / max(r["N_sig"], 1.0), 0.0, 5.0))
        r = fit_data(m, lo, hi, m0, shape, fc, tmpl, rr)
    else:
        r = fit_data(m, lo, hi, m0, shape, fc, tmpl, rr)
    info["refl_ratio_used"] = rr
    r.update(info)
    return r
