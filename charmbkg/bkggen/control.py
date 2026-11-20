"""Klasyczne próbki kontrolne tła dostępne w h5 v1.1 (bez nowych ntupli).

  sideband    — bliskie / dalekie pasma masy (D0, D+)
  reflection  — D0: zamiana hipotez K↔π; kandydaci, których masa „lustrzana” trafia w pik → szablon refleksji z danych
  dstar_wrong — D*: (a) pasmo Δm poza pikiem, (b) D0 + losowy tor z innego kandydata tego samego zdarzenia (π_s losowy)
  mix         — pary D–D z różnych zdarzeń (pule: liczba kandydatów w zdarzeniu × bin pT), porównanie z parami z tego samego zdarzenia

Każda funkcja zwraca ramkę (kopię wierszy wejścia + kolumny pomocnicze) z kolumną `ctrl` = nazwa próbki.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .data import M_PI, four_vec, inv_mass, swapped_mass

# ----------------------------------------------------------------------------- regiony masy (D0 domyślnie)
REGIONS = {
    "D0": {"sr": (0.0, 0.055), "inner": (0.040, 0.075), "near": (0.075, 0.100), "far": (0.100, 0.200)},
    "Dplus": {"sr": (0.0, 0.045), "inner": (0.045, 0.070), "near": (0.070, 0.110), "far": (0.110, 0.400)},
}


def region_mask(L, sp, species, name):
    lo, hi = REGIONS[species][name]
    d = np.abs(L["m"].to_numpy(dtype=float) - sp["m0"])
    return (d >= lo) & (d < hi)


def sideband(L, sp, species="D0", which="near"):
    """which: near | far | both."""
    if which == "both":
        m = region_mask(L, sp, species, "near") | region_mask(L, sp, species, "far")
    else:
        m = region_mask(L, sp, species, which)
    out = L[m].copy()
    out["ctrl"] = f"sideband_{which}"
    out["side"] = np.where(out["m"] < sp["m0"], "L", "R")
    return out


# ----------------------------------------------------------------------------- refleksja K↔π (D0)
def reflection(L, sp, win=0.025, sb=(0.050, 0.100)):
    """Kandydaci z m_swap w oknie piku (|m_swap−m0|<win) = wzbogacone w prawdziwe D0 o błędnej hipotezie.
    Waga −(win)/(szerokość SB) dla kandydatów z m_swap w pasmach bocznych sb → odejmowanie kombinatoryki
    (szablon refleksji = Σ wag w binach m). Zwraca wiersze z kolumnami m_swap i w_ctrl."""
    ms = swapped_mass(L, sp.get("kaon_index", 0), (sp.get("pion_indices") or [1])[0])
    d = np.abs(ms - sp["m0"])
    pk = d < win
    s = (d >= sb[0]) & (d < sb[1])
    out = L[pk | s].copy()
    out["m_swap"] = ms[pk | s]
    out["w_ctrl"] = np.where(pk[pk | s], 1.0, -win / (sb[1] - sb[0]))
    out["ctrl"] = "reflection"
    return out


def reflection_template(R, lo, hi, nbins=100):
    """Histogram szablonu refleksji (z odjęciem SB w m_swap), ujemne biny → 0."""
    h, e = np.histogram(R["m"], nbins, (lo, hi), weights=R["w_ctrl"])
    return np.clip(h, 0, None), e


def reflection_template_var(R, lo, hi, nbins=100, clip=False):
    """AA (S7 §2.2): szablon refleksji BEZ obcinania ujemnych binów + wariancja Σw² per bin (błędy Poissona
    ważone: pik − SB). Obcinanie (clip=True) tworzy szablon sztucznie przy pustej statystyce — używać tylko do rysowania."""
    h, e = np.histogram(R["m"], nbins, (lo, hi), weights=R["w_ctrl"])
    v, _ = np.histogram(R["m"], nbins, (lo, hi), weights=np.asarray(R["w_ctrl"], float) ** 2)
    return (np.clip(h, 0, None) if clip else h), v, e


# ----------------------------------------------------------------------------- D*: zły Δm / losowy π_s
def dstar_wrong_dm(Ls, dm_sb=(0.150, 0.165)):
    """Ls: kandydaci D* (decayType 1). Pasmo Δm poza pikiem."""
    dm = Ls["dm"].to_numpy(dtype=float)
    out = Ls[(dm >= dm_sb[0]) & (dm < dm_sb[1])].copy()
    out["ctrl"] = "dstar_wrong_dm"
    return out


def dstar_random_pis(L, sp, n_per=3, seed=1, pt_min=0.0):
    """D0 (decayType 10) + tor z INNEGO kandydata tego samego zdarzenia (nie córka tego D0) w hipotezie π.
    Daje Δm = m(Kππ_x) − m(Kπ) dla losowego „miękkiego pionu” — kształt tła kombinatorycznego Δm.
    Pula torów = córki wszystkich kandydatów D0 w zdarzeniu (unikalne po pT). Brak ładunków córek w v1.1 → bez warunku znaku."""
    rng = np.random.default_rng(seed)
    ev = L["run"].astype(np.int64).to_numpy() * 10**10 + L["event"].astype(np.int64).to_numpy()
    tr = pd.DataFrame({"ev": np.r_[ev, ev],
                       "pt": np.r_[L["d_pt0"], L["d_pt1"]], "eta": np.r_[L["d_eta0"], L["d_eta1"]],
                       "phi": np.r_[L["d_phi0"], L["d_phi1"]]})
    tr["key"] = np.round(tr["pt"].to_numpy(), 5)
    tr = tr.drop_duplicates(["ev", "key"]).sort_values("ev").reset_index(drop=True)
    tr = tr[tr["pt"] > pt_min]
    evs = tr["ev"].to_numpy()
    starts = np.searchsorted(evs, ev, "left")
    stops = np.searchsorted(evs, ev, "right")
    n = stops - starts
    rows, pick = [], []
    k0 = np.round(L["d_pt0"].to_numpy(), 5)
    k1 = np.round(L["d_pt1"].to_numpy(), 5)
    keys = tr["key"].to_numpy()
    for j in range(n_per):
        r = starts + (rng.random(len(L)) * np.maximum(n, 1)).astype(np.int64)
        ok = (n > 2) & (r < stops)
        r = np.where(ok, r, 0)
        ok &= (keys[r] != k0) & (keys[r] != k1)
        rows.append(np.flatnonzero(ok))
        pick.append(r[ok])
    idx, tix = np.concatenate(rows), np.concatenate(pick)
    S = L.iloc[idx].reset_index(drop=True).copy()
    T = tr.iloc[tix].reset_index(drop=True)
    g = lambda q, i: S[f"{q}{i}"].to_numpy(dtype=float)  # noqa: E731
    from .data import M_K
    a = four_vec(g("d_pt", 0), g("d_eta", 0), g("d_phi", 0), M_K)
    b = four_vec(g("d_pt", 1), g("d_eta", 1), g("d_phi", 1), M_PI)
    c = four_vec(T["pt"].to_numpy(), T["eta"].to_numpy(), T["phi"].to_numpy(), M_PI)
    S["dm_rand"] = inv_mass(a, b, c) - inv_mass(a, b)
    S["pis_pt"] = T["pt"].to_numpy()
    S["ctrl"] = "dstar_random_pis"
    return S


# ----------------------------------------------------------------------------- mixing par D–D
def _event_key(L):
    return L["run"].astype(np.int64).to_numpy() * 10**10 + L["event"].astype(np.int64).to_numpy()


def same_event_pairs(L, max_per_event=20):
    """Pary kandydatów (i<j) z tego samego zdarzenia, bez wspólnych torów (po pT córek)."""
    ev = _event_key(L)
    o = np.argsort(ev, kind="stable")
    evs = ev[o]
    b = np.flatnonzero(np.r_[True, evs[1:] != evs[:-1], True])
    I, J = [], []
    for s, e in zip(b[:-1], b[1:]):
        k = min(e - s, max_per_event)
        if k < 2:
            continue
        ii, jj = np.triu_indices(k, 1)
        I.append(o[s + ii])
        J.append(o[s + jj])
    if not I:
        return np.array([], int), np.array([], int)
    I, J = np.concatenate(I), np.concatenate(J)
    p = lambda c: np.round(L[c].to_numpy(), 5)  # noqa: E731
    a0, a1, = p("d_pt0"), p("d_pt1")
    share = (a0[I] == a0[J]) | (a0[I] == a1[J]) | (a1[I] == a0[J]) | (a1[I] == a1[J])
    return I[~share], J[~share]


def pool_index(L, mult_bins=(1, 3, 6, 12, 25, 10**6), pt_bins=(3, 5, 8, 30)):
    """Pula = (krotność kandydatów w zdarzeniu, bin pT). Brak z_vtx i krotności torów w h5 v1.1 → krotność kandydatów jako proxy."""
    ev = _event_key(L)
    _, inv, cnt = np.unique(ev, return_inverse=True, return_counts=True)
    mult = cnt[inv]
    im = np.clip(np.digitize(mult, mult_bins) - 1, 0, len(mult_bins) - 2)
    ip = np.clip(np.digitize(L["pt"].to_numpy(), pt_bins) - 1, 0, len(pt_bins) - 2)
    return im * 10 + ip, mult


def mixed_pairs(L, I, J, seed=2, factor=1):
    """Dla każdej pary z tego samego zdarzenia (I,J): drugi kandydat zastąpiony losowym kandydatem z INNEGO zdarzenia
    z tej samej puli (krotność × pT drugiej nogi). Zwraca (I, Jmix)."""
    rng = np.random.default_rng(seed)
    pool, _ = pool_index(L)
    ev = _event_key(L)
    Is, Js = [], []
    by = {k: np.flatnonzero(pool == k) for k in np.unique(pool)}
    for _ in range(factor):
        Jm = np.empty_like(J)
        for k, members in by.items():
            sel = np.flatnonzero(pool[J] == k)
            if not len(sel):
                continue
            Jm[sel] = members[rng.integers(0, len(members), len(sel))]
        bad = ev[Jm] == ev[I]
        Jm[bad] = J[bad]  # rzadkie trafienie w to samo zdarzenie → odrzuć niżej
        keep = ~bad
        Is.append(I[keep])
        Js.append(Jm[keep])
    return np.concatenate(Is), np.concatenate(Js)


def pair_frame(L, I, J, tag):
    d = lambda c: L[c].to_numpy()  # noqa: E731
    dphi = np.abs((d("phi")[I] - d("phi")[J] + np.pi) % (2 * np.pi) - np.pi)
    return pd.DataFrame({"m1": d("m")[I], "m2": d("m")[J], "pt1": d("pt")[I], "pt2": d("pt")[J],
                         "dphi": dphi, "deta": np.abs(d("eta")[I] - d("eta")[J]),
                         "t1": d("truth")[I], "t2": d("truth")[J], "ctrl": tag})
