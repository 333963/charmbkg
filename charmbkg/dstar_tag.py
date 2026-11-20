"""HYP T4 — tag D*⁺→D0π⁺_soft dla wierszy D0 (obie hipotezy K/π): poprawna hipoteza na DANYCH z ładunku miękkiego π.

Format h5 v1.2 (konwerter `Q_lxplus/h2_pdg/charmpp_h5_v12.py`, README_H2):
  * D* = decayType 1, córki 0 = K (hipoteza), 1 = π, 2 = π_soft; `Dmeson_delta_mass` = m(Kππ) − m(Kπ) [GeV];
    `Dmeson_charge` = znak pdgId (dla D* = ładunek D*, czyli π_soft); `Dmeson_fit_charge` dla D* = 0 (fit D0) — nie używać.
  * D0 = decayType 10, córka 0 = K (hipoteza), 1 = π; obie hipotezy tej samej pary torów to dwa wiersze.
  * `Dmeson_daughter{i}_trackId` / `_charge` (int32, −999 = brak córki).
  * Kolejność wierszy w zdarzeniu: D⁺, D*, D0, Λc (D* przed D0 → przy cięciu porcji wystarczy przenieść D* ostatniego zdarzenia).

Fizyka: D*⁺ → D0 π⁺_s, D0 → K⁻π⁺ (D*⁻ sprzężone). Ładunek K = −q(π_s). Wiersz D0 jest zgodny z zapachem (poprawna hipoteza),
gdy jego córka 0 (hipoteza K) ma ładunek −q(π_s); inaczej to refleksja K↔π. DCS D0→K⁺π⁻ (R ≈ 0,38 %) i mieszanie dają ułamek
~0,4 % „fałszywych refleksji” — poziom tła tej metody.

Dopasowanie D0↔D*: to samo zdarzenie (blok `evblk`: ciągłe wiersze o tym samym (plik, run, event)) i ten sam NIEUPORZĄDKOWANY zbiór
trackId córek {t0, t1} D0 = {t0, t1} D*. Kilka D* na jeden D0 (różne π_soft) → wybór D* z Δm najbliżej 145,43 MeV, licznik
`dstar_n` i flaga `dstar_conflict` (różne ładunki π_soft wśród dopasowanych).

Kolumny wyjściowe (wiersz D0):
  dstar_tag       int8    1 = jest dopasowany D* z 0 < Δm < dm_max, 0 = brak
  dm              float32 Δm wybranego D* [GeV] (NaN bez tagu)
  soft_pi_charge  int8    ±1 (0 bez tagu)
  hypo_from_tag   int8    0 = wiersz zgodny z zapachem, 1 = refleksja, −1 = brak tagu / brak ładunków
  dstar_n         int8    liczba dopasowanych D* (≤ 127)
  dstar_conflict  int8    1 = dopasowane D* mają różne ładunki π_soft
  dstar_row       int64   wiersz globalny wybranego D* (−1)
  pis_pt          float32 pT π_soft [GeV] (NaN)
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .constants import DM_DSTAR

DSTAR_DT = 1
DM_MAX = 0.165          # GeV, jak EV_dstar (ev_extract.DM_MAX)
TAG_COLS = ["dstar_tag", "dm", "soft_pi_charge", "hypo_from_tag", "dstar_n", "dstar_conflict", "dstar_row", "pis_pt"]
# kolumny h5 v1.2 potrzebne tylko dla wierszy D* (córki 0..2); nazwy logiczne → kandydaci w pliku
DSTAR_V12 = {
    "ds_trk": ["Dmeson_daughter{i}_trackId"],
    "ds_q": ["Dmeson_daughter{i}_charge"],
}


def _valid_q(q):
    q = np.asarray(q, dtype=np.float64)
    return np.isfinite(q) & (np.abs(q) == 1)


def dstar_table(evblk, row, trk0, trk1, trk2, dm, q_soft, q0=None, q1=None, pis_pt=None, dm_max: float = DM_MAX) -> pd.DataFrame:
    """Tabela kandydatów D* (jeden wiersz = jeden D*). q_soft: ładunek π_soft (córka 2; brak → Dmeson_charge).
    q0/q1: ładunki córek 0/1 D* (do ustalenia toru K, gdy wiersz D0 nie ma ładunków). Odrzuca: brak trackId, Δm poza (0, dm_max),
    ładunek ≠ ±1."""
    n = len(np.asarray(evblk))
    a = lambda v, fill=np.nan: np.full(n, fill, np.float64) if v is None else np.asarray(v, np.float64)  # noqa: E731
    t0, t1, t2 = a(trk0, -1), a(trk1, -1), a(trk2, -1)
    dm, qs = a(dm), a(q_soft)
    q0, q1 = a(q0), a(q1)
    q0[~_valid_q(q0)] = np.nan
    q1[~_valid_q(q1)] = np.nan
    ok = (t0 >= 0) & (t1 >= 0) & np.isfinite(dm) & (dm > 0) & (dm < dm_max) & _valid_q(qs)
    T = pd.DataFrame({
        "evblk": np.asarray(evblk, np.int64)[ok],
        "lo": np.minimum(t0, t1)[ok].astype(np.int64), "hi": np.maximum(t0, t1)[ok].astype(np.int64),
        "ds_t0": t0[ok].astype(np.int64), "ds_t1": t1[ok].astype(np.int64), "ds_t2": t2[ok].astype(np.int64),
        "ds_dm": dm[ok], "ds_qs": qs[ok].astype(np.int8),
        "ds_q0": q0[ok], "ds_q1": q1[ok],
        "ds_row": np.asarray(row, np.int64)[ok], "ds_pis_pt": a(pis_pt)[ok],
    })
    return T


def empty_tags(n: int) -> dict:
    return {"dstar_tag": np.zeros(n, np.int8), "dm": np.full(n, np.nan, np.float32), "soft_pi_charge": np.zeros(n, np.int8),
            "hypo_from_tag": np.full(n, -1, np.int8), "dstar_n": np.zeros(n, np.int8), "dstar_conflict": np.zeros(n, np.int8),
            "dstar_row": np.full(n, -1, np.int64), "pis_pt": np.full(n, np.nan, np.float32)}


def tag_d0(evblk, trk0, trk1, T: pd.DataFrame, q0=None) -> dict:
    """Tagi dla wierszy D0 (evblk, trackId córek 0/1, ładunek córki 0 = hipoteza K). T = `dstar_table`. Zwraca dict kolumn TAG_COLS."""
    evblk = np.asarray(evblk, np.int64)
    n = len(evblk)
    out = empty_tags(n)
    if n == 0 or T is None or not len(T):
        return out
    t0, t1 = np.asarray(trk0, np.float64), np.asarray(trk1, np.float64)
    okt = (t0 >= 0) & (t1 >= 0)
    if not okt.any():
        return out
    idx = np.flatnonzero(okt)
    D = pd.DataFrame({"_i": idx, "evblk": evblk[idx], "lo": np.minimum(t0, t1)[idx].astype(np.int64),
                      "hi": np.maximum(t0, t1)[idx].astype(np.int64), "_t0": t0[idx].astype(np.int64)})
    J = D.merge(T, on=["evblk", "lo", "hi"], how="inner")
    if not len(J):
        return out
    J["_dd"] = np.abs(J["ds_dm"].to_numpy() - DM_DSTAR)
    g = J.groupby("_i")
    nmatch = g.size()
    nq = g["ds_qs"].nunique()
    B = J.sort_values(["_i", "_dd", "ds_row"], kind="mergesort").drop_duplicates("_i", keep="first")
    i = B["_i"].to_numpy()
    qs = B["ds_qs"].to_numpy().astype(np.int8)
    # hipoteza: ładunek córki 0 wiersza D0 (gdy znany), inaczej tor K z ładunków córek D*, inaczej córka 0 D* (konwencja charmpp)
    hyp = np.full(len(B), -1, np.int8)
    if q0 is not None:
        qd = np.asarray(q0, np.float64)[i]
        v = _valid_q(qd)
        hyp[v] = np.where(qd[v] == -qs[v], 0, 1)
    else:
        v = np.zeros(len(B), bool)
    rest = ~v
    if rest.any():
        dq0, dq1 = B["ds_q0"].to_numpy(), B["ds_q1"].to_numpy()
        kid = np.where(_valid_q(dq0) & (dq0 == -qs), B["ds_t0"].to_numpy(),
                       np.where(_valid_q(dq1) & (dq1 == -qs), B["ds_t1"].to_numpy(), B["ds_t0"].to_numpy()))
        hyp[rest] = np.where(B["_t0"].to_numpy()[rest] == kid[rest], 0, 1)
    out["dstar_tag"][i] = 1
    out["dm"][i] = B["ds_dm"].to_numpy(np.float32)
    out["soft_pi_charge"][i] = qs
    out["hypo_from_tag"][i] = hyp
    out["dstar_n"][i] = np.minimum(nmatch.loc[i].to_numpy(), 127).astype(np.int8)
    out["dstar_conflict"][i] = (nq.loc[i].to_numpy() > 1).astype(np.int8)
    out["dstar_row"][i] = B["ds_row"].to_numpy(np.int64)
    out["pis_pt"][i] = B["ds_pis_pt"].to_numpy(np.float32)
    return out


def carry_last_event(T: pd.DataFrame, last_evblk) -> pd.DataFrame:
    """D* z ostatniego zdarzenia porcji — przenoszone do następnej (zdarzenie przecięte granicą porcji)."""
    if T is None or not len(T) or last_evblk is None:
        return T.iloc[0:0] if T is not None else None
    return T[T["evblk"].to_numpy() == int(last_evblk)].reset_index(drop=True)


def confusion(hypo_from_tag, hypo_swap, mask=None) -> dict:
    """MC: hypo_from_tag vs prawda hypo_swap (tylko wiersze z oboma ≥ 0). Klucze 'tag{a}_true{b}'."""
    h, s = np.asarray(hypo_from_tag), np.asarray(hypo_swap)
    m = (h >= 0) & (s >= 0)
    if mask is not None:
        m &= np.asarray(mask, bool)
    return {f"tag{a}_true{b}": int(((h == a) & (s == b) & m).sum()) for a in (0, 1) for b in (0, 1)}


def merge_counts(acc: dict, add: dict) -> dict:
    for k, v in add.items():
        acc[k] = acc.get(k, 0) + int(v)
    return acc
