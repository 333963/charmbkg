"""B_trig — emulacja triggera w MC (mb_nd / fd_nd nie mają bitów triggera: trigPassed tylko w danych v2.0).

Zmienne x = (nSelBL, lead_pt [GeV]) — obie w ROOT/h5 v2.0 także w MC. Dwie metody:

A. MAPA WYDAJNOŚCI (nominał, „fizyka z MC, trigger z danych”):
   ε_k(x) = P(łańcuch k odpala bez prescale | x) mierzone w 427929 na próbce odniesienia MB427929 (MBTS_1, PS = 1 → bez
   obciążenia, o ile k jest niezależny od MBTS_1 przy danym x); k ∈ {J12, j20MBTS, sptrk4} (+ mbts ≈ 1).
   Klasy 428770 składane z prescale’ami runu 428770 (PS_k, niezależne decyzje prescale):
       P_MBTS(x)   = 1 − Π_k∈MBTS (1 − ε_k(x)/PS_k)          (j20MBTS, sptrk4, mbts, …)
       P_J12only(x) = ε_J12(x) · Π_k∈MBTS (1 − ε_k(x)/PS_k)
   Waga MC w = P_klasa(x) (normalizacja nieistotna dla ε). Dla 427929: P_MB(x) = ε_MBTS1(x) ≈ 1 (waga płaska, tylko GRL).
B. GĘSTOŚĆ (wariant, jak w2): w(x) = p_dane,klasa(x) / p_MC(x) na zdarzeniach z ≥ 1 kandydatem D0 — wchłania też
   różnice fizyki dane/MC w krotności (nie tylko trigger), więc ε_MC,klasa(B) − ε_MC,klasa(A) = miara modelowania krotności.

Ograniczenia (PROCEDURA §4): L1_J12 to kalorymetr (neutralne, E_T jetu) — x torowe opisuje go tylko częściowo; sam kandydat D
(pT 8–30 GeV) może odpalić J12 → mapa per bin pT(D) (opcja pt_bins); ε_k z 427929 (inny run: μ, warunki L1) przeniesione do
428770; korelacje decyzji przy danym x pominięte; MC bez pile-up (dane 5,9 % nPileupVertices > 0 → weto npu == 0 w obu).
"""
from __future__ import annotations

import json

import numpy as np

NSEL_EDGES = np.array([0, 5, 10, 15, 20, 30, 40, 60, 1e9])
LEAD_EDGES = np.array([0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 1e9])
CHAINS_EMUL = {"J12": "HLT_noalg_L1J12", "j20MBTS": "HLT_j20_L1MBTS_2", "sptrk4": "HLT_mb_sptrk_pt4_L1MBTS_2",
               "mbts": "HLT_mb_mbts_L1MBTS_1"}
PS_428770 = {"J12": 1.0, "j20MBTS": 50.0, "sptrk4": 50.0, "mbts": 50006.6}   # LB 21–323 (LB 4–20: 25/25) — z trigPS gdy dostępne


def _bin2(nsel, lead, ne=NSEL_EDGES, le=LEAD_EDGES):
    i = np.clip(np.searchsorted(ne, np.asarray(nsel, float), side="right") - 1, 0, len(ne) - 2)
    j = np.clip(np.searchsorted(le, np.asarray(lead, float), side="right") - 1, 0, len(le) - 2)
    return i, j


def effmap(nsel, lead, ref_mask, fired: dict, prior=0.5, min_ref=20):
    """ε_k(x) na próbce odniesienia: (n_fired + prior·ε̄_k) / (n_ref + prior) per komórka; komórki z n_ref < min_ref → ε̄_k
    (globalne) i flaga. Zwraca dict k → {eps (2D), n_ref (2D), eps_global, low_stat (2D bool)}."""
    i, j = _bin2(nsel, lead)
    r = np.asarray(ref_mask, bool)
    shape = (len(NSEL_EDGES) - 1, len(LEAD_EDGES) - 1)
    nref = np.zeros(shape)
    np.add.at(nref, (i[r], j[r]), 1.0)
    out = {}
    for k, f in fired.items():
        f = np.asarray(f, bool) & r
        nf = np.zeros(shape)
        np.add.at(nf, (i[f], j[f]), 1.0)
        g = nf.sum() / max(nref.sum(), 1.0)
        eps = (nf + prior * g) / (nref + prior)
        low = nref < min_ref
        eps[low] = g
        out[k] = {"eps": eps, "n_ref": nref, "eps_global": float(g), "low_stat": low}
    return out


def class_prob(em: dict, nsel, lead, cls: str, ps: dict | None = None):
    """P(klasa | x) z map ε_k i prescale’ów (metoda A). cls ∈ {J12only, MBTSbased, j20MBTS, sptrk4, all, MB427929, …}."""
    ps = {**PS_428770, **(ps or {})}
    i, j = _bin2(nsel, lead)
    e = {k: em[k]["eps"][i, j] if k in em else np.zeros(len(i)) for k in CHAINS_EMUL}
    none_mbts = np.ones(len(i))
    for k in ("j20MBTS", "sptrk4", "mbts"):
        none_mbts *= 1.0 - np.clip(e[k] / ps[k], 0, 1)
    if cls == "J12only":
        return e["J12"] * none_mbts
    if cls == "MBTSbased":
        return 1.0 - none_mbts
    if cls in ("j20MBTS", "sptrk4"):
        return e[cls] / ps[cls]
    if cls == "other":
        return np.clip((1 - e["J12"]) * none_mbts, 0, 1)      # przybliżenie (AFP/RD0 nieemulowane)
    if cls == "sptrk4_427929":
        return e["sptrk4"]
    return np.ones(len(i))                                    # all, all427929, MB427929 (MBTS_1 ≈ 1 przy D)


def density_weights(nsel_mc, lead_mc, nsel_d, lead_d, min_mc=20):
    """Metoda B: w(x) = p_dane(x)/p_MC(x) (histogramy 2D, normalizowane); komórki z < min_mc MC → 1."""
    im, jm = _bin2(nsel_mc, lead_mc)
    idd, jd = _bin2(nsel_d, lead_d)
    shape = (len(NSEL_EDGES) - 1, len(LEAD_EDGES) - 1)
    hm, hd = np.zeros(shape), np.zeros(shape)
    np.add.at(hm, (im, jm), 1.0)
    np.add.at(hd, (idd, jd), 1.0)
    pm, pd_ = hm / max(hm.sum(), 1), hd / max(hd.sum(), 1)
    r = np.where(hm >= min_mc, pd_ / np.where(pm > 0, pm, 1), 1.0)
    return r[im, jm], {"ratio": r.tolist(), "n_mc": hm.tolist(), "n_data": hd.tolist()}


def save_map(em: dict, path: str, meta: dict | None = None):
    js = {"nsel_edges": NSEL_EDGES.tolist(), "lead_edges": LEAD_EDGES.tolist(), "meta": meta or {},
          "chains": {k: {"eps": v["eps"].tolist(), "n_ref": v["n_ref"].tolist(), "eps_global": v["eps_global"]} for k, v in em.items()}}
    json.dump(js, open(path, "w"), indent=0)


def load_map(path: str) -> dict:
    js = json.load(open(path))
    return {k: {"eps": np.array(v["eps"]), "n_ref": np.array(v["n_ref"]), "eps_global": v["eps_global"]} for k, v in js["chains"].items()}


def n_eff(w):
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / max((w ** 2).sum(), 1e-300))
