"""S6: próbki kontrolne na poziomie kandydata z kolumn h5 v1.1 (bez Athena, bez refitu wierzchołka).

Transformacje jednej córki kandydata (ramka logiczna charmbkg, GeV):
  rot_pi     — φ córki += π            (klasyczna „rotacja toru”: niszczy korelację z pikiem)
  rot_rand   — φ córki += U(π/2, 3π/2) (losowy kąt z dala od 0)
  flip_eta   — η córki → −η
  swap       — zamiana hipotez mas K↔π (D0: refleksja; D⁺: K↔π_j = „pseudo-WS”, bo ładunków córek brak w h5)
  mix        — córka zastąpiona torem (córką kandydata) z INNEGO zdarzenia o podobnym pT (event mixing toru)
Po transformacji przeliczane: m, pT, η, φ D, cosθ*, m(Kπ)_1,2 (D⁺), m_swap (D0), Dalitz (D⁺) i cechy pochodne
charmbkg (features.compute). Zmienne wierzchołka (Lxy, χ², IP, IPsig, z0 D, cosθ_P) NIE są przeliczane —
`vertex_mode`: keep (z kandydata źródłowego) | resample (losowanie z pasm bocznych warunkowo w binach pT×|η| D).
Zmienne torowe (d0, z0, dE/dx, npts) wędrują z torem (przy mix — z toru z innego zdarzenia).

API: make(L, sp, species, method, ...) → ramka kontrolna (kolumny charmbkg + ctrl, src_row, w_ctrl).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .. import features as F
from . import kin3
from .data import M_K, M_PI

MASSES = {"D0": (M_K, M_PI), "Dplus": (M_K, M_PI, M_PI)}
PER_D = list(F.PER_D)
VERTEX_COLS = ["Lxy", "chi2", "impact", "impact_sig", "impact_z0", "cos_theta_p"]
METHODS = ("rot_pi", "rot_rand", "flip_eta", "swap", "mix")
ETA_MAX = {"D0": None, "Dplus": 2.2}  # |η_D| z preselekcji v1.1 (D⁺ ucięte na 2,2; D0 bez cięcia)
PT_RANGE = (3.0, 30.0)


# ----------------------------------------------------------------------------- kinematyka
def _p4(L, nd, masses):
    g = lambda q, i: L[f"{q}{i}"].to_numpy(dtype=float)  # noqa: E731
    out = []
    for i in range(nd):
        pt, eta, phi = g("d_pt", i), g("d_eta", i), g("d_phi", i)
        px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        out.append(np.stack([np.sqrt(px ** 2 + py ** 2 + pz ** 2 + masses[i] ** 2), px, py, pz], -1))
    return out


def _m(*p):
    s = sum(p)
    return np.sqrt(np.maximum(s[:, 0] ** 2 - (s[:, 1:] ** 2).sum(-1), 0))


def cos_theta_star(p4D, p4K):
    """cosθ* = kąt K w spoczynku D względem kierunku lotu D (definicja jak Dmeson_cosThetaStar, sprawdzone)."""
    P, E = p4D[:, 1:], p4D[:, 0]
    k, Ek = p4K[:, 1:], p4K[:, 0]
    kst, _ = kin3._boost(k, Ek, P / E[:, None], -1)
    return (kst * P).sum(-1) / (np.linalg.norm(kst, axis=-1) * np.linalg.norm(P, axis=-1))


def recompute(L: pd.DataFrame, sp: dict, species: str) -> pd.DataFrame:
    """Przelicza kinematykę D z córek (hipotezy MASSES, K = kaon_index) i cechy pochodne. Modyfikuje L."""
    nd = int(sp["n_daughters"])
    ms = MASSES[species]
    p = _p4(L, nd, ms)
    P = sum(p)
    pt = np.hypot(P[:, 1], P[:, 2])
    L["m"] = _m(*p)
    L["pt"] = pt
    L["eta"] = np.arcsinh(P[:, 3] / np.maximum(pt, 1e-12))
    L["phi"] = np.arctan2(P[:, 2], P[:, 1])
    L["cos_theta_star"] = cos_theta_star(P, p[0])
    for i in range(nd):
        L[f"d_p{i}"] = np.nan
    if species == "D0":
        a = p[0].copy(); b = p[1].copy()  # noqa: E702
        a[:, 0] = np.sqrt((a[:, 1:] ** 2).sum(-1) + M_PI ** 2)
        b[:, 0] = np.sqrt((b[:, 1:] ** 2).sum(-1) + M_K ** 2)
        L["m_swap"] = _m(a, b)
    else:
        L["m_kpi1"], L["m_kpi2"] = _m(p[0], p[1]), _m(p[0], p[2])
        L["s12"], L["s13"] = L["m_kpi1"] ** 2, L["m_kpi2"] ** 2
    F.compute(L, sp, units={})
    return L


# ----------------------------------------------------------------------------- transformacje
def _which(n, nd, daughter, rng):
    if daughter == "random":
        return rng.integers(0, nd, n)
    return np.full(n, int(daughter))


def rotate_phi(L, idx, angle):
    for i in np.unique(idx):
        s = idx == i
        c = f"d_phi{i}"
        v = L[c].to_numpy(dtype=float).copy()
        v[s] = (v[s] + angle[s] + np.pi) % (2 * np.pi) - np.pi
        L[c] = v


def flip_eta(L, idx):
    for i in np.unique(idx):
        s = idx == i
        v = L[f"d_eta{i}"].to_numpy(dtype=float).copy()
        v[s] = -v[s]
        L[f"d_eta{i}"] = v


def swap_hyp(L, i, j, rows=None):
    """Zamiana wszystkich kolumn córek i↔j (indeks K pozostaje stały → K dostaje tor j); rows = maska wierszy."""
    s = np.ones(len(L), bool) if rows is None else np.asarray(rows)
    for q in PER_D:
        a, b = L[f"{q}{i}"].to_numpy().astype(float), L[f"{q}{j}"].to_numpy().astype(float)
        na, nb = np.where(s, b, a), np.where(s, a, b)
        L[f"{q}{i}"], L[f"{q}{j}"] = na, nb


def _event_key(L):
    return L["run"].fillna(0).astype(np.int64).to_numpy() * 10 ** 10 + L["event"].astype(np.int64).to_numpy()


def track_pool(L, nd):
    """Pula torów = córki kandydatów (unikalne w zdarzeniu po pT) z wszystkimi kolumnami toru."""
    ev = _event_key(L)
    parts = []
    for i in range(nd):
        d = {q: L[f"{q}{i}"].to_numpy() for q in PER_D}
        d["ev"] = ev
        parts.append(pd.DataFrame(d))
    T = pd.concat(parts, ignore_index=True)
    T["key"] = np.round(T["d_pt"].to_numpy(dtype=float), 5)
    T = T.drop_duplicates(["ev", "key"]).reset_index(drop=True)
    T["lpt"] = np.log(T["d_pt"].to_numpy(dtype=float))
    return T.sort_values("lpt").reset_index(drop=True)


def mix_track(L, idx, pool, rng, window=50):
    """Córka idx ← tor z innego zdarzenia o zbliżonym pT (losowy z ±window sąsiadów w rankingu log pT)."""
    ev = _event_key(L)
    lp = pool["lpt"].to_numpy()
    pev = pool["ev"].to_numpy()
    n = len(L)
    nd = sum(1 for c in L.columns if c.startswith("d_pt") and c[4:].isdigit())
    pts = np.column_stack([L[f"d_pt{i}"].to_numpy(dtype=float) for i in range(nd)])
    cur = pts[np.arange(n), idx]
    pos = np.searchsorted(lp, np.log(cur))
    pick = np.clip(pos + rng.integers(-window, window + 1, n), 0, len(pool) - 1)
    for _ in range(5):  # ponów, gdy tor z tego samego zdarzenia
        bad = pev[pick] == ev
        if not bad.any():
            break
        pick[bad] = np.clip(pos[bad] + rng.integers(-window, window + 1, bad.sum()), 0, len(pool) - 1)
    ok = pev[pick] != ev
    for i in np.unique(idx):
        s = idx == i
        for q in PER_D:
            v = L[f"{q}{i}"].to_numpy().astype(float).copy()
            v[s] = pool[q].to_numpy(dtype=float)[pick[s]]
            L[f"{q}{i}"] = v
    return ok


# ----------------------------------------------------------------------------- preselekcja / wierzchołek
def presel_mask(L, sp, species, pt_range=PT_RANGE):
    lo, hi = sp["mass_range"]
    m, pt = L["m"].to_numpy(dtype=float), L["pt"].to_numpy(dtype=float)
    ok = (m >= lo) & (m < hi) & (pt >= pt_range[0]) & (pt < pt_range[1])
    em = ETA_MAX.get(species)
    if em is not None:
        ok &= np.abs(L["eta"].to_numpy(dtype=float)) < em
    return ok


def vertex_resample(C, ref, rng, pt_bins=(3, 4, 5, 6, 8, 10, 14, 30), eta_bins=(0, 0.5, 1.0, 1.5, 2.0, 3.0)):
    """Hot-deck: zmienne wierzchołka C ← losowy kandydat ref (np. pasma boczne) z tego samego binu (pT, |η|) D."""
    def cell(X):
        ip = np.clip(np.digitize(X["pt"].to_numpy(dtype=float), pt_bins) - 1, 0, len(pt_bins) - 2)
        ie = np.clip(np.digitize(np.abs(X["eta"].to_numpy(dtype=float)), eta_bins) - 1, 0, len(eta_bins) - 2)
        return ip * 10 + ie
    cc, cr = cell(C), cell(ref)
    R = ref[VERTEX_COLS].to_numpy(dtype=float)
    out = C[VERTEX_COLS].to_numpy(dtype=float).copy()
    for k in np.unique(cc):
        s = np.flatnonzero(cc == k)
        pool = np.flatnonzero(cr == k)
        if not len(pool):
            pool = np.arange(len(ref))
        out[s] = R[pool[rng.integers(0, len(pool), len(s))]]
    for j, c in enumerate(VERTEX_COLS):
        C[c] = out[:, j]
    return C


# ----------------------------------------------------------------------------- główne API
def make(L: pd.DataFrame, sp: dict, species: str, method: str, n_rep: int = 1, daughter="random", seed: int = 0,
         source_mask=None, vertex_mode: str = "keep", vertex_ref: pd.DataFrame | None = None, presel: bool = True,
         pool: pd.DataFrame | None = None, swap_pion: int | None = None) -> pd.DataFrame:
    """Ramka kontrolna: n_rep kopii każdego kandydata źródłowego (source_mask) po transformacji `method`.
    w_ctrl = 1/n_rep (suma wag = liczba źródeł, które trafiły w okno preselekcji)."""
    if method not in METHODS:
        raise ValueError(method)
    rng = np.random.default_rng(seed)
    nd = int(sp["n_daughters"])
    src = np.flatnonzero(np.ones(len(L), bool) if source_mask is None else np.asarray(source_mask))
    if method == "mix" and pool is None:
        pool = track_pool(L, nd)
    reps = 1 if method in ("flip_eta", "swap") and daughter != "random" else n_rep
    outs = []
    base = L.iloc[src].reset_index(drop=True)
    for r in range(reps):
        C = base.copy()
        n = len(C)
        ok = np.ones(n, bool)
        if method == "swap":
            j = (rng.integers(1, nd, n) if swap_pion is None else np.full(n, swap_pion)) if nd > 2 else np.ones(n, int)
            for jj in np.unique(j):
                swap_hyp(C, int(sp.get("kaon_index", 0)), int(jj), rows=(j == jj))
        else:
            idx = _which(n, nd, daughter, rng)
            if method == "rot_pi":
                rotate_phi(C, idx, np.full(n, np.pi))
            elif method == "rot_rand":
                rotate_phi(C, idx, rng.uniform(0.5 * np.pi, 1.5 * np.pi, n))
            elif method == "flip_eta":
                flip_eta(C, idx)
            elif method == "mix":
                ok = mix_track(C, idx, pool, rng)
        C["src_m"] = C["m"].to_numpy()
        C["src_truth"] = C["truth"].to_numpy()
        recompute(C, sp, species)
        if presel:
            ok &= presel_mask(C, sp, species)
        C = C[ok].copy()
        C["src_row"] = src[ok]
        C["rep"] = r
        outs.append(C)
    C = pd.concat(outs, ignore_index=True)
    C["w_ctrl"] = 1.0 / reps
    C["ctrl"] = method
    if vertex_mode == "resample":
        if vertex_ref is None:
            raise ValueError("vertex_mode=resample wymaga vertex_ref")
        vertex_resample(C, vertex_ref, rng)
        F.compute(C, sp, units={})
    elif vertex_mode != "keep":
        raise ValueError(vertex_mode)
    return C


# ----------------------------------------------------------------------------- CLI (podpróbki npz → npz kontrolne)
def main(argv=None):
    """python3 -m charmbkg.bkggen.candrot --npz 'sub/dp_data_p*.npz' --species Dplus --method swap --source sb -o out.npz"""
    import argparse
    import glob

    from .. import config as Cf
    from .data import load_npz
    ap = argparse.ArgumentParser(description="S6: próbki kontrolne na poziomie kandydata (bez Athena)")
    ap.add_argument("--npz", required=True)
    ap.add_argument("-c", "--config", default=None)
    ap.add_argument("--species", default="D0", choices=["D0", "Dplus"])
    ap.add_argument("--method", default="swap", choices=METHODS)
    ap.add_argument("--source", default="sb", choices=["all", "sb"])
    ap.add_argument("--n-rep", type=int, default=3)
    ap.add_argument("--vertex", default="keep", choices=["keep", "resample"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args(argv)
    cfg = Cf.load(a.config)
    sp = cfg["species"][a.species]
    L = pd.concat([load_npz(cfg, p, a.species) for p in sorted(glob.glob(a.npz))], ignore_index=True)
    recompute(L, sp, a.species)
    sb = np.abs(L["m"].to_numpy(float) - sp["m0"]) > sp["sideband_min"]
    src = sb if a.source == "sb" else None
    ref = L[sb].reset_index(drop=True) if a.vertex == "resample" else None
    C = make(L, sp, a.species, a.method, n_rep=a.n_rep, seed=a.seed, source_mask=src, vertex_mode=a.vertex,
             vertex_ref=ref)
    num = {c: C[c].to_numpy() for c in C.columns if C[c].dtype.kind in "fiub"}
    np.savez_compressed(a.out, **num)
    print(f"{a.method}/{a.source}/{a.vertex}: {len(L)} kandydatów → {len(C)} obiektów kontrolnych (Σw={C['w_ctrl'].sum():.0f}) → {a.out}")


if __name__ == "__main__":
    main()
