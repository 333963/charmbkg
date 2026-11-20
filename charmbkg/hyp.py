"""HYP (2026-11-08): arbiter hipotezy masowej D0→Kπ (para wierszy tych samych torów, swap K↔π).

ntupla zapisuje obie hipotezy K/π tej samej pary torów jako osobne wiersze (córka 0 = hipoteza K). Arbiter przewiduje,
który wiersz pary ma poprawną hipotezę — BEZ masy i bez zmiennych, z których masa daje się odtworzyć:

* zestaw `pid` (domyślny): dE/dx, npts i log p obu torów + pulle Bethe-Blocha (K/π) + Δχ² hipotez.
  Z samych |p| obu torów masy nie da się policzyć (brak kąta otwarcia) → wybór nie rzeźbi tła.
* zestaw `pid_kin`: + α Armenteros, cosθ*, pT córek. UWAGA: (p0, p1, α) wyznaczają p_L obu torów, q_T i kąt otwarcia,
  czyli PEŁNĄ masę obu hipotez → klasyfikator uczy się po cichu „closest_mass" i rzeźbi tło (test `sculpt_test`).
  Tylko do porównań / pokazania efektu.

Model wierszowy h(f), antysymetryzowany: P(A poprawny) = ½·[h(f_A) + 1 − h(f_B)], f_B = cechy bliźniaka (zestaw `pid`:
f_B = swap(f_A) dokładnie z jednego wiersza → liczalne w `iter_species` bez grupowania). Para: zostaje wiersz z P > ½
(P_A + P_B = 1); pojedyncze wiersze bez zmian.
"""
from __future__ import annotations

import json
import logging
import os
import pickle

import numpy as np
import pandas as pd

from .constants import M_D0, M_K, M_PI

log = logging.getLogger("charmbkg")

DEDX_SCALE = 0.85      # dE/dx danych [MIP-jedn. ntupli] / plateau BB (π przy p ≫ 1 GeV ≈ 0,87 w ntupli, BB plateau ≈ 1,03)
FEATS = {
    "pid": ["dx0", "dx1", "n0", "n1", "lp0", "lp1", "zK0", "zP0", "zK1", "zP1", "dchi"],
    "pid_kin": ["dx0", "dx1", "n0", "n1", "lp0", "lp1", "zK0", "zP0", "zK1", "zP1", "dchi", "alpha", "cts", "pt0", "pt1"],
}
_SWAP = {"dx0": "dx1", "dx1": "dx0", "n0": "n1", "n1": "n0", "lp0": "lp1", "lp1": "lp0", "zK0": "zK1", "zK1": "zK0",
         "zP0": "zP1", "zP1": "zP0", "pt0": "pt1", "pt1": "pt0"}
_NEG = {"dchi", "alpha"}
METHODS = ("keep_both", "hash", "closest_mass", "score", "arbiter")


# ----------------------------------------------------------------------------- cechy
def _a(L, c):
    return L[c].to_numpy(dtype=float) if c in L else np.full(len(L), np.nan)


def track_p(L, i):
    p = _a(L, f"d_p{i}")
    alt = _a(L, f"d_pt{i}") * np.cosh(_a(L, f"d_eta{i}"))
    return np.where(np.isfinite(p) & (p > 0), p, alt)


def armenteros_alpha(L) -> np.ndarray:
    """α = (p_L0 − p_L1)/(p_L0 + p_L1) względem kierunku p0+p1 (córka 0 = hipoteza K). Zmienia znak przy swapie."""
    pt0, pt1 = _a(L, "d_pt0"), _a(L, "d_pt1")
    e0, e1, f0, f1 = _a(L, "d_eta0"), _a(L, "d_eta1"), _a(L, "d_phi0"), _a(L, "d_phi1")
    p0 = np.stack([pt0 * np.cos(f0), pt0 * np.sin(f0), pt0 * np.sinh(e0)], 1)
    p1 = np.stack([pt1 * np.cos(f1), pt1 * np.sin(f1), pt1 * np.sinh(e1)], 1)
    d = p0 + p1
    n = np.linalg.norm(d, axis=1)
    with np.errstate(all="ignore"):
        u = d / n[:, None]
        l0, l1 = (p0 * u).sum(1), (p1 * u).sum(1)
        return (l0 - l1) / (l0 + l1)


def row_features(L: pd.DataFrame, dedx_scale: float = DEDX_SCALE) -> pd.DataFrame:
    """Cechy arbitra dla wiersza (córka 0 = K). Bez masy. Zwraca nową ramkę (te same indeksy)."""
    from .bethe_bloch import pull
    out = {}
    for i in (0, 1):
        p = track_p(L, i)
        dx = _a(L, f"d_dedx{i}")
        n = _a(L, f"d_npts{i}")
        dx = np.where((dx > 0) & ~(n == 0), dx, np.nan)
        out[f"dx{i}"], out[f"n{i}"], out[f"lp{i}"] = dx, n, np.log(np.clip(p, 1e-3, None))
        out[f"zK{i}"] = pull(dx, p, M_K, scale=dedx_scale)
        out[f"zP{i}"] = pull(dx, p, M_PI, scale=dedx_scale)
        out[f"pt{i}"] = _a(L, f"d_pt{i}")
    z = lambda x: np.where(np.isfinite(x), x, 0.0) ** 2  # noqa: E731  brak pomiaru → 0 (brak informacji)
    out["dchi"] = (z(out["zK0"]) + z(out["zP1"])) - (z(out["zP0"]) + z(out["zK1"]))   # > 0 → hipoteza wiersza gorsza
    out["alpha"] = armenteros_alpha(L)
    out["cts"] = _a(L, "cos_theta_star")
    return pd.DataFrame(out, index=L.index)


def swap_features(F: pd.DataFrame) -> pd.DataFrame:
    """Cechy bliźniaka (K↔π) z jednego wiersza: zamiana torów, α i Δχ² zmieniają znak. cosθ* NIE do odtworzenia → NaN."""
    G = pd.DataFrame(index=F.index)
    for c in F.columns:
        if c in _SWAP:
            G[c] = F[_SWAP[c]].to_numpy()
        elif c in _NEG:
            G[c] = -F[c].to_numpy()
        else:
            G[c] = np.nan
    return G


# ----------------------------------------------------------------------------- pary
def track_keys(L: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """ID torów córek 0/1: trackId v1.2 (`d_trk*` ≥ 0), inaczej bity (d0, z0·sinθ_PV) (PF), inaczej bity pT (zaokr.)."""
    out = []
    for i in (0, 1):
        t = _a(L, f"d_trk{i}")
        if np.isfinite(t).all() and (t >= 0).all():
            out.append(t.astype(np.int64))
            continue
        if f"d_d0{i}" in L and f"d_z0{i}" in L and np.isfinite(_a(L, f"d_d0{i}")).all():
            from .tracks import trk_id_d0z0 as _trk_id_d0z0
            alt = _trk_id_d0z0(np.ones(len(L), np.float32), _a(L, f"d_d0{i}"), _a(L, f"d_z0{i}"))
        else:
            alt = np.round(_a(L, f"d_pt{i}") * 1e4).astype(np.int64)
        out.append(np.where(np.isfinite(t) & (t >= 0), np.nan_to_num(t).astype(np.int64), alt))
    return out[0], out[1]


def event_blocks(run, event, start: int = 0, prev=None):
    """Id bloku zdarzenia = numer ciągłego bloku wierszy o tym samym (run, event) w kolejności pliku (fd: event_number NIE jest
    unikalny między plikami ani blokami → klucz pary musi brać blok, nie numer). prev = (run, event, blk) ostatniego wiersza
    poprzedniej porcji (ciągłość między porcjami); zwraca (blk, nowy prev)."""
    run = np.asarray(run, np.int64)
    event = np.asarray(event, np.int64)
    if not len(event):
        return np.zeros(0, np.int64), prev
    ch = np.ones(len(event), bool)
    ch[1:] = (run[1:] != run[:-1]) | (event[1:] != event[:-1])
    if prev is None:
        blk = start + np.cumsum(ch) - 1
    else:
        ch[0] = not (run[0] == prev[0] and event[0] == prev[1])
        blk = prev[2] + np.cumsum(ch)
    return blk.astype(np.int64), (int(run[-1]), int(event[-1]), int(blk[-1]))


def pair_index(L: pd.DataFrame, sample_col: str | None = "sample"):
    """(gid, n, rank): gid = id grupy (blok zdarzenia `evblk` gdy jest, inaczej (run, event); nieuporządkowana para torów),
    n = liczność, rank = 0/1 w grupie (po `row`)."""
    a, b = track_keys(L)
    ev = L["evblk"].to_numpy().astype(np.int64) if "evblk" in L else L["event"].to_numpy().astype(np.int64)
    key = pd.DataFrame({"s": L[sample_col].astype(str).to_numpy() if sample_col and sample_col in L else "",
                        "run": _a(L, "run").astype(np.int64) if "run" in L else 0, "ev": ev,
                        "lo": np.minimum(a, b), "hi": np.maximum(a, b), "row": _a(L, "row") if "row" in L else np.arange(len(L))})
    g = key.groupby(["s", "run", "ev", "lo", "hi"], sort=False)
    gid = g.ngroup().to_numpy()
    n = g["row"].transform("size").to_numpy()
    rank = g["row"].rank(method="first").to_numpy().astype(int) - 1
    return gid, n, rank


def twin_index(gid: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Indeks bliźniaka w parze (n == 2), inaczej −1."""
    tw = np.full(len(gid), -1, np.int64)
    m = n == 2
    idx = np.flatnonzero(m)
    order = idx[np.argsort(gid[idx], kind="stable")]
    a, b = order[0::2], order[1::2]
    tw[a], tw[b] = b, a
    return tw


def _mix64(x):
    x = np.asarray(x, np.uint64)
    with np.errstate(over="ignore"):
        x = x + np.uint64(0x9E3779B97F4A7C15)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return x ^ (x >> np.uint64(31))


def choose(L: pd.DataFrame, method: str, p_hyp=None, m0: float = M_D0, gid=None, n=None, rank=None) -> np.ndarray:
    """Maska wierszy pozostawionych: pojedyncze zawsze; para → jeden wiersz wg `method`:
    keep_both | hash (ślepo na masę, z klucza pary) | closest_mass (RZEŹBI tło — tylko odniesienie) | score (wyższy score D0_H)
    | arbiter (wyższe P(poprawny), `p_hyp`)."""
    if gid is None:
        gid, n, rank = pair_index(L)
    keep = np.ones(len(L), bool)
    if method == "keep_both":
        return keep
    tw = twin_index(gid, n)
    pr = tw >= 0
    i = np.flatnonzero(pr)
    j = tw[i]
    if method == "hash":
        a, b = track_keys(L)
        k = _mix64((_a(L, "run").astype(np.int64).astype(np.uint64) * np.uint64(1000003)) ^ L["event"].to_numpy().astype(np.int64).astype(np.uint64)
                   ^ _mix64(np.minimum(a, b).astype(np.uint64)) ^ (_mix64(np.maximum(a, b).astype(np.uint64)) >> np.uint64(1)))
        pick = (k % np.uint64(2)).astype(int)
        keep[i] = rank[i] == pick[i]
    elif method == "closest_mass":
        d = np.abs(_a(L, "m") - m0)
        keep[i] = (d[i] < d[j]) | ((d[i] == d[j]) & (rank[i] == 0))
    elif method in ("score", "arbiter"):
        s = _a(L, "score") if method == "score" else np.asarray(p_hyp, float)
        s = np.nan_to_num(s, nan=-1.0)
        keep[i] = (s[i] > s[j]) | ((s[i] == s[j]) & (rank[i] == 0))
    else:
        raise ValueError(f"nieznana metoda {method!r} ({METHODS})")
    if (n > 2).any():   # > 2 wiersze na parę torów (nie powinno się zdarzać) → zostaw wszystkie, raportuj
        log.warning("hyp.choose: %d grup z > 2 wierszami (pozostawione bez arbitrażu)", int((n > 2).sum()))
    return keep


# ----------------------------------------------------------------------------- model
def _backend():
    try:
        import xgboost  # noqa: F401
        return "xgboost"
    except Exception:
        return "sklearn"


def _new(backend, seed=1, n_estimators=400):
    if backend == "xgboost":
        import xgboost as xgb
        return xgb.XGBClassifier(n_estimators=n_estimators, learning_rate=0.08, max_depth=5, min_child_weight=5, subsample=0.8,
                                 colsample_bytree=0.9, tree_method="hist", objective="binary:logistic", eval_metric="logloss",
                                 early_stopping_rounds=30, n_jobs=0, random_state=seed)
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(max_iter=n_estimators, learning_rate=0.08, max_depth=5, early_stopping=True,
                                          validation_fraction=0.1, n_iter_no_change=30, random_state=seed)


def decor_weights(L: pd.DataFrame, y: np.ndarray, nb: int = 30) -> np.ndarray:
    """Wagi dekorelacyjne par treningowych (T1, 09.11): w komórkach 30 kwantyli |ln p_K/p_π| × p_min {<1,1-2,2-4,>4} GeV
    znak ln(p_K/p_π) jest 50/50 dla etykiety ⇒ reguła „K = twardszy tor” (= wybór niższej masy, ρ(Δm, ln p_K/p_π) = −0,997)
    przestaje się opłacać. Bez tego arbiter `pid` rzeźbi tło bez cięcia (Z = +21)."""
    y = np.asarray(y, int)
    lr = np.log(np.clip(track_p(L, 0), 1e-3, None)) - np.log(np.clip(track_p(L, 1), 1e-3, None))
    x = np.where(y == 1, lr, -lr)                                  # prawdziwe ln(p_K/p_π)
    pm = np.minimum(track_p(L, 0), track_p(L, 1))
    qa = np.quantile(np.abs(x), np.linspace(0, 1, nb + 1))
    cell = np.clip(np.searchsorted(qa, np.abs(x), side="right") - 1, 0, nb - 1) * 4 + np.searchsorted([1.0, 2.0, 4.0], pm)
    s = x > 0
    w = np.ones(len(L))
    for c in np.unique(cell):
        m = cell == c
        f = s[m].mean()
        if 0 < f < 1:
            w[m & s], w[m & ~s] = 0.5 / f, 0.5 / (1 - f)
    return w


# -- PID jednotorowy (mode="pertrack", T1): l = LLR_K/π(dE/dx | npts, ln p, |η|, pulle) toru, P(A) = σ(l0' − l1')
_PT_PB = np.array([0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.35, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8.5, 10, 13, 1e9])
_PT_EB = np.array([0, 0.8, 1.6, 2.6])
_PT_EBW = np.array([0, 0.6, 1.2, 1.8, 2.6])


def _pt_cell(p, n, aeta):
    ip = np.clip(np.searchsorted(_PT_PB, p, side="right") - 1, 0, len(_PT_PB) - 2)
    inn = np.clip(np.nan_to_num(n, nan=1), 1, 4).astype(int) - 1
    ie = np.clip(np.searchsorted(_PT_EB, aeta, side="right") - 1, 0, len(_PT_EB) - 2)
    return ip, inn, ie


def _pt_feats(x, n, p, aeta, dedx_scale):
    from .bethe_bloch import pull
    return np.stack([x, n, np.log(np.clip(p, 1e-3, None)), aeta, pull(x, p, M_K, scale=dedx_scale),
                     pull(x, p, M_PI, scale=dedx_scale)], 1).astype(np.float32)


def pertrack_table(L: pd.DataFrame) -> pd.DataFrame:
    """Unikalne tory (blok/run, event, trackId) obu córek z prawdą K/π (`d_tpdg`) i pomiarem dE/dx."""
    a, b = track_keys(L)
    run = _a(L, "run")
    # fd MC: event_number powtarza się między blokami/plikami → klucz zdarzenia = blok `evblk` (globalnie unikalny w próbce),
    # inaczej dedup po (run, event, trackId) zlewa tory różnych zdarzeń (audyt 16.11: −40 % torów fd)
    ev = L["evblk"].to_numpy() if "evblk" in L else L["event"].to_numpy()
    if "evblk" in L and "sample" in L:   # evblk koduje (plik, run, event) w obrębie próbki
        run = pd.factorize(L["sample"])[0].astype(np.int64)
    parts = [pd.DataFrame({"run": run, "event": ev, "key": k, "x": _a(L, f"d_dedx{i}"), "n": _a(L, f"d_npts{i}"),
                           "p": track_p(L, i), "aeta": np.abs(_a(L, f"d_eta{i}")), "pdg": np.abs(_a(L, f"d_tpdg{i}"))})
             for i, k in ((0, a), (1, b))]
    T = pd.concat(parts, ignore_index=True).drop_duplicates(["run", "event", "key"])
    return T[(T.x > 0) & (T.n > 0) & np.isfinite(T.x) & np.isfinite(T.p) & T.pdg.isin([211, 321])].reset_index(drop=True)


def _pt_flat_weights(T):
    """Płaskie wagi K/π w komórkach p × |η| ⇒ logit klasyfikatora = LLR (bez prioru gatunku w kinematyce)."""
    c = (np.clip(np.searchsorted(_PT_PB, T.p, side="right") - 1, 0, len(_PT_PB) - 2) * 10
         + np.clip(np.searchsorted(_PT_EBW, T.aeta, side="right") - 1, 0, len(_PT_EBW) - 2))
    k = (T.pdg == 321).to_numpy()
    w = np.zeros(len(T))
    for cc in np.unique(c):
        m = c == cc
        nk, npi = (m & k).sum(), (m & ~k).sum()
        if nk and npi:
            w[m & k], w[m & ~k] = 0.5 / nk, 0.5 / npi
    return w * len(T) / max(w.sum(), 1e-12)


class Arbiter:
    """k-fold (po zdarzeniu, splitmix jak charmbkg) model wierszowy + antysymetryzacja. `predict` używa modelu foldu,
    w którym zdarzenie NIE było trenowane (OOF dla MC treningowego, dowolny fold dla danych).

    mode="pairs" (domyślnie, jak dotąd): model na wierszach par sygnału; `decor=True` → wagi `decor_weights` (zalecane,
    bez nich zestaw `pid` rzeźbi tło bez cięcia). mode="pertrack" (v2, T1): model K/π per tor uczony na torach z prawdą
    (`d_tpdg`), płaskie wagi K/π w p×|η|, centrowanie na pionach (l' = l − E_π[l | p, npts, |η|]); P(A) = σ(l0' − l1').
    Ślepy na masę z konstrukcji (para ππ ma E[l0' − l1'] = 0 przy każdym p0, p1), kalibrowalny per tor na D*-tag."""

    def __init__(self, feats="pid", n_folds=5, dedx_scale=DEDX_SCALE, seed=12345, mode="pairs", center=True):
        self.feat_set = feats
        self.feats = list(FEATS[feats]) if isinstance(feats, str) else list(feats)
        self.n_folds, self.dedx_scale, self.seed = int(n_folds), float(dedx_scale), int(seed)
        self.mode, self.center = mode, bool(center)
        if mode not in ("pairs", "pertrack"):
            raise ValueError(f"Arbiter: mode {mode!r} (pairs|pertrack)")
        self.models, self.backend, self.info, self.centers = [], _backend(), [], []

    def fold(self, L):
        from .samples import fold_of
        return fold_of(L["run"].fillna(0).astype(np.int64), L["event"].astype(np.int64), self.n_folds)

    def fit(self, L: pd.DataFrame, y: np.ndarray, n_estimators=400, sample_weight=None, decor=False) -> np.ndarray:
        """L: wiersze treningowe (obie hipotezy par sygnału), y = 1 gdy hipoteza wiersza poprawna. Zwraca OOF h(f)
        (mode="pairs") albo OOF P (mode="pertrack"; tory brane z wierszy L, y nieużywane). sample_weight: wagi wierszy;
        decor=True i brak sample_weight → `decor_weights(L, y)`."""
        if self.mode == "pertrack":
            return self._fit_pertrack(L)
        F = row_features(L, self.dedx_scale)[self.feats].to_numpy(np.float32)
        y = np.asarray(y, int)
        w = np.asarray(sample_weight, float) if sample_weight is not None else (decor_weights(L, y) if decor else np.ones(len(L)))
        fo = self.fold(L)
        rng = np.random.default_rng(self.seed)
        oof = np.full(len(L), np.nan, np.float32)
        self.models = []
        for k in range(self.n_folds):
            te = fo == k
            tr = ~te
            va = tr & (rng.random(len(L)) < 0.1)
            ft = tr & ~va
            m = _new(self.backend, self.seed + k, n_estimators)
            if self.backend == "xgboost":
                m.fit(F[ft], y[ft], sample_weight=w[ft], eval_set=[(F[va], y[va])], sample_weight_eval_set=[w[va]], verbose=False)
            else:
                m.fit(F[tr], y[tr], sample_weight=w[tr])
            oof[te] = m.predict_proba(F[te])[:, 1]
            self.models.append(m)
            self.info.append({"fold": k, "n_train": int(ft.sum()), "n_test": int(te.sum()), "decor": bool(decor or sample_weight is not None)})
        return oof

    def _fit_pertrack(self, L, n_estimators=600):
        if self.backend != "xgboost":
            raise SystemExit("hyp: mode=pertrack wymaga xgboost")
        import xgboost as xgb
        T = pertrack_table(L)
        if not len(T):
            raise SystemExit("hyp: mode=pertrack — brak torów z prawdą K/π (kolumny d_tpdg0/1) i pomiarem dE/dx")
        T["fold"] = self.fold(T)
        y = (T.pdg == 321).to_numpy().astype(int)
        F = _pt_feats(T.x.to_numpy(), T.n.to_numpy(), T.p.to_numpy(), T.aeta.to_numpy(), self.dedx_scale)
        rng = np.random.default_rng(self.seed)
        self.models, self.centers = [], []
        for k in range(self.n_folds):
            tr = (T.fold != k).to_numpy()
            w = np.zeros(len(T))
            w[tr] = _pt_flat_weights(T[tr])
            va = tr & (rng.random(len(T)) < 0.1)
            ft = tr & ~va
            m = xgb.XGBClassifier(n_estimators=n_estimators, learning_rate=0.08, max_depth=5, min_child_weight=20, subsample=0.8,
                                  tree_method="hist", eval_metric="logloss", early_stopping_rounds=30, n_jobs=0, random_state=self.seed + k)
            m.fit(F[ft], y[ft], sample_weight=w[ft], eval_set=[(F[va], y[va])], sample_weight_eval_set=[w[va]], verbose=False)
            c = np.zeros((len(_PT_PB) - 1, 4, len(_PT_EB) - 1))
            if self.center:
                pi = tr & (y == 0)
                l = m.predict(F[pi], output_margin=True).astype(float)
                idx = _pt_cell(T.p.to_numpy()[pi], T.n.to_numpy()[pi], T.aeta.to_numpy()[pi])
                s, cnt = np.zeros_like(c), np.zeros_like(c)
                np.add.at(s, idx, l)
                np.add.at(cnt, idx, 1.0)
                c = np.where(cnt > 20, s / np.maximum(cnt, 1), 0.0)
            self.models.append(m)
            self.centers.append(c)
            self.info.append({"fold": k, "n_train_tracks": int(ft.sum()), "n_K": int((ft & (y == 1)).sum())})
        return self.predict(L)

    def _llr(self, L, i, fo):
        x, n = _a(L, f"d_dedx{i}"), _a(L, f"d_npts{i}")
        p, ae = track_p(L, i), np.abs(_a(L, f"d_eta{i}"))
        ok = (x > 0) & (n > 0) & np.isfinite(x) & np.isfinite(p)
        l = np.zeros(len(L))                       # brak pomiaru dE/dx → 0 (brak informacji)
        for k, m in enumerate(self.models):
            s = ok & (fo == k)
            if s.any():
                l[s] = m.predict(_pt_feats(x[s], n[s], p[s], ae[s], self.dedx_scale), output_margin=True) - self.centers[k][_pt_cell(p[s], n[s], ae[s])]
        return l

    def _h(self, F: np.ndarray, fo: np.ndarray) -> np.ndarray:
        out = np.full(len(F), np.nan, np.float32)
        for k, m in enumerate(self.models):
            s = fo == k
            if s.any():
                out[s] = m.predict_proba(F[s])[:, 1]
        return out

    def predict(self, L: pd.DataFrame, twin: np.ndarray | None = None) -> np.ndarray:
        """P(hipoteza wiersza poprawna), antysymetryczne w parze. twin = indeks bliźniaka (−1: brak) — zestaw `pid` go nie
        potrzebuje (cechy bliźniaka = swap), `pid_kin` wymaga (cosθ*), inaczej bliźniak syntetyczny z NaN cosθ*."""
        fo = self.fold(L)
        if self.mode == "pertrack":
            d = self._llr(L, 0, fo) - self._llr(L, 1, fo)
            return (1.0 / (1.0 + np.exp(-np.clip(d, -50, 50)))).astype(np.float32)
        R = row_features(L, self.dedx_scale)
        hA = self._h(R[self.feats].to_numpy(np.float32), fo)
        S = swap_features(R)
        if twin is not None and "cts" in self.feats:
            ok = twin >= 0
            S.loc[S.index[ok], "cts"] = R["cts"].to_numpy()[twin[ok]]
        hB = self._h(S[self.feats].to_numpy(np.float32), fo)
        return 0.5 * (hA + 1.0 - hB)

    # -- zapis / odczyt
    def save(self, d: str):
        os.makedirs(d, exist_ok=True)
        for k, m in enumerate(self.models):
            if self.backend == "xgboost":
                m.save_model(os.path.join(d, f"fold{k}.json"))
            else:
                with open(os.path.join(d, f"fold{k}.pkl"), "wb") as f:
                    pickle.dump(m, f)
        if self.mode == "pertrack":
            np.save(os.path.join(d, "centers.npy"), np.stack(self.centers))
        with open(os.path.join(d, "meta.json"), "w") as f:
            json.dump({"kind": "hyp_arbiter", "feats": self.feats, "feat_set": self.feat_set if isinstance(self.feat_set, str) else None,
                       "n_folds": self.n_folds, "dedx_scale": self.dedx_scale, "backend": self.backend, "seed": self.seed,
                       "mode": self.mode, "center": self.center, "folds": self.info}, f, indent=1)

    @classmethod
    def load(cls, d: str) -> "Arbiter":
        meta = json.load(open(os.path.join(d, "meta.json")))
        a = cls(meta["feats"], meta["n_folds"], meta["dedx_scale"], meta.get("seed", 12345), meta.get("mode", "pairs"),
                meta.get("center", True))
        a.backend = meta["backend"]
        if a.mode == "pertrack":
            a.centers = list(np.load(os.path.join(d, "centers.npy")))
        for k in range(a.n_folds):
            if a.backend == "xgboost":
                import xgboost as xgb
                m = xgb.XGBClassifier()
                m.load_model(os.path.join(d, f"fold{k}.json"))
            else:
                with open(os.path.join(d, f"fold{k}.pkl"), "rb") as f:
                    m = pickle.load(f)
            a.models.append(m)
        return a


# ----------------------------------------------------------------------------- integracja z charmbkg (addytywnie)
_CACHE: dict = {}


def apply_arbiter(L: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """`signal.swap_resolution: arbiter` + `signal.arbiter.model_dir` → kolumny `p_hyp` (P poprawnej hipotezy wiersza) i
    `hyp_keep` (wiersz zostaje po arbitrażu: pojedynczy albo wygrany w parze). Grupowanie par w porcji `iter_species`
    (pary na granicy porcji → traktowane jak pojedyncze, ≤ 1 zdarzenie na 2 M wierszy)."""
    a = (cfg.get("signal") or {}).get("arbiter") or {}
    d = a.get("model_dir")
    if not d:
        return L
    d = os.path.expanduser(d)
    if d not in _CACHE:
        _CACHE[d] = Arbiter.load(d)
    arb = _CACHE[d]
    if "cts" in arb.feats:
        raise SystemExit("hyp: model z cechami kinematycznymi (pid_kin) rzeźbi tło — w charmbkg dozwolony tylko zestaw `pid`")
    if not len(L):
        return L.assign(p_hyp=np.float32(0), hyp_keep=np.bool_(True))
    p = arb.predict(L)
    gid, n, rank = pair_index(L, None)
    keep = choose(L.assign(score=p), "arbiter", p, gid=gid, n=n, rank=rank)
    return L.assign(p_hyp=p.astype(np.float32), hyp_keep=keep)


def truth_from_hypo_swap(L: pd.DataFrame, cfg: dict) -> bool:
    """`swap_resolution: arbiter`: prawda swapu z `Dmeson_hypo_swap` (v1.2): truth += 100 dla hipotezy zamienionej.
    Zwraca False gdy kolumny brak / same −1 (wtedy wołający wraca do closest_mass)."""
    if "hypo_swap" not in L:
        return False
    hs = L["hypo_swap"].to_numpy(dtype=float)
    if not (np.isfinite(hs) & (hs >= 0)).any():
        return False
    t = L["truth"].to_numpy(dtype=float)
    sw = np.isin(t, cfg["signal"]["truth_values"]) & (hs == 1)
    L.loc[sw, "truth"] = L.loc[sw, "truth"] + 100
    return True


# ----------------------------------------------------------------------------- testy rzeźbienia i szablon
def sculpt_test(m: np.ndarray, keep: np.ndarray, gid: np.ndarray, n: np.ndarray, m0: float = M_D0, half: float = 0.030,
                controls=(-0.090, -0.060, 0.060, 0.090)) -> dict:
    """Test rzeźbienia na TLE (pary, n == 2): wśród par z dokładnie jednym wierszem w oknie |m − c| < half, ułamek, w którym
    metoda wybiera wiersz w oknie. H0 (brak informacji o masie): 0,5 dokładnie (rozkład dwumianowy). Okno w m0 vs okna
    kontrolne (c = m0 + δ): pik = f(m0) − średnia f(kontrola). Zwraca f, błędy, Z."""
    tw = twin_index(gid, n)
    i = np.flatnonzero((tw >= 0) & (np.arange(len(m)) < tw))
    j = tw[i]
    res = {}
    for name, c in [("m0", m0)] + [(f"c{d:+.3f}", m0 + d) for d in controls]:
        wi, wj = np.abs(m[i] - c) < half, np.abs(m[j] - c) < half
        one = wi ^ wj
        k = np.where(wi[one], keep[i][one], keep[j][one]).sum()
        N = int(one.sum())
        f = k / N if N else float("nan")
        res[name] = {"n_pairs_one_in_window": N, "f_pick_in_window": f, "err": float(0.5 / np.sqrt(N)) if N else float("nan"),
                     "z_vs_half": float((k - N / 2) / np.sqrt(N / 4)) if N else float("nan")}
    fc = [v["f_pick_in_window"] for k_, v in res.items() if k_ != "m0" and np.isfinite(v["f_pick_in_window"])]
    ec = [v["err"] for k_, v in res.items() if k_ != "m0" and np.isfinite(v["err"])]
    if fc and np.isfinite(res["m0"]["f_pick_in_window"]):
        ctrl = float(np.mean(fc))
        e = float(np.hypot(res["m0"]["err"], np.sqrt(np.sum(np.square(ec))) / len(ec)))
        res["peak"] = {"f_m0_minus_ctrl": res["m0"]["f_pick_in_window"] - ctrl, "err": e,
                       "z": (res["m0"]["f_pick_in_window"] - ctrl) / e if e > 0 else float("nan")}
    return res


def mass_hist_ratio(m, keep, gid, n, lo, hi, nb=25) -> dict:
    """Histogram m wybranych wierszy par tła vs odniesienie „ślepe" (każdy wiersz pary z wagą ½)."""
    pr = n == 2
    h_sel, xe = np.histogram(m[pr & keep], nb, (lo, hi))
    h_ref, _ = np.histogram(m[pr], nb, (lo, hi), weights=np.full(int(pr.sum()), 0.5))
    with np.errstate(all="ignore"):
        r = h_sel / h_ref
    return {"edges": xe.tolist(), "sel": h_sel.tolist(), "ref_half": h_ref.tolist(), "ratio": np.nan_to_num(r, nan=0).tolist()}
