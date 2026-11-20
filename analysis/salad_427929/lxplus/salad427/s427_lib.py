"""SALAD427 — wspólne definicje paczki Q_lxplus/salad427 (projekt B, 17.11.2026).

SALAD v1.0 (zamrożony yaml `charmbkg/data/salad_v1.yaml`, wersja pliku 1.1: D0 truth >= 2) powtórzony na runie 427929
(physics_MinBias, MBTS) i w klasach trygera runu 428770. Tylko pasma boczne i tło: żadnych fitów piku sygnału w danych.
Używane w s427_extract.py (lxplus), s427_study.py (lxplus) i SALAD_427929/analyze.py (lokalnie, tylko stałe).
Pary D-D poza zakresem. charmbkg używany jako biblioteka (kopia w src/, bez edycji).
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VERSION = "salad427-1.0"
SPECIES = ("D0", "Dplus")
PT_BINS = [3.0, 5.0, 8.0, 30.0]
BIN_NAMES = ["3-5", "5-8", "8-30", "all"]
WPS = (0.95, 0.80, 0.70, 0.50)
K_FOLDS = 5
REFL_VETO = 0.025                 # D0: |m_swap - m0| > 25 MeV wszędzie (B1-D0 / S4)
LXY_MIN = 0.25                    # mm, nominał v1.1/v1.2 (h5 v2.0 ma 0,10)
MASS_BW = 0.005                   # biny masy 5 MeV do interpolacji pasm bocznych
# budżet systematyki kształtu tła (METHODS §5.2): |SALAD/dane - 1| przy WP80/70/50; WP95 bez budżetu
BUDGET = {"Dplus": {0.80: 0.05, 0.70: 0.10, 0.50: 0.25}, "D0": {0.80: 0.05, 0.70: 0.10, 0.50: 0.15}}
BUDGET_D0_T8_WP50 = 0.20          # D0 z MC truth >= 8 (wariant kontrolny)

MBTS_CHAINS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1",
               "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
J12 = "HLT_noalg_L1J12"
AFP_CHAINS = ["HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"]
CHAINS = MBTS_CHAINS + [J12] + AFP_CHAINS            # kolejność = numer bitu w `tbits`
GRL_COL = "passGRL_MinBiasAFP"
UID_COL = "event_uid"
# klasy (rozłączne w 428770; definicje jak B_trig/btrig_common i B_trig/AFPSTREAM):
#   MB427929  = run 427929, GRL i HLT_mb_mbts_L1MBTS_1
#   J12only   = 428770, J12 i żaden łańcuch MBTS
#   MBTSbased = 428770, dowolny łańcuch L1MBTS
#   AFPonly   = 428770, bit AFP (noalg lub sptrk_pt4 L1AFP_A_OR_C), bez J12 i bez MBTS („strumień AFP”)
CLASSES = ["MB427929", "J12only", "MBTSbased", "AFPonly"]
CLASS_SAMPLE = {"MB427929": "data_427929", "J12only": "data_428770", "MBTSbased": "data_428770", "AFPonly": "data_428770"}


def add_src():
    src = os.environ.get("CHARMBKG_SRC") or os.path.join(HERE, "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    return src


# ----------------------------------------------------------------------------------------------------- hash / foldy
def mix64(x, salt=0):
    """splitmix64 -> uint64 (deterministyczny hash zdarzenia)."""
    with np.errstate(over="ignore"):
        x = np.asarray(x, np.uint64) + np.uint64(0x9E3779B97F4A7C15) * np.uint64(salt + 1)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return x ^ (x >> np.uint64(31))


def unit_hash(x, salt=0):
    return (mix64(x, salt) >> np.uint64(11)).astype(np.float64) / float(1 << 53)


def fold_of(uid, k=K_FOLDS, salt=7):
    return (mix64(uid, salt) % np.uint64(k)).astype(np.int8)


# ----------------------------------------------------------------------------------------------------- klasy trygera
def bit(tbits, chain):
    return (np.asarray(tbits).astype(np.int64) >> CHAINS.index(chain)) & 1 == 1


def class_masks(tbits, grl, sample):
    """{klasa: maska} dla próbki danych (`data_427929` albo `data_428770`)."""
    m = np.zeros(len(tbits), bool)
    for c in MBTS_CHAINS:
        m |= bit(tbits, c)
    j = bit(tbits, J12)
    a = bit(tbits, AFP_CHAINS[0]) | bit(tbits, AFP_CHAINS[1])
    if sample == "data_427929":
        return {"MB427929": bit(tbits, "HLT_mb_mbts_L1MBTS_1") & (np.asarray(grl) == 1)}
    return {"J12only": j & ~m, "MBTSbased": m, "AFPonly": a & ~j & ~m}


# ----------------------------------------------------------------------------------------------------- regiony / biny
def ptbin(pt):
    return np.clip(np.digitize(pt, PT_BINS) - 1, -1, len(PT_BINS) - 1)


class MassGrid:
    """Biny masy 5 MeV w oknie gatunku; SB = środek binu >= sb_min od m0, SR = środek < sr_half (jak BKG_427929/bkg_common)."""

    def __init__(self, sp):
        self.m0, (self.lo, self.hi) = sp["m0"], sp["mass_range"]
        self.sr_half, self.sb_min = sp["sr_half"], sp["sb_min"]
        self.nb = int(round((self.hi - self.lo) / MASS_BW))
        self.edges = np.linspace(self.lo, self.hi, self.nb + 1)
        self.cen = 0.5 * (self.edges[1:] + self.edges[:-1])
        self.t = (self.cen - self.m0) / 0.1
        self.sb = np.abs(self.cen - self.m0) >= self.sb_min
        self.bw = self.edges[1] - self.edges[0]
        self.xf = np.linspace(self.m0 - self.sr_half, self.m0 + self.sr_half, 401)
        self.tf = (self.xf - self.m0) / 0.1
        self.sb_max = min(self.m0 - self.lo, self.hi - self.m0)      # zasięg „symetrycznego” SB

    def cat_hist(self, m, ncat_pass, w=None):
        """Histogram [kategoria = liczba spełnionych WP (0..nWP), bin masy]; pass(WP j) = suma kategorii >= j+1."""
        im = np.clip(((np.asarray(m) - self.lo) / self.bw).astype(int), 0, self.nb - 1)
        H = np.zeros((len(WPS) + 1, self.nb))
        np.add.at(H, (np.asarray(ncat_pass, int), im), 1.0 if w is None else w)
        return H

    def glm_pois(self, h, order):
        """Poisson ML: N(bin) = exp(sum p_k t^k) na binach SB (Newton). None gdy brak zliczeń."""
        X = np.vander(self.t[self.sb], order + 1, increasing=True)
        y = h[self.sb]
        if y.sum() <= 0:
            return None
        p = np.zeros(order + 1)
        p[0] = np.log(max(y.mean(), 1e-9))
        for _ in range(100):
            mu = np.exp(np.clip(X @ p, -50, 50))
            g = X.T @ (y - mu)
            Hs = X.T @ (X * mu[:, None])
            try:
                dp = np.linalg.solve(Hs + 1e-10 * np.eye(order + 1), g)
            except np.linalg.LinAlgError:
                return None
            st = np.abs(dp).max()
            if st > 5:
                dp *= 5 / st
            p = p + dp
            if st < 1e-9:
                break
        return p

    def pred_sr(self, p):
        if p is None:
            return np.nan
        dens = np.exp(np.clip(np.vander(self.tf, len(p), increasing=True) @ p, -50, 50)) / self.bw
        return float(np.sum(0.5 * (dens[1:] + dens[:-1]) * np.diff(self.xf)))

    def n_sr(self, h, order=2):
        return self.pred_sr(self.glm_pois(h, order))

    def eps_expk(self, h_all, h_pass, order):
        a, b = self.n_sr(h_all, order), self.n_sr(h_pass, order)
        return b / a if a > 0 else np.nan

    def eps_sym(self, h_all, h_pass):
        s = self.sb & (np.abs(self.cen - self.m0) <= self.sb_max)
        return h_pass[s].sum() / max(h_all[s].sum(), 1e-12)

    def eps_lr(self, h_all, h_pass):
        L, R = self.sb & (self.cen < self.m0), self.sb & (self.cen > self.m0)
        if h_all[L].sum() <= 0 or h_all[R].sum() <= 0:
            return np.nan
        eL, eR = h_pass[L].sum() / h_all[L].sum(), h_pass[R].sum() / h_all[R].sum()
        mL, mR = np.average(self.cen[L], weights=h_all[L]), np.average(self.cen[R], weights=h_all[R])
        return eL + (eR - eL) * (self.m0 - mL) / (mR - mL)

    def eps_all(self, H, method):
        """H = cat_hist -> wektor eps_B(SR) po WP dla metody pol2|pol1|sym|lr."""
        ha = H.sum(0)
        out = []
        for j in range(len(WPS)):
            hp = H[j + 1:].sum(0)
            out.append({"pol2": lambda: self.eps_expk(ha, hp, 2), "pol1": lambda: self.eps_expk(ha, hp, 1),
                        "sym": lambda: self.eps_sym(ha, hp), "lr": lambda: self.eps_lr(ha, hp)}[method]())
        return np.array(out, float)


def neff(w):
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / max((w ** 2).sum(), 1e-300)) if len(w) else 0.0


def weff(w, k):
    """eps = sum(w k)/sum(w) i błąd (ważony dwumianowy): sqrt(sum w^2 (k - eps)^2)/sum w."""
    w = np.asarray(w, float)
    k = np.asarray(k, float)
    s = w.sum()
    if s <= 0:
        return float("nan"), float("nan")
    e = float((w * k).sum() / s)
    return e, float(np.sqrt((w ** 2 * (k - e) ** 2).sum()) / s)
