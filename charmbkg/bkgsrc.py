"""Źródło tła treningowego XGB (etap S5, addytywnie). Konfiguracja: `training.bkg_source` + `training.bkg_options`.

  sideband           — DOMYŚLNE, bez zmian: sidebandy data_none |m−m0| > sideband_min (samples.training_set)
  sideband_far_veto  — dalekie sidebandy |m−m0| > far_min (D0 100, D⁺ 110 MeV) + weto refleksji D0 |m_swap−m0| > refl_veto
                       (wynik S: najlepiej oddaje tło pod pikiem, ε_B ratio 0,91/0,93/1,06 @WP80/70/50)
  mc                 — tło MC: incl_mc truth ≥ 8 (bez sygnału/swapów), cały zakres mas (mc_region: all) albo tylko SB
  synthetic          — kopuła Gaussa p(x|m) uczona na sidebandach danych (+weto refleksji), próbkowana w |m−m0| < sideband_min
                       (masy z exp(pol2) dopasowanego w SB); synthetic_model: flow → charmbkg.bkggen.gen (S/S2), gdy dostępny

Zwraca ramkę z kolumnami jak load_species (cechy, m, pt, abs_eta, weight, run, event); cap/reweight robi training_set.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger("charmbkg")

SOURCES = ("sideband", "sideband_far_veto", "mc", "synthetic")
from .constants import M_K, M_PI  # noqa: F401 (re-eksport dla refldata)
from .kinematics import four_vec as _four_vec, swapped_mass_sp as _swapped_mass_sp
from .samples import pt_range
DEFAULT_OPTS = {
    "far_min": {"D0": 0.100, "Dplus": 0.110},
    "refl_veto": 0.025,            # D0: odrzuć |m_swap − m0| < 25 MeV (0 → bez weta)
    "mc_sample": "incl_mc",
    "mc_region": "all",            # all | sideband
    "mc_truth_min": 8,
    "synthetic_model": "copula",   # copula | flow (flow wymaga charmbkg.bkggen.gen + zuko) | salad (S8: MC tła ważone, bez torcha)
    "synthetic_factor": 1.0,       # N_syn = factor × N(SB) (przed cap max_bkg_train)
    "synthetic_seed": 7,
    "synthetic_veto": True,        # weto refleksji w SB uczącym generator (D0)
}


def source(cfg) -> str:
    s = (cfg.get("training") or {}).get("bkg_source") or "sideband"
    if s not in SOURCES:
        raise ValueError(f"training.bkg_source={s!r}; dozwolone {SOURCES}")
    return s


def options(cfg) -> dict:
    o = dict(DEFAULT_OPTS)
    o.update((cfg.get("training") or {}).get("bkg_options") or {})
    return o


def _fv(pt, eta, phi, m):
    """Alias zgodności → charmbkg.kinematics.four_vec."""
    return _four_vec(pt, eta, phi, m)


def swapped_mass(L: pd.DataFrame, sp: dict) -> np.ndarray:
    """D0: masa przy zamianie hipotez K↔π (d_pt w GeV, po features.compute). Implementacja: charmbkg.kinematics."""
    return _swapped_mass_sp(L, sp)


def refl_veto_mask(L, sp, species, opts) -> np.ndarray:
    """True = przechodzi weto refleksji (tylko D0 z reflection: true i refl_veto > 0)."""
    w = float(opts.get("refl_veto") or 0)
    if species != "D0" or not sp.get("reflection") or w <= 0 or not len(L):
        return np.ones(len(L), bool)
    ms = swapped_mass(L, sp)
    return ~(np.abs(ms - sp["m0"]) < w)


def _far(opts, species, sp):
    f = opts.get("far_min")
    return float(f.get(species, sp["sideband_min"]) if isinstance(f, dict) else (f or sp["sideband_min"]))


# ----------------------------------------------------------------------------- źródła
def _data_sb(cfg, species, dmin, veto, tagx):
    from . import samples as SA
    sp = cfg["species"][species]
    o = options(cfg)
    _, n_sb = SA.count_species(cfg, "data_none", species)
    cap = cfg["io"]["max_bkg_train"]
    frac = min(1.0, 1.6 * cap / max(n_sb, 1)) if cap else 1.0

    def sel(L):
        m = SA.in_range(L, sp).to_numpy() & (np.abs(L["m"].to_numpy() - sp["m0"]) > dmin) & SA.dstar_veto_mask(L, sp)
        if veto:
            m &= refl_veto_mask(L, sp, species, o)
        return m
    return SA.load_species(cfg, "data_none", species, keep_frac=frac, tagx=tagx, select=sel)


def bkg_far_veto(cfg, species):
    sp = cfg["species"][species]
    o = options(cfg)
    dmin = _far(o, species, sp)
    B = _data_sb(cfg, species, dmin, True, f"sbfar{dmin:g}v{o['refl_veto']:g}")
    log.info("bkg_source=sideband_far_veto %s: |m−m0|>%.0f MeV, weto refl. %.0f MeV → %d kandydatów",
             species, 1e3 * dmin, 1e3 * float(o["refl_veto"] or 0), len(B))
    return B


def bkg_mc(cfg, species):
    from . import samples as SA
    sp = cfg["species"][species]
    o = options(cfg)
    key = o["mc_sample"]
    tmin = float(o["mc_truth_min"])

    def sel(L):
        t = L["truth"].to_numpy(dtype=float)
        m = SA.in_range(L, sp).to_numpy() & (t >= tmin) & (t < 100) & SA.dstar_veto_mask(L, sp)
        if o["mc_region"] == "sideband":
            m &= SA.sideband_mask(L, sp)
        return m
    B = SA.load_species(cfg, key, species, tagx=f"mcbkg{tmin:g}{o['mc_region']}", select=sel)
    log.info("bkg_source=mc %s: %s truth≥%g region=%s → %d kandydatów", species, key, tmin, o["mc_region"], len(B))
    return B


class Copula:
    """Kopuła Gaussa z|m: Gaussianizacja kwantylami (drganie atomów), μ(m), Σ(m) liniowe w (m−m0) z binów masy SB."""

    def __init__(self, m0, n_q=2001, nbins_m=6, seed=0):
        self.m0, self.n_q, self.nb, self.rng = m0, n_q, nbins_m, np.random.default_rng(seed)

    def _to_u(self, X):
        U = np.empty_like(X)
        for j in range(X.shape[1]):
            q = self.q[:, j]
            lo = np.searchsorted(q, X[:, j], "left")
            hi = np.searchsorted(q, X[:, j], "right")
            U[:, j] = (lo + self.rng.random(len(X)) * np.maximum(hi - lo, 1)) / len(q)
        return np.clip(U, 0.5 / self.n_q, 1 - 0.5 / self.n_q)

    def fit(self, X, m):
        from scipy.special import ndtri
        X = np.asarray(X, float)
        self.med = np.nanmedian(X, 0)
        self.nanf = np.mean(~np.isfinite(X), 0)
        X = np.where(np.isfinite(X), X, self.med)
        self.p = np.linspace(0, 1, self.n_q)
        self.q = np.quantile(X, self.p, axis=0)
        Z = ndtri(self._to_u(X))
        t = (np.asarray(m, float) - self.m0) / 0.1
        e = np.unique(np.quantile(t, np.linspace(0, 1, self.nb + 1)))
        tc, mus, covs = [], [], []
        for a, b in zip(e[:-1], e[1:]):
            s = (t >= a) & (t <= b)
            if s.sum() < 50:
                continue
            tc.append(t[s].mean())
            mus.append(Z[s].mean(0))
            covs.append(np.cov(Z[s], rowvar=False))
        tc = np.array(tc)
        A = np.c_[np.ones_like(tc), tc] if len(tc) > 1 else np.c_[np.ones_like(tc), np.zeros_like(tc)]
        self.mu = np.linalg.lstsq(A, np.array(mus), rcond=None)[0]
        self.cv = np.linalg.lstsq(A, np.array(covs).reshape(len(tc), -1), rcond=None)[0]
        self.d = X.shape[1]
        return self

    def sample(self, m):
        from scipy.special import ndtr
        t = (np.asarray(m, float) - self.m0) / 0.1
        Z = np.empty((len(t), self.d))
        g = np.linspace(t.min(), t.max() + 1e-9, 21)
        ib = np.clip(np.digitize(t, g) - 1, 0, 19)
        for k in np.unique(ib):
            s = ib == k
            tc = 0.5 * (g[k] + g[k + 1])
            mu = self.mu[0] + self.mu[1] * tc
            C = (self.cv[0] + self.cv[1] * tc).reshape(self.d, self.d)
            w, V = np.linalg.eigh(0.5 * (C + C.T))
            Z[s] = mu + (self.rng.standard_normal((s.sum(), self.d)) * np.sqrt(np.clip(w, 1e-4, None))) @ V.T
        U = ndtr(Z)
        X = np.column_stack([np.interp(U[:, j], self.p, self.q[:, j]) for j in range(self.d)])
        miss = self.rng.random(X.shape) < self.nanf  # zachowaj udział braków (np. dE/dx)
        return np.where(miss, np.nan, X)


def _mass_sampler(m_sb, lo, hi, m0, dmin, rng):
    """exp(a0 + a1 x + a2 x²) dopasowane (Poisson) do histogramu SB; próbki w |m−m0| < dmin."""
    from iminuit import Minuit
    h, e = np.histogram(m_sb, 50, (lo, hi))
    c = 0.5 * (e[1:] + e[:-1])
    use = np.abs(c - m0) > dmin
    x = (c - lo) / (hi - lo)

    def nll(a0, a1, a2):
        mu = np.exp(a0 + a1 * x[use] + a2 * x[use] ** 2)
        return float(np.sum(mu - h[use] * np.log(mu + 1e-300)))
    mi = Minuit(nll, a0=np.log(max(h[use].mean(), 1)), a1=0.0, a2=0.0)
    mi.errordef = Minuit.LIKELIHOOD
    mi.migrad()
    a = np.array(mi.values)
    g = np.linspace(max(lo, m0 - dmin), min(hi, m0 + dmin), 4001)
    xg = (g - lo) / (hi - lo)
    cdf = np.cumsum(np.exp(a[0] + a[1] * xg + a[2] * xg ** 2))
    cdf /= cdf[-1]
    n_exp = float(np.trapezoid(np.exp(a[0] + a[1] * xg + a[2] * xg ** 2), g) / (e[1] - e[0]))
    return (lambda n: np.interp(rng.random(n), cdf, g)), n_exp, a.tolist()


def bkg_synthetic(cfg, species, feats=None):
    from . import features as F
    sp = cfg["species"][species]
    o = options(cfg)
    feats = list(feats or F.feature_list(cfg, cfg["training"].get("variant", "H")))
    if o["synthetic_model"] == "salad":          # S8: MC tła ważone SALAD (nie generacja) — charmbkg/s8salad.py
        from .s8salad import training_background_salad
        return training_background_salad(cfg, species, feats)
    SB = _data_sb(cfg, species, sp["sideband_min"], bool(o["synthetic_veto"]), f"synsb{o['refl_veto']:g}")
    feats = [f for f in feats if f in SB and np.isfinite(SB[f].to_numpy(dtype=float)).any()]
    gen_cols = list(dict.fromkeys(feats + ["abs_eta", "log_pt"]))
    rng = np.random.default_rng(int(o["synthetic_seed"]))
    lo, hi = sp["mass_range"]
    msamp, n_exp, coef = _mass_sampler(SB["m"].to_numpy(float), lo, hi, sp["m0"], sp["sideband_min"], rng)
    n = int(round(float(o["synthetic_factor"]) * len(SB)))
    cap = cfg["io"].get("max_bkg_train")
    n = min(n, int(1.3 * cap)) if cap else n
    m = msamp(n)
    if o["synthetic_model"] == "flow":
        from .bkggen.gen import BkgGenerator  # moduł S/S2 (opcjonalny)
        G = BkgGenerator(gen_cols, sp["m0"], lambda x: np.abs(x - sp["m0"]) > sp["sideband_min"], (lo, hi), model="flow",
                         seed=int(o["synthetic_seed"]))
        G.fit(SB)
        X = G.sample(n, None, m=m)[gen_cols].to_numpy(float)
    else:
        X = Copula(sp["m0"], seed=int(o["synthetic_seed"])).fit(SB[gen_cols].to_numpy(float), SB["m"].to_numpy(float)).sample(m)
    B = pd.DataFrame(X, columns=gen_cols)
    B["m"] = m
    B["pt"] = np.expm1(B["log_pt"].to_numpy())
    B["abs_eta"] = np.abs(B["abs_eta"].to_numpy())
    B["eta"] = B["abs_eta"]
    B["weight"] = 1.0
    B["run"] = 0  # syntetyczne: własne „zdarzenia” (fold z hash event)
    B["event"] = rng.integers(0, 2**62, n)
    B["truth"] = -1.0
    B["sample"] = "synthetic"
    plo, phi = pt_range()
    B = B[(B["pt"] >= plo) & (B["pt"] < phi)].reset_index(drop=True)
    B.attrs["synthetic"] = {"n_sb": int(len(SB)), "n_syn": int(len(B)), "n_exp_window_from_sb_fit": n_exp, "mass_coef": coef,
                            "model": o["synthetic_model"], "features": gen_cols}
    log.info("bkg_source=synthetic %s: model=%s, SB=%d → syn=%d w |m−m0|<%.0f MeV (oczek. tło w oknie z fitu SB %.0f)",
             species, o["synthetic_model"], len(SB), len(B), 1e3 * sp["sideband_min"], n_exp)
    return B


def training_background(cfg, species):
    s = source(cfg)
    if s == "sideband_far_veto":
        return bkg_far_veto(cfg, species)
    if s == "mc":
        return bkg_mc(cfg, species)
    if s == "synthetic":
        return bkg_synthetic(cfg, species)
    raise ValueError(s)
