"""charmbkg.bkggen.salad_v1 — SALAD v1.0 (etap B1): zamrożony, produkcyjny estymator ε_B / kształtu tła pod pikiem.

Idea (SALAD, 2001.05001; S4/S8): klasyfikator dane(1) vs MC tła(0) uczony w sidebandach z masą jako wejściem
(t = (m−m0)/t_scale), w(x, m) = p/(1−p) interpolowane do SR. Zwraca PRAWDZIWE zdarzenia MC z wagą — nie generuje.

API (wszystko z jednego pliku konfiguracji `config/bkggen/salad_v1.yaml`):
    cfg = load_cfg()                                       # lub load_cfg(path)
    rw  = fit_salad_v1(cfg, "Dplus", data_sb, mc)          # DataFrame'y z cechami H + "m" (+ "truth" w mc)
    w   = salad_weights(rw, mc)                            # wagi dla WSZYSTKICH wierszy mc (średnia 1 w całym zakresie)
    path = run(cfg, "Dplus", data_sb, mc, out_dir, sample="data")   # → salad_w_Dplus_data_v1.npz (+ diagnostyka)
    w, meta = load_weights(path)
Diagnostyka w rw.diag / meta: val_auc_SB, neff_over_N_SR, c2st_band (jeśli data_sb zawiera pasmo), sha256 wejść.
Założenie (sprawdzone w B1): skład tła zmienia się gładko w m; składnik zlokalizowany pod pikiem nie jest korygowany.
Addytywne: nie zmienia S4/S8; rejestruje generator api "salad_v1" (fit(SB_data, mc=MC) / sample(m) → MC w SR z wagą `w`).
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import time

import numpy as np
import pandas as pd

from . import flow as F
from .feta import SaladReweighter

log = logging.getLogger("charmbkg.bkggen.salad_v1")
VERSION = "1.1"   # 1.1 (B1-D0/B4): mc_truth_min per gatunek (D0: 2 — truth 2–7 = 1 % tła, nieskorygowane przy ≥ 8)
HERE = os.path.dirname(os.path.abspath(__file__))
_PKG = os.path.dirname(HERE)
_CFG_CANDIDATES = (  # kopia w pakiecie (package-data; REPO) → config/ (AA_merge) → configs/ (layout repo)
    os.path.join(_PKG, "data", "salad_v1.yaml"),
    os.path.join(os.path.dirname(_PKG), "config", "bkggen", "salad_v1.yaml"),
    os.path.join(os.path.dirname(_PKG), "configs", "bkggen", "salad_v1.yaml"),
)
DEFAULT_CFG = next((c for c in _CFG_CANDIDATES if os.path.exists(c)), _CFG_CANDIDATES[0])
KINDS = ("hgb", "sk_mlp")


# --------------------------------------------------------------------------- konfiguracja / regiony
def load_cfg(path: str | None = None) -> dict:
    import yaml
    with open(path or DEFAULT_CFG) as f:
        cfg = yaml.safe_load(f)
    cfg.setdefault("version", VERSION)
    return cfg


def species_params(cfg: dict, species: str) -> dict:
    r = cfg["regions"][species]
    return dict(m0=float(r["m0"]), mass_range=tuple(r["mass_range"]), sr_half=float(r["sr_half"]),
                band=tuple(r["band"]), sb_min=float(r["sb_min"]), t_scale=float(cfg.get("t_scale", 0.1)))


def regions(m, sp: dict) -> np.ndarray:
    """'SR' | 'band' | 'SB' | 'out' dla każdej masy."""
    m = np.asarray(m, float)
    dm = np.abs(m - sp["m0"])
    lo, hi = sp["mass_range"]
    out = np.where(dm < sp["sr_half"], "SR", np.where(dm < sp["sb_min"], "band", "SB"))
    out = np.where((m < lo) | (m > hi), "out", out)
    return out


def input_hash(df: pd.DataFrame, cols) -> str:
    """sha256 kolumn (float32, kolejność cols) — identyfikacja wejścia w artefakcie wag."""
    h = hashlib.sha256()
    for c in cols:
        h.update(np.ascontiguousarray(df[c].to_numpy(dtype=np.float32)).tobytes())
    return h.hexdigest()


# --------------------------------------------------------------------------- reweighter
class SaladV1(SaladReweighter):
    """kind='hgb' (HistGB z t) | 'sk_mlp' (MLP na z-PIT + t). Wagi = p/(1−p) na próbce zbalansowanej, clip, średnia 1."""

    def __init__(self, feats, sp, cfg, kind=None, seed=None):
        kind = kind or cfg["kind"]
        if kind not in KINDS:
            raise ValueError(f"kind={kind!r}; dozwolone {KINDS}")
        super().__init__(feats, sp, kind=kind, clip=tuple(cfg.get("clip", (0.05, 20.0))),
                         seed=int(cfg["seed"] if seed is None else seed))
        self.cfg = cfg
        self.diag = {"kind": kind, "version": cfg.get("version", VERSION)}

    def _inp(self, X, m):
        Z = self.pre.forward(X)
        t = (np.asarray(m, float) - self.sp["m0"]) / self.sp.get("t_scale", 0.1)
        return np.column_stack([Z, t])

    def fit(self, X_mc, m_mc, X_data, m_data, max_n=None, **_):
        from sklearn.metrics import roc_auc_score
        t0 = time.time()
        max_n = int(max_n or self.cfg.get("max_n", 400_000))
        rng = np.random.default_rng(self.seed)
        X_mc, X_data = np.asarray(X_mc, float), np.asarray(X_data, float)
        self.pre = F.Preproc(self.feats, seed=self.seed).fit(X_mc)
        i_mc = rng.choice(len(X_mc), min(len(X_mc), max_n), replace=False)
        i_d = rng.choice(len(X_data), min(len(X_data), max_n), replace=False)
        A, B = self._inp(X_mc[i_mc], np.asarray(m_mc)[i_mc]), self._inp(X_data[i_d], np.asarray(m_data)[i_d])
        X = np.vstack([A, B]).astype(np.float32)
        y = np.r_[np.zeros(len(A)), np.ones(len(B))]
        w = np.r_[np.full(len(A), len(B) / len(A)), np.ones(len(B))]
        perm = rng.permutation(len(X)); nv = len(X) // 5
        iv, it = perm[:nv], perm[nv:]
        if self.kind == "sk_mlp":
            from sklearn.neural_network import MLPClassifier
            p = self.cfg.get("sk_mlp", {})
            take = it[rng.random(len(it)) < np.clip(w[it] / w[it].max(), 0, 1)]     # balans klas losowaniem
            self.clf = MLPClassifier(hidden_layer_sizes=tuple(p.get("hidden", (128, 128, 64))), activation="relu",
                                     alpha=float(p.get("alpha", 1e-4)), batch_size=int(p.get("batch_size", 2048)),
                                     learning_rate_init=1e-3, max_iter=int(p.get("max_iter", 60)), early_stopping=True,
                                     validation_fraction=0.15, n_iter_no_change=8, random_state=self.seed)
            self.clf.fit(X[take], y[take])
        else:
            from sklearn.ensemble import HistGradientBoostingClassifier
            p = self.cfg.get("hgb", {})
            self.clf = HistGradientBoostingClassifier(max_iter=int(p.get("max_iter", 300)),
                                                      learning_rate=float(p.get("learning_rate", 0.08)),
                                                      max_leaf_nodes=int(p.get("max_leaf_nodes", 31)),
                                                      l2_regularization=float(p.get("l2_regularization", 1.0)),
                                                      early_stopping=True, random_state=self.seed)
            self.clf.fit(X[it], y[it], sample_weight=w[it])
        pv = self.clf.predict_proba(X[iv])[:, 1]
        self.val_auc = float(roc_auc_score(y[iv], pv, sample_weight=w[iv]))
        self.diag.update(val_auc_SB=self.val_auc, n_fit_mc=int(len(A)), n_fit_data=int(len(B)), fit_s=round(time.time() - t0, 1))
        log.info("SALAD v1 %s: val AUC dane vs MC (SB) = %.4f (%.0f s)", self.kind, self.val_auc, time.time() - t0)
        return self

    def logit(self, X, m):
        Z = self._inp(np.asarray(X, float), m).astype(np.float32)
        p = np.clip(self.clf.predict_proba(Z)[:, 1].astype(np.float64), 1e-6, 1 - 1e-6)
        return np.log(p / (1 - p))


# --------------------------------------------------------------------------- API wysokiego poziomu
def _feats(cfg, *dfs):
    feats = [f for f in cfg["features"] if all(f in d for d in dfs)]
    feats = [f for f in feats if all(np.isfinite(d[f].to_numpy(dtype=float)).any() for d in dfs)]
    miss = [f for f in cfg["features"] if f not in feats]
    if miss:
        log.warning("SALAD v1: pominięte cechy (brak/NaN): %s", miss)
    return feats


def mc_truth_min(cfg: dict, species: str | None = None) -> float:
    """Próg truth tła MC: `mc_truth_min_species[species]` (v1.1), inaczej `mc_truth_min` (domyślnie 8)."""
    per = cfg.get("mc_truth_min_species") or {}
    return float(per.get(species, cfg.get("mc_truth_min", 8)) if species else cfg.get("mc_truth_min", 8))


def _mc_bkg(cfg, mc, species=None):
    if "truth" in mc:
        return mc[mc["truth"].to_numpy(dtype=float) >= mc_truth_min(cfg, species)].reset_index(drop=True)
    return mc.reset_index(drop=True)


def fit_salad_v1(cfg: dict, species: str, data_sb: pd.DataFrame, mc: pd.DataFrame, kind=None, seed=None,
                 c2st_band=True) -> SaladV1:
    """Fit na |m−m0| > sb_min obu próbek (reszta wierszy data_sb/mc ignorowana w ficie). Jeśli data_sb zawiera pasmo
    → C2ST(dane pasmo, MC pasmo ważone) jako jedyny test bez zaglądania w SR. mc może zawierać `truth` (filtr ≥ mc_truth_min)."""
    sp = species_params(cfg, species)
    mcb = _mc_bkg(cfg, mc, species)
    feats = _feats(cfg, data_sb, mcb)
    rD, rM = regions(data_sb["m"], sp), regions(mcb["m"], sp)
    nD = int((rD == "SB").sum())
    if nD < int(cfg.get("min_n_data_sb", 0)):
        log.warning("SALAD v1: mało danych w SB: %d < %d", nD, cfg.get("min_n_data_sb"))
    rw = SaladV1(feats, sp, cfg, kind=kind, seed=seed).fit(
        mcb.loc[rM == "SB", feats].to_numpy(float), mcb.loc[rM == "SB", "m"].to_numpy(float),
        data_sb.loc[rD == "SB", feats].to_numpy(float), data_sb.loc[rD == "SB", "m"].to_numpy(float))
    rw.feats_fit = feats
    rw.species = species
    w = rw.weights(mcb[feats].to_numpy(float), mcb["m"].to_numpy(float))
    sr = rM == "SR"
    rw.diag.update(species=species, n_data_SB=nD, n_mc=int(len(mcb)), n_mc_SR=int(sr.sum()),
                   neff_over_N_SR=float(w[sr].sum() ** 2 / (w[sr] ** 2).sum() / max(sr.sum(), 1)) if sr.any() else None,
                   w_quantiles_SR=[float(x) for x in np.quantile(w[sr], [0.01, 0.1, 0.5, 0.9, 0.99])] if sr.any() else None,
                   sha256_mc=input_hash(mcb, feats + ["m"]), sha256_data_sb=input_hash(data_sb, feats + ["m"]))
    if c2st_band and (rD == "band").sum() > 1000 and (rM == "band").sum() > 1000:
        from .feta import c2st
        c = cfg.get("c2st", {})
        rw.diag["c2st_band"] = c2st(data_sb.loc[rD == "band", feats].to_numpy(float), mcb.loc[rM == "band", feats].to_numpy(float),
                                    None, w[rM == "band"], n_folds=int(c.get("n_folds", 3)), max_n=int(c.get("max_n", 120_000)),
                                    seed=rw.seed)
    rw.diag["ok"] = bool((rw.diag.get("neff_over_N_SR") or 0) >= float(cfg.get("min_neff_over_N_SR", 0))
                         and nD >= int(cfg.get("min_n_data_sb", 0))
                         and rw.diag.get("c2st_band", 0.5) <= float(cfg.get("max_c2st_band", 1.0)))
    return rw


def salad_weights(rw: SaladV1, mc: pd.DataFrame | np.ndarray, m=None) -> np.ndarray:
    """Wagi dla wierszy mc (DataFrame z cechami + 'm', lub X ndarray + m). Średnia 1 na podanej próbce."""
    if isinstance(mc, pd.DataFrame):
        return rw.weights(mc[rw.feats_fit].to_numpy(float), mc["m"].to_numpy(float))
    return rw.weights(np.asarray(mc, float), np.asarray(m, float))


def save_weights(path: str, w: np.ndarray, meta: dict, row=None) -> str:
    meta = {**meta, "version": meta.get("version", VERSION), "saved": time.strftime("%Y-%m-%dT%H:%M:%S")}
    kw = dict(w=np.asarray(w, np.float32), meta=json.dumps(meta))
    if row is not None:
        kw["row"] = np.asarray(row)
    np.savez_compressed(path, **kw)
    return path


def load_weights(path: str):
    z = np.load(path, allow_pickle=False)
    return z["w"].astype(np.float64), json.loads(str(z["meta"]))


def run(cfg: dict, species: str, data_sb: pd.DataFrame, mc: pd.DataFrame, out_dir: str, sample: str = "data",
        kind=None, seed=None) -> str:
    """Fit + wagi dla mc (po filtrze truth) → <out_dir>/salad_w_<species>_<sample>_v1.npz; kolumna `row` mc (jeśli jest)."""
    os.makedirs(out_dir, exist_ok=True)
    rw = fit_salad_v1(cfg, species, data_sb, mc, kind=kind, seed=seed)
    mcb = _mc_bkg(cfg, mc, species)
    w = salad_weights(rw, mcb)
    meta = {**rw.diag, "feats": rw.feats_fit, "clip": list(rw.clip), "seed": rw.seed, "sample": sample,
            "cfg": {**{k: cfg[k] for k in ("kind", "t_scale", "clip", "max_n", "mc_truth_min") if k in cfg},
                    "mc_truth_min_used": mc_truth_min(cfg, species)}}
    path = os.path.join(out_dir, f"salad_w_{species}_{sample}_v{str(cfg.get('version', VERSION)).split('.')[0]}.npz")
    save_weights(path, w, meta, row=mcb["row"].to_numpy() if "row" in mcb else None)
    log.info("SALAD v1 %s/%s: %s  n_eff/N(SR)=%.2f  AUC_SB=%.3f  C2ST pasmo=%s  ok=%s", species, sample, path,
             meta.get("neff_over_N_SR") or -1, meta.get("val_auc_SB") or -1, meta.get("c2st_band"), meta["ok"])
    return path


# --------------------------------------------------------------------------- generator api "salad_v1"
from . import api as _api  # noqa: E402


@_api.register("salad_v1")
class SaladV1Generator(_api.Salad):
    """fit(SB_data, mc=MC_all[, cfg=..., kind=...]) → sample(m): zdarzenia MC z |m−m0|<sr_half z wagą `w` (SALAD v1.0)."""

    def _fit(self, SB, mc=None, cfg=None, kind=None, **kw):
        if mc is None:
            raise ValueError("salad_v1.fit wymaga mc=DataFrame (MC tła, cały zakres mas)")
        cfg = cfg or load_cfg()
        self.rw = fit_salad_v1(cfg, self.species, SB, mc, kind=kind, seed=self.seed, c2st_band=False)
        self._feats_fit = self.feats_out = list(self.rw.feats_fit)
        self.mc = _mc_bkg(cfg, mc, self.species)
        self.mass_model = F.MassModel(self.sp).fit(SB["m"].to_numpy(float))
