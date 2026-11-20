"""charmbkg.bkggen.salad8 — SALAD (2001.05001) bez torcha: reweighter MC tła → dane, warunkowany masą (etap S8).

Rozszerza `feta.SaladReweighter` (S4) o:
  * kind='sk_mlp'  — sklearn MLPClassifier (SiLU brak → ReLU 128-128-64), wejście z-PIT + t=(m−m0)/0,1; gładki w m
                     (interpolacja SB→SR), działa bez torch/GPU (VM Maca, lxplus bez venv_s2).
  * kind='hgb' | 'hgb_nom' — jak w S4.
  * val_auc        — AUC klasyfikatora dane vs MC na 20 % walidacji SB (diagnostyka: 0,5 = MC już zgodne).
  * gatunek dowolny: wymaga tylko sp['m0'] (D0, D⁺, D*, Λc) — cechy H z charmbkg.features.compute.
  * dostępne przez bkggen.api.make_generator('salad', kind='sk_mlp') (api.Salad przekazuje kind do SaladReweighter8).
Wagi: w = p/(1−p) na próbce zbalansowanej → znormalizowane do średniej 1 (normalizację tła w SR daje fit masy).
"""
from __future__ import annotations

import logging

import numpy as np

from . import flow as F
from .feta import SaladReweighter

log = logging.getLogger("charmbkg.bkggen.salad8")


class SaladReweighter8(SaladReweighter):
    KINDS = ("sk_mlp", "hgb", "hgb_nom", "mlp")

    def __init__(self, feats, sp, kind="sk_mlp", clip=(0.05, 20.0), seed=0, hidden=(128, 128, 64), max_iter=60):
        if kind not in self.KINDS:
            raise ValueError(f"kind={kind!r}; dozwolone {self.KINDS}")
        super().__init__(feats, sp, kind=kind, clip=clip, seed=seed)
        self.hidden, self.max_iter = tuple(hidden), int(max_iter)
        self.val_auc = None

    def fit(self, X_mc, m_mc, X_data, m_data, epochs=12, max_n=400_000):
        if self.kind == "mlp":                      # torch (S4) — gdy dostępny
            return super().fit(X_mc, m_mc, X_data, m_data, epochs=epochs, max_n=max_n)
        from sklearn.metrics import roc_auc_score
        rng = np.random.default_rng(self.seed)
        self.pre = F.Preproc(self.feats, seed=self.seed).fit(X_mc)
        i_mc = rng.choice(len(X_mc), min(len(X_mc), max_n), replace=False)
        i_d = rng.choice(len(X_data), min(len(X_data), max_n), replace=False)
        A, B = self._inp(X_mc[i_mc], m_mc[i_mc]), self._inp(X_data[i_d], m_data[i_d])
        X = np.vstack([A, B]).astype(np.float32)
        y = np.r_[np.zeros(len(A)), np.ones(len(B))]
        w = np.r_[np.full(len(A), len(B) / len(A)), np.ones(len(B))]
        perm = rng.permutation(len(X)); nv = len(X) // 5
        iv, it = perm[:nv], perm[nv:]
        if self.kind == "sk_mlp":
            from sklearn.neural_network import MLPClassifier
            # równoważenie klas przez powielanie wag nie jest wspierane w MLPClassifier → losowanie z wagami
            take = it[rng.random(len(it)) < np.clip(w[it] / w[it].max(), 0, 1)]
            self.clf = MLPClassifier(hidden_layer_sizes=self.hidden, activation="relu", alpha=1e-4, batch_size=2048,
                                     learning_rate_init=1e-3, max_iter=self.max_iter, early_stopping=True,
                                     validation_fraction=0.15, n_iter_no_change=8, random_state=self.seed)
            self.clf.fit(X[take], y[take])
            self._bal = (len(A), len(B))
        else:
            from sklearn.ensemble import HistGradientBoostingClassifier
            self.clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08, max_leaf_nodes=31,
                                                      l2_regularization=1.0, early_stopping=True, random_state=self.seed)
            self.clf.fit(X[it], y[it], sample_weight=w[it])
        pv = self.clf.predict_proba(X[iv])[:, 1]
        self.val_auc = float(roc_auc_score(y[iv], pv, sample_weight=w[iv]))
        log.info("SALAD %s: val AUC dane vs MC (SB) = %.4f", self.kind, self.val_auc)
        return self

    def logit(self, X, m):
        if self.kind == "mlp":
            return super().logit(X, m)
        Z = self._inp(X, m).astype(np.float32)
        p = np.clip(self.clf.predict_proba(Z)[:, 1].astype(np.float64), 1e-6, 1 - 1e-6)
        return np.log(p / (1 - p))                  # sk_mlp: próbka `take` zbalansowana losowaniem z wagami
