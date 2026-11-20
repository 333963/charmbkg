"""charmbkg.bkggen.c2st — JEDEN test dwóch próbek (etap AA).

Kanoniczna implementacja: feta.c2st (S4) — HistGB OOF, wagi, klasy zbalansowane. Wcześniej były trzy:
validate.c2st_auc (S, XGBoost, bez wag), flow.auc_two_samples (S2, HistGB, bez wag), feta.c2st (S4).
Liczby S i S2 z doców NIE są porównywalne 1:1 z tym testem (S7 §3); do odtworzenia starych liczb: engine="xgb".
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _X(A, feats):
    if isinstance(A, pd.DataFrame):
        feats = list(feats) if feats is not None else [c for c in A.columns if c not in ("m", "w", "pt")]
        return A[feats].to_numpy(np.float64), feats
    return np.asarray(A, np.float64), feats


def c2st(A, B, feats=None, wa=None, wb=None, n_folds=3, max_n=150_000, seed=0, engine="hgb"):
    """AUC klasyfikatora A(1) vs B(0), out-of-fold; 0,5 = nierozróżnialne.
    A, B: DataFrame (kolumny feats) lub ndarray. wa, wb: wagi (opcjonalnie; engine="xgb" je ignoruje)."""
    XA, feats = _X(A, feats)
    XB, _ = _X(B, feats)
    if engine == "xgb":
        from .validate import _c2st_auc_xgb
        fa = feats or [f"f{i}" for i in range(XA.shape[1])]
        return _c2st_auc_xgb(pd.DataFrame(XA, columns=fa), pd.DataFrame(XB, columns=fa), fa, n_max=max_n, seed=seed, folds=n_folds)
    from .feta import c2st as _c2st
    return _c2st(XA, XB, wa=wa, wb=wb, n_folds=n_folds, max_n=max_n, seed=seed)
