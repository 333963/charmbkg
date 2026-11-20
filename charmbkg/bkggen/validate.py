"""Metryki zgodności kształtów: χ²/ndf (znormalizowane histogramy, ważone), KS (ważony), AUC testu dwupróbkowego (C2ST)."""
from __future__ import annotations

import numpy as np


def _w(x, w):
    return np.ones(len(x)) if w is None else np.asarray(w, float)


def chi2_ndf(a, b, wa=None, wb=None, bins=30, rng=None):
    """χ²/ndf porównania kształtów (oba znormalizowane do 1), błędy Σw² w obu. Biny = kwantyle a (równa statystyka)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    wa, wb = _w(a, wa), _w(b, wb)
    fa, fb = np.isfinite(a), np.isfinite(b)
    a, wa, b, wb = a[fa], wa[fa], b[fb], wb[fb]
    if len(a) < 50 or len(b) < 50:
        return float("nan")
    if rng is None:
        e = np.unique(np.quantile(a, np.linspace(0, 1, bins + 1)))
    else:
        e = np.linspace(*rng, bins + 1)
    if len(e) < 3:
        return float("nan")
    ha, _ = np.histogram(np.clip(a, e[0], e[-1]), e, weights=wa)
    hb, _ = np.histogram(np.clip(b, e[0], e[-1]), e, weights=wb)
    va, _ = np.histogram(np.clip(a, e[0], e[-1]), e, weights=wa ** 2)
    vb, _ = np.histogram(np.clip(b, e[0], e[-1]), e, weights=wb ** 2)
    sa, sb = ha.sum(), hb.sum()
    if sa <= 0 or sb <= 0:
        return float("nan")
    ha, hb, va, vb = ha / sa, hb / sb, va / sa ** 2, vb / sb ** 2
    m = (va + vb) > 0
    return float(np.sum((ha[m] - hb[m]) ** 2 / (va[m] + vb[m])) / max(m.sum() - 1, 1))


def chi2_ndf_hist(ha, hb, va=None, vb=None, normalize=True, min_var=1.0):
    """AA (S7 §2.2): χ²/ndf dwóch HISTOGRAMÓW z błędami Poissona. va/vb = wariancje per bin (domyślnie = h, Poisson);
    dla szablonu z odejmowaniem podać Σw² (reflection_template_var). normalize: oba kształty do 1 (ndf = n−1).
    Zastępuje błędne wywołanie chi2_ndf(centra, centra, weights=h) (σ = h → χ² zaniżony ~4×)."""
    ha, hb = np.asarray(ha, float), np.asarray(hb, float)
    va = np.maximum(ha, min_var) if va is None else np.asarray(va, float)
    vb = np.maximum(hb, min_var) if vb is None else np.asarray(vb, float)
    sa, sb = ha.sum(), hb.sum()
    if sa <= 0 or sb <= 0:
        return float("nan")
    if normalize:
        ha, va, hb, vb = ha / sa, va / sa ** 2, hb / sb, vb / sb ** 2
    m = (va + vb) > 0
    ndf = max(int(m.sum()) - (1 if normalize else 0), 1)
    return float(np.sum((ha[m] - hb[m]) ** 2 / (va[m] + vb[m])) / ndf)


def ks(a, b, wa=None, wb=None):
    a, b = np.asarray(a, float), np.asarray(b, float)
    wa, wb = _w(a, wa), _w(b, wb)
    fa, fb = np.isfinite(a), np.isfinite(b)
    a, wa, b, wb = a[fa], wa[fa], b[fb], wb[fb]
    if not len(a) or not len(b):
        return float("nan")
    x = np.sort(np.r_[a, b])
    oa, ob = np.argsort(a), np.argsort(b)
    ca = np.r_[0, np.cumsum(wa[oa])] / wa.sum()
    cb = np.r_[0, np.cumsum(wb[ob])] / wb.sum()
    Fa = ca[np.searchsorted(a[oa], x, "right")]
    Fb = cb[np.searchsorted(b[ob], x, "right")]
    return float(np.max(np.abs(Fa - Fb)))


def compare_features(ref, ctrl, feats, w_ref=None, w_ctrl=None):
    """Tabela χ²/ndf i KS per cecha (ref = prawda, np. tło truth≥8 w SR)."""
    out = {}
    for f in feats:
        a, b = ref[f].to_numpy(dtype=float), ctrl[f].to_numpy(dtype=float)
        out[f] = {"chi2ndf": chi2_ndf(a, b, w_ref, w_ctrl), "ks": ks(a, b, w_ref, w_ctrl)}
    return out


def c2st_auc(A, B, feats, n_max=200000, seed=0, folds=3, engine="hgb"):
    """AA: JEDEN test dwóch próbek → charmbkg.bkggen.c2st (HistGB, wagi). engine="xgb" = stara implementacja S."""
    if engine == "xgb":
        return _c2st_auc_xgb(A, B, feats, n_max=n_max, seed=seed, folds=folds)
    from .c2st import c2st
    return c2st(A, B, feats, n_folds=folds, max_n=n_max, seed=seed)


def _c2st_auc_xgb(A, B, feats, n_max=200000, seed=0, folds=3):
    """AUC klasyfikatora A vs B (XGBoost, k-fold OOF). 0.5 = nierozróżnialne. (oryginał S, zachowany)"""
    import xgboost as xgb
    rng = np.random.default_rng(seed)
    XA, XB = A[feats].to_numpy(np.float32), B[feats].to_numpy(np.float32)
    if len(XA) > n_max:
        XA = XA[rng.choice(len(XA), n_max, replace=False)]
    if len(XB) > n_max:
        XB = XB[rng.choice(len(XB), n_max, replace=False)]
    n = min(len(XA), len(XB))
    XA, XB = XA[:n], XB[:n]
    X = np.r_[XA, XB]
    y = np.r_[np.ones(n), np.zeros(n)]
    f = rng.integers(0, folds, len(y))
    s = np.zeros(len(y))
    for k in range(folds):
        m = xgb.XGBClassifier(n_estimators=200, max_depth=4, learning_rate=0.1, subsample=0.8, n_jobs=2,
                              tree_method="hist")
        m.fit(X[f != k], y[f != k])
        s[f == k] = m.predict_proba(X[f == k])[:, 1]
    return auc(s, y)


def auc(s, y, w=None):
    s, y = np.asarray(s), np.asarray(y)
    w = np.ones(len(s)) if w is None else np.asarray(w)
    o = np.argsort(s, kind="stable")
    s, y, w = s[o], y[o], w[o]
    ws, wb = np.where(y == 1, w, 0), np.where(y == 0, w, 0)
    cb = np.cumsum(wb) - wb / 2
    return float(np.sum(ws * cb) / max(ws.sum() * wb.sum(), 1e-12))
