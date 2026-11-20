"""SIG_constr — narzędzia ML: XGB z wczesnym stopem, klasyfikator z ujemnymi wagami (własna funkcja celu), progi przy zadanym ε_B."""
import numpy as np
import xgboost as xgb

NTH = 1
PAR = dict(max_bin=64, max_depth=4, eta=0.25, subsample=0.8, colsample_bytree=0.8, min_child_weight=5.0, reg_lambda=2.0, tree_method="hist", nthread=NTH, verbosity=0)


def fit_xgb(X, y, w, seed=0, rounds=150, par=None, val_frac=0.25, base_margin=None, es=10):
    """Klasyfikator binarny; wagi klas zrównoważone (Σw_1 = Σw_0); wczesny stop na losowej części treningu. Wagi ≥ 0."""
    p = dict(PAR, objective="binary:logistic", eval_metric="logloss", seed=seed); p.update(par or {})
    w = np.asarray(w, float).copy()
    s1, s0 = w[y == 1].sum(), w[y == 0].sum()
    w[y == 1] *= 0.5 * (s1 + s0) / max(s1, 1e-9); w[y == 0] *= 0.5 * (s1 + s0) / max(s0, 1e-9)
    w *= len(w) / w.sum()
    rng = np.random.default_rng(seed + 17)
    v = rng.random(len(y)) < val_frac
    dt = xgb.DMatrix(X[~v], y[~v], weight=w[~v], nthread=NTH); dv = xgb.DMatrix(X[v], y[v], weight=w[v], nthread=NTH)
    if base_margin is not None:
        dt.set_base_margin(base_margin[~v]); dv.set_base_margin(base_margin[v])
    return xgb.train(p, dt, rounds, evals=[(dv, "v")], early_stopping_rounds=es, verbose_eval=False)


def pred(bst, X, base_margin=None, margin=False):
    d = xgb.DMatrix(X, nthread=NTH)
    if base_margin is not None:
        d.set_base_margin(base_margin)
    return bst.predict(d, iteration_range=(0, bst.best_iteration + 1) if hasattr(bst, "best_iteration") else (0, 0), output_margin=margin)


def fit_ratio(X1, w1, X0, w0, seed=0, rounds=80, depth=3, eta=0.2, mcw=10.0, val_seed=None):
    """Stosunek gęstości p1/p0 klasyfikatorem z wagami dowolnego znaku (w1 może zawierać ujemne wagi odejmowania pasm).
    Własna funkcja celu: grad = w(p − y), hess = w·p(1 − p); liście z Σhess < mcw są zabronione (naturalny regularyzator ujemnych wag).
    Klasy normowane do równych sum wag → wynik exp(margin) = p1/p0 (unormowane gęstości). Wczesny stop: ważona logloss na 25 %."""
    X = np.vstack([X1, X0]); y = np.r_[np.ones(len(X1)), np.zeros(len(X0))]
    w = np.r_[np.asarray(w1, float) / np.sum(w1), np.asarray(w0, float) / np.sum(w0)] * 0.5 * (len(X1) + len(X0))
    rng = np.random.default_rng((seed if val_seed is None else val_seed) + 5)
    v = rng.random(len(y)) < 0.25
    # normalizacja klas osobno w części treningowej i walidacyjnej
    def norm(m):
        ww = w[m].copy(); yy = y[m]
        ww[yy == 1] *= 0.5 * m.sum() / ww[yy == 1].sum(); ww[yy == 0] *= 0.5 * m.sum() / ww[yy == 0].sum()
        return ww
    wt, wv = norm(~v), norm(v)
    dt = xgb.DMatrix(X[~v], y[~v], nthread=NTH); dv = xgb.DMatrix(X[v], y[v], nthread=NTH)

    def obj(pm, d):
        p = 1 / (1 + np.exp(-pm))
        return wt * (p - y[~v]), wt * np.maximum(p * (1 - p), 1e-6)
    par = dict(max_depth=depth, eta=eta, subsample=0.8, colsample_bytree=0.9, min_child_weight=mcw, reg_lambda=5.0, max_bin=64, tree_method="hist", nthread=NTH, verbosity=0,
               seed=seed, base_score=0.5)
    best, best_it, bst, bad = np.inf, 0, None, 0
    bst = xgb.Booster(par, [dt])
    yv = y[v]
    for it in range(rounds):
        bst.update(dt, it, fobj=obj)
        if it % 5 == 4:
            p = np.clip(1 / (1 + np.exp(-bst.predict(dv, output_margin=True))), 1e-6, 1 - 1e-6)
            ll = -np.sum(wv * (yv * np.log(p) + (1 - yv) * np.log(1 - p))) / wv.sum()
            if ll < best - 1e-5:
                best, best_it, bad = ll, it + 1, 0
            else:
                bad += 1
                if bad >= 4:
                    break
    return bst, max(best_it, 5)


def ratio_pred(model, X, clip=(0.05, 20.0)):
    bst, nit = model
    m = bst.predict(xgb.DMatrix(X, nthread=NTH), output_margin=True, iteration_range=(0, nit))
    return np.clip(np.exp(m), *clip)


def auc_w(s1, s0, w1=None, w0=None):
    """AUC z wagami (dowolnego znaku po stronie sygnału)."""
    w1 = np.ones(len(s1)) if w1 is None else w1; w0 = np.ones(len(s0)) if w0 is None else w0
    s = np.r_[s1, s0]; o = np.argsort(s, kind="mergesort")
    y = np.r_[np.ones(len(s1)), np.zeros(len(s0))][o]; w = np.r_[w1, w0][o]
    c0 = np.cumsum(w * (1 - y))
    return float(np.sum(w * y * c0) / (w1.sum() * w0.sum()))


def thr_at_epsb(sb, epsb, wb=None):
    """Próg t: ułamek tła ze score > t równy epsb."""
    if wb is None:
        return float(np.quantile(sb, 1 - epsb))
    o = np.argsort(sb); c = np.cumsum(wb[o]) / wb.sum()
    return float(sb[o][np.searchsorted(c, 1 - epsb)])
