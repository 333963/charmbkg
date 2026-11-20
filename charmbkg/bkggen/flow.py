"""charmbkg.bkggen.flow — model generatywny tła warunkowany masą kandydata D (etap S2).

p(x | m) uczone na sidebandach masy (|m − m0| > sb_min), próbkowanie w oknie sygnału.
Modele:
  * FlowBkg    — warunkowy neural spline flow (zuko.flows.NSF), jak CATHODE (2109.00546).
  * CopulaBkg  — baseline: kopuła Gaussa; z = Φ⁻¹(F(x)), z | m ~ N(a + b·m, Σ) (liniowe w m).
Wspólne: Preproc (kwantyle → N(0,1), dekwantyzacja cech całkowitych, imputacja NaN dE/dx),
MassModel (kształt m z sidebandów: exp·pol1 lub pol2, próbkowanie w SR), walidacja (KS, χ², AUC),
CATHODE (klasyfikator dane-w-oknie vs syntetyczne).
Izolacja (ptcone*) zakazana jako cecha — sprawdzane.
"""
from __future__ import annotations

import json
import logging

import numpy as np

from .. import constants as _C

log = logging.getLogger("charmbkg.bkggen.flow")

FEATS_H = ["log_pt", "abs_eta", "chi2", "Lxy", "impact_sig", "abs_impact_z0", "d_max_abs_d0",
           "pointing", "abs_cos_theta_star", "abs_impact_sig", "d_min_abs_d0", "d_max_abs_z0", "d_max_abs_eta",
           "dedx_K", "log_p_K", "npts_K", "dedx_pi_max", "log_p_pi_min"]
# Zmienne GENEROWANE (bez ścisłych więzów): kinematyka 2-ciał przez (pT, η, cosθ*, φ*) przy danym m,
# d0 córek przez (max, min/max); |IPsig| = IPsig (IPsig ≥ 0 w v1.1).
FEATS_GEN = ["log_pt", "eta", "cts", "phis", "chi2", "Lxy", "impact_sig", "abs_impact_z0", "d_max_abs_d0", "r_d0",
             "pointing", "d_max_abs_z0", "dedx_K", "npts_K", "dedx_pi_max", "npts_pi"]
NAN_IF_ZERO = {"dedx_K": "npts_K", "dedx_pi_max": "npts_pi"}         # dE/dx brak ⇔ brak punktów (≈0,6% kandydatów)


def derive_H(G, m, gen=FEATS_GEN):
    """Macierz G (kolumny FEATS_GEN) + m → macierz cech H (kolumny FEATS_H), więzy ściśle spełnione."""
    from . import kin
    c = {f: G[:, j] for j, f in enumerate(gen)}
    pt = np.expm1(c["log_pt"])
    pK, eK, pP, eP = kin.daughters(np.asarray(m, float), pt, c["eta"], np.clip(c["cts"], -1, 1), c["phis"])
    h = dict(log_pt=c["log_pt"], abs_eta=np.abs(c["eta"]), chi2=c["chi2"], Lxy=c["Lxy"], impact_sig=c["impact_sig"],
             abs_impact_z0=c["abs_impact_z0"], d_max_abs_d0=c["d_max_abs_d0"], pointing=c["pointing"],
             abs_cos_theta_star=np.abs(np.clip(c["cts"], -1, 1)), abs_impact_sig=np.abs(c["impact_sig"]),
             d_min_abs_d0=c["d_max_abs_d0"] * np.clip(c["r_d0"], 0, 1), d_max_abs_z0=c["d_max_abs_z0"],
             d_max_abs_eta=np.maximum(np.abs(eK), np.abs(eP)), dedx_K=c["dedx_K"], log_p_K=np.log(pK),
             npts_K=c["npts_K"], dedx_pi_max=c["dedx_pi_max"], log_p_pi_min=np.log(pP))
    return np.column_stack([h[f] for f in FEATS_H])


FORBIDDEN = ("ptcone", "iso")

SPECIES = {"D0": dict(m0=_C.M_D0, mass_range=(1.75, 2.00), sb_min=0.075, sr_half=0.055),
           "Dplus": dict(m0=_C.M_DPLUS, mass_range=(1.70, 2.20), sb_min=0.070, sr_half=0.045)}


def check_features(feats):
    bad = [f for f in feats if any(b in f for b in FORBIDDEN)]
    if bad:
        raise ValueError(f"cechy zakazane (izolacja): {bad}")
    return list(feats)


def regions(m, sp):
    dm = np.abs(np.asarray(m) - sp["m0"])
    lo, hi = sp["mass_range"]
    inr = (m >= lo) & (m < hi)
    return inr & (dm > sp["sb_min"]), inr & (dm < sp["sr_half"])


# --------------------------------------------------------------------------- preprocessing
class Preproc:
    """x → z: losowany PIT względem posortowanej próbki odniesienia (remisy/spajki, np. npts, pointing=7,
    dostają losowe u w swoim przedziale → ciągłe z), imputacja NaN z rozkładu brzegowego, Φ⁻¹.
    Odwrotność: interpolacja kwantyli — wartości z remisów odtwarzane dokładnie, brak ekstrapolacji poza zakres."""

    def __init__(self, feats, n_ref=200_000, seed=1):
        self.feats = check_features(feats)
        self.n_ref = n_ref
        self.rng = np.random.default_rng(seed)

    def _impute(self, X):
        X = np.array(X, dtype=np.float64, copy=True)
        for j in range(X.shape[1]):
            bad = ~np.isfinite(X[:, j])
            if bad.any():
                pool = X[~bad, j]
                X[bad, j] = self.rng.choice(pool, bad.sum()) if len(pool) else 0.0
        return X

    def fit(self, X):
        X = self._impute(X)
        i = self.rng.choice(len(X), min(len(X), self.n_ref), replace=False)
        self.ref = np.sort(X[i], axis=0)
        n = len(self.ref)
        self.uc = (np.arange(n) + 0.5) / n
        return self

    def forward(self, X):
        from scipy.special import ndtri
        X = self._impute(X)
        n = len(self.ref)
        Z = np.empty_like(X)
        for j in range(X.shape[1]):
            r = self.ref[:, j]
            lo = np.searchsorted(r, X[:, j], "left") / n
            hi = np.searchsorted(r, X[:, j], "right") / n
            u = np.interp(X[:, j], r, self.uc)
            tie = (hi - lo) > 1.5 / n
            u[tie] = lo[tie] + self.rng.random(tie.sum()) * (hi[tie] - lo[tie])
            Z[:, j] = ndtri(np.clip(u, 0.25 / n, 1 - 0.25 / n))
        return Z

    def inverse(self, Z):
        from scipy.special import ndtr
        X = np.empty_like(Z)
        for j in range(Z.shape[1]):
            X[:, j] = np.interp(ndtr(Z[:, j]), self.uc, self.ref[:, j])
        for f, g in NAN_IF_ZERO.items():
            if f in self.feats and g in self.feats:
                X[X[:, self.feats.index(g)] <= 0, self.feats.index(f)] = np.nan
        return X


def m_context(m, sp):
    return ((np.asarray(m, dtype=np.float64) - sp["m0"]) / 0.1)[:, None]


# --------------------------------------------------------------------------- mass model
class MassModel:
    """Gęstość tła w m z sidebandów: exp(a·t + b·t²) (t = (m−m0)/0.1), ML bez binowania (scipy)."""

    def __init__(self, sp):
        self.sp = sp

    def fit(self, m):
        from scipy.optimize import minimize
        sb, _ = regions(m, self.sp)
        t = (m[sb] - self.sp["m0"]) / 0.1
        grid = np.linspace(*[(x - self.sp["m0"]) / 0.1 for x in self.sp["mass_range"]], 4001)
        sbg = np.abs(grid) > self.sp["sb_min"] / 0.1

        def nll(p):
            lg = p[0] * grid + p[1] * grid ** 2
            norm = np.trapezoid(np.exp(lg) * sbg, grid)
            return -(np.sum(p[0] * t + p[1] * t ** 2) - len(t) * np.log(norm))

        self.p = minimize(nll, [0.0, 0.0], method="Nelder-Mead").x
        self.grid = grid
        self.n_sb = int(sb.sum())
        dens = np.exp(self.p[0] * grid + self.p[1] * grid ** 2)
        self.frac_sr_over_sb = float(np.trapezoid(dens * (np.abs(grid) < self.sp["sr_half"] / 0.1), grid)
                                     / np.trapezoid(dens * sbg, grid))
        return self

    def expected_sr(self):
        return self.n_sb * self.frac_sr_over_sb

    def sample(self, n, region="SR", rng=None):
        rng = rng or np.random.default_rng()
        g = self.grid
        dens = np.exp(self.p[0] * g + self.p[1] * g ** 2)
        if region == "SR":
            dens = dens * (np.abs(g) < self.sp["sr_half"] / 0.1)
        cdf = np.cumsum(dens)
        cdf /= cdf[-1]
        return self.sp["m0"] + 0.1 * np.interp(rng.random(n), cdf, g)


# --------------------------------------------------------------------------- models
class CopulaBkg:
    """Baseline: z | m ~ N(A·[1, c, c²] , Σ) w przestrzeni kwantyli Gaussa (kopuła Gaussa, średnie kwadratowe w m)."""

    name = "copula"

    def __init__(self, feats, sp, deg=1, seed=1):
        self.pre = Preproc(feats, seed=seed)
        self.sp, self.deg = sp, deg
        self.rng = np.random.default_rng(seed)

    def _B(self, m):
        c = m_context(m, self.sp)[:, 0]
        return np.column_stack([c ** k for k in range(self.deg + 1)])

    def fit(self, X, m, **kw):
        Z = self.pre.fit(X).forward(X)
        B = self._B(m)
        self.A, *_ = np.linalg.lstsq(B, Z, rcond=None)
        R = Z - B @ self.A
        self.L = np.linalg.cholesky(np.cov(R.T) + 1e-9 * np.eye(R.shape[1]))
        return self

    def sample_z(self, m):
        return self._B(m) @ self.A + self.rng.standard_normal((len(m), self.L.shape[0])) @ self.L.T

    def sample(self, m):
        return self.pre.inverse(self.sample_z(np.asarray(m)))


class FlowBkg:
    """Warunkowy NSF (zuko) na z = Preproc(x), kontekst c = (m − m0)/0.1. Ensemble ostatnich/najlepszych epok."""

    name = "flow"

    def __init__(self, feats, sp, transforms=6, hidden=(256, 256), bins=8, seed=1, device=None):
        import torch
        import zuko
        self.torch = torch
        torch.manual_seed(seed)
        self.pre = Preproc(feats, seed=seed)
        self.sp = sp
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.flow = zuko.flows.NSF(features=len(feats), context=1, transforms=transforms,
                                   hidden_features=list(hidden), bins=bins).to(self.device)
        self.history = []

    def fit(self, X, m, epochs=30, batch=1024, lr=1e-3, val_frac=0.2, patience=6, max_minutes=None, verbose=True):
        import time
        torch = self.torch
        Z = self.pre.fit(X).forward(X).astype(np.float32)
        C = m_context(m, self.sp).astype(np.float32)
        rng = np.random.default_rng(0)
        idx = rng.permutation(len(Z))
        nv = int(val_frac * len(Z))
        tv = lambda a, i: torch.from_numpy(a[i]).to(self.device)  # noqa: E731
        Zt, Ct, Zv, Cv = tv(Z, idx[nv:]), tv(C, idx[nv:]), tv(Z, idx[:nv]), tv(C, idx[:nv])
        opt = torch.optim.AdamW(self.flow.parameters(), lr=lr, weight_decay=1e-5)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
        best, best_state, bad, t0 = np.inf, None, 0, time.time()
        for ep in range(epochs):
            self.flow.train()
            perm = torch.randperm(len(Zt), device=self.device)
            tot = 0.0
            for i in range(0, len(Zt), batch):
                b = perm[i:i + batch]
                loss = -self.flow(Ct[b]).log_prob(Zt[b]).mean()
                opt.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.flow.parameters(), 5.0)
                opt.step()
                tot += float(loss.detach()) * len(b)
            sched.step()
            self.flow.eval()
            with torch.no_grad():
                vl = float(np.mean([-self.flow(Cv[i:i + 8192]).log_prob(Zv[i:i + 8192]).mean().item()
                                    for i in range(0, len(Zv), 8192)]))
            self.history.append((ep, tot / len(Zt), vl, time.time() - t0))
            if verbose:
                log.info("epoka %d: train %.4f val %.4f (%.0f s)", ep, tot / len(Zt), vl, time.time() - t0)
            if vl < best - 1e-4:
                best, bad = vl, 0
                best_state = {k: v.detach().clone() for k, v in self.flow.state_dict().items()}
            else:
                bad += 1
                if bad >= patience:
                    break
            if max_minutes and (time.time() - t0) > 60 * max_minutes:
                break
        if best_state is not None:
            self.flow.load_state_dict(best_state)
        self.best_val = best
        return self

    def sample_z(self, m, chunk=4_000):
        torch = self.torch
        C = torch.from_numpy(m_context(m, self.sp).astype(np.float32)).to(self.device)
        out = []
        self.flow.eval()
        with torch.no_grad():
            for i in range(0, len(C), chunk):
                c = C[i:i + chunk]
                out.append(_sample_ctx(self.flow, c).cpu().numpy())
        Z = np.concatenate(out).astype(np.float64)
        return np.clip(Z, -6, 6)

    def sample(self, m):
        return self.pre.inverse(self.sample_z(np.asarray(m)))

    def save(self, path):
        torch = self.torch
        torch.save({"state": self.flow.state_dict(), "ref": self.pre.ref, "feats": self.pre.feats,
                    "sp": self.sp, "history": self.history}, path)


def load_flow(path, transforms=6, hidden=(256, 256), bins=8, device=None):
    import torch
    ck = torch.load(path, map_location="cpu", weights_only=False)
    fb = FlowBkg(ck["feats"], ck["sp"], transforms=transforms, hidden=hidden, bins=bins, device=device)
    fb.flow.load_state_dict(ck["state"])
    fb.pre.ref = ck["ref"]
    n = len(fb.pre.ref)
    fb.pre.uc = (np.arange(n) + 0.5) / n
    fb.history = ck.get("history", [])
    return fb


def _sample_ctx(flow, c):
    """Jedna próbka na każdy wiersz kontekstu c (N,1)."""
    d = flow(c)
    return d.sample()


# --------------------------------------------------------------------------- walidacja
def ks_chi2(Xa, Xb, feats, nb=40):
    from scipy.stats import ks_2samp
    res = {}
    for j, f in enumerate(feats):
        a, b = Xa[:, j], Xb[:, j]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        ks = ks_2samp(a, b)
        edges = np.unique(np.quantile(a, np.linspace(0, 1, nb + 1)))
        ha, _ = np.histogram(a, edges)
        hb, _ = np.histogram(b, edges)
        hb = hb * (ha.sum() / max(hb.sum(), 1))
        ok = (ha + hb) > 0
        chi2 = float(np.sum((ha[ok] - hb[ok]) ** 2 / (ha[ok] + hb[ok] * (len(a) / max(len(b), 1)))))
        res[f] = dict(ks=float(ks.statistic), ks_p=float(ks.pvalue), chi2_ndf=chi2 / max(ok.sum() - 1, 1))
    return res


def auc_two_samples(Xa, Xb, n_folds=3, max_n=200_000, seed=0, return_scores=False):
    """AUC klasyfikatora HistGB a vs b (out-of-fold, zbalansowane); cel ≈ 0,5.
    AA: bez return_scores deleguje do charmbkg.bkggen.c2st (jedna implementacja)."""
    if not return_scores:
        from .c2st import c2st
        return c2st(np.asarray(Xa), np.asarray(Xb), n_folds=n_folds, max_n=max_n, seed=seed)
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    rng = np.random.default_rng(seed)
    n = min(len(Xa), len(Xb), max_n)
    A = Xa[rng.choice(len(Xa), n, replace=False)]
    B = Xb[rng.choice(len(Xb), n, replace=False)]
    X = np.vstack([A, B])
    y = np.r_[np.ones(n), np.zeros(n)]
    s = np.zeros(len(y))
    for tr, te in StratifiedKFold(n_folds, shuffle=True, random_state=seed).split(X, y):
        clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, max_leaf_nodes=31,
                                             early_stopping=True, random_state=seed)
        clf.fit(X[tr], y[tr])
        s[te] = clf.predict_proba(X[te])[:, 1]
    auc = float(roc_auc_score(y, s))
    return (auc, s, y) if return_scores else auc


# --------------------------------------------------------------------------- CATHODE
def cathode_scores(X_sr, X_syn, n_folds=5, seed=0, clf_params=None):
    """Klasyfikator SR (y=1) vs syntetyczne (y=0); zwraca wynik OOF dla każdego kandydata SR.
    Dla dużej liczby syntetycznych waga klasy 0 = N_SR/N_syn (stosunek prawdopodobieństw ∝ p_data/p_bkg)."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.model_selection import KFold
    p = dict(max_iter=400, learning_rate=0.05, max_leaf_nodes=31, early_stopping=True, l2_regularization=1.0,
             random_state=seed)
    p.update(clf_params or {})
    rng = np.random.default_rng(seed)
    s_sr = np.zeros(len(X_sr))
    syn_fold = rng.integers(0, n_folds, len(X_syn))
    for k, (tr, te) in enumerate(KFold(n_folds, shuffle=True, random_state=seed).split(X_sr)):
        Xs = X_syn[syn_fold != k]
        X = np.vstack([X_sr[tr], Xs])
        y = np.r_[np.ones(len(tr)), np.zeros(len(Xs))]
        w = np.r_[np.ones(len(tr)), np.full(len(Xs), len(tr) / len(Xs))]
        clf = HistGradientBoostingClassifier(**p)
        clf.fit(X, y, sample_weight=w)
        s_sr[te] = clf.predict_proba(X_sr[te])[:, 1]
    return s_sr


def save_json(obj, path):
    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=conv)
