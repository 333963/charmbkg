"""R3 — silnik: jednoczesny fit Δm w kategoriach score (D0<-D*), model jak dstar_refit.fit_passfail (próg × pol2, b = 0).
Kategoria c: lam_c = S_c * G2(dm; mu, s*r_c*sc, k, f) + B_c * u^a_c (1 + c1_c t + c2_c t^2), u = dm - m_pi, t = u/20.
Parametry: mu, s, potem per kategoria (S, B, a, c1, c2). NLL Poissona (postać nasycona), gradient analityczny.
Pass/fail = to samo z 2 kategoriami. Projekt B; nic nie zapisuje poza katalogiem R3.
"""
from __future__ import annotations

import json
import os

import numpy as np
from iminuit import Minuit
from scipy.special import ndtr

HERE = os.path.dirname(os.path.abspath(__file__))
BT = os.path.expanduser("~/Documents/charm/charm_work/B_trig/runs/lxplus")
SFJSON = os.path.join(BT, "out", "btrig_sf.json")
M0_D0, NOM_WIN, MPI = 1.86484, 0.055, 139.57
BINS = {"3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0)}
RUNS = (428770, 427929)
E = np.round(np.arange(140.5, 160.0001, 0.1), 6)          # MeV; 195 binów (jak _rebin_pf na E_DM)
NB = len(E) - 1
CEN = 0.5 * (E[1:] + E[:-1])
WP_MAIN = (0.95, 0.80, 0.70, 0.50, 0.25)                    # dolne granice kategorii 1..5 (skumulowana eps MC)
WP_EQ = (0.95, 0.76, 0.57, 0.38, 0.19)                      # wariant: równe ułamki sygnału
SIDE_GEN = 3.5                                              # MeV, wyłączenie piku w fitach tła do generatora
SIDE_START = 2.5

_GLX, _GLW = np.polynomial.legendre.leggauss(6)
_X = CEN[:, None] + 0.05 * _GLX[None, :]
_U = _X - MPI
_LU = np.log(_U)
_T = _U / 20.0
_W = 0.05 * _GLW[None, :]
_SQ2PI = np.sqrt(2 * np.pi)


def bkg_shape(a, c1, c2, grad=False):
    """Znormalizowany kształt tła w binach (suma = 1) [+ pochodne po a, c1, c2]."""
    ua = np.exp(a * _LU)
    P = 1 + c1 * _T + c2 * _T * _T
    ok = P > 1e-12
    Pc = np.where(ok, P, 1e-12)
    I = (ua * Pc * _W).sum(1)
    Tt = I.sum()
    B = I / Tt
    if not grad:
        return B
    dIa = (ua * Pc * _LU * _W).sum(1)
    dI1 = (ua * _T * ok * _W).sum(1)
    dI2 = (ua * _T * _T * ok * _W).sum(1)
    return B, [(d - B * d.sum()) / Tt for d in (dIa, dI1, dI2)]


def sig_shape(mu, sig1, k, f, grad=False):
    """Znormalizowany na [E0, E-1] kształt sygnału (2 Gaussy) w binach [+ pochodne po mu i po ln(sig1)]."""
    z1 = (E - mu) / sig1
    z2 = z1 / k
    C = f * ndtr(z1) + (1 - f) * ndtr(z2)
    D = np.diff(C)
    Tt = C[-1] - C[0]
    S = D / Tt
    if not grad:
        return S
    p1 = np.exp(-0.5 * z1 * z1) / _SQ2PI
    p2 = np.exp(-0.5 * z2 * z2) / _SQ2PI
    dCmu = -(f * p1 / sig1 + (1 - f) * p2 / (k * sig1))
    dCls = -(f * p1 * z1 + (1 - f) * p2 * z2)
    out = []
    for dC in (dCmu, dCls):
        dD = np.diff(dC)
        out.append((dD - S * (dC[-1] - dC[0])) / Tt)
    return S, out


NSB_FULL, NSB_LIN = 400, 100


def ncoef_rule(H, mu0):
    """Złożoność tła per kategoria wg liczby wpisów poza pikiem (|dm - mu| > 2,5 MeV): >= 400 -> a, c1, c2; 100-400 -> c2 = 0;
    < 100 -> c1 = c2 = 0 (sam próg u^a). Reguła ustalona na pseudo-danych, przed fitami danych (KRYTERIA §1 uzupełnienie)."""
    nsb = np.asarray(H, float)[:, np.abs(CEN - mu0) > SIDE_START].sum(1)
    return np.where(nsb >= NSB_FULL, 2, np.where(nsb >= NSB_LIN, 1, 0))


class CatFit:
    """Fit jednoczesny. H: (ncat, NB) zliczenia; shape: dict(mu, sc, k, f); r: (ncat,) stosunki szerokości (z MC)."""

    def __init__(self, H, shape, r=None, ncoef=None):
        self.H = np.asarray(H, float)
        self.ncoef = ncoef_rule(self.H, shape["mu"]) if ncoef is None else np.asarray(ncoef, int)
        self.nc = self.H.shape[0]
        self.shape = shape
        self.r = np.ones(self.nc) if r is None else np.asarray(r, float)
        self.pos = self.H > 0
        self.names = ["mu", "s"] + [f"{p}{c}" for c in range(self.nc) for p in ("S", "B", "a", "c1_", "c2_")]

    def lam(self, x):
        mu, s = x[0], x[1]
        L = np.empty_like(self.H)
        for c in range(self.nc):
            S, B, a, c1, c2 = x[2 + 5 * c: 7 + 5 * c]
            L[c] = S * sig_shape(mu, self.shape["sc"] * s * self.r[c], self.shape["k"], self.shape["f"]) + B * bkg_shape(a, c1, c2)
        return L

    def nll(self, x):
        L = np.clip(self.lam(x), 1e-9, None)
        n = self.H
        d = np.where(self.pos, L / np.where(self.pos, n, 1.0) - 1.0, 0.0)
        return float(np.sum(np.where(self.pos, n * (d - np.log1p(d)), L)))

    def grad(self, x):
        mu, s = x[0], x[1]
        g = np.zeros(len(x))
        for c in range(self.nc):
            S, B, a, c1, c2 = x[2 + 5 * c: 7 + 5 * c]
            sg, (dmu, dls) = sig_shape(mu, self.shape["sc"] * s * self.r[c], self.shape["k"], self.shape["f"], grad=True)
            bg, (da, d1, d2) = bkg_shape(a, c1, c2, grad=True)
            L = S * sg + B * bg
            w = np.where(L > 1e-9, 1.0 - self.H[c] / np.clip(L, 1e-9, None), 0.0)
            g[0] += S * (w @ dmu)
            g[1] += S * (w @ dls) / s
            g[2 + 5 * c: 7 + 5 * c] = (w @ sg, w @ bg, B * (w @ da), B * (w @ d1), B * (w @ d2))
        return g

    # ---- starty: tło z pasm bocznych per kategoria, sygnał z nadwyżki w oknie
    def bkg_only(self, c, side=SIDE_START, mu0=None):
        mu0 = self.shape["mu"] if mu0 is None else mu0
        sb = np.abs(CEN - mu0) > side
        h = self.H[c][sb]
        pos = h > 0

        def f(B, a, c1, c2):
            L = np.clip(B * bkg_shape(a, c1, c2)[sb], 1e-9, None)
            d = np.where(pos, L / np.where(pos, h, 1.0) - 1.0, 0.0)
            return float(np.sum(np.where(pos, h * (d - np.log1p(d)), L)))
        best = None
        for a0 in (0.5, 1.5):
            m = Minuit(f, B=max(self.H[c].sum(), 1.0), a=a0, c1=0.0, c2=0.0)
            m.errordef = 0.5
            m.errors["B"] = max(np.sqrt(self.H[c].sum()), 1.0)
            m.limits["B"] = (0, None); m.limits["a"] = (0.01, 3); m.limits["c1"] = (-5, 20); m.limits["c2"] = (-10, 40)
            m.fixed["c2"] = self.ncoef[c] < 2; m.fixed["c1"] = self.ncoef[c] < 1
            m.migrad(ncall=5000)
            if best is None or m.fval < best.fval:
                best = m
        v = [float(best.values[k]) for k in ("B", "a", "c1", "c2")]
        return v

    def start(self):
        x = [self.shape["mu"], 1.0]
        win = np.abs(CEN - self.shape["mu"]) <= SIDE_START
        for c in range(self.nc):
            B, a, c1, c2 = self.bkg_only(c)
            exc = self.H[c][win].sum() - B * bkg_shape(a, c1, c2)[win].sum()
            x += [float(exc), B, a, c1, c2]
        return np.array(x)

    def minuit(self, x0, strategy=1):
        m = Minuit(self.nll, np.asarray(x0, float), grad=self.grad, name=self.names)
        m.errordef = 0.5
        m.strategy = strategy
        m.tol = 0.01
        m.limits["mu"] = (self.shape["mu"] - 0.5, self.shape["mu"] + 0.5)
        m.limits["s"] = (0.8, 1.6)
        m.errors["mu"] = 0.05; m.errors["s"] = 0.05
        for c in range(self.nc):
            n = max(self.H[c].sum(), 1.0)
            m.errors[f"S{c}"] = max(np.sqrt(n) * 0.5, 5.0)
            m.errors[f"B{c}"] = max(np.sqrt(n), 1.0)
            m.limits[f"B{c}"] = (0, None)
            m.limits[f"a{c}"] = (0.01, 3); m.limits[f"c1_{c}"] = (-5, 20); m.limits[f"c2_{c}"] = (-10, 40)
            m.errors[f"a{c}"] = 0.05; m.errors[f"c1_{c}"] = 0.1; m.errors[f"c2_{c}"] = 0.1
            if self.ncoef[c] < 2:
                m.values[f"c2_{c}"] = 0.0; m.fixed[f"c2_{c}"] = True
            if self.ncoef[c] < 1:
                m.values[f"c1_{c}"] = 0.0; m.fixed[f"c1_{c}"] = True
        return m

    def fit(self, x0=None, hesse=True):
        """Zwraca dict: x, S (ncat), covS (ncat×ncat) lub None, valid, posdef, good, nll, at_limit, mu, s."""
        x0 = self.start() if x0 is None else x0
        m = self.minuit(x0, 1)
        m.migrad(ncall=40000)
        if not m.valid:                                   # druga próba: simplex + strategia 2
            m2 = self.minuit(np.array(m.values), 2)
            m2.simplex(ncall=20000); m2.migrad(ncall=60000)
            if m2.valid or m2.fval < m.fval:
                m = m2
        covS, posdef, nfix = None, False, 0
        if hesse:
            # HESSE przelicza EDM z dokładną kowariancją; gdy EDM > próg (tol 0,01 → 1e-5), Migrad kontynuuje z tej kowariancji
            for _ in range(4):
                try:
                    m.hesse()
                except Exception:  # noqa: BLE001
                    break
                if m.valid:
                    break
                m.migrad(ncall=40000)
            posdef = bool(m.fmin.has_posdef_covar) and m.covariance is not None
            if not (posdef and m.valid):
                # parametry tła na granicy (typowo a = 0,01 przy małej statystyce kategorii) psują HESSE: ustalamy je na granicy
                # i powtarzamy migrad + hesse dla pozostałych (kowariancja warunkowa — standard dla parametrów na granicy)
                lim = [p for p in m.parameters if p[0] in "ac" and not m.fixed[p] and any(abs(m.values[p] - b) < 1e-3 * max(1.0, abs(b)) for b in m.limits[p])]
                if lim:
                    for p in lim:
                        m.fixed[p] = True
                    nfix = len(lim)
                    for _ in range(4):
                        m.migrad(ncall=40000)
                        try:
                            m.hesse()
                        except Exception:  # noqa: BLE001
                            break
                        if m.valid:
                            break
                    posdef = bool(m.fmin.has_posdef_covar) and m.covariance is not None
            if m.covariance is not None:
                idx = [2 + 5 * c for c in range(self.nc)]
                covS = np.array(m.covariance)[np.ix_(idx, idx)]
        x = np.array(m.values)
        atl = [p for p in m.parameters if any(np.isfinite(b) and abs(m.values[p] - b) < 1e-3 * max(1.0, abs(b))
                                              for b in m.limits[p] if b is not None) and not p.startswith("B") and not m.fixed[p]]
        S = x[2::5].copy()
        return {"x": x, "S": S, "covS": covS, "valid": bool(m.valid), "posdef": posdef, "good": bool(m.valid and (posdef or not hesse)),
                "nll": float(m.fval), "at_limit": atl, "n_fixed_at_limit": nfix, "mu": float(x[0]), "s": float(x[1]),
                "mu_err": float(m.errors["mu"]), "s_err": float(m.errors["s"])}


# ---------------------------------------------------------------------------------------------- pochodne z wydajków
def cum_eff(S, covS=None):
    """Kategorie uporządkowane rosnąco w score (0 = odrzucone). Zwraca N_tot, f (ncat), E_j = sum_{c>=j} f_c dla j = 1..ncat-1
    oraz (jeśli covS) błąd N_tot, cov(f), cov(E)."""
    S = np.asarray(S, float)
    N = S.sum()
    f = S / N
    nc = len(S)
    A = np.array([[1.0 if c >= j else 0.0 for c in range(nc)] for j in range(1, nc)])     # E = A f
    Ev = A @ f
    if covS is None:
        return N, f, Ev
    Jf = (np.eye(nc) - f[:, None]) / N                    # df_c/dS_k = (delta_ck - f_c)/N
    covf = Jf @ covS @ Jf.T
    covE = A @ covf @ A.T
    return N, f, Ev, float(np.sqrt(covS.sum())), covf, covE


def mc_cov_f(f, n):
    """Kowariancja wielomianowa ułamków MC (n zdarzeń)."""
    f = np.asarray(f, float)
    return (np.diag(f) - np.outer(f, f)) / n


# ---------------------------------------------------------------------------------------------- wejścia
def load_npz(name, keys):
    z = np.load(os.path.join(BT, f"bt_{name}.npz"), allow_pickle=False)
    return {k: z[k] for k in keys}


def official():
    return json.load(open(SFJSON))


def prep_inputs():
    """Zwraca dict: mc[bin] (score, dm_mev), data[run][bin] (score, dm_mev, evt_inv, n_evt), thr[bin][wp-str], shapes[bin]."""
    off = official()
    M = load_npz("fd_nd", ["m", "pt", "score", "dm", "truth", "truth_ds"])
    fs = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - M0_D0) < NOM_WIN)
    out = {"mc": {}, "data": {}, "thr": off["thresholds"], "shapes": off["shapes"]}
    for b, (lo, hi) in BINS.items():
        s = fs & (M["pt"] >= lo) & (M["pt"] < hi)
        out["mc"][b] = {"score": M["score"][s].astype(float), "dm": M["dm"][s].astype(float) * 1e3}
    for run in RUNS:
        D = load_npz(f"data_{run}", ["m", "pt", "score", "dm", "run", "event", "e_grl"])
        base = np.abs(D["m"] - M0_D0) < NOM_WIN
        if run == 427929:
            base &= D["e_grl"] == 1
        out["data"][run] = {}
        for b, (lo, hi) in BINS.items():
            s = base & (D["pt"] >= lo) & (D["pt"] < hi)
            _, inv = np.unique(D["event"][s], return_inverse=True)
            out["data"][run][b] = {"score": D["score"][s].astype(float), "dm": D["dm"][s].astype(float) * 1e3, "inv": inv,
                                   "nev": int(inv.max()) + 1 if len(inv) else 0}
    return out


def edges_for(inp, b, wps):
    """Progi score (rosnąco) dla listy skumulowanych wydajności (malejąco). Oficjalne progi tam, gdzie istnieją; inaczej kwantyl MC par."""
    t = []
    for w in wps:
        key = f"{w:.2f}"
        if key in inp["thr"][b]:
            t.append(inp["thr"][b][key])
        else:
            t.append(float(np.quantile(inp["mc"][b]["score"], 1 - w)))
    return np.array(t)


def cat_index(score, thr):
    """Kategoria 0 = score < thr[0]; c = liczba progów <= score."""
    return np.searchsorted(thr, score, side="right")


def hist_cats(dm, cat, ncat, w=None):
    ib = np.floor((dm - E[0]) / 0.1 + 1e-9).astype(int)
    ok = (ib >= 0) & (ib < NB)
    flat = cat[ok] * NB + ib[ok]
    h = np.bincount(flat, weights=None if w is None else w[ok], minlength=ncat * NB)
    return h.reshape(ncat, NB).astype(float)
