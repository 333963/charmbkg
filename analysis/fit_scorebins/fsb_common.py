"""FIT_scorebins — silnik: binowany fit Poissona masy w kategoriach score (Fisher scoring + LM).

Model w kategorii c, binie masy i:  mu_ci = a_c(theta_s) * U_ci + exp(X_i . beta_c)
  U_c  — szablon sygnału (+ refleksja) na jednostkę prawdziwego sygnału w kategorii (suma części sygnałowej = 1)
  a_c  — amplituda sygnału w kategorii:
     'fixed'  a_c = N * f_c                      (ułamki z MC)
     'free'   a_c = N_c                          (N = suma)
     'constr' a_c = N * phi_c(delta), phi_c = E_c(1+d_c) - E_{c+1}(1+d_{c+1}), d_0 = 0, d_k ~ Gauss(0, sig)
              (E_k = skumulowana wydajność MC powyżej dolnej granicy kategorii k — niepewność typu SF na każdej granicy)
Kategorie uporządkowane od najniższego score (c = 0) do najwyższego.
"""
import numpy as np

EFF_EDGES = np.array([1.0, .95, .90, .85, .80, .75, .70, .65, .60, .55, .50, .45, .40, .35, .30, .20, .10, 0.0])
NCELL = len(EFF_EDGES) - 1          # komórka j: eps_S skumulowana w (EFF_EDGES[j+1], EFF_EDGES[j]]; j = 0 najniższy score
WPS = [.95, .90, .85, .80, .75, .70, .65, .60, .55, .50, .45, .40, .35, .30]
PTB = [(3, 5), (5, 8), (8, 30)]


def xbasis(nm):
    x = np.linspace(-1, 1, nm)
    return np.stack([np.ones(nm), x, x * x], 1)


def merge(arr, bounds):
    """arr [NCELL, M] -> [C, M]; bounds = lista początków kategorii (indeksy komórek), np. [0, 5, 12]; koniec = NCELL
    lub bounds[-1] gdy podano krotkę (starts, end)."""
    starts, end = bounds if isinstance(bounds, tuple) else (bounds, arr.shape[0])
    ed = list(starts) + [end]
    return np.stack([arr[ed[i]:ed[i + 1]].sum(0) for i in range(len(starts))])


def init_beta(b, X, mask=None):
    """start: regresja log(b) na bazie X (b > 0)."""
    w = np.ones(len(b)) if mask is None else mask.astype(float)
    y = np.log(np.maximum(b, 1e-3))
    return np.linalg.lstsq(X * w[:, None], y * w, rcond=None)[0]


class Model:
    def __init__(self, U, X, mode, f=None, sig=None, mask=None):
        self.U, self.X, self.mode = U, X, mode
        self.C, self.M = U.shape
        self.f = None if f is None else np.asarray(f, float)
        self.sig = sig
        self.mask = np.ones((self.C, self.M)) if mask is None else np.broadcast_to(mask, (self.C, self.M)).astype(float)
        self.ns = {"fixed": 1, "free": self.C, "constr": self.C, "none": 0}[mode]
        if mode == "constr":
            self.E = np.cumsum(self.f[::-1])[::-1]          # E_c = suma f od c w górę; E_0 = pokrycie
        self.np_ = self.ns + 3 * self.C

    def amps(self, ts):
        C = self.C
        if self.mode == "none":
            return np.zeros(C), np.zeros((C, 0))
        if self.mode == "fixed":
            return ts[0] * self.f, self.f[:, None].copy()
        if self.mode == "free":
            return ts.copy(), np.eye(C)
        N, d = ts[0], np.r_[0.0, ts[1:]]
        Ee = self.E * (1 + d)
        phi = Ee - np.r_[Ee[1:], 0.0]
        A = np.zeros((C, C))
        A[:, 0] = phi
        for k in range(1, C):
            A[k, k] += N * self.E[k]
            A[k - 1, k] -= N * self.E[k]
        return N * phi, A

    def mu(self, th):
        a, A = self.amps(th[:self.ns])
        b = np.exp(np.clip(self.X @ th[self.ns:].reshape(self.C, 3).T, -50, 50)).T      # [C, M]
        return a[:, None] * self.U + b, b, A

    def nll(self, th, n):
        mu, _, _ = self.mu(th)
        if np.any(mu[self.mask > 0] <= 0):
            return np.inf
        mu = np.maximum(mu, 1e-300)
        v = (self.mask * (mu - n * np.log(mu))).sum()
        if self.mode == "constr":
            v += 0.5 * (th[1:self.ns] ** 2).sum() / self.sig ** 2
        return v

    def grad_info(self, th, n):
        mu, b, A = self.mu(th)
        ns, C = self.ns, self.C
        r = self.mask * (1 - n / mu)
        w = self.mask / mu
        g = np.zeros(self.np_)
        H = np.zeros((self.np_, self.np_))
        X = self.X
        for c in range(C):
            J = np.empty((self.M, ns + 3))
            J[:, :ns] = self.U[c][:, None] * A[c][None, :]
            J[:, ns:] = b[c][:, None] * X
            idx = np.r_[np.arange(ns), ns + 3 * c + np.arange(3)]
            g[idx] += J.T @ r[c]
            H[np.ix_(idx, idx)] += (J * w[c][:, None]).T @ J
        if self.mode == "constr":
            g[1:ns] += th[1:ns] / self.sig ** 2
            H[np.arange(1, ns), np.arange(1, ns)] += 1 / self.sig ** 2
        return g, H

    def fit(self, n, th0, maxit=200, tol=1e-7):
        th = np.array(th0, float)
        f0 = self.nll(th, n)
        lam = 1e-3
        ok = False
        small = 0
        for _ in range(maxit):
            g, H = self.grad_info(th, n)
            dH = np.diag(H).copy()
            while True:
                try:
                    step = -np.linalg.solve(H + lam * np.diag(dH + 1e-12), g)
                except np.linalg.LinAlgError:
                    step = None
                f1 = self.nll(th + step, n) if step is not None else np.inf
                if f1 <= f0 + 1e-12:
                    break
                lam *= 10
                if lam > 1e10:
                    break
            if lam > 1e10:
                break
            th = th + step
            small = small + 1 if (f0 - f1) < tol else 0
            f0 = f1
            lam = max(lam / 10, 1e-9)
            if small >= 2:
                ok = True
                break
        g, H = self.grad_info(th, n)
        try:
            cov = np.linalg.inv(H)
        except np.linalg.LinAlgError:
            return th, np.nan, np.nan, False, f0
        if self.mode == "free":
            N, var = th[:self.ns].sum(), cov[:self.ns, :self.ns].sum()
        elif self.mode == "none":
            N, var = 0.0, 0.0
        else:
            N, var = th[0], cov[0, 0]
        return th, N, np.sqrt(max(var, 0)), ok, f0

    def theta0(self, N, bkg):
        beta = np.concatenate([init_beta(bkg[c], self.X) for c in range(self.C)])
        if self.mode == "fixed":
            s = [N]
        elif self.mode == "free":
            s = list(N * self.f)
        elif self.mode == "constr":
            s = [N] + [0.0] * (self.C - 1)
        else:
            s = []
        return np.r_[s, beta]


def build(inp, bounds, mode, sig=None, f_fit=None):
    """inp: dict z prep (f[NCELL], Us[NCELL, M] sygnał na jednostkę sygnału w komórce, Ur — refleksja na jednostkę sygnału).
    Zwraca Model oraz ułamki f_c; szablon kategorii = średnia ważona f komórek (jak histogram MC regionu)."""
    f = inp["f"] if f_fit is None else f_fit
    T = f[:, None] * (inp["Us"] + inp["Ur"])                # oczekiwany sygnał+refl. na jednostkę N_S
    Tm = merge(T, bounds)
    fc = merge(f[:, None], bounds)[:, 0]
    U = Tm / fc[:, None]
    return Model(U, inp["X"], mode, f=fc, sig=sig), fc
