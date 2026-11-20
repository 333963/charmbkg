"""Generator tła warunkowany masą (styl CATHODE / FETA-lite): p(x | m) uczony na sidebandach → próbkowanie w oknie sygnału.

Wspólne przetwarzanie: każda cecha → N(0,1) przez globalny empiryczny CDF (z drganiem dla wartości dyskretnych/obciętych),
model warunkowy w przestrzeni z, powrót przez kwantyle. Dwa modele:
  CopulaGen — z | m ~ N(μ(m), Σ(m)), μ i Σ liniowe w m (dopasowane na binach masy sidebandów) — baseline, sekundy
  FlowGen   — neural spline flow (zuko NSF) z kontekstem m — nieparametryczny (CATHODE)
Masa próbkowana z gęstości tła dopasowanej w sidebandach (exp(pol2)) i interpolowanej pod pik.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri


class Gaussianizer:
    """Cecha → z = Φ⁻¹(F̂(x)); odwrotnie przez interpolację kwantyli. NaN → mediana (flaga w nan_frac)."""

    def __init__(self, n_q=2001, seed=0):
        self.n_q, self.rng = n_q, np.random.default_rng(seed)

    def fit(self, X: np.ndarray):
        X = np.asarray(X, float)
        self.med = np.nanmedian(X, 0)
        self.nan_frac = np.mean(~np.isfinite(X), 0)
        X = np.where(np.isfinite(X), X, self.med)
        p = np.linspace(0, 1, self.n_q)
        self.q = np.quantile(X, p, axis=0)  # (n_q, d)
        self.p = p
        return self

    def to_z(self, X):
        X = np.asarray(X, float)
        X = np.where(np.isfinite(X), X, self.med)
        U = np.empty_like(X)
        for j in range(X.shape[1]):
            q = self.q[:, j]
            lo = np.searchsorted(q, X[:, j], "left")
            hi = np.searchsorted(q, X[:, j], "right")
            # drganie w przedziale wartości równych (atomy: npts, obcięcia chi2/impact)
            u = (lo + self.rng.random(len(X)) * np.maximum(hi - lo, 1)) / (len(q))
            U[:, j] = u
        U = np.clip(U, 1e-6, 1 - 1e-6)
        return ndtri(U)

    def from_z(self, Z):
        U = ndtr(np.asarray(Z, float))
        X = np.empty_like(U)
        for j in range(U.shape[1]):
            X[:, j] = np.interp(U[:, j], self.p, self.q[:, j])
        return X


class MassDensity:
    """Gęstość tła w m: exp(pol2) dopasowana do histogramu sidebandów (Poisson, iminuit) → próbkowanie w zakresie."""

    def __init__(self, lo, hi, sb_mask_fn, nbins=50):
        self.lo, self.hi, self.sb, self.nbins = lo, hi, sb_mask_fn, nbins

    def fit(self, m):
        from iminuit import Minuit
        h, e = np.histogram(m, self.nbins, (self.lo, self.hi))
        c = 0.5 * (e[1:] + e[:-1])
        use = self.sb(c)
        x = (c - self.lo) / (self.hi - self.lo)

        def nll(a0, a1, a2):
            mu = np.exp(a0 + a1 * x[use] + a2 * x[use] ** 2)
            return float(np.sum(mu - h[use] * np.log(mu + 1e-300)))
        mi = Minuit(nll, a0=np.log(max(h[use].mean(), 1)), a1=0.0, a2=0.0)
        mi.errordef = Minuit.LIKELIHOOD
        mi.migrad()
        self.a = np.array(mi.values)
        self.edges, self.h = e, h
        return self

    def pdf(self, m):
        x = (np.asarray(m) - self.lo) / (self.hi - self.lo)
        return np.exp(self.a[0] + self.a[1] * x + self.a[2] * x ** 2)

    def expected(self, a, b, n=2001):
        g = np.linspace(a, b, n)
        return float(np.trapezoid(self.pdf(g), g) / (self.edges[1] - self.edges[0]))

    def sample(self, n, a, b, seed=0):
        rng = np.random.default_rng(seed)
        g = np.linspace(a, b, 4001)
        cdf = np.cumsum(self.pdf(g))
        cdf /= cdf[-1]
        return np.interp(rng.random(n), cdf, g)


def _mstd(m, m0):
    return (np.asarray(m, float) - m0) / 0.1


class CopulaGen:
    """z | m ~ N(μ(m), Σ(m)); μ(m), Σ(m) liniowe w m (2 punkty: lewy/prawy sideband ⇒ interpolacja liniowa)."""

    def __init__(self, m0, nbins_m=6, seed=0):
        self.m0, self.nb, self.seed = m0, nbins_m, seed

    def fit(self, Z, m):
        t = _mstd(m, self.m0)
        e = np.unique(np.quantile(t, np.linspace(0, 1, self.nb + 1)))
        tc, mus, covs = [], [], []
        for a, b in zip(e[:-1], e[1:]):
            s = (t >= a) & (t <= b)
            if s.sum() < 200:
                continue
            tc.append(t[s].mean())
            mus.append(Z[s].mean(0))
            covs.append(np.cov(Z[s], rowvar=False))
        tc, mus, covs = np.array(tc), np.array(mus), np.array(covs)
        A = np.c_[np.ones_like(tc), tc]
        self.mu_coef = np.linalg.lstsq(A, mus, rcond=None)[0]           # (2, d)
        self.cov_coef = np.linalg.lstsq(A, covs.reshape(len(tc), -1), rcond=None)[0]
        self.d = Z.shape[1]
        return self

    def sample(self, m, seed=None):
        rng = np.random.default_rng(self.seed if seed is None else seed)
        t = _mstd(m, self.m0)
        out = np.empty((len(t), self.d))
        # Σ(m) w siatce 50 punktów (szybko), Cholesky z regularyzacją
        g = np.linspace(t.min(), t.max() + 1e-9, 51)
        ib = np.clip(np.digitize(t, g) - 1, 0, 49)
        for k in np.unique(ib):
            s = ib == k
            tc = 0.5 * (g[k] + g[k + 1])
            mu = self.mu_coef[0] + self.mu_coef[1] * tc
            C = (self.cov_coef[0] + self.cov_coef[1] * tc).reshape(self.d, self.d)
            C = 0.5 * (C + C.T)
            w, V = np.linalg.eigh(C)
            L = V * np.sqrt(np.clip(w, 1e-4, None))
            out[s] = mu + rng.standard_normal((s.sum(), self.d)) @ L.T
        return out


class FlowGen:
    """Neural spline flow (zuko NSF), kontekst = m. Uczenie na sidebandach w przestrzeni z."""

    def __init__(self, m0, transforms=4, hidden=(128, 128), bins=8, epochs=30, batch=2048, lr=2e-3, seed=0,
                 val_frac=0.1, patience=5, verbose=True):
        self.m0, self.kw = m0, dict(transforms=transforms, hidden_features=hidden, bins=bins)
        self.epochs, self.batch, self.lr, self.seed = epochs, batch, lr, seed
        self.val_frac, self.patience, self.verbose = val_frac, patience, verbose

    def fit(self, Z, m):
        import torch
        import zuko
        torch.manual_seed(self.seed)
        d = Z.shape[1]
        self.flow = zuko.flows.NSF(features=d, context=1, **self.kw)
        X = torch.tensor(Z, dtype=torch.float32)
        C = torch.tensor(_mstd(m, self.m0)[:, None], dtype=torch.float32)
        g = torch.Generator().manual_seed(self.seed)
        perm = torch.randperm(len(X), generator=g)
        nv = int(len(X) * self.val_frac)
        iv, it = perm[:nv], perm[nv:]
        opt = torch.optim.AdamW(self.flow.parameters(), lr=self.lr, weight_decay=1e-5)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, self.epochs)
        best, bad, self.hist = np.inf, 0, []
        best_state = None
        for ep in range(self.epochs):
            self.flow.train()
            p = it[torch.randperm(len(it), generator=g)]
            tl = 0.0
            for k in range(0, len(p), self.batch):
                b = p[k:k + self.batch]
                loss = -self.flow(C[b]).log_prob(X[b]).mean()
                opt.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.flow.parameters(), 5.0)
                opt.step()
                tl += float(loss.detach()) * len(b)
            sched.step()
            self.flow.eval()
            with torch.no_grad():
                vl = float(-self.flow(C[iv]).log_prob(X[iv]).mean())
            self.hist.append((ep, tl / len(it), vl))
            if self.verbose:
                print(f"  flow ep {ep}: train {tl/len(it):.4f} val {vl:.4f}", flush=True)
            if vl < best - 1e-4:
                best, bad = vl, 0
                best_state = {k: v.clone() for k, v in self.flow.state_dict().items()}
            else:
                bad += 1
                if bad >= self.patience:
                    break
        if best_state is not None:
            self.flow.load_state_dict(best_state)
        self.best_val = best
        return self

    def sample(self, m, seed=None, chunk=5000):  # NSF: ~110 MB / 1k próbek przy odwracaniu
        import torch
        torch.manual_seed(self.seed + 1 if seed is None else seed)
        t = torch.tensor(_mstd(m, self.m0)[:, None], dtype=torch.float32)
        out = []
        with torch.no_grad():
            for k in range(0, len(t), chunk):
                out.append(self.flow(t[k:k + chunk]).sample().numpy())
        Z = np.concatenate(out)
        return np.clip(Z, -6, 6)

    def log_prob(self, Z, m, chunk=20000):
        import torch
        out = []
        with torch.no_grad():
            for k in range(0, len(Z), chunk):
                x = torch.tensor(Z[k:k + chunk], dtype=torch.float32)
                c = torch.tensor(_mstd(m[k:k + chunk], self.m0)[:, None], dtype=torch.float32)
                out.append(self.flow(c).log_prob(x).numpy())
        return np.concatenate(out)


class BkgGenerator:
    """Całość: Gaussianizer + model warunkowy + gęstość masy. fit(sidebandy) → sample(n, okno)."""

    def __init__(self, feats, m0, sb_fn, mass_range, model="flow", seed=0, **kw):
        self.feats, self.m0, self.sb_fn, self.mass_range = list(feats), m0, sb_fn, mass_range
        self.model_name, self.seed, self.kw = model, seed, kw

    def fit(self, SB: pd.DataFrame, m_all_for_density=None):
        X = SB[self.feats].to_numpy(float)
        m = SB["m"].to_numpy(float)
        self.G = Gaussianizer(seed=self.seed).fit(X)
        Z = self.G.to_z(X)
        if self.model_name == "copula":
            self.M = CopulaGen(self.m0, seed=self.seed).fit(Z, m)
        else:
            self.M = FlowGen(self.m0, seed=self.seed, **self.kw).fit(Z, m)
        lo, hi = self.mass_range
        self.mass = MassDensity(lo, hi, self.sb_fn).fit(m if m_all_for_density is None else m_all_for_density)
        self.mass.sb = self.sb_fn = None  # lambda nie jest picklowalna
        return self

    def sample(self, n, window, seed=None, m=None):
        m = self.mass.sample(n, *window, seed=self.seed + 7 if seed is None else seed) if m is None else np.asarray(m)
        Z = self.M.sample(m, seed=seed)
        X = self.G.from_z(Z)
        df = pd.DataFrame(X, columns=self.feats)
        df["m"] = m
        if "log_pt" in df:
            df["pt"] = np.expm1(df["log_pt"])
        return df
