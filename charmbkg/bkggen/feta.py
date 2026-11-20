"""charmbkg.bkggen.feta — tło pod pikiem D z MC skorygowanego do danych w sidebandach (etap S4).

Metody (warunkowane masą m kandydata, korekta uczona w SB, stosowana w SR):
  * FETA (2212.11285): flow bazowy p_MC(z|m) na MC tła (cały zakres m, MC nie ma piku), zamrożony;
    flow transportujący T(·|m) (NSF, inicjalizacja = tożsamość) uczony ML na danych SB:
        log p_dane(z|m) = log p_MC(T(z)|m) + log|det ∂T/∂z|.
    Tło syntetyczne w SR = T⁻¹(z_MC | m_MC) dla PRAWDZIWYCH zdarzeń MC w SR (m bez zmian → więzy kinematyczne
    spełnione, cechy H z derive_H).
  * SALAD (2001.05001): klasyfikator dane-SB vs MC-SB z m na wejściu → w(x,m) = p/(1−p)·N_MC/N_dane,
    wagi stosowane do MC w SR. Warianty: 'mlp' (gładki w m → interpolacja), 'hgb' (drzewa, m jako cecha),
    'hgb_nom' (bez m: korekta niezależna od m).
Wspólne: Preproc z flow.py (PIT→N(0,1)) dopasowany na MC; cechy generowane FEATS_GEN (kin.py, prep.py — S2).
"""
from __future__ import annotations

import logging
import time

import numpy as np

from . import flow as F

log = logging.getLogger("charmbkg.bkggen.feta")
SWAP_VETO = 0.025   # |m_swap − m0| > 25 MeV (weto refleksji K↔π, etap S)


def swap_mass_kin(m, pt, eta, cts, phis):
    """m(K↔π zamienione) z (m, pT, η, cosθ*, φ*) — te same zmienne, które generuje/transportuje model,
    więc weto działa identycznie na danych, MC i tle syntetycznym."""
    from . import kin
    m, pt, eta = (np.asarray(v, float) for v in (m, pt, eta))
    cts = np.clip(np.asarray(cts, float), -1, 1)
    n = len(m)
    pD = kin._vec(pt, eta, np.zeros(n))
    E = np.sqrt((pD ** 2).sum(-1) + m ** 2)
    ps = np.sqrt(np.clip((m ** 2 - (kin.MK + kin.MPI) ** 2) * (m ** 2 - (kin.MK - kin.MPI) ** 2), 0, None)) / (2 * m)
    x, y, z = kin._frame(pD)
    st = np.sqrt(1 - cts ** 2)
    u = (st * np.cos(phis))[:, None] * x + (st * np.sin(phis))[:, None] * y + cts[:, None] * z
    b = pD / E[:, None]
    b2 = (b ** 2).sum(-1)
    g = 1 / np.sqrt(1 - b2)
    P = []
    for sgn, mm in ((1, kin.MK), (-1, kin.MPI)):
        q = sgn * ps[:, None] * u
        Eq = np.sqrt(ps ** 2 + mm ** 2)
        bq = (b * q).sum(-1)
        P.append(q + (((g - 1) * bq / b2) + g * Eq)[:, None] * b)
    Ea = np.sqrt((P[0] ** 2).sum(-1) + kin.MPI ** 2)      # tor „K” jako π
    Eb = np.sqrt((P[1] ** 2).sum(-1) + kin.MK ** 2)       # tor „π” jako K
    return np.sqrt(np.clip((Ea + Eb) ** 2 - ((P[0] + P[1]) ** 2).sum(-1), 0, None))


def swap_mass(df):
    return swap_mass_kin(df.m.values, df.pt.values, df.eta.values, df.cts.values, df.phis.values)


def swap_mass_G(G, m, gen=F.FEATS_GEN):
    c = {f: G[:, j] for j, f in enumerate(gen)}
    return swap_mass_kin(m, np.expm1(c["log_pt"]), c["eta"], c["cts"], c["phis"])


def veto_ok(df, sp, cut=SWAP_VETO):
    return np.abs(swap_mass(df) - sp["m0"]) > cut


def veto_ok_G(G, m, sp, cut=SWAP_VETO):
    return np.abs(swap_mass_G(G, m) - sp["m0"]) > cut


def band(m, sp, lo, hi):
    dm = np.abs(np.asarray(m) - sp["m0"])
    return (dm >= lo) & (dm < hi)


def _zero_last(nsf):
    import torch
    for t in nsf.transform.transforms:
        lins = [mm for mm in t.hyper.modules() if isinstance(mm, torch.nn.Linear)]
        torch.nn.init.zeros_(lins[-1].weight)
        torch.nn.init.zeros_(lins[-1].bias)


# --------------------------------------------------------------------------- FETA
class FetaTransport:
    """Flow bazowy na MC + flow transportujący MC↔dane (w przestrzeni z Preproc dopasowanego na MC)."""

    def __init__(self, feats, sp, base_transforms=6, base_hidden=(192, 192), top_transforms=3, top_hidden=(96, 96),
                 bins=8, seed=1, device=None):
        import torch
        import zuko
        self.torch = torch
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        torch.manual_seed(seed)
        self.feats, self.sp = list(feats), sp
        self.pre = F.Preproc(feats, seed=seed)
        self.base = zuko.flows.NSF(len(feats), 1, transforms=base_transforms, hidden_features=list(base_hidden), bins=bins)
        self.top = zuko.flows.NSF(len(feats), 1, transforms=top_transforms, hidden_features=list(top_hidden), bins=bins)
        _zero_last(self.top)          # T = tożsamość na starcie: MC bez korekty
        self.base.to(self.device); self.top.to(self.device)
        self.hist = {"base": [], "top": []}

    # -- pomocnicze
    def _ctx(self, m):
        return self.torch.from_numpy(F.m_context(m, self.sp).astype(np.float32)).to(self.device)

    def _z(self, X):
        return self.torch.from_numpy(self.pre.forward(X).astype(np.float32)).to(self.device)

    def _train(self, which, Z, C, loss_fn, params, epochs, batch, lr, max_minutes, patience=5, val_frac=0.2):
        torch = self.torch
        idx = np.random.default_rng(0).permutation(len(Z))
        nv = int(val_frac * len(Z))
        Zv, Cv, Zt, Ct = Z[idx[:nv]], C[idx[:nv]], Z[idx[nv:]], C[idx[nv:]]
        opt = torch.optim.AdamW(params, lr=lr, weight_decay=1e-5)
        sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
        mod = self.base if which == "base" else self.top
        best, best_state, bad, t0 = np.inf, None, 0, time.time()
        for ep in range(epochs):
            mod.train()
            perm = torch.randperm(len(Zt), device=self.device)
            tot = 0.0
            for i in range(0, len(Zt), batch):
                b = perm[i:i + batch]
                loss = loss_fn(Zt[b], Ct[b])
                opt.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(params, 5.0)
                opt.step()
                tot += float(loss.detach()) * len(b)
            sch.step()
            mod.eval()
            with torch.no_grad():
                vl = float(np.mean([loss_fn(Zv[i:i + 8192], Cv[i:i + 8192]).item() for i in range(0, len(Zv), 8192)]))
            self.hist[which].append((ep, tot / len(Zt), vl, time.time() - t0))
            log.info("%s epoka %d: train %.4f val %.4f (%.0f s)", which, ep, tot / len(Zt), vl, time.time() - t0)
            if vl < best - 1e-4:
                best, bad = vl, 0
                best_state = {k: v.detach().clone() for k, v in mod.state_dict().items()}
            else:
                bad += 1
                if bad >= patience:
                    break
            if max_minutes and time.time() - t0 > 60 * max_minutes:
                break
        if best_state is not None:
            mod.load_state_dict(best_state)
        return best

    def fit_base(self, X_mc, m_mc, epochs=40, batch=1024, lr=1e-3, max_minutes=20):
        """Flow bazowy p_MC(z|m) — na MC tła z całego zakresu m (SB+SR)."""
        self.pre.fit(X_mc)
        Z, C = self._z(X_mc), self._ctx(m_mc)
        nll = lambda z, c: -self.base(c).log_prob(z).mean()  # noqa: E731
        self.base_val = self._train("base", Z, C, nll, list(self.base.parameters()), epochs, batch, lr, max_minutes)
        for p in self.base.parameters():
            p.requires_grad_(False)
        return self

    def _logp_data(self, z, c):
        y, ladj = self.top(c).transform.call_and_ladj(z)
        return self.base(c).log_prob(y) + ladj

    def fit_top(self, X_data, m_data, epochs=40, batch=1024, lr=5e-4, max_minutes=20):
        """Transport T uczony na danych (tylko SB!)."""
        Z, C = self._z(X_data), self._ctx(m_data)
        with self.torch.no_grad():
            self.top_val0 = float(-self._logp_data(Z[:50_000], C[:50_000]).mean())   # T = tożsamość
        nll = lambda z, c: -self._logp_data(z, c).mean()  # noqa: E731
        self.top_val = self._train("top", Z, C, nll, list(self.top.parameters()), epochs, batch, lr, max_minutes)
        log.info("FETA top: NLL dane-SB tożsamość %.4f → po transporcie (val) %.4f", self.top_val0, self.top_val)
        return self

    def transport(self, X_mc, m_mc, chunk=4_000):
        """Prawdziwe zdarzenia MC → „dane”: z_dane = T⁻¹(z_MC | m)."""
        torch = self.torch
        Z, C = self._z(X_mc), self._ctx(m_mc)
        out = []
        self.top.eval()
        with torch.no_grad():
            for i in range(0, len(Z), chunk):
                out.append(self.top(C[i:i + chunk]).transform.inv(Z[i:i + chunk]).cpu().numpy())
        Zd = np.clip(np.concatenate(out).astype(np.float64), -6, 6)
        return self.pre.inverse(Zd)

    def sample_base(self, m, chunk=4000):
        """Próbki z flow bazowego (MC) — diagnostyka."""
        torch = self.torch
        C = self._ctx(np.asarray(m))
        with torch.no_grad():
            Z = np.concatenate([self.base(C[i:i + chunk]).sample().cpu().numpy() for i in range(0, len(C), chunk)])
        return self.pre.inverse(np.clip(Z.astype(np.float64), -6, 6))

    def save(self, path):
        self.torch.save({"base": self.base.state_dict(), "top": self.top.state_dict(), "ref": self.pre.ref,
                         "feats": self.feats, "sp": self.sp, "hist": self.hist}, path)

    def load(self, path, top=True):
        ck = self.torch.load(path, map_location=self.device, weights_only=False)
        self.base.load_state_dict(ck["base"])
        if top:
            self.top.load_state_dict(ck["top"])
        self.pre.ref = ck["ref"]
        n = len(self.pre.ref)
        self.pre.uc = (np.arange(n) + 0.5) / n
        self.hist = ck.get("hist", self.hist)
        for p in self.base.parameters():
            p.requires_grad_(False)
        return self

    def reset_top(self, top_transforms=3, top_hidden=(96, 96), bins=8):
        import zuko
        self.top = zuko.flows.NSF(len(self.feats), 1, transforms=top_transforms, hidden_features=list(top_hidden), bins=bins)
        _zero_last(self.top)
        self.top.to(self.device)
        self.hist["top"] = []
        return self


# --------------------------------------------------------------------------- SALAD
class SaladReweighter:
    """w(x, m) = p/(1−p)·N_MC/N_dane, klasyfikator dane(1) vs MC(0) uczony w SB.
    kind: 'mlp' (torch, wejście z-PIT + t=(m−m0)/0.1, gładki w m), 'hgb' (HistGB z m), 'hgb_nom' (HistGB bez m)."""

    def __init__(self, feats, sp, kind="mlp", clip=(0.05, 20.0), seed=0):
        self.feats, self.sp, self.kind, self.clip, self.seed = list(feats), sp, kind, clip, seed

    def _inp(self, X, m):
        Z = self.pre.forward(X)
        if self.kind == "hgb_nom":
            return Z
        return np.column_stack([Z, F.m_context(m, self.sp)[:, 0]])

    def fit(self, X_mc, m_mc, X_data, m_data, epochs=12, max_n=400_000):
        rng = np.random.default_rng(self.seed)
        self.pre = F.Preproc(self.feats, seed=self.seed).fit(X_mc)
        i_mc = rng.choice(len(X_mc), min(len(X_mc), max_n), replace=False)
        i_d = rng.choice(len(X_data), min(len(X_data), max_n), replace=False)
        A, B = self._inp(X_mc[i_mc], m_mc[i_mc]), self._inp(X_data[i_d], m_data[i_d])
        X = np.vstack([A, B]).astype(np.float32)
        y = np.r_[np.zeros(len(A)), np.ones(len(B))].astype(np.float32)
        w = np.r_[np.full(len(A), len(B) / len(A)), np.ones(len(B))].astype(np.float32)   # zbalansowane klasy
        if self.kind == "mlp":
            import torch
            torch.manual_seed(self.seed)
            d = X.shape[1]
            self.net = torch.nn.Sequential(torch.nn.Linear(d, 128), torch.nn.SiLU(), torch.nn.Linear(128, 128),
                                           torch.nn.SiLU(), torch.nn.Linear(128, 64), torch.nn.SiLU(), torch.nn.Linear(64, 1))
            opt = torch.optim.AdamW(self.net.parameters(), lr=1e-3, weight_decay=1e-4)
            Xt, yt, wt = map(torch.from_numpy, (X, y, w))
            perm0 = torch.randperm(len(Xt))
            nv = len(Xt) // 5
            iv, it = perm0[:nv], perm0[nv:]
            lossf = torch.nn.BCEWithLogitsLoss(reduction="none")
            best, best_state = np.inf, None
            for ep in range(epochs):
                self.net.train()
                p = it[torch.randperm(len(it))]
                for i in range(0, len(p), 2048):
                    b = p[i:i + 2048]
                    loss = (lossf(self.net(Xt[b])[:, 0], yt[b]) * wt[b]).mean()
                    opt.zero_grad(); loss.backward(); opt.step()
                self.net.eval()
                with torch.no_grad():
                    vl = float((lossf(self.net(Xt[iv])[:, 0], yt[iv]) * wt[iv]).mean())
                log.info("SALAD mlp epoka %d val %.5f", ep, vl)
                if vl < best:
                    best, best_state = vl, {k: v.clone() for k, v in self.net.state_dict().items()}
            self.net.load_state_dict(best_state)
            self.val_loss = best
        else:
            from sklearn.ensemble import HistGradientBoostingClassifier
            self.clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08, max_leaf_nodes=31,
                                                      l2_regularization=1.0, early_stopping=True, random_state=self.seed)
            self.clf.fit(X, y, sample_weight=w)
        return self

    def logit(self, X, m):
        Z = self._inp(X, m).astype(np.float32)
        if self.kind == "mlp":
            import torch
            with torch.no_grad():
                return np.concatenate([self.net(torch.from_numpy(Z[i:i + 65536]))[:, 0].numpy()
                                       for i in range(0, len(Z), 65536)]).astype(np.float64)
        p = np.clip(self.clf.predict_proba(Z)[:, 1], 1e-6, 1 - 1e-6)
        return np.log(p / (1 - p))

    def weights(self, X, m):
        """Wagi znormalizowane do średniej 1 (normalizację tła w SR daje fit masy / MassModel)."""
        w = np.clip(np.exp(self.logit(X, m)), *self.clip)
        return w / w.mean()


# --------------------------------------------------------------------------- walidacja
def c2st(Xa, Xb, wa=None, wb=None, n_folds=3, max_n=150_000, seed=0):
    """AUC klasyfikatora (HistGB, OOF) a(1) vs b(0) z wagami; cel 0,5. Klasy zbalansowane wagami."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    rng = np.random.default_rng(seed)
    wa = np.ones(len(Xa)) if wa is None else np.asarray(wa, float)
    wb = np.ones(len(Xb)) if wb is None else np.asarray(wb, float)
    ia = rng.choice(len(Xa), min(len(Xa), max_n), replace=False)
    ib = rng.choice(len(Xb), min(len(Xb), max_n), replace=False)
    X = np.vstack([Xa[ia], Xb[ib]])
    X[~np.isfinite(X)] = -1.0
    y = np.r_[np.ones(len(ia)), np.zeros(len(ib))]
    w = np.r_[wa[ia] / wa[ia].sum(), wb[ib] / wb[ib].sum()] * len(y) / 2
    s = np.zeros(len(y))
    for tr, te in StratifiedKFold(n_folds, shuffle=True, random_state=seed).split(X, y):
        clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, max_leaf_nodes=31,
                                             early_stopping=True, random_state=seed)
        clf.fit(X[tr], y[tr], sample_weight=w[tr])
        s[te] = clf.predict_proba(X[te])[:, 1]
    return float(roc_auc_score(y, s, sample_weight=w))


def wmean_ks(a, b, wa=None, wb=None):
    """Ważona statystyka KS (max |ΔCDF|)."""
    wa = np.ones(len(a)) if wa is None else wa
    wb = np.ones(len(b)) if wb is None else wb
    ka, kb = np.isfinite(a), np.isfinite(b)
    a, wa, b, wb = a[ka], wa[ka], b[kb], wb[kb]
    g = np.unique(np.quantile(np.r_[a, b], np.linspace(0, 1, 401)))
    ca = np.array([wa[a <= x].sum() for x in g]) / wa.sum()
    cb = np.array([wb[b <= x].sum() for x in g]) / wb.sum()
    return float(np.max(np.abs(ca - cb)))
