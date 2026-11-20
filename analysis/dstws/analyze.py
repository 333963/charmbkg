#!/usr/bin/env python3
"""DSTWS — analiza: tło D* przez konstrukcję (D0 × tor z Tracks: RS / WS / ROT / MIX) jako szablon kształtu tła w Δm (projekt B, 17.11.2026).

KRYTERIA (zapisane PRZED obejrzeniem danych; liczby progów nie są dostrajane po wynikach):
 (a) Walidacja RS vs nominalne D* z ntupla:
     K-a1  pary nominalne z pT(π_s) ≥ 0,2 GeV odnalezione wśród naszych RS: ≥ 99 % (per bin pT);
     K-a2  ⟨Δm_własne − Δm_nominalne⟩: |średnia| < 0,05 MeV, RMS < 0,30 MeV;
     K-a3  N_sig z fitu nominalnego (próg × pol2) na naszych RS vs na parach nominalnych (ten sam próg pT π_s ≥ 0,2): zgodne w 2σ(różnicy niezal.) —
           informacyjnie, bo zbiory silnie skorelowane.
     Raportowane bez progu: udział RS spoza nominału (luźniejsza selekcja wolnego pionu), strata nominału przez pT(π_s) < 0,2 GeV.
 (b) MC (mb_nd = skład tła jak w danych; fd_nd = kontrola): kształt Δm szablonu T ∈ {WS, ROT, MIX} vs prawdziwe tło RS
     (RS i NIE [D0 prawda ∈ {0, 1} i tor = prawdziwy wolny pion D*]) na [140,47; 160] MeV, biny 0,5 MeV, norma swobodna:
     K-b1  χ²/ndf < 2 i p > 0,01 w każdym binie pT × {wszystkie, pass WP80, fail WP80} → „szablon odtwarza tło”;
     K-b2  to samo po poprawce liniowej (1 + c·x) → „odtwarza po poprawce pol1”; inaczej „nie odtwarza”;
     K-b3  przeciek sygnału do WS (tor = prawdziwy wolny pion D*, bliźniak K↔π): N_WS,przeciek / N_RS,sygnał < 5 % (inaczej korekta/odrzucenie WS).
 (c) Dane, per run (428770 wszystkie; 427929 GRL ∧ mb_mbts_L1MBTS_1), bin pT, WP 95/80/50: fit Δm pass/fail z tłem = szablon (kształt ustalony,
     norma swobodna; wariant +pol1) vs nominał `charmtag.dstar_refit.fit_passfail` (próg × pol2) na TYCH SAMYCH histogramach RS:
     K-c1  bin „mierzalny”: fit określony (valid, Hesse dodatnio określona, N_sig > 3σ, ε poza granicami, MINOS ε zamknięty) i σ(SF) < 0,15
           (σ = stat ⊕ statystyka szablonu ⊕ stat MC);
     K-c2  zysk: σ(N_sig)_szablon / σ(N_sig)_nominał < 0,8;
     K-c3  systematyka wyboru szablonu: max |SF_T − SF_T'| po T, T' ∈ {WS, ROT, MIX} (tylko te, które zdały K-b) — podawana obok σ_stat;
           „stabilny”, gdy < σ_stat.
     K-c4  funkcja progowa: fit samego tła (próg × pol2) do WS w danych — χ²/ndf < 2 i p > 0,01 → „funkcja progowa właściwa”.
 (d) Domknięcie MC (mb_nd, fd_nd): N_sig z fitu szablonowego vs prawda (liczba RS-sygnał w zakresie fitu):
     K-d1  |N_fit − N_true| < 2σ_fit i |N_fit/N_true − 1| < 10 %; ε: |ε_fit − ε_true| < 2σ.

Wejście: --in KATALOG z dw_<próbka>.npz (dstws_merge.py). Progi WP: kwantyle score prawdziwych D0 w fd (bt_fd_nd.npz z B_trig, jak btrig_analyze).
Wyjście: --out KATALOG → dstws.json + RAPORT.md (--report ścieżka).
Użycie: .venv_h5/bin/python charm_work/DSTWS/analyze.py [--in charm_work/DSTWS/runs/lxplus] [--fast] [--nproc 8]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CW, "AA_merge", "src"))

M0 = 1.86484
MPI = 139.57
LO, HI = MPI + 0.9, 160.0
E = np.arange(139.6, 160.0001, 0.1)                       # MeV, siatka jak B_trig/EV (do 160 — zakres fitu)
BINS = {"3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0)}
WPS = ["0.95", "0.80", "0.50"]
TYPES = {"RS": 0, "WS": 1, "ROT": 2, "MIX": 3}
TPL = ["WS", "ROT", "MIX"]
NOM_WIN = 0.055
RUNS = {"428770": ("data_428770", 0, 0), "427929": ("data_427929", 36, 36)}      # (próbka, maska bitów trig, wartość)
SAMPLES = ["data_428770", "data_427929", "fd_nd", "mb_nd"]


# ----------------------------------------------------------------------------------------------------------- wejście
class Sample:
    def __init__(self, path):
        z = np.load(path, allow_pickle=False)
        self.meta = json.loads(str(z["meta"]))
        self.z = {k: z[k] for k in z.files if k != "meta"}
        i = self.z["c_d0"]
        g = lambda k: self.z["d_" + k][i]  # noqa: E731
        self.dm = self.z["c_dm"].astype(np.float64) * 1e3
        self.typ = self.z["c_typ"]; self.k = self.z["c_k"]; self.flag = self.z["c_flag"]; self.pis = self.z["c_pis_pt"]
        self.pt = g("pt"); self.m = g("m"); self.score = g("score"); self.truth = g("truth"); self.trig = g("trig"); self.hs = g("hypo_swap")
        self.nom = np.abs(self.m - M0) < NOM_WIN
        self.d0true = np.isin(self.truth, [0, 1])
        self.sig = (self.typ == 0) & self.d0true & ((self.flag & 2) != 0)        # MC: prawdziwe D*→D0π, D0→Kπ z poprawną hipotezą

    def bin(self, b):
        lo, hi = BINS[b]
        return (self.pt >= lo) & (self.pt < hi)

    def hist(self, mask):
        return np.histogram(self.dm[mask], E)[0].astype(float)


def thresholds(bt_fd, fallback):
    """Progi WP per bin pT: kwantyle score prawdziwych D0 (fd, okno ±55 MeV) — jak btrig_analyze.thresholds."""
    if bt_fd and os.path.exists(bt_fd):
        z = np.load(bt_fd, allow_pickle=False)
        sc, pt, tr, m = z["all_score"], z["all_pt"], z["all_truth"], z["all_m"]
        sig = np.isin(tr, [0, 1]) & (np.abs(m - M0) < NOM_WIN)
        return {b: {wp: float(np.quantile(sc[sig & (pt >= lo) & (pt < hi)], 1 - float(wp))) for wp in WPS} for b, (lo, hi) in BINS.items()}, bt_fd
    return {b: {wp: float(fallback[b][wp]) for wp in WPS if wp in fallback[b]} for b in BINS}, "meta modelu (brak bt_fd_nd.npz; bez WP 0.95)"


# ----------------------------------------------------------------------------------------------------------- fity
def _sel(e):
    s = (e[:-1] >= LO - 1e-9) & (e[1:] <= HI + 1e-9)
    return s, np.r_[e[:-1][s], e[1:][s][-1]]


def tpl_density(h, e, group=5):
    """Szablon → gładka gęstość na binach fitu: grupy po `group` binów (0,5 MeV), interpolacja liniowa między środkami grup, norma 1."""
    s, ee = _sel(e)
    hh = np.asarray(h, float)[s]
    n = len(hh) // group * group
    g = hh[:n].reshape(-1, group).sum(1) / group
    cg = 0.5 * (ee[:-1][:n].reshape(-1, group)[:, 0] + ee[1:][:n].reshape(-1, group)[:, -1])
    if n < len(hh):
        g = np.r_[g, hh[n:].mean()]; cg = np.r_[cg, 0.5 * (ee[n] + ee[-1])]
    c = 0.5 * (ee[1:] + ee[:-1])
    d = np.clip(np.interp(c, cg, g), 1e-6, None)
    return d / d.sum()


def _nll(lam, n):
    lam = np.clip(lam, 1e-9, None)
    pos = n > 0
    out = lam.sum() - n.sum()
    out -= np.sum(n[pos] * np.log(lam[pos] / n[pos]))
    return float(out)


def fit_tpl(hp, hf, tp, tf, e, shape, pol1=False, nboot=0, seed=1):
    """Jednoczesny fit pass/fail: sygnał 2 Gaussy (kształt MC, wspólne μ i skala σ), tło = nb · T(Δm)[·(1 + c·x)], T z szablonu (ustalony).
    Zwraca dict (eps, eps_err, N, N_err, chi2, ndf, determined, …). nboot > 0: rozrzut ε i N przy fluktuacji Poissona szablonu (stat. szablonu)."""
    from iminuit import Minuit
    from scipy.stats import norm
    s, ee = _sel(e)
    hp_, hf_ = np.asarray(hp, float)[s], np.asarray(hf, float)[s]
    c = 0.5 * (ee[1:] + ee[:-1]); x = (c - c.mean()) / 10.0
    base = {"n_pass": int(hp_.sum()), "n_fail": int(hf_.sum()), "n_tpl_pass": float(np.asarray(tp)[s].sum()), "n_tpl_fail": float(np.asarray(tf)[s].sum())}
    if hp_.sum() + hf_.sum() < 50 or base["n_tpl_pass"] < 20 or base["n_tpl_fail"] < 20:
        return {**base, "eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False, "determined": False, "reason": "za mało wpisów"}
    sc, k, f = shape["sc"], shape["k"], shape["f"]

    def sig(mu, sg):
        cdf = f * norm.cdf(ee, mu, sc * sg) + (1 - f) * norm.cdf(ee, mu, sc * sg * k)
        return np.diff(cdf) / (cdf[-1] - cdf[0])

    def one(dp, df_, full=True):
        def lam(Nt, eps, mu, sg, nbp, nbf, cp, cf):
            S = sig(mu, sg)
            bp = dp * np.clip(1 + cp * x, 1e-3, None); bf = df_ * np.clip(1 + cf * x, 1e-3, None)
            return Nt * eps * S + nbp * bp / bp.sum(), Nt * (1 - eps) * S + nbf * bf / bf.sum()

        def nll(Nt, eps, mu, sg, nbp, nbf, cp, cf):
            lp, lf = lam(Nt, eps, mu, sg, nbp, nbf, cp, cf)
            return _nll(lp, hp_) + _nll(lf, hf_)
        win = np.abs(c - shape["mu"]) < 2.5
        exc = lambda h, d: max(h[win].sum() - h[~win].sum() * d[win].sum() / max(d[~win].sum(), 1e-9), 1.0)  # noqa: E731
        sp0, sf0 = exc(hp_, dp), exc(hf_, df_)
        best = None
        for mult in (1.0, 0.3, 3.0):
            n0 = max((sp0 + sf0) * mult, 10.0)
            m = Minuit(nll, Nt=n0, eps=float(np.clip(sp0 / (sp0 + sf0), 0.05, 0.95)), mu=shape["mu"], sg=1.0,
                       nbp=max(hp_.sum() - sp0 * mult, 1.0), nbf=max(hf_.sum() - sf0 * mult, 1.0), cp=0.0, cf=0.0)
            m.errordef = 0.5; m.strategy = 2
            m.limits["Nt"] = (0, None); m.limits["eps"] = (-0.2, 1.5); m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
            m.limits["sg"] = (0.8, 1.6); m.limits["nbp"] = (0, None); m.limits["nbf"] = (0, None)
            m.limits["cp"] = (-0.9, 0.9); m.limits["cf"] = (-0.9, 0.9)
            if not pol1:
                m.fixed["cp"] = True; m.fixed["cf"] = True
            try:
                m.migrad(ncall=20000)
                if full:
                    m.hesse()
            except (RuntimeError, ValueError):
                continue
            if best is None or (m.valid and not best.valid) or (m.valid == best.valid and m.fval < best.fval):
                best = m
            if not full and best is not None and best.valid:
                break
        return best, lam

    m, lam = one(tpl_density(tp, e), tpl_density(tf, e))
    if m is None:
        return {**base, "eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False, "determined": False, "reason": "brak zbieżności"}
    v = {p: float(m.values[p]) for p in m.parameters}
    lp, lf = lam(*[v[p] for p in m.parameters])
    chi2 = float(np.sum((hp_ - lp) ** 2 / np.maximum(lp, 1)) + np.sum((hf_ - lf) ** 2 / np.maximum(lf, 1)))
    lo = hi = np.nan
    try:
        m.minos("eps")
        me = m.merrors["eps"]
        if me.is_valid and not me.at_lower_limit and not me.at_upper_limit:
            lo, hi = v["eps"] + me.lower, v["eps"] + me.upper
    except Exception:  # noqa: BLE001
        pass
    posdef = bool(m.fmin.has_posdef_covar); at_lim = bool(v["eps"] + 0.2 < 1e-3 or 1.5 - v["eps"] < 1e-3)
    Nerr = float(m.errors["Nt"])
    why = [w for w, bad in (("invalid", not m.valid), ("not_posdef", not posdef), ("eps_at_limit", at_lim), ("interval_open", not np.isfinite(lo + hi)),
                            ("low_signal", not v["Nt"] > 3.0 * Nerr)) if bad]
    out = {**base, "eps": v["eps"], "eps_err": float(m.errors["eps"]), "eps_lo": lo, "eps_hi": hi, "N": v["Nt"], "N_err": Nerr, "mu": v["mu"],
           "scale": v["sg"], "cp": v["cp"], "cf": v["cf"], "nbp": v["nbp"], "nbf": v["nbf"], "valid": bool(m.valid), "posdef": posdef, "chi2": chi2,
           "ndf": int(2 * len(hp_) - m.nfit), "determined": not why, "reason": ",".join(why), "eps_tpl": 0.0, "N_tpl": 0.0}
    if nboot and m.valid:
        rng = np.random.default_rng(seed); ev, nv = [], []
        for _ in range(nboot):
            mb, _l = one(tpl_density(rng.poisson(np.asarray(tp, float)), e), tpl_density(rng.poisson(np.asarray(tf, float)), e), full=False)
            if mb is not None and mb.valid:
                ev.append(float(mb.values["eps"])); nv.append(float(mb.values["Nt"]))
        if len(ev) >= 3:      # fluktuowany szablon = szablon + szum → rozrzut zawiera stat. szablonu raz
            out["eps_tpl"], out["N_tpl"] = float(np.std(ev)), float(np.std(nv))
    return out


def sym_err(f):
    lo, hi = f.get("eps_lo", np.nan), f.get("eps_hi", np.nan)
    e = 0.5 * (hi - lo) if np.isfinite(lo) and np.isfinite(hi) else f.get("eps_err", np.nan)
    return float(np.hypot(e, f.get("eps_tpl", 0.0))) if np.isfinite(e) else np.nan


def fit_nom(hp, hf, shape, interval):
    """Nominał: charmtag.dstar_refit.fit_passfail (próg × pol2) — jak btrig_common.fit_eps."""
    from charmtag import dstar_refit as DR
    base = {"n_pass": int(hp.sum()), "n_fail": int(hf.sum())}
    if hp.sum() + hf.sum() < 50:
        return {**base, "eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False, "determined": False, "reason": "za mało wpisów"}
    try:
        r = DR.fit_passfail(hp, hf, E, shape, interval=interval, model_check=False)
    except Exception as ex:  # noqa: BLE001
        return {**base, "eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False, "determined": False, "reason": str(ex)}
    g = lambda k: float(r.get(k, np.nan)) if r.get(k) is not None else np.nan  # noqa: E731
    return {**base, "eps": g("eps"), "eps_err": g("eps_err"), "eps_lo": g("eps_lo"), "eps_hi": g("eps_hi"), "N": g("N_tot"), "N_err": g("N_tot_err"),
            "valid": bool(r.get("valid", False)), "determined": bool(r.get("determined", False)), "reason": r.get("undetermined_reason"),
            "chi2": g("chi2"), "ndf": int(r.get("ndf", 0)), "bkg_model": r.get("bkg_model")}


def fit_bkg_only(h):
    """Fit samego tła u^a(1 + c1·t + c2·t²) (funkcja progowa nominału) do histogramu bez sygnału → χ², ndf, p."""
    from charmtag import dstar_refit as DR
    from iminuit import Minuit
    from scipy.stats import chi2 as C2
    s, ee = _sel(E)
    hh = np.asarray(h, float)[s]
    if hh.sum() < 200:
        return {"n": float(hh.sum()), "chi2": np.nan, "ndf": 0, "p": np.nan}
    best = None
    for a0, c10, c20 in ((0.5, 0.0, 0.0), (1.0, 5.0, -3.0), (0.3, -2.0, 1.0), (2.0, 0.0, 0.0)):
        def nll(a, c1, c2):
            return _nll(hh.sum() * np.diff(DR._bkg_cdf_t(ee, a, 0.0, c1, c2, False)), hh)
        m = Minuit(nll, a=a0, c1=c10, c2=c20); m.errordef = 0.5
        m.limits["a"] = (0.01, 3); m.limits["c1"] = (-5, 20); m.limits["c2"] = (-10, 40)
        try:
            m.migrad(ncall=20000)
        except (RuntimeError, ValueError):
            continue
        if best is None or m.fval < best.fval:
            best = m
    if best is None:
        return {"n": float(hh.sum()), "chi2": np.nan, "ndf": 0, "p": np.nan}
    lam = hh.sum() * np.diff(DR._bkg_cdf_t(ee, best.values["a"], 0.0, best.values["c1"], best.values["c2"], False))
    g = 5; n = len(hh) // g * g
    o = hh[:n].reshape(-1, g).sum(1); p = lam[:n].reshape(-1, g).sum(1)
    x2 = float(np.sum((o - p) ** 2 / np.maximum(p, 1))); ndf = len(o) - 3
    return {"n": float(hh.sum()), "chi2": x2, "ndf": ndf, "p": float(C2.sf(x2, ndf)), "a": float(best.values["a"]), "c1": float(best.values["c1"]),
            "c2": float(best.values["c2"])}


def shape_chi2(ha, hb, group=5):
    """Zgodność kształtu dwóch histogramów (norma swobodna) na zakresie fitu, biny 0,5 MeV; także z poprawką liniową (1 + c·x) do hb."""
    from scipy.stats import chi2 as C2
    s, ee = _sel(E)
    a, b = np.asarray(ha, float)[s], np.asarray(hb, float)[s]
    n = len(a) // group * group
    a = a[:n].reshape(-1, group).sum(1); b = b[:n].reshape(-1, group).sum(1)
    if a.sum() < 100 or b.sum() < 100:
        return {"na": float(a.sum()), "nb": float(b.sum()), "chi2": np.nan, "ndf": 0, "p": np.nan, "chi2_pol1": np.nan, "p_pol1": np.nan, "c": np.nan}
    x = np.linspace(-1, 1, len(a))
    r = a.sum() / b.sum()
    var = a + r * r * b
    x2 = float(np.sum((a - r * b) ** 2 / np.maximum(var, 1)))
    w = 1 / np.maximum(var, 1)                      # a ≈ b·(p0 + p1·x): ważone najmniejsze kwadraty (liniowe)
    A = np.column_stack([b, b * x])
    p = np.linalg.lstsq(A * np.sqrt(w)[:, None], a * np.sqrt(w), rcond=None)[0]
    pred = A @ p
    var2 = a + (pred / np.maximum(b, 1e-9)) ** 2 * b
    x2p = float(np.sum((a - pred) ** 2 / np.maximum(var2, 1)))
    return {"na": float(a.sum()), "nb": float(b.sum()), "chi2": x2, "ndf": len(a) - 1, "p": float(C2.sf(x2, len(a) - 1)),
            "chi2_pol1": x2p, "ndf_pol1": len(a) - 2, "p_pol1": float(C2.sf(x2p, len(a) - 2)), "c": float(p[1] / p[0]) if p[0] else np.nan}


def _job(job):
    key, kind, args = job
    try:
        if kind == "nom":
            return key, fit_nom(*args)
        if kind == "tpl":
            return key, fit_tpl(*args[:6], pol1=args[6], nboot=args[7])
    except Exception as ex:  # noqa: BLE001
        return key, {"eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False, "determined": False, "reason": f"wyjątek: {ex}"}
    raise ValueError(kind)


# ----------------------------------------------------------------------------------------------------------- części analizy
def part_a(S, name, pmin):
    """Walidacja RS vs pary nominalne (tablica n_*)."""
    z = S.z; out = {}
    if not len(z["n_dm"]):
        return out
    npt, nm, ndm, nppt, found = z["n_pt"], z["n_m"], z["n_dm"], z["n_pis_pt"], z["n_found"]
    nomw = np.abs(nm - M0) < NOM_WIN
    for b, (lo, hi) in BINS.items():
        k = nomw & (npt >= lo) & (npt < hi) & (ndm < S.meta["row_dm"])
        hi_pt = k & (nppt >= pmin)
        d = (z["n_dm_own"] - ndm)[hi_pt & (found == 1)] * 1e3
        ours = (S.typ == 0) & S.nom & S.bin(b) & (S.pis >= pmin)
        inpk = np.abs(S.dm - 145.43) < 1.5
        r = {"n_nominal": int(k.sum()), "frac_nominal_pis_lt_pmin": float(np.mean(nppt[k] < pmin)) if k.any() else np.nan,
             "found_rs": float(np.mean(found[hi_pt] == 1)) if hi_pt.any() else np.nan,
             "found_ws_charge": float(np.mean(found[hi_pt] == 2)) if hi_pt.any() else np.nan,
             "d_dm_mean_MeV": float(d.mean()) if len(d) else np.nan, "d_dm_rms_MeV": float(d.std()) if len(d) else np.nan,
             "d_dm_frac_gt_0p5": float(np.mean(np.abs(d) > 0.5)) if len(d) else np.nan,
             "n_rs_ours": int(ours.sum()), "frac_rs_ours_nominal": float(np.mean((S.flag[ours] & 1) != 0)) if ours.any() else np.nan,
             "n_rs_ours_peak": int((ours & inpk).sum()),
             "frac_rs_ours_nominal_peak": float(np.mean((S.flag[ours & inpk] & 1) != 0)) if (ours & inpk).any() else np.nan,
             "frac_daughters_matched": float(np.mean(S.z["d_nmatch"][(S.z["d_pt"] >= lo) & (S.z["d_pt"] < hi)] == 2))}
        if name in ("fd_nd", "mb_nd"):
            sg = k & np.isin(z["n_truth"], [0, 1]) & np.isin(z["n_truth_ds"], [0, 1])
            r["mc_signal_nominal"] = int(sg.sum()); r["mc_signal_frac_pis_lt_pmin"] = float(np.mean(nppt[sg] < pmin)) if sg.any() else np.nan
            so = ours & S.sig
            r["mc_signal_ours"] = int(so.sum()); r["mc_signal_ours_frac_nominal"] = float(np.mean((S.flag[so] & 1) != 0)) if so.any() else np.nan
        out[b] = r
    return out


def part_b(S, thr, pmin):
    """MC: kształt szablonów vs prawdziwe tło RS; przeciek sygnału do WS; tło pikujące w RS."""
    out = {}
    for b in BINS:
        base = S.nom & S.bin(b) & (S.pis >= pmin)
        t80 = thr[b]["0.80"]
        cats = {"all": np.ones(len(S.dm), bool), "pass80": S.score >= t80, "fail80": S.score < t80}
        rb = {}
        for cn, cm in cats.items():
            k = base & cm
            hbkg = S.hist(k & (S.typ == 0) & ~S.sig)
            r = {"n_true_bkg": float(hbkg.sum()), "n_signal": float(S.hist(k & S.sig).sum())}
            for t in TPL:
                r[t] = shape_chi2(hbkg, S.hist(k & (S.typ == TYPES[t])))
            r["threshold_fn_on_true_bkg"] = fit_bkg_only(hbkg)
            rb[cn] = r
        k = base
        nsig = int((k & S.sig).sum())
        leak = k & (S.typ == 1) & ((S.flag & 2) != 0)
        pk_rs = k & (S.typ == 0) & ((S.flag & 2) != 0) & ~S.d0true
        inpk = np.abs(S.dm - 145.43) < 2.5
        rb["leak"] = {"n_signal": nsig, "n_ws_true_pion": int(leak.sum()), "ratio": float(leak.sum() / nsig) if nsig else np.nan,
                      "n_ws_true_pion_swap": int((leak & np.isin(S.truth, [100, 101])).sum()),
                      "ws_true_pion_in_peak": float(np.mean(inpk[leak])) if leak.any() else np.nan,
                      "n_rs_true_pion_d0_not_true": int(pk_rs.sum()), "ratio_rs_peaking": float(pk_rs.sum() / nsig) if nsig else np.nan,
                      "rs_peaking_in_peak": float(np.mean(inpk[pk_rs])) if pk_rs.any() else np.nan,
                      "n_rs_peak_window": int((k & (S.typ == 0) & inpk).sum()), "n_ws_peak_window": int((k & (S.typ == 1) & inpk).sum())}
        out[b] = rb
    return out


def build_jobs(tag, S, mask, shapes, thr, pmin, fast, nboot, wps=WPS, pol1=(False, True)):
    """Zadania fitów (nominał + szablony) dla zbioru (próbka, maska wierszy)."""
    jobs = []
    for b in BINS:
        base = mask & S.nom & S.bin(b) & (S.pis >= pmin)
        for wp in wps:
            if wp not in thr[b]:
                continue
            pas = S.score >= thr[b][wp]
            hp, hf = S.hist(base & (S.typ == 0) & pas), S.hist(base & (S.typ == 0) & ~pas)
            jobs.append(((tag, b, wp, "nominal"), "nom", (hp, hf, shapes[b], not fast)))
            for t in TPL:
                tp, tf = S.hist(base & (S.typ == TYPES[t]) & pas), S.hist(base & (S.typ == TYPES[t]) & ~pas)
                for p1 in pol1:
                    jobs.append(((tag, b, wp, t + ("+pol1" if p1 else "")), "tpl", (hp, hf, tp, tf, E, shapes[b], p1, nboot)))
    return jobs


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=os.path.join(HERE, "runs", "lxplus"))
    ap.add_argument("--out", default=None); ap.add_argument("--report", default=None)
    ap.add_argument("--bt-fd", default=os.path.join(CW, "B_trig", "runs", "lxplus", "bt_fd_nd.npz"))
    ap.add_argument("--pis-min", type=float, default=0.2, help="próg pT wolnego pionu [GeV] (Tracks: ≥ 0,2)")
    ap.add_argument("--pis-scan", default="0.25,0.30,0.40")
    ap.add_argument("--fast", action="store_true", help="nominał bez skanu profilu (HESSE)")
    ap.add_argument("--nboot", type=int, default=20); ap.add_argument("--nproc", type=int, default=6)
    ap.add_argument("--label", default="")
    a = ap.parse_args(argv)
    out_dir = a.out or os.path.join(a.inp, "out")
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    from charmtag import dstar_refit as DR
    S = {s: Sample(os.path.join(a.inp, f"dw_{s}.npz")) for s in SAMPLES if os.path.exists(os.path.join(a.inp, f"dw_{s}.npz"))}
    if "fd_nd" not in S:
        sys.exit("STOP: brak dw_fd_nd.npz (kształt sygnału, ε_MC)")
    thr, thr_src = thresholds(a.bt_fd, S["fd_nd"].meta.get("thresholds") or {})
    fd = S["fd_nd"]
    shapes, eps_mc = {}, {}
    for b in BINS:
        k = fd.sig & fd.nom & fd.bin(b) & (fd.pis >= a.pis_min)
        shapes[b] = DR.shape_from_mc(fd.dm[k])
        shapes[b]["n"] = int(k.sum())
        for wp in thr[b]:
            p = float(np.mean(fd.score[k] >= thr[b][wp])) if k.any() else np.nan
            eps_mc[(b, wp)] = (p, float(np.sqrt(p * (1 - p) / max(k.sum(), 1))))
    res = {"label": a.label, "input": a.inp, "pis_min": a.pis_min, "thresholds": thr, "thresholds_src": thr_src, "shapes": shapes,
           "eps_mc": {f"{b}@{wp}": v for (b, wp), v in eps_mc.items()}, "summary": {s: {k: v for k, v in S[s].meta.items() if k != "files"} for s in S}}
    res["a"] = {s: part_a(S[s], s, a.pis_min) for s in S}
    res["b"] = {s: part_b(S[s], thr, a.pis_min) for s in ("mb_nd", "fd_nd") if s in S}
    # (c) dane + (d) domknięcie MC
    jobs = []
    for run, (s, bm, bv) in RUNS.items():
        if s in S:
            jobs += build_jobs(run, S[s], (S[s].trig & bm) == bv, shapes, thr, a.pis_min, a.fast, a.nboot)
            for pm in [float(x) for x in a.pis_scan.split(",") if x]:
                jobs += [((f"{run}|pis{pm:.2f}",) + j[0][1:], j[1], j[2]) for j in
                         build_jobs(run, S[s], (S[s].trig & bm) == bv, shapes, thr, pm, True, 0, wps=["0.80"], pol1=(False,))
                         if j[0][1] == "3-5"]
    for s in ("mb_nd", "fd_nd"):
        if s in S:
            jobs += build_jobs("MC:" + s, S[s], np.ones(len(S[s].dm), bool), shapes, thr, a.pis_min, True, a.nboot, wps=["0.80"])
    print(f"fitów: {len(jobs)} (nproc {a.nproc})", flush=True)
    with Pool(a.nproc) as P:
        fits = dict(P.map(_job, jobs, chunksize=1))
    rows = []
    for (tag, b, wp, model), f in fits.items():
        e, ee = eps_mc.get((b, wp), (np.nan, np.nan))
        se = sym_err(f)
        sf = f["eps"] / e if e and np.isfinite(f.get("eps", np.nan)) else np.nan
        r = {"set": tag, "bin": b, "wp": wp, "model": model, "N_sig": f.get("N"), "N_err": float(np.hypot(f.get("N_err", np.nan), f.get("N_tpl", 0.0))),
             "N_err_fit": f.get("N_err"), "N_tpl": f.get("N_tpl", 0.0), "eps": f.get("eps"), "eps_err": se, "eps_mc": e, "SF": sf,
             "SF_err": float(np.hypot(se / e, sf * ee / e)) if e and np.isfinite(se) else np.nan, "determined": bool(f.get("determined")),
             "reason": f.get("reason"), "chi2": f.get("chi2"), "ndf": f.get("ndf"), "n_pass": f.get("n_pass"), "n_fail": f.get("n_fail"),
             "cp": f.get("cp"), "cf": f.get("cf")}
        r["measurable"] = bool(r["determined"] and np.isfinite(r["SF_err"]) and r["SF_err"] < 0.15)
        rows.append(r)
    res["rows"] = rows
    # prawda MC do domknięcia
    res["d"] = {}
    for s in ("mb_nd", "fd_nd"):
        if s not in S:
            continue
        X = S[s]; res["d"][s] = {}
        for b in BINS:
            k = X.nom & X.bin(b) & (X.pis >= a.pis_min) & (X.dm >= LO) & (X.dm < HI)
            sg = k & X.sig
            pk = k & (X.typ == 0) & ((X.flag & 2) != 0) & ~X.d0true
            res["d"][s][b] = {"N_true": int(sg.sum()), "eps_true": float(np.mean(X.score[sg] >= thr[b]["0.80"])) if sg.any() else np.nan,
                              "N_peaking_other": int(pk.sum())}
    # (c4) funkcja progowa na szablonach w danych
    res["c4"] = {}
    for run, (s, bm, bv) in RUNS.items():
        if s not in S:
            continue
        X = S[s]; mk = (X.trig & bm) == bv
        for b in BINS:
            base = mk & X.nom & X.bin(b) & (X.pis >= a.pis_min)
            for t in TPL:
                res["c4"][f"{run}|{b}|{t}"] = fit_bkg_only(X.hist(base & (X.typ == TYPES[t])))
            for t in ("ROT", "MIX"):
                res["c4"][f"{run}|{b}|{t}_vs_WS"] = shape_chi2(X.hist(base & (X.typ == 1)), X.hist(base & (X.typ == TYPES[t])))
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(out_dir, "dstws.json"), "w"), indent=1, default=lambda o: None if isinstance(o, float) and not np.isfinite(o) else o)
    rep = a.report or os.path.join(HERE, "RAPORT.md")
    open(rep, "w").write(report(res))
    print(f"zapisano {os.path.join(out_dir, 'dstws.json')} i {rep} ({res['seconds']} s)")


# ----------------------------------------------------------------------------------------------------------- raport
def _f(x, n=2):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


def _row(rows, s, b, wp, model):
    for r in rows:
        if r["set"] == s and r["bin"] == b and r["wp"] == wp and r["model"] == model:
            return r
    return None


def verdict_b(x):
    if not x or not x.get("ndf"):
        return "brak stat."
    if x["chi2"] / x["ndf"] < 2 and x["p"] > 0.01:
        return "TAK"
    if x["chi2_pol1"] / x["ndf_pol1"] < 2 and x["p_pol1"] > 0.01:
        return "po pol1"
    return "NIE"


def report(R):
    L = []
    w = L.append
    partial = any(v.get("partial") for v in R["summary"].values())
    w("# DSTWS — tło D* przez konstrukcję: D0 × tor z drzewa Tracks (RS / WS / ROT / MIX) jako szablon tła w Δm\n")
    w(f"Wygenerowano {time.strftime('%Y-%m-%d %H:%M')} przez `DSTWS/analyze.py` (projekt B, pojedyncze D, nieślepe). Wejście `{R['input']}`"
      f"{' — **PRÓBKA TESTOWA (części plików), liczby tylko do sprawdzenia kodu**' if partial else ''}. {R.get('label', '')}\n")
    w("## 0. Konstrukcja i kryteria (zapisane przed danymi — docstring `analyze.py`)\n")
    w("- D0: wiersze D0 h5 v2.0, L_xy ≥ 0,25 mm, |m − m0| < 75 MeV (analiza: ±55 MeV), score H modelem nominalnym — selekcja D0 par B_trig, ale wszystkie D0, "
      "nie tylko z partnerem D* w ntuplu. Wolny pion: dowolny tor z drzewa `Tracks` (pT > 0,2 GeV; **bez d0/z0 i jakości toru**), niebędący córką tego D0.")
    w("- Δm = m(D0 + π) − m(D0) z 4-pędu kandydata D0 i toru w PV (nominał: ponowne dopasowanie D*). Typy: RS (q_π = −q_K), WS (q_π = q_K), "
      "ROT (RS obrócony w φ o k·60°, k = 1…5), MIX (tor RS-ładunku z 5 zdarzeń o najbliższej krotności torów w tej samej porcji pliku ≈ ten sam LB).")
    w("- Różnica względem KF (`KF_kontrfakt`: „realny szablon WS gorzej, rotacja k≈9: 0…+13 %”, S6 candrot χ²/ndf 8–127): tam obracano/mieszano tory "
      "WCHODZĄCE do wierzchołka istniejącego kandydata (wierzchołek przestaje pasować do kinematyki) albo robiono pseudo-WS przez zamianę K↔π, i celem był "
      "kształt masy m(D) pojedynczych D / par. Tu wierzchołek D0 jest nietknięty; zmieniany jest tylko wolny pion, który w ogóle nie wchodzi do wierzchołka, "
      "a zmienną jest Δm.")
    w("- Kryteria: K-a (odtworzenie nominału ≥ 99 %, |⟨δΔm⟩| < 0,05 MeV, RMS < 0,30 MeV), K-b (χ²/ndf < 2 i p > 0,01 kształtu szablonu vs prawdziwe tło RS w MC; "
      "przeciek WS < 5 %), K-c (bin mierzalny: fit określony i σ(SF) < 0,15; zysk σ(N) < 0,8 nominału; stabilność WS/ROT/MIX; funkcja progowa χ²/ndf < 2), "
      "K-d (domknięcie MC: |ΔN| < 2σ i < 10 %).")
    w(f"- Próg pT(π_s) ≥ {R['pis_min']} GeV. Progi WP: {R['thresholds_src']}. Bin 2–3 GeV: **brak** — h5 v2.0 i ta paczka obejmują tylko produkcję dpt3to150 "
      "(dpt2to3 wymaga osobnej konwersji h5 i modelu H dla pT < 3).\n")
    w("| próbka | plików | D0 w oknie ±75 MeV | kombinacji (Δm < 160) | RS | WS | ROT | MIX | par nominalnych | córki D0 dopasowane w Tracks |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for s, m in R["summary"].items():
        nt = m["n_by_type"]
        w(f"| {s} | {m['n_files']} | {m['n_d0_window']} | {m['n_comb']} | {nt['RS']} | {nt['WS']} | {nt['ROT']} | {nt['MIX']} | {m['n_nominal_pairs']} | "
          f"{_f(m.get('frac_daughters_matched'), 4)} |")
    w("\n## (a) Walidacja: nasze RS vs nominalne pary D0←D* z ntupla\n")
    w("| próbka | bin | par nom. (Δm < 160) | z pT(π_s) < próg | odnalezione w RS (pT ≥ próg) | ⟨δΔm⟩ [MeV] | RMS [MeV] | nasze RS | z nich nominalne | w piku ±1,5 MeV: nominalne | K-a1 | K-a2 |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s, A in R["a"].items():
        for b, r in A.items():
            k1 = "TAK" if r["found_rs"] is not None and r["found_rs"] >= 0.99 else "NIE"
            k2 = "TAK" if abs(r["d_dm_mean_MeV"] or 9) < 0.05 and (r["d_dm_rms_MeV"] or 9) < 0.30 else "NIE"
            w(f"| {s} | {b} | {r['n_nominal']} | {_f(r['frac_nominal_pis_lt_pmin'], 3)} | {_f(r['found_rs'], 4)} | {_f(r['d_dm_mean_MeV'], 3)} | "
              f"{_f(r['d_dm_rms_MeV'], 3)} | {r['n_rs_ours']} | {_f(r['frac_rs_ours_nominal'], 3)} | {_f(r['frac_rs_ours_nominal_peak'], 3)} | {k1} | {k2} |")
    w("\nMC, sygnał prawdziwy (D0 ∈ {0, 1} i pion z D*):\n")
    w("| próbka | bin | sygnał w parach nominalnych | z pT(π_s) < próg (strata Tracks) | sygnał w naszych RS | z nich w nominale |")
    w("|---|---|---|---|---|---|")
    for s, A in R["a"].items():
        for b, r in A.items():
            if "mc_signal_nominal" in r:
                w(f"| {s} | {b} | {r['mc_signal_nominal']} | {_f(r['mc_signal_frac_pis_lt_pmin'], 3)} | {r['mc_signal_ours']} | {_f(r['mc_signal_ours_frac_nominal'], 3)} |")
    w("\nK-a3 (N_sig z fitu nominalnego na naszych RS) — tabela (c), kolumna „nominał”; porównanie z B_trig `RAPORT_DANE.md` (tam pT(π_s) ≥ 0,1 GeV).\n")
    w("## (b) MC: czy szablon ma kształt prawdziwego tła RS? (χ²/ndf, biny 0,5 MeV, [140,47; 160] MeV)\n")
    w("| próbka | bin | kategoria | tło RS prawda | WS: χ²/ndf (p) → po pol1; c | ROT | MIX | funkcja progowa na tle prawdziwym χ²/ndf (p) | WS | ROT | MIX |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for s, B in R["b"].items():
        for b, rb in B.items():
            for cn in ("all", "pass80", "fail80"):
                r = rb[cn]
                cell = lambda x: (f"{_f(x['chi2'], 1)}/{x['ndf']} ({_f(x['p'], 3)}) → {_f(x['chi2_pol1'], 1)}/{x.get('ndf_pol1', 0)}; {_f(x['c'], 3)}"  # noqa: E731
                                  if x.get("ndf") else "–")
                tf = r["threshold_fn_on_true_bkg"]
                w(f"| {s} | {b} | {cn} | {int(r['n_true_bkg'])} | {cell(r['WS'])} | {cell(r['ROT'])} | {cell(r['MIX'])} | "
                  f"{_f(tf['chi2'], 1)}/{tf['ndf']} ({_f(tf['p'], 3)}) | {verdict_b(r['WS'])} | {verdict_b(r['ROT'])} | {verdict_b(r['MIX'])} |")
    w("\nPrzeciek sygnału i tło pikujące (okno ±55 MeV, cały zakres Δm < 160):\n")
    w("| próbka | bin | sygnał RS | WS z prawdziwym pionem D* | przeciek WS/sygnał | z tego bliźniak K↔π | w piku ±2,5 MeV | RS: pion z D*, D0 nie-sygnał | / sygnał | K-b3 |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for s, B in R["b"].items():
        for b, rb in B.items():
            r = rb["leak"]
            w(f"| {s} | {b} | {r['n_signal']} | {r['n_ws_true_pion']} | {_f(r['ratio'], 3)} | {r['n_ws_true_pion_swap']} | {_f(r['ws_true_pion_in_peak'], 2)} | "
              f"{r['n_rs_true_pion_d0_not_true']} | {_f(r['ratio_rs_peaking'], 3)} | {'TAK' if (r['ratio'] or 9) < 0.05 else 'NIE'} |")
    w("\nUwaga: „prawdziwy pion D*” = tor dopasowany do wolnego pionu dowolnego nominalnego kandydata D* z prawdą |211| ← |413|; piony D*, które nie mają "
      "żadnego kandydata D* w ntuplu, są w „tle prawdziwym” (resztkowy pik w RS-tło zawyża χ² szablonów).\n")
    w("## (c) Dane: fit Δm pass/fail — nominał (próg × pol2) vs szablon tła\n")
    rows = R["rows"]
    for run in RUNS:
        if not any(r["set"] == run for r in rows):
            continue
        for wp in WPS:
            w(f"\n### run {run}, WP {wp}\n")
            w("| bin | model tła | N_sig ± σ (w tym stat. szablonu) | σ(N)/σ(N)_nom | ε_dane | SF ± σ | χ²/ndf | określony | mierzalny (K-c1) |")
            w("|---|---|---|---|---|---|---|---|---|")
            for b in BINS:
                n0 = _row(rows, run, b, wp, "nominal")
                for model in ["nominal"] + [t + p for t in TPL for p in ("", "+pol1")]:
                    r = _row(rows, run, b, wp, model)
                    if r is None:
                        continue
                    ratio = r["N_err"] / n0["N_err"] if n0 and n0["N_err"] and np.isfinite(n0["N_err"]) and r["N_err"] is not None else np.nan
                    w(f"| {b} | {model} | {_f(r['N_sig'], 0)} ± {_f(r['N_err'], 0)} ({_f(r['N_tpl'], 0)}) | {_f(ratio, 2)} | {_f(r['eps'], 3)} ± {_f(r['eps_err'], 3)} | "
                      f"{_f(r['SF'], 2)} ± {_f(r['SF_err'], 2)} | {_f(r['chi2'], 0)}/{r['ndf']} | {'tak' if r['determined'] else 'NIE (' + str(r['reason']) + ')'} | "
                      f"{'**TAK**' if r['measurable'] else 'nie'} |")
    w("\n### Stabilność na wybór szablonu (K-c3) i werdykt binu 3–5 GeV\n")
    w("| run | WP | bin | SF nominał | SF WS | SF ROT | SF MIX | max ΔSF (szablony) | min σ_stat | stabilny | mierzalny: nominał / WS / ROT / MIX |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for run in RUNS:
        for wp in WPS:
            for b in BINS:
                rr = {mn: _row(rows, run, b, wp, mn) for mn in ["nominal"] + TPL}
                if rr["nominal"] is None:
                    continue
                sfs = [rr[t]["SF"] for t in TPL if rr[t] and rr[t]["SF"] is not None and np.isfinite(rr[t]["SF"])]
                ers = [rr[t]["SF_err"] for t in TPL if rr[t] and rr[t]["SF_err"] is not None and np.isfinite(rr[t]["SF_err"])]
                d = max(sfs) - min(sfs) if len(sfs) > 1 else np.nan
                c = lambda r: f"{_f(r['SF'], 2)} ± {_f(r['SF_err'], 2)}" if r else "–"  # noqa: E731
                w(f"| {run} | {wp} | {b} | {c(rr['nominal'])} | {c(rr['WS'])} | {c(rr['ROT'])} | {c(rr['MIX'])} | {_f(d, 2)} | {_f(min(ers) if ers else np.nan, 2)} | "
                  f"{'tak' if ers and np.isfinite(d) and d < min(ers) else 'NIE'} | " + " / ".join('TAK' if rr[mn] and rr[mn]['measurable'] else 'nie' for mn in ['nominal'] + TPL) + " |")
    scan = [r for r in rows if "|pis" in r["set"]]
    if scan:
        w("\n### Skan progu pT(π_s), bin 3–5 GeV, WP 0.80 (nominał HESSE; szablony bez stat. szablonu)\n")
        w("| run | próg | model | N_sig ± σ | SF ± σ | określony |")
        w("|---|---|---|---|---|---|")
        for r in sorted(scan, key=lambda r: (r["set"], r["model"])):
            run, p = r["set"].split("|pis")
            w(f"| {run} | {p} | {r['model']} | {_f(r['N_sig'], 0)} ± {_f(r['N_err'], 0)} | {_f(r['SF'], 2)} ± {_f(r['SF_err'], 2)} | {'tak' if r['determined'] else 'NIE'} |")
    w("\n### Funkcja progowa nominału na szablonach w danych (K-c4) i zgodność szablonów między sobą\n")
    w("| run | bin | fit próg × pol2 do WS: χ²/ndf (p) | do ROT | do MIX | kształt ROT vs WS χ²/ndf (p) | MIX vs WS | K-c4 |")
    w("|---|---|---|---|---|---|---|---|")
    c4 = R["c4"]
    for run in RUNS:
        for b in BINS:
            if f"{run}|{b}|WS" not in c4:
                continue
            f1 = lambda x: f"{_f(x['chi2'], 1)}/{x['ndf']} ({_f(x['p'], 3)})" if x.get("ndf") else "–"  # noqa: E731
            x = c4[f"{run}|{b}|WS"]
            ok = "TAK" if x.get("ndf") and x["chi2"] / x["ndf"] < 2 and x["p"] > 0.01 else ("NIE" if x.get("ndf") else "–")
            w(f"| {run} | {b} | {f1(x)} | {f1(c4[f'{run}|{b}|ROT'])} | {f1(c4[f'{run}|{b}|MIX'])} | {f1(c4[f'{run}|{b}|ROT_vs_WS'])} | {f1(c4[f'{run}|{b}|MIX_vs_WS'])} | {ok} |")
    w("\n## (d) Domknięcie na MC (WP 0.80): N_sig i ε z fitu vs prawda\n")
    w("| próbka | bin | N prawda (RS sygnał w zakresie fitu) | inne pikujące RS | model | N_fit ± σ | (N_fit − N)/σ | N_fit/N | ε_fit ± σ | ε prawda | K-d1 |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for s, Dd in R["d"].items():
        for b, t in Dd.items():
            for model in ["nominal"] + [x + p for x in TPL for p in ("", "+pol1")]:
                r = _row(rows, "MC:" + s, b, "0.80", model)
                if r is None or not t["N_true"]:
                    continue
                pull = (r["N_sig"] - t["N_true"]) / r["N_err"] if r["N_err"] else np.nan
                ratio = r["N_sig"] / t["N_true"]
                ok = np.isfinite(pull) and abs(pull) < 2 and abs(ratio - 1) < 0.10 and r["eps_err"] is not None and np.isfinite(r["eps_err"]) \
                    and abs(r["eps"] - t["eps_true"]) < 2 * r["eps_err"]
                w(f"| {s} | {b} | {t['N_true']} | {t['N_peaking_other']} | {model} | {_f(r['N_sig'], 0)} ± {_f(r['N_err'], 0)} | {_f(pull, 2)} | {_f(ratio, 3)} | "
                  f"{_f(r['eps'], 3)} ± {_f(r['eps_err'], 3)} | {_f(t['eps_true'], 3)} | {'TAK' if ok else 'NIE'} |")
    w("\n## Zastrzeżenia\n")
    w("- Tracks: pT > 0,2 GeV, brak d0/z0/jakości → RS luźniejsze niż nominał w selekcji toru, a ciaśniejsze w pT (nominał: pT(π_s) > 0,1 GeV) — "
      "strata sygnału w binie 3–5 GeV w tabeli (a).")
    w("- WS zawiera przeciek prawdziwych D* przez bliźniaka K↔π (tabela (b)) — pik w szablonie odejmuje sygnał; ROT i MIX go nie mają.")
    w("- ε_MC i kształt sygnału z fd_nd z NASZEJ definicji Δm i sygnału (D0 prawda ∈ {0, 1}, pion dopasowany do prawdziwego wolnego pionu D*).")
    w("- Błąd szablonu: rozrzut ε i N przy fluktuacji Poissona szablonu (bootstrap), dodany w kwadraturze; szablon wygładzany w grupach 0,5 MeV.")
    w("- MC bez wag (jak B_trig bez emulacji triggera). Klasy triggera 428770 (J12only/MBTSbased) nie są rozdzielane w tej analizie (bity są w `d_trig`).")
    if partial:
        w("- **Wynik testowy na 1–2 plikach/próbkę** — werdykty K-* nie są wynikiem analizy.")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
