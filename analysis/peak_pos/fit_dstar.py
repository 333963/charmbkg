"""PEAK_pos: D*±-tag (pełna lumi 428770, HYP_d0/runs/lxplus_dstar). (1) Δm D*+ vs D*− (dane vs MC),
(2) m(Kπ) D0 vs D̄0 z tagu (zapach z π_s, przypisanie K/π z tagu → bez refleksji). WP70 modelu D0_H.
Wynik: out/dstar_rows.csv"""
import numpy as np, pandas as pd
from iminuit import Minuit
from common import THR, PTB, M0, OUT, pass_wp
from fitlib import fit_mc, fit_data, dg_binfrac, seff

MPI = 0.13957039
DM0 = 0.14543


def load(tag, wp=0.7, mwin=0.04):
    L = pd.read_parquet(f"{OUT}/dstar_{tag}.parquet")
    ok = (L.hypo_from_tag == 0) & (L.dstar_conflict == 0) & pass_wp(L.score.to_numpy(), L.pt.to_numpy(), THR["D0"], wp)
    if tag == "fd_nd":
        ok &= np.isin(L.truth, [0, 1]) & (L.hypo_swap == 0)
    L = L[ok].reset_index(drop=True)
    L["q"] = L.soft_pi_charge.astype(int)
    return L


def fit_dm(x, lo=0.1396, hi=0.160, bw=0.0002, shape=None, nb=1):
    """Δm: 2 Gaussy (wolny kształt gdy shape=None, inaczej μ, k wolne) + tło u^a·(1+c1·u[+c2·u²]), u = Δm − m_π."""
    e = np.arange(lo, hi + 1e-12, bw); n, _ = np.histogram(x, e); c = 0.5 * (e[1:] + e[:-1]); u = (c - MPI) * 1e3

    def bk(a, cs):
        b = np.power(np.maximum(u, 1e-6), a) * (1 + sum(ci * u ** (i + 1) for i, ci in enumerate(cs)))
        return b
    if shape is None:
        names = ["Ns", "mu", "s1", "s2", "f", "Nb", "a"] + [f"c{i}" for i in range(nb)]

        def nll(*p):
            Ns, mu, s1, s2, f, Nb, a = p[:7]; b = bk(a, p[7:])
            if np.any(b <= 0): return 1e12
            lam = Ns * dg_binfrac(e, mu, s1, s2, f) + Nb * b / b.sum() + 1e-12
            return float(np.sum(lam - n * np.log(lam)))
        st = [n.sum() * 0.3, DM0, 0.0004, 0.001, 0.6, n.sum() * 0.7, 0.5] + [0.0] * nb
    else:
        names = ["Ns", "mu", "k", "Nb", "a"] + [f"c{i}" for i in range(nb)]

        def nll(*p):
            Ns, mu, k, Nb, a = p[:5]; b = bk(a, p[5:])
            if np.any(b <= 0): return 1e12
            lam = Ns * dg_binfrac(e, mu, k * shape["s1"], k * shape["s2"], shape["f"]) + Nb * b / b.sum() + 1e-12
            if np.any(lam <= 0): return 1e12
            return float(np.sum(lam - n * np.log(lam)))
        st = [n.sum() * 0.3, DM0, 1.0, n.sum() * 0.7, 0.5] + [0.0] * nb
    M = Minuit(nll, *st, name=names); M.errordef = 0.5
    M.limits["Ns"] = (0, n.sum()); M.limits["Nb"] = (0, 2 * n.sum()); M.limits["mu"] = (DM0 - 0.002, DM0 + 0.002); M.limits["a"] = (0, 3)
    if shape is None:
        M.limits["s1"] = (0.0001, 0.003); M.limits["s2"] = (0.0002, 0.006); M.limits["f"] = (0.05, 0.98)
    else:
        M.limits["k"] = (0.4, 2.5)
    M.migrad(); M.migrad(); M.hesse()
    v, er = M.values.to_dict(), M.errors.to_dict()
    r = {**{k: float(x) for k, x in v.items()}, **{"e_" + k: float(x) for k, x in er.items()}, "ok": bool(M.valid), "n": int(n.sum())}
    if shape is None:
        if r["s1"] > r["s2"]:
            r["s1"], r["s2"], r["f"] = r["s2"], r["s1"], 1 - r["f"]
        r["seff"] = seff(r["s1"], r["s2"], r["f"])
    else:
        r["seff"] = r["k"] * seff(shape["s1"], shape["s2"], shape["f"]); r["e_seff"] = r["e_k"] * seff(shape["s1"], shape["s2"], shape["f"])
    return r


rows = []
BINS = [(f"pT{lo}-{hi}", lambda L, lo=lo, hi=hi: (L.pt >= lo) & (L.pt < hi)) for lo, hi in PTB] + \
       [("pT3-30", lambda L: np.ones(len(L), bool)), ("abseta<1", lambda L: np.abs(L.eta) < 1), ("abseta>1", lambda L: np.abs(L.eta) >= 1),
        ("eta<-1", lambda L: L.eta < -1), ("eta-1..0", lambda L: (L.eta >= -1) & (L.eta < 0)), ("eta0..1", lambda L: (L.eta >= 0) & (L.eta < 1)),
        ("eta>1", lambda L: L.eta >= 1)]
for wp in (0.7, 0.5):
    D, S = load("data_428770", wp), load("fd_nd", wp)
    for bname, sel in BINS:
        sb, db = np.asarray(sel(S)), np.asarray(sel(D))
        # (1) Δm, okno m(Kπ) ±40 MeV (wariant ±60 MeV)
        for mw in (0.04, 0.06):
            ms, md = np.abs(S.m - M0["D0"]).to_numpy() < mw, np.abs(D.m - M0["D0"]).to_numpy() < mw
            shp = fit_dm(S.dm.to_numpy()[sb & ms])
            for q in (0, 1, -1):
                qs = np.ones(len(S), bool) if q == 0 else (S.q == q).to_numpy()
                qd = np.ones(len(D), bool) if q == 0 else (D.q == q).to_numpy()
                mc = fit_dm(S.dm.to_numpy()[sb & ms & qs])
                for var, nb, hi in (("nom", 1, 0.160), ("bkg2", 2, 0.160), ("range165", 1, 0.165)):
                    if (var != "nom" or mw != 0.04) and not bname.startswith("pT"):
                        continue
                    r = fit_dm(D.dm.to_numpy()[db & md & qd], hi=hi, shape=shp, nb=nb)
                    rows.append(dict(obs="dm", wp=wp, bin=bname, q=q, var=var if mw == 0.04 else f"mwin{int(mw*1e3)}_{var}",
                                     mu=r["mu"], e_mu=r["e_mu"], k=r["k"], e_k=r["e_k"], seff=r["seff"], e_seff=r["e_seff"],
                                     Ns=r["Ns"], e_Ns=r["e_Ns"], ok=r["ok"], mu_mc=mc["mu"], e_mu_mc=mc["e_mu"], seff_mc=mc["seff"],
                                     n_mc=int((sb & ms & qs).sum()), shape_seff=shp["seff"]))
        # (2) m(Kπ) z tagu, |Δm − 145,43| < 1,2 MeV (wariant 0,8 MeV)
        for dw in (0.0012, 0.0008):
            ds_, dd_ = np.abs(S.dm - DM0).to_numpy() < dw, np.abs(D.dm - DM0).to_numpy() < dw
            shp = fit_mc(S.m.to_numpy()[sb & ds_], 1.76, 1.98)
            for q in (0, 1, -1):
                qs = np.ones(len(S), bool) if q == 0 else (S.q == q).to_numpy()
                qd = np.ones(len(D), bool) if q == 0 else (D.q == q).to_numpy()
                mc = fit_mc(S.m.to_numpy()[sb & ds_ & qs], 1.76, 1.98)
                for var, bkg, (lo, hi) in (("nom", "cheb2", (1.76, 1.98)), ("cheb3", "cheb3", (1.76, 1.98)),
                                           ("exp1", "exp1", (1.76, 1.98)), ("narrow", "cheb2", (1.79, 1.94))):
                    if (var != "nom" or dw != 0.0012) and not bname.startswith("pT"):
                        continue
                    r = fit_data(D.m.to_numpy()[db & dd_ & qd], lo, hi, shp, bkg=bkg, mu0=M0["D0"])
                    rows.append(dict(obs="mKpi_tag", wp=wp, bin=bname, q=q, var=var if dw == 0.0012 else f"dmwin08_{var}",
                                     mu=r["mu"], e_mu=r["e_mu"], k=r["k"], e_k=r["e_k"], seff=r["seff"], e_seff=r["e_seff"],
                                     Ns=r["Ns"], e_Ns=r["e_Ns"], ok=r["ok"], mu_mc=mc["mu"], e_mu_mc=mc["e_mu"], seff_mc=mc["seff"],
                                     n_mc=int((sb & ds_ & qs).sum()), shape_seff=shp["seff"]))
        print(wp, bname, flush=True)
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/dstar_rows.csv", index=False)
T["dmu"] = (T.mu - T.mu_mc) * 1e3; T["e_dmu"] = np.hypot(T.e_mu, T.e_mu_mc) * 1e3; T["rs"] = T.seff / T.seff_mc
print(T[T["var"] == "nom"][["obs", "wp", "bin", "q", "dmu", "e_dmu", "rs", "e_k", "Ns", "e_Ns", "ok", "n_mc"]].round(4).to_string())
