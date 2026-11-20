"""PEAK_pos: tabele do RAPORT.md (out/tables.md) + liczby zbiorcze (out/summary.json)."""
import json
import numpy as np, pandas as pd
from common import load, OUT
from fitlib import fit_mc

P = pd.read_csv(f"{OUT}/peaks_rows.csv"); P["obs"] = np.where(P.sp == "D0", np.where(P.veto, "D0_veto", "D0_keepboth"), "Dplus")
Dd = pd.read_csv(f"{OUT}/dstar_rows.csv")
A = pd.concat([P, Dd], ignore_index=True)
SAG = json.load(open(f"{OUT}/sagitta.json"))
BINS = ["pT3-5", "pT5-8", "pT8-30", "pT3-30", "abseta<1", "abseta>1"]
out, summ = [], {}


def g(T, obs, wp, b, q, var="nom"):
    x = T[(T.obs == obs) & (T.wp == wp) & (T.bin == b) & (T.q == q) & (T["var"] == var)]
    return None if len(x) == 0 else x.iloc[0]


def wavg(v, e):
    v, e = np.array(v), np.array(e); w = 1 / e**2
    m = np.sum(w * v) / np.sum(w); return float(m), float(1 / np.sqrt(np.sum(w))), float(np.sum(w * (v - m) ** 2))


for obs, unit in [("Dplus", "D⁺→Kππ, m(Kππ)"), ("D0_veto", "D⁰→Kπ inkluzywnie, weto refl."), ("D0_keepboth", "D⁰ keep_both + szablon refl. z MC"),
                  ("mKpi_tag", "D⁰ z tagu D* (π_s → zapach, bez refl.)"), ("dm", "Δm(D*) = m(Kππ_s) − m(Kπ)")]:
    for wp in (0.7, 0.5):
        rows = []
        for b in BINS:
            r = {q: g(A, obs, wp, b, q) for q in (1, -1, 0)}
            if r[0] is None:
                continue
            f = lambda x: (x.mu - x.mu_mc) * 1e3, lambda x: np.hypot(x.e_mu, x.e_mu_mc) * 1e3
            dp, ep = f[0](r[1]), f[1](r[1]); dn, en = f[0](r[-1]), f[1](r[-1]); d0, e0 = f[0](r[0]), f[1](r[0])
            asy_d = (r[1].mu - r[-1].mu) * 1e3; e_asy_d = np.hypot(r[1].e_mu, r[-1].e_mu) * 1e3
            asy_mc = (r[1].mu_mc - r[-1].mu_mc) * 1e3; e_asy_mc = np.hypot(r[1].e_mu_mc, r[-1].e_mu_mc) * 1e3
            rs = lambda x: (x.seff / x.seff_mc, x.e_seff / x.seff_mc)
            vs = [g(A, obs, wp, b, 0, v) for v in ("cheb3", "exp1", "narrow", "bw1", "bkg2", "range165", "mwin60_nom", "dmwin08_nom")]
            vs = [v for v in vs if v is not None]
            spread_mu = max([abs(v.mu - r[0].mu) * 1e3 for v in vs], default=np.nan)
            spread_s = max([abs(v.seff / r[0].seff - 1) for v in vs], default=np.nan)
            rows.append(dict(bin=b, dmu_D=dp, e_D=ep, dmu_Dbar=dn, e_Dbar=en, dmu_all=d0, e_all=e0, asy_data=asy_d, e_asy_data=e_asy_d,
                             asy_mc=asy_mc, e_asy_mc=e_asy_mc, rs_D=rs(r[1])[0], e_rs_D=rs(r[1])[1], rs_Dbar=rs(r[-1])[0],
                             e_rs_Dbar=rs(r[-1])[1], rs_all=rs(r[0])[0], e_rs_all=rs(r[0])[1], spread_mu=spread_mu, spread_s=spread_s,
                             Ns=r[0].Ns, e_Ns=r[0].e_Ns, mu_mc=r[0].mu_mc * 1e3, seff_mc=r[0].seff_mc * 1e3,
                             fitok=bool(r[1].ok and r[-1].ok and r[0].ok and 0.41 < r[0].k < 2.49 and 0.41 < r[1].k < 2.49 and 0.41 < r[-1].k < 2.49)))
        R = pd.DataFrame(rows)
        if len(R) == 0:
            continue
        sc = 1e3 if obs == "dm" else 1  # Δm w keV
        u = "keV" if obs == "dm" else "MeV"
        out.append(f"\n### {unit} — WP{int(wp*100)}  (Δμ = μ_dane − μ_MC [{u}]; σ = σ_eff 2 Gaussów)\n")
        out.append(f"| bin | N_S dane | μ_MC [{u}] | σ_MC [{u}] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |")
        out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for x in R.itertuples():
            flag = "" if x.fitok else " ⚠"
            out.append(f"| {x.bin}{flag} | {x.Ns:.0f}±{x.e_Ns:.0f} | {x.mu_mc*sc:.2f} | {x.seff_mc*sc:.2f} | {x.dmu_D*sc:+.2f}±{x.e_D*sc:.2f} | "
                       f"{x.dmu_Dbar*sc:+.2f}±{x.e_Dbar*sc:.2f} | {x.dmu_all*sc:+.2f}±{x.e_all*sc:.2f} | {x.asy_data*sc:+.2f}±{x.e_asy_data*sc:.2f} | "
                       f"{x.asy_mc*sc:+.2f}±{x.e_asy_mc*sc:.2f} | {x.rs_D:.2f}±{x.e_rs_D:.2f} | {x.rs_Dbar:.2f}±{x.e_rs_Dbar:.2f} | "
                       f"{x.rs_all:.2f}±{x.e_rs_all:.2f} | {x.spread_mu*sc:.2f} / {100*x.spread_s:.0f} % |")
        pt = R[R.bin.isin(["pT3-5", "pT5-8", "pT8-30"]) & R.fitok]
        if len(pt):
            m, e, chi = wavg(pt.dmu_all, pt.e_all)
            ma, ea, _ = wavg(pt.asy_data - pt.asy_mc, np.hypot(pt.e_asy_data, pt.e_asy_mc))
            ms, es, chis = wavg(pt.rs_all, pt.e_rs_all)
            summ[f"{obs}_WP{int(wp*100)}"] = dict(dmu=m * sc, e_dmu=e * sc, chi2_dmu=chi, ndf=len(pt) - 1, asy=ma * sc, e_asy=ea * sc,
                                                  rs=ms, e_rs=es, chi2_rs=chis, bins_used=list(pt.bin), unit=u,
                                                  spread_mu_max=float(np.nanmax(pt.spread_mu)) * sc, spread_s_max=float(np.nanmax(pt.spread_s)))
            out.append(f"\nŚrednia ważona binów pT (tylko fity bez ⚠): Δμ = {m*sc:+.2f} ± {e*sc:.2f} {u} (χ²={chi:.1f}/{len(pt)-1}); "
                       f"(D−D̄)_dane − (D−D̄)_MC = {ma*sc:+.2f} ± {ea*sc:.2f} {u}; σ_d/σ_MC = {ms:.3f} ± {es:.3f} (χ²={chis:.1f}/{len(pt)-1}).")
# sagitta: D+ (czuły) i D0
for obs, sp in (("Dplus", "Dplus"), ("D0_veto", "D0"), ("mKpi_tag", "D0")):
    rr = []
    for b in BINS:
        r1, r2 = g(A, obs, 0.7, b, 1), g(A, obs, 0.7, b, -1)
        if r1 is None:
            continue
        asy = ((r1.mu - r2.mu) - (r1.mu_mc - r2.mu_mc)) * 1e3
        ea = np.hypot(np.hypot(r1.e_mu, r2.e_mu), np.hypot(r1.e_mu_mc, r2.e_mu_mc)) * 1e3
        k = SAG[sp][b]["dm_ddelta_MeV_per_TeVinv"]
        rr.append(dict(bin=b, asy=asy, e=ea, sens=2 * k, delta=asy / (2 * k), e_delta=ea / abs(2 * k)))
    summ[f"sagitta_{obs}"] = rr
json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1)
# MC D+: model FEAT (fd) vs model analizy (mb_nd_dp prawda)
mc_cmp = []
for nm in ("mc", "mbmc"):
    S = load("Dplus", nm)
    for lo, hi in [(3, 5), (5, 8), (8, 30)]:
        s = S[(S.pt >= lo) & (S.pt < hi)]
        r = fit_mc(s.m.to_numpy(), 1.74, 2.04)
        mc_cmp.append(f"| {nm} | {lo}-{hi} | {len(s)} | {r['mu']*1e3:.2f}±{r['e_mu']*1e3:.2f} | {r['seff']*1e3:.2f} |")
out.append("\n### Kontrola MC D⁺: fd (model FEAT H, własne progi) vs mb_nd (model analizy Dplus_H), WP70\n")
out.append("| próbka | pT | N | μ [MeV] | σ_eff [MeV] |\n|---|---|---|---|---|")
out += mc_cmp
open(f"{OUT}/tables.md", "w").write("\n".join(out) + "\n")
print("\n".join(out)); print(json.dumps(summ, indent=1))
