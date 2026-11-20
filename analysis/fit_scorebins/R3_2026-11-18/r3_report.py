"""R3 — tabele, kryteria, verdict.json, rysunki z out/toys.npz, out/data.json, out/boot.npz."""
import json, os
import numpy as np
from scipy.stats import chi2 as C2
import r3_common as R

OUT = os.path.join(R.HERE, "out")
WPL = ["95", "80", "70", "50", "25"]
NC = 6
A = np.array([[1.0 if c >= j else 0.0 for c in range(NC)] for j in range(1, NC)])
L = []
P = L.append
V = {"toys": {}, "data": {}, "criteria": {}}


def unpack(rows, n, off=0):
    S = rows[:, off:off + n]
    good = rows[:, off + n] > 0.5
    cov = rows[:, off + n + 6:off + n + 6 + n * n].reshape(-1, n, n)
    return S, good, cov, rows[:, off + n + 3], rows[:, off + n + 4]


def rstd(x):
    return float(np.std(x, ddof=1))


# ------------------------------------------------------------------------------------------------ toys
T = np.load(os.path.join(OUT, "toys.npz"))
P("## T. Test wstrzyknięcia (pseudo-dane)\n")
P("| run | pT | SF70 prawda | N_S; B w oknie ±2,5 MeV (kat. 0…5) | dobre fity kształt / p-f | ⟨SF70⟩ kształt | obciążenie [%] (w σ) | σ toy kształt | szer. pulla | pokrycie 1σ | ⟨SF70⟩ p-f (obc. %) | σ toy p-f | σ_k/σ_pf | (a) |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
toy_other = []
for run in R.RUNS:
    for b in R.BINS:
        for c in ("1", "0"):
            k = f"{run}_{b}_{c}"
            rows = T[k + "_rows"]; fm, ft = T[k + "_f_mc"], T[k + "_f_true"]
            S, good, cov, mu, s = unpack(rows, NC)
            Sp, goodp, covp, _, _ = unpack(rows, 2, off=NC + 6 + NC * NC)
            Emc, Etr = A @ fm, A @ ft
            sft = Etr / Emc
            res = {"n": int(len(rows)), "good": float(good.mean()), "good_pf": float(goodp.mean()), "sf_true": sft.tolist()}
            Sg, cg = S[good], cov[good]
            N = Sg.sum(1); f = Sg / N[:, None]
            Ef = f @ A.T
            sf = Ef / Emc
            errs = np.empty_like(sf)
            for i in range(len(Sg)):
                Jf = (np.eye(NC) - f[i][:, None]) / N[i]
                errs[i] = np.sqrt(np.clip(np.diag(A @ Jf @ cg[i] @ Jf.T @ A.T), 0, None)) / Emc
            pull = (sf - sft) / errs
            res["wp"] = {}
            for j, w in enumerate(WPL):
                res["wp"][w] = {"mean": float(sf[:, j].mean()), "std": rstd(sf[:, j]), "bias_pct": float(100 * (sf[:, j].mean() / sft[j] - 1)),
                                "bias_sigma": float((sf[:, j].mean() - sft[j]) / rstd(sf[:, j])), "bias_err_pct": float(100 * rstd(sf[:, j]) / np.sqrt(len(sf)) / sft[j]),
                                "pull_width": rstd(pull[:, j]), "pull_mean": float(pull[:, j].mean()), "cover": float((np.abs(pull[:, j]) < 1).mean()),
                                "med_err": float(np.median(errs[:, j]))}
            both = good & goodp
            sfp = Sp[both, 1] / Sp[both].sum(1) / Emc[2]
            res["pf70"] = {"mean": float(sfp.mean()), "std": rstd(sfp), "bias_pct": float(100 * (sfp.mean() / sft[2] - 1))}
            sfk_both = (S[both] / S[both].sum(1)[:, None]) @ A.T[:, 2] / Emc[2]
            res["ratio_toy"] = rstd(sfk_both) / rstd(sfp)
            res["Ntot_mean"] = float(N.mean()); res["NS"] = float(T[k + "_NS"])
            w70 = res["wp"]["70"]
            res["a_pass"] = bool(abs(w70["bias_pct"]) <= 3 and abs(w70["bias_sigma"]) <= 1 and 0.8 <= w70["pull_width"] <= 1.25 and res["good"] >= 0.9)
            V["toys"][k] = res
            P(f"| {run} | {b} | {sft[2]:.3f} | {res['NS']:.0f}; {'/'.join(f'{x:.0f}' for x in T[k + '_B_win'])} | {100 * res['good']:.0f} % / {100 * res['good_pf']:.0f} % | "
              f"{w70['mean']:.3f} | {w70['bias_pct']:+.1f} ± {w70['bias_err_pct']:.1f} ({w70['bias_sigma']:+.2f}σ) | {w70['std']:.3f} | {w70['pull_width']:.2f} | "
              f"{100 * w70['cover']:.0f} % | {res['pf70']['mean']:.3f} ({res['pf70']['bias_pct']:+.1f}) | {res['pf70']['std']:.3f} | {res['ratio_toy']:.2f} | {'TAK' if res['a_pass'] else 'NIE'} |")
            toy_other.append(f"| {run} | {b} | {sft[2]:.2f} | " + " | ".join(
                f"{res['wp'][w]['bias_pct']:+.1f} % ({res['wp'][w]['bias_sigma']:+.2f}σ), σ {res['wp'][w]['std']:.3f}, pull {res['wp'][w]['pull_width']:.2f}" for w in ("95", "80", "50", "25")) + " |")
P("\nPozostałe WP (obciążenie względem prawdy, rozrzut toyów, szerokość pulla HESSE):\n")
P("| run | pT | SF70 prawda | WP95 | WP80 | WP50 | WP25 |"); P("|---|---|---|---|---|---|---|")
L.extend(toy_other)

# ------------------------------------------------------------------------------------------------ dane
D = json.load(open(os.path.join(OUT, "data.json")))
Bt = np.load(os.path.join(OUT, "boot.npz"))
off = {(r["run"], r["bin"], r["wp"]): r for r in R.official()["rows"] if r["cls"] in ("all", "all427929")}
P("\n## D. Dane\n")
P("### D1. Fit nominalny (6 kategorii score) i mierzalność\n")
P("| run | pT | wpisy w kat. 0…5 | param. tła (c1,c2 wolne: 2/1/0) | N_tot ± σ (Z) | μ [MeV] | skala s | fit dobry | param. na granicy | bootstrap: dobre | mierzalny |")
P("|---|---|---|---|---|---|---|---|---|---|---|")
meas = {}
tab_sf, tab_f, tab_corr, tab_eq, tab_cl = [], [], [], [], []
for run in R.RUNS:
    for b in R.BINS:
        k = f"{run}_{b}"; d = D[k]; nom = d["nom"]
        S = np.array(nom["S"]); covS = np.array(nom["covS"]); fm = np.array(d["f_mc"]); nmc = d["n_mc"]
        N, f, Ev, Nerr, covf, covE = R.cum_eff(S, covS)
        Emc = A @ fm; cmc = R.mc_cov_f(fm, nmc); cEmc = A @ cmc @ A.T
        rows = Bt[k]
        Sb, gb, _, _, _ = unpack(rows, NC)
        o = NC + 6 + NC * NC
        pfb = {}
        for i, j in enumerate((1, 2, 3, 4)):
            Sp, gp, _, _, _ = unpack(rows, 2, off=o + i * (2 + 6 + 4))
            pfb[j] = (Sp, gp)
        Z = N / Nerr
        ok = bool(Z >= 3 and nom["good"])
        meas[k] = ok
        fb = Sb[gb] / Sb[gb].sum(1)[:, None]
        Eb = fb @ A.T
        cov_fb = np.cov(fb.T); cov_Eb = np.cov(Eb.T)
        sf = Ev / Emc
        cov_sf_stat = cov_Eb / np.outer(Emc, Emc)
        cov_sf_mc = np.outer(sf, sf) * cEmc / np.outer(Emc, Emc)
        cov_sf = cov_sf_stat + cov_sf_mc
        e_boot = np.sqrt(np.diag(cov_sf_stat)); e_hes = np.sqrt(np.clip(np.diag(covE), 0, None)) / Emc; e_mc = np.sqrt(np.diag(cov_sf_mc))
        e_tot = np.sqrt(np.diag(cov_sf))
        corr = cov_sf / np.outer(e_tot, e_tot)
        res = {"N_tot": float(N), "N_err": Nerr, "Z": float(Z), "good": nom["good"], "measurable": ok, "mu": nom["mu"], "s": nom["s"], "at_limit": nom["at_limit"],
               "boot_good": float(gb.mean()), "n_boot": int(len(rows)), "f_data": f.tolist(), "f_mc": fm.tolist(), "E_mc": Emc.tolist(),
               "SF": sf.tolist(), "SF_err_boot": e_boot.tolist(), "SF_err_hesse": e_hes.tolist(), "SF_err_mc": e_mc.tolist(), "SF_err_tot": e_tot.tolist(),
               "SF_cov_tot": cov_sf.tolist(), "SF_corr_tot": corr.tolist(), "wps": WPL, "ncoef": d["ncoef"], "r": d["r"], "thr": d["thr"]}
        P(f"| {run} | {b} | {'/'.join(f'{x:.0f}' for x in d['n_cat'])} | {'/'.join(str(x) for x in d['ncoef'])} | {N:.0f} ± {Nerr:.0f} ({Z:.1f}) | {nom['mu']:.3f} ± {nom['mu_err']:.3f} | "
          f"{nom['s']:.2f} ± {nom['s_err']:.2f} | {'tak' if nom['good'] else 'NIE'} | {', '.join(nom['at_limit']) or '–'} | {100 * gb.mean():.0f} % z {len(rows)} | {'TAK' if ok else 'NIE'} |")
        # pass/fail tym samym silnikiem
        res["pf"] = {}
        for i, j in enumerate((1, 2, 3, 4)):
            w = WPL[i]; p = d["pf"][str(j)]
            Sp = np.array(p["S"]); cp = np.array(p["covS"])
            Np = Sp.sum(); ep = Sp[1] / Np
            Jp = np.array([-Sp[1], Sp[0]]) / Np ** 2
            ehp = float(np.sqrt(max(Jp @ cp @ Jp, 0))) / Emc[i]
            Sbp, gp = pfb[j]
            both = gb & gp
            sfp_b = Sbp[both, 1] / Sbp[both].sum(1) / Emc[i]
            sfk_b = (Sb[both] / Sb[both].sum(1)[:, None]) @ A.T[:, i] / Emc[i]
            ebp = rstd(sfp_b); ebk = rstd(sfk_b); ed = rstd(sfk_b - sfp_b)
            sfp = ep / Emc[i]
            o_ = off.get((run, b, f"0.{w}"))
            res["pf"][w] = {"SF": float(sfp), "err_boot": ebp, "err_hesse": ehp, "good": p["good"], "boot_good": float(gp.mean()), "N_tot": float(Np),
                            "shape_err_boot_same": ebk, "ratio_boot": ebk / ebp, "ratio_hesse": float(e_hes[i] / ehp), "diff": float(sf[i] - sfp), "diff_err": ed,
                            "diff_pull": float((sf[i] - sfp) / ed), "corr_shape_pf": float(np.corrcoef(sfk_b, sfp_b)[0, 1]),
                            "official": None if o_ is None else {"SF": o_["SF"], "err": o_["SF_err"], "determined": o_["determined"], "eps_mc": o_["eps_mc"]}}
            x = res["pf"][w]
            offs = "–" if o_ is None else f"{o_['SF']:.3f} ± {o_['SF_err']:.3f}{'' if o_['determined'] else ' (nieokr.)'}"
            tab_sf.append(f"| {run} | {b} | {w} | {sf[i]:.3f} ± {e_boot[i]:.3f} (H {e_hes[i]:.3f}; MC {e_mc[i]:.3f}) | {sfp:.3f} ± {ebp:.3f} (H {ehp:.3f}){'' if p['good'] else ' [fit zły]'} | {offs} | "
                          f"{x['ratio_boot']:.2f} (H {x['ratio_hesse']:.2f}) | {x['diff']:+.3f} ± {ed:.3f} ({x['diff_pull']:+.1f}σ) | {x['corr_shape_pf']:.2f} |")
            if o_ is not None:
                tab_cl.append(f"| {run} | {b} | {w} | {sfp:.3f} ± {ehp:.3f} | {o_['SF']:.3f} ± {o_['SF_err']:.3f} | {(sfp - o_['SF']) / o_['SF_err']:+.2f}σ | {Np:.0f} vs {o_['N_sig']:.0f} |")
        # ułamki i R_c
        Rc = f / fm
        cov_d = cov_fb + cmc
        ef = np.sqrt(np.diag(cov_fb)); eR = np.sqrt(np.diag(cov_fb / np.outer(fm, fm) + np.outer(Rc, Rc) * cmc / np.outer(fm, fm)))
        dd = (f - fm)[1:]
        chi = float(dd @ np.linalg.solve(cov_d[1:, 1:], dd))
        covH = covf + cmc
        chiH = float(dd @ np.linalg.solve(covH[1:, 1:], dd))
        res.update({"R": Rc.tolist(), "R_err": eR.tolist(), "f_err_boot": ef.tolist(), "chi2_f": chi, "p_f": float(C2.sf(chi, 5)), "chi2_f_hesse": chiH,
                    "p_f_hesse": float(C2.sf(chiH, 5)), "pull_R": ((Rc - 1) / eR).tolist(), "S_data": S.tolist(), "S_err_hesse": np.sqrt(np.diag(covS)).tolist()})
        tab_f.append(f"| {run} | {b} | " + " | ".join(f"{Rc[c]:.2f} ± {eR[c]:.2f} ({(Rc[c] - 1) / eR[c]:+.1f})" for c in range(NC)) + f" | {chi:.1f}/5 (p = {res['p_f']:.3f}); HESSE {chiH:.1f} (p = {res['p_f_hesse']:.3f}) |")
        tab_f.append(f"| | f_dane [%] (S_c ± σ_H) | " + " | ".join(f"{100 * f[c]:.1f} ({S[c]:.0f} ± {np.sqrt(covS[c, c]):.0f})" for c in range(NC)) + f" | f_MC [%]: {'/'.join(f'{100 * x:.1f}' for x in fm)} |")
        tab_corr.append(f"**{run}, {b} GeV** — SF ± σ_całk: " + ", ".join(f"WP{w} {sf[i]:.3f} ± {e_tot[i]:.3f}" for i, w in enumerate(WPL)) + "\n\n| ρ | " + " | ".join("WP" + w for w in WPL) + " |\n|---|" + "---|" * 5 + "\n" +
                        "\n".join(f"| WP{w} | " + " | ".join(f"{corr[i, j]:+.2f}" for j in range(5)) + " |" for i, w in enumerate(WPL)) + "\n")
        # wariant równych ułamków (HESSE)
        e = d["eq"]; Se = np.array(e["S"]); ce = np.array(e["covS"]); fme = np.array(e["f_mc"])
        Ne, fe, Ee, Nee, cfe, cEe = R.cum_eff(Se, ce)
        Re = fe / fme; eRe = np.sqrt(np.clip(np.diag(cfe), 0, None)) / fme
        de = (fe - fme)[1:]; chie = float(de @ np.linalg.solve((cfe + R.mc_cov_f(fme, nmc))[1:, 1:], de))
        res["eq"] = {"good": e["good"], "N_tot": float(Ne), "R": Re.tolist(), "R_err_hesse": eRe.tolist(), "chi2": chie, "p": float(C2.sf(chie, 5)),
                     "SF": (Ee / (A @ fme)).tolist(), "SF_err_hesse": (np.sqrt(np.clip(np.diag(cEe), 0, None)) / (A @ fme)).tolist()}
        tab_eq.append(f"| {run} | {b} | {'tak' if e['good'] else 'NIE'} | {Ne:.0f} ± {Nee:.0f} | " + " | ".join(f"{Re[c]:.2f} ± {eRe[c]:.2f}" for c in range(NC)) + f" | {chie:.1f}/5 (p = {res['eq']['p']:.3f}) |")
        # kontrola szerokości wspólnej (3–5)
        if d["nom_r1"]["x"] != nom["x"]:
            S1 = np.array(d["nom_r1"]["S"]); res["SF_r1"] = ((A @ (S1 / S1.sum())) / Emc).tolist()
        V["data"][k] = res

P("\n### D2. SF: kształt (6 kategorii) vs pass/fail — błąd = bootstrap sparowany po zdarzeniach (H = HESSE; MC = statystyka MC, wspólna dla obu metod)\n")
P("| run | pT | WP | SF kształt | SF pass/fail (ten silnik) | SF pass/fail oficjalny (btrig_sf) | σ_kształt/σ_pf | kształt − pf (sparowane) | ρ(kształt, pf) |")
P("|---|---|---|---|---|---|---|---|---|")
L.extend(tab_sf)
P("\n### D3. Kontrola silnika: pass/fail tego silnika (HESSE) vs oficjalny `fit_passfail`\n")
P("| run | pT | WP | ten silnik | oficjalny | różnica / σ_oficj. | N_tot |"); P("|---|---|---|---|---|---|---|")
L.extend(tab_cl)
P("\n### D4. Ułamki sygnału w kategoriach score: R_c = f_c(dane)/f_c(MC) (błąd bootstrap ⊕ MC; w nawiasie pull względem 1)\n")
P("| run | pT | kat. 0 (< WP95) | kat. 1 (95–80) | kat. 2 (80–70) | kat. 3 (70–50) | kat. 4 (50–25) | kat. 5 (> WP25) | χ²(dane = MC) |")
P("|---|---|---|---|---|---|---|---|---|")
L.extend(tab_f)
P("\n### D5. Macierze korelacji SF między WP (bootstrap ⊕ MC)\n")
L.extend(tab_corr)
P("\n### D6. Wariant kontrolny: równe ułamki sygnału (odrzucone < WP95, potem 5 × 19 %), błędy HESSE\n")
P("| run | pT | fit dobry | N_tot | R_0 | R_1 | R_2 | R_3 | R_4 | R_5 | χ² |"); P("|---|---|---|---|---|---|---|---|---|---|---|")
L.extend(tab_eq)

# ------------------------------------------------------------------------------------------------ kryteria
P("\n## K. Kryteria (KRYTERIA.md §3)\n")
P("| run | pT | mierzalny | (a) c=1 / c=0 | (b) σ_k/σ_pf przy WP70 (≤ 0,8) | (c) |Δ|/σ_Δ przy WP70 (≤ 2) | werdykt binu |"); P("|---|---|---|---|---|---|---|")
worst = 3
rank = {"UŻYTECZNA": 3, "UŻYTECZNA CZĘŚCIOWO": 2, "NIEUŻYTECZNA": 1}
for run in R.RUNS:
    for b in R.BINS:
        k = f"{run}_{b}"; r = V["data"][k]
        a1, a0 = V["toys"][k + "_1"]["a_pass"], V["toys"][k + "_0"]["a_pass"]
        x = r["pf"]["70"]
        bb = x["ratio_boot"] <= 0.8; cc = abs(x["diff_pull"]) <= 2
        if not meas[k]:
            v = "pominięty (niemierzalny)"
        else:
            v = "NIEUŻYTECZNA" if not (a1 and a0) else ("UŻYTECZNA" if (bb and cc) else ("UŻYTECZNA CZĘŚCIOWO" if cc else "NIEUŻYTECZNA CZĘŚCIOWO: (a) tak, (c) nie"))
            worst = min(worst, rank.get(v, 1))
        V["criteria"][k] = {"measurable": meas[k], "a_c1": a1, "a_c0": a0, "b": bool(bb), "b_ratio": x["ratio_boot"], "c": bool(cc), "c_pull": x["diff_pull"], "verdict": v}
        P(f"| {run} | {b} | {'tak' if meas[k] else 'NIE'} | {'TAK' if a1 else 'NIE'} / {'TAK' if a0 else 'NIE'} | {x['ratio_boot']:.2f} → {'TAK' if bb else 'NIE'} | {abs(x['diff_pull']):.1f} → {'TAK' if cc else 'NIE'} | {v} |")
V["verdict"] = {3: "UŻYTECZNA", 2: "UŻYTECZNA CZĘŚCIOWO", 1: "NIEUŻYTECZNA"}[worst]
P(f"\n**Werdykt łączny (najgorszy z mierzalnych): {V['verdict']}**\n")

json.dump(V, open(os.path.join(R.HERE, "verdict.json"), "w"), indent=1)
open(os.path.join(R.HERE, "tables.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))

# ------------------------------------------------------------------------------------------------ rysunki
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
INK, MUTED, C1, C2c, SURF = "#0b0b0b", "#52514e", "#2a78d6", "#eb6834", "#fcfcfb"
names = ["kat. 0: < WP95", "kat. 1: WP95–80", "kat. 2: WP80–70", "kat. 3: WP70–50", "kat. 4: WP50–25", "kat. 5: > WP25"]
for run in R.RUNS:
    fig, ax = plt.subplots(3, 6, figsize=(21, 9.5), facecolor=SURF)
    for i, b in enumerate(R.BINS):
        d = D[f"{run}_{b}"]; H = np.array(d["H"]); lam = np.array(d["lam"]); lb = np.array(d["lam_bkg"])
        rb = 5
        for c in range(6):
            a = ax[i, c]; a.set_facecolor(SURF)
            h = H[c].reshape(-1, rb).sum(1); x = R.CEN.reshape(-1, rb).mean(1)
            a.errorbar(x, h, np.sqrt(h), fmt="o", ms=2.5, lw=0.8, color=INK, label="dane", zorder=3)
            a.plot(x, lam[c].reshape(-1, rb).sum(1), color=C1, lw=2, label="fit: sygnał + tło")
            a.plot(x, lb[c].reshape(-1, rb).sum(1), color=C2c, lw=2, ls="--", label="fit: tło")
            a.set_xlim(140.5, 160); a.grid(alpha=0.25, lw=0.5); a.tick_params(colors=MUTED, labelsize=8)
            for sp in ("top", "right"):
                a.spines[sp].set_visible(False)
            S_ = d["nom"]["S"][c]; eS = np.sqrt(np.array(d["nom"]["covS"])[c, c])
            a.set_title(f"pT {b} GeV, {names[c]}\nS = {S_:.0f} ± {eS:.0f}", fontsize=9, color=INK)
            if i == 2:
                a.set_xlabel("Δm = m(Kππ) − m(Kπ) [MeV]", fontsize=9, color=MUTED)
            if c == 0:
                a.set_ylabel("pary / 0,5 MeV", fontsize=9, color=MUTED)
    ax[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle(f"Run {run}: Δm kandydatów D⁰←D* w kategoriach score H, fit jednoczesny (R3)", fontsize=12, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(R.HERE, f"fig_dm_{run}.png"), dpi=110)
    plt.close(fig)
