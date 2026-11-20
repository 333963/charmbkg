"""DLZ_dplus/DATA: tabele (tables.md), werdykt wg KRYTERIA.md (verdict.json) i rysunki PNG z wyników węzła `ana`
(runs/lxplus/ana/res_<run>.json, mc_checks.json) lub lokalnych (out/).  Użycie: python report.py [KATALOG_Z_res_*.json]"""
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "runs", "lxplus", "ana")
FIG = os.path.join(BASE, "fig")
os.makedirs(FIG, exist_ok=True)
BINS = ["3-5", "5-8", "8-30"]
WPS = [0.95, 0.9, 0.8, 0.7, 0.5]
H, D, DAL = "H_s12345", "HDLZ_s12345", "Hdal_s12345"
RUNS = [r for r in ("427929", "428770", "428770loc") if os.path.exists(os.path.join(SRC, f"res_{r}.json"))]
R = {r: json.load(open(os.path.join(SRC, f"res_{r}.json"))) for r in RUNS}
MC = json.load(open(os.path.join(SRC, "mc_checks.json"))) if os.path.exists(os.path.join(SRC, "mc_checks.json")) else {}
LO, HI, NB, M0 = 1.70, 2.20, 100, 1.86966
C = LO + (np.arange(NB) + 0.5) * 0.005
out = []
P = out.append
pm = lambda v, e, f="{:.0f}": f"{f.format(v)} ± {f.format(e)}"  # noqa: E731


def atlas(ax, txt=""):
    ax.text(0.04, 0.95, "ATLAS", transform=ax.transAxes, fontweight="bold", fontstyle="italic", va="top", fontsize=11)
    ax.text(0.19, 0.95, "", transform=ax.transAxes, va="top", fontsize=11)
    if txt:
        ax.text(0.04, 0.87, txt, transform=ax.transAxes, va="top", fontsize=8)


verdict = {}
for run in RUNS:
    r = R[run]
    F, DR, SC, SP = r["fits"], r["derived"], r["sculpt"], r["splot"]
    DX, FX, AL = r.get("derived_fixshape", {}), r.get("fits_fixshape", {}), r.get("at_limit", {})
    P(f"\n## Run {run}\n")
    P(f"Kandydatów D⁺ po weto D* (pT 3–30 GeV, 1,70–2,20 GeV, L_xy ≥ 0,25 mm): {r['n_cand_all']}; plików h5: {r['n_files']}; "
      f"wierszy „luźnych” (≥ WP95 dowolnego modelu): {r['n_loose']} w {r['n_events_loose']} zdarzeniach.\n")
    P("### Fit masy: N_S, błąd względny, S/B (±45 MeV), istotność\n")
    P("(† = parametr kształtu sygnału (dμ lub skala σ) na granicy w fitcie swobodnym → dla tego wiersza wiarygodny jest wariant „kształt ustalony na MC”.)\n")
    P("| bin pT | WP | H: N_S | σ/N [%] | S/B | Z | H+DLZ: N_S | σ/N [%] | S/B | Z | zysk σ/N [%] (hesjan; ± bootstrap) | zysk (rozrzut replik) [%] | χ²/ndf H / DLZ | kształt ustalony: N_H / N_DLZ | zysk [%] |")
    P("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for b in BINS:
        for wp in WPS:
            a, c, d = F[f"{H}|{b}|{wp}"], F[f"{D}|{b}|{wp}"], DR[f"{D}|{b}|{wp}"]
            P(f"| {b} | {int(wp*100)} | {pm(a['N'], a['N_err'])} | {100*d['rel_H']:.1f} | {a['SB']:.3f} | {a['Z']:.1f} | {pm(c['N'], c['N_err'])} | {100*d['rel_D']:.1f} | {c['SB']:.3f} | {c['Z']:.1f} | "
              f"{100*d['gain']:+.1f} ± {100*d['gain_err']:.1f}{'†' if AL.get(f'{H}|{b}|{wp}') or AL.get(f'{D}|{b}|{wp}') else ''} | {100*d['gain_boot']:+.1f} | {a['chi2_ndf']:.2f} / {c['chi2_ndf']:.2f} | "
              + (f"{pm(FX[f'{H}|{b}|{wp}']['N'], FX[f'{H}|{b}|{wp}']['N_err'])} / {pm(FX[f'{D}|{b}|{wp}']['N'], FX[f'{D}|{b}|{wp}']['N_err'])} | {100*DX[f'{D}|{b}|{wp}']['gain']:+.1f} ± {100*DX[f'{D}|{b}|{wp}']['gain_err']:.1f} |" if DX else "— | — |"))
    P("\n### Stosunek uzysków, pik, SF względny (ΔSF_rel = SF_rel(DLZ)/SF_rel(H); błędy: bootstrap po zdarzeniach, sparowany)\n")
    P("| bin pT | WP | N_DLZ/N_H | Δμ (DLZ−H) [MeV] | skala σ DLZ/H | μ_H − μ_MC [MeV] | skala σ_H | SF_rel H (ref WP90) | SF_rel DLZ | ΔSF_rel (ref WP90) | ΔSF_rel (ref WP95) | połowy danych ΔSF (ref 90) |")
    P("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for b in BINS:
        for wp in WPS:
            a, d = F[f"{H}|{b}|{wp}"], DR[f"{D}|{b}|{wp}"]
            g = lambda k, f="{:.3f}": (f"{f.format(d[k])} ± {f.format(d[k + '_err'])}" if k in d else "—")  # noqa: E731
            hv = " / ".join(f"{x:.3f}" for x in d["dSF_ref90_halves"]) if "dSF_ref90_halves" in d else "—"
            P(f"| {b} | {int(wp*100)} | {d['Nratio']:.3f} ± {d['Nratio_err']:.3f} | {1e3*d['ddmu']:+.2f} ± {1e3*d['ddmu_err']:.2f} | {d['scratio']:.3f} ± {d['scratio_err']:.3f} | "
              f"{1e3*a['dmu']:+.2f} ± {1e3*a['dmu_err']:.2f} | {a['sc']:.3f} ± {a['sc_err']:.3f} | {g('SFrel_H_ref90')} | {g('SFrel_D_ref90')} | {g('dSF_ref90')} | {g('dSF_ref95')} | {hv} |")
    if DX:
        P("\nΔSF_rel z fitów o kształcie ustalonym na MC (ref WP90): " + "; ".join(
            f"{b} WP{int(wp*100)}: {DX[f'{D}|{b}|{wp}']['dSF_ref90']:.3f} ± {DX[f'{D}|{b}|{wp}']['dSF_ref90_err']:.3f}" for b in BINS for wp in (0.8, 0.7, 0.5)) + ".")
    if f"{DAL}|3-5|0.7" in DR:
        P("\n### Wariant zapasowy H+dalitz (WP70)\n")
        P("| bin pT | N_S | σ/N [%] | zysk σ/N vs H [%] | ΔSF_rel (ref WP90) |")
        P("|---|---|---|---|---|")
        for b in BINS:
            c, d = F[f"{DAL}|{b}|0.7"], DR[f"{DAL}|{b}|0.7"]
            P(f"| {b} | {pm(c['N'], c['N_err'])} | {100*d['rel_D']:.1f} | {100*d['gain']:+.1f} ± {100*d['gain_err']:.1f} | {d['dSF_ref90']:.3f} ± {d['dSF_ref90_err']:.3f} |")
    P("\n### Warianty modelu tła / okna (WP70): N_S\n")
    P("| bin pT | klasyfikator | nominał exp(pol2) | exp(pol1) | exp(pol3) | okno 1,72–2,10 | kształt MC ustalony |")
    P("|---|---|---|---|---|---|---|")
    for b in BINS:
        for t, nm in ((H, "H"), (D, "H+DLZ")):
            a, al = F[f"{t}|{b}|0.7"], r["alt"][f"{t}|{b}|0.7"]
            P(f"| {b} | {nm} | {pm(a['N'], a['N_err'])} | " + " | ".join(f"{pm(al[k]['N'], al[k]['N_err'])} (χ² {al[k]['chi2_ndf']:.2f})" for k in ("pol1", "pol3", "win", "fixshape")) + " |")
    P("\n### Rzeźbienie w danych: ε_B(m) w SB; R = ε_B(pasmo 45–70 MeV) / przewidywanie z liniowego fitu w |Δm| > 70 MeV\n")
    P("| bin pT | WP | R (H) | R (H+DLZ) | nachylenie H [/100 MeV] | nachylenie DLZ | χ²/ndf fitu lin. H / DLZ | χ²_SB/n_SB fitu masy H / DLZ |")
    P("|---|---|---|---|---|---|---|---|")
    for b in BINS + ["3-30"]:
        for wp in (0.8, 0.7):
            a, c = SC[f"{H}|{b}|{wp}"], SC[f"{D}|{b}|{wp}"]
            fs = "—"
            if b != "3-30":
                fa, fc = F[f"{H}|{b}|{wp}"], F[f"{D}|{b}|{wp}"]
                fs = f"{fa['chi2_sb']/(fa['n_sb']-fa['npar_b']):.2f} / {fc['chi2_sb']/(fc['n_sb']-fc['npar_b']):.2f}"
            P(f"| {b} | {int(wp*100)} | {a['R']:.3f} ± {a['R_err']:.3f} | {c['R']:.3f} ± {c['R_err']:.3f} | {a['slope']:+.3f} ± {a['slope_err']:.3f} | {c['slope']:+.3f} ± {c['slope_err']:.3f} | "
              f"{a['chi2_lin']/a['ndf_lin']:.2f} / {c['chi2_lin']/c['ndf_lin']:.2f} | {fs} |")
    P("\n### sPlot przy H WP90: ε warunkowe sygnału dane vs MC\n")
    P("| bin pT | Σ sWeights (N_fit) | cięcie | ε dane | ε MC | SF warunkowy |")
    P("|---|---|---|---|---|---|")
    for b in BINS:
        for t, nm in ((H, "H"), (D, "H+DLZ"), (DAL, "H+dalitz")):
            for wp in (80, 70, 50):
                k = f"{t}_WP{wp}"
                if k in SP[b]:
                    x = SP[b][k]
                    P(f"| {b} | {SP[b]['sum_sw']:.0f} ({SP[b]['N_fit']:.0f}) | {nm} WP{wp} | {x['eps_data']:.3f} ± {x['eps_data_err']:.3f} | {x['eps_mc']:.3f} ± {x['eps_mc_err']:.3f} | {x['sf']:.3f} ± {x['sf_err']:.3f} |")
    if "dMinPt_bins" in SP:
        P(f"\nW tercylach dMinPt (granice {SP['dMinPt_bins']['edges'][0]:.2f}, {SP['dMinPt_bins']['edges'][1]:.2f} GeV; pT-integralnie), WP70 | H WP90:\n")
        P("| klasyfikator | tercyl | udział dane / MC | ε dane | ε MC | SF |")
        P("|---|---|---|---|---|---|")
        for t, nm in ((D, "H+DLZ"), (H, "H")):
            for q, x in enumerate(SP["dMinPt_bins"].get(t, [])):
                P(f"| {nm} | {q+1} | {x['frac_data']:.3f} / {x['frac_mc']:.3f} | {x['eps_data']:.3f} ± {x['eps_data_err']:.3f} | {x['eps_mc']:.3f} | {x['sf']:.3f} ± {x['sf_err']:.3f} |")
    if "feat1d" in SP:
        P("\nRozkłady 1D cech sygnału (sWeights, H WP90) dane vs MC, χ²/9 w 10 kwantylach MC:\n")
        P("| cecha | χ²/ndf | średnia dane | średnia MC |")
        P("|---|---|---|---|")
        for f, x in SP["feat1d"].items():
            if f.endswith("_cn"):
                P(f"| {f} | {x['chi2']/x['ndf']:.2f} | {x['mean_data']:.4f} | {x['mean_mc']:.4f} |")
        hh = [x['chi2'] / x['ndf'] for f, x in SP["feat1d"].items() if not f.endswith("_cn")]
        P(f"\n18 cech H: χ²/ndf mediana {np.median(hh):.2f}, zakres {min(hh):.2f}–{max(hh):.2f}.")
    if "c2st_signal" in SP:
        c = SP["c2st_signal"]
        P("\nC2ST sygnału (dane sWeighted vs MC, OOF, 3 ziarna; błąd = bootstrap po zdarzeniach):\n")
        P("| zestaw | AUC (ziarna) | błąd bootstrap |")
        P("|---|---|---|")
        for nm in ("H", "HDLZ", "Hdal"):
            if nm in c:
                P(f"| {nm} | {' / '.join(f'{x:.4f}' for x in c[nm])} | {c[nm + '_booterr']:.4f} |")
        for nm in ("HDLZ", "Hdal"):
            if f"d_{nm}_minus_H" in c:
                P(f"\nΔC2ST({nm} − H) = {c[f'd_{nm}_minus_H']:+.4f} ± {c[f'd_{nm}_minus_H_err']:.4f}")
    # ziarna
    if f"H_s23456|3-5|0.7" in F:
        P("\n### Ziarna treningu (WP70): N_S seed 12345 vs 23456\n")
        P("| bin pT | H: N₁ / N₂ (Δ/σ) | H+DLZ: N₁ / N₂ (Δ/σ) |")
        P("|---|---|---|")
        for b in BINS:
            row = []
            for t in ("H", "HDLZ"):
                a, c = F[f"{t}_s12345|{b}|0.7"], F[f"{t}_s23456|{b}|0.7"]
                row.append(f"{a['N']:.0f} / {c['N']:.0f} ({(c['N']-a['N'])/a['N_err']:+.2f})")
            P(f"| {b} | {row[0]} | {row[1]} |")
    # --- werdykt
    v = {}
    meas = [b for b in BINS if F[f"{H}|{b}|0.7"]["N"] / max(F[f"{H}|{b}|0.7"]["N_err"], 1e-9) >= 3]
    v["mierzalne_biny"] = meas
    lim = {b: bool(any(AL.get(f"{t}|{b}|{wp}") for t in (H, D) for wp in (0.7, 0.9))) for b in BINS}
    v["fit_swobodny_na_granicy"] = lim
    DRn = DR
    DR = {k: (DX[k] if (DX and lim.get(k.split("|")[1])) else x) for k, x in DR.items()}   # bin z parametrem na granicy → wariant z kształtem ustalonym
    g = [DR[f"{D}|{b}|0.7"]["gain"] for b in meas]
    v["K1"] = dict(gain={b: DR[f"{D}|{b}|0.7"]["gain"] for b in meas}, gain_err={b: DR[f"{D}|{b}|0.7"]["gain_err"] for b in meas},
                   mean=float(np.mean(g)) if g else None, PASS=bool(g and all(x >= 0 for x in g) and np.mean(g) >= 0.08))
    k2 = True
    for b in BINS + ["3-30"]:
        for wp in (0.8, 0.7):
            a, c = SC[f"{H}|{b}|{wp}"], SC[f"{D}|{b}|{wp}"]
            okR = abs(c["R"] - a["R"]) < 2 * np.hypot(a["R_err"], c["R_err"]) or abs(c["R"] - 1) <= abs(a["R"] - 1)
            okS = abs(c["slope"] - a["slope"]) <= 0.10
            k2 &= bool(okR and okS)
    chi_sb = {b: F[f"{D}|{b}|0.7"]["chi2_sb"] / (F[f"{D}|{b}|0.7"]["n_sb"] - F[f"{D}|{b}|0.7"]["npar_b"]) for b in BINS}
    v["K2"] = dict(PASS_R_slope=bool(k2), chi2_sb_fit_DLZ_WP70=chi_sb, PASS_chi2=bool(all(x < 1.5 for x in chi_sb.values())), PASS=bool(k2 and all(x < 1.5 for x in chi_sb.values())))
    k3 = {b: (DR[f"{D}|{b}|0.7"]["dSF_ref90"], DR[f"{D}|{b}|0.7"]["dSF_ref90_err"]) for b in meas}
    v["K3"] = dict(dSF=k3, PASS=bool(k3 and all(abs(x - 1) < 2 * e and abs(x - 1) < 0.10 for x, e in k3.values())),
                   pass_2sigma=bool(all(abs(x - 1) < 2 * e for x, e in k3.values())), pass_010=bool(all(abs(x - 1) < 0.10 for x, e in k3.values())))
    if "c2st_signal" in SP and "d_HDLZ_minus_H" in SP["c2st_signal"]:
        v["K4"] = dict(d=SP["c2st_signal"]["d_HDLZ_minus_H"], err=SP["c2st_signal"]["d_HDLZ_minus_H_err"], PASS=(None if SP["c2st_signal"]["d_HDLZ_minus_H_err"] > 0.02 else bool(SP["c2st_signal"]["d_HDLZ_minus_H"] <= 0.02)))
    k5 = {b: (DR[f"{D}|{b}|0.7"]["ddmu"], DR[f"{D}|{b}|0.7"]["ddmu_err"], DR[f"{D}|{b}|0.7"]["scratio"], DR[f"{D}|{b}|0.7"]["scratio_err"]) for b in meas}
    v["K5"] = dict(vals=k5, PASS=bool(all(abs(a) < 2 * max(e, 1e-9) and abs(a) < 0.002 and abs(s - 1) < 2 * max(se, 1e-9) and abs(s - 1) < 0.10 for a, e, s, se in k5.values())))
    if f"HDLZ_s23456|3-5|0.7" in F:
        v["K7_N"] = dict(d={b: (F[f"HDLZ_s23456|{b}|0.7"]["N"] - F[f"{D}|{b}|0.7"]["N"]) / F[f"{D}|{b}|0.7"]["N_err"] for b in meas})
        v["K7_N"]["PASS"] = bool(all(abs(x) < 1 for x in v["K7_N"]["d"].values()))
    verdict[run] = v
    DR = DRn
    # --- rysunki
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for j, b in enumerate(BINS):
        for i, (t, nm) in enumerate(((H, "H"), (D, "H+DLZ"))):
            ax = axs[i, j]
            k = f"{t}|{b}|0.7"
            h, lam, f = np.array(r["hists"][k]), np.array(r["lam"][k]), F[k]
            ax.errorbar(C, h, np.sqrt(np.maximum(h, 1)), fmt="k.", ms=3, lw=0.6, label="dane")
            ax.plot(C, lam, "C3-", lw=1.2, label="fit")
            ax.set_ylim(0, 1.45 * max(h.max(), 1))
            atlas(ax, f"run {run}, D$^+\\to K\\pi\\pi$, {b} GeV, {nm} WP70\n$N_S$ = {f['N']:.0f} ± {f['N_err']:.0f} ({100*f['N_err']/max(f['N'],1e-9):.1f} %), S/B = {f['SB']:.3f}, Z = {f['Z']:.1f}\n"
                      f"$\\chi^2$/ndf = {f['chi2_ndf']:.2f}")
            if i == 1:
                ax.set_xlabel("m(K$\\pi\\pi$) [GeV]")
            if j == 0:
                ax.set_ylabel("kandydaci / 5 MeV")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"fig1_massfit_WP70_{run}.png"), dpi=130)
    plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    e20 = LO + (np.arange(25) + 0.5) * 0.02
    for j, b in enumerate(BINS):
        ax = axs[j]
        for t, nm, col in ((H, "H", "C0"), (D, "H+DLZ", "C3")):
            s = SC[f"{t}|{b}|0.7"]
            e, ee = np.array(s["eps"]), np.array(s["eps_err"])
            sb = np.abs(e20 - M0) > 0.045
            nrm = np.mean(e[sb])
            ax.errorbar(e20[sb], e[sb] / nrm, ee[sb] / nrm, fmt="o", ms=3, color=col, label=f"{nm}: R = {s['R']:.3f} ± {s['R_err']:.3f}, nachyl. {s['slope']:+.2f}/100 MeV")
        ax.axvspan(M0 - 0.045, M0 + 0.045, color="0.9")
        ax.set_ylim(0.4, 1.9)
        ax.legend(fontsize=7, loc="lower left")
        ax.set_xlabel("m(K$\\pi\\pi$) [GeV]")
        ax.set_ylabel("$\\varepsilon_B(m)$ / średnia w SB")
        atlas(ax, f"run {run}, {b} GeV, WP70, pasma boczne danych")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"fig3_sculpt_{run}.png"), dpi=130)
    plt.close(fig)
    if "feat1d" in SP:
        fs = [f for f in SP["feat1d"] if f.endswith("_cn")]
        fig, axs = plt.subplots(2, 4, figsize=(16, 6.5))
        for ax, f in zip(axs.ravel(), fs):
            x = SP["feat1d"][f]
            ax.errorbar(np.arange(10) + 0.5, x["ratio"], np.abs(np.nan_to_num(x["ratio_err"])), fmt="ko", ms=3)
            ax.axhline(1, color="C3", lw=0.8)
            ax.set_ylim(0.3, 1.9)
            ax.set_xlabel(f"{f} (decyl MC)")
            ax.set_ylabel("dane (sWeights) / MC")
            atlas(ax, f"run {run}, sygnał D$^+$, H WP90\n$\\chi^2$/ndf = {x['chi2']/x['ndf']:.2f}")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, f"fig5_splot_features_{run}.png"), dpi=130)
        plt.close(fig)

# rysunki wspólne: zysk i ΔSF
if RUNS:
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    for j, b in enumerate(BINS):
        ax = axs[j]
        for run, col, off in zip(RUNS, ("C0", "C3", "C2"), (-0.6, 0.6, 0)):
            d = [R[run]["derived"][f"{D}|{b}|{wp}"] for wp in WPS]
            ax.errorbar(np.array(WPS) * 100 + off, [100 * x["gain"] for x in d], [100 * x["gain_err"] for x in d], fmt="o", color=col, capsize=2, label=f"run {run}")
        ax.axhline(0, color="k", lw=0.6)
        ax.axhline(8, color="0.5", lw=0.6, ls="--")
        ax.set_xlabel("WP ($\\varepsilon_S^{MC}$) [%]")
        if j == 0:
            ax.set_ylabel("zysk $1-(\\sigma_N/N)_{H+DLZ}/(\\sigma_N/N)_{H}$ [%]")
        ax.legend(fontsize=8, loc="lower left")
        atlas(ax, f"D$^+$, {b} GeV; linia przerywana: próg 8 %")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_gain.png"), dpi=130)
    plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    for j, b in enumerate(BINS):
        ax = axs[j]
        for run, col, off in zip(RUNS, ("C0", "C3", "C2"), (-0.6, 0.6, 0)):
            w = [wp for wp in WPS if wp < 0.9]
            d = [R[run]["derived"][f"{D}|{b}|{wp}"] for wp in w]
            ax.errorbar(np.array(w) * 100 + off, [x["dSF_ref90"] for x in d], [x["dSF_ref90_err"] for x in d], fmt="o", color=col, capsize=2, label=f"run {run}")
        ax.axhline(1, color="k", lw=0.6)
        ax.axhspan(0.9, 1.1, color="0.92")
        ax.set_ylim(0.5, 1.5)
        ax.set_xlabel("WP [%]")
        if j == 0:
            ax.set_ylabel("$\\Delta SF_{rel}$ = SF$_{rel}$(H+DLZ)/SF$_{rel}$(H), ref. WP90")
        ax.legend(fontsize=8, loc="lower left")
        atlas(ax, f"D$^+$, {b} GeV; pas: ±0,10")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_dSF.png"), dpi=130)
    plt.close(fig)

if MC:
    P("\n## Kontrole na MC (bez danych)\n")
    P("| model | AUC trening − OOF | best_it | naruszenia foldów |")
    P("|---|---|---|---|")
    for t, x in MC["auc"].items():
        P(f"| {t} | {x['diff']:+.4f} (trening {np.mean(x['train']):.4f}, OOF {np.mean(x['oof']):.4f}) | {x['best_it']} | {x['fold_violations']} |")
    P("\n| model | WP | R_sculpt MC truth ≥ 8 (okno ±45 MeV) | nachylenie ε_B [/100 MeV] |")
    P("|---|---|---|---|")
    for k, x in MC["sculpt_mc"].items():
        t, wp = k.split("|")
        P(f"| {t} | {int(float(wp)*100)} | {x['R']:.3f} ± {x['R_err']:.3f} | {x['slope']:+.3f} ± {x['slope_err']:.3f} |")
    P("\nClosure fitu na mb_nd (ε_fit(WP)/ε_fit(WP90) względem zliczenia prawdy):\n")
    P("| model | bin pT | WP | N_fit(WP) | N_prawda(WP) | closure ε_fit/ε_prawda |")
    P("|---|---|---|---|---|---|")
    for k, x in MC["closure"].items():
        t, b, wp = k.split("|")
        if t in (H, D) and float(wp) in (0.7, 0.5):
            P(f"| {t} | {b} | {int(float(wp)*100)} | {x['N_fit']:.0f} ± {x['N_fit_err']:.0f} | {x['N_truth']} | {x['ratio']:.3f} ± {x['ratio_err']:.3f} |")
open(os.path.join(BASE, "tables.md"), "w").write("\n".join(out) + "\n")
json.dump(verdict, open(os.path.join(BASE, "verdict.json"), "w"), indent=1)
print(json.dumps(verdict, indent=1))
