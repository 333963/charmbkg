#!/usr/bin/env python3
"""SIG_constr — tabele Markdown z out/*.json → out/tables.md (T1–T7)."""
import json, os, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
J = lambda n: json.load(open(os.path.join(OUT, n))) if os.path.exists(os.path.join(OUT, n)) else None  # noqa: E731
f = lambda x, n=3: "–" if x is None or not np.isfinite(x) else f"{x:.{n}f}".replace(".", ",")  # noqa: E731
pm = lambda x, e, n=3: f"{f(x, n)} ± {f(e, n)}"  # noqa: E731
L = ["# SIG_constr — tabele", ""]
d = J("s1_neff.json")
if d:
    L += ["## T1. Sygnał z tagu D*: N_sig i n_eff (okno |m − m₀| < 40 MeV)", "",
          "| klasa | obszar | pT | par | w oknie Δm | czystość okna | N_sig (SBS) | n_eff (SBS) | N_sig (sW) | n_eff (sW) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in d["rows"]:
        sw = (pm(r["Ns_sw"], r["Ns_sw_err"], 0), f(r["neff_sw"], 0)) if r["Ns_sw"] else ("–", "–")
        L.append(f"| {r['cls']} | {r['reg']} | {r['bin']} | {r['n_pairs']} | {r['n_win']} | {f(r['pur_win'], 2)} | {pm(r['Ns_sbs'], r['Ns_sbs_err'], 0)} | {f(r['neff_sbs'], 0)} | {sw[0]} | {sw[1]} |")
    L += ["", f"MC fd_nd, prawdziwe D⁰←D*: {d['mc_sig']}", ""]
d = J("s2_learn.json")
if d:
    L += ["## T2. Krzywa uczenia (sygnał MC, N zdarzeń; tło: pasma Δm 9:MB)", "", "| zadanie | N | AUC | Z_rel WP95 | Z_rel WP80 | Z_rel WP50 |", "|---|---|---|---|---|---|"]
    for t in d:
        for s in d[t]["summary"]:
            L.append(f"| {t} | {s['N']} | {f(s['auc'], 4)} | {f(s['Zrel0.95'])} | {pm(s['Zrel0.8'], s['Zrel0.8_sd'])} | {f(s['Zrel0.5'])} |")
        L.append(f"| {t} | N*(2/5/10 %) | {d[t]['Nstar_2']} / {d[t]['Nstar_5']} / {d[t]['Nstar_10']} | | | |")
    L.append("")
rs = []
for p in sorted(glob.glob(os.path.join(OUT, "s3_variants_*.json"))):
    rs += json.load(open(p))
if rs:
    L += ["## T3. Warianty treningu na danych (OOF): ε_S(dane) i Z_rel względem modelu MC przy tym samym ε_B", "",
          "| klasa | zadanie | WP (ε_B) | pT | ε_S MC-model: dane / MC | RW: ε_S, Z_rel | DATA: ε_S, Z_rel | STACK: ε_S, Z_rel | H: ε_S, Z_rel |", "|---|---|---|---|---|---|---|---|---|"]
    for x in rs:
        for wp, res in x["wp"].items():
            for b, rb in res.items():
                if b == "epsB_target":
                    continue
                c = lambda m: f"{pm(rb[m]['epsS'], rb[m]['epsS_err'])}; {pm(rb[m]['Zrel_vs_MC'], rb[m]['Zrel_vs_MC_err'])}"  # noqa: E731
                L.append(f"| {x['cls']} | {x['task']} | {wp} ({f(res['epsB_target'], 4)}) | {b} | {pm(rb['MC']['epsS'], rb['MC']['epsS_err'])} / {f(rb['MC']['epsS_mc'])} | {c('RW')} | {c('DATA')} | {c('STACK')} | {c('H')} |")
    L += ["", "Wyciek (R² score vs Δm / m(Kπ) na tle, pasmo górne) i wagi RW:", "", "| klasa | zadanie | N_sig (SBS) | n_eff | model: R²(Δm), R²(m) | wagi RW: sd, n_eff/N |", "|---|---|---|---|---|---|"]
    for x in rs:
        L.append(f"| {x['cls']} | {x['task']} | {f(x['Ns_sbs'], 0)} | {f(x['neff_sbs'], 0)} | " + "; ".join(f"{m}: {f(l['R2_dm_bkg'], 4)}, {f(l['R2_m_bkg'], 4)}" for m, l in x["leak"].items()) +
                 f" | {f(x['rw_weights']['sd'], 2)}, {f(x['rw_weights']['neff_frac'], 2)} |")
    L.append("")
d = J("s4_closure.json")
if d:
    for key, ttl in (("iso", "T4. Domknięcie wag sygnałowych — izolacja (H06)"), ("H", "T5. Domknięcie wag sygnałowych — score H (PRE)")):
        L += [f"## {ttl}", "", "| klasa | pT | wagi | cięcie | N_sig | ε_dane | ε_MC | ε_MC,rw ± MC ± wagi | dane/MC_rw | pull | n_eff/N |", "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in d[key]:
            L.append(f"| {r['cls']} | {r['bin']} | {r['set']} | {r['cut']} | {f(r['Ns'], 0)} | {pm(r['eps_data'], r['eps_data_err'])} ({r['method']}) | {f(r['eps_mc_raw'])} | "
                     f"{f(r['eps_mc_rw'])} ± {f(r['err_mcstat'])} ± {f(r['err_w_jk'])} | {pm(r['ratio'], r['ratio_err'])} | {f(r['pull'], 1)} | {f(r['neff_frac'], 2)} |")
        L.append("")
    L += ["## T6. Przenoszenie wag między klasami (ptcone40 < 2 GeV, H06)", "", "| wagi z klasy | do klasy | zestaw | ε_dane(B) | ε przewidziane | dane/przew. | pull |", "|---|---|---|---|---|---|---|"]
    for r in d["transfer"]:
        L.append(f"| {r['from']} | {r['to']} | {r['set']} | {pm(r['eps_data_B'], r['eps_data_B_err'])} | {f(r['eps_pred'])} | {pm(r['ratio'], r['ratio_err'])} | {f(r['pull'], 1)} |")
    L.append("")
d = J("s5_transfer.json")
if d:
    L += ["## T7. Przenoszenie D⁰←D* → inkluzywne D⁰ (MC fd_nd, prawda)", "", "| cięcie | pT | N inkl. | N tag | ε inkl. | ε tag | ε tag (pT-rw) | Δ | Δ/ε [%] | pull |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in d["rows"]:
        L.append(f"| {r['test']} | {r['bin']} | {r['n_incl']} | {r['n_tag']} | {f(r['eps_incl'], 4)} | {f(r['eps_tag_raw'], 4)} | {f(r['eps_tag_ptrw'], 4)} | {pm(r['diff'], r['err'], 4)} | {f(100 * r['rel'], 1)} | {f(r['pull'], 1)} |")
    L.append("")
d = J("s6_perp.json")
if d:
    L += ["## T8. Stożki prostopadłe / losowe wokół kandydatów (TRKJET2) — ε(Σ pT < 2 GeV), H06", "",
          "| klasa | pT | N_sig | własny (bez π_s): dane / MC | prostopadły: dane / MC | losowy: dane / MC | nośniki-tło: prostop. / losowy / własny | Δ_corr dane | Δ_corr MC |", "|---|---|---|---|---|---|---|---|---|"]
    for r in d["rows"]:
        L.append(f"| {r['cls']} | {r['bin']} | {f(r['Ns'], 0)} | {pm(r['own'], r['own_err'])} / {f(r['mc_own'])} | {pm(r['perp1'], r['perp1_err'])} / {f(r['mc_perp1'])} | {pm(r['r_s2'], r['r_s2_err'])} / {f(r['mc_r_s2'])} | "
                 f"{f(r['bkg_perp1'])} / {f(r['bkg_r_s2'])} / {f(r['bkg_own'])} | {pm(r['dcorr'], r['dcorr_err'])} | {f(r['dcorr_mc'])} |")
    L.append("")
open(os.path.join(OUT, "tables.md"), "w").write("\n".join(L) + "\n")
print("→ out/tables.md", len(L), "wierszy")
