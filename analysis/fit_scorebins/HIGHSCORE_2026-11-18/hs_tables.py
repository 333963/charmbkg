"""HIGHSCORE — tabele zbiorcze z out/*.json -> tables.md (+ kontrola przetrenowania na mb_nd). EKSPLORACYJNE, PO FAKCIE."""
import json, os
import numpy as np
from scipy.stats import chi2 as CHI2
import hs_common as C
R = C.R
J = lambda n: json.load(open(os.path.join(C.OUT, n)))
L = ["# HIGHSCORE — tabele (analiza eksploracyjna po fakcie, 18.11.2026)\n"]
rob = J("robust.json")
L += ["## T1. Odporność (zliczenie w oknie Δm − tło z fitu pasm; błąd = bootstrap po zdarzeniach ⊕ MC)\n",
      "R_top = f(>WP25) dane/MC (zależy od całkowitego N, czyli od tła w brudnych kategoriach); ρ50 = [N(>WP25)/N(>WP50)] dane/MC, ρ70 analogicznie (tylko czyste kategorie).\n",
      "| run | wariant | N_sig | R_top | R(WP50–25) | ρ50 | ρ70 | SF50 | f_MC(>WP25) | n_MC |", "|---|---|---|---|---|---|---|---|---|---|"]
for r in rob:
    big = r["R_err"][5] > 1
    L.append(f"| {r['run']} | {r['name']} | {r['N']:.0f} ± {r['N_err']:.0f} | " + ("nieokr." if big else f"{r['R'][5]:.2f} ± {r['R_err'][5]:.2f}") + " | "
             + ("nieokr." if big else f"{r['R'][4]:.2f} ± {r['R_err'][4]:.2f}") + f" | {r['rho50']:.3f} ± {r['rho50_err']:.3f} | {r['rho70']:.3f} ± {r['rho70_err']:.3f} | "
             + ("nieokr." if big else f"{r['SF50']:.2f} ± {r['SF50_err']:.2f}") + f" | {r['f_mc'][5]:.3f} | {r['n_mc']} |")
ptw = J("ptw.json")
L += ["\n## T2. Samo widmo pT(D0) w binie (MC przeważone do widma danych) i reszta\n",
      "| run | wagi z | biny pT | zakres | wielkość | dane/MC | przewidywanie z samego widma pT | dane/MC(pT-waż.) |", "|---|---|---|---|---|---|---|---|"]
for r in ptw:
    for k in ("8-30", "8-12", "12-30"):
        for q in ("Rtop", "rho50", "rho70", "R4", "SF50"):
            a, b, c = r[f"{k}|{q}|MC"], r[f"{k}|{q}|pred"], r[f"{k}|{q}|MCptw"]
            L.append(f"| {r['run']} | ≥WP{r['wp']} | {r['ename']} | {k} | {q} | {a[0]:.3f} ± {a[1]:.3f} | {b[0]:.3f} ± {b[1]:.3f} | {c[0]:.3f} ± {c[1]:.3f} |")
L.append("\nWagi pT (428770, ≥WP70, biny " + str(ptw[0]["edges"]) + "): " + ", ".join(f"{a:.2f}±{b:.2f}" for a, b in zip(*ptw[0]["wbin"])))
sh = J("shape.json")
L += ["\n## T3. Kształt: R_c = f_c(dane)/f_c(MC) w kwantylach MC (kategorie wg skumulowanej wydajności MC) + udział nieprompt i swap w MC\n"]
for r in sh:
    S = np.array(r["S"]); N = S.sum(); f = np.array(r["f_mc"]); fw = np.array(r["f_mc_ptw"]); w = r["wps"]
    L += [f"\n**{r['run']}, pT {r['bin'][0]}–{r['bin'][1]}** (N = {N:.0f}, n_MC = {r['n_mc']}, n_eff po wagach pT = {r['neff_ptw']:.0f}, swap/sygnał w MC = {r['x_swap_mc']:.3f})\n",
          "| kategoria (ε skum.) | S | R vs MC | R vs MC pT-waż. | udział nieprompt (MC) | rozkład swap (MC) |", "|---|---|---|---|---|---|"]
    for c in range(len(S)):
        lab = f"< {w[0]}" if c == 0 else (f"> {w[-1]}" if c == len(S) - 1 else f"{w[c-1]}–{w[c]}")
        L.append(f"| {lab} | {S[c]:.0f} ± {r['S_err'][c]:.0f} | {S[c]/N/f[c]:.2f} ± {r['S_err'][c]/N/f[c]:.2f} | {S[c]/N/fw[c]:.2f} | {r['np_frac'][c]:.3f} | {r['f_swap'][c]:.3f} |")
    L += ["\n| model (transformacja score MC) | χ²/ndf | p | parametry |", "|---|---|---|---|"]
    for k, m in r["models"].items():
        L.append(f"| {k} | {m['chi2']:.1f}/{m['ndf']} | {m['p']:.2g} | {np.round(m['par'], 3).tolist()} |")
ft = J("feats.json")
L += ["\n## T4. Cechy: sygnał danych (okno ±2 MeV − k·pasma) vs MC; χ² w 10 kwantylach MC, z_rank = przesunięcie średniej rangi [σ]; „domkn.” = przewidywane ρ50 (ρ70) po przeważeniu MC w tej jednej cesze (w regionie), względem MC nieważonego\n"]
for r in ft:
    for rg in ("all|MC", "geWP70|MC", "geWP70|MCptw", "gtWP25|MCptw"):
        v = r["regions"][rg]
        L += [f"\n**{r['run']}, pT {r['bin'][0]}–{r['bin'][1]}, region {rg}** (N_sig = {v['N_sig']:.0f}; same wagi MC → ρ50 {v['mc_w_only']['rho50']:.3f}, ρ70 {v['mc_w_only']['rho70']:.3f})\n",
              "| cecha | H? | χ²/ndf | z_rank | KS | domkn. ρ50 | domkn. ρ70 |", "|---|---|---|---|---|---|---|"]
        for f in sorted(v["feats"], key=lambda f: -v["feats"][f]["chi2"] / max(v["feats"][f]["ndf"], 1)):
            x = v["feats"][f]
            L.append(f"| {f} | {'H' if f in C.FEATS else '–'} | {x['chi2']:.1f}/{x['ndf']} | {x['z_rank']:+.1f} | {x['ks']:.3f} | {x['closure']['rho50']:.3f} | {x['closure']['rho70']:.3f} ± {x['closure_err']['rho70']:.3f} |")
        L.append(f"\nraking 2 najgorszych cech H {v['pair']}: ρ50 = {v['rake']['rho50']:.3f} ± {v['rake_err']['rho50']:.3f}, ρ70 = {v['rake']['rho70']:.3f} ± {v['rake_err']['rho70']:.3f}")
env = J("env.json")
L += ["\n## T5. MC: ułamek > WP25 / > WP50 i udział nieprompt vs ptcone40 przy ustalonym pT\n", "| pT | ptcone40 [GeV] | n | f(>WP25) | f(>WP50) | nieprompt | f(>WP25) prompt |", "|---|---|---|---|---|---|---|"]
for t in env["mc_table"]:
    L.append(f"| {t['pt'][0]}–{t['pt'][1]} | {t['cone'][0]:.0f}–{t['cone'][1]:.0f} | {t['n']} | {t['f25']:.3f} ± {t['f25_err']:.3f} | {t['f50']:.3f} | {t['np']:.3f} | {t['f25_prompt']:.3f} |")
L.append("\nNieprompt: " + json.dumps(env["np_solve"]))
L += ["\n## T6. Test: MC przeważone w pT i/lub ptcone40 vs wydajki danych w 6 kategoriach (N profilowane; wagi ustalone z nominału)\n",
      "| run | pT | MC | χ²/5 (p) | n_eff | przewid. R_top | przewid. ρ50 | ρ50 dane/przewid. | ρ70 dane/przewid. |", "|---|---|---|---|---|---|---|---|---|"]
for t in env["tests"]:
    L.append(f"| {t['run']} | {t['bin']} | {t['mc']} | {t['chi2']:.1f} ({t['p']:.2g}) | {t['neff']:.0f} | {t['pred_Rtop']:.3f} | {t['pred_rho50']:.3f} | {t['rho50_data_over_pred'][0]:.3f} ± {t['rho50_data_over_pred'][1]:.3f} | {t['rho70_data_over_pred'][0]:.3f} ± {t['rho70_data_over_pred'][1]:.3f} |")
# --- kontrola przetrenowania: prawdziwe D0<-D* z mb_nd (nieużywane w treningu H) vs fd_nd (sygnał treningowy, score poza foldem)
z = np.load(os.path.join(R.BT, "bt_mb_nd.npz"))
M = C.load("fd_nd")
L += ["\n## T7. Kontrola przetrenowania: ułamek > WP25 i > WP50 prawdziwych par w mb_nd (poza treningiem) vs fd_nd (trening, OOF)\n", "| pT | n(mb_nd) | f25 mb_nd | f25 fd_nd | f50 mb_nd | f50 fd_nd |", "|---|---|---|---|---|---|"]
for b in ("5-8", "8-30"):
    thr = C.thr_for(b, R.WP_MAIN); lo, hi = R.BINS[b]
    s = np.isin(z["truth"], [0, 1]) & np.isin(z["truth_ds"], [0, 1]) & (np.abs(z["m"] - R.M0_D0) < R.NOM_WIN) & (z["pt"] >= lo) & (z["pt"] < hi)
    sf = C.base_mask(M, b, isdata=False)
    n = s.sum(); a = (z["score"][s] >= thr[4]).mean(); a5 = (z["score"][s] >= thr[3]).mean()
    L.append(f"| {b} | {n} | {a:.3f} ± {np.sqrt(a*(1-a)/n):.3f} | {(M['score'][sf] >= thr[4]).mean():.3f} | {a5:.3f} ± {np.sqrt(a5*(1-a5)/n):.3f} | {(M['score'][sf] >= thr[3]).mean():.3f} |")
    print(L[-1], " <pT> mb, fd:", z["pt"][s].mean(), M["pt"][sf].mean())
# chi2 wąskich binów pT (rho50 vs 1)
for run in R.RUNS:
    rr = [r for r in rob if r["run"] == run and r["name"] in ("pT 8-10", "pT 10-12", "pT 12-15", "pT 15-20", "pT 20-30") and r["rho50_err"] > 0]
    c2 = sum(((r["rho50"] - 1) / r["rho50_err"]) ** 2 for r in rr)
    L.append(f"\nWąskie biny pT ({run}): χ²(ρ50 = 1) = {c2:.1f}/{len(rr)} (p = {CHI2.sf(c2, len(rr)):.2g}); wartości: " + ", ".join(f"{r['name']}: {r['rho50']:.2f} ± {r['rho50_err']:.2f}" for r in rr))
    print(L[-1])
open(os.path.join(C.HERE, "tables.md"), "w").write("\n".join(L) + "\n")
