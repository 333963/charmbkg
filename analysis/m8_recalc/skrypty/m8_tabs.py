#!/usr/bin/env python3
"""Tabele pkt 2 i 3 z out/d23_m8.json, out/R_m8.json → tabele/d23_m8.md, tabele/R_m8.md"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "..", "out", "d23_m8.json")))
HEAD = {("427929", "drzewo"): 0.1, ("427929", "siec"): 0.2, ("428770", "drzewo"): 0.1, ("428770", "siec"): 0.1}
LAB = {"S4": "S4 etap 4 (kształt stały, 1,75–2,00) = STARY", "M0": "M0 (kształt wolny, 1,75–2,00)", "M8n": "M8 bez składowej DMM+CM (Gauss, 1,70–2,05)",
       "M8a": "M8 + szablon 3–5 GeV, κ wg przepisu (1,154 × 1,364) — WARIANT", "M8b": "M8 + szablon 3–5 GeV, κ z prawdy mb_nd 3–5 przy tym ε (× 1,364) — WARIANT"}
L = ["# M8_przeliczenie — pkt 2: pik D⁰→Kπ w 2–3 GeV; post hoc", "",
     "W MC nie ma tła dla 2–3 GeV (mb_nd tylko pT ≥ 3) → **szablon DMM+CM dla tego binu nie istnieje**. M8a/M8b używają szablonu i κ z 3–5 GeV — warianty, nie nominał.",
     f"κ_T(3–5 GeV, ε ≤ 0,5) = {d['kappa_T_3-5_eps0.5']:.3f}; κ z prawdy mb_nd 3–5 GeV: " + "; ".join(f"{k.replace('|', ' WP').replace('0.', '')}0: {v['kap']:.2f}" for k, v in d["kappa_truth_3-5"].items()) + ".",
     "σ = max(HESSE, bootstrap po zdarzeniach; 100 replik w punktach z werdyktu etapu 4). Z_LR = √Δ(−2 ln L) względem N = 0 (ten sam zakres i tło).", "",
     "## D1. Punkty z werdyktu etapu 4", "", "| run | klasyfikator, WP | model | N ± σ | Δ wzgl. S4 [σ_S4] | N/σ | Z_LR | p(χ²) | dmu [MeV]; sc | N/ε_MC | uwagi |", "|---|---|---|---|---|---|---|---|---|---|---|"]
for (r, n), wp in HEAD.items():
    o = d["wp"][f"{r}|{n}|{wp}"]
    for m, v in o.items():
        dl = (v["N"] - o["S4"]["N"]) / o["S4"]["err"]
        L.append(f"| {r} | {n} WP{int(wp * 100)} | {LAB[m]} | {v['N']:.0f} ± {v['err']:.0f} | {dl:+.2f} | {v['Z']:.1f} | {v['ZLR']:.1f} | {v['p']:.2g} | {1e3 * v['dmu']:+.1f}; {v['sc']:.2f} | {v['y']:.0f} ± {v['yerr']:.0f} | "
                 f"{'' if v['valid'] else 'fit formalnie nieważny (EDM); '}{'granica' if v['at_limit'] else ''} |")
L += ["", "## D2. N/ε_MC per punkt pracy (skumulowane WP30 / WP20 / WP10)", "", "| run | klasyfikator | model | WP30 | WP20 | WP10 |", "|---|---|---|---|---|---|"]
for r in ("427929", "428770"):
    for n in ("drzewo", "siec"):
        for m in LAB:
            L.append(f"| {r} | {n} | {m} | " + " | ".join(f"{d['wp'][f'{r}|{n}|{w}'][m]['y']:.0f} ± {d['wp'][f'{r}|{n}|{w}'][m]['yerr']:.0f}" for w in (0.3, 0.2, 0.1)) + " |")
L += ["", "## D3. ⟨N/ε⟩ z trzech rozłącznych wycinków ε < 0,3 (jak test 5 D23) i stosunek do y(D⁺) = 3441 ± 260 (427929), 3055 ± 363 (428770); oczekiwane 0,82–1,03", "",
      "| run | klasyfikator | model | ⟨y⟩ | χ²/2 | y(D⁰)/y(D⁺) | nadwyżka wobec środka oczekiwania 0,925 |", "|---|---|---|---|---|---|---|"]
for r in ("427929", "428770"):
    for n in ("drzewo", "siec"):
        for m in LAB:
            v = d["slice"][f"{r}|{n}|mean|{m}"]
            L.append(f"| {r} | {n} | {m} | {v['y']:.0f} ± {v['yerr']:.0f} | {v['chi2']:.1f} | {v['ratio_Dp']:.2f} ± {v['ratio_Dp_err']:.2f} | × {v['ratio_Dp'] / 0.925:.2f} |")
open(os.path.join(HERE, "..", "tabele", "d23_m8.md"), "w").write("\n".join(L) + "\n")
R = json.load(open(os.path.join(HERE, "..", "out", "R_m8.json")))["res"]
L = ["# M8_przeliczenie — pkt 3: metoda R (sieć × drzewo, kategorie b/n/t) dla D⁰ — domknięcie na mb_nd; post hoc", "",
     "Histogramów łącznych lokalnie nie było; zbudowane z plików score per kandydat pobranych z EOS (`dnet4/app/d150`, tylko fd_nd i mb_nd, 241 MB; bez joba).",
     "M0 odtwarza etap 4 co do cyfry (R, HESSE, R_prawda); rms bootstrapu różni się w granicach szumu 100 replik (inne ziarno).",
     "M8 dla kategorii łącznych = rozszerzenie: κ_c z prawdy mb_nd per kategoria (przepis M8 definiuje κ tylko dla drzewa); M8n = to samo bez składowej DMM+CM.",
     "Bramka etapu 4: |R_fit − R_prawda| ≤ 2,5σ i ≤ 0,04 oraz σ ≤ 0,05; σ = max(HESSE, bootstrap 100 replik po zdarzeniach).", "",
     "| ziarno | pT | WP | R_prawda | κ b / n / t | M0: R_fit ± σ (odchylenie) [bramka] | M8n | M8 | G_fit(MC): M0 / M8 | p(χ²) M8 |", "|---|---|---|---|---|---|---|---|---|---|"]
for r in R:
    f = lambda m: f"{r[m]['R']:.3f} ± {r[m]['R_err']:.3f} ({r[m]['bias']:+.3f}) [{'TAK' if r[m]['pass'] else 'NIE'}]"  # noqa: E731
    L.append(f"| {r['tag'][-2:]} | {r['bin']} | {int(r['wp'] * 100)}{'' if r['wp'] == 0.7 else ' (opisowo, HESSE)'} | {r['R_truth']:.3f} ± {r['R_truth_err']:.3f} | " + " / ".join(f"{k:.2f}" for k in r["kap"]) +
             f" | {f('M0')} | {f('M8n')} | {f('M8')} | {r['M0']['G_fit']:.2f} / {r['M8']['G_fit']:.2f} | {r['M8']['p']:.2g} |")
open(os.path.join(HERE, "..", "tabele", "R_m8.md"), "w").write("\n".join(L) + "\n")
print(open(os.path.join(HERE, "..", "tabele", "d23_m8.md")).read()[:3000])
