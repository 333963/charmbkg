"""SCAN_2026-11-17: tabele z out/scan.json → tables.md. Użycie: python tables.py"""
import json
import os

import numpy as np

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(os.path.dirname(OUT), "DATA_2026-11-17")
import sys
SUF = sys.argv[1] if len(sys.argv) > 1 else ""   # "_posthoc" → wariant z odjętymi ogonami sygnału
S = json.load(open(os.path.join(OUT, "out", f"scan{SUF}.json")))
WPS = S["wps"]
WPX = [0.5, 0.7, 0.8, 0.9, 0.95]
L = []
p = L.append
pm = lambda v, e, f="{:.3f}": (f + " ± " + f).format(v, e)  # noqa: E731

p(f"# SCAN_2026-11-17 — tabele{' POST HOC (tło prawdy R2 z odjętymi ogonami sygnału)' if SUF else ''} (R1 skan cięcia, R2 pseudo-dane; {S['ntoy']} toyów/komórkę)\n")
p("Wariant fitu per komórka (reguła granicy, KRYTERIA §2): " + ", ".join(f"{k}: **{c['fitvar']}**" for k, c in S["cells"].items()) + "\n")

p("## 1. Kryteria per komórka (wariant fitu komórki)\n")
p("| run | pT | fit | b_A(WP70) H | b_A(WP70) DLZ | czuła (c) | R surowe WP70 | R_corr WP70 | nσ | (a1) | ⟨R_corr⟩_WP | Δχ²_flat surowe | Δχ²_flat skor. | (b) | rozstrz. |")
p("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for k, c in S["cells"].items():
    run, b = k.split("|")
    q = c["crit"]
    p(f"| {run} | {b} | {c['fitvar']} | {pm(*q['b70']['H'][:2])} | {pm(*q['b70']['DLZ'][:2])} | {'tak' if q['c_sensitive'] else 'NIE'} | {pm(*q['raw_R70'])} | {pm(*q['a1_R70'])} | "
      f"{q['a1_nsig']:+.1f} | {'PASS' if q['a1_pass'] else 'FAIL'} | {pm(*q['a2_cell'])} | {q['raw_dchi2']:+.1f} | {q['b_dchi2']:+.1f} | {'PASS' if q['b_pass'] else 'FAIL'} | {'tak' if q['c_decisive'] else 'nie'} |")
p("")
for run, v in S["verdict"].items():
    p(f"- **{run}: {v['verdict']}**; (a2) ⟨R_corr⟩ po WP i binach = {pm(*v['a2_mean'])} → {'PASS' if v['a2_pass'] else 'FAIL'}; czułe: {v['sensitive']}; rozstrzygające: {v['decisive']}; FAIL w: {v['fail_cells']}")
p("")

for v in ("free", "fix"):
    p(f"## 2{'a' if v == 'free' else 'b'}. R1: N_fit/ε^MC vs WP — fit {'swobodny' if v == 'free' else 'o kształcie ustalonym na MC'}\n")
    p("| run | pT | klas. | " + " | ".join(f"WP{int(w*100)}" for w in WPS) + " | WP95 (opis.) | średnia | χ²_flat/3 | nachylenie [1/ε] | po korekcie R2: χ²_flat/3 | nachylenie | repliki ok |")
    p("|---|---|---|" + "---|" * (len(WPS) + 7))
    for k, c in S["cells"].items():
        run, b = k.split("|")
        d = c["R1"][v]
        for t in ("H", "DLZ"):
            x, xc = d["der"][t], d["corr"][t]
            p(f"| {run} | {b} | {t} | " + " | ".join(f"{y:.0f} ± {e:.0f}" for y, e in zip(x["y"], x["y_err"])) + f" | {d['der']['wp95'][t]['y']:.0f} ± {d['der']['wp95'][t]['y_err']:.0f} | "
              f"{x['mean']:.0f} ± {x['mean_err']:.0f} | {x['chi2_flat']:.1f} | {pm(x['slope'], x['slope_err'], '{:+.2f}')} | {xc['chi2_flat']:.1f} | {pm(xc['slope'], xc['slope_err'], '{:+.2f}')} | {d['n_rep_ok']} |")
    p("")
    p("Stosunek R = (N/ε)_DLZ/(N/ε)_H, surowy → po korekcie R2 (wariant tła A; błąd zawiera stat. toyów i |A − C|):\n")
    p("| run | pT | " + " | ".join(f"WP{int(w*100)}" for w in WPS) + " | ⟨R⟩ surowe | ⟨R_corr⟩ |")
    p("|---|---|" + "---|" * (len(WPS) + 2))
    for k, c in S["cells"].items():
        run, b = k.split("|")
        r, rc = c["R1"][v]["der"]["R"], c["R1"][v]["corr"]["R"]
        p(f"| {run} | {b} | " + " | ".join(f"{a:.3f} ± {ae:.3f} → {x:.3f} ± {xe:.3f}" for a, ae, x, xe in zip(r["R"], r["R_err"], rc["R"], rc["R_err"])) + f" | {pm(r['mean'], r['mean_err'])} | {pm(rc['mean'], rc['mean_err'])} |")
    p("")

p("## 3. Jakość fitów w danych (fit swobodny / kształt ust.): χ²/ndf, na granicy, S/B (okno ±45 MeV)\n")
p("| run | pT | klas. | " + " | ".join(f"WP{int(w*100)}" for w in WPX) + " |")
p("|---|---|---|" + "---|" * len(WPX))
for k, c in S["cells"].items():
    run, b = k.split("|")
    for ti, t in enumerate(("H", "DLZ")):
        fr, fx = c["R1"]["free"], c["R1"]["fix"]
        cells = []
        for i, w in enumerate(WPX):
            j = ti * len(WPX) + i
            cells.append(f"{fr['chi2_ndf'][j]:.2f}{'*' if c['at_limit_free'][t][str(w)] else ''}{'' if fr['valid'][j] else '!'} / {fx['chi2_ndf'][j]:.2f}{'' if fx['valid'][j] else '!'}; S/B {fr['SB'][j]:.3f}; rep ok {fr['boot_ok_frac'][j]:.2f}")
        p(f"| {run} | {b} | {t} | " + " | ".join(cells) + " |")
p("\n`*` = fit swobodny na granicy parametru kształtu; `!` = MIGRAD nieważny; „rep ok” = udział replik bootstrapu ze zbieżnym fitem swobodnym nie na granicy.\n")

p("## 4. Kontrola błędu różnicy selekcji zagnieżdżonych, N(WP90) − N(WP70) (wariant fitu komórki)\n")
p("| run | pT | klas. | σ bootstrap | √\\|σ²_90 − σ²_70\\| | √(σ²_90 + σ²_70) (nieskorel.) |")
p("|---|---|---|---|---|---|")
for k, c in S["cells"].items():
    run, b = k.split("|")
    for t in ("H", "DLZ"):
        n = c["R1"][c["fitvar"]]["der"][t]["nested_check_90_70"]
        p(f"| {run} | {b} | {t} | {n['sig_diff_boot']:.0f} | {n['sig_diff_approx']:.0f} | {n['sig_diff_uncorr']:.0f} |")
p("")

for v in ("free", "fix"):
    p(f"## 5{'a' if v == 'free' else 'b'}. R2: obciążenie b = ⟨N_fit/N_true⟩ — fit {'swobodny' if v == 'free' else 'o kształcie ustalonym'}\n")
    p("Format: tło A: średnia ± błąd (mediana; rozrzut; udział udanych; średni pull) // tło C: średnia // Asimov A // P (tło = exp(pol2) z fitu danych, kontrola maszynerii).\n")
    p("| run | pT | klas. | S/B prawdy WP70 | " + " | ".join(f"WP{int(w*100)}" for w in WPX) + " |")
    p("|---|---|---|---|" + "---|" * len(WPX))
    for k, c in S["cells"].items():
        run, b = k.split("|")
        for t in ("H", "DLZ"):
            cells = []
            for w in WPX:
                A, C_, P = (c["R2"][f"{t}|{w}|{x}"][v] for x in "ACP")
                cells.append(f"{pm(A['b'], A['b_err'])} ({A['med']:.3f}; {A['rms']:.3f}; {A['f_ok']:.2f}; {A['pull_mean']:+.2f}) // {C_['b']:.3f} // {A['asimov']:.3f} // {P['b']:.3f}")
            p(f"| {run} | {b} | {t} | {c['R2'][f'{t}|0.7|A']['SB_true']:.3f} | " + " | ".join(cells) + " |")
    p("")

p("## 6. Tło prawdy R2: jakość exp(pol4) w pasmach bocznych i tło pod pikiem (±45 MeV)\n")
p("| run | pT | klas. | WP | χ²/ndf SB: exp(pol2) z fitu danych | A: exp(pol4), SB>45 | C: exp(pol4), SB>70 | B_win A | B_win C | B_win z fitu danych | (A−fit)/fit | S0·ε |")
p("|---|---|---|---|---|---|---|---|---|---|---|---|")
for k, c in S["cells"].items():
    run, b = k.split("|")
    for t in ("H", "DLZ"):
        for w in WPX:
            q = c["R2"][f"{t}|{w}|truthbkg"]
            p(f"| {run} | {b} | {t} | {int(w*100)} | {q['pol2_chi2sb_data']:.2f} | {q['A_chi2sb']:.2f}{'' if q['A_valid'] else '!'} | {q['C_chi2sb']:.2f}{'' if q['C_valid'] else '!'} | {q['B_win_A']:.0f} | {q['B_win_C']:.0f} | {q['B_win_fit']:.0f} | "
              f"{(q['B_win_A'] - q['B_win_fit']) / q['B_win_fit']:+.4f} | {c['R2'][f'{t}|{w}|A']['N_true']:.0f} |")
p("")

# domknięcie MC mb (istniejące, opisowo): N_fit/N_truth per WP
mc = json.load(open(os.path.join(DATA, "runs", "lxplus", "ana", "mc_checks.json")))["closure"]
p("## 7. Opisowo: domknięcie bezwzględne na MC mb z istniejącego `mc_checks.json` (prawdziwi kandydaci sygnału, tło MC): N_fit/N_truth\n")
p("| klas. | pT | WP80 | WP70 | WP50 | WP90 (odniesienie) |")
p("|---|---|---|---|---|---|")
for t, tn in (("H", "H_s12345"), ("DLZ", "HDLZ_s12345")):
    for b in ("3-5", "5-8", "8-30"):
        row = []
        for w in (0.8, 0.7, 0.5):
            q = mc[f"{tn}|{b}|{w}"]
            row.append(f"{q['N_fit'] / max(q['N_truth'], 1):.2f} ± {q['N_fit_err'] / max(q['N_truth'], 1):.2f} (N_truth {q['N_truth']})")
        q = mc[f"{tn}|{b}|0.7"]
        row.append(f"{q['N_fit90'] / max(q['N_truth90'], 1):.2f} (N_truth {q['N_truth90']})")
        p(f"| {t} | {b} | " + " | ".join(row) + " |")
open(os.path.join(OUT, f"tables{SUF}.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L[:14]))
