#!/usr/bin/env python3
"""TAGGED — tabele wynikowe (markdown) z out/unblind.json, out/prep.json, out/dst_unblind.json. Bez fitów."""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tg_lib as L  # noqa: E402

U = json.load(open(os.path.join(L.OUT, "unblind.json")))["items"]
D = json.load(open(os.path.join(L.OUT, "dst_unblind.json")))
o = []
pc = lambda x: f"{100 * x:.3f}"  # noqa: E731
o.append("## T1. Uzyski N_tag (kształt stały z fitu inkluzywnego; przedział 68 % z profilu; Z = √q₀ lokalne)\n")
o.append("| gat. | run | region | pT | WP | sel. | rząd tła | kand. w zakresie / w oknie | tło w oknie (fit) | N_tag | −stat | +stat | syst kształtu | przedział 95 % | Z | p (zabawki) | skala błędu | N/L [nb] |")
o.append("|" + "---|" * 18)
for k, r in U.items():
    f = r["fit"]
    o.append(f"| {r['sp']} | {r['run']} | {r['reg']} | {k.split('|')[3]} | {int(100 * r['eps'])} | {r['xsel']} | {r['nbkg']} | {r['n_tot']:.0f} / {r['n_win']:.0f} | {r['bkg_win']:.1f} | **{f['N']:.1f}** | {r['stat_lo']:.1f} | {r['stat_hi']:.1f} | "
             f"{(r.get('syst') or {}).get('razem', float('nan')):.1f} | [{f['lo95']:.1f}; {f['hi95']:.1f}] | {f['Z']:.2f} | {r.get('p_toys', float('nan')):.4f} | {r['err_scale']:.2f} | {r.get('N_per_nb', float('nan')):.1f} |")
o.append("\n## T2. Ułamek tagowany (tylko regiony z mianownikiem inkluzywnym)\n")
o.append("| gat. | run | region | pT | WP | sel. | N_incl ± (HESSE × 1,3) | P_tag,SB [%] | f_raw [%] | N_acc | f_true [%] | −stat | +stat | syst (kształt ⊕ f_acc ⊕ N_incl) | przedział 68 % | przedział 95 % |")
o.append("|" + "---|" * 16)
for k, r in U.items():
    if "f_raw" not in r:
        continue
    o.append(f"| {r['sp']} | {r['run']} | {r['reg']} | {k.split('|')[3]} | {int(100 * r['eps'])} | {r['xsel']} | {r['N_incl']:.0f} ± {r['N_incl_err']:.0f} | {pc(r['P_tag_SB'])} | {pc(r['f_raw'])} | {r['N_acc']:.2f} ± {r['N_acc_err']:.2f} | **{pc(r['f_true'])}** | "
             f"{pc(r['f_stat_lo'])} | {pc(r['f_stat_hi'])} | {pc(r['f_syst'])} | [{pc(r['f_true'] - r['f_true_lo68'])}; {pc(r['f_true'] + r['f_true_hi68'])}] | [{pc(r['f_true'] - r['f_true_lo95'])}; {pc(r['f_true'] + r['f_true_hi95'])}] |")
o.append("\n## T3. Porównanie z przewidywaniami (427929 MBTS, T2)\n")
o.append("| pozycja | przewidywanie | wartość [%] (zakres) | werdykt wg kryteriów | (f_true − P)/σ | P centralne w przedziale 68 % |")
o.append("|---|---|---|---|---|---|")
for k, r in U.items():
    for nm, p in (r.get("pred") or {}).items():
        o.append(f"| {k} | {'pasma (S)' if nm == 'S' else 'Pythia (M)'} | {pc(p['central'])} ({pc(p['lo'])}–{pc(p['hi'])}) | {p['werdykt']} | {p['pull_central']:+.2f} | {'tak' if p['central_in_68'] else 'nie'} |")
o.append("\n## T4. Systematyka kształtu, ΔN per źródło\n")
o.append("| pozycja | sc, dmu | wagi pT | rząd tła + 1 | zakres alt. | odbicie, κ | razem |")
o.append("|---|---|---|---|---|---|---|")
for k, r in U.items():
    s = r.get("syst")
    if s:
        g = lambda n: f"{s[n]:.1f}" if n in s else "–"  # noqa: E731
        o.append(f"| {k} | {g('ksztalt_sc_dmu')} | {g('wagi_pT')} | {g('rzad_tla')} | {g('zakres')} | {g('odbicie_kappa')} | {s['razem']:.1f} |")
o.append("\n## T5. Zabawki (przed odczytem) i K4\n")
o.append("| pozycja | N wstrz. | n | obciążenie | rms N̂ | ⟨pull⟩ / rms | pokrycie 68 % | mediana Z | udział Z ≥ 3 |")
o.append("|---|---|---|---|---|---|---|---|---|")
for k, r in U.items():
    for t in r.get("toys") or []:
        o.append(f"| {k} | {t['ninj']:.0f} | {t['n']} | {t['bias']:+.2f} ± {t['bias_err']:.2f} | {t['rmsN']:.1f} | {t['pull_mean']:+.2f} / {t['pull_rms']:.2f} | {t['cov68']:.3f} | {t['medZ']:.2f} | {t['fracZ3']:.3f} |")
    if "K4" in r:
        t = r["K4"]
        o.append(f"| {k} **K4 (tag losowy, oczekiwane {t['expected']:.1f})** | {t['expected']:.1f} | {t['n']} | {t['bias']:+.2f} ± {t['bias_err']:.2f} | {t['rmsN']:.1f} | {t['pull_mean']:+.2f} / {t['pull_rms']:.2f} | {t['cov68']:.3f} | {t['medZ']:.2f} | {t['fracZ3']:.3f} |")
o.append("\n## T6. D* (K2; kandydaci slim bez klasyfikatora)\n")
o.append("| pozycja | N_incl ± HESSE | pasma tag / incl | P_tag,SB [%] | N_tag | przedział 68 % | przedział 95 % | Z | f_raw [%] |")
o.append("|---|---|---|---|---|---|---|---|---|")
for k, r in D.items():
    f = r["fit"]
    o.append(f"| {k} | {r['N_incl']:.0f} ± {r['N_incl_err']:.0f} | {r['nsb_tag']:.0f} / {r['nsb_incl']:.0f} | {pc(r['P_tag_SB'])} | {f['N']:.1f} | [{f['lo68']:.1f}; {f['hi68']:.1f}] | [{f['lo95']:.1f}; {f['hi95']:.1f}] | {f['Z']:.2f} | {pc(f['N'] / r['N_incl'])} |")
# kontrolki pochodne
M = lambda sp, x="3-30", e="0.7", reg="MBTS", s="T2": U[f"{sp}|427929|{reg}|{x}|{e}|{s}"]  # noqa: E731
o.append("\n## T7. Kontrolki pochodne\n")
a, b = M("D0"), M("Dplus")
sa = (a["f_true_lo68"] + a["f_true_hi68"]) / 2; sb = (b["f_true_lo68"] + b["f_true_hi68"]) / 2
rat = b["f_true"] / a["f_true"]; er = rat * math.hypot(sa / a["f_true"], sb / b["f_true"])
w = (a["f_true"] / sa ** 2 + b["f_true"] / sb ** 2) / (1 / sa ** 2 + 1 / sb ** 2); ew = 1 / math.sqrt(1 / sa ** 2 + 1 / sb ** 2)
o.append(f"- K1: f_true(D⁺)/f_true(D⁰) = {rat:.2f} ± {er:.2f} (odchylenie od 1: {(rat - 1) / er:+.2f} σ); różnica f(D⁰) − f(D⁺) = {pc(a['f_true'] - b['f_true'])} ± {pc(math.hypot(sa, sb))} % ({(a['f_true'] - b['f_true']) / math.hypot(sa, sb):+.2f} σ); "
         f"ułamek wspólny D⁰ + D⁺ = {pc(w)} ± {pc(ew)} %.")
for sp in L.SPN:
    n = M(sp)
    sn = (n["f_stat_lo"] + n["f_stat_hi"]) / 2
    for e in ("0.5", "0.8", "0.9"):
        v = M(sp, e=e)
        sv = (v["f_stat_lo"] + v["f_stat_hi"]) / 2
        sd = math.sqrt(abs(sv ** 2 - sn ** 2))
        o.append(f"- K3 {sp} WP{int(float(e) * 100)}: f_raw = {pc(v['f_raw'])} ± {pc(sv)} % wobec WP70 {pc(n['f_raw'])} ± {pc(sn)} %; Δ = {pc(v['f_raw'] - n['f_raw'])}, σ_Δ = {pc(sd)} → {abs(v['f_raw'] - n['f_raw']) / sd:.2f} σ_Δ ({'OK' if abs(v['f_raw'] - n['f_raw']) <= 2 * sd else 'NIESPEŁNIONE'}).")
for sp in L.SPN:
    for run, reg in (("427929", "MBTS"), ("428770", "J12")):
        c = U[f"{sp}|{run}|{reg}|3-30|0.7|CR"]
        f = c["fit"]
        exp = c["N_incl"] * c["P_tag_SB"]
        o.append(f"- K5 {sp} {run} {reg} (CR, tag przypadkowy z konstrukcji): N = {f['N']:.1f} [{f['lo68']:.1f}; {f['hi68']:.1f}], oczekiwane N_incl · P_CR,SB = {exp:.2f} (pasma CR: {c['nsb_tag']:.0f} kandydatów).")
for sp in L.SPN:
    x, y = U[f"{sp}|427929|MBAFP|3-30|0.7|T2"], U[f"{sp}|428770|R1|3-30|0.7|T2"]
    sx, sy = (x["stat_lo"] + x["stat_hi"]) / 2 / L.LUMI["427929"], (y["stat_lo"] + y["stat_hi"]) / 2 / L.LUMI["R1"]
    d = x["N_per_nb"] - y["N_per_nb"]
    o.append(f"- K6 {sp}: N_tag/L 427929 (MBTS ∧ bit noalg_L1AFP) = {x['N_per_nb']:.0f} ± {sx:.0f} /nb; 428770 R1 = {y['N_per_nb']:.0f} ± {sy:.0f} /nb; różnica {d:.0f} ± {math.hypot(sx, sy):.0f} ({d / math.hypot(sx, sy):+.2f} σ, tylko stat).")
x, y = U["D0|427929|MBTS_A|3-30|0.7|T2"], U["D0|427929|MBTS_C|3-30|0.7|T2"]
o.append(f"- K7 D⁰ strony: tylko A {x['fit']['N']:.1f} −{x['stat_lo']:.1f} +{x['stat_hi']:.1f}; tylko C {y['fit']['N']:.1f} −{y['stat_lo']:.1f} +{y['stat_hi']:.1f}.")
open(os.path.join(L.TOP, "tabele", "wyniki.md"), "w").write("# TAGGED — tabele wynikowe (generowane przez skrypty/t4_tables.py)\n\n" + "\n".join(o) + "\n")
print("\n".join(o))
