#!/usr/bin/env python3
"""ISO_origin: tabele markdown z out/s1.json, out/s2.json → out/tables.md."""
import json, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
s1 = json.load(open(f"{H}/out/s1.json")); s2 = json.load(open(f"{H}/out/s2.json"))
L = []
f = lambda x, n=3: "–" if x is None or not np.isfinite(x) else f"{x:.{n}f}"  # noqa
def pm(r, a="eps", b="err"): return "–" if r is None or r.get(a) is None else f"{f(r[a])} ± {f(r.get(b))}"
def flag(r):
    S = r["S"]; sig = S / r["S_err"] if r["S_err"] else 0
    return "niemiarodajny" if (S < 50 or sig < 3) else ("słaby" if S < 150 else "")
L += ["## T1. Klasy triggera (ptcone40 < 2 GeV; score_H > 0,6; |m−m0| < 40 MeV)", "",
      "| próbka | pT | N_kand | N_sig (SBS) | ε SBS | ε fit pass/fail (profil) | N_sig fit | det | ε_MC bez wag | tło f (score>0,6) | tło f (bez score), N | uwaga |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for k, r in s1["classes"].items():
    nm, pb = k.split("|"); ft = r.get("fit")
    L.append(f"| {nm} | {pb} | {r['ncand']} | {r['S']:.0f} ± {r['S_err']:.0f} | {pm(r)} | {pm(ft) if ft else '–'} | {ft['N']:.0f} ± {ft['N_err']:.0f} |" .replace("| – | ", "| – | ") if ft else
             f"| {nm} | {pb} | {r['ncand']} | {r['S']:.0f} ± {r['S_err']:.0f} | {pm(r)} | – | – |")
    L[-1] += f" {('tak' if ft['det'] else 'NIE') if ft else '–'} | {f(r['mc_eps'])} | {f(r['bkg']['f'])} ± {f(r['bkg']['err'])} | {f(r['bkg_noscore']['f'],4)} ± {f(r['bkg_noscore']['err'],4)} ({r['bkg_noscore']['n']}) | {flag(r)} |"
L += ["", "## T2. ε_iso różniczkowo (dane; SBS, fit gdzie N_sig ≥ 150)", "", "| próbka | zmienna | bin | N_sig | ε SBS | ε fit | ε_MC (bin, bez wag) | tło f | uwaga |", "|---|---|---|---|---|---|---|---|---|"]
for k, r in s1["diff"].items():
    nm, var, lo, hi = k.split("|"); ft = r.get("fit"); mc = s1["mc_diff"].get(f"{var}|{lo}|{hi}")
    L.append(f"| {nm} | {var} | {lo}–{hi} | {r['S']:.0f} ± {r['S_err']:.0f} | {pm(r)} | {pm(ft) + ('' if ft['det'] else ' (nieokr.)') if ft else '–'} | {f(mc['eps']) if mc else '–'} | {f(r['bkg_noscore']['f'])} | {flag(r)} |")
L += ["", "## T3. Pile-up na tle (pT 3–5, Δm > 150 MeV, bez score)", "", "| próbka | zmienna | bin | N | f(ptcone40<2) | ⟨μ⟩ | med nSelBL | udział GRL |", "|---|---|---|---|---|---|---|---|"]
for k, r in s1["pileup_bkg"].items():
    a, var, lo, hi = k.split("|"); L.append(f"| {a} | {var} | {lo}–{hi} | {r['n']} | {r['f']:.4f} ± {r['err']:.4f} | {r['mu_mean']:.4g} | {r['nSelBL_med']:.0f} | {r['grl_frac']:.2f} |")
fo = s1["forensic"]
L += ["", "## T4. Bity triggera w tabelach zdarzeń", "", f"428770: {fo['n_events_428770']} zdarzeń, bez żadnego zapisanego łańcucha: {fo['none_of_stored']}; P(bit | J12) w 428770: {fo['P_in_J12_events_428770']}; ⟨1/PS_j20⟩ = {fo['mean_inv_ps_j20']:.4f}", "",
      "| łańcuch | N (428770) | sam (bez innego bitu) | PS | ⟨nSelBL⟩ 428770 | med lead 428770 | N (427929 GRL) | 427929: ułamek bez innego twardego łańcucha | ⟨nSelBL⟩ 427929 | med lead 427929 |", "|---|---|---|---|---|---|---|---|---|---|"]
for k, d in fo["chains"].items():
    L.append(f"| {k} | {d['n8']} | {d['alone8']} ({d['alone8']/d['n8']:.3f}) | {d['ps8']:.0f} | {d['nSelBL_mean8']:.1f} | {d['lead_med8']:.2f} | {d.get('n9','–')} | {f(d.get('frac_wo_other_hard9'))} | {f(d.get('nSelBL_mean9'),1)} | {f(d.get('lead_med9'),2)} |")
L += ["", "Domknięcie emulacji (ułamek MBTSbased w J12only ∪ MBTSbased per komórka (nSelBL, lead_pt)):", ""]
for k, d in fo["closure"].items():
    L.append(f"- {k}: χ²/ndf = {d['chi2']:.0f}/{d['ndf']}, max |pred/obs − 1| = {d['max_rel_dev']:.2f}, ułamek przewidziany {d['pred_frac_MBTS']:.3f} vs obserwowany {d['obs_frac_MBTS']:.3f}")
n = s2["nominal"]
L += ["", "## T5. Uniwersalność ε(zmienna): χ² między próbkami (sygnał, ptcone40 < 2)", "", "| zmienna | 427929-MB vs J12only | 427929-MB vs MBTSbased | J12only vs MBTSbased | 427929-MB vs MC | J12only vs MC | MBTSbased vs MC | TŁO: 9MB vs J12 | TŁO: 9MB vs MBTSb. | TŁO: J12 vs MBTSb. |", "|---|---|---|---|---|---|---|---|---|---|"]
for var, u in n["uni"].items():
    c = u["chi2"]; g = lambda k: f"{c[k]['chi2']:.1f}/{c[k]['ndf']}" + (f" (p={c[k]['p']:.2g})" if c[k].get("p") is not None else "")  # noqa
    L.append(f"| {var} | " + " | ".join(g(k) for k in ("9:MB vs 8:J12only", "9:MB vs 8:MBTSbased", "8:J12only vs 8:MBTSbased", "9:MB vs MC", "8:J12only vs MC", "8:MBTSbased vs MC", "TŁO 9:MB vs 8:J12only", "TŁO 9:MB vs 8:MBTSbased", "TŁO 8:J12only vs 8:MBTSbased")) + " |")
for var, u in n["uni"].items():
    L += ["", f"ε({var}):", "", "| bin | MC | 427929 MB (N_sig) | 428770 J12only (N_sig) | 428770 MBTSbased (N_sig) | tło 9MB / J12 / MBTSb. |", "|---|---|---|---|---|---|"]
    t = u["table"]
    for i, r in enumerate(t["MC"]):
        cell = lambda x: f"{pm(x)} ({x['S']:.0f})"  # noqa
        L.append(f"| {r['lo']:g}–{r['hi']:g} | {f(r['eps'])} | {cell(t['9:MB'][i])} | {cell(t['8:J12only'][i])} | {cell(t['8:MBTSbased'][i])} | {f(t['9:MB'][i]['bkg_f'])} / {f(t['8:J12only'][i]['bkg_f'])} / {f(t['8:MBTSbased'][i]['bkg_f'])} |")
L += ["", "## T6. Ten sam trigger w obu runach", ""]
for k, rows in n["same_trig"].items():
    a, b = k.split(" || ")
    L += [f"**{a}  vs  {b}**", "", "| pT | 427929: ε sygn. (N_sig) | 428770: ε sygn. (N_sig) | Δ sygn. [σ] | 427929: tło f (N) | 428770: tło f (N) | Δ tła [σ] |", "|---|---|---|---|---|---|---|"]
    for r in rows:
        x, y = r[a], r[b]
        ds = (x["eps"] - y["eps"]) / np.hypot(x["err"], y["err"]) if x["eps"] is not None and y["eps"] is not None and x["err"] and y["err"] else np.nan
        db = (x["bkg_f"] - y["bkg_f"]) / np.hypot(x["bkg_err"], y["bkg_err"]) if x["bkg_n"] and y["bkg_n"] else np.nan
        L.append(f"| {r['pt']} | {pm(x)} ({x['S']:.0f}) | {pm(y)} ({y['S']:.0f}) | {f(ds,1)} | {f(x['bkg_f'],4)} ± {f(x['bkg_err'],4)} ({x['bkg_n']}) | {f(y['bkg_f'],4)} ± {f(y['bkg_err'],4)} ({y['bkg_n']}) | {f(db,1)} |")
    L.append("")
SCH = list(dict.fromkeys(k.split("|")[3] for k in n["sf"]))
CL = ["9:MB", "9:MB&j20&!J12", "9:MB&!j20&!J12", "8:all", "8:J12only", "8:MBTSbased", "8:other"]
for cut in ("ptcone40<1", "ptcone40<2", "ptcone40<4", "rel40<0.5"):
    L += [f"## T7. SF po przeważeniu MC do rozkładu sygnału w danych — {cut}, pT 3–30", "", "| schemat | " + " | ".join(CL) + " | χ² vs 1 (4 kl. rozł.) | χ² spójności (p) |", "|---|" + "---|" * (len(CL) + 2)]
    for s in SCH:
        cl = n["closure"].get(f"{cut}|{s}")
        L.append(f"| {s} | " + " | ".join(pm(n["sf"].get(f"{cut}|{c}|3-30|{s}"), "SF") for c in CL) + f" | {cl['chi2_vs1']:.1f}/{cl['ndf']} | {cl['chi2_spojnosc']:.1f}/3 ({cl['p_spojnosc']:.2g}) |")
    L.append("")
L += ["## T8. SF per bin pT (ptcone40 < 2): brak → pT×lead_pt×nSelBL; ułamek sygnału w komórkach z cofnięciem", "", "| klasa | pT | N_sig | ε_dane | ε_MC bez wag | SF bez wag | ε_MC rw | SF rw | cofnięcie |", "|---|---|---|---|---|---|---|---|---|"]
for c in CL:
    for pb in ("3-30", "3-5", "5-8", "8-30"):
        a, b = n["sf"].get(f"ptcone40<2|{c}|{pb}|brak"), n["sf"].get(f"ptcone40<2|{c}|{pb}|pT×lead_pt×nSelBL")
        if a and b:
            L.append(f"| {c} | {pb} | {a['S']:.0f} | {f(a['eps_data'])} | {f(a['eps_mc_rw'])} | {pm(a,'SF')} | {f(b['eps_mc_rw'])} | {pm(b,'SF')} | {b['frac_fallback']:.2f} |")
L += ["", "## T9. Sprawdzenia krzyżowe (ptcone40 < 2, pT 3–30): ε_dane (SBS) i SF", "", "| wariant | klasa | N_sig | ε_dane | SF bez wag | SF nSelBL | SF pT×lead_pt | SF pT×lead_pt×nSelBL |", "|---|---|---|---|---|---|---|---|"]
for tag in ("nominal", "LB_parzyste", "LB_nieparzyste", "okno_m_25MeV", "okno_m_55MeV", "score>0.0", "score>0.3", "score>0.8"):
    for c in ("9:MB", "8:J12only", "8:MBTSbased", "8:all"):
        g = lambda s: s2[tag]["sf"].get(f"ptcone40<2|{c}|3-30|{s}")  # noqa
        a = g("brak")
        if a:
            L.append(f"| {tag} | {c} | {a['S']:.0f} | {f(a['eps_data'])} ± {f(a['err_data']*a['eps_mc_rw'])} | {pm(a,'SF')} | {pm(g('nSelBL'),'SF')} | {pm(g('pT×lead_pt'),'SF')} | {pm(g('pT×lead_pt×nSelBL'),'SF')} |")
L += ["", "Ten sam trigger (j20 ∧ ¬J12), pT 3–30, w wariantach:", "", "| wariant | 427929 ε sygn. | 428770 ε sygn. | 427929 tło | 428770 tło |", "|---|---|---|---|---|"]
for tag in ("nominal", "LB_parzyste", "LB_nieparzyste", "okno_m_25MeV", "okno_m_55MeV", "score>0.0", "score>0.3", "score>0.8"):
    r = s2[tag]["same_trig"]["9:MB&j20&!J12 || 8:j20&!J12"][0]; x, y = r["9:MB&j20&!J12"], r["8:j20&!J12"]
    L.append(f"| {tag} | {pm(x)} ({x['S']:.0f}) | {pm(y)} ({y['S']:.0f}) | {f(x['bkg_f'],4)} | {f(y['bkg_f'],4)} |")
d = s2["sig_dist"]
L += ["", "## T10. Rozkład sygnału (ułamki, wagi SBS) vs MC", ""]
for var, dd in d.items():
    ed = dd["edges"]; L += [f"{var}:", "", "| próbka | " + " | ".join(f"{a:g}–{b:g}" for a, b in zip(ed[:-1], ed[1:])) + " |", "|---|" + "---|" * (len(ed) - 1)]
    for nm in ("MC", "9:MB", "8:J12only", "8:MBTSbased", "8:other"):
        L.append(f"| {nm} | " + " | ".join(f"{x:.3f}" + (f" ± {e:.3f}" if "err" in dd[nm] else "") for x, e in zip(dd[nm]["frac"], dd[nm].get("err", dd[nm]["frac"]))) + " |")
    L.append("")
open(f"{H}/out/tables.md", "w").write("\n".join(L) + "\n"); print(len(L), "wierszy → out/tables.md")
