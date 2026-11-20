"""EST_spread — tabele markdown do raportu (out/tables.md) + podsumowania domknięcia. Wejście: out/est_new.json, out/est_new_rep.npz, out/collect.json."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
N = json.load(open(f"{OUT}/est_new.json")); C = json.load(open(f"{OUT}/collect.json")); REP = np.load(f"{OUT}/est_new_rep.npz")
S, ROWS = C["spread"], C["rows"]
WPS = N["meta"]["wps"]; BINS = ["3-5", "5-8", "8-30"]; SBE = N["meta"]["sb_est"]
L = []
pc = lambda x, p=1: f"{100*x:.{p}f}".replace(".", ",")
fx = lambda x, p=3: f"{x:.{p}f}".replace(".", ",")

# T1: nowe ε_B, dane, per run
for run, nm in (("d427", "427929"), ("d428", "428770")):
    L += [f"### T1-{nm}. D⁰, dane {nm}: ε_B(SR) × 100 per estymator (NOWE)", "",
          "| bin | WP | exp(pol2) | exp(pol1) | exp(pol3) | zliczenie sym | liniowa L/P | dalekie ≥ 90 | MC V2 (typ SALAD) | MC V1 | MC surowe t≥2 | MC surowe t≥8 |", "|" + "---|" * 12]
    for b in BINS:
        d, m = N["data"][f"{run}|{b}"], N["mc_est"][b]
        for j, wp in enumerate(WPS):
            c = [f"{100*d[e]['val'][j]:.3f} ± {100*d[e]['err'][j]:.3f}".replace(".", ",") for e in ("pol2", "pol1", "pol3", "sym", "lr", "far")]
            c += [f"{100*m[k]['val'][j]:.3f} ± {100*m[k]['err'][j]:.3f}".replace(".", ",") for k in (f"mcV2_t2_{run}", f"mcV1_t2_{run}", "mcraw_t2", "mcraw_t8")]
            L.append(f"| {b} | {int(float(wp)*100)} | " + " | ".join(c) + " |")
    L.append("")

# T2: rozrzut dane D0 (nowe)
def sp(group, pt, wp, st):
    return next(r for r in S if r["group"] == group and r["pt"] == pt and r["wp"] == wp and r["set"].startswith(st))
L += ["### T2. D⁰, dane: rozrzut estymatorów [%] (RMS / (max−min)/2), χ²/ndf zgodności (pełna kowariancja), błąd stat. nominału", "",
      "| run | bin | WP | A: 5 interpolacji | χ²/ndf A | B: pol2, sym, lin., MC V2 | χ²/ndf B | E: pol2, sym, MC V2 | χ²/ndf E | C: B + MC surowe | σ_stat(pol2) | budżet D⁰ |", "|" + "---|" * 12]
BUD0 = {95: "—", 90: "—", 80: "±5", 70: "±10", 50: "±15–20"}
for g, nm in (("G1_D0_dane_427", "427929"), ("G2_D0_dane_428", "428770")):
    for b in BINS:
        for wp in (95, 90, 80, 70, 50):
            a, bb, c, e3 = sp(g, b, wp, "A_"), sp(g, b, wp, "B_"), sp(g, b, wp, "C_"), sp(g, b, wp, "E_")
            L.append(f"| {nm} | {b} | {wp} | {pc(a['rms'])} / {pc(a['half'])} | {a['chi2']:.0f}/{a['ndf']} | {pc(bb['rms'])} / {pc(bb['half'])} | {bb['chi2']:.0f}/{bb['ndf']} | {pc(e3['rms'])} / {pc(e3['half'])} | {e3['chi2']:.1f}/{e3['ndf']} | "
                     f"{pc(c['rms'])} / {pc(c['half'])} | {pc(a['stat_nom'])} | {BUD0[wp]} |")
L.append("")

# pary: istotność różnic (kontrola χ²)
L += ["### T2b. Kontrola: różnice par estymatorów względem exp(pol2), 427929, 3–5 GeV (błąd różnicy z replik)", "", "| WP | sym − pol2 [%] | lin. − pol2 [%] | pol1 − pol2 [%] | MC V2 − pol2 [%] |", "|---|---|---|---|---|"]
Rd, Rm = REP["d427|3-5"], REP["mcest|3-5"]
MCN = ["mcraw_t2", "mcraw_t8", "mcV1_t2_d427", "mcV2_t2_d427", "mcV1_t2_d428", "mcV2_t2_d428"]
for j, wp in enumerate(WPS):
    d = N["data"]["d427|3-5"]; ref = d["pol2"]["val"][j]; c = []
    for e in ("sym", "lr", "pol1"):
        dv = d[e]["val"][j] - ref; de = np.std(Rd[:, SBE.index(e), j] - Rd[:, 0, j])
        c.append(f"{pc(dv/ref)} ± {pc(de/ref, 2)}")
    mv = N["mc_est"]["3-5"]["mcV2_t2_d427"]; dv = mv["val"][j] - ref; de = np.hypot(mv["err"][j], d["pol2"]["err"][j])
    c.append(f"{pc(dv/ref)} ± {pc(de/ref)}")
    L.append(f"| {int(float(wp)*100)} | " + " | ".join(c) + " |")
L.append("")

# T3: domknięcie MC D0 bt
L += ["### T3. D⁰, MC (bt_mb_nd, tło truth ≥ 2): estymator/prawda (NOWE)", "",
      "| bin | WP | exp(pol2) | exp(pol1) | exp(pol3) | zliczenie sym | liniowa L/P | dalekie ≥ 90 | RMS / półzakres 5 interp. [%] | najbliższy prawdy | prawda w [min,max] |", "|" + "---|" * 11]
NAME = {"pol2": "exp(pol2)", "pol1": "exp(pol1)", "pol3": "exp(pol3)", "sym": "zliczenie sym", "lr": "liniowa L/P", "far": "dalekie"}
for b in BINS:
    d = N["mc_closure"][f"t2|{b}"]
    for j, wp in enumerate(WPS):
        r = sp("G3_D0_MC_bt_t2", b, int(float(wp) * 100), "A_")
        c = [f"{fx(d[e]['ratio'][j])} ± {fx(d[e]['ratio_err'][j])}" for e in ("pol2", "pol1", "pol3", "sym", "lr", "far")]
        L.append(f"| {b} | {int(float(wp)*100)} | " + " | ".join(c) + f" | {pc(r['rms'])} / {pc(r['half'])} | {NAME[r['best']]} ({pc(r['best_dev'])} %) | {'tak' if r['truth_in_range'] else 'nie'} |")
L.append("")

# T3b: per estymator — χ² wobec 1 po 15 komórkach, średnie ważone odchylenie, pokrycie
L += ["### T3b. D⁰, MC: jakość estymatorów po 15 komórkach (bin × WP); χ² bez korelacji między WP (przybliżenie)", "",
      "| tło MC | estymator | χ²/15 wobec 1 | średnie ważone odchylenie [%] | średnie \\|odchylenie\\| [%] (3–5 / 5–8 / 8–30) |", "|---|---|---|---|---|"]
summ = {}
for tn in ("t2", "t8"):
    for e in SBE:
        dv, er, per = [], [], []
        for b in BINS:
            d = N["mc_closure"][f"{tn}|{b}"][e]
            v, s_ = np.array(d["ratio"]) - 1, np.array(d["ratio_err"])
            dv += list(v); er += list(s_); per.append(np.abs(v).mean())
        dv, er = np.array(dv), np.array(er); w = 1 / er ** 2
        summ[(tn, e)] = (float((dv ** 2 * w).sum()), float((w * dv).sum() / w.sum()), float(1 / np.sqrt(w.sum())))
        L.append(f"| truth ≥ {tn[1:]} | {NAME[e]} | {(dv**2*w).sum():.1f} | {pc((w*dv).sum()/w.sum(), 2)} ± {pc(1/np.sqrt(w.sum()), 2)} | " + " / ".join(pc(x) for x in per) + " |")
L.append("")
# pokrycie: czy |odchylenie nominału| <= rozrzut ; pull z rozrzutem i statystyką
L += ["### T3c. D⁰, MC (truth ≥ 2): czy rozrzut 5 interpolacji pokrywa odchylenie nominału exp(pol2) od prawdy", "",
      "| bin | WP | odchylenie pol2 [%] | σ_stat [%] | RMS [%] | półzakres [%] | \\|odch.\\| ≤ RMS | \\|odch.\\| ≤ półzakres | \\|odch.\\|/√(RMS²+σ²) | \\|odch.\\|/√(budżet²+σ²) |", "|" + "---|" * 10]
cnt = {"rms": 0, "half": 0, "n": 0, "p_rms": [], "p_bud": []}
BUDV = {95: 0.05, 90: 0.05, 80: 0.05, 70: 0.10, 50: 0.15}
for b in BINS:
    for wp in (95, 90, 80, 70, 50):
        r = sp("G3_D0_MC_bt_t2", b, wp, "A_")
        dv, st = r["nom_dev"], r["nom_err"]
        p1, p2 = abs(dv) / np.hypot(r["rms"], st), abs(dv) / np.hypot(BUDV[wp], st)
        cnt["n"] += 1; cnt["rms"] += abs(dv) <= r["rms"]; cnt["half"] += abs(dv) <= r["half"]; cnt["p_rms"].append(p1); cnt["p_bud"].append(p2)
        L.append(f"| {b} | {wp} | {pc(dv)} | {pc(st)} | {pc(r['rms'])} | {pc(r['half'])} | {'tak' if abs(dv)<=r['rms'] else 'nie'} | {'tak' if abs(dv)<=r['half'] else 'nie'} | {fx(p1,2)} | {fx(p2,2)} |")
L += ["", f"Pokrycie: |odch.| ≤ RMS w {cnt['rms']}/{cnt['n']} komórek, ≤ półzakres w {cnt['half']}/{cnt['n']}. "
      f"Σ pull² z (RMS ⊕ stat) = {sum(x*x for x in cnt['p_rms']):.1f}/15; z (budżet ⊕ stat; WP95/90 liczone z ±5 %) = {sum(x*x for x in cnt['p_bud']):.1f}/15; "
      f"z samą statystyką = {summ[('t2','pol2')][0]:.1f}/15.", ""]

# T4: grupy z literatury projektu — wartości / nominał lub prawda
def tab_group(title, group, ests, ref=None, pts=("all",), wps=(80, 70, 50), names=None, pcnt=False):
    L.extend([title, "", "| bin | WP | " + " | ".join((names or {}).get(e, e) for e in ests) + " |", "|---|---|" + "---|" * len(ests)])
    for pt in pts:
        for wp in wps:
            rr = {r["estimator"]: r for r in ROWS if r["group"] == group and r["pt"] == pt and r["wp"] == wp}
            if not rr:
                continue
            rv = rr[ref]["value"] if ref else 1.0
            c = []
            for e in ests:
                if e not in rr:
                    c.append("—"); continue
                v, er = rr[e]["value"] / rv, rr[e]["err"]
                c.append(fx(v) + (f" ± {fx(er / rv)}" if er != "" else ""))
            L.append(f"| {pt} | {wp} | " + " | ".join(c) + " |")
    L.append("")

tab_group("### T4a. D⁺, dane (podpróbka S8): ε_B(SR) estymatora / ε_B z fitu pol2 (błąd podany tylko dla fitu)", "G4_Dp_dane_S8sub",
          ["pol2_fit", "sb_sym", "salad_hgb", "salad_sk_mlp", "mcraw_t8", "sb_asym", "sb_left", "sb_right"], "pol2_fit", ("3-5", "5-8", "8-30", "all"),
          names={"pol2_fit": "fit pol2", "sb_sym": "SB sym", "salad_hgb": "SALAD hgb", "salad_sk_mlp": "SALAD sk_mlp", "mcraw_t8": "MC surowe", "sb_asym": "SB całe (asym.)", "sb_left": "SB lewe", "sb_right": "SB prawe"})
tab_group("### T4b. D⁰, dane (podpróbka B1, 3–30 GeV): estymator / fit pol2", "G7_D0_dane_B1",
          ["pol2_fit", "sb_copy", "salad_hgb_t2", "salad_sk_mlp_t2", "salad_hgb_t8", "salad_sk_mlp_t8", "mcraw_t2", "mcraw_t8"], "pol2_fit",
          names={"pol2_fit": "fit pol2", "sb_copy": "SB zliczenie", "salad_hgb_t2": "SALAD hgb t≥2", "salad_sk_mlp_t2": "SALAD sk_mlp t≥2", "salad_hgb_t8": "SALAD hgb t≥8", "salad_sk_mlp_t8": "SALAD sk_mlp t≥8", "mcraw_t2": "MC surowe t≥2", "mcraw_t8": "MC surowe t≥8"})
tab_group("### T4c. D⁰, dane (S4, bieg lxplus full, 3–30 GeV): estymator / fit", "G8_D0_dane_S4", ["pol2_fit", "sb_copy", "salad_mlp", "feta"], "pol2_fit",
          names={"pol2_fit": "fit Gauss+exp(pol2)", "sb_copy": "SB zliczenie", "salad_mlp": "SALAD mlp", "feta": "FETA"})
tab_group("### T4d. D⁰, MC (S4, bieg lxplus full, 3–30 GeV): ε_B przewidziane / prawda", "G9_D0_MC_S4",
          ["sb_copy", "salad_mlp", "salad_hgb", "salad_hgb_nom", "feta", "flow", "mcraw_t8", "ideal"],
          names={"sb_copy": "SB zliczenie", "salad_mlp": "SALAD mlp", "salad_hgb": "SALAD hgb", "salad_hgb_nom": "SALAD hgb bez m", "feta": "FETA", "flow": "flow p(x\\|m)", "mcraw_t8": "MC surowe", "ideal": "druga kopia prawdy"})
for g, t in (("G10_D0_MC_S6", "D⁰"), ("G11_Dp_MC_S6", "D⁺")):
    tab_group(f"### T4e-{t}. {t}, MC (S6, 3–30 GeV): N_B(SR) z szablonu masy unormowanego w pasmach bocznych / prawda", g,
              ["pol2", "rot_rand", "rot_pi", "flip_eta", "mix", "swap"],
              names={"pol2": "fit pol2", "rot_rand": "rotacja losowa", "rot_pi": "rotacja π", "flip_eta": "odbicie η", "mix": "tor z innego zdarzenia", "swap": "zamiana K↔π"})
tab_group("### T4f. D⁺, MC (VAL smoke_vm): estymator / prawda", "G12_Dp_MC_VAL", ["pol2", "sb_sym", "sb_asym", "inner"], None, ("3-5", "5-8", "8-30", "all"),
          names={"pol2": "fit pol2", "sb_sym": "SB sym", "sb_asym": "SB całe (asym.)", "inner": "pasmo wewn."})
tab_group("### T4g. D⁺, MC (S8, gładkie zniekształcenie): ε_B przewidziane / prawda", "G13_Dp_MC_S8smooth", ["salad_hgb", "salad_sk_mlp", "sb_copy", "mcraw_t8"], None,
          ("3-5", "5-8", "8-30", "all"), names={"salad_hgb": "SALAD hgb", "salad_sk_mlp": "SALAD sk_mlp", "sb_copy": "SB całe (asym.)", "mcraw_t8": "MC surowe"})
tab_group("### T4h. D⁰, MC (S, 3–30 GeV, bez weta): ε_B kontroli / prawda", "G15_D0_MC_S2gen", ["sb_both", "sb_near", "copula", "flow"],
          names={"sb_both": "SB oba", "sb_near": "SB bliskie", "copula": "kopuła", "flow": "flow"})
tab_group("### T4i. D⁰, MC (S, z wetem refleksji): ε_B kontroli / prawda", "G15b_D0_MC_S2veto", ["sb_both_veto", "sb_far_veto", "sb_near_veto", "copula_veto"],
          names={"sb_both_veto": "SB oba + weto", "sb_far_veto": "SB dalekie + weto", "sb_near_veto": "SB bliskie + weto", "copula_veto": "kopuła + weto"})
tab_group("### T4j. D⁰, MC (S2 cloud_v1): ε_B generatora w SR / prawda", "G16_D0_MC_S2flow", ["flow", "copula"])
for g, t in (("G5_Dp_dane_S8full", "D⁺"), ("G6_D0_dane_S8full", "D⁰")):
    tab_group(f"### T4k-{t}. {t}, dane 428770 pełna lumi (S8 lxplus): ε_B × 1 (wartości bezwzględne; brak fitu pol2 w źródle)", g, ["salad_sk_mlp", "mcraw_t8", "sb_asym"], None,
              ("3-5", "5-8", "8-30", "all"), names={"salad_sk_mlp": "SALAD sk_mlp", "mcraw_t8": "MC surowe", "sb_asym": "SB całe (asym.)"})

# T5: rozrzuty pozostałych grup
L += ["### T5. Rozrzut w pozostałych grupach (RMS / półzakres [%]; na MC: odchylenie średniej i najlepszego od prawdy)", "",
      "| grupa | bin | WP | zestaw | n | RMS | półzakres | średnia − prawda | najlepszy (odch.) | prawda w [min,max] |", "|" + "---|" * 10]
for r in S:
    if r["group"].startswith(("G1_", "G2_", "G3_")):
        continue
    if r["pt"] != "all" and not r["group"].startswith(("G4_", "G12_")):
        continue
    mc = r["mode"] != "dane"
    L.append(f"| {r['group']} | {r['pt']} | {r['wp']} | {r['set']} | {r['n']} | {pc(r['rms'])} | {pc(r['half'])} | " +
             (f"{pc(r['mean_dev'])} | {r['best']} ({pc(r['best_dev'])}) | {'tak' if r['truth_in_range'] else 'nie'} |" if mc else "— | — | — |"))
L.append("")
open(f"{OUT}/tables.md", "w").write("\n".join(L))
print("\n".join(L[-90:]))
