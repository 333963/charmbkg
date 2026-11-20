#!/usr/bin/env python3
"""XSEC_atlas — tabele wynikowe i budzet (wlasna kopia rachunku; wiersze XSEC przepisane z XSEC_2026-11-19/out/tabele.md, T1; nic tam nie zmieniono)."""
import json, math, os
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = lambda n: json.load(open(os.path.join(H, "out", n)))  # noqa: E731
ds, cz, fo, dx, gm, ng = J("dstar.json"), J("ciecia_D0_Dplus.json"), J("fonll.json")["K"], J("dedx.json"), J("gen_model.json"), J("truth_ngen.json")["bins"]
q = lambda *a: math.sqrt(sum(x * x for x in a))  # noqa: E731
FND = {"3-5": 0.9712, "3.5-5": 0.9712, "5-8": 0.9768, "8-30": 0.9797}
L = ["# XSEC_atlas — tabele wynikowe (427929, L = 0,4648 nb⁻¹; C ekskluzywne; POST HOC względem XSEC)", ""]
# ---------- A
AT = {"3.5-5": (145 * 1.5, 1.5 * q(15, 14)), "5-8": ((43.4 + 20.8) * 1.5, 1.5 * q(4.2, 1.9, 3.6 + 1.7)), "8-30": (6.34 * 4 + 0.757 * 8 + 0.0788 * 10, q(0.50 * 4, 0.101 * 8, 0.0056 * 10, 0.51 * 4 + 0.065 * 8 + 0.0064 * 10))}
STAB = {"3.5-5": 0.0, "5-8": 0.14, "8-30": 0.0}
L += ["## A. σ(D*± → D⁰π± → K∓π±π±), fit Δm", "", "| bin pT(D*) | kandydaci | N_fit | χ²/ndf | skala σ | f_pk | C_ekskl [%] | C(decayType) [%] | domknięcie mb_nd | σ [µb] ± stat ± syst ± L | warianty selekcji σ [µb]: L_xy ≥ 0,10 / pT > 0,6 / okno 55 MeV |", "|---|---|---|---|---|---|---|---|---|---|---|"]
A = {}
for bn in ("3-5", "3.5-5", "5-8", "8-30"):
    r = ds[f"nom|{bn}"]; f = r["fits"]
    fit = max(abs(f[k]["N"] / r["N"] - 1) for k in ("ksztalt_fd", "zakres152", "bin05", "tlo_expo"))
    g = gm.get(f"Dstar|Δm|{bn if bn != '3.5-5' else '3-5'}")["tot"]
    rows = {"stat MC (C)": r["C_err"] / r["C"], "fit (kształt, zakres, bin, tło)": fit, "składowa pikująca (±20 % f_pk)": 0.2 * r["f_pk"], "śledzenie (3 tory × 2,6 %)": 0.077, "rozdzielczość IP / materiał": 0.015,
            "model generatora (pT, η ⊕ krotność 2,3 % lit.)": q(g, 0.023), "nieprompt": 0.017, "stabilność selekcji (L_xy ≥ 0,10; > 2σ)": STAB.get(bn, 0.0), "BR": 0.0106}
    sy = q(*rows.values()); A[bn] = (r, rows, sy)
    v = " / ".join(f"{ds[f'{vn}|{bn}']['sigma_ub']:.0f} ± {ds[f'{vn}|{bn}']['stat_ub']:.0f}" for vn in ("Lxy010", "ptmin06", "okno55"))
    c = r["closure"]
    L.append(f"| {bn} | {r['n_cand']} | {r['N']:.0f} ± {r['N_err']:.0f} | {f['nom']['chi2']:.0f}/{f['nom']['ndf']} | {f['nom']['scale']:.2f} | {r['f_pk']:.3f} | {100 * r['C']:.2f} ± {100 * r['C_err']:.2f} | {100 * r['C_dty']:.2f} | {c['ratio_to_sig_pk']:.2f} ± {c['ratio_err']:.2f} | "
             f"**{r['sigma_ub']:.1f} ± {r['stat_ub']:.1f} ± {r['sigma_ub'] * sy:.1f} ± {r['sigma_ub'] * 0.022:.1f}** | {v} |")
L += ["", "Budżet systematyczny σ(D*) [%]:", "", "| źródło | 3,5–5 | 5–8 | 8–30 |", "|---|---|---|---|"]
for k in A["5-8"][1]:
    L.append(f"| {k} | " + " | ".join(f"{100 * A[b][1][k]:.1f}" for b in ("3.5-5", "5-8", "8-30")) + " |")
L.append("| **razem syst** | " + " | ".join(f"**{100 * A[b][2]:.1f}**" for b in ("3.5-5", "5-8", "8-30")) + " |")
L += ["", "Porównania σ(D*±) [µb]:", "", "| bin | dane | ATLAS 7 TeV (tab. 4, \\|η\\| < 2,1) | K FONLL 7→13,6 | ATLAS przeskalowane | dane / ATLAS skal. | FONLL 13,6 TeV ×2 (c + c̄), f(c→D*) = 0,2429 (min–max) | Pythia A3 (÷ f_ND) | dane / Pythia |", "|---|---|---|---|---|---|---|---|---|"]
sd = sa = 0; vd = va = 0
for bn in ("3.5-5", "5-8", "8-30"):
    r, rows, sy = A[bn]; s, e = AT[bn]; K = fo[bn]; tot = r["sigma_ub"] * q(r["stat_ub"] / r["sigma_ub"], sy, 0.022)
    sc = s * K["K"]; ec = sc * q(e / s, abs(K["K_min"] / K["K"] - 1))
    py = r["pythia_ub"] / FND[bn]
    L.append(f"| {bn} | {r['sigma_ub']:.1f} ± {tot:.1f} | {s:.1f} ± {e:.1f} | {K['K']:.3f} (min/min {K['K_min']:.3f}) | {sc:.1f} ± {ec:.1f} | {r['sigma_ub'] / sc:.2f} ± {r['sigma_ub'] / sc * q(tot / r['sigma_ub'], ec / sc):.2f} | "
             f"{2 * K['sigma13p6'][0]:.0f} ({2 * K['sigma13p6'][1]:.0f}–{2 * K['sigma13p6'][2]:.0f}) | {py:.1f} | {r['sigma_ub'] / py:.2f} ± {tot / py:.2f} |")
# ---------- B
XS = {"D0": {"3-5": (1365, 4105), "5-8": (376, 1530), "8-30": (132, 681)}, "Dplus": {"3-5": (542, 3547), "5-8": (176, 1651), "8-30": (54.8, 687)}}
L += ["", "## B. Selekcja cięciowa (bez klasyfikatora) wobec WP70 drzewa H — te same tabele, ten sam kod fitu (D⁺ M0, D⁰ M8, κ × 1)", "",
      "| gatunek | bin | selekcja | kandydaci | N_fit | p(χ²) | C_ekskl [%] | σ_popr/σ_XSEC = C_XSEC/C_ekskl = N_ekskl/N_decayType | domknięcie mb_nd | σ [µb] ± stat | σ_XSEC × czynnik [µb] (N_XSEC) |", "|---|---|---|---|---|---|---|---|---|---|---|"]
for sp in ("D0", "Dplus"):
    for bn in ("3-5", "5-8", "8-30"):
        for sn in ("WP70", "ciecia"):
            r = cz[f"{sp}|{bn}"][sn]; n = ng[f"{sp}|{bn}"]
            L.append(f"| {sp} | {bn} | {sn} | {r['n_cand']} | {r['N']:.0f} ± {r['N_err']:.0f} | {r['p']:.3f} | {100 * r['C']:.2f} ± {100 * r['C_err']:.2f} | {n['ratio']:.4f} ± {n['ratio_err']:.4f} | {r['closure']:.2f} ± {r['closure_err']:.2f} | "
                     f"{r['sigma_ub']:.1f} ± {r['stat_ub']:.1f} | " + (f"{XS[sp][bn][0] * n['ratio']:.0f} ({XS[sp][bn][1]})" if sn == "WP70" else "—") + " |")
L += ["", "Stosunek R = σ_cięcia/σ_WP70 (SF_klas = 1 po stronie WP70; SF_klas = 1/R przy założeniu, że MC opisuje wydajność cięć prostokątnych):", "",
      "| gatunek | bin | R | ± stat (bootstrap po zdarzeniach, 200 replik) | ρ(N_c, N_k) | ± stat MC | R po poprawce domknięć | warianty R | SF_klas = 1/R ± razem |", "|---|---|---|---|---|---|---|---|---|"]
comb = {}
for sp in ("D0", "Dplus"):
    for bn in ("3-5", "5-8", "8-30"):
        z = cz[f"{sp}|{bn}"]; c = z["ciecia"]; b = c["boot"]
        e = c["R"] * q(b["R_rel_sd"], c["R_mcstat_rel"])
        v = ", ".join(f"{k} {z[k]['R']:.2f}" for k in z if k not in ("WP70", "ciecia"))
        L.append(f"| {sp} | {bn} | **{c['R']:.3f}** | {c['R'] * b['R_rel_sd']:.3f} | {b['rho']:.2f} | {c['R'] * c['R_mcstat_rel']:.3f} | {c['R_closure_corr']:.3f} | {v} | {1 / c['R']:.2f} ± {e / c['R'] ** 2:.2f} |")
        if bn != "3-5":
            comb.setdefault(sp, []).append((c["R"], e))
for sp, v in comb.items():
    w = [1 / e ** 2 for _, e in v]; m = sum(r * x for (r, _), x in zip(v, w)) / sum(w); e = 1 / math.sqrt(sum(w))
    L.append(f"| {sp} | 5–30 (średnia ważona) | **{m:.3f}** | {e:.3f} (stat ⊕ stat MC) | | | | χ² = {sum((r - m) ** 2 * x for (r, _), x in zip(v, w)):.1f}/1 | **{1 / m:.2f} ± {e / m ** 2:.2f}** |")
# ---------- C
L += ["", "## C. Budżet systematyczny — własna kopia (wiersze XSEC z `XSEC_2026-11-19/out/tabele.md` T1 + wiersze nowe); wartości względne [%]", ""]
OLD = {"D0": {"5-8": {"fit: model tła": 1.2, "fit: obciążenie": 8.3, "SF": 12.5, "stat. MC": 2.0, "BR": 0.8, "odbicie": 0.1, "nieprompt": 1.7, "wybór zdarzeń + tryger": 0.2},
              "8-30": {"fit: model tła": 11.3, "fit: obciążenie": 9.4, "SF": 10.1, "stat. MC": 2.9, "BR": 0.8, "odbicie": 0.1, "nieprompt": 0.4, "wybór zdarzeń + tryger": 0.2}},
       "Dplus": {"3-5": {"fit: model tła": 3.9, "fit: obciążenie": 7.8, "SF": 10.0, "stat. MC": 1.6, "BR": 1.5, "odbicie": 0.0, "nieprompt": 1.6, "wybór zdarzeń + tryger": 0.2},
                 "5-8": {"fit: model tła": 2.9, "fit: obciążenie": 12.4, "SF": 10.0, "stat. MC": 2.1, "BR": 1.5, "odbicie": 0.0, "nieprompt": 0.6, "wybór zdarzeń + tryger": 0.2},
                 "8-30": {"fit: model tła": 7.7, "fit: obciążenie": 5.6, "SF": 10.0, "stat. MC": 2.9, "BR": 1.5, "odbicie": 0.0, "nieprompt": 0.6, "wybór zdarzeń + tryger": 0.2}}}
NT = {"D0": 2, "Dplus": 3}
bud = {}
for sp in ("D0", "Dplus"):
    for bn, old in OLD[sp].items():
        z = cz[f"{sp}|{bn}"]
        new = {"śledzenie (2,6 %/tor, 1512.02913; skorelowane)": 2.6 * NT[sp], "rozdzielczość d₀/z₀ + materiał (lit.)": 1.5, "rekonstrukcja wierzchołka (lit., tylko D⁺)": 3.0 if sp == "Dplus" else 0.0,
               "dE/dx (przesunięcie ∓0,023 na MC)": 100 * dx[f"{sp}|{bn}"]["syst"],
               "norm. składowej pikującej (±5,8 %; M8)": (abs(z["WP70"]["N_kap1364"] / z["WP70"]["N"] - 1) / 0.364 * 5.8 if sp == "D0" else 0.0),
               "model generatora (pT, η z MC ⊕ krotność 2,0 % lit.)": q(100 * gm[f"{sp}|WP70|{bn}"]["tot"], 2.0)}
        o, n = q(*old.values()), q(*new.values())
        bud[(sp, bn)] = (old, new, o, n, q(o, n))
keys_o, keys_n = list(OLD["D0"]["5-8"]), None
cols = [("D0", "5-8"), ("D0", "8-30"), ("Dplus", "3-5"), ("Dplus", "5-8"), ("Dplus", "8-30")]
L += ["| źródło | " + " | ".join(f"{sp} {bn}" + (" (scen.)" if sp == "Dplus" else "") for sp, bn in cols) + " |", "|---|" + "---|" * len(cols)]
for k in keys_o:
    L.append(f"| {k} (XSEC) | " + " | ".join(f"{bud[c][0][k]:.1f}" for c in cols) + " |")
L.append("| razem XSEC | " + " | ".join(f"{bud[c][2]:.1f}" for c in cols) + " |")
for k in bud[cols[0]][1]:
    L.append(f"| **nowy:** {k} | " + " | ".join(f"{bud[c][1][k]:.1f}" for c in cols) + " |")
L.append("| razem nowe wiersze | " + " | ".join(f"{bud[c][3]:.1f}" for c in cols) + " |")
L.append("| **razem syst (nowy budżet)** | " + " | ".join(f"**{bud[c][4]:.1f}**" for c in cols) + " |")
alt = []
for c in cols:
    sp, bn = c; z = cz[f"{sp}|{bn}"]; o, n = dict(bud[c][0]), dict(bud[c][1])
    n["śledzenie (2,6 %/tor, 1512.02913; skorelowane)"] = 1.35 * NT[sp]
    alt.append(q(*o.values(), *n.values()))
L.append("| wariant: śledzenie Run 2 (1,35 %/tor, 2302.00336) | " + " | ".join(f"{x:.1f}" for x in alt) + " |")
alt = []
for c in cols:
    sp, bn = c; z = cz[f"{sp}|{bn}"]; o, n = dict(bud[c][0]), dict(bud[c][1])
    if sp == "D0":
        n["norm. składowej pikującej (±5,8 %; M8)"] = abs(z["WP70"]["N_kap1364"] / z["WP70"]["N"] - 1) / 0.364 * 20
    alt.append(q(*o.values(), *n.values()))
L.append("| wariant: składowa pikująca ±20 % | " + " | ".join(f"{x:.1f}" for x in alt) + " |")
alt = []
for c in cols:
    sp, bn = c; o, n = dict(bud[c][0]), dict(bud[c][1]); zc = cz[f"{sp}|{bn}"]["ciecia"]
    o["SF"] = 100 * q(zc["boot"]["R_rel_sd"], zc["R_mcstat_rel"])
    alt.append((o["SF"], q(*o.values(), *n.values())))
L.append("| wariant: SF z ograniczenia cięciowego (błąd 1/R per bin zamiast wiersza SF) | " + " | ".join(f"{t:.1f} (SF {s:.1f})" for s, t in alt) + " |")
L += ["", "Stosunek D⁺/D⁰ (wiersze nowe, część nieskorelowana): śledzenie 2,6 (jeden tor więcej); rozdzielczość/materiał ≈ 1,0; wierzchołek 3 torów 3,0; dE/dx ≤ 0,2; pikująca D⁰ jak wyżej; model generatora ≈ 1,5.", "",
      "| bin | razem bez SF (XSEC T2) | nowe wiersze | razem bez SF (nowy) |", "|---|---|---|---|"]
for bn, o in (("5-8", 18.3), ("8-30", 21.6)):
    n = q(2.6, 1.0, 3.0, 0.2, bud[("D0", bn)][1]["norm. składowej pikującej (±5,8 %; M8)"], 1.5)
    L.append(f"| {bn} | {o:.1f} | {n:.1f} | {q(o, n):.1f} |")
open(os.path.join(H, "tabele/wyniki.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
