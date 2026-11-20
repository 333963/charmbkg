#!/usr/bin/env python3
"""AFPTRIG/report_md.py — RAPORT.md z facc_afptrig.json i plików pomocniczych (liczby wstawiane automatycznie)."""
import json, math, time
import lib
R = json.load(open("facc_afptrig.json")); DESC = json.load(open("describe.json")); EX = json.load(open("extra.json")); LK = json.load(open("leak.json")); SV = json.load(open("strat_var.json"))
RAW = {r: json.load(open(f"results_raw_{r}.json"))[r] for r in ("427929", "428770")}
KR = open("KRYTERIA_WALIDACJI.txt").read()
RUNS = ("427929", "428770")
def n(x, k=0): return "—" if x is None else f"{x:,.{k}f}".replace(",", " ")
def f4(x): return "—" if x is None else f"{x:.4f}"
def pc(x, k=1): return "—" if x is None else f"{100 * x:.{k}f} %"
SN = {"afp_any": "AFP (dowolny łańcuch)", "pt4_L1AFP": "sptrk_pt4_L1AFP", "noalg_L1AFP": "noalg_L1AFP", "nonAFP": "nie-AFP"}
PN = {"twin": "bliźniak sptrk_pt4_L1MBTS_2", "mb": "MB (noalg_mb/sptrk L1MBTS)", "mbts": "bliźniak ∪ MB", "nonAFP": "wszystkie nie-AFP", "j12": "L1J12", "self": "sam strumień"}
L = []
A = L.append
A("# Tło przypadkowego tagu AFP w strumieniu trygera AFP — metoda kształtu Δ = ξ_cal − ξ_AFP (WSTĘPNE)")
A("")
A(f"Wygenerowano {R['meta']['generated']} przez `AFPTRIG/report_md.py`. Projekt B (ATLAS minbias 2022, runy 427929 i 428770, h5 v2.0, kandydaci D po luźnej selekcji, bez score).")
A("Praca lokalna, tylko odczyt danych z `AFP_mix/runs/lxplus/`; kod metody v2 (`charmbkg`) i wyniki `out/`, `out_v2/`, `TEMPLATE/` nietknięte.")
A("**Czego tu nie ma:** uzysków sygnału D w zdarzeniach z tagiem, fitów masy, przekrojów. N_tag to liczba zdarzeń albo kandydatów (sygnał + kombinatoryka) w oknie.")
A("")
A("## 0. Wynik w skrócie")
A("")
A("| run | strumień | N_tag(SR) zdarzeń z D (A + C) | N_acc | **f_acc ± stat + syst − syst** | − syst z efektem lustrzanym | f_acc przy normalizacji w Δ > 0,10 | K1 | K4 | K5 | używalna (K6) |")
A("|---|---|---|---|---|---|---|---|---|---|---|")
for run in RUNS:
    for st in ("afp_any", "pt4_L1AFP", "noalg_L1AFP", "nonAFP"):
        if run == "427929" and st == "noalg_L1AFP":
            continue
        S = R["runs"][run]["streams"][st]; K = S["validation"]
        A(f"| {run} | {SN[st]} | {n(S['n_tag_SR'])} | {n(S['n_acc'], 1)} ± {n(S['n_acc_err_stat'], 1)} | **{f4(S['f_acc'])} ± {f4(S['f_acc_err_stat'])} + {f4(S['f_acc_err_syst_up'])} − {f4(S['f_acc_err_syst_down_no_mirror'])}** | − {f4(S['f_acc_err_syst_down'])} | {f4(S['f_acc_far_norm'])} ± {f4(S['f_acc_far_norm_err_stat'])} | {K['K1']} | {K['K4']['status']} | {K['K5']['status']} | {'TAK' if S['K6_usable'] else 'NIE'} |")
A("")
S8 = R["runs"]["428770"]["streams"]["afp_any"]; S7 = R["runs"]["427929"]["streams"]["afp_any"]; b8 = S8["budget_rel"]
A(f"- **428770, strumień AFP:** f_acc = {f4(S8['f_acc'])} ± {f4(S8['f_acc_err_stat'])} (stat) + {f4(S8['f_acc_err_syst_up'])} − {f4(S8['f_acc_err_syst_down_no_mirror'])} (syst); z nierozstrzygniętym efektem lustrzanym dół rośnie do − {f4(S8['f_acc_err_syst_down'])}. "
  f"Wartość jest bliska f_acc strumienia nie-AFP (v2: 0,0477 ± 0,0135). **Walidacja K1 (domknięcie w obszarze kontrolnym) i K5 (stabilność względem marginesu) nie przechodzą**: "
  f"dane mają nad szablonem nadwyżkę {pc(S8['validation']['norm>0.10_pred_0.045-0.10']['ratio'] - 1, 0)} w 0,045 < Δ < 0,10, układającą się wzdłuż ξ_cal ≈ ξ_AFP + stała (rys. 2) i rosnącą z μ. "
  f"To przeciek własnych protonów zdarzenia ponad cięcie, nie tło przypadkowe; zawyża normalizację. Normalizacja dalej od cięcia (Δ > 0,10, wybór po obejrzeniu walidacji) daje f_acc = {f4(S8['f_acc_far_norm'])} ± {f4(S8['f_acc_far_norm_err_stat'])}; różnica jest w systematyce „w dół”.")
A(f"- **427929, strumień AFP** (99,98 % zdarzeń ma też tryger MBTS): f_acc = {f4(S7['f_acc'])} ± {f4(S7['f_acc_err_stat'])} + {f4(S7['f_acc_err_syst_up'])} − {f4(S7['f_acc_err_syst_down_no_mirror'])}; wszystkie kryteria przechodzą, błąd statystyczny {pc(S7['budget_rel']['stat'], 0)} (413 zdarzeń w obszarze kontrolnym).")
A("- **Łańcuch HLT_noalg_L1AFP_A_OR_C w 428770:** metoda nie przechodzi walidacji (K1: 1,53; kształt K4 FAIL) — pula centralna z łańcuchów L1MBTS nie opisuje tych zdarzeń. Liczba w tabeli jest tylko orientacyjna.")
A("- Metoda kształtu w strumieniu nie-AFP odtwarza v2 (0,0087 i 0,0495 wobec 0,0087 i 0,0477), ale to w dużej części ta sama informacja (§4.2) — sprawdza implementację, nie fizykę.")
A("- Ograniczenia, których metoda nie usuwa: obszar ξ_cal ≲ 0,075 po stronie protonu pozostaje nietestowalny (toy: skok częstości ×2 poniżej 0,07 daje bias −35 % bez żadnego sygnału w walidacji), efekt lustrzany (−19 %) jest nierozstrzygnięty.")
A("")
A("## 1. Kryteria walidacji (zapisane przed obejrzeniem rozkładów Δ; plik `KRYTERIA_WALIDACJI.txt`)")
A("")
A("```"); A(KR.strip()); A("```")
A("")
A("Odstępstwa od zapisu, wszystkie po obejrzeniu wyniku i tak oznaczone w tekście: (1) wariant „normalizacja w Δ > 0,10” (margines 0,08) jako liczba pomocnicza i źródło systematyki „przeciek”; "
  "(2) błąd różnicy wariantów marginesu liczony wzorem dla obszarów zagnieżdżonych (`make_report.py::sd_nested`), bo pierwotny |σ_a² − σ_b²| dawał zero przy równych błędach; "
  "(3) w K4 zakres 0,045–0,30 z łączeniem przedziałów do ≥ 10 oczekiwanych, szablon przeskalowany do danych w tym zakresie (test samego kształtu).")
A("")
A("## 2. Opis próbki")
A("")
A("Selekcja wspólna: GRL ∧ nPileupVertices == 0 ∧ ≥ 1 zapisany łańcuch. „Z D” = ≥ 1 kandydat D⁺, D⁰ lub D* (slim). Tag po stronie = dokładnie jeden proton po stronie, 0,03 ≤ ξ ≤ 0,10. "
  "SR: Δ ≤ 0,02 (cięcie analizy); przerwa: 0,02 < Δ ≤ 0,045; CR: Δ > 0,045 („proton niemożliwy”).")
A("")
A("| run | strumień / łańcuch | zdarzeń | z D | kandydatów w oknach | ⟨nSelBL⟩ z D | strona | tagów z D | SR | przerwa | CR | w tym ξ_cal > 0,14 | CR/SR |")
A("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for run in RUNS:
    for nm, lab in (("noalg_L1AFP", "HLT_noalg_L1AFP_A_OR_C"), ("sptrk_pt4_L1AFP", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"), ("AFP_any", "AFP: dowolny z obu"), ("AFP_only", "AFP i żaden nie-AFP"), ("nonAFP", "nie-AFP (dowolny)"),
                    ("sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2 (bliźniak)"), ("MB_noalg_L1MBTS", "HLT_noalg_mb_L1MBTS_1/2"), ("L1J12", "HLT_noalg_L1J12")):
        r = DESC[run]["streams"][nm]
        for s in "AC":
            sr, cr = r[f"hasD_SR_{s}"], r[f"hasD_CR_{s}"]
            A(f"| {run} | {lab} | {n(r['n_ev'])} | {n(r['n_ev_hasD'])} | {n(r['n_cand_windows'])} | {n(r['mean_nsel_hasD'], 1)} | {s} | {n(r[f'hasD_tag_{s}'])} | {n(sr)} | {n(r[f'hasD_gap_{s}'])} | {n(cr)} | {n(r[f'hasD_CR_xical>0.14_{s}'])} | {cr / sr:.3f} |" if sr else f"| {run} | {lab} | {n(r['n_ev'])} | {n(r['n_ev_hasD'])} | — | — | {s} | 0 | 0 | 0 | 0 | 0 | — |")
A("")
p7, p8 = DESC["427929"]["prescale"], DESC["428770"]["prescale"]
i8 = DESC["428770"]["info"]; i7 = DESC["427929"]["info"]
A("Preskale (nie ma ich w ntupli; oszacowanie z danych: wśród zdarzeń zapisanych przez łańcuch nie-AFP z bitem L1_AFP_A_OR_C „before prescale”, a dla pt4 dodatkowo z odpalonym bliźniakiem pt4_L1MBTS_2, odwrotność ułamka z odpalonym łańcuchem AFP):")
A("")
A("| run | łańcuch | zdarzeń odniesienia | efektywna preskala | w 8 blokach LB |")
A("|---|---|---|---|---|")
for run, p in (("427929", p7), ("428770", p8)):
    for c in ("noalg_L1AFP", "sptrk_pt4_L1AFP"):
        A(f"| {run} | {c} | {n(p[c]['n_den'])} | {p[c]['prescale_eff']:.2f} | {', '.join(f'{x:.1f}' for x in p[c]['per_lb_block'])} |")
A("")
A(f"- 428770: zdarzenia z D w strumieniu AFP pochodzą w 87 % z łańcucha sptrk_pt4 (preskala ≈ 1, na początku runu 1,8); noalg ma preskalę ≈ 23. "
  f"{n(i8['n_none'])} zdarzeń po GRL nie ma żadnego z zapisanych bitów trygera (inne łańcuchy; 41 % z tagiem) — poza zakresem.")
A(f"- 427929: łańcuchy AFP bez preskali, ale MBTS też — tylko {i7['n_afp_only_grl']} zdarzeń wyłącznie z AFP. Strumień AFP jest tu podzbiorem strumienia MBTS, więc służy do walidacji.")
ea, ec = p8["L1AFP_TBP_eff_vs_xi_A"], p8["L1AFP_TBP_eff_vs_xi_C"]; ea7, ec7 = p7["L1AFP_TBP_eff_vs_xi_A"], p7["L1AFP_TBP_eff_vs_xi_C"]
A(f"- Wydajność L1_AFP_A_OR_C (bit TBP) dla zdarzeń nie-AFP z dokładnie jednym protonem po jednej stronie, w przedziałach ξ 0,03–0,04 … 0,09–0,10: strona A {', '.join(f'{x:.2f}' for x in ea7['eff'][2:9])} (427929), {', '.join(f'{x:.2f}' for x in ea['eff'][2:9])} (428770); "
  f"strona C {', '.join(f'{x:.2f}' for x in ec7['eff'][2:9])} (427929), {', '.join(f'{x:.2f}' for x in ec['eff'][2:9])} (428770). Strona C ma wydajność o połowę mniejszą — stąd asymetria liczby tagów A/C w strumieniu AFP.")
A("")
A("Rozkłady ξ_AFP (tagi z D, strumień AFP, udział w przedziałach 0,03–0,04 … 0,09–0,10):")
A("")
A("| run | strona | próbka | N | udziały | ⟨ξ⟩ |")
A("|---|---|---|---|---|---|")
for run in RUNS:
    for s in "AC":
        x = EX[run]["afp_any"][s]["xi_spectra"]
        for k, lab in (("SR_stream_hasD", "SR (Δ ≤ 0,02)"), ("CR_stream_hasD", "CR (Δ > 0,045)"), ("impossible_stream", "ξ_cal > 0,14, strumień AFP (pula protonów)"), ("impossible_nonAFP", "ξ_cal > 0,14, strumień nie-AFP")):
            h = x[k]["h"]; t = max(sum(h), 1)
            A(f"| {run} | {s} | {lab} | {n(x[k]['n'])} | {', '.join(f'{v / t:.2f}' for v in h)} | {x[k]['mean']:.4f} ± {x[k]['mean_err']:.4f} |")
A("")
A("Protony z SR mają widmo rosnące z ξ, protony „niemożliwe” — malejące. Pula protonów do szablonu to te drugie (ξ_cal > 0,14: nie zależy od ξ, więc widmo nieobciążone cięciem). Rozkład Δ: rys. 1 i histogramy w `results_raw_*.json` (`hist_delta`), ξ_cal: `describe.json`.")
A("")
A("![Rozkład Δ: dane i szablon](fig1_delta.png)")
A("")
A("Rys. 1. Rozkład Δ zdarzeń z D i tagiem w strumieniu AFP (strony A + C) i szablon tła przypadkowego znormalizowany w Δ > 0,045. Na dole stosunek dane/szablon dla Δ > 0,02. W 428770 nadwyżka danych do Δ ≈ 0,07.")
A("")
A("## 3. Metoda")
A("")
A("Zdarzenie z tagiem przypadkowym = dowolne oddziaływanie z D bez własnego protonu po stronie s ⊕ niezależny proton. Tryger AFP zapisuje je, bo odpala od protonu przypadkowego, więc liczby „wystawionych” zdarzeń centralnych nie znamy i częstości nie da się użyć. Znamy natomiast kształt:")
A("")
A("    N_acc(SR) = Σ_c N_CR,c(dane) · R_c,     R_c = T_c(Δ ≤ 0,02) / T_c(Δ > 0,045)")
A("")
A("- T_c: wszystkie pary (zdarzenie centralne z puli, proton z puli protonów) w komórce c; Δ = ξ_cal,s zdarzenia centralnego − ξ protonu. Liczone dokładnie (bez losowania).")
A("- **Pula protonów:** tagi strumienia AFP po stronie s ze zdarzeń z ξ_cal,s > 0,14 (z D i bez D), ten sam blok LB (8 bloków o równej liczbie tagów).")
A("- **Pula centralna (nominał):** zdarzenia z D z łańcuchów nie-AFP o tym samym warunku HLT, bez protonu po stronie s: dla sptrk_pt4_L1AFP bliźniak HLT_mb_sptrk_pt4_L1MBTS_2; dla noalg_L1AFP łańcuchy MB (noalg_mb_L1MBTS, sptrk(_pt2)_L1MBTS); dla sumy łańcuchów AFP bliźniak ∪ MB; w 427929 wszystkie nie-AFP (tam to MBTS bez preskali).")
A("- **Komórki:** blok LB × klasa nSelBL (0, 5, 10, 15, 20, 25, 30, 40, 50, 80, ∞). Normalizacja w CR osobno w każdej komórce jest równoważna przeważeniu puli do rozkładu aktywności zdarzeń przypadkowych widzianego w CR, więc usuwa zarówno różnicę aktywności pula ↔ strumień, jak i zależność prawdopodobieństwa protonu przypadkowego od nSelBL.")
A("- Błąd stat.: Poisson N_CR ⊕ skończona pula centralna ⊕ skończona liczba protonów (estymator ilorazowy). Kategorie D: wagi = liczba kandydatów zdarzenia w kategorii, po stronie danych i puli.")
A("- f_acc(A + C) = (N_acc,A + N_acc,C)/(N_SR,A + N_SR,C). Zdarzenia z tagiem po obu stronach (428770, z D: 1 032 z 99 959 tagów) wchodzą po obu stronach; regionu „A∨C z cięciem na każdej stronie” nie odtwarzano.")
A("")
A("### 3.1 Warianty puli centralnej (N_acc, strony A + C; nominał pogrubiony)")
A("")
A("| run | strumień | " + " | ".join(PN[p] for p in ("twin", "mb", "mbts", "nonAFP", "j12", "self")) + " | nominał bez warstw nSelBL |")
A("|---|---|---|---|---|---|---|---|---|")
for run in RUNS:
    for st in ("afp_any", "pt4_L1AFP", "noalg_L1AFP", "nonAFP"):
        if run == "427929" and st == "noalg_L1AFP":
            continue
        o = RAW[run][st]; cells = []
        for p in ("twin", "mb", "mbts", "nonAFP", "j12", "self"):
            if p in o["A"]["pool"]:
                v = o["A"]["pool"][p]["n_acc"] + o["C"]["pool"][p]["n_acc"]; e = math.hypot(o["A"]["pool"][p]["err"], o["C"]["pool"][p]["err"])
                t = f"{n(v, 0)} ± {n(e, 0)}"; cells.append(f"**{t}**" if p == o["nominal_pool"] else t)
            else:
                cells.append("—")
        v = o["A"]["strat"]["no_nsel"]["n_acc"] + o["C"]["strat"]["no_nsel"]["n_acc"]
        A(f"| {run} | {SN[st]} | " + " | ".join(cells) + f" | {n(v, 0)} |")
A("")
A("- **Tasowanie wewnątrz strumienia AFP (pula „sam strumień”) nie nadaje się**: strumień to w ~95 % prawdziwe tagi z małym ξ_cal, więc szablon przewiduje kilkanaście razy za dużo zdarzeń w SR (toy: ×8–20). Wzorzec CMS-TOTEM działa, bo tam próbka centralna (dijety) nie jest wybrana protonem.")
A("- Pula nie-AFP zdominowana przez L1J12 daje po warstwowaniu +7…8 % względem bliźniaka; bez warstw nSelBL wynik spada o 19–23 % — warstwowanie jest konieczne (toy: pula obciążona aktywnością bez warstw daje 0,43 prawdy, z warstwami 1,00).")
A("")
A("## 4. Walidacja")
A("")
A("### 4.1 K1, K4, K5: obszar kontrolny (zdarzenia z D, A + C)")
A("")
A("| run | strumień | norm. Δ > 0,10 → 0,045–0,10: obs / pred | stosunek (z) | 0,045–0,07 | 0,07–0,10 | norm. Δ > 0,20 → 0,10–0,20 | K1 | K4: χ²/ndf (p), maks. odchylenie | K4 | margines 0,01 / 0,04 / 0,06 / 0,08 / 0,12: N_acc/nominał | K5 |")
A("|---|---|---|---|---|---|---|---|---|---|---|---|")
for run in RUNS:
    for st in ("afp_any", "pt4_L1AFP", "noalg_L1AFP", "nonAFP"):
        if run == "427929" and st == "noalg_L1AFP":
            continue
        K = R["runs"][run]["streams"][st]["validation"]; a = K["norm>0.10_pred_0.045-0.10"]; k4 = K["K4"]; k5 = K["K5"]["rows"]
        def rr(k): return f"{K[k]['ratio']:.2f} ± {K[k]['ratio_err']:.2f}"
        A(f"| {run} | {SN[st]} | {n(a['obs'])} / {n(a['pred'], 0)} ± {n(a['pred_err'], 0)} | {a['ratio']:.3f} ± {a['ratio_err']:.3f} ({a['z']:+.1f}) | {rr('norm>0.10_pred_0.045-0.07')} | {rr('norm>0.10_pred_0.07-0.10')} | {rr('norm>0.20_pred_0.10-0.20')} | **{K['K1']}** | "
          f"{k4['chi2']:.0f}/{k4['ndf']} ({k4['p']:.2g}), {pc(k4['max_grouped_dev'])} | **{k4['status']}** | " + " / ".join(f"{k5[m]['ratio']:.2f}" for m in ("0.01", "0.04", "0.06", "0.08", "0.12")) + f" | **{K['K5']['status']}** |")
A("")
A("- 427929: wszystko przechodzi w granicach (dużych) błędów. 428770, strumień nie-AFP: przechodzi.")
A("- **428770, strumień AFP (suma i sptrk_pt4): K1 i K5 FAIL.** Dane leżą nad szablonem o 12–16 % w 0,045 < Δ < 0,10 i o 5–8 % w 0,10–0,20; N_acc maleje monotonicznie z marginesem (−13 % przy normalizacji w Δ > 0,10). K4 formalnie PASS tylko dzięki klauzuli „odchylenie zgrupowane < 10 %” (p ≪ 0,01).")
A("- 428770, noalg_L1AFP: K1 1,53, K4 FAIL — szablon z puli MB nie opisuje danych także daleko od cięcia (norm. > 0,20 → 0,10–0,20: 1,78). Przyczyna nieustalona (pula L1MBTS a zdarzenia bez wymogu MBTS; możliwy większy udział przecieku przy SR/CR = 6,8). Dla tego łańcucha metoda w obecnej postaci nie jest zwalidowana.")
A("")
A("![Diagnostyka i skan marginesu](fig2_walidacja.png)")
A("")
A("Rys. 2. Lewy i środkowy: dane/szablon w funkcji ξ_cal dla czterech przedziałów ξ_AFP (szablon znormalizowany w ξ_cal > 0,20). Niedopasowanie puli dałoby jedną krzywą dla wszystkich ξ_AFP; nadwyżka zaczyna się zaś tuż nad cięciem każdego przedziału ξ_AFP (ξ_cal ≈ ξ_AFP + 0,02…0,07), czyli jest skorelowana z protonem. Prawy: N_acc względem nominału w funkcji dolnej granicy obszaru normalizacji.")
A("")
A("### 4.2 K2: zgodność z metodą częstościową v2 w strumieniu nie-AFP")
A("")
A("| run | strona | N_acc kształt | N_acc v2 | f_acc kształt (A + C) | f_acc v2 (A∨C) | z | K2 |")
A("|---|---|---|---|---|---|---|---|")
for run in RUNS:
    k = R["runs"][run]["K2"]
    for s in "AC":
        A(f"| {run} | {s} | {n(k['per_side'][s]['shape'], 1)} ± {n(k['per_side'][s]['shape_err'], 1)} | {n(k['per_side'][s]['v2'], 1)} ± {n(k['per_side'][s]['v2_err'], 1)} | | | | |")
    A(f"| {run} | A + C | {n(k['n_acc_shape'], 1)} | {n(k['n_acc_v2'], 1)} | {f4(k['f_shape'])} ± {f4(k['f_shape_err'])} | {f4(k['f_v2'])} ± {f4(k['f_v2_err'])} | {k['z']:+.2f} | **{k['status']}** |")
A("")
A("**Zastrzeżenie:** zgodność jest w dużej części z konstrukcji. v2 liczy N_acc = (otagowani dawcy z ξ_cal > 0,14)/(dawcy) × Σ_biorcy P(para w SR); metoda kształtu z CR ograniczonym do ξ_cal > 0,14 byłaby algebraicznie tym samym. "
  "Różnią się tylko dodatkiem obszaru Δ > 0,045 ∧ ξ_cal ≤ 0,14 (ok. 20 % CR) i komórkami. K2 sprawdza więc implementację i to, że normalizacja w CR nie psuje wyniku tam, gdzie v2 działa — nie jest niezależnym pomiarem. "
  f"Dodatkowy test z 427929: w podpróbce z odpalonym łańcuchem AFP (wydajność L1 ≈ 0,9 / 0,45) f_acc = {f4(S7['f_acc'])} ± {f4(S7['f_acc_err_stat'])} wobec 0,0087 ± 0,0044 (v2) — przerzedzenie trygerem protonowym nie zmienia frakcji.")
A("")
A("### 4.3 K3: toy z trygerem wymagającym protonu")
A("")
T = R["K3_toy"]
A(f"Generator `charmbkg.afp_mix_toy.generate` bez zmian (3·10⁶ zdarzeń na seed, {T['n_seeds']} seedów) + dodany w `toy.py` tryger: L1 odpala od każdego widzianego protonu z ε_A = 0,9, ε_C = 0,45; strumień AFP = zdarzenia z odpalonym L1; pula centralna = wszystkie zdarzenia (jak MBTS w 427929) albo klasa trygera 1 (obciążona aktywnością, ≈ J12). Prawda: tag w SR, który nie jest własnym protonem zdarzenia. Stosunki N_acc(est)/N_acc(prawda), A + C, średnia ± błąd średniej:")
A("")
A("| scenariusz | f_acc prawda | **nominał** | bez warstw nSelBL | pula „J12” z warstwami | pula „J12” bez warstw | tasowanie w strumieniu | protony: wszystkie tagi | margines 0,01 | margines 0,06 | domknięcie K1 (obs/pred) |")
A("|---|---|---|---|---|---|---|---|---|---|---|")
for sc, r in T["table"].items():
    def c(k): return f"{r[k][0]:.2f} ± {r[k][1]:.2f}"
    A(f"| {sc} | {r['f_acc_true']:.3f} | **{c('nominal')}** | {c('no_nsel')} | {c('pool_J12_strat')} | {c('pool_J12_nostrat')} | {r['pool_self'][0]:.1f} | {c('prot_alltags')} | {c('m0.01')} | {c('m0.06')} | {c('closure')} |")
A("")
A(f"- **K3 {T['status']}**: nominał {T['nominal_bias']:.3f} ± {T['nominal_bias_err']:.3f}. Zależność częstości od nSelBL, większa częstość i wydajność trygera zależna od ξ nie psują estymatora (1,01–1,04).")
A("- **Bias, gdy widmo ξ protonów w puli różni się od widma w obszarze użycia** (`halo_xi_tilt_low`: inne widmo halo przy ξ_cal < 0,14): −19 %; domknięcie K1 to wykrywa (1,23). Pula „wszystkie tagi strumienia” (skażona protonami sygnału) daje −9 %.")
A("- **Martwe pole jak w v2:** częstość przypadkowych ×2 tylko przy ξ_cal < 0,07 → estymator 0,65 prawdy, domknięcie 0,97 (nic nie widać). Skok poniżej 0,12 (×1,5): 0,72, domknięcie 1,11.")
A("- **Przeciek własnych protonów** (`xical_smear_0.5`): +32 %, domknięcie 1,59, N_acc maleje z marginesem (1,51 → 1,15) — ten sam obraz co w danych 428770, w danych o mniejszej amplitudzie.")
A("")
A("![Toy](fig3_toy.png)")
A("")
A("Rys. 3. Toy: N_acc(estymator)/N_acc(prawda). Szary pas: ±5 % (kryterium K3).")
A("")
A("## 5. Wynik: f_acc w strumieniu AFP")
A("")
A("### 5.1 Per strona (zdarzenia z D)")
A("")
A("| run | strumień | strona | N_tag(SR) | N_CR | N_acc ± stat | f_acc ± stat | pula centralna (zdarzeń) | protonów w puli |")
A("|---|---|---|---|---|---|---|---|---|")
for run in RUNS:
    for st in ("afp_any", "pt4_L1AFP", "noalg_L1AFP"):
        if run == "427929" and st == "noalg_L1AFP":
            continue
        S = R["runs"][run]["streams"][st]
        for s in "AC":
            x = S["sides"][s]
            A(f"| {run} | {SN[st]} | {s} | {n(x['n_sr'])} | {n(x['n_cr'])} | {n(x['n_acc'], 1)} ± {n(x['err'], 1)} | {f4(x['f_acc'])} ± {f4(x['f_acc_err_stat'])} | {PN[S['central_pool']]} ({n(x['n_pool'])}) | {n(x['n_protons'])} |")
A("")
A("### 5.2 Rozbicie błędu względnego N_acc (A + C)")
A("")
A("| run | strumień | stat | pula centralna (½Δ, alternatywa) | pula protonów | LB / μ | ziarnistość nSelBL | niedomknięcie daleko od cięcia (↑) | przeciek: norm. Δ > 0,10 (↓) | lustrzany κ (↓, nierozstrzygnięty) | **razem ↑** | **razem ↓ bez lustrzanego** | razem ↓ z lustrzanym |")
A("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for run in RUNS:
    for st in ("afp_any", "pt4_L1AFP", "noalg_L1AFP", "nonAFP"):
        if run == "427929" and st == "noalg_L1AFP":
            continue
        b = R["runs"][run]["streams"][st]["budget_rel"]
        A(f"| {run} | {SN[st]} | {pc(b['stat'])} | {pc(b['pool_half_diff'])} ({PN[b['pool_alt']]}: ×{b['pool_alt_ratio']:.3f}) | {pc(b['protons'])} | {pc(b['lb_mu'])} | {pc(b['nsel_gran'])} | {pc(b['closure_up'])} ({b['closure_far_ratio']:.3f} ± {b['closure_far_err']:.3f}) | "
          f"{pc(b['leak_down'])} ({b['far_norm_ratio']:.3f} ± {b['far_norm_ratio_err']:.3f}) | " + (f"{pc(b['mirror_down'])} ({b['mirror_kappa']:.2f} ± {b['mirror_kappa_err']:.2f})" if b['mirror_kappa'] else "nie liczono") + f" | **{pc(b['total_up'])}** | **{pc(b['total_down_no_mirror'])}** | {pc(b['total_down'])} |")
A("")
A("Definicje składników:")
A("- pula centralna: połowa różnicy nominał ↔ pula alternatywna (reguła CMS-TOTEM 2002.12146); pula protonów: |nominał − protony z weta strumienia nie-AFP| (te nie przeszły przez tryger AFP);")
A("- LB / μ: największa zmiana przy 1, 4, 16 blokach LB i 8 blokach μ zamiast LB; ziarnistość nSelBL: największa zmiana przy klasach zgrubnych (5) i drobnych (24, 44). W 427929/sptrk_pt4 (140 zdarzeń CR) te składniki to głównie szum statystyczny;")
A("- niedomknięcie (↑): istotna część |obs/pred − 1| w teście „norm. Δ > 0,20 → 0,10–0,20”, czyli poza obszarem przecieku — dopuszcza, że gęstość przypadkowych rośnie ku małym ξ_cal szybciej niż w puli;")
A("- przeciek (↓): istotna część różnicy nominał ↔ normalizacja w Δ > 0,10. Znak znany (przeciek zawyża N_CR), wartość centralna nieprzesunięta, bo nominał był ustalony z góry;")
A("- lustrzany (↓): √((1 − κ)² + σ_κ²), zob. §6.3. Obowiązuje tylko, jeśli efekt nie zależy od strony; podany osobno.")
A("")
A("### 5.3 Zależność od LB, μ, nSelBL (428770, strumień AFP, f_acc ± stat)")
A("")
for dep, lab in (("lb", "blok LB"), ("mu", "μ"), ("nsel", "nSelBL")):
    A(f"| {lab} | " + " | ".join(f"{r['lo']:.3g}–{r['hi']:.3g}" if r['hi'] < 9000 else f"≥ {r['lo']:.3g}" for r in RAW["428770"]["afp_any"]["A"]["dep"][dep]) + " |")
    A("|---|" + "---|" * len(RAW["428770"]["afp_any"]["A"]["dep"][dep]))
    for s in "AC":
        A(f"| strona {s} | " + " | ".join(f"{r['f_acc']:.4f} ± {r['f_acc_err_stat']:.4f}" for r in RAW["428770"]["afp_any"][s]["dep"][dep]) + " |")
    A("")
A("f_acc nie zależy od bloku LB, rośnie z μ (0,035 → 0,072–0,077; składnik pile-up) i z nSelBL (0,04–0,05 → 0,085–0,10). Wynik zbiorczy jest ważony składem próbki; przy innej selekcji (score, pT) trzeba go przeliczyć, a nie przenosić.")
A("")
A("### 5.4 Gatunki × pT × okno masy (strumień AFP, A + C; kandydaci, bez fitu)")
A("")
A("N_tag = kandydaci w oknie w zdarzeniach z tagiem w SR; N_CR = kandydaci w zdarzeniach z CR; N_acc = Σ_c N_CR,c R_c z wagami = liczba kandydatów. Błąd: stat ⊕ systematyka względna strumienia z §5.2 (↑ / ↓ bez lustrzanego). „—”: N_CR < 5, estymata bez znaczenia. Okna: D⁺ 1,84–1,90 (SB 1,75–1,80, 1,94–2,00), D⁰ 1,835–1,895 (SB 1,75–1,80, 1,93–2,00) GeV, D* Δm 0,144–0,147 (SB 0,150–0,160) GeV. Per strona: `facc_afptrig.json`.")
for run in RUNS:
    S = R["runs"][run]["streams"]["afp_any"]; C = S["categories"]; b = S["budget_rel"]
    A("")
    A(f"**Run {run}**")
    A("")
    A("| gatunek | pT [GeV] | N_tag pik | N_CR pik | N_acc pik ± stat | f_acc pik ± stat + syst − syst | N_tag SB | N_CR SB | f_acc SB ± stat | pik/SB | pula alternatywna / nominał (pik) |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for sp in lib.SPECIES:
        for lo, hi in lib.PT_BINS:
            pk = C[f"{sp}|{lo:g}-{hi:g}|peak"]; sb = C[f"{sp}|{lo:g}-{hi:g}|sb"]
            def ff(c, syst):
                if c["n_CR"] < 5 or not c["f_acc"]:
                    return "—"
                t = f"{f4(c['f_acc'])} ± {f4(c['f_acc_err_stat'])}"
                return t + (f" + {f4(c['f_acc'] * b['syst_up'])} − {f4(c['f_acc'] * b['syst_down_no_mirror'])}" if syst else "")
            ps = "—"
            if pk["n_CR"] >= 5 and sb["n_CR"] >= 5 and pk["f_acc"] and sb["f_acc"]:
                r = pk["f_acc"] / sb["f_acc"]; ps = f"{r:.2f} ± {r * math.hypot(pk['f_acc_err_stat'] / pk['f_acc'], sb['f_acc_err_stat'] / sb['f_acc']):.2f}"
            A(f"| {sp} | {lo:g}–{hi:g} | {n(pk['n_tag_SR'])} | {n(pk['n_CR'])} | {n(pk['n_acc'], 1)} ± {n(pk['n_acc_err_stat'], 1)} | {ff(pk, True)} | {n(sb['n_tag_SR'])} | {n(sb['n_CR'])} | {ff(sb, False)} | {ps} | {(pk['n_acc_alt_pool'] / pk['n_acc']):.2f} |" if pk["n_acc"] > 0 else f"| {sp} | {lo:g}–{hi:g} | {n(pk['n_tag_SR'])} | {n(pk['n_CR'])} | 0 | — | {n(sb['n_tag_SR'])} | {n(sb['n_CR'])} | {ff(sb, False)} | — | — |")
A("")
A("- 428770: f_acc w kategoriach jest w granicach błędów stałe (0,05–0,06) i zgodne ze zbiorczym 0,051; pik i sidebandy zgodne (pik/SB 0,98–1,07 tam, gdzie błąd < 10 %). Odchylenia D⁰ 8–30 GeV (0,039) i D* 8–30 GeV to 2–3σ stat. Pula alternatywna (wszystkie nie-AFP) daje w kategoriach +9…22 %, więcej niż zbiorczo (+8 %).")
A("- 427929: policzalne tylko 3–5 GeV dla D⁺ i D⁰ (N_CR 102–236); pozostałe kategorie mają N_CR ≤ 26 i błąd rzędu samej wartości.")
A("- Sidebandy nie są próbką bez prawdziwych tagów; stosunek pik/SB jest informacyjny.")
A("")
A("## 6. Czego metoda nie daje")
A("")
A("### 6.1 Przeciek własnych protonów zdarzenia do obszaru kontrolnego")
A("")
A("Lokalnie nie ma MC SD z ξ_cal (`DIFF_model` ma tylko prawdę generatora bez ξ; `afp_ext` nie ma MC), więc przeciek oszacowano z danych: nadwyżka danych nad szablonem znormalizowanym w Δ > 0,10, jako ułamek prawdziwych tagów N_SR − N_acc (strony A + C), w tercylach μ:")
A("")
A("| run | strumień | ⟨μ⟩ | N_tag(SR) | nadwyżka w przerwie 0,02 < Δ ≤ 0,045 | nadwyżka w CR 0,045 < Δ ≤ 0,10 |")
A("|---|---|---|---|---|---|")
for run in RUNS:
    for st in ("afp_any", "nonAFP"):
        for i, r in enumerate(LK[run][st]):
            A(f"| {run} | {SN[st]} | {'wszystkie' if i == 0 else ''} {r['mu_mean']:.4f} | {n(r['n_sr'])} | {pc(r['gap']['frac_of_true'], 2)} ± {pc(r['gap']['frac_err'], 2)} | {pc(r['crlow']['frac_of_true'], 2)} ± {pc(r['crlow']['frac_err'], 2)} |")
A("")
l8 = LK["428770"]["afp_any"]
A(f"- 428770, strumień AFP: {pc(l8[0]['gap']['frac_of_true'] + l8[0]['crlow']['frac_of_true'], 1)} prawdziwych tagów ma 0,02 < Δ ≤ 0,10, z tego {pc(l8[0]['crlow']['frac_of_true'], 1)} w nominalnym CR. Ułamek rośnie z μ ({pc(l8[1]['gap']['frac_of_true'] + l8[1]['crlow']['frac_of_true'], 1)} → {pc(l8[3]['gap']['frac_of_true'] + l8[3]['crlow']['frac_of_true'], 1)} między skrajnymi tercylami), a w 427929 (μ dziesięć razy mniejsze) jest zgodny z zerem.")
A("- Hipoteza (niepotwierdzona): prawdziwy tag ⊕ depozyt kalorymetryczny oddziaływania pile-up bez zrekonstruowanego wierzchołka, który podnosi ξ_cal. Sama rozdzielczość ξ_cal nie zależałaby od μ. Po stronie ujemnej Δ nie ma z czego zbudować kształtu ogona: rozkład prawdziwych tagów jest tam szeroki i niesymetryczny (rys. 1).")
A("- Skutek dla metody: przeciek zawyża N_CR, więc N_acc. Przeciek poza Δ > 0,10 nie jest ograniczony tym oszacowaniem (szablon jest tam normalizowany). Te same zdarzenia to zarazem strata prawdziwych tagów na cięciu (ok. 2 % w 428770) — do uwzględnienia w wydajności, nie w tle.")
A("")
A("### 6.2 Obciążenie trygerem AFP")
A("")
x8 = EX["428770"]["afp_any"]
A(f"- **Zależność od ξ:** wydajność L1 rośnie z ξ o kilka procent (§2). Protony „niemożliwe” ze strumienia AFP i z weta strumienia nie-AFP mają to samo widmo: ⟨ξ⟩ = {x8['A']['xi_spectra']['impossible_stream']['mean']:.4f} ± {x8['A']['xi_spectra']['impossible_stream']['mean_err']:.4f} wobec {x8['A']['xi_spectra']['impossible_nonAFP']['mean']:.4f} ± {x8['A']['xi_spectra']['impossible_nonAFP']['mean_err']:.4f} (A), "
  f"{x8['C']['xi_spectra']['impossible_stream']['mean']:.4f} ± {x8['C']['xi_spectra']['impossible_stream']['mean_err']:.4f} wobec {x8['C']['xi_spectra']['impossible_nonAFP']['mean']:.4f} ± {x8['C']['xi_spectra']['impossible_nonAFP']['mean_err']:.4f} (C); zamiana puli protonów zmienia N_acc o {pc(b8['protons'])}. Metoda bierze protony z tego samego strumienia, więc zależność ε(ξ) znosi się z konstrukcji (toy: 1,01 ± 0,02). "
  "Nie znosi się natomiast dla prawdziwych tagów: wydajność trygera dla sygnału trzeba zmierzyć osobno (strumień centralny, wzorzec ALFA 2212.00664).")
A(f"- **Topologia:** zdarzenie z protonem po stronie przeciwnej odpala tryger samo, więc jest w strumieniu nadreprezentowane. Udział zdarzeń z protonem po stronie przeciwnej: CR {pc(x8['A']['opp_proton_frac']['CR'])} (A) i {pc(x8['C']['opp_proton_frac']['CR'])} (C) wobec {pc(x8['A']['opp_proton_frac']['pool'], 2)} w puli. "
  f"Usunięcie ich z danych i puli zmienia f_acc z {f4(x8['A']['est_nominal']['f_acc'])} na {f4(x8['A']['est_no_opp_proton']['f_acc'])} (A) i z {f4(x8['C']['est_nominal']['f_acc'])} na {f4(x8['C']['est_no_opp_proton']['f_acc'])} (C), czyli −2…3 % względnie; nieuwzględnione w budżecie osobno (mniejsze od składnika puli).")
A("- **Strona C** ma wydajność L1 ≈ 0,45: połowa tagów C jest zapisana tylko wtedy, gdy odpali coś innego. Frakcje f_acc stron A i C są mimo to zgodne (0,0513 i 0,0514).")
A("- Zdarzeń bez żadnego zapisanego bitu trygera (1,36 mln w 428770) i innych łańcuchów AFP ta analiza nie obejmuje.")
A("")
A("### 6.3 Zależność prawdopodobieństwa protonu przypadkowego od topologii zdarzenia (test lustrzany w strumieniu AFP)")
A("")
A("Zdarzenia z D, ξ_cal,s > 0,14 (proton po stronie s na pewno przypadkowy), bez protonu po stronie przeciwnej o. Rozkład ξ_cal,o danych wobec puli w tych samych komórkach, normalizacja w ξ_cal,o > 0,14:")
A("")
m8 = EX["428770"]["afp_any"]; ed = m8["A"]["mirror"]["edges"]
A("| ξ_cal,o | " + " | ".join(f"{ed[i]:g}–{ed[i + 1]:g}" for i in range(7)) + " | razem < 0,14 |")
A("|---|" + "---|" * 8)
for run in RUNS:
    mm = EX[run]["afp_any"]
    ob = [mm["A"]["mirror"]["obs"][i] + mm["C"]["mirror"]["obs"][i] for i in range(7)]; pr = [mm["A"]["mirror"]["pred"][i] + mm["C"]["mirror"]["pred"][i] for i in range(7)]
    bb = R["runs"][run]["streams"]["afp_any"]["budget_rel"]
    A(f"| {run}: obs/pred | " + " | ".join(f"{o_:.0f}/{p_:.0f} = {o_ / p_:.2f}" for o_, p_ in zip(ob, pr)) + f" | **{bb['mirror_kappa']:.2f} ± {bb['mirror_kappa_err']:.2f}** |")
A("")
A("- W strumieniu AFP efekt jest ten sam co w strumieniu nie-AFP (`TEMPLATE/RAPORT.md`: κ_L = 0,79 ± 0,04): zdarzenia z małym ξ_cal po stronie przeciwnej mają o ok. 20 % mniej protonów przypadkowych, a przy ξ_cal,o < 0,02 trzy razy mniej.")
A("- Jeśli efekt nie zależy od strony, szablon zawyża N_acc w SR o ok. 20 % (stąd kolumna „z lustrzanym”). Hipoteza alternatywna: proton przypadkowy pochodzi z dyfrakcyjnego oddziaływania pile-up, którego układ X zostawia energię po stronie przeciwnej do protonu (więc takie zdarzenie nie może mieć małego ξ_cal,o), a po stronie protonu prawie nic (przegląd: ⟨ξ_cal | tag⟩ − ⟨ξ_cal | brak⟩ = −0,0001 ± 0,0020). Wtedy SR nie jest obciążony. Dane tego nie rozstrzygają.")
A("")
A("### 6.4 Pozostałe")
A("")
A("- **Obszar ξ_cal,s ≲ 0,075 jest nietestowalny** (prawdziwe i przypadkowe tagi nieodróżnialne). Metoda kształtu przenosi to ograniczenie z v2 bez zmian: zakłada, że rozkład ξ_cal zdarzeń centralnych tagów przypadkowych jest poniżej 0,075 taki jak w puli. Toy: naruszenie daje bias −28…−35 % przy poprawnym domknięciu.")
A("- **Pula centralna** to zdarzenia z L1MBTS; strumień AFP nie wymaga MBTS. Dla łańcucha sptrk_pt4 walidacja daleko od cięcia daje 1,05 ± 0,02, dla noalg 1,78 ± 0,13.")
A("- **Zdarzenia z tagiem po obu stronach** (1 % tagów z D w 428770) liczone po obu stronach; regionu A∨C z cięciem na każdej otagowanej stronie nie odtworzono.")
A("- **Kandydaci bez score**: f_acc zależy od krotności i pT (§5.3, §5.4), po selekcji score trzeba przeliczyć.")
A("- Preskale oszacowane z danych, nie z bazy trygera; do samej frakcji f_acc nie wchodzą.")
A("- K2 nie jest niezależnym testem (§4.2). Jedyne testy niezależne od v2 to domknięcia w CR, toy i podpróbka 427929.")
A("")
A("## 7. Ocena")
A("")
A(f"- **Czy metoda jest używalna?** Dla łańcucha HLT_mb_sptrk_pt4_L1AFP_A_OR_C i sumy łańcuchów AFP — **warunkowo**. Według kryteriów zapisanych z góry K6 nie jest spełnione w 428770 (K1, K5 FAIL), bo obszar kontrolny tuż nad cięciem jest skażony przeciekiem własnych protonów. "
  f"Rozrzut wyniku między normalizacją w Δ > 0,045 i Δ > 0,10 to {pc(1 - b8['far_norm_ratio'], 0)} N_acc, czyli {f4(S8['f_acc'] - S8['f_acc_far_norm'])} w f_acc — przy f_acc ≈ 0,05 przesuwa to liczbę prawdziwych tagów o ≈ 0,7 %.")
A("- **Zalecenie (po obejrzeniu walidacji, do decyzji):** normalizować w Δ > 0,10 zamiast Δ > 0,045 i przyjąć połowę różnicy względem 0,045 jako systematykę; nominał z tego raportu zostawiono zgodny z zapisem sprzed analizy.")
A("- Dla HLT_noalg_L1AFP_A_OR_C w 428770 metoda w obecnej postaci nie jest zwalidowana; potrzebna pula centralna bez wymogu MBTS (tryger losowy — w tych danych 0 zdarzeń HLT_noalg_L1RD0_FILLED).")
A("- Do zrobienia przed użyciem w pomiarze: MC SD z pełną symulacją (kształt ogona Δ prawdziwych tagów, wpływ pile-up na ξ_cal), fit dwuskładnikowy szablon przypadkowych + ogon sygnału zamiast samej normalizacji w CR (wzorzec ALFA 2212.00664), rozstrzygnięcie efektu lustrzanego, próbka po score.")
A("")
A("## 8. Pliki")
A("")
A("| plik | zawartość |")
A("|---|---|")
A("| `facc_afptrig.json` | wynik: per run × strumień (afp_any, pt4_L1AFP, noalg_L1AFP, nonAFP): N_tag, N_CR, N_acc, f_acc, budżet błędu, walidacja K1/K4/K5/K6, strony, kategorie gatunek × pT × okno (A + C i per strona); K2, K3 |")
A("| `KRYTERIA_WALIDACJI.txt` | kryteria zapisane przed analizą |")
A("| `lib.py`, `core.py`, `analysis.py` | wczytanie danych (cache `cache/ev_*.npz`, 3,3 GB — do odtworzenia przez `prep.py`), estymator, bieg na danych → `results_raw_<run>.json` |")
A("| `describe.py`, `extra.py`, `leak.py`, `diag2d.py`, `strat_var.py` | opis próbki, test lustrzany / widma ξ / topologia, przeciek, diagnostyka (ξ_AFP, ξ_cal), ziarnistość nSelBL → `*.json`, `*.log` |")
A("| `toy.py` | toy z trygerem protonowym → `toy_results.json`, `toy.log` |")
A("| `make_report.py`, `report_md.py` | budżet błędu, kryteria, rysunki `fig1_delta.png`, `fig2_walidacja.png`, `fig3_toy.png`, ten raport |")
A("| `explore0.py`, `peek.py`, `prep.py` | pomocnicze |")
A("")
A("Odtworzenie: `prep.py` → `describe.py` → `analysis.py 427929`, `analysis.py 428770` → `extra.py`, `leak.py`, `diag2d.py`, `strat_var.py`, `toy.py` → `make_report.py` → `report_md.py` (Python `~/Documents/charm/.venv_h5/bin/python`, katalog roboczy `AFPTRIG/`).")
open("RAPORT.md", "w").write("\n".join(L) + "\n")
print("RAPORT.md", len(L), "linii")
