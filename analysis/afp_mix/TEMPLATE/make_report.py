#!/usr/bin/env python3
"""RAPORT.md z template.json + jump_tests.json (przez derive.py). Tylko formatowanie — żadnych nowych obliczeń na danych."""
from __future__ import annotations
import json, os, time
import numpy as np
from derive import derive, HERE

T = json.load(open(os.path.join(HERE, "template.json")))
S = derive()
SP = ("Dplus", "D0", "Dstar")
PTB = ("3-5", "5-8", "8-30", "3-30")
L = []
A = L.append


def f(x, n=4):
    return "—" if x is None or x != x else f"{x:.{n}f}"


def pc(x, n=1):
    return "—" if x is None or x != x else f"{100 * x:.{n}f} %"


def Rs(d, n=2):
    return "—" if d["R"] is None else f"{d['k']:.0f} / {d['pred']:.1f} = {d['R']:.{n}f} ± {d['err']:.{n}f}"


A("# Szablon tła przypadkowego tagu AFP (AFP-MIX v2) i test skoku częstości poniżej progu ξ_cal — WSTĘPNE")
A("")
A(f"Wygenerowano {time.strftime('%Y-%m-%d %H:%M')} przez `TEMPLATE/make_report.py` z `template.json` i `jump_tests.json`. Projekt B (ATLAS minbias 2022, runy 427929 i 428770, h5 v2.0, kandydaci bez score).")
A("Metoda bez zmian: `charmbkg/afp_mix.py` v2 (estymator A), systematyki z `runs/lxplus/out_v2/result_v2.json`. Kod metody i wyniki `out/`, `out_v2/` nietknięte.")
A("")
A("**Zakres ważności:** tylko zdarzenia z trygerów nie-AFP (GRL ∧ ≥ 1 zapisany łańcuch bez L1AFP ∧ nPileupVertices == 0); w 428770 to 4,4 % wszystkich tagów. Kandydaci po luźnej selekcji v2.0 (pT ≥ 3 GeV).")
A("**Czego tu nie ma:** uzysków sygnału D w zdarzeniach z tagiem, fitów masy, przekrojów, histogramów masy próbki tagowanej. N_tag to liczba kandydatów (sygnał + kombinatoryka) w oknie masy.")
A("")
A("## 1. Szablon (ZADANIE 1)")
A("")
A("Pliki: `template.json` (wszystko), `template.csv` (liczności i f_acc, wiersz = run × strona × gatunek × pT), `template_shapes.csv` (kształty masy, wiersz = bin masy). Schemat pól: `README.md`.")
A("")
A("- N_tag: kandydaci w oknie masy (D⁺ 1,84–1,90; D⁰ 1,835–1,895 GeV; D* Δm 0,1440–0,1470 GeV) w zdarzeniach z tagiem w regionie po cięciu ξ_cal ≤ ξ_AFP + 0,02.")
A("- N_acc = Σ p_acc po kandydatach w oknie (wszyscy biorcy, mieszanie nominalne v2); f_acc = N_acc / N_tag.")
A("- Błędy: stat = bootstrap otagowanych dawców (B = 500) ⊕ Poisson N_tag; pula, ekstrapolacja, ½Δ metod = wspólne dla runu błędy względne N_acc z v2 (tabela niżej), stosowane do każdego wiersza i każdej strony.")
A("- Przedziału pT 2–3 GeV **nie ma w danych** (slim: pT ≥ 3,0; D* ≥ 3,11 GeV). Przedziały: 3–5, 5–8, 8–30 oraz 3–30 (suma). Kandydaci z pT ≥ 30 GeV pominięci (0,005 % i 0,12 %).")
A("")
A("Systematyki względne runu (z v2) i kontrola zgodności z `out_v2` (zdarzenia, A∨C):")
A("")
A("| run | N_tag (zdarzenia) | N_acc | stat (B = 500) | pula | ekstrapolacja | ½Δ metod | syst razem |")
A("|---|---|---|---|---|---|---|---|")
for s in T["samples"]:
    c, y = s["check_vs_v2"]["AorC"], s["syst_rel_run"]
    A(f"| {s['run']} | {c['n_tag']} | {c['n_acc']:.1f} ± {c['n_acc_err_stat']:.1f} | {pc(c['n_acc_err_stat'] / c['n_acc'])} | {pc(y['pool'])} | {pc(y['extrap'])} | {pc(y['half'])} | {pc(y['syst'])} |")
A("")
A("N_tag i N_acc zgadzają się z `out_v2` co do zdarzenia (kontrola w kodzie); wiersz kontrolny pT 3–150 w `template.csv` odtwarza tabele gatunków z `RAPORT_DANE_v2.md`.")
A("")
for s in T["samples"]:
    rows = {(r["region"], r["species"], r["pt_bin"]): r for r in s["rows"]}
    A(f"### Run {s['run']}: N_tag, N_acc, f_acc")
    A("")
    A("f_acc ± stat ± pula ± ekstrapolacja ± ½Δ metod (błędy bezwzględne).")
    A("")
    A("| gatunek | pT [GeV] | strona | N_tag | N_acc ± stat | f_acc | ± stat | ± pula | ± ekstrap. | ± ½Δ | ± razem |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for sp in SP:
        for pt in PTB:
            for reg in ("A", "C", "AorC"):
                r = rows[(reg, sp, pt)]
                nd = 2 if r["n_acc"] < 10 else 1
                A(f"| {sp} | {pt} | {reg.replace('AorC', 'A∨C')} | {r['n_tag']:.0f} | {r['n_acc']:.{nd}f} ± {r['n_acc_err_stat']:.{nd}f} | {f(r['f_acc'])} | {f(r['f_acc_err_stat'])} | {f(r['f_acc_err_pool'])} | "
                  f"{f(r['f_acc_err_extrap'])} | {f(r['f_acc_err_half'])} | {f(r['f_acc_err_tot'])} |")
    A("")
A("### Kształt masy składowej przypadkowej a kształt wszystkich kandydatów bez wymogu tagu")
A("")
A("Składowa przypadkowa: histogram Σ p_acc po kandydatach (biny 10 MeV dla D⁺ 1,72–2,02 i D⁰ 1,75–2,00; 1 MeV dla Δm D* 0,139–0,160). Odniesienie: wszyscy kandydaci biorców bez wymogu tagu "
  "(zdarzenia z tagiem to < 1 % tej próbki). Test: d_i = h_acc,i − (Σh_acc)·s_i, χ² = dᵀC⁻¹d, ndf = liczba binów − 1.")
A("")
nul = [r["shape_chi2_null_perm"] / r["shape_ndf"] for s in T["samples"] for r in s["rows"]]
ov = {s["run"]: np.mean([r["shape_overlap_frac"] for r in s["rows"]]) for s in T["samples"]}
A("Kowariancja C = C_boot + C_samp − C_overlap:")
A("- C_boot: bootstrap otagowanych dawców przy ustalonych biorcach (B = 500). To jest korelacja wag z jednej puli dawców: jeden dawca daje wagę wszystkim biorcom swojej komórki (LB × nSelBL), więc szum wag przenosi się między komórkami o różnym kształcie masy.")
A("- C_samp: próbkowanie zdarzeń-biorców przy zrealizowanych wagach, Σ_t (p_t − p̄)²·g_t g_tᵀ (g_t = wektor kandydatów zdarzenia w binach minus oczekiwanie).")
A(f"- C_overlap = Σ_t Var_boot(p_t)·g_t g_tᵀ: część C_samp, którą C_boot już zawiera (E[(p − p̄)²] = (π − π̄)² + Var p). Bez tej poprawki wariancja jest liczona podwójnie: C_overlap/C_samp ≈ {ov[427929]:.2f} (427929, mało dawców) i {ov[428770]:.2f} (428770).")
A("")
A(f"Warianty: χ²_cons z C_boot + C_samp (wariancja zawyżona → dolna granica istotności); χ²_samp z samym C_samp (H0 ścisłe: masa niezależna od wag), kalibrowane permutacją p_acc między biorcami: χ²/ndf = {np.mean(nul):.2f} ± {np.std(nul):.2f} w {len(nul)} kategoriach. "
  "Kontrola χ² nominalnego: w 427929 dla D⁰ i D* (gdzie efektu nie widać) średnie χ²/ndf = "
  f"{np.mean([r['shape_chi2'] / r['shape_ndf'] for s in T['samples'] if s['run'] == 427929 for r in s['rows'] if r['species'] != 'Dplus' and r['pt_bin'] != '3-150_v2check']):.2f} (cons: "
  f"{np.mean([r['shape_chi2_cons'] / r['shape_ndf'] for s in T['samples'] if s['run'] == 427929 for r in s['rows'] if r['species'] != 'Dplus' and r['pt_bin'] != '3-150_v2check']):.2f}). "
  "W kategoriach z kilkoma kandydatami z tagiem (D* 8–30 GeV w 427929) przybliżenie gaussowskie χ² nie jest wiarygodne.")
A("")
A("Kolumny: udział okna piku = (część h_acc w oknie piku)/(część próbki nietagowanej w oknie piku); nachylenie = liniowe nachylenie stosunku kształtów acc/nietagowane względem masy [1/GeV].")
A("")
A("| run | gatunek | pT | strona | χ²/ndf | p | χ²_cons/ndf | p_cons | χ²_samp/ndf | udział okna piku | nachylenie [1/GeV] |")
A("|---|---|---|---|---|---|---|---|---|---|---|")
for s in T["samples"]:
    rows = {(r["region"], r["species"], r["pt_bin"]): r for r in s["rows"]}
    for sp in SP:
        for pt in PTB:
            for reg in ("AorC",) if pt != "3-30" else ("A", "C", "AorC"):
                r = rows[(reg, sp, pt)]
                ns = 2 if sp != "Dstar" else 1
                A(f"| {s['run']} | {sp} | {pt} | {reg.replace('AorC', 'A∨C')} | {r['shape_chi2']:.1f}/{r['shape_ndf']} | {r['shape_p']:.2g} | {r['shape_chi2_cons']:.1f}/{r['shape_ndf']} | {r['shape_p_cons']:.2g} | {r['shape_chi2_samp']:.1f}/{r['shape_ndf']} | "
                  f"{r['peakfrac_ratio']:.3f} ± {r['peakfrac_ratio_err']:.3f} | {r['shape_slope_per_GeV']:.{ns}f} ± {r['shape_slope_err']:.{ns}f} |")
A("")
r1 = {(s["run"], r["region"], r["species"], r["pt_bin"]): r for s in T["samples"] for r in s["rows"]}
a, b, c0 = r1[(428770, "AorC", "Dplus", "3-30")], r1[(427929, "AorC", "Dplus", "3-30")], r1[(428770, "AorC", "D0", "3-30")]
ds = r1[(428770, "AorC", "Dstar", "3-30")]
A("**Wynik.**")
A(f"- D⁺: kształt składowej przypadkowej **różni się** od kształtu próbki nietagowanej. 428770: χ²/ndf = {a['shape_chi2']:.0f}/{a['shape_ndf']} (zachowawczo {a['shape_chi2_cons']:.0f}/{a['shape_ndf']}), "
  f"nachylenie {a['shape_slope_per_GeV']:.3f} ± {a['shape_slope_err']:.3f} /GeV, czyli około ±{100 * abs(a['shape_slope_per_GeV']) * 0.15:.1f} % na krańcach zakresu 1,72–2,02 GeV (więcej przy małej masie); to samo w każdym przedziale pT i po obu stronach. "
  f"427929: nachylenie {b['shape_slope_per_GeV']:.2f} ± {b['shape_slope_err']:.2f} /GeV ({abs(b['shape_slope_per_GeV']) / b['shape_slope_err']:.1f}σ, ten sam znak i wielkość), χ²/ndf = {b['shape_chi2']:.1f}/{b['shape_ndf']} (p = {b['shape_p']:.2f}; zachowawczo p = {b['shape_p_cons']:.2f}).")
A(f"- D⁰: w 428770 mniejsze nachylenie tego samego znaku ({c0['shape_slope_per_GeV']:.3f} ± {c0['shape_slope_err']:.3f} /GeV; χ²/ndf = {c0['shape_chi2']:.0f}/{c0['shape_ndf']}, p = {c0['shape_p']:.1g}; zachowawczo p = {c0['shape_p_cons']:.2g}), rosnące z pT; w 427929 zgodność.")
A(f"- D* (Δm): χ² zgodne w obu runach (A∨C: najmniejsze p = "
  f"{min(r['shape_p'] for s in T['samples'] for r in s['rows'] if r['species'] == 'Dstar' and r['pt_bin'] != '3-150_v2check' and r['region'] == 'AorC'):.2f}); "
  f"w 428770 słaba wskazówka nachylenia ({ds['shape_slope_per_GeV']:.1f} ± {ds['shape_slope_err']:.1f} /GeV, {abs(ds['shape_slope_per_GeV']) / ds['shape_slope_err']:.1f}σ) i udziału okna piku {ds['peakfrac_ratio']:.3f} ± {ds['peakfrac_ratio_err']:.3f}.")
A(f"- Udział okna piku dla D⁺ i D⁰ jest zgodny z 1 (428770, 3–30 GeV: {a['peakfrac_ratio']:.4f} ± {a['peakfrac_ratio_err']:.4f} i {c0['peakfrac_ratio']:.4f} ± {c0['peakfrac_ratio_err']:.4f}): okno leży w środku zakresu, więc nachylenie nie zmienia "
  "f_acc w oknie piku. Zmienia natomiast przewidywanie, jeśli tło przypadkowe normalizuje się w sidebandach po jednej stronie piku (do kilku procent).")
A("- Zalecenie: jako kształt masy tła przypadkowego brać histogram ważony p_acc z `template_shapes.csv`, nie kształt próbki nietagowanej.")
A("")
A("Hipoteza (niesprawdzana tutaj): p_acc faworyzuje zdarzenia o małym ξ_cal i małej krotności, a kształt masy kombinatoryki D⁺ (trzy tory) zależy od krotności.")
A("")

# =============================================================================================== ZADANIE 2
A("## 2. Test skoku częstości przypadkowych poniżej progu ξ_cal dawców (ZADANIE 2)")
A("")
A("Pytanie: czy częstość protonów przypadkowych u biorców z małym ξ_cal po stronie tagu (obszar użycia, ξ_cal ≲ 0,12) jest taka jak u dawców (ξ_cal > 0,14). Obszaru użycia nie da się sprawdzić wprost, bo tagi prawdziwe i przypadkowe są tam nieodróżnialne.")
A("")
A("Profil obszaru użycia (wagi par „acc” z mieszania nominalnego; s = strona tagu, o = strona przeciwna):")
A("")
A("| run | strona | ⟨ln ξ_cal,s⟩ | część wag z ξ_cal,s < 0,06 | część wag z ξ_cal,o > 0,14 | ⟨nSelBL⟩ użycie | ⟨nSelBL⟩ pula z D |")
A("|---|---|---|---|---|---|---|")
J = {r["run"]: r for r in json.load(open(os.path.join(HERE, "jump_tests.json")))["results"]}
for s in T["samples"]:
    for sd, u in s["usage_profile"].items():
        H = np.array(u["w2d"])
        A(f"| {s['run']} | {sd} | {u['mean_ln_xical_same']:.2f} | {pc(H[:6].sum() / H.sum(), 0)} | {pc(H[:, 7:].sum() / H.sum(), 0)} | {u['mean_nsel']:.1f} | {J[s['run']]['test_i']['sides'][sd]['mean_nsel_pool_hasD']:.1f} |")
A("")
A("Wszystkie testy niżej porównują liczbę obserwowanych tagów z przewidywaniem z tych samych warstw co nominał (LB × nSelBL × „z D”), R = obs/pred. „z D” = zdarzenia z ≥ 1 kandydatem (klasa biorców); „wszystkie” = z D i bez D (każde względem własnej klasy).")
A("")
A("### 2.1 Test (i): strona przeciwna i estymator „lustrzany”")
A("")
A("**Co mierzy.** W zdarzeniach z ξ_cal,s > 0,14 (własny proton po stronie s wykluczony) częstość tagu po stronie s w funkcji ξ_cal,o strony przeciwnej, do najmniejszych ξ_cal,o. "
  "To zależność częstości przypadkowej od typu zdarzenia (jedna strona z małym ξ_cal) przy ustalonych LB i nSelBL.")
A("**Czego nie mierzy.** Efektów związanych ze stroną protonu przy małym ξ_cal,s. Przeniesienie na obszar użycia wymaga założenia, że częstość nie zależy od tego, po której stronie jest małe ξ_cal („ślepota na stronę”); "
  "ponadto tylko 29–44 % wag obszaru użycia ma ξ_cal,o > 0,14 (dokładne odbicie lustrzane zdarzeń testowych), reszta ma małe ξ_cal po obu stronach.")
A("")
for run, x in S.items():
    pr = x["i_profile_hasD"]
    ed = pr["edges"]
    A(f"Run {run}, zdarzenia z D, strony A + C: R w przedziałach ξ_cal,o")
    A("")
    A("| ξ_cal,o | " + " | ".join(f"{ed[i]:g}–{ed[i + 1] if ed[i + 1] is not None else '∞'}" for i in range(len(ed) - 1)) + " |")
    A("|---|" + "---|" * (len(ed) - 1))
    A("| obs | " + " | ".join(f"{b_['k']:.0f}" for b_ in pr["bins"]) + " |")
    A("| R | " + " | ".join("—" if b_["R"] is None else f"{b_['R']:.2f} ± {b_['err']:.2f}" for b_ in pr["bins"]) + " |")
    A("")
A("| run | próbka | R(ξ_cal,o < 0,14) | R(ξ_cal,o > 0,14) | stosunek mały/duży |")
A("|---|---|---|---|---|")
for run, x in S.items():
    for lab, nm in (("hasD", "z D"), ("hasD_npo0", "z D, bez protonu po stronie o"), ("all", "wszystkie")):
        pr = x[f"i_profile_{lab}"]
        A(f"| {run} | {nm} | {Rs(pr['low'], 3)} | {Rs(pr['high'], 3)} | {pr['low_over_high'][0]:.3f} ± {pr['low_over_high'][1]:.3f} |")
A("")
A("Estymator lustrzany: nominał z dodatkową klasą warstw — dawca (ξ_cal,s > 0,14) ma ξ_cal,o w tym samym przedziale (0; 0,02; …; 0,14) co ξ_cal,s biorcy. "
  "Błąd stosunku ze wspólnego bootstrapu dawców. „bez p_o” = dawcy bez protonu po stronie przeciwnej (u biorców proton po stronie o ma < 0,12 % wag).")
A("")
A("| run | N_acc nominał | nominał, dawcy bez p_o / nominał | lustrzany / nominał | **lustrzany bez p_o / nominał bez p_o (κ_L)** | strona A | strona C | górna granica κ_L (95 % CL) |")
A("|---|---|---|---|---|---|---|---|")
for run, x in S.items():
    m, n0 = x["mirror_npo0"], x["nominal_npo0"]
    A(f"| {run} | {x['nominal']['n_acc']:.1f} ± {x['nominal']['err']:.1f} | {n0['ratio']:.3f} ± {n0['ratio_err_corr']:.3f} | {x['mirror_all']['ratio']:.3f} ± {x['mirror_all']['ratio_err_corr']:.3f} | "
      f"**{m['ratio_to_nominal_npo0']:.3f} ± {m['ratio_to_nominal_npo0_err_corr']:.3f}** | {m['ratio_A_to_nominal_npo0']:.2f} ± {m['ratio_A_to_nominal_npo0_err_corr']:.2f} | "
      f"{m['ratio_C_to_nominal_npo0']:.2f} ± {m['ratio_C_to_nominal_npo0_err_corr']:.2f} | {x['proposal']['kappa_upper95']:.2f} |")
A("")
p9, p8 = S["427929"], S["428770"]
A("**Wynik.** Częstość tagu w zdarzeniach z małym ξ_cal po stronie przeciwnej jest **mniejsza**, nie większa: "
  f"stosunek mały/duży = {p8['i_profile_hasD']['low_over_high'][0]:.2f} ± {p8['i_profile_hasD']['low_over_high'][1]:.2f} (428770) i {p9['i_profile_hasD']['low_over_high'][0]:.2f} ± {p9['i_profile_hasD']['low_over_high'][1]:.2f} (427929); "
  f"κ_L = {p8['proposal']['kappa_mirror']:.2f} ± {p8['proposal']['kappa_mirror_err']:.2f} i {p9['proposal']['kappa_mirror']:.2f} ± {p9['proposal']['kappa_mirror_err']:.2f}. "
  "Skoku w górę o mechanizmie ślepym na stronę dane nie dopuszczają (górne granice κ_L w tabeli; scenariusz toy ×2 dałby R ≈ 2 przy małym ξ_cal,o).")
A("")
A("Test symetrii stron (czy zależność jest ślepa na stronę): zdarzenia puli z ξ_cal > 0,14 po obu stronach, nachylenie dR/d ln ξ_cal względem strony tagu i strony przeciwnej (przedziały 0,14–0,2–0,3–0,5; te same zdarzenia, błąd różnicy w przybliżeniu niezależnych nachyleń).")
A("")
A("| run | próbka | nachylenie, strona tagu | nachylenie, strona przeciwna | różnica | z |")
A("|---|---|---|---|---|---|")
for run, x in S.items():
    for lab, nm in (("hasD", "z D"), ("all", "wszystkie")):
        y = x[f"ib_symmetry_{lab}"]
        A(f"| {run} | {nm} | {y['same']['slope']:.2f} ± {y['same']['slope_err']:.2f} | {y['opp']['slope']:.2f} ± {y['opp']['slope_err']:.2f} | {y['diff']:.2f} ± {y['diff_err']:.2f} | {y['z']:.1f} |")
A("")
y = p8["ib_symmetry_hasD"]
A(f"W 428770 częstość rośnie z ξ_cal strony przeciwnej ({y['opp']['slope']:.2f} ± {y['opp']['slope_err']:.2f}), a nie zależy od ξ_cal strony tagu ({y['same']['slope']:.2f} ± {y['same']['slope_err']:.2f}); różnica {y['z']:.1f}σ. "
  "To wskazówka (nie dowód), że efekt jest związany ze stroną przeciwną do protonu, a więc że założenie ślepoty na stronę nie jest spełnione i κ_L nie przenosi się wprost na obszar użycia. W 427929 test nie ma mocy.")
A("")
A("### 2.2 Test (ii): protony kinematycznie niemożliwe")
A("")
A("**Co mierzy.** Dokładnie jeden proton po stronie s z ξ_p w wycinku [a, b] (0,02–0,03 poniżej okna tagu oraz 0,03–0,04, 0,04–0,05, 0,05–0,06, 0,06–0,08, 0,08–0,10) i ξ_cal,s > b + 0,02 + m. "
  "Taki proton nie pochodzi z tego oddziaływania, jeśli m przewyższa rozdzielczość ξ_cal. Klasy marginesu m: bliski 0,010–0,025, średni 0,025–0,040, daleki ≥ 0,040. Zakres ξ_cal,s: 0,06–0,14.")
A("**Czego nie mierzy.** ξ_cal,s < 0,06 (najmniejsze ξ_p w danych to ≈ 0,02; 55 % wag obszaru użycia leży poniżej 0,06) oraz, przy dalekim marginesie, ξ_cal,s < 0,08. Margines bliski jest zanieczyszczony własnymi protonami (rozdzielczość ξ_cal).")
A("")
A("| run | próbka | bliski | średni | daleki | bliski − daleki |")
A("|---|---|---|---|---|---|")
for run, x in S.items():
    for lab, nm in (("hasD", "z D"), ("all", "wszystkie")):
        t = x["ii"][lab]
        A(f"| {run} | {nm} | " + " | ".join(Rs(b_) for b_ in t["by_margin"]) + f" | {t['near_minus_far'][0]:.2f} ± {t['near_minus_far'][1]:.2f} |")
A("")
A("R w przedziałach ξ_cal,s (zdarzenia z D; obs i R ± błąd):")
A("")
for run, x in S.items():
    t = x["ii"]["hasD"]
    e = x["ii"]["edges_xical_same"]
    A(f"| {run}: ξ_cal,s | bliski | średni | daleki |")
    A("|---|---|---|---|")
    for i in range(len(e) - 1):
        A(f"| {e[i]:.2f}–{e[i + 1]:.2f} | " + " | ".join("—" if b_["R"] is None or b_["k"] + b_["pred"] < 1 else f"{b_['k']:.0f}: {b_['R']:.2f} ± {b_['err']:.2f}" for b_ in t["by_xical_margin"][i]) + " |")
    A("")
A("R w funkcji ξ_cal strony przeciwnej (zdarzenia z D):")
A("")
A("| run | ξ_cal,o | bliski | średni | daleki |")
A("|---|---|---|---|---|")
for run, x in S.items():
    t = x["ii"]["hasD"]
    for q, nm in enumerate(("< 0,06", "0,06–0,14", "> 0,14")):
        A(f"| {run} | {nm} | " + " | ".join(f"{b_['k']:.0f}: {b_['R']:.2f} ± {b_['err']:.2f}" for b_ in t["by_opp_margin"][q]) + " |")
A("")
t8, t9 = p8["ii"]["hasD"], p9["ii"]["hasD"]
A(f"**Wynik.** Przy dalekim marginesie (ξ_cal,s 0,09–0,14) R = {t8['by_margin'][2]['R']:.2f} ± {t8['by_margin'][2]['err']:.2f} (428770) i {t9['by_margin'][2]['R']:.2f} ± {t9['by_margin'][2]['err']:.2f} (427929): "
  f"brak skoku w tym zakresie; górne granice 95 % CL: {p8['proposal']['R_far_upper95']:.2f} i {p9['proposal']['R_far_upper95']:.2f}.")
A(f"Nadwyżka 428770 tuż pod progiem (przegląd: 1,14 ± 0,06; B(ii): 1,12 ± 0,06) siedzi w marginesie bliskim ({t8['by_margin'][0]['R']:.2f} ± {t8['by_margin'][0]['err']:.2f}) "
  f"i tam rośnie z ξ_cal strony przeciwnej ({t8['by_opp_margin'][0][0]['R']:.2f} ± {t8['by_opp_margin'][0][0]['err']:.2f} → {t8['by_opp_margin'][2][0]['R']:.2f} ± {t8['by_opp_margin'][2][0]['err']:.2f}), "
  "a przy dalekim marginesie jej nie ma przy tym samym ξ_cal,s. Skok częstości przypadkowych dawałby to samo R dla każdego ξ_p; "
  "obraz pasuje do przecieku własnych protonów przez rozdzielczość ξ_cal (prawdziwa dyfrakcja po stronie s ma układ X po stronie o, czyli duże ξ_cal,o). To interpretacja: różnica bliski − daleki ma "
  f"{t8['near_minus_far'][0] / t8['near_minus_far'][1]:.1f}σ.")
A("")
A("### 2.3 Test (iii): zależność od nSelBL, μ, LB")
A("")
A("**Co mierzy.** (a) Czy niedobór z testu (i) jest skupiony w jakiejś klasie nSelBL, przedziale μ lub części runu. (b) Czy warstwy nSelBL nominału są dość drobne (użycie ma ⟨nSelBL⟩ ≈ 20–27, pula 38–47).")
A("**Czego nie mierzy.** Niczego poniżej progu po stronie tagu.")
A("")
A("Stosunek R(ξ_cal,o < 0,12)/R(ξ_cal,o ≥ 0,12), zdarzenia z D, strony razem:")
A("")
for run, x in S.items():
    for var, nm in (("nsel", "nSelBL"), ("mu", "μ (tercyle)"), ("lb", "LB (tercyle)")):
        d = x["iii_dep"][var]
        e = d["edges"]
        cells = []
        for i, (v, lo_) in enumerate(zip(d["hasD_low_over_high"], d["hasD_low"])):
            if v[0] is None or lo_["k"] < 5:
                continue
            last = i == len(e) - 2
            lab = (f"≥ {e[i]:.0f}" if last and var != "nsel" else f"{e[i]:.0f}–{e[i + 1]:.0f}") if var != "mu" else (f"≥ {e[i]:.3f}" if last else f"{e[i]:.3f}–{e[i + 1]:.3f}")
            cells.append(f"{lab}: {v[0]:.2f} ± {v[1]:.2f}")
        A(f"- {run}, {nm}: " + "; ".join(cells))
A("")
A("| run | drobne warstwy nSelBL (19 brzegów) / nominał | dawcy bez protonu po stronie o / nominał | v2: warstwy μ zamiast LB |")
A("|---|---|---|---|")
V2 = {r["run_number"]: r for r in json.load(open(os.path.join(os.path.dirname(HERE), "runs", "lxplus", "out_v2", "result_v2.json")))["results"]}
for run, x in S.items():
    v = V2[int(run)]["variants"]["strata_mu"]
    A(f"| {run} | {x['act_fine']['ratio']:.3f} ± {x['act_fine']['ratio_err_corr']:.3f} | {x['nominal_npo0']['ratio']:.3f} ± {x['nominal_npo0']['ratio_err_corr']:.3f} | {v['ratio']:.3f} ± {v['delta_err']:.3f} |")
A("")
k1 = sum(sum(J[428770]["test_i"]["sides"][sd]["tables"]["hasD_npo1"]["k"]) for sd in "AC")
q1 = sum(sum(J[428770]["test_i"]["sides"][sd]["tables"]["hasD_npo1"]["pred"]) for sd in "AC")
A("**Wynik.** Niedobór z testu (i) jest podobny w tercylach μ i LB. W klasach nSelBL: 427929 bez trendu w granicach błędu; 428770 silniejszy przy nSelBL < 20 (0,33 ± 0,10 i 0,64 ± 0,12), dla nSelBL ≥ 20 w zakresie 0,78–0,95. "
  "Drobniejsze warstwy nSelBL zmieniają N_acc o +1,1 ± 0,7 % (428770) i +3 ± 5 % (427929), czyli resztkowa zależność od aktywności wewnątrz klas jest mała. "
  f"Dawcy z protonem po stronie przeciwnej mają częstość tagu wielokrotnie większą od przewidywanej (428770, z D, strony razem: {k1:.0f} tagów wobec {q1:.1f} przewidywanych — ten sam efekt co G9); ich usunięcie obniża N_acc o 1,8 ± 0,5 % i 3,3 ± 2,8 %.")
A("Wynik poboczny: przy ustalonych warstwach LB R rośnie z μ zdarzenia o czynnik około 1,6–2 między skrajnymi tercylami (zarówno przy małym, jak i dużym ξ_cal,o). "
  "actual_mu zmienia się głównie w obrębie LB (rozrzut względny w LB: 0,86 w 427929 i 0,53 w 428770; między LB: 0,09 i 0,03), więc warstwy LB nie warstwują μ. "
  "N_acc jest na to niewrażliwe tylko dlatego, że biorcy i dawcy mają w LB ten sam rozkład μ (wariant v2 „warstwy μ”: 1,007 ± 0,003 i 1,035 ± 0,017).")
A("")
A("### 2.4 Wniosek i propozycja systematyki")
A("")
A("Wyniki (pomiary):")
A(f"1. Mechanizm ślepy na stronę: skok w górę wykluczony; κ_L = {p8['proposal']['kappa_mirror']:.2f} ± {p8['proposal']['kappa_mirror_err']:.2f} (428770), {p9['proposal']['kappa_mirror']:.2f} ± {p9['proposal']['kappa_mirror_err']:.2f} (427929); "
  f"granice 95 % CL: κ_L < {p8['proposal']['kappa_upper95']:.2f} i < {p9['proposal']['kappa_upper95']:.2f}.")
A(f"2. Strona tagu, ξ_cal 0,09–0,14: R = {p8['proposal']['R_far_hasD']:.2f} ± {p8['proposal']['R_far_hasD_err']:.2f} i {p9['proposal']['R_far_hasD']:.2f} ± {p9['proposal']['R_far_hasD_err']:.2f}; skok < {100 * (p8['proposal']['R_far_upper95'] - 1):.0f} % i < {100 * (p9['proposal']['R_far_upper95'] - 1):.0f} % (95 % CL).")
A("3. Nadwyżka 1,14 ± 0,06 w 428770 nie wygląda na skok częstości przypadkowych (2.2).")
A("")
A("Hipotezy (niepotwierdzone):")
A("- Jeśli zależność z testu (i) jest ślepa na stronę, nominał **zawyża** N_acc o około 20 % w obu runach. Test symetrii (2,4σ w 428770, zdarzenia z D) przemawia przeciw, ale nie rozstrzyga.")
A("- Jeśli efekt jest związany ze stroną przeciwną (np. proton przypadkowy z oddziaływania, które zostawia energię w kalorymetrze po stronie o i wypycha zdarzenie z małych ξ_cal,o), nominał jest obciążony mniej; "
  "znak ewentualnej poprawki (depozyt także po stronie tagu → mniejsza wydajność cięcia) jest ten sam. Zależność niedoboru od μ, która by to potwierdziła, nie jest widoczna w granicach błędu.")
A("- Mechanizm działający tylko po stronie tagu i tylko dla ξ_cal,s < 0,06 pozostaje **nieograniczony przez dane** (bez obszaru sygnału nie ma próbki). Za jego brakiem przemawia tylko argument: proton niezależny nie „wie”, po której stronie jest przerwa.")
A("")
A("**PROPOZYCJA** (kod metody bez zmian): składnik „ekstrapolacja” zastąpić asymetryczną obwiednią z pomiarów — w dół √((1 − κ_L)² + σ²), w górę √((R_daleki − 1)² + σ²); bez korekty wartości centralnej. Pozostałe składniki (stat, pula, ½Δ metod) bez zmian.")
A("")
A("| run | f_acc | v2: ekstrapolacja | v2: razem | **propozycja: zależność od ξ_cal** | **propozycja: razem** | **f_acc, propozycja** | v2 |")
A("|---|---|---|---|---|---|---|---|")
for run, x in S.items():
    p = x["proposal"]
    A(f"| {run} | {p['f_acc']:.4f} | ± {pc(p['v2']['extrap'])} | ± {pc(p['v2']['total'])} | +{pc(p['xical_syst_up'])} / −{pc(p['xical_syst_down'])} | +{pc(p['total_up'])} / −{pc(p['total_down'])} | "
      f"**{p['f_acc']:.4f} +{p['f_err_up']:.4f} / −{p['f_err_down']:.4f}** | {p['f_acc']:.4f} ± {p['v2']['f_err']:.4f} |")
A("")
A(f"- 427929: dane **pozwalają zawęzić** błąd: ± {pc(p9['proposal']['v2']['total'])} → +{pc(p9['proposal']['total_up'])} / −{pc(p9['proposal']['total_down'])}. Ekstrapolacja v2 (47 %) była zdominowana błędem statystycznym nachylenia z ~440 tagów puli na ramieniu 1,45 w ln ξ_cal; estymator lustrzany mierzy tę samą wielkość bez ekstrapolacji.")
A(f"- 428770: dane **nie pozwalają zawęzić** błędu. Składnik w dół rośnie z 10,4 % do {pc(p8['proposal']['xical_syst_down'])}, bo test (i) pokazuje istotny ({(1 - p8['proposal']['kappa_mirror']) / p8['proposal']['kappa_mirror_err']:.1f}σ) efekt, którego ekstrapolacja liniowa po stronie tagu nie widzi; w górę maleje do {pc(p8['proposal']['xical_syst_up'])}. "
  f"Błąd całkowity dalej zdominowany przez ½Δ metod (B(i), 25,8 %): +{pc(p8['proposal']['total_up'])} / −{pc(p8['proposal']['total_down'])} wobec ± {pc(p8['proposal']['v2']['total'])}.")
A("- Wariant symetryczny (większa z obu stron): ± "
  f"{pc(p9['proposal']['total_sym'])} (427929), ± {pc(p8['proposal']['total_sym'])} (428770).")
A("")
A("## 3. Zastrzeżenia")
A("")
A("- Wszystko WSTĘPNE; zakres tylko trygery nie-AFP; kandydaci bez score. Po selekcji score f_acc i kształty trzeba przeliczyć.")
A("- Systematyki w szablonie są wspólne dla runu (z regionu A∨C, inkluzywnie) i przypisane także stronom A, C i przedziałom pT; ekstrapolacja per strona w v2 różni się (427929: κ_A = 0,94 ± 0,37, κ_C = 0,22 ± 0,43).")
A("- Kolumny błędów w `template.csv` to systematyki v2, nie propozycja z 2.4 (propozycja: `jump_summary.json`, pole `proposal`).")
A("- Propozycja w górę opiera się na teście (ii) w zakresie 0,09–0,14 i na założeniu braku mechanizmu właściwego stronie tagu poniżej 0,06; to założenie, nie pomiar.")
A("- Testy (i)–(iii) i warianty mieszania dzielą dawców z nominałem; błędy stosunków pochodzą ze wspólnego bootstrapu (B = 500), błędy R z Poissona obs ⊕ 1/√(tagi puli).")
A("- Test symetrii stron: błąd różnicy nachyleń liczony jak dla wielkości niezależnych (te same tagi, inne binowanie) — przybliżenie.")
A("- B = 500 zamiast 200 w v2: błędy stat. N_acc różnią się od v2 o kilka procent względnie (np. 428770 A∨C: 27,2 wobec 25,3).")
A("")
A("## 4. Pliki")
A("")
A("| plik | zawartość |")
A("|---|---|")
A("| `template.json`, `template.csv`, `template_shapes.csv` | szablon (schemat: `README.md`) |")
A("| `jump_tests.json`, `jump_summary.json` | surowe liczności testów (i)–(iii) i wielkości pochodne z propozycją |")
A("| `common.py`, `make_template.py`, `jump_tests.py`, `derive.py`, `make_report.py` | kod (selekcja = kopia logiki `run_data.py` v2; metoda importowana bez zmian) |")
A("| `make_template.log`, `jump_tests.log` | logi biegów |")
open(os.path.join(HERE, "RAPORT.md"), "w").write("\n".join(L) + "\n")
print("RAPORT.md", len(L), "wierszy")
