#!/usr/bin/env python3
"""ZMIANY.md: kazda liczba v1 -> v2 (out/final.json obu wersji). python -B zmiany.py"""
import json, os, math
H = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(H, "..", "out", "final.json"))); B = json.load(open(os.path.join(H, "out", "final.json")))
q = lambda *x: math.sqrt(sum(v * v for v in x))
NM = {"D0": "D⁰", "Dplus": "D⁺"}
def f(x, n): return f"{x:.{n}f}".replace(".", ",")
def s(x, n=1): return f"{x:+.{n}f}".replace(".", ",").replace("-", "−")
def nd(x): return 0 if x >= 10 else 1
T = ["# XSEC_final v2 — zmiany względem v1 (generowane przez `zmiany.py` z `../out/final.json` i `out/final.json`)", "",
     "Run 427929, WP70 drzewa H. PRELIMINARY, POST HOC. Pliki v1 (`../xsec_final.py`, `../out/`, `../tabele/`) nie zostały zmienione.", "",
     "## Co zmieniono w rachunku", "",
     "1. **D⁺: model fitu A2** (`FITBIAS_Dplus_2026-11-20/out/dplus_yields.json`) zamiast M0 z `XSEC_2026-11-19`: uzysk, błąd statystyczny (teraz max(1,3 × HESSE, bootstrap), jak D⁰; v1: max(HESSE, bootstrap)), wiersz „fit: model” (max(|N − N_A0|, |N − N_A4|); v1: 4 warianty tła M0), wiersz „fit: obciążenie”.",
     "   **Domknięcie — spójność z D⁰:** v1 liczyło dla D⁰ wiersz |b| ⊕ σ_b z domknięcia mb_nd, bez poprawiania uzysku; dla D⁺ v1 miało wiersz tej samej definicji z domknięcia M0 (7,8 / 12,4 / 5,6 %). v2 **zastępuje** go domknięciem A2 (9,0 / 11,1 / 7,3 %) w tej samej definicji — jeden wiersz na gatunek, bez dublowania, bez poprawki uzysku w obu gatunkach.",
     "2. **SF D⁺ 3–5 GeV: 0,949 ± 0,080** (decyzja; v1 ± 0,048). W D⁺/D⁰ (3–5): do członu SF⁺/SF⁰ = 0,981 ± 0,042 dodany kwadratowo nadmiar błędu SF⁺ (jak AUDYT_noty pkt 5).",
     f"3. **SF D⁰ z tagu D* (5–8, 8–30 GeV):** do błędu dodana kwadratowo połowa rozpiętości czterech estymatorów Δm (AUDYT_noty pkt 6 b; `AUDYT_noty_2026-11-20/out/przeliczenia.json`): {f(B['cells']['D0|5-8']['SF_spread'],3)} i {f(B['cells']['D0|8-30']['SF_spread'],3)} → "
     f"± {f(A['cells']['D0|5-8']['SF_err'],3)} → ± {f(B['cells']['D0|5-8']['SF_err'],3)}; ± {f(A['cells']['D0|8-30']['SF_err'],3)} → ± {f(B['cells']['D0|8-30']['SF_err'],3)}. Wchodzi też do D⁺/D⁰ i D*/D⁰ (kowariancja z SF_ref liczona nadal ze statystycznej części błędu).",
     f"4. **L_rel = 0,0103 / 0,4648 = {f(100*B['L_rel'],3)} %** (v1: {f(100*A['L_rel'],3)} %). Wyjaśnienie 2,24 %: `XSEC_2026-11-19/xsec.py` (`load_lumi`) składał kalibrację 2,2 % ⊕ live 0,3 % ⊕ **kompletność plików 0,3 %** (`L_COMPL_REL`); `lumi.json` podaje `błąd_wzgl` 0,0222 = kalibracja ⊕ live (0,0103 nb⁻¹), "
     "a DEFINICJA §1 — 0,4648 ± 0,0103 nb⁻¹. v2 używa zapisu DEFINICJI; wiersz kompletności plików (0,3 %) wypada z budżetu (różnica 0,02 punktu procentowego).",
     "", "D⁰: uzyski, C, BR, wartości SF, σ centralne — bez zmian. σ(D*±) bez zmian (XSEC_atlas); zmienia się tylko jej błąd świetlności.", "",
     "## Przekroje [µb] (σ ± stat ± syst ± L)", "", "| | bin | v1 | **v2** | Δσ | razem v1 → v2 [%] |", "|---|---|---|---|---|---|"]
for k in A["cells"]:
    a, b = A["cells"][k], B["cells"][k]; sp, bn = k.split("|"); n = nd(b["stat_ub"])
    g = lambda c: f"{f(c['sigma_ub'], n)} ± {f(c['stat_ub'], n)} ± {f(c['syst_ub'], n)} ± {f(c['lumi_ub'], n)}"
    T.append(f"| {NM[sp]} | {bn.replace('-', '–')} | {g(a)} | **{g(b)}** | {s(100 * (b['sigma_ub'] / a['sigma_ub'] - 1))} % | {f(100 * a['tot_rel'], 1)} → {f(100 * b['tot_rel'], 1)} |")
for bn in ("5-8", "8-30"):
    a, b = A["dstar"][bn], B["dstar"][bn]
    g = lambda c: f"{f(c['sigma_ub'], 1)} ± {f(c['stat_ub'], 1)} ± {f(c['syst_ub'], 1)} ± {f(c['lumi_ub'], 2)}"
    T.append(f"| D*± | {bn.replace('-', '–')} | {g(a)} | **{g(b)}** | 0 | — |")
T += ["", "## Wejścia i budżet, które się zmieniły", "", "| wielkość | komórka | v1 | **v2** |", "|---|---|---|---|"]
for k in A["cells"]:
    a, b = A["cells"][k], B["cells"][k]; sp, bn = k.split("|"); lab = f"{NM[sp]} {bn.replace('-', '–')}"
    if abs(a["N"] - b["N"]) > 0.5 or abs(a["rel"]["stat"] - b["rel"]["stat"]) > 1e-4:
        T.append(f"| uzysk N ± stat | {lab} | {a['N']:.0f} ± {a['N'] * a['rel']['stat']:.0f} | **{b['N']:.0f} ± {b['N'] * b['rel']['stat']:.0f}** |")
    if abs(a["SF_err"] - b["SF_err"]) > 1e-4:
        T.append(f"| SF | {lab} | {f(a['SF'], 3)} ± {f(a['SF_err'], 3)} | **{f(b['SF'], 3)} ± {f(b['SF_err'], 3)}** |")
    for key, name in (("stat", "stat [%]"), ("fit_model", "fit: model [%]"), ("fit_bias", "fit: obciążenie (domknięcie) [%]"), ("SF", "SF [%]")):
        if abs(a["rel"][key] - b["rel"][key]) > 5e-5:
            T.append(f"| {name} | {lab} | {f(100 * a['rel'][key], 1)} | **{f(100 * b['rel'][key], 1)}** |")
    if abs(a["syst_rel"] - b["syst_rel"]) > 5e-5:
        T.append(f"| razem syst [%] | {lab} | {f(100 * a['syst_rel'], 1)} | **{f(100 * b['syst_rel'], 1)}** |")
    if abs(a["closure_b"] - b["closure_b"]) > 1e-4:
        T.append(f"| domknięcie b [%] | {lab} | {s(100 * a['closure_b'])} ± {f(100 * a['closure_b_err'], 1)} | **{s(100 * b['closure_b'])} ± {f(100 * b['closure_b_err'], 1)}** |")
T.append(f"| świetlność [%] | wszystkie | {f(100 * A['L_rel'], 2)} | **{f(100 * B['L_rel'], 2)}** |")
T += ["", "## Stosunki", "", "| | bin | v1 | **v2** | razem v1 → v2 [%] | odchylenia v1 → v2 [σ] |", "|---|---|---|---|---|---|"]
for bn in ("3-5", "5-8", "8-30"):
    a, b = A["ratio"][bn], B["ratio"][bn]
    z = "; ".join(f"{k} {s(a['cmp'][k]['Z'])} → {s(b['cmp'][k]['Z'])}" for k in a["cmp"])
    T.append(f"| D⁺/D⁰ | {bn.replace('-', '–')} | {f(a['R'], 3)} ± {f(a['stat'], 3)} ± {f(a['syst'], 3)} | **{f(b['R'], 3)} ± {f(b['stat'], 3)} ± {f(b['syst'], 3)}** | {f(100 * a['tot'] / a['R'], 1)} → {f(100 * b['tot'] / b['R'], 1)} | {z} |")
for bn in ("5-8", "8-30"):
    a, b = A["dstar"][bn], B["dstar"][bn]
    z = "; ".join(f"{k} {s(a['cmp'][k]['Z'])} → {s(b['cmp'][k]['Z'])}" for k in a["cmp"])
    T.append(f"| D*/D⁰ | {bn.replace('-', '–')} | {f(a['R'], 3)} ± {f(a['stat'], 3)} ± {f(a['syst'], 3)} | **{f(b['R'], 3)} ± {f(b['stat'], 3)} ± {f(b['syst'], 3)}** | {f(100 * a['tot'] / a['R'], 1)} → {f(100 * b['tot'] / b['R'], 1)} | {z} |")
T += ["", "Budżet D⁺/D⁰ [%] (v1 → v2), wiersze zmienione:", "", "| bin | stat | fit model | fit obciążenie | SF⁺/SF⁰ |", "|---|---|---|---|---|"]
for bn in ("3-5", "5-8", "8-30"):
    a, b = A["ratio"][bn]["rel"], B["ratio"][bn]["rel"]
    T.append(f"| {bn.replace('-', '–')} | " + " | ".join(f"{f(100 * a[k], 1)} → {f(100 * b[k], 1)}" for k in ("stat", "fit_model", "fit_bias", "SF")) + " |")
T += ["", "## Porównania (stosunek ± błąd (Z))", "", "| | bin | odniesienie | v1 | **v2** |", "|---|---|---|---|---|"]
def r(o, zk="Z"): return f"{f(o['ratio'], 2)} ± {f(o['ratio_err'], 2)} ({s(o[zk])}σ)"
for k in A["cmp"]:
    sp, bn = k.split("|"); lab = f"| {NM[sp]} | {bn.replace('-', '–')} |"
    a, b = A["cmp"][k], B["cmp"][k]
    T.append(f"{lab} Pythia A3 | {r(a['pythia'])} | **{r(b['pythia'])}** |")
    fa, fb = a["fonll"], b["fonll"]
    T.append(f"{lab} FONLL prompt (centr.; Z do krawędzi pasma) | {f(fa['ratio_centr'], 2)} ± {f(fa['ratio_err'], 2)} ({s(fa['Z_band'])}σ) | **{f(fb['ratio_centr'], 2)} ± {f(fb['ratio_err'], 2)} ({s(fb['Z_band'])}σ)** |")
    if "cms" in a: T.append(f"{lab} CMS 13 TeV | {r(a['cms'])} | **{r(b['cms'])}** |")
    if "atlas7" in a: T.append(f"{lab} ATLAS 7 TeV × K_FONLL | {r(a['atlas7'])} | **{r(b['atlas7'])}** |")
    T.append(f"{lab} ALICE 13 TeV (POST HOC) | {r(a['alice'])} | **{r(b['alice'])}** |")
for bn in ("5-8", "8-30"):
    a, b = A["dstar"][bn], B["dstar"][bn]; lab = f"| D*± | {bn.replace('-', '–')} |"
    for key, name in (("pythia", "Pythia A3"), ("cms", "CMS 13 TeV"), ("atlas7", "ATLAS 7 TeV × K_FONLL")):
        T.append(f"{lab} {name} | {r(a[key])} | **{r(b[key])}** |")
T += ["", "Dane prompt / w definicji odniesienia [µb], v1 → v2:", "", "| | bin | prompt (do FONLL) | w def. CMS | w def. ALICE dσ/dy |", "|---|---|---|---|---|"]
for k in A["cmp"]:
    sp, bn = k.split("|"); a, b = A["cmp"][k], B["cmp"][k]; n = 1
    cm = f"{f(a['cms']['data_conv'], n)} ± {f(a['cms']['data_conv_err'], n)} → {f(b['cms']['data_conv'], n)} ± {f(b['cms']['data_conv_err'], n)}" if "cms" in a else "—"
    T.append(f"| {NM[sp]} | {bn.replace('-', '–')} | {f(a['fonll']['prompt_data'], n)} ± {f(a['fonll']['err'], n)} → {f(b['fonll']['prompt_data'], n)} ± {f(b['fonll']['err'], n)} | {cm} | "
             f"{f(a['alice']['data_conv'], 1)} ± {f(a['alice']['data_conv_err'], 1)} → {f(b['alice']['data_conv'], 1)} ± {f(b['alice']['data_conv_err'], 1)} |")
T += ["", "## Uwagi", "",
      "- D⁺ 8–30 GeV wobec Pythii: Z rośnie z −3,8 do " + s(B['cmp']['Dplus|8-30']['pythia']['Z']) + " (uzysk −5,2 %; błąd tylko danych, bez niepewności generatora).",
      "- D⁺ 8–30 GeV wobec FONLL: dane prompt wchodzą w pasmo (v1: +0,3σ nad górną krawędzią).",
      "- D⁰ 5–8 GeV wobec CMS: " + s(A['cmp']['D0|5-8']['cms']['Z']) + " → " + s(B['cmp']['D0|5-8']['cms']['Z']) + "σ (tylko większy błąd SF).",
      "- Odniesienia (Pythia, FONLL, CMS, ATLAS 7 TeV, ALICE, LHCb, e⁺e⁻) i czynniki przeliczeń — bez zmian.",
      "- Kryterium domknięcia fitu D⁺ nie jest spełnione w 5–8 GeV przez żadne z 9 ustawień (b = +9,2 ± 6,2 % dla A2); wybór A2 wg rankingu, post hoc."]
open(os.path.join(H, "ZMIANY.md"), "w").write("\n".join(T) + "\n"); print("\n".join(T))
