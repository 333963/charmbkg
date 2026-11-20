# Szablon tła przypadkowego tagu AFP (AFP-MIX v2) i test skoku częstości poniżej progu ξ_cal — WSTĘPNE

Wygenerowano 2026-11-17 06:19 przez `TEMPLATE/make_report.py` z `template.json` i `jump_tests.json`. Projekt B (ATLAS minbias 2022, runy 427929 i 428770, h5 v2.0, kandydaci bez score).
Metoda bez zmian: `charmbkg/afp_mix.py` v2 (estymator A), systematyki z `runs/lxplus/out_v2/result_v2.json`. Kod metody i wyniki `out/`, `out_v2/` nietknięte.

**Zakres ważności:** tylko zdarzenia z trygerów nie-AFP (GRL ∧ ≥ 1 zapisany łańcuch bez L1AFP ∧ nPileupVertices == 0); w 428770 to 4,4 % wszystkich tagów. Kandydaci po luźnej selekcji v2.0 (pT ≥ 3 GeV).
**Czego tu nie ma:** uzysków sygnału D w zdarzeniach z tagiem, fitów masy, przekrojów, histogramów masy próbki tagowanej. N_tag to liczba kandydatów (sygnał + kombinatoryka) w oknie masy.

## 1. Szablon (ZADANIE 1)

Pliki: `template.json` (wszystko), `template.csv` (liczności i f_acc, wiersz = run × strona × gatunek × pT), `template_shapes.csv` (kształty masy, wiersz = bin masy). Schemat pól: `README.md`.

- N_tag: kandydaci w oknie masy (D⁺ 1,84–1,90; D⁰ 1,835–1,895 GeV; D* Δm 0,1440–0,1470 GeV) w zdarzeniach z tagiem w regionie po cięciu ξ_cal ≤ ξ_AFP + 0,02.
- N_acc = Σ p_acc po kandydatach w oknie (wszyscy biorcy, mieszanie nominalne v2); f_acc = N_acc / N_tag.
- Błędy: stat = bootstrap otagowanych dawców (B = 500) ⊕ Poisson N_tag; pula, ekstrapolacja, ½Δ metod = wspólne dla runu błędy względne N_acc z v2 (tabela niżej), stosowane do każdego wiersza i każdej strony.
- Przedziału pT 2–3 GeV **nie ma w danych** (slim: pT ≥ 3,0; D* ≥ 3,11 GeV). Przedziały: 3–5, 5–8, 8–30 oraz 3–30 (suma). Kandydaci z pT ≥ 30 GeV pominięci (0,005 % i 0,12 %).

Systematyki względne runu (z v2) i kontrola zgodności z `out_v2` (zdarzenia, A∨C):

| run | N_tag (zdarzenia) | N_acc | stat (B = 500) | pula | ekstrapolacja | ½Δ metod | syst razem |
|---|---|---|---|---|---|---|---|
| 427929 | 18871 | 164.7 ± 20.6 | 12.5 % | 2.9 % | 47.1 % | 9.9 % | 48.2 % |
| 428770 | 16065 | 766.1 ± 27.2 | 3.5 % | 4.8 % | 10.4 % | 25.8 % | 28.2 % |

N_tag i N_acc zgadzają się z `out_v2` co do zdarzenia (kontrola w kodzie); wiersz kontrolny pT 3–150 w `template.csv` odtwarza tabele gatunków z `RAPORT_DANE_v2.md`.

### Run 427929: N_tag, N_acc, f_acc

f_acc ± stat ± pula ± ekstrapolacja ± ½Δ metod (błędy bezwzględne).

| gatunek | pT [GeV] | strona | N_tag | N_acc ± stat | f_acc | ± stat | ± pula | ± ekstrap. | ± ½Δ | ± razem |
|---|---|---|---|---|---|---|---|---|---|---|
| Dplus | 3-5 | A | 1923 | 21.1 ± 2.4 | 0.0110 | 0.0013 | 0.0003 | 0.0052 | 0.0011 | 0.0054 |
| Dplus | 3-5 | C | 1795 | 16.9 ± 2.2 | 0.0094 | 0.0012 | 0.0003 | 0.0044 | 0.0009 | 0.0047 |
| Dplus | 3-5 | A∨C | 3705 | 37.8 ± 3.2 | 0.0102 | 0.0009 | 0.0003 | 0.0048 | 0.0010 | 0.0050 |
| Dplus | 5-8 | A | 144 | 1.73 ± 0.18 | 0.0120 | 0.0017 | 0.0003 | 0.0057 | 0.0012 | 0.0060 |
| Dplus | 5-8 | C | 166 | 1.32 ± 0.15 | 0.0080 | 0.0011 | 0.0002 | 0.0038 | 0.0008 | 0.0040 |
| Dplus | 5-8 | A∨C | 309 | 3.04 ± 0.24 | 0.0098 | 0.0010 | 0.0003 | 0.0046 | 0.0010 | 0.0048 |
| Dplus | 8-30 | A | 25 | 0.26 ± 0.03 | 0.0104 | 0.0028 | 0.0003 | 0.0049 | 0.0010 | 0.0057 |
| Dplus | 8-30 | C | 37 | 0.23 ± 0.04 | 0.0062 | 0.0016 | 0.0002 | 0.0029 | 0.0006 | 0.0034 |
| Dplus | 8-30 | A∨C | 62 | 0.49 ± 0.05 | 0.0078 | 0.0014 | 0.0002 | 0.0037 | 0.0008 | 0.0040 |
| Dplus | 3-30 | A | 2092 | 23.1 ± 2.6 | 0.0110 | 0.0013 | 0.0003 | 0.0052 | 0.0011 | 0.0055 |
| Dplus | 3-30 | C | 1998 | 18.4 ± 2.3 | 0.0092 | 0.0012 | 0.0003 | 0.0043 | 0.0009 | 0.0046 |
| Dplus | 3-30 | A∨C | 4076 | 41.3 ± 3.5 | 0.0101 | 0.0009 | 0.0003 | 0.0048 | 0.0010 | 0.0050 |
| D0 | 3-5 | A | 1630 | 16.3 ± 2.5 | 0.0100 | 0.0016 | 0.0003 | 0.0047 | 0.0010 | 0.0051 |
| D0 | 3-5 | C | 1597 | 13.1 ± 2.2 | 0.0082 | 0.0014 | 0.0002 | 0.0039 | 0.0008 | 0.0042 |
| D0 | 3-5 | A∨C | 3219 | 29.2 ± 3.3 | 0.0091 | 0.0010 | 0.0003 | 0.0043 | 0.0009 | 0.0045 |
| D0 | 5-8 | A | 156 | 1.65 ± 0.22 | 0.0106 | 0.0017 | 0.0003 | 0.0050 | 0.0010 | 0.0054 |
| D0 | 5-8 | C | 169 | 1.31 ± 0.19 | 0.0078 | 0.0013 | 0.0002 | 0.0037 | 0.0008 | 0.0040 |
| D0 | 5-8 | A∨C | 324 | 2.95 ± 0.29 | 0.0091 | 0.0010 | 0.0003 | 0.0043 | 0.0009 | 0.0045 |
| D0 | 8-30 | A | 25 | 0.22 ± 0.03 | 0.0087 | 0.0023 | 0.0002 | 0.0041 | 0.0009 | 0.0047 |
| D0 | 8-30 | C | 26 | 0.18 ± 0.03 | 0.0071 | 0.0017 | 0.0002 | 0.0034 | 0.0007 | 0.0038 |
| D0 | 8-30 | A∨C | 49 | 0.40 ± 0.04 | 0.0082 | 0.0015 | 0.0002 | 0.0039 | 0.0008 | 0.0042 |
| D0 | 3-30 | A | 1811 | 18.2 ± 2.7 | 0.0100 | 0.0015 | 0.0003 | 0.0047 | 0.0010 | 0.0051 |
| D0 | 3-30 | C | 1792 | 14.6 ± 2.4 | 0.0081 | 0.0013 | 0.0002 | 0.0038 | 0.0008 | 0.0042 |
| D0 | 3-30 | A∨C | 3592 | 32.6 ± 3.6 | 0.0091 | 0.0010 | 0.0003 | 0.0043 | 0.0009 | 0.0045 |
| Dstar | 3-5 | A | 121 | 1.17 ± 0.14 | 0.0097 | 0.0015 | 0.0003 | 0.0046 | 0.0010 | 0.0049 |
| Dstar | 3-5 | C | 102 | 1.00 ± 0.14 | 0.0098 | 0.0017 | 0.0003 | 0.0046 | 0.0010 | 0.0050 |
| Dstar | 3-5 | A∨C | 223 | 2.16 ± 0.19 | 0.0097 | 0.0011 | 0.0003 | 0.0046 | 0.0010 | 0.0048 |
| Dstar | 5-8 | A | 15 | 0.17 ± 0.03 | 0.0114 | 0.0034 | 0.0003 | 0.0054 | 0.0011 | 0.0065 |
| Dstar | 5-8 | C | 15 | 0.11 ± 0.02 | 0.0076 | 0.0023 | 0.0002 | 0.0036 | 0.0007 | 0.0043 |
| Dstar | 5-8 | A∨C | 30 | 0.28 ± 0.03 | 0.0094 | 0.0020 | 0.0003 | 0.0044 | 0.0009 | 0.0050 |
| Dstar | 8-30 | A | 4 | 0.02 ± 0.00 | 0.0042 | 0.0024 | 0.0001 | 0.0020 | 0.0004 | 0.0032 |
| Dstar | 8-30 | C | 0 | 0.01 ± 0.00 | — | — | — | — | — | — |
| Dstar | 8-30 | A∨C | 4 | 0.03 ± 0.01 | 0.0076 | 0.0041 | 0.0002 | 0.0036 | 0.0008 | 0.0055 |
| Dstar | 3-30 | A | 140 | 1.36 ± 0.16 | 0.0097 | 0.0014 | 0.0003 | 0.0046 | 0.0010 | 0.0049 |
| Dstar | 3-30 | C | 117 | 1.13 ± 0.15 | 0.0096 | 0.0016 | 0.0003 | 0.0045 | 0.0009 | 0.0049 |
| Dstar | 3-30 | A∨C | 257 | 2.48 ± 0.22 | 0.0096 | 0.0011 | 0.0003 | 0.0045 | 0.0010 | 0.0048 |

### Run 428770: N_tag, N_acc, f_acc

f_acc ± stat ± pula ± ekstrapolacja ± ½Δ metod (błędy bezwzględne).

| gatunek | pT [GeV] | strona | N_tag | N_acc ± stat | f_acc | ± stat | ± pula | ± ekstrap. | ± ½Δ | ± razem |
|---|---|---|---|---|---|---|---|---|---|---|
| Dplus | 3-5 | A | 3100 | 157.5 ± 5.3 | 0.0508 | 0.0022 | 0.0024 | 0.0053 | 0.0131 | 0.0145 |
| Dplus | 3-5 | C | 2793 | 155.7 ± 5.2 | 0.0557 | 0.0025 | 0.0027 | 0.0058 | 0.0144 | 0.0159 |
| Dplus | 3-5 | A∨C | 5863 | 312.2 ± 7.5 | 0.0533 | 0.0016 | 0.0026 | 0.0055 | 0.0137 | 0.0151 |
| Dplus | 5-8 | A | 939 | 42.2 ± 1.5 | 0.0449 | 0.0024 | 0.0022 | 0.0047 | 0.0116 | 0.0129 |
| Dplus | 5-8 | C | 858 | 42.4 ± 1.5 | 0.0494 | 0.0028 | 0.0024 | 0.0051 | 0.0127 | 0.0142 |
| Dplus | 5-8 | A∨C | 1789 | 84.2 ± 2.1 | 0.0471 | 0.0018 | 0.0023 | 0.0049 | 0.0121 | 0.0134 |
| Dplus | 8-30 | A | 569 | 24.4 ± 0.9 | 0.0428 | 0.0028 | 0.0021 | 0.0044 | 0.0110 | 0.0124 |
| Dplus | 8-30 | C | 563 | 24.8 ± 0.9 | 0.0441 | 0.0029 | 0.0021 | 0.0046 | 0.0113 | 0.0128 |
| Dplus | 8-30 | A∨C | 1128 | 49.0 ± 1.3 | 0.0434 | 0.0020 | 0.0021 | 0.0045 | 0.0112 | 0.0124 |
| Dplus | 3-30 | A | 4608 | 224.1 ± 7.7 | 0.0486 | 0.0020 | 0.0023 | 0.0050 | 0.0125 | 0.0138 |
| Dplus | 3-30 | C | 4214 | 222.9 ± 7.6 | 0.0529 | 0.0023 | 0.0025 | 0.0055 | 0.0136 | 0.0151 |
| Dplus | 3-30 | A∨C | 8780 | 445.4 ± 10.8 | 0.0507 | 0.0015 | 0.0024 | 0.0053 | 0.0131 | 0.0144 |
| D0 | 3-5 | A | 2253 | 108.3 ± 4.1 | 0.0481 | 0.0022 | 0.0023 | 0.0050 | 0.0124 | 0.0137 |
| D0 | 3-5 | C | 2043 | 110.0 ± 4.1 | 0.0539 | 0.0025 | 0.0026 | 0.0056 | 0.0139 | 0.0154 |
| D0 | 3-5 | A∨C | 4281 | 217.5 ± 5.9 | 0.0508 | 0.0017 | 0.0024 | 0.0053 | 0.0131 | 0.0144 |
| D0 | 5-8 | A | 684 | 35.2 ± 1.4 | 0.0515 | 0.0030 | 0.0025 | 0.0053 | 0.0133 | 0.0148 |
| D0 | 5-8 | C | 671 | 35.1 ± 1.4 | 0.0524 | 0.0030 | 0.0025 | 0.0054 | 0.0135 | 0.0151 |
| D0 | 5-8 | A∨C | 1352 | 70.1 ± 2.0 | 0.0518 | 0.0021 | 0.0025 | 0.0054 | 0.0133 | 0.0148 |
| D0 | 8-30 | A | 437 | 18.0 ± 0.8 | 0.0411 | 0.0027 | 0.0020 | 0.0043 | 0.0106 | 0.0119 |
| D0 | 8-30 | C | 389 | 18.7 ± 0.8 | 0.0481 | 0.0033 | 0.0023 | 0.0050 | 0.0124 | 0.0139 |
| D0 | 8-30 | A∨C | 820 | 36.5 ± 1.1 | 0.0446 | 0.0021 | 0.0021 | 0.0046 | 0.0115 | 0.0127 |
| D0 | 3-30 | A | 3374 | 161.5 ± 6.3 | 0.0479 | 0.0021 | 0.0023 | 0.0050 | 0.0123 | 0.0137 |
| D0 | 3-30 | C | 3103 | 163.9 ± 6.3 | 0.0528 | 0.0024 | 0.0025 | 0.0055 | 0.0136 | 0.0151 |
| D0 | 3-30 | A∨C | 6453 | 324.1 ± 8.9 | 0.0502 | 0.0016 | 0.0024 | 0.0052 | 0.0129 | 0.0142 |
| Dstar | 3-5 | A | 209 | 9.53 ± 0.34 | 0.0456 | 0.0038 | 0.0022 | 0.0047 | 0.0117 | 0.0134 |
| Dstar | 3-5 | C | 185 | 9.46 ± 0.34 | 0.0511 | 0.0044 | 0.0025 | 0.0053 | 0.0132 | 0.0151 |
| Dstar | 3-5 | A∨C | 392 | 18.9 ± 0.5 | 0.0483 | 0.0029 | 0.0023 | 0.0050 | 0.0124 | 0.0139 |
| Dstar | 5-8 | A | 66 | 3.12 ± 0.12 | 0.0472 | 0.0064 | 0.0023 | 0.0049 | 0.0122 | 0.0148 |
| Dstar | 5-8 | C | 66 | 3.19 ± 0.12 | 0.0484 | 0.0063 | 0.0023 | 0.0050 | 0.0125 | 0.0150 |
| Dstar | 5-8 | A∨C | 131 | 6.27 ± 0.17 | 0.0478 | 0.0045 | 0.0023 | 0.0050 | 0.0123 | 0.0142 |
| Dstar | 8-30 | A | 24 | 1.14 ± 0.05 | 0.0474 | 0.0106 | 0.0023 | 0.0049 | 0.0122 | 0.0171 |
| Dstar | 8-30 | C | 24 | 1.22 ± 0.06 | 0.0507 | 0.0110 | 0.0024 | 0.0053 | 0.0131 | 0.0181 |
| Dstar | 8-30 | A∨C | 47 | 2.35 ± 0.08 | 0.0499 | 0.0079 | 0.0024 | 0.0052 | 0.0129 | 0.0161 |
| Dstar | 3-30 | A | 299 | 13.8 ± 0.5 | 0.0461 | 0.0033 | 0.0022 | 0.0048 | 0.0119 | 0.0134 |
| Dstar | 3-30 | C | 275 | 13.9 ± 0.5 | 0.0504 | 0.0037 | 0.0024 | 0.0052 | 0.0130 | 0.0147 |
| Dstar | 3-30 | A∨C | 570 | 27.5 ± 0.7 | 0.0483 | 0.0025 | 0.0023 | 0.0050 | 0.0124 | 0.0138 |

### Kształt masy składowej przypadkowej a kształt wszystkich kandydatów bez wymogu tagu

Składowa przypadkowa: histogram Σ p_acc po kandydatach (biny 10 MeV dla D⁺ 1,72–2,02 i D⁰ 1,75–2,00; 1 MeV dla Δm D* 0,139–0,160). Odniesienie: wszyscy kandydaci biorców bez wymogu tagu (zdarzenia z tagiem to < 1 % tej próbki). Test: d_i = h_acc,i − (Σh_acc)·s_i, χ² = dᵀC⁻¹d, ndf = liczba binów − 1.

Kowariancja C = C_boot + C_samp − C_overlap:
- C_boot: bootstrap otagowanych dawców przy ustalonych biorcach (B = 500). To jest korelacja wag z jednej puli dawców: jeden dawca daje wagę wszystkim biorcom swojej komórki (LB × nSelBL), więc szum wag przenosi się między komórkami o różnym kształcie masy.
- C_samp: próbkowanie zdarzeń-biorców przy zrealizowanych wagach, Σ_t (p_t − p̄)²·g_t g_tᵀ (g_t = wektor kandydatów zdarzenia w binach minus oczekiwanie).
- C_overlap = Σ_t Var_boot(p_t)·g_t g_tᵀ: część C_samp, którą C_boot już zawiera (E[(p − p̄)²] = (π − π̄)² + Var p). Bez tej poprawki wariancja jest liczona podwójnie: C_overlap/C_samp ≈ 0.95 (427929, mało dawców) i 0.68 (428770).

Warianty: χ²_cons z C_boot + C_samp (wariancja zawyżona → dolna granica istotności); χ²_samp z samym C_samp (H0 ścisłe: masa niezależna od wag), kalibrowane permutacją p_acc między biorcami: χ²/ndf = 1.04 ± 0.33 w 90 kategoriach. Kontrola χ² nominalnego: w 427929 dla D⁰ i D* (gdzie efektu nie widać) średnie χ²/ndf = 1.16 (cons: 0.51). W kategoriach z kilkoma kandydatami z tagiem (D* 8–30 GeV w 427929) przybliżenie gaussowskie χ² nie jest wiarygodne.

Kolumny: udział okna piku = (część h_acc w oknie piku)/(część próbki nietagowanej w oknie piku); nachylenie = liniowe nachylenie stosunku kształtów acc/nietagowane względem masy [1/GeV].

| run | gatunek | pT | strona | χ²/ndf | p | χ²_cons/ndf | p_cons | χ²_samp/ndf | udział okna piku | nachylenie [1/GeV] |
|---|---|---|---|---|---|---|---|---|---|---|
| 427929 | Dplus | 3-5 | A∨C | 42.6/29 | 0.049 | 23.1/29 | 0.77 | 64.3/29 | 1.001 ± 0.005 | -0.20 ± 0.04 |
| 427929 | Dplus | 5-8 | A∨C | 41.9/29 | 0.057 | 17.9/29 | 0.95 | 36.6/29 | 1.009 ± 0.021 | -0.34 ± 0.13 |
| 427929 | Dplus | 8-30 | A∨C | 37.8/29 | 0.13 | 15.3/29 | 0.98 | 29.6/29 | 1.068 ± 0.056 | -0.50 ± 0.25 |
| 427929 | Dplus | 3-30 | A | 32.6/29 | 0.29 | 17.6/29 | 0.95 | 48.6/29 | 0.999 ± 0.006 | -0.19 ± 0.05 |
| 427929 | Dplus | 3-30 | C | 43.0/29 | 0.045 | 17.8/29 | 0.95 | 54.4/29 | 1.005 ± 0.008 | -0.25 ± 0.07 |
| 427929 | Dplus | 3-30 | A∨C | 43.4/29 | 0.042 | 24.0/29 | 0.73 | 68.4/29 | 1.002 ± 0.005 | -0.21 ± 0.04 |
| 427929 | D0 | 3-5 | A∨C | 22.7/24 | 0.54 | 8.8/24 | 1 | 16.2/24 | 1.003 ± 0.007 | -0.06 ± 0.05 |
| 427929 | D0 | 5-8 | A∨C | 29.6/24 | 0.2 | 11.8/24 | 0.98 | 22.8/24 | 1.000 ± 0.025 | -0.30 ± 0.17 |
| 427929 | D0 | 8-30 | A∨C | 27.7/24 | 0.27 | 11.5/24 | 0.98 | 22.0/24 | 1.087 ± 0.063 | 0.06 ± 0.41 |
| 427929 | D0 | 3-30 | A | 21.2/24 | 0.63 | 9.4/24 | 1 | 20.8/24 | 1.003 ± 0.011 | -0.15 ± 0.08 |
| 427929 | D0 | 3-30 | C | 30.2/24 | 0.18 | 9.8/24 | 1 | 20.7/24 | 1.006 ± 0.009 | 0.01 ± 0.06 |
| 427929 | D0 | 3-30 | A∨C | 23.6/24 | 0.48 | 9.2/24 | 1 | 17.6/24 | 1.004 ± 0.007 | -0.08 ± 0.05 |
| 427929 | Dstar | 3-5 | A∨C | 13.5/20 | 0.85 | 6.1/20 | 1 | 12.3/20 | 1.011 ± 0.024 | 0.8 ± 2.1 |
| 427929 | Dstar | 5-8 | A∨C | 28.1/20 | 0.11 | 12.1/20 | 0.91 | 25.1/20 | 1.188 ± 0.080 | -1.9 ± 6.1 |
| 427929 | Dstar | 8-30 | A∨C | 21.3/20 | 0.38 | 14.1/20 | 0.83 | 21.7/20 | 1.146 ± 0.199 | -7.7 ± 15.3 |
| 427929 | Dstar | 3-30 | A | 23.8/20 | 0.25 | 12.1/20 | 0.91 | 27.3/20 | 1.009 ± 0.025 | -0.4 ± 2.7 |
| 427929 | Dstar | 3-30 | C | 18.0/20 | 0.59 | 7.4/20 | 1 | 16.2/20 | 1.056 ± 0.039 | 1.2 ± 2.8 |
| 427929 | Dstar | 3-30 | A∨C | 15.7/20 | 0.73 | 7.0/20 | 1 | 14.3/20 | 1.030 ± 0.022 | 0.5 ± 1.9 |
| 428770 | Dplus | 3-5 | A∨C | 364.4/29 | 1.2e-59 | 249.7/29 | 5.7e-37 | 500.2/29 | 1.001 ± 0.002 | -0.24 ± 0.01 |
| 428770 | Dplus | 5-8 | A∨C | 196.5/29 | 8.3e-27 | 124.0/29 | 1e-13 | 224.2/29 | 0.997 ± 0.003 | -0.28 ± 0.02 |
| 428770 | Dplus | 8-30 | A∨C | 117.9/29 | 1.1e-12 | 70.6/29 | 2.5e-05 | 120.5/29 | 1.013 ± 0.005 | -0.28 ± 0.03 |
| 428770 | Dplus | 3-30 | A | 309.8/29 | 9.4e-49 | 221.2/29 | 1.8e-31 | 564.3/29 | 1.001 ± 0.002 | -0.29 ± 0.02 |
| 428770 | Dplus | 3-30 | C | 328.1/29 | 2.2e-52 | 231.7/29 | 1.7e-33 | 597.0/29 | 1.001 ± 0.002 | -0.30 ± 0.02 |
| 428770 | Dplus | 3-30 | A∨C | 578.5/29 | 1.8e-103 | 423.9/29 | 1.1e-71 | 1050.9/29 | 1.001 ± 0.001 | -0.29 ± 0.01 |
| 428770 | D0 | 3-5 | A∨C | 49.2/24 | 0.0018 | 30.6/24 | 0.17 | 55.5/24 | 1.003 ± 0.002 | -0.10 ± 0.02 |
| 428770 | D0 | 5-8 | A∨C | 66.6/24 | 7e-06 | 37.1/24 | 0.043 | 60.3/24 | 1.003 ± 0.004 | -0.20 ± 0.03 |
| 428770 | D0 | 8-30 | A∨C | 71.5/24 | 1.3e-06 | 40.2/24 | 0.02 | 67.7/24 | 1.006 ± 0.007 | -0.31 ± 0.05 |
| 428770 | D0 | 3-30 | A | 48.1/24 | 0.0025 | 30.3/24 | 0.17 | 60.9/24 | 1.005 ± 0.003 | -0.13 ± 0.02 |
| 428770 | D0 | 3-30 | C | 31.4/24 | 0.14 | 19.2/24 | 0.74 | 39.9/24 | 1.002 ± 0.003 | -0.10 ± 0.02 |
| 428770 | D0 | 3-30 | A∨C | 69.9/24 | 2.3e-06 | 43.0/24 | 0.01 | 78.9/24 | 1.003 ± 0.002 | -0.12 ± 0.02 |
| 428770 | Dstar | 3-5 | A∨C | 24.6/20 | 0.22 | 15.1/20 | 0.77 | 24.9/20 | 1.006 ± 0.008 | -2.2 ± 0.6 |
| 428770 | Dstar | 5-8 | A∨C | 15.5/20 | 0.75 | 9.1/20 | 0.98 | 14.3/20 | 1.016 ± 0.014 | 0.1 ± 1.0 |
| 428770 | Dstar | 8-30 | A∨C | 29.9/20 | 0.072 | 17.5/20 | 0.62 | 29.7/20 | 1.063 ± 0.025 | -2.8 ± 1.8 |
| 428770 | Dstar | 3-30 | A | 30.0/20 | 0.07 | 16.2/20 | 0.71 | 25.0/20 | 1.016 ± 0.008 | -1.3 ± 0.6 |
| 428770 | Dstar | 3-30 | C | 18.6/20 | 0.55 | 10.5/20 | 0.96 | 17.4/20 | 1.011 ± 0.009 | -1.5 ± 0.6 |
| 428770 | Dstar | 3-30 | A∨C | 27.2/20 | 0.13 | 15.3/20 | 0.76 | 23.1/20 | 1.013 ± 0.006 | -1.3 ± 0.5 |

**Wynik.**
- D⁺: kształt składowej przypadkowej **różni się** od kształtu próbki nietagowanej. 428770: χ²/ndf = 579/29 (zachowawczo 424/29), nachylenie -0.290 ± 0.013 /GeV, czyli około ±4.4 % na krańcach zakresu 1,72–2,02 GeV (więcej przy małej masie); to samo w każdym przedziale pT i po obu stronach. 427929: nachylenie -0.21 ± 0.04 /GeV (5.2σ, ten sam znak i wielkość), χ²/ndf = 43.4/29 (p = 0.04; zachowawczo p = 0.73).
- D⁰: w 428770 mniejsze nachylenie tego samego znaku (-0.117 ± 0.015 /GeV; χ²/ndf = 70/24, p = 2e-06; zachowawczo p = 0.01), rosnące z pT; w 427929 zgodność.
- D* (Δm): χ² zgodne w obu runach (A∨C: najmniejsze p = 0.07); w 428770 słaba wskazówka nachylenia (-1.3 ± 0.5 /GeV, 2.9σ) i udziału okna piku 1.013 ± 0.006.
- Udział okna piku dla D⁺ i D⁰ jest zgodny z 1 (428770, 3–30 GeV: 1.0008 ± 0.0015 i 1.0034 ± 0.0020): okno leży w środku zakresu, więc nachylenie nie zmienia f_acc w oknie piku. Zmienia natomiast przewidywanie, jeśli tło przypadkowe normalizuje się w sidebandach po jednej stronie piku (do kilku procent).
- Zalecenie: jako kształt masy tła przypadkowego brać histogram ważony p_acc z `template_shapes.csv`, nie kształt próbki nietagowanej.

Hipoteza (niesprawdzana tutaj): p_acc faworyzuje zdarzenia o małym ξ_cal i małej krotności, a kształt masy kombinatoryki D⁺ (trzy tory) zależy od krotności.

## 2. Test skoku częstości przypadkowych poniżej progu ξ_cal dawców (ZADANIE 2)

Pytanie: czy częstość protonów przypadkowych u biorców z małym ξ_cal po stronie tagu (obszar użycia, ξ_cal ≲ 0,12) jest taka jak u dawców (ξ_cal > 0,14). Obszaru użycia nie da się sprawdzić wprost, bo tagi prawdziwe i przypadkowe są tam nieodróżnialne.

Profil obszaru użycia (wagi par „acc” z mieszania nominalnego; s = strona tagu, o = strona przeciwna):

| run | strona | ⟨ln ξ_cal,s⟩ | część wag z ξ_cal,s < 0,06 | część wag z ξ_cal,o > 0,14 | ⟨nSelBL⟩ użycie | ⟨nSelBL⟩ pula z D |
|---|---|---|---|---|---|---|
| 427929 | A | -3.00 | 55 % | 33 % | 20.5 | 38.3 |
| 427929 | C | -3.02 | 56 % | 29 % | 19.4 | 38.1 |
| 428770 | A | -2.96 | 54 % | 44 % | 26.8 | 47.2 |
| 428770 | C | -2.98 | 55 % | 42 % | 25.9 | 47.0 |

Wszystkie testy niżej porównują liczbę obserwowanych tagów z przewidywaniem z tych samych warstw co nominał (LB × nSelBL × „z D”), R = obs/pred. „z D” = zdarzenia z ≥ 1 kandydatem (klasa biorców); „wszystkie” = z D i bez D (każde względem własnej klasy).

### 2.1 Test (i): strona przeciwna i estymator „lustrzany”

**Co mierzy.** W zdarzeniach z ξ_cal,s > 0,14 (własny proton po stronie s wykluczony) częstość tagu po stronie s w funkcji ξ_cal,o strony przeciwnej, do najmniejszych ξ_cal,o. To zależność częstości przypadkowej od typu zdarzenia (jedna strona z małym ξ_cal) przy ustalonych LB i nSelBL.
**Czego nie mierzy.** Efektów związanych ze stroną protonu przy małym ξ_cal,s. Przeniesienie na obszar użycia wymaga założenia, że częstość nie zależy od tego, po której stronie jest małe ξ_cal („ślepota na stronę”); ponadto tylko 29–44 % wag obszaru użycia ma ξ_cal,o > 0,14 (dokładne odbicie lustrzane zdarzeń testowych), reszta ma małe ξ_cal po obu stronach.

Run 427929, zdarzenia z D, strony A + C: R w przedziałach ξ_cal,o

| ξ_cal,o | 0–0.02 | 0.02–0.04 | 0.04–0.06 | 0.06–0.08 | 0.08–0.1 | 0.1–0.12 | 0.12–0.14 | 0.14–0.2 | 0.2–0.3 | 0.3–0.5 | 0.5–∞ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| obs | 1 | 9 | 16 | 21 | 19 | 30 | 43 | 128 | 129 | 41 | 1 |
| R | 0.37 ± 0.61 | 1.01 ± 0.34 | 1.01 ± 0.25 | 0.90 ± 0.21 | 0.61 ± 0.18 | 0.83 ± 0.17 | 1.08 ± 0.16 | 1.08 ± 0.10 | 1.06 ± 0.09 | 1.06 ± 0.17 | 1.63 ± 1.63 |

Run 428770, zdarzenia z D, strony A + C: R w przedziałach ξ_cal,o

| ξ_cal,o | 0–0.02 | 0.02–0.04 | 0.04–0.06 | 0.06–0.08 | 0.08–0.1 | 0.1–0.12 | 0.12–0.14 | 0.14–0.2 | 0.2–0.3 | 0.3–0.5 | 0.5–∞ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| obs | 13 | 47 | 74 | 132 | 159 | 201 | 242 | 762 | 1058 | 595 | 32 |
| R | 0.78 ± 0.24 | 0.81 ± 0.13 | 0.73 ± 0.10 | 0.91 ± 0.08 | 0.86 ± 0.07 | 0.91 ± 0.07 | 0.98 ± 0.06 | 0.97 ± 0.04 | 1.06 ± 0.03 | 1.12 ± 0.05 | 1.28 ± 0.23 |

| run | próbka | R(ξ_cal,o < 0,14) | R(ξ_cal,o > 0,14) | stosunek mały/duży |
|---|---|---|---|---|
| 427929 | z D | 139 / 158.3 = 0.878 ± 0.079 | 299 / 279.2 = 1.071 ± 0.062 | 0.820 ± 0.084 |
| 427929 | z D, bez protonu po stronie o | 134 / 157.6 = 0.850 ± 0.080 | 293 / 279.2 = 1.049 ± 0.061 | 0.810 ± 0.084 |
| 427929 | wszystkie | 455 / 484.2 = 0.940 ± 0.045 | 447 / 417.4 = 1.071 ± 0.051 | 0.877 ± 0.058 |
| 428770 | z D | 868 / 973.6 = 0.892 ± 0.032 | 2447 / 2341.6 = 1.045 ± 0.021 | 0.853 ± 0.034 |
| 428770 | z D, bez protonu po stronie o | 849 / 968.4 = 0.877 ± 0.032 | 2412 / 2339.9 = 1.031 ± 0.021 | 0.850 ± 0.034 |
| 428770 | wszystkie | 1252 / 1514.3 = 0.827 ± 0.026 | 2928 / 2665.5 = 1.098 ± 0.020 | 0.753 ± 0.025 |

Estymator lustrzany: nominał z dodatkową klasą warstw — dawca (ξ_cal,s > 0,14) ma ξ_cal,o w tym samym przedziale (0; 0,02; …; 0,14) co ξ_cal,s biorcy. Błąd stosunku ze wspólnego bootstrapu dawców. „bez p_o” = dawcy bez protonu po stronie przeciwnej (u biorców proton po stronie o ma < 0,12 % wag).

| run | N_acc nominał | nominał, dawcy bez p_o / nominał | lustrzany / nominał | **lustrzany bez p_o / nominał bez p_o (κ_L)** | strona A | strona C | górna granica κ_L (95 % CL) |
|---|---|---|---|---|---|---|---|
| 427929 | 164.7 ± 20.6 | 0.967 ± 0.028 | 0.881 ± 0.126 | **0.825 ± 0.124** | 0.78 ± 0.17 | 0.89 ± 0.17 | 1.03 |
| 428770 | 766.1 ± 27.2 | 0.982 ± 0.005 | 0.798 ± 0.042 | **0.790 ± 0.044** | 0.81 ± 0.06 | 0.77 ± 0.06 | 0.86 |

**Wynik.** Częstość tagu w zdarzeniach z małym ξ_cal po stronie przeciwnej jest **mniejsza**, nie większa: stosunek mały/duży = 0.85 ± 0.03 (428770) i 0.82 ± 0.08 (427929); κ_L = 0.79 ± 0.04 i 0.83 ± 0.12. Skoku w górę o mechanizmie ślepym na stronę dane nie dopuszczają (górne granice κ_L w tabeli; scenariusz toy ×2 dałby R ≈ 2 przy małym ξ_cal,o).

Test symetrii stron (czy zależność jest ślepa na stronę): zdarzenia puli z ξ_cal > 0,14 po obu stronach, nachylenie dR/d ln ξ_cal względem strony tagu i strony przeciwnej (przedziały 0,14–0,2–0,3–0,5; te same zdarzenia, błąd różnicy w przybliżeniu niezależnych nachyleń).

| run | próbka | nachylenie, strona tagu | nachylenie, strona przeciwna | różnica | z |
|---|---|---|---|---|---|
| 427929 | z D | 0.07 ± 0.25 | -0.03 ± 0.25 | -0.10 ± 0.35 | -0.3 |
| 427929 | wszystkie | 0.18 ± 0.21 | -0.03 ± 0.21 | -0.21 ± 0.30 | -0.7 |
| 428770 | z D | -0.05 ± 0.08 | 0.22 ± 0.08 | 0.26 ± 0.11 | 2.4 |
| 428770 | wszystkie | -0.21 ± 0.07 | 0.13 ± 0.07 | 0.35 ± 0.10 | 3.3 |

W 428770 częstość rośnie z ξ_cal strony przeciwnej (0.22 ± 0.08), a nie zależy od ξ_cal strony tagu (-0.05 ± 0.08); różnica 2.4σ. To wskazówka (nie dowód), że efekt jest związany ze stroną przeciwną do protonu, a więc że założenie ślepoty na stronę nie jest spełnione i κ_L nie przenosi się wprost na obszar użycia. W 427929 test nie ma mocy.

### 2.2 Test (ii): protony kinematycznie niemożliwe

**Co mierzy.** Dokładnie jeden proton po stronie s z ξ_p w wycinku [a, b] (0,02–0,03 poniżej okna tagu oraz 0,03–0,04, 0,04–0,05, 0,05–0,06, 0,06–0,08, 0,08–0,10) i ξ_cal,s > b + 0,02 + m. Taki proton nie pochodzi z tego oddziaływania, jeśli m przewyższa rozdzielczość ξ_cal. Klasy marginesu m: bliski 0,010–0,025, średni 0,025–0,040, daleki ≥ 0,040. Zakres ξ_cal,s: 0,06–0,14.
**Czego nie mierzy.** ξ_cal,s < 0,06 (najmniejsze ξ_p w danych to ≈ 0,02; 55 % wag obszaru użycia leży poniżej 0,06) oraz, przy dalekim marginesie, ξ_cal,s < 0,08. Margines bliski jest zanieczyszczony własnymi protonami (rozdzielczość ξ_cal).

| run | próbka | bliski | średni | daleki | bliski − daleki |
|---|---|---|---|---|---|
| 427929 | z D | 55 / 59.4 = 0.93 ± 0.14 | 50 / 55.7 = 0.90 ± 0.14 | 122 / 107.1 = 1.14 ± 0.12 | -0.21 ± 0.17 |
| 427929 | wszystkie | 375 / 391.5 = 0.96 ± 0.06 | 280 / 306.3 = 0.91 ± 0.07 | 496 / 457.6 = 1.08 ± 0.07 | -0.13 ± 0.07 |
| 428770 | z D | 345 / 299.0 = 1.15 ± 0.06 | 290 / 280.3 = 1.03 ± 0.06 | 504 / 505.0 = 1.00 ± 0.05 | 0.16 ± 0.08 |
| 428770 | wszystkie | 582 / 490.9 = 1.19 ± 0.05 | 499 / 437.7 = 1.14 ± 0.05 | 798 / 735.5 = 1.09 ± 0.04 | 0.10 ± 0.06 |

R w przedziałach ξ_cal,s (zdarzenia z D; obs i R ± błąd):

| 427929: ξ_cal,s | bliski | średni | daleki |
|---|---|---|---|
| 0.06–0.07 | 8: 1.54 ± 0.55 | — | — |
| 0.07–0.08 | 10: 0.96 ± 0.32 | 3: 0.96 ± 0.58 | — |
| 0.08–0.09 | 11: 1.13 ± 0.35 | 11: 1.00 ± 0.31 | — |
| 0.09–0.10 | 8: 1.02 ± 0.37 | 9: 0.74 ± 0.30 | 9: 1.28 ± 0.44 |
| 0.10–0.12 | 8: 0.60 ± 0.29 | 12: 0.84 ± 0.28 | 42: 1.03 ± 0.17 |
| 0.12–0.14 | 10: 0.78 ± 0.29 | 15: 1.00 ± 0.28 | 71: 1.19 ± 0.15 |

| 428770: ξ_cal,s | bliski | średni | daleki |
|---|---|---|---|
| 0.06–0.07 | 20: 0.93 ± 0.22 | — | — |
| 0.07–0.08 | 57: 1.23 ± 0.17 | 14: 1.14 ± 0.31 | — |
| 0.08–0.09 | 64: 1.32 ± 0.17 | 49: 1.07 ± 0.16 | — |
| 0.09–0.10 | 49: 1.11 ± 0.16 | 66: 1.16 ± 0.15 | 23: 0.83 ± 0.19 |
| 0.10–0.12 | 83: 1.19 ± 0.13 | 93: 1.16 ± 0.12 | 190: 1.04 ± 0.08 |
| 0.12–0.14 | 72: 1.05 ± 0.13 | 68: 0.80 ± 0.11 | 291: 0.99 ± 0.06 |

R w funkcji ξ_cal strony przeciwnej (zdarzenia z D):

| run | ξ_cal,o | bliski | średni | daleki |
|---|---|---|---|---|
| 427929 | < 0,06 | 10: 1.26 ± 0.40 | 9: 1.34 ± 0.45 | 14: 1.26 ± 0.34 |
| 427929 | 0,06–0,14 | 25: 0.99 ± 0.20 | 11: 0.48 ± 0.21 | 53: 1.26 ± 0.18 |
| 427929 | > 0,14 | 20: 0.76 ± 0.20 | 30: 1.15 ± 0.22 | 55: 1.02 ± 0.15 |
| 428770 | < 0,06 | 22: 0.61 ± 0.17 | 27: 0.86 ± 0.18 | 48: 0.98 ± 0.14 |
| 428770 | 0,06–0,14 | 110: 1.04 ± 0.10 | 85: 0.87 ± 0.10 | 162: 0.97 ± 0.08 |
| 428770 | > 0,14 | 213: 1.36 ± 0.09 | 178: 1.17 ± 0.09 | 294: 1.02 ± 0.06 |

**Wynik.** Przy dalekim marginesie (ξ_cal,s 0,09–0,14) R = 1.00 ± 0.05 (428770) i 1.14 ± 0.12 (427929): brak skoku w tym zakresie; górne granice 95 % CL: 1.08 i 1.33.
Nadwyżka 428770 tuż pod progiem (przegląd: 1,14 ± 0,06; B(ii): 1,12 ± 0,06) siedzi w marginesie bliskim (1.15 ± 0.06) i tam rośnie z ξ_cal strony przeciwnej (0.61 ± 0.17 → 1.36 ± 0.09), a przy dalekim marginesie jej nie ma przy tym samym ξ_cal,s. Skok częstości przypadkowych dawałby to samo R dla każdego ξ_p; obraz pasuje do przecieku własnych protonów przez rozdzielczość ξ_cal (prawdziwa dyfrakcja po stronie s ma układ X po stronie o, czyli duże ξ_cal,o). To interpretacja: różnica bliski − daleki ma 2.0σ.

### 2.3 Test (iii): zależność od nSelBL, μ, LB

**Co mierzy.** (a) Czy niedobór z testu (i) jest skupiony w jakiejś klasie nSelBL, przedziale μ lub części runu. (b) Czy warstwy nSelBL nominału są dość drobne (użycie ma ⟨nSelBL⟩ ≈ 20–27, pula 38–47).
**Czego nie mierzy.** Niczego poniżej progu po stronie tagu.

Stosunek R(ξ_cal,o < 0,12)/R(ξ_cal,o ≥ 0,12), zdarzenia z D, strony razem:

- 427929, nSelBL: 10–15: 1.74 ± 1.88; 15–20: 0.69 ± 0.28; 20–25: 0.72 ± 0.22; 25–30: 0.67 ± 0.18; 30–40: 0.70 ± 0.17; 40–50: 0.82 ± 0.25
- 427929, μ (tercyle): 0.001–0.004: 0.85 ± 0.20; 0.004–0.008: 0.66 ± 0.15; ≥ 0.008: 0.76 ± 0.13
- 427929, LB (tercyle): 119–193: 0.84 ± 0.16; 193–258: 0.67 ± 0.14; ≥ 258: 0.77 ± 0.16
- 428770, nSelBL: 10–15: 0.33 ± 0.10; 15–20: 0.64 ± 0.12; 20–25: 0.87 ± 0.13; 25–30: 0.84 ± 0.11; 30–40: 0.78 ± 0.07; 40–50: 0.81 ± 0.09; 50–80: 0.95 ± 0.11
- 428770, μ (tercyle): 0.010–0.045: 0.76 ± 0.07; 0.045–0.058: 0.89 ± 0.07; ≥ 0.058: 0.82 ± 0.06
- 428770, LB (tercyle): 4–232: 0.87 ± 0.07; 232–278: 0.84 ± 0.07; ≥ 278: 0.78 ± 0.06

| run | drobne warstwy nSelBL (19 brzegów) / nominał | dawcy bez protonu po stronie o / nominał | v2: warstwy μ zamiast LB |
|---|---|---|---|
| 427929 | 1.030 ± 0.047 | 0.967 ± 0.028 | 1.035 ± 0.017 |
| 428770 | 1.011 ± 0.007 | 0.982 ± 0.005 | 1.007 ± 0.003 |

**Wynik.** Niedobór z testu (i) jest podobny w tercylach μ i LB. W klasach nSelBL: 427929 bez trendu w granicach błędu; 428770 silniejszy przy nSelBL < 20 (0,33 ± 0,10 i 0,64 ± 0,12), dla nSelBL ≥ 20 w zakresie 0,78–0,95. Drobniejsze warstwy nSelBL zmieniają N_acc o +1,1 ± 0,7 % (428770) i +3 ± 5 % (427929), czyli resztkowa zależność od aktywności wewnątrz klas jest mała. Dawcy z protonem po stronie przeciwnej mają częstość tagu wielokrotnie większą od przewidywanej (428770, z D, strony razem: 54 tagów wobec 6.9 przewidywanych — ten sam efekt co G9); ich usunięcie obniża N_acc o 1,8 ± 0,5 % i 3,3 ± 2,8 %.
Wynik poboczny: przy ustalonych warstwach LB R rośnie z μ zdarzenia o czynnik około 1,6–2 między skrajnymi tercylami (zarówno przy małym, jak i dużym ξ_cal,o). actual_mu zmienia się głównie w obrębie LB (rozrzut względny w LB: 0,86 w 427929 i 0,53 w 428770; między LB: 0,09 i 0,03), więc warstwy LB nie warstwują μ. N_acc jest na to niewrażliwe tylko dlatego, że biorcy i dawcy mają w LB ten sam rozkład μ (wariant v2 „warstwy μ”: 1,007 ± 0,003 i 1,035 ± 0,017).

### 2.4 Wniosek i propozycja systematyki

Wyniki (pomiary):
1. Mechanizm ślepy na stronę: skok w górę wykluczony; κ_L = 0.79 ± 0.04 (428770), 0.83 ± 0.12 (427929); granice 95 % CL: κ_L < 0.86 i < 1.03.
2. Strona tagu, ξ_cal 0,09–0,14: R = 1.00 ± 0.05 i 1.14 ± 0.12; skok < 8 % i < 33 % (95 % CL).
3. Nadwyżka 1,14 ± 0,06 w 428770 nie wygląda na skok częstości przypadkowych (2.2).

Hipotezy (niepotwierdzone):
- Jeśli zależność z testu (i) jest ślepa na stronę, nominał **zawyża** N_acc o około 20 % w obu runach. Test symetrii (2,4σ w 428770, zdarzenia z D) przemawia przeciw, ale nie rozstrzyga.
- Jeśli efekt jest związany ze stroną przeciwną (np. proton przypadkowy z oddziaływania, które zostawia energię w kalorymetrze po stronie o i wypycha zdarzenie z małych ξ_cal,o), nominał jest obciążony mniej; znak ewentualnej poprawki (depozyt także po stronie tagu → mniejsza wydajność cięcia) jest ten sam. Zależność niedoboru od μ, która by to potwierdziła, nie jest widoczna w granicach błędu.
- Mechanizm działający tylko po stronie tagu i tylko dla ξ_cal,s < 0,06 pozostaje **nieograniczony przez dane** (bez obszaru sygnału nie ma próbki). Za jego brakiem przemawia tylko argument: proton niezależny nie „wie”, po której stronie jest przerwa.

**PROPOZYCJA** (kod metody bez zmian): składnik „ekstrapolacja” zastąpić asymetryczną obwiednią z pomiarów — w dół √((1 − κ_L)² + σ²), w górę √((R_daleki − 1)² + σ²); bez korekty wartości centralnej. Pozostałe składniki (stat, pula, ½Δ metod) bez zmian.

| run | f_acc | v2: ekstrapolacja | v2: razem | **propozycja: zależność od ξ_cal** | **propozycja: razem** | **f_acc, propozycja** | v2 |
|---|---|---|---|---|---|---|---|
| 427929 | 0.0087 | ± 47.1 % | ± 49.9 % | +18.3 % / −21.4 % | +24.5 % / −26.9 % | **0.0087 +0.0021 / −0.0024** | 0.0087 ± 0.0044 |
| 428770 | 0.0477 | ± 10.4 % | ± 28.4 % | +4.9 % / −21.5 % | +26.9 % / −34.0 % | **0.0477 +0.0128 / −0.0162** | 0.0477 ± 0.0135 |

- 427929: dane **pozwalają zawęzić** błąd: ± 49.9 % → +24.5 % / −26.9 %. Ekstrapolacja v2 (47 %) była zdominowana błędem statystycznym nachylenia z ~440 tagów puli na ramieniu 1,45 w ln ξ_cal; estymator lustrzany mierzy tę samą wielkość bez ekstrapolacji.
- 428770: dane **nie pozwalają zawęzić** błędu. Składnik w dół rośnie z 10,4 % do 21.5 %, bo test (i) pokazuje istotny (4.8σ) efekt, którego ekstrapolacja liniowa po stronie tagu nie widzi; w górę maleje do 4.9 %. Błąd całkowity dalej zdominowany przez ½Δ metod (B(i), 25,8 %): +26.9 % / −34.0 % wobec ± 28.4 %.
- Wariant symetryczny (większa z obu stron): ± 26.9 % (427929), ± 34.0 % (428770).

## 3. Zastrzeżenia

- Wszystko WSTĘPNE; zakres tylko trygery nie-AFP; kandydaci bez score. Po selekcji score f_acc i kształty trzeba przeliczyć.
- Systematyki w szablonie są wspólne dla runu (z regionu A∨C, inkluzywnie) i przypisane także stronom A, C i przedziałom pT; ekstrapolacja per strona w v2 różni się (427929: κ_A = 0,94 ± 0,37, κ_C = 0,22 ± 0,43).
- Kolumny błędów w `template.csv` to systematyki v2, nie propozycja z 2.4 (propozycja: `jump_summary.json`, pole `proposal`).
- Propozycja w górę opiera się na teście (ii) w zakresie 0,09–0,14 i na założeniu braku mechanizmu właściwego stronie tagu poniżej 0,06; to założenie, nie pomiar.
- Testy (i)–(iii) i warianty mieszania dzielą dawców z nominałem; błędy stosunków pochodzą ze wspólnego bootstrapu (B = 500), błędy R z Poissona obs ⊕ 1/√(tagi puli).
- Test symetrii stron: błąd różnicy nachyleń liczony jak dla wielkości niezależnych (te same tagi, inne binowanie) — przybliżenie.
- B = 500 zamiast 200 w v2: błędy stat. N_acc różnią się od v2 o kilka procent względnie (np. 428770 A∨C: 27,2 wobec 25,3).

## 4. Pliki

| plik | zawartość |
|---|---|
| `template.json`, `template.csv`, `template_shapes.csv` | szablon (schemat: `README.md`) |
| `jump_tests.json`, `jump_summary.json` | surowe liczności testów (i)–(iii) i wielkości pochodne z propozycją |
| `common.py`, `make_template.py`, `jump_tests.py`, `derive.py`, `make_report.py` | kod (selekcja = kopia logiki `run_data.py` v2; metoda importowana bez zmian) |
| `make_template.log`, `jump_tests.log` | logi biegów |
