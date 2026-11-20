# AFP-MIX v2 na danych — WSTĘPNE

**STATUS: WSTĘPNE.** Wygenerowano 2026-11-18 19:17 przez `AFP_mix/run_data.py` v2 (czas 239 s). Wejście: `/home/k/Documents/charm/charm_work/AFP_mix/runs/lxplus`.
Metoda: `charmbkg/afp_mix.py` v2, `docs/afp_mixing.md` (charmbkg); przegląd `REV_afp_mix/RAPORT.md`, literatura `RECON_lit_B/RAPORT.md`. Kandydaci **bez score** (luźna selekcja v2.0) — liczby pod pikiem to wszyscy kandydaci w oknie (sygnał + kombinatoryka), bez fitu. Raport zawiera wyłącznie frakcje tła przypadkowego tagu i liczności tagów; uzysków sygnału dyfrakcyjnego ani przekrojów nie liczono.

## 0. Zakres ważności

**Wynik dotyczy wyłącznie zdarzeń wybranych trygerami nie-AFP** (GRL ∧ ≥ 1 zapisany łańcuch bez L1AFP ∧ nPileupVertices == 0). Zdarzenia wybrane wyłącznie przez łańcuchy L1AFP (tryger wymaga protonu) są wykluczone u biorców i dawców; f_acc z tej metody ich NIE opisuje (inna aktywność, inne ξ_cal, tag wymuszony trygerem). Liczby per próbka poniżej („Wejścia”). Próbek nie łączono.

## 1. Wynik: f_acc (region A∨C, po cięciu ξ_cal) z rozbiciem błędu

f_acc = N_acc/N_tag, estymator A (mieszanie w warstwach LB × nSelBL × „z D/bez D”, pule jednostronne ξ_cal > 0,14). Składniki błędu względnego N_acc: stat = bootstrap otagowanych dawców; pula = największa istotna różnica wariantów puli (|Δ| − σ_Δ, pełna statystyka, wspólny bootstrap); ekstrapolacja = √((κ−1)² + σ_κ²) z liniowego modelu R(ln ξ_cal) w puli; ½Δ metod = połowa największej różnicy A–B(i), A–B(ii).

| próbka | N_tag | N_acc (A) | **f_acc** | stat | pula | ekstrapolacja | ½Δ metod | razem (wzgl.) | **f_acc ± razem** | v1 f_acc | werdykt bramek v2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| data_428770 | 15833 | 754.0 ± 24.6 | **0.0476** | 3.3 % | 4.5 % | 9.3 % | 7.0 % | 12.9 % | **0.0476 ± 0.0061** | 0.0459 | **PASS**; N/A: G9_opposite_side |

Uwaga: składniki „ekstrapolacja” i „½Δ metod (B(ii))” dotyczą tego samego pytania (zależność częstości od ξ_cal) — suma kwadratowa jest konserwatywna.

## data_428770 (run 428770)

### Wejścia, selekcja, zakres

- zdarzeń 7,989,496; GRL 7,989,496; LB 4–323 (136 bloków); μ (5/50/95 %) = 0.0224/0.0512/0.0849
- **zakres (W2):** zdarzenia wyłącznie z łańcuchów L1AFP po GRL: 890,219, w tym z tagiem 738,763 — **wykluczone**; zdarzenia z tagiem w wybranej próbce (nie-AFP, npu0): 32,020 → metoda opisuje 4.2 % wszystkich tagów po GRL; bez żadnego zapisanego łańcucha: 1,469,105
- po GRL ∧ tryger: 5,630,172; z nPileupVertices == 0: 5,350,292; dawcy 5,350,292; biorcy (≥ 1 kandydat slim) 4,403,767
- tryger losowy HLT_noalg_L1RD0_FILLED: **0 zdarzeń** → absolutne P(tag) z trygera losowego niedostępne; B(i) używa: łańcuchy L1MBTS quasi-minbias (z D i bez D) (48,341 zdarzeń; łańcuchy L1MBTS to 0.9 % wybranych)
- klasy trygera (3 = L1J12, 2 = j20, 1 = sptrk_pt4_L1MBTS_2, 0 = pozostałe L1MBTS): biorcy 2: 0.172, 3: 0.828; dawcy 2: 0.203, 3: 0.797
- liczności tagów (przed cięciem ξ_cal): dawcy A 16,408 / C 15,760 / A∧C 148 z 5,350,292; biorcy A 10,421 / C 10,029 / A∧C 84 z 4,403,767

### Estymatory (N_acc, region A∨C)

| estymator | N_acc | N/N_A | uwagi |
|---|---|---|---|
| **A** — mieszanie, warstwy LB × nSelBL × „z D”, pule jednostronne ξ_cal > 0.14 (nominał) | 754.0 ± 24.6 | 1 | otagowanych dawców 4069; ε_ξcal = 0.147; f̄ = 0.0012 |
| informacyjnie: warstwy LB × nSelBL, dawcy z D i bez D | 869.3 | 1.153 ± 0.027 | nie w systematyce — efekt dawców bez D (W1; zob. I2) |
| informacyjnie: tylko warstwy LB, wszyscy dawcy (≈ v1, pule jednostronne) | 725.8 | 0.963 ± 0.027 | nie w systematyce — wielkość efektu warstwowania |
| v1 (raport `out/RAPORT_DANE.md`) | 736.9 ± 12.2 | 0.977 | f_acc v1 = 0.0459; v1 syst. dekl. 8.3 % |
| **B(i)** — druga pula: łańcuchy L1MBTS quasi-minbias (z D i bez D), weto ξ_cal po obu stronach, nSelBL ≥ 20, absolutne P(tag) w 6 przedziałach μ | 859.9 ± 205.1 | 1.140 ± 0.272 | otagowanych dawców 25 z 17,966 |
| **B(ii)** — A × N_obs/N_pred w obszarze kontrolnym (ξ_cal ≤ 0.12, Δ > 0.045) | 838.8 ± 52.1 | 1.113 ± 0.061 | CR: obs 394, pred 354.1 |

Rozbicie błędu względnego N_acc: stat 3.3 %, pula 4.5 %, ekstrapolacja 9.3 % (κ = 1.017 ± 0.092), ½Δ metod 7.0 % (Bi: 7.0 %, Bii: 5.6 %) → **razem 12.9 %**.

Warianty puli (pełna statystyka; błąd = błąd różnicy względem nominału ze wspólnego bootstrapu):

| wariant | N_acc / nominał | istotna część |
|---|---|---|
| margin_0 | 1.000 ± 0.015 | 0.0 % |
| margin_0.05 | 0.973 ± 0.024 | 0.2 % |
| margin_0.1 | 0.908 ± 0.048 | 4.5 % |
| strata_mu | 1.007 ± 0.003 | 0.4 % |
| act_coarse | 0.989 ± 0.011 | 0.0 % |
| plus_trig_class | 1.008 ± 0.005 | 0.3 % |

Profil częstości tagu w puli po warstwowaniu (R = obs/pred; G1 i ekstrapolacja):

| strona | przedziały ξ_cal | tagi k | pred | nachylenie dR/d ln ξ_cal | ⟨ln ξ_cal⟩ pula → użycie | κ |
|---|---|---|---|---|---|---|
| A | 0.14, 0.2, 0.3, 0.5, 1, inf | 605, 671, 316, 15, 0 | 602.4, 670.3, 319.8, 14.6, 0.0 | -0.016 ± 0.087 | -1.47 → -2.96 | 1.023 ± 0.129 |
| C | 0.14, 0.2, 0.3, 0.5, 1, inf | 610, 670, 328, 17, 0 | 605.4, 675.4, 328.4, 15.9, 0.0 | -0.007 ± 0.086 | -1.47 → -2.98 | 1.010 ± 0.129 |

Obszar kontrolny — skan dolnej granicy Δ − δ (spadek R z marginesem = przeciek własnych protonów przez rozdzielczość ξ_cal):

| margines | N_obs | N_pred | R |
|---|---|---|---|
| 0.010 | 678 | 575.1 | 1.179 ± 0.049 |
| 0.025 | 394 | 354.1 | 1.113 ± 0.061 |
| 0.040 | 177 | 175.1 | 1.011 ± 0.083 |
| 0.060 | 20 | 23.3 | 0.860 ± 0.215 |

### Bramki v2

| bramka | status | kluczowe | co testuje | czego NIE testuje |
|---|---|---|---|---|
| G0_inputs | **PASS** | coverage=1, frac_no_pool=0, median_pool=4557, n_afp_only_excluded=8.902e+05; zdarzenia wyłącznie z trygera L1AFP wykluczone u biorców i dawców — dla nich wynik NIE obowiązuje (W2) | pokrycie warstw (≥ 99 % biorców w komórce zachowującej klasę aktywności: LB×klasa → μ×klasa → run×klasa; 0 bez puli), brak zdarzeń wyłącznie z trygera AFP | poprawności fizycznej; zakresu poza trygerami nie-AFP (W2) — tam metoda nie obowiązuje |
| G1_pool_homogeneity | **PASS** | lr=0.2199, ndf=6, p=0.9998, n_pool_tags=3232 | czy po warstwowaniu częstość tagu w puli nie zależy od ξ_cal dawcy POWYŻEJ progu (Poisson LR obs vs pred w przedziałach ξ_cal) | obszaru poniżej progu (obszar użycia); zależności, które znoszą się w przedziale |
| G2_disjoint_closure | **PASS** | n_obs=2662, n_pred=2764, n_pred_err=124.3, ratio=0.963, ratio_err=0.04884, z=-0.7576; strony: A: 1305/1474.5, C: 1357/1289.7 | domknięcie na zbiorach ROZŁĄCZNYCH w obszarze anty-ξ_cal: tagi zdarzeń puli z najliczniejszej klasy trygera vs przewidywanie z dawców pozostałych klas (te same warstwy LB × aktywność × D) | obszaru użycia (ξ_cal małe); zależności od ξ_cal; przy jednej klasie trygera N/A (v1: G2 była w ~90 % tożsamością) |
| G5_cr_shape | **PASS** | lr=9.677, n_bins=16, p=0.8436, n_obs=394, n_pred=354.1 | kształt (Δ, ξ_AFP) tagów w obszarze kontrolnym Δ > δ + 0,025 u biorców SPOZA puli (ξ_cal ≤ 0,12) vs mieszanie; Poisson LR, p z pseudoeksperymentów z bootstrapem dawców | normalizacji (to G7); obszaru Δ ≤ δ; przy < ~50 zdarzeniach nie ma mocy |
| G7_methods_A_vs_Bii | **PASS** | n_obs=394, n_pred=354.1, n_pred_err=10.6, ratio=1.113, ratio_err=0.06099, z=1.846; strony: A: 205/177.2, C: 189/176.9 | GŁÓWNA: częstość przypadkowych najbliżej obszaru użycia — N_obs/N_pred w obszarze kontrolnym (ξ_cal 0,075–0,12, Δ > 0,045; „przesuwane weto”) = stosunek metod B(ii)/A | obszaru ξ_cal < 0,075 (tam prawdziwe tagi są nieodróżnialne); przeciek własnych protonów przez rozdzielczość ξ_cal zawyża R (skan marginesu) |
| G8_methods_A_vs_Bi | **PASS** | ratio=1.14, ratio_err=0.2718, z=0.5166, n_acc_A=754, n_acc_Bi=859.9, n_tagged_donors=25 | GŁÓWNA: zgodność z drugą metodą — absolutne P(tag) z próbki quasi-minbias o pełnej aktywności (weto ξ_cal po obu stronach, nSelBL ≥ 20) w przedziałach μ, bez warstw aktywności (styl ALFA: pełna niezależność) | obszaru użycia (obie pule mają ξ_cal > próg); przy tej samej próbce trygerowej (427929) pule częściowo się pokrywają |
| G9_opposite_side | **N/A** | n_obs=10, n_pred=3.971, n_pred_err=0.06986, ratio=2.518, ratio_err=0.5021, z=2.418; strony: A: 5/1.9, C: 5/2.1 | obszar „dwuprotonowy”: zdarzenia z prawdziwym tagiem (po cięciu) po stronie przeciwnej i ξ_cal > próg po stronie badanej — częstość tagu po stronie badanej vs mieszanie (czy zdarzenia dyfrakcyjne mają tę samą częstość przypadkowych) | zależności od ξ_cal po stronie badanej; CD z protonem poza oknem nie szkodzi, ale CD z rozdzielczością ξ_cal zawyża obs |
| I2_hasD_vs_noD | **INFO (informacyjna)** | ndf=9, p=4.794e-06, n_obs=3232, n_pred=3257 | (informacyjna) skład puli W1: tagi zdarzeń z D w obszarze anty vs przewidywanie z dawców BEZ D, per klasa aktywności — uzasadnia klasę „z D/bez D” w warstwach | niczego o estymacie nominalnej (dawcy bez D nie są w niej używani) — nie wchodzi do werdyktu |
| I3_mu_slope | **INFO (informacyjna)** | slope_data=0.0005216, slope_data_err=0.001424, slope_mix=7.556e-05, z_slope=0.3133 | (informacyjna) nachylenie f_tag(μ) dane vs mieszanie | niczego przy σ_b(dane) ≫ b_mix — nie wchodzi do werdyktu |
| I6_eta_asymmetry | **INFO (informacyjna)** | z=-2.828, A_obs=0.2665, A_err=0.003565, f_acc_eta=-0.3326, f_acc_eta_err=0.1344, f_acc_mix=0.04762, n_cand=7.314e+04 | (informacyjna) f_acc z asymetrii η* kandydatów przy założonym A_diff | A_diff z Pythii nie pasuje do kandydatów inkluzywnych — nie wchodzi do werdyktu |

### f_acc inkluzywnie per region

| region | N_tag | N_acc (A) | f_acc ± stat | ε_ξcal(przyp.) | f̄ | v1 f_acc |
|---|---|---|---|---|---|---|
| AorC | 15833 | 754.0 ± 24.6 | 0.0476 ± 0.0016 | 0.147 | 0.0012 | 0.0459 |
| A | 8140 | 367.6 ± 18.2 | 0.0452 ± 0.0023 | 0.144 | 0.0008 | 0.0456 |
| C | 7756 | 390.2 ± 19.0 | 0.0503 ± 0.0025 | 0.151 | 0.0008 | 0.0463 |

### f_acc(strona, μ (actual_mu)) — błąd tylko stat.

| μ (actual_mu) | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C |
|---|---|---|---|---|---|---|---|
| 0.01037–0.03488 | 0.02555 | 1358 | 0.0440 ± 0.0025 | 1317 | 0.0497 ± 0.0038 | 2663 | 0.0468 ± 0.0022 |
| 0.03488–0.04455 | 0.04051 | 1351 | 0.0445 ± 0.0026 | 1289 | 0.0506 ± 0.0029 | 2630 | 0.0474 ± 0.0018 |
| 0.04455–0.05091 | 0.04786 | 1290 | 0.0474 ± 0.0027 | 1305 | 0.0506 ± 0.0032 | 2587 | 0.0489 ± 0.0020 |
| 0.05091–0.05803 | 0.05421 | 1355 | 0.0454 ± 0.0026 | 1262 | 0.0511 ± 0.0031 | 2607 | 0.0481 ± 0.0019 |
| 0.05803–0.06924 | 0.0634 | 1406 | 0.0437 ± 0.0025 | 1273 | 0.0497 ± 0.0028 | 2671 | 0.0465 ± 0.0018 |
| 0.06924–0.2772 | 0.08923 | 1380 | 0.0461 ± 0.0026 | 1310 | 0.0501 ± 0.0033 | 2675 | 0.0481 ± 0.0020 |

### f_acc(strona, lumiBlock) — błąd tylko stat.

| lumiBlock | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C |
|---|---|---|---|---|---|---|---|
| 4–20 | 11.89 | 999 | 0.0459 ± 0.0068 | 932 | 0.0536 ± 0.0072 | 1923 | 0.0495 ± 0.0048 |
| 20–221 | 147.2 | 1014 | 0.0478 ± 0.0067 | 991 | 0.0494 ± 0.0070 | 1998 | 0.0485 ± 0.0049 |
| 221–238 | 229 | 1000 | 0.0417 ± 0.0058 | 987 | 0.0585 ± 0.0076 | 1982 | 0.0499 ± 0.0048 |
| 238–256 | 246.5 | 1048 | 0.0510 ± 0.0063 | 963 | 0.0470 ± 0.0059 | 1999 | 0.0491 ± 0.0043 |
| 256–273 | 264 | 994 | 0.0361 ± 0.0056 | 979 | 0.0464 ± 0.0060 | 1967 | 0.0412 ± 0.0042 |
| 273–290 | 281 | 1077 | 0.0406 ± 0.0058 | 978 | 0.0453 ± 0.0053 | 2048 | 0.0428 ± 0.0038 |
| 290–306 | 297.5 | 972 | 0.0451 ± 0.0069 | 880 | 0.0527 ± 0.0083 | 1848 | 0.0486 ± 0.0052 |
| 306–323 | 314.4 | 1036 | 0.0528 ± 0.0067 | 1046 | 0.0499 ± 0.0070 | 2068 | 0.0514 ± 0.0052 |

### Gatunki × pT: tag przypadkowy pod pikiem i w sidebandach (A∨C; kandydaci, bez fitu)

N_tag = kandydaci w oknie w zdarzeniach z tagiem po cięciu ξ_cal; N_acc = Σ p_acc po kandydatach. f_acc ± stat ± syst (syst wspólna, względna 12.5 %: pula ⊕ ekstrapolacja ⊕ ½Δ metod). Sidebandy NIE są próbką „bez prawdziwych tagów” (kombinatoryka masowa w zdarzeniach dyfrakcyjnych niesie prawdziwy proton) — stosunek pik/SB jest tylko informacyjny.

**Dplus** — Dmeson_mass: pik (1.84, 1.9), SB ((1.75, 1.8), (1.94, 2.0))

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik ± stat ± syst | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | v1 f_acc pik |
|---|---|---|---|---|---|---|---|---|
| 3–150 | 8757 | 441.2 ± 10.1 | 0.0504 ± 0.0015 ± 0.0063 | 0.0481 ± 0.0020 | 0.0528 ± 0.0023 | 0.0507 ± 0.0014 | 0.993 ± 0.040 | 0.0537 |
| 3–4 | 4022 | 218.8 ± 5.0 | 0.0544 ± 0.0017 ± 0.0068 | 0.0522 ± 0.0023 | 0.0566 ± 0.0026 | 0.0542 ± 0.0016 | 1.004 ± 0.043 | 0.0586 |
| 4–5 | 1811 | 89.7 ± 2.1 | 0.0495 ± 0.0019 ± 0.0062 | 0.0459 ± 0.0024 | 0.0534 ± 0.0029 | 0.0522 ± 0.0018 | 0.948 ± 0.048 | 0.0527 |
| 5–6 | 927 | 43.0 ± 1.0 | 0.0463 ± 0.0021 ± 0.0058 | 0.0438 ± 0.0027 | 0.0492 ± 0.0032 | 0.0466 ± 0.0018 | 0.995 ± 0.059 | 0.0490 |
| 6–8 | 853 | 40.3 ± 1.0 | 0.0473 ± 0.0022 ± 0.0059 | 0.0452 ± 0.0029 | 0.0496 ± 0.0034 | 0.0459 ± 0.0018 | 1.031 ± 0.063 | 0.0501 |
| 8–12 | 672 | 30.3 ± 0.8 | 0.0451 ± 0.0023 ± 0.0056 | 0.0449 ± 0.0033 | 0.0455 ± 0.0034 | 0.0436 ± 0.0019 | 1.033 ± 0.071 | 0.0467 |
| 12–150 | 472 | 19.2 ± 0.5 | 0.0406 ± 0.0026 ± 0.0051 | 0.0395 ± 0.0035 | 0.0416 ± 0.0038 | 0.0429 ± 0.0024 | 0.946 ± 0.080 | 0.0419 |

**D0** — Dmeson_mass: pik (1.835, 1.895), SB ((1.75, 1.8), (1.93, 2.0))

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik ± stat ± syst | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | v1 f_acc pik |
|---|---|---|---|---|---|---|---|---|
| 3–150 | 6404 | 320.6 ± 8.2 | 0.0501 ± 0.0015 ± 0.0062 | 0.0475 ± 0.0021 | 0.0529 ± 0.0024 | 0.0495 ± 0.0014 | 1.011 ± 0.042 | 0.0515 |
| 3–4 | 2763 | 141.2 ± 3.5 | 0.0511 ± 0.0017 ± 0.0064 | 0.0491 ± 0.0024 | 0.0533 ± 0.0026 | 0.0510 ± 0.0016 | 1.001 ± 0.046 | 0.0530 |
| 4–5 | 1473 | 73.2 ± 1.9 | 0.0497 ± 0.0019 ± 0.0062 | 0.0452 ± 0.0025 | 0.0549 ± 0.0032 | 0.0504 ± 0.0017 | 0.987 ± 0.051 | 0.0510 |
| 5–6 | 692 | 36.5 ± 1.0 | 0.0527 ± 0.0026 ± 0.0066 | 0.0522 ± 0.0036 | 0.0535 ± 0.0038 | 0.0488 ± 0.0019 | 1.081 ± 0.068 | 0.0539 |
| 6–8 | 639 | 32.7 ± 0.9 | 0.0512 ± 0.0025 ± 0.0064 | 0.0503 ± 0.0036 | 0.0522 ± 0.0037 | 0.0485 ± 0.0019 | 1.057 ± 0.067 | 0.0520 |
| 8–12 | 499 | 22.6 ± 0.6 | 0.0454 ± 0.0025 ± 0.0057 | 0.0401 ± 0.0031 | 0.0516 ± 0.0042 | 0.0464 ± 0.0021 | 0.978 ± 0.069 | 0.0461 |
| 12–150 | 338 | 14.4 ± 0.4 | 0.0425 ± 0.0027 ± 0.0053 | 0.0412 ± 0.0038 | 0.0434 ± 0.0039 | 0.0411 ± 0.0021 | 1.035 ± 0.085 | 0.0428 |

**Dstar** — Dmeson_delta_mass: pik (0.144, 0.147), SB ((0.15, 0.16),)

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik ± stat ± syst | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | v1 f_acc pik |
|---|---|---|---|---|---|---|---|---|
| 3–150 | 570 | 27.2 ± 0.7 | 0.0478 ± 0.0024 ± 0.0060 | 0.0451 ± 0.0032 | 0.0505 ± 0.0037 | 0.0520 ± 0.0018 | 0.919 ± 0.056 | 0.0507 |
| 3–4 | 224 | 11.3 ± 0.3 | 0.0504 ± 0.0037 ± 0.0063 | 0.0488 ± 0.0051 | 0.0525 ± 0.0055 | 0.0550 ± 0.0025 | 0.916 ± 0.079 | 0.0537 |
| 4–5 | 166 | 7.4 ± 0.2 | 0.0447 ± 0.0039 ± 0.0056 | 0.0401 ± 0.0047 | 0.0495 ± 0.0064 | 0.0526 ± 0.0026 | 0.849 ± 0.085 | 0.0471 |
| 5–6 | 71 | 3.5 ± 0.1 | 0.0486 ± 0.0062 ± 0.0061 | 0.0438 ± 0.0077 | 0.0552 ± 0.0100 | 0.0471 ± 0.0029 | 1.033 ± 0.146 | 0.0517 |
| 6–8 | 60 | 2.7 ± 0.1 | 0.0453 ± 0.0061 ± 0.0057 | 0.0497 ± 0.0098 | 0.0410 ± 0.0075 | 0.0501 ± 0.0036 | 0.905 ± 0.138 | 0.0480 |
| 8–12 | 26 | 1.5 ± 0.1 | 0.0558 ± 0.0115 ± 0.0070 | 0.0467 ± 0.0130 | 0.0629 ± 0.0185 | 0.0521 ± 0.0051 | 1.072 ± 0.246 | 0.0624 |
| 12–150 | 23 | 0.9 ± 0.0 | 0.0394 ± 0.0091 ± 0.0049 | 0.0393 ± 0.0131 | 0.0399 ± 0.0127 | 0.0392 ± 0.0050 | 1.005 ± 0.266 | 0.0385 |

## Zmiany względem v1

- Estymator: warstwy (LB × nSelBL × „z D/bez D”) zamiast samych LB; region A∨C z pul jednostronnych (v1: pula dwustronna); overlay „union” na parach (bez iteracyjnej poprawki własnego okna); pełna pula na poziomach zapasowych.
- Błąd: bootstrap B = 200 wspólny dla wariantów; warianty puli na pełnej statystyce z błędem różnicy (v1: podpróbka → szum, np. 32 % w 427929); nowe składniki: ekstrapolacja w ξ_cal, połowa różnicy metod A–B(i), A–B(ii).
- Bramki: G1 → Poisson LR po warstwowaniu (odporna na zera); G2 → zbiory rozłączne (klasy trygera); I2 (informacyjna): z D vs bez D per klasa nSelBL; G4 usunięta; G5 → kształt w obszarze kontrolnym poza pulą, Poisson LR z pseudoeksperymentami, próg nominalny p ≥ 0,01; nowe G7 (obszar kontrolny / „przesuwane weto” = A vs B(ii)), G8 (A vs B(i)), G9 (obszar dwuprotonowy); G3/G6 → informacyjne I3/I6.

## Ograniczenia i sprawy otwarte (WSTĘPNE)

- **Zakres (W2):** tylko trygery nie-AFP. Dla próbki z trygera L1AFP potrzebni są biorcy reprezentatywni dla tej klasy (decyzja WP2).
- **Obszar ξ_cal < 0,075 pozostaje nietestowalny bezpośrednio** (prawdziwe i przypadkowe tagi nieodróżnialne). G7 sięga do ξ_cal ≈ 0,075, G9 testuje zdarzenia dyfrakcyjne po drugiej stronie; toy: skok częstości poniżej 0,07 nie jest wykrywany przez żadną bramkę (docs/afp_mixing.md, tabela mocy).
- Brak trygera losowego w zapisanych łańcuchach (0 zdarzeń RD0_FILLED) — B(i) jest drugą pulą z łańcuchów L1MBTS, nie absolutnym P(tag) z pustych przecięć.
- Brak score CharmNet/XGB w h5 v2.0 → kandydaci po luźnej selekcji; f_acc pod pikiem dotyczy mieszanki sygnał + kombinatoryka.
- Fluktuacje mianownika (liczba dawców w komórce) pominięte w bootstrapie (względnie ~ częstość tagu, ≲ 10⁻³).
