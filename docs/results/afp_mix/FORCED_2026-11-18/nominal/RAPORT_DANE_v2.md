# AFP-MIX v2 na danych — WSTĘPNE

**STATUS: WSTĘPNE.** Wygenerowano 2026-11-18 19:14 przez `AFP_mix/run_data.py` v2 (czas 142 s). Wejście: `runs/lxplus`.
Metoda: `charmbkg/afp_mix.py` v2, `docs/afp_mixing.md` (charmbkg); przegląd `REV_afp_mix/RAPORT.md`, literatura `RECON_lit_B/RAPORT.md`. Kandydaci **bez score** (luźna selekcja v2.0) — liczby pod pikiem to wszyscy kandydaci w oknie (sygnał + kombinatoryka), bez fitu. Raport zawiera wyłącznie frakcje tła przypadkowego tagu i liczności tagów; uzysków sygnału dyfrakcyjnego ani przekrojów nie liczono.

## 0. Zakres ważności

**Wynik dotyczy wyłącznie zdarzeń wybranych trygerami nie-AFP** (GRL ∧ ≥ 1 zapisany łańcuch bez L1AFP ∧ nPileupVertices == 0). Zdarzenia wybrane wyłącznie przez łańcuchy L1AFP (tryger wymaga protonu) są wykluczone u biorców i dawców; f_acc z tej metody ich NIE opisuje (inna aktywność, inne ξ_cal, tag wymuszony trygerem). Liczby per próbka poniżej („Wejścia”). Próbek nie łączono.

## 1. Wynik: f_acc (region A∨C, po cięciu ξ_cal) z rozbiciem błędu

f_acc = N_acc/N_tag, estymator A (mieszanie w warstwach LB × nSelBL × „z D/bez D”, pule jednostronne ξ_cal > 0,14). Składniki błędu względnego N_acc: stat = bootstrap otagowanych dawców; pula = największa istotna różnica wariantów puli (|Δ| − σ_Δ, pełna statystyka, wspólny bootstrap); ekstrapolacja = √((κ−1)² + σ_κ²) z liniowego modelu R(ln ξ_cal) w puli; ½Δ metod = połowa największej różnicy A–B(i), A–B(ii).

| próbka | N_tag | N_acc (A) | **f_acc** | stat | pula | ekstrapolacja | ½Δ metod | razem (wzgl.) | **f_acc ± razem** | v1 f_acc | werdykt bramek v2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| data_428770 | 16065 | 766.1 ± 25.3 | **0.0477** | 3.3 % | 4.8 % | 10.4 % | 25.8 % | 28.4 % | **0.0477 ± 0.0135** | 0.0459 | **PASS**; N/A: G9_opposite_side |

Uwaga: składniki „ekstrapolacja” i „½Δ metod (B(ii))” dotyczą tego samego pytania (zależność częstości od ξ_cal) — suma kwadratowa jest konserwatywna.

## data_428770 (run 428770)

### Wejścia, selekcja, zakres

- zdarzeń 7,989,496; GRL 7,989,496; LB 4–323 (136 bloków); μ (5/50/95 %) = 0.0224/0.0512/0.0849
- **zakres (W2):** zdarzenia wyłącznie z łańcuchów L1AFP po GRL: 888,801, w tym z tagiem 737,591 — **wykluczone**; zdarzenia z tagiem w wybranej próbce (nie-AFP, npu0): 33,592 → metoda opisuje 4.4 % wszystkich tagów po GRL; bez żadnego zapisanego łańcucha: 1,364,623
- po GRL ∧ tryger: 5,736,072; z nPileupVertices == 0: 5,449,549; dawcy 5,449,549; biorcy (≥ 1 kandydat slim) 4,475,334
- tryger losowy HLT_noalg_L1RD0_FILLED: **0 zdarzeń** → absolutne P(tag) z trygera losowego niedostępne; B(i) używa: łańcuchy L1MBTS quasi-minbias (z D i bez D) (63,632 zdarzeń; łańcuchy L1MBTS to 1.2 % wybranych)
- klasy trygera (3 = L1J12, 2 = j20, 1 = sptrk_pt4_L1MBTS_2, 0 = pozostałe L1MBTS): biorcy 0: 0.001, 1: 0.015, 2: 0.170, 3: 0.814; dawcy 0: 0.002, 1: 0.016, 2: 0.199, 3: 0.782
- liczności tagów (przed cięciem ξ_cal): dawcy A 17,445 / C 16,558 / A∧C 411 z 5,449,549; biorcy A 10,609 / C 10,188 / A∧C 93 z 4,475,334

### Estymatory (N_acc, region A∨C)

| estymator | N_acc | N/N_A | uwagi |
|---|---|---|---|
| **A** — mieszanie, warstwy LB × nSelBL × „z D”, pule jednostronne ξ_cal > 0.14 (nominał) | 766.1 ± 25.3 | 1 | otagowanych dawców 4180; ε_ξcal = 0.146; f̄ = 0.0012 |
| informacyjnie: warstwy LB × nSelBL, dawcy z D i bez D | 881.9 | 1.151 ± 0.027 | nie w systematyce — efekt dawców bez D (W1; zob. I2) |
| informacyjnie: tylko warstwy LB, wszyscy dawcy (≈ v1, pule jednostronne) | 737.0 | 0.962 ± 0.027 | nie w systematyce — wielkość efektu warstwowania |
| v1 (raport `out/RAPORT_DANE.md`) | 736.9 ± 12.2 | 0.962 | f_acc v1 = 0.0459; v1 syst. dekl. 8.3 % |
| **B(i)** — druga pula: łańcuchy L1MBTS quasi-minbias (z D i bez D), weto ξ_cal po obu stronach, nSelBL ≥ 20, absolutne P(tag) w 6 przedziałach μ | 1160.7 ± 218.2 | 1.515 ± 0.288 | otagowanych dawców 43 z 23,239 |
| **B(ii)** — A × N_obs/N_pred w obszarze kontrolnym (ξ_cal ≤ 0.12, Δ > 0.045) | 857.7 ± 50.0 | 1.120 ± 0.059 | CR: obs 404, pred 360.9 |

Rozbicie błędu względnego N_acc: stat 3.3 %, pula 4.8 %, ekstrapolacja 10.4 % (κ = 1.051 ± 0.090), ½Δ metod 25.8 % (Bi: 25.8 %, Bii: 6.0 %) → **razem 28.4 %**.

Warianty puli (pełna statystyka; błąd = błąd różnicy względem nominału ze wspólnego bootstrapu):

| wariant | N_acc / nominał | istotna część |
|---|---|---|
| margin_0 | 1.000 ± 0.015 | 0.0 % |
| margin_0.05 | 0.972 ± 0.023 | 0.5 % |
| margin_0.1 | 0.904 ± 0.048 | 4.8 % |
| strata_mu | 1.007 ± 0.003 | 0.4 % |
| act_coarse | 0.990 ± 0.011 | 0.0 % |
| plus_trig_class | 1.004 ± 0.007 | 0.0 % |

Profil częstości tagu w puli po warstwowaniu (R = obs/pred; G1 i ekstrapolacja):

| strona | przedziały ξ_cal | tagi k | pred | nachylenie dR/d ln ξ_cal | ⟨ln ξ_cal⟩ pula → użycie | κ |
|---|---|---|---|---|---|---|
| A | 0.14, 0.2, 0.3, 0.5, 1, inf | 627, 683, 325, 16, 0 | 616.6, 688.8, 330.6, 15.0, 0.0 | -0.037 ± 0.086 | -1.47 → -2.96 | 1.056 ± 0.127 |
| C | 0.14, 0.2, 0.3, 0.5, 1, inf | 624, 689, 334, 17, 0 | 616.5, 692.1, 339.2, 16.4, 0.0 | -0.030 ± 0.085 | -1.46 → -2.98 | 1.046 ± 0.128 |

Obszar kontrolny — skan dolnej granicy Δ − δ (spadek R z marginesem = przeciek własnych protonów przez rozdzielczość ξ_cal):

| margines | N_obs | N_pred | R |
|---|---|---|---|
| 0.010 | 692 | 585.5 | 1.182 ± 0.048 |
| 0.025 | 404 | 360.9 | 1.120 ± 0.059 |
| 0.040 | 182 | 178.9 | 1.017 ± 0.080 |
| 0.060 | 21 | 23.9 | 0.878 ± 0.211 |

### Bramki v2

| bramka | status | kluczowe | co testuje | czego NIE testuje |
|---|---|---|---|---|
| G0_inputs | **PASS** | coverage=1, frac_no_pool=0, median_pool=4633, n_afp_only_excluded=8.888e+05; zdarzenia wyłącznie z trygera L1AFP wykluczone u biorców i dawców — dla nich wynik NIE obowiązuje (W2) | pokrycie warstw (≥ 99 % biorców w komórce zachowującej klasę aktywności: LB×klasa → μ×klasa → run×klasa; 0 bez puli), brak zdarzeń wyłącznie z trygera AFP | poprawności fizycznej; zakresu poza trygerami nie-AFP (W2) — tam metoda nie obowiązuje |
| G1_pool_homogeneity | **PASS** | lr=0.5981, ndf=6, p=0.9964, n_pool_tags=3315 | czy po warstwowaniu częstość tagu w puli nie zależy od ξ_cal dawcy POWYŻEJ progu (Poisson LR obs vs pred w przedziałach ξ_cal) | obszaru poniżej progu (obszar użycia); zależności, które znoszą się w przedziale |
| G2_disjoint_closure | **PASS** | n_obs=2662, n_pred=2859, n_pred_err=104.9, ratio=0.9312, ratio_err=0.04118, z=-1.67; strony: A: 1305/1524.8, C: 1357/1333.9 | domknięcie na zbiorach ROZŁĄCZNYCH w obszarze anty-ξ_cal: tagi zdarzeń puli z najliczniejszej klasy trygera vs przewidywanie z dawców pozostałych klas (te same warstwy LB × aktywność × D) | obszaru użycia (ξ_cal małe); zależności od ξ_cal; przy jednej klasie trygera N/A (v1: G2 była w ~90 % tożsamością) |
| G5_cr_shape | **PASS** | lr=10.99, n_bins=16, p=0.7556, n_obs=404, n_pred=360.9 | kształt (Δ, ξ_AFP) tagów w obszarze kontrolnym Δ > δ + 0,025 u biorców SPOZA puli (ξ_cal ≤ 0,12) vs mieszanie; Poisson LR, p z pseudoeksperymentów z bootstrapem dawców | normalizacji (to G7); obszaru Δ ≤ δ; przy < ~50 zdarzeniach nie ma mocy |
| G7_methods_A_vs_Bii | **PASS** | n_obs=404, n_pred=360.9, n_pred_err=9.375, ratio=1.12, ratio_err=0.0587, z=2.037; strony: A: 211/180.6, C: 193/180.3 | GŁÓWNA: częstość przypadkowych najbliżej obszaru użycia — N_obs/N_pred w obszarze kontrolnym (ξ_cal 0,075–0,12, Δ > 0,045; „przesuwane weto”) = stosunek metod B(ii)/A | obszaru ξ_cal < 0,075 (tam prawdziwe tagi są nieodróżnialne); przeciek własnych protonów przez rozdzielczość ξ_cal zawyża R (skan marginesu) |
| G8_methods_A_vs_Bi | **PASS** | ratio=1.515, ratio_err=0.2884, z=1.786, n_acc_A=766.1, n_acc_Bi=1161, n_tagged_donors=43 | GŁÓWNA: zgodność z drugą metodą — absolutne P(tag) z próbki quasi-minbias o pełnej aktywności (weto ξ_cal po obu stronach, nSelBL ≥ 20) w przedziałach μ, bez warstw aktywności (styl ALFA: pełna niezależność) | obszaru użycia (obie pule mają ξ_cal > próg); przy tej samej próbce trygerowej (427929) pule częściowo się pokrywają |
| G9_opposite_side | **N/A** | n_obs=10, n_pred=4.076, n_pred_err=0.0708, ratio=2.453, ratio_err=0.4956, z=2.356; strony: A: 5/1.9, C: 5/2.1 | obszar „dwuprotonowy”: zdarzenia z prawdziwym tagiem (po cięciu) po stronie przeciwnej i ξ_cal > próg po stronie badanej — częstość tagu po stronie badanej vs mieszanie (czy zdarzenia dyfrakcyjne mają tę samą częstość przypadkowych) | zależności od ξ_cal po stronie badanej; CD z protonem poza oknem nie szkodzi, ale CD z rozdzielczością ξ_cal zawyża obs |
| I2_hasD_vs_noD | **INFO (informacyjna)** | ndf=9, p=5.545e-06, n_obs=3315, n_pred=3297 | (informacyjna) skład puli W1: tagi zdarzeń z D w obszarze anty vs przewidywanie z dawców BEZ D, per klasa aktywności — uzasadnia klasę „z D/bez D” w warstwach | niczego o estymacie nominalnej (dawcy bez D nie są w niej używani) — nie wchodzi do werdyktu |
| I3_mu_slope | **INFO (informacyjna)** | slope_data=0.0005679, slope_data_err=0.001413, slope_mix=7.132e-05, z_slope=0.3515 | (informacyjna) nachylenie f_tag(μ) dane vs mieszanie | niczego przy σ_b(dane) ≫ b_mix — nie wchodzi do werdyktu |
| I6_eta_asymmetry | **INFO (informacyjna)** | z=-2.825, A_obs=0.2664, A_err=0.003553, f_acc_eta=-0.3319, f_acc_eta_err=0.1344, f_acc_mix=0.04769, n_cand=7.363e+04 | (informacyjna) f_acc z asymetrii η* kandydatów przy założonym A_diff | A_diff z Pythii nie pasuje do kandydatów inkluzywnych — nie wchodzi do werdyktu |

### f_acc inkluzywnie per region

| region | N_tag | N_acc (A) | f_acc ± stat | ε_ξcal(przyp.) | f̄ | v1 f_acc |
|---|---|---|---|---|---|---|
| AorC | 16065 | 766.1 ± 25.3 | 0.0477 ± 0.0016 | 0.146 | 0.0012 | 0.0459 |
| A | 8271 | 376.5 ± 18.0 | 0.0455 ± 0.0022 | 0.144 | 0.0008 | 0.0456 |
| C | 7863 | 393.4 ± 19.6 | 0.0500 ± 0.0025 | 0.149 | 0.0008 | 0.0463 |

### f_acc(strona, μ (actual_mu)) — błąd tylko stat.

| μ (actual_mu) | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C |
|---|---|---|---|---|---|---|---|
| 0.01037–0.03488 | 0.02555 | 1380 | 0.0444 ± 0.0024 | 1343 | 0.0492 ± 0.0038 | 2709 | 0.0468 ± 0.0022 |
| 0.03488–0.04455 | 0.04052 | 1371 | 0.0449 ± 0.0025 | 1309 | 0.0502 ± 0.0029 | 2669 | 0.0475 ± 0.0018 |
| 0.04455–0.05092 | 0.04786 | 1314 | 0.0476 ± 0.0027 | 1318 | 0.0505 ± 0.0031 | 2623 | 0.0490 ± 0.0019 |
| 0.05092–0.05803 | 0.05421 | 1375 | 0.0459 ± 0.0026 | 1274 | 0.0511 ± 0.0031 | 2638 | 0.0484 ± 0.0019 |
| 0.05803–0.06925 | 0.0634 | 1427 | 0.0441 ± 0.0024 | 1288 | 0.0495 ± 0.0029 | 2707 | 0.0466 ± 0.0018 |
| 0.06925–0.2772 | 0.08923 | 1404 | 0.0463 ± 0.0026 | 1331 | 0.0496 ± 0.0033 | 2719 | 0.0480 ± 0.0021 |

### f_acc(strona, lumiBlock) — błąd tylko stat.

| lumiBlock | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C |
|---|---|---|---|---|---|---|---|
| 4–20 | 11.89 | 1021 | 0.0449 ± 0.0061 | 949 | 0.0531 ± 0.0074 | 1962 | 0.0488 ± 0.0048 |
| 20–221 | 147.1 | 1030 | 0.0486 ± 0.0069 | 1003 | 0.0488 ± 0.0075 | 2026 | 0.0486 ± 0.0050 |
| 221–238 | 229 | 1017 | 0.0420 ± 0.0057 | 996 | 0.0608 ± 0.0078 | 2007 | 0.0512 ± 0.0046 |
| 238–256 | 246.5 | 1072 | 0.0515 ± 0.0066 | 973 | 0.0468 ± 0.0052 | 2032 | 0.0493 ± 0.0042 |
| 256–273 | 264 | 1006 | 0.0361 ± 0.0051 | 999 | 0.0449 ± 0.0055 | 1999 | 0.0404 ± 0.0038 |
| 273–289 | 280.5 | 1035 | 0.0435 ± 0.0066 | 934 | 0.0473 ± 0.0056 | 1960 | 0.0453 ± 0.0045 |
| 289–306 | 297 | 1043 | 0.0447 ± 0.0064 | 948 | 0.0488 ± 0.0074 | 1987 | 0.0465 ± 0.0047 |
| 306–323 | 314.4 | 1047 | 0.0523 ± 0.0070 | 1061 | 0.0497 ± 0.0062 | 2092 | 0.0511 ± 0.0047 |

### Gatunki × pT: tag przypadkowy pod pikiem i w sidebandach (A∨C; kandydaci, bez fitu)

N_tag = kandydaci w oknie w zdarzeniach z tagiem po cięciu ξ_cal; N_acc = Σ p_acc po kandydatach. f_acc ± stat ± syst (syst wspólna, względna 28.2 %: pula ⊕ ekstrapolacja ⊕ ½Δ metod). Sidebandy NIE są próbką „bez prawdziwych tagów” (kombinatoryka masowa w zdarzeniach dyfrakcyjnych niesie prawdziwy proton) — stosunek pik/SB jest tylko informacyjny.

**Dplus** — Dmeson_mass: pik (1.84, 1.9), SB ((1.75, 1.8), (1.94, 2.0))

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik ± stat ± syst | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | v1 f_acc pik |
|---|---|---|---|---|---|---|---|---|
| 3–150 | 8797 | 446.2 ± 10.3 | 0.0507 ± 0.0015 ± 0.0143 | 0.0487 ± 0.0020 | 0.0528 ± 0.0022 | 0.0511 ± 0.0014 | 0.993 ± 0.040 | 0.0537 |
| 3–4 | 4045 | 221.6 ± 5.1 | 0.0548 ± 0.0017 ± 0.0154 | 0.0529 ± 0.0023 | 0.0567 ± 0.0026 | 0.0546 ± 0.0016 | 1.004 ± 0.043 | 0.0586 |
| 4–5 | 1818 | 90.6 ± 2.1 | 0.0498 ± 0.0019 ± 0.0140 | 0.0464 ± 0.0024 | 0.0535 ± 0.0029 | 0.0524 ± 0.0018 | 0.951 ± 0.048 | 0.0527 |
| 5–6 | 932 | 43.4 ± 1.0 | 0.0466 ± 0.0021 ± 0.0131 | 0.0443 ± 0.0027 | 0.0492 ± 0.0032 | 0.0469 ± 0.0018 | 0.994 ± 0.059 | 0.0490 |
| 6–8 | 857 | 40.8 ± 1.0 | 0.0476 ± 0.0022 ± 0.0134 | 0.0456 ± 0.0029 | 0.0496 ± 0.0034 | 0.0463 ± 0.0018 | 1.028 ± 0.063 | 0.0501 |
| 8–12 | 673 | 30.6 ± 0.8 | 0.0455 ± 0.0024 ± 0.0128 | 0.0455 ± 0.0033 | 0.0456 ± 0.0034 | 0.0440 ± 0.0020 | 1.034 ± 0.071 | 0.0467 |
| 12–150 | 472 | 19.2 ± 0.5 | 0.0407 ± 0.0026 ± 0.0115 | 0.0397 ± 0.0036 | 0.0416 ± 0.0038 | 0.0433 ± 0.0024 | 0.940 ± 0.080 | 0.0419 |

**D0** — Dmeson_mass: pik (1.835, 1.895), SB ((1.75, 1.8), (1.93, 2.0))

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik ± stat ± syst | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | v1 f_acc pik |
|---|---|---|---|---|---|---|---|---|
| 3–150 | 6471 | 324.8 ± 8.3 | 0.0502 ± 0.0015 ± 0.0141 | 0.0479 ± 0.0020 | 0.0528 ± 0.0023 | 0.0497 ± 0.0014 | 1.010 ± 0.042 | 0.0515 |
| 3–4 | 2789 | 143.0 ± 3.6 | 0.0513 ± 0.0017 ± 0.0144 | 0.0495 ± 0.0023 | 0.0533 ± 0.0026 | 0.0513 ± 0.0016 | 0.999 ± 0.045 | 0.0530 |
| 4–5 | 1492 | 74.5 ± 1.9 | 0.0499 ± 0.0019 ± 0.0141 | 0.0456 ± 0.0025 | 0.0550 ± 0.0032 | 0.0504 ± 0.0017 | 0.991 ± 0.051 | 0.0510 |
| 5–6 | 703 | 37.0 ± 1.0 | 0.0526 ± 0.0025 ± 0.0148 | 0.0521 ± 0.0035 | 0.0534 ± 0.0037 | 0.0490 ± 0.0019 | 1.075 ± 0.067 | 0.0539 |
| 6–8 | 649 | 33.1 ± 0.9 | 0.0510 ± 0.0025 ± 0.0144 | 0.0508 ± 0.0036 | 0.0512 ± 0.0036 | 0.0487 ± 0.0020 | 1.046 ± 0.067 | 0.0520 |
| 8–12 | 500 | 22.8 ± 0.6 | 0.0456 ± 0.0025 ± 0.0129 | 0.0407 ± 0.0031 | 0.0513 ± 0.0042 | 0.0467 ± 0.0021 | 0.977 ± 0.070 | 0.0461 |
| 12–150 | 338 | 14.4 ± 0.4 | 0.0427 ± 0.0027 ± 0.0120 | 0.0418 ± 0.0038 | 0.0433 ± 0.0039 | 0.0414 ± 0.0021 | 1.031 ± 0.085 | 0.0428 |

**Dstar** — Dmeson_delta_mass: pik (0.144, 0.147), SB ((0.15, 0.16),)

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik ± stat ± syst | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | v1 f_acc pik |
|---|---|---|---|---|---|---|---|---|
| 3–150 | 572 | 27.6 ± 0.7 | 0.0482 ± 0.0025 ± 0.0136 | 0.0458 ± 0.0033 | 0.0505 ± 0.0037 | 0.0524 ± 0.0018 | 0.920 ± 0.057 | 0.0507 |
| 3–4 | 226 | 11.4 ± 0.3 | 0.0506 ± 0.0037 ± 0.0142 | 0.0495 ± 0.0052 | 0.0520 ± 0.0054 | 0.0555 ± 0.0025 | 0.911 ± 0.079 | 0.0537 |
| 4–5 | 166 | 7.5 ± 0.2 | 0.0452 ± 0.0039 ± 0.0127 | 0.0407 ± 0.0048 | 0.0499 ± 0.0064 | 0.0530 ± 0.0026 | 0.853 ± 0.085 | 0.0471 |
| 5–6 | 71 | 3.5 ± 0.1 | 0.0495 ± 0.0063 ± 0.0139 | 0.0449 ± 0.0079 | 0.0559 ± 0.0102 | 0.0475 ± 0.0030 | 1.042 ± 0.147 | 0.0517 |
| 6–8 | 60 | 2.8 ± 0.1 | 0.0458 ± 0.0062 ± 0.0129 | 0.0505 ± 0.0100 | 0.0413 ± 0.0075 | 0.0500 ± 0.0036 | 0.917 ± 0.140 | 0.0480 |
| 8–12 | 26 | 1.5 ± 0.0 | 0.0567 ± 0.0117 ± 0.0160 | 0.0481 ± 0.0134 | 0.0631 ± 0.0185 | 0.0528 ± 0.0052 | 1.074 ± 0.246 | 0.0624 |
| 12–150 | 23 | 0.9 ± 0.0 | 0.0398 ± 0.0092 ± 0.0112 | 0.0395 ± 0.0132 | 0.0404 ± 0.0129 | 0.0396 ± 0.0051 | 1.006 ± 0.265 | 0.0385 |

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
