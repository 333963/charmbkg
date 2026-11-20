# MC_templates — WYNIKI: test szablonów tła z MC (zad. #7; projekt B)

Wygenerowane 2026-11-17 08:04 przez `mctmpl_analyze.py` z `charm_work/MC_templates/runs/test`. Pliki: data_428770 1, data_427929 4, mb_nd 2, fd_nd 1. Kryteria: `KRYTERIA.md` (zapisane przed wynikami). Selekcja i klasyfikator H są nasze (nie CharmNet). Dane: GRL, bez zdarzeń regionu sygnału głównego łańcucha (tag AFP ∧ npu = 0 ∧ ξ_cal).

## 1. Skład MC wg kategorii głównego łańcucha (TRExFitter) i stabilność szablonów (K1)

Kategorie = próbki TRExFitter z `HistogramManager.cxx:555-556` (Sig 0,1; Reflection = Sig z zamienionymi K/π, tylko D⁰; CharmMismatch 2,3; DecayModeMismatch 6,7; PartComb 8–11; Comb 4,5,12–15). Bez podziału Sig/SigNonFid (wymaga tagu prawdy AFP — poza zakresem). Udziały w pełnym oknie (D⁰ 1,70–2,05; D⁺ 1,70–2,20 GeV, po wecie D*), mb_nd bez wag; κ = N/N_Sig.

### D0

| pT | WP | N_Sig | κ Reflection | κ CharmMismatch | κ DecayModeMismatch | κ PartComb | κ Comb | Sig w oknie [%] |
|---|---|---|---|---|---|---|---|---|
| 3-5 | all | 124 | 0.621 ± 0.071 | 0.637 ± 0.072 | 1.677 ± 0.116 | 39.492 ± 0.564 | 699.855 ± 2.376 | 0.13 |
| 3-5 | 95 | 116 | 0.647 ± 0.075 | 0.345 ± 0.055 | 1.095 ± 0.097 | 16.948 ± 0.382 | 236.310 ± 1.427 | 0.39 |
| 3-5 | 80 | 98 | 0.663 ± 0.082 | 0.265 ± 0.052 | 0.908 ± 0.096 | 7.408 ± 0.275 | 69.684 ± 0.843 | 1.25 |
| 3-5 | 50 | 54 | 0.870 ± 0.127 | 0.296 ± 0.074 | 0.981 ± 0.135 | 2.815 ± 0.228 | 8.648 ± 0.400 | 6.84 |
| 5-8 | all | 60 | 0.650 ± 0.104 | 0.583 ± 0.099 | 1.533 ± 0.160 | 13.100 ± 0.467 | 166.950 ± 1.668 | 0.54 |
| 5-8 | 95 | 54 | 0.704 ± 0.114 | 0.444 ± 0.091 | 1.259 ± 0.153 | 5.389 ± 0.316 | 53.019 ± 0.991 | 1.62 |
| 5-8 | 80 | 51 | 0.725 ± 0.119 | 0.333 ± 0.081 | 1.000 ± 0.140 | 1.882 ± 0.192 | 14.255 ± 0.529 | 5.21 |
| 5-8 | 50 | 31 | 0.774 ± 0.158 | 0.323 ± 0.102 | 0.935 ± 0.174 | 0.968 ± 0.177 | 1.645 ± 0.230 | 17.71 |
| 8-30 | all | 17 | 0.706 ± 0.204 | 0.882 ± 0.228 | 1.765 ± 0.322 | 8.294 ± 0.698 | 76.529 ± 2.122 | 1.12 |
| 8-30 | 95 | 16 | 0.750 ± 0.217 | 0.750 ± 0.217 | 1.438 ± 0.300 | 4.062 ± 0.504 | 22.875 ± 1.196 | 3.24 |
| 8-30 | 80 | 14 | 0.857 ± 0.247 | 0.714 ± 0.226 | 1.429 ± 0.319 | 2.071 ± 0.385 | 5.643 ± 0.635 | 8.54 |
| 8-30 | 50 | 8 | 0.875 ± 0.331 | 0.625 ± 0.280 | 1.375 ± 0.415 | 1.000 ± 0.354 | 1.000 ± 0.354 | 17.02 |

Podskład (WP95, pT 3–30, zdarzenia w oknie): CM_411 55, CM_431 18, CM_4122 2, CM_other 1, DMM_C 185, DMM_B 33, PC_411 655, PC_421 1269, PC_413 14, PC_431 299, PC_4122 71, PC_other 14.

### Dplus

| pT | WP | N_Sig | κ CharmMismatch | κ DecayModeMismatch | κ PartComb | κ Comb | Sig w oknie [%] |
|---|---|---|---|---|---|---|---|
| 3-5 | all | 133 | 1.053 ± 0.089 | 0.346 ± 0.051 | 146.744 ± 1.050 | 1703.714 ± 3.579 | 0.05 |
| 3-5 | 95 | 128 | 0.539 ± 0.065 | 0.250 ± 0.044 | 39.750 ± 0.557 | 213.148 ± 1.290 | 0.39 |
| 3-5 | 80 | 106 | 0.396 ± 0.061 | 0.236 ± 0.047 | 14.302 ± 0.367 | 31.425 ± 0.544 | 2.11 |
| 3-5 | 50 | 68 | 0.191 ± 0.053 | 0.279 ± 0.064 | 4.279 ± 0.251 | 2.324 ± 0.185 | 12.39 |
| 5-8 | all | 66 | 1.167 ± 0.133 | 0.364 ± 0.074 | 45.864 ± 0.834 | 286.455 ± 2.083 | 0.30 |
| 5-8 | 95 | 64 | 0.703 ± 0.105 | 0.219 ± 0.058 | 15.219 ± 0.488 | 29.672 ± 0.681 | 2.14 |
| 5-8 | 80 | 53 | 0.434 ± 0.090 | 0.208 ± 0.063 | 6.170 ± 0.341 | 4.226 ± 0.282 | 8.31 |
| 5-8 | 50 | 33 | 0.364 ± 0.105 | 0.212 ± 0.080 | 2.242 ± 0.261 | 0.515 ± 0.125 | 23.08 |
| 8-30 | all | 39 | 0.821 ± 0.145 | 0.154 ± 0.063 | 17.333 ± 0.667 | 66.103 ± 1.302 | 1.17 |
| 8-30 | 95 | 36 | 0.528 ± 0.121 | 0.167 ± 0.068 | 6.000 ± 0.408 | 8.417 ± 0.484 | 6.21 |
| 8-30 | 80 | 28 | 0.214 ± 0.087 | 0.036 ± 0.036 | 2.821 ± 0.317 | 1.286 ± 0.214 | 18.67 |
| 8-30 | 50 | 20 | 0.150 ± 0.087 | 0.000 ± 0.000 | 0.550 ± 0.166 | 0.100 ± 0.071 | 55.56 |

Podskład (WP95, pT 3–30, zdarzenia w oknie): CM_421 22, CM_431 56, CM_413 16, CM_423 34, CM_4122 4, CM_other 1, DMM_C 44, DMM_B 8, PC_411 1983, PC_421 1155, PC_423 1256, PC_413 1012, PC_431 716, PC_4122 134, PC_other 22.

### Stabilność przy przeważeniu MC do aktywności klasy trygera (K1)

Wagi w(nSelBL, pT toru wiodącego) = dane/MC w pasmach bocznych (|m − m₀| > 75 MeV, bez cięcia H), per gatunek i bin pT. Δκ/κ i odległość KS kształtu (szablon przeważony vs nieprzeważony), WP80. n_eff/N < 0,2 → „nieokreślony”.

| próbka danych | gatunek | pT | n_eff/N | dane poza MC [%] | Reflection: Δκ/κ [%] (KS) | CharmMismatch: Δκ/κ [%] (KS) | DecayModeMismatch: Δκ/κ [%] (KS) | PartComb: Δκ/κ [%] (KS) | Comb: Δκ/κ [%] (KS) | werdykt |
|---|---|---|---|---|---|---|---|---|---|---|
| 427929:all | D0 | 3-5 | 0.95 | 0.0 | -0.2 (0.018) | +11.7 (0.057) | -3.5 (0.011) | -1.3 (0.007) | -0.5 (0.004) | stabilne |
| 427929:all | D0 | 5-8 | 0.87 | 0.4 | -2.3 (0.061) | -4.0 (0.073) | +12.4 (0.054) | +18.6 (0.030) | +18.1 (0.007) | stabilne |
| 427929:all | D0 | 8-30 | 0.72 | 0.7 | -6.5 (0.167) | +7.6 (0.059) | +93.0 (0.259)! | +53.0 (0.024) | +64.6 (0.110)! | NIESTABILNE: DecayModeMismatch, Comb |
| 427929:all | Dplus | 3-5 | 0.90 | 0.0 | — | +0.2 (0.023) | +3.2 (0.020) | +4.6 (0.004) | +8.1 (0.003) | stabilne |
| 427929:all | Dplus | 5-8 | 0.78 | 0.0 | — | +0.5 (0.019) | -5.0 (0.029) | +7.6 (0.012) | +9.5 (0.014) | stabilne |
| 427929:all | Dplus | 8-30 | 0.75 | 0.1 | — | — | — | +8.8 (0.038) | +17.2 (0.081) | stabilne |
| 427929:MB | D0 | 3-5 | 0.95 | 0.0 | -0.2 (0.018) | +11.7 (0.057) | -3.5 (0.011) | -1.3 (0.007) | -0.5 (0.004) | stabilne |
| 427929:MB | D0 | 5-8 | 0.87 | 0.4 | -2.3 (0.061) | -4.0 (0.073) | +12.4 (0.054) | +18.6 (0.030) | +18.1 (0.007) | stabilne |
| 427929:MB | D0 | 8-30 | 0.72 | 0.7 | -6.5 (0.167) | +7.6 (0.059) | +93.0 (0.259)! | +53.0 (0.024) | +64.6 (0.110)! | NIESTABILNE: DecayModeMismatch, Comb |
| 427929:MB | Dplus | 3-5 | 0.90 | 0.0 | — | +0.2 (0.023) | +3.2 (0.020) | +4.7 (0.004) | +8.2 (0.003) | stabilne |
| 427929:MB | Dplus | 5-8 | 0.78 | 0.0 | — | +0.5 (0.019) | -5.0 (0.029) | +7.6 (0.012) | +9.5 (0.014) | stabilne |
| 427929:MB | Dplus | 8-30 | 0.75 | 0.1 | — | — | — | +8.8 (0.038) | +17.2 (0.081) | stabilne |
| 428770:all | D0 | 3-5 | 0.14 | 0.0 | -7.1 (0.154) | +179.5 (0.287)! | -7.2 (0.147) | +139.1 (0.059)! | +93.6 (0.024)! | nieokreślony (n_eff/N < 0,2); NIESTABILNE: CharmMismatch, PartComb, Comb |
| 428770:all | D0 | 5-8 | 0.25 | 0.1 | -4.0 (0.140) | +105.8 (0.314)! | -11.2 (0.135) | +74.0 (0.175)! | +55.0 (0.051)! | NIESTABILNE: CharmMismatch, PartComb, Comb |
| 428770:all | D0 | 8-30 | 0.49 | 0.4 | -14.9 (0.054) | +10.3 (0.174) | +71.5 (0.076)! | +76.3 (0.088)! | +59.6 (0.151)! | NIESTABILNE: DecayModeMismatch, PartComb, Comb |
| 428770:all | Dplus | 3-5 | 0.14 | 0.0 | — | +117.4 (0.369)! | +28.5 (0.245) | +131.8 (0.051)! | +149.0 (0.043)! | nieokreślony (n_eff/N < 0,2); NIESTABILNE: CharmMismatch, PartComb, Comb |
| 428770:all | Dplus | 5-8 | 0.22 | 0.0 | — | +62.9 (0.206)! | -27.6 (0.303) | +3.6 (0.057) | +16.6 (0.059) | NIESTABILNE: CharmMismatch |
| 428770:all | Dplus | 8-30 | 0.56 | 0.0 | — | — | — | -2.0 (0.057) | +20.2 (0.226) | stabilne |
| 428770:J12only | D0 | 3-5 | 0.11 | 0.0 | -8.2 (0.218) | +211.7 (0.306)! | -3.7 (0.183) | +184.5 (0.074)! | +118.8 (0.026)! | nieokreślony (n_eff/N < 0,2); NIESTABILNE: CharmMismatch, PartComb, Comb |
| 428770:J12only | D0 | 5-8 | 0.22 | 0.0 | -4.3 (0.169) | +131.2 (0.354)! | -16.0 (0.155) | +84.7 (0.193)! | +58.7 (0.061)! | NIESTABILNE: CharmMismatch, PartComb, Comb |
| 428770:J12only | D0 | 8-30 | 0.46 | 0.3 | -16.0 (0.055) | +13.9 (0.179) | +80.1 (0.079)! | +82.8 (0.099)! | +66.1 (0.161)! | NIESTABILNE: DecayModeMismatch, PartComb, Comb |
| 428770:J12only | Dplus | 3-5 | 0.11 | 0.0 | — | +158.8 (0.427)! | +36.5 (0.288) | +159.2 (0.060)! | +174.2 (0.051)! | nieokreślony (n_eff/N < 0,2); NIESTABILNE: CharmMismatch, PartComb, Comb |
| 428770:J12only | Dplus | 5-8 | 0.20 | 0.0 | — | +70.9 (0.218)! | -30.1 (0.339) | +0.3 (0.067) | +14.8 (0.066) | NIESTABILNE: CharmMismatch |
| 428770:J12only | Dplus | 8-30 | 0.53 | 0.0 | — | — | — | -3.2 (0.060) | +20.7 (0.246) | stabilne |
| 428770:MBTSbased | D0 | 3-5 | 0.41 | 0.0 | -4.1 (0.048) | +116.1 (0.241)! | -16.0 (0.093) | +47.2 (0.024)! | +42.1 (0.019)! | NIESTABILNE: CharmMismatch, PartComb, Comb |
| 428770:MBTSbased | D0 | 5-8 | 0.62 | 0.1 | -6.3 (0.082) | +19.3 (0.093) | +11.0 (0.085) | +44.4 (0.093)! | +40.2 (0.014)! | NIESTABILNE: PartComb, Comb |
| 428770:MBTSbased | D0 | 8-30 | 0.72 | 1.2 | -3.1 (0.193) | -22.6 (0.116) | -0.6 (0.156) | +22.0 (0.058) | +4.9 (0.085) | stabilne |
| 428770:MBTSbased | Dplus | 3-5 | 0.39 | 0.0 | — | +25.0 (0.126) | +6.7 (0.144) | +78.3 (0.027)! | +104.0 (0.016)! | NIESTABILNE: PartComb, Comb |
| 428770:MBTSbased | Dplus | 5-8 | 0.57 | 0.0 | — | +18.6 (0.108) | -14.9 (0.160) | +23.9 (0.027) | +32.0 (0.027)! | NIESTABILNE: Comb |
| 428770:MBTSbased | Dplus | 8-30 | 0.85 | 0.1 | — | — | — | +5.7 (0.032) | +17.2 (0.079) | stabilne |
| 428770:AFPstream | D0 | 3-5 | 0.29 | 0.0 | -5.4 (0.134) | +192.4 (0.304)! | -22.0 (0.106) | +43.3 (0.053)! | +27.2 (0.026)! | NIESTABILNE: CharmMismatch, PartComb, Comb |
| 428770:AFPstream | D0 | 5-8 | 0.60 | 0.0 | -11.0 (0.103) | +12.5 (0.146) | +9.6 (0.157) | +36.9 (0.040)! | +28.2 (0.020) | NIESTABILNE: PartComb |
| 428770:AFPstream | D0 | 8-30 | 0.49 | 1.7 | +0.3 (0.247) | -54.0 (0.400) | +27.9 (0.075) | +49.4 (0.059) | +22.2 (0.084) | stabilne |
| 428770:AFPstream | Dplus | 3-5 | 0.26 | 0.0 | — | -8.9 (0.204) | +0.5 (0.202) | +69.7 (0.040)! | +102.1 (0.019)! | NIESTABILNE: PartComb, Comb |
| 428770:AFPstream | Dplus | 5-8 | 0.44 | 0.0 | — | +26.6 (0.158) | +21.5 (0.217) | +49.8 (0.060)! | +49.4 (0.037)! | NIESTABILNE: PartComb, Comb |
| 428770:AFPstream | Dplus | 8-30 | 0.51 | 0.0 | — | — | — | +35.5 (0.082) | +43.3 (0.115) | stabilne |

## 2. Testy w obszarach wzbogaconych (K2)

### 2a. Refleksja K↔π D⁰ z tagu D* (bliźniak o zamienionej hipotezie)

Wiersze D* (znak pionu miękkiego ustala hipotezę) dopasowane po torach do wierszy D⁰: R = bliźniak (zamieniona hipoteza) przechodzi WP, S = poprawna hipoteza przechodzi WP; tło pod pikiem Δm odjęte (funkcja progowa dopasowana poza pikiem). MC: ta sama procedura na mb_nd + fd_nd; „prawda” = N(Reflection)/N(Sig) z kategorii. r = (R/S)_dane / (R/S)_MC.

| próbka danych | pT | WP | N_S(D*) dane | (R/S) dane | (R/S) MC tag | (R/S) MC prawda | r = dane/MC | K2 | χ²/ndf kształtu (p) |
|---|---|---|---|---|---|---|---|---|---|
| 427929:all | 3-5 | 95 | 6 ± 5 | 0.383 ± 0.575 | 0.646 ± 0.034 | 0.647 | 0.593 ± 0.891 | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 3-5 | 80 | 6 ± 3 | 0.580 ± 0.543 | 0.701 ± 0.037 | 0.702 | 0.827 ± 0.776 | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 3-5 | 50 | — | — | 0.818 ± 0.052 | 0.798 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 5-8 | 95 | — | — | 0.565 ± 0.039 | 0.619 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 5-8 | 80 | — | — | 0.656 ± 0.047 | 0.690 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 5-8 | 50 | — | — | 0.763 ± 0.067 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 8-30 | 95 | — | — | 0.557 ± 0.056 | 0.618 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 8-30 | 80 | — | — | 0.618 ± 0.065 | 0.668 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 8-30 | 50 | — | — | — | 0.725 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 3-30 | 95 | 9 ± 5 | 0.395 ± 0.436 | 0.608 ± 0.023 | 0.635 | 0.649 ± 0.719 | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 3-30 | 80 | 9 ± 4 | 0.539 ± 0.388 | 0.675 ± 0.027 | 0.694 | 0.798 ± 0.575 | niemierzalne (δ = 0,5 a priori) | — |
| 427929:all | 3-30 | 50 | — | — | 0.772 ± 0.037 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 3-5 | 95 | 6 ± 5 | 0.383 ± 0.575 | 0.646 ± 0.034 | 0.647 | 0.593 ± 0.891 | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 3-5 | 80 | 6 ± 3 | 0.580 ± 0.543 | 0.701 ± 0.037 | 0.702 | 0.827 ± 0.776 | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 3-5 | 50 | — | — | 0.818 ± 0.052 | 0.798 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 5-8 | 95 | — | — | 0.565 ± 0.039 | 0.619 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 5-8 | 80 | — | — | 0.656 ± 0.047 | 0.690 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 5-8 | 50 | — | — | 0.763 ± 0.067 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 8-30 | 95 | — | — | 0.557 ± 0.056 | 0.618 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 8-30 | 80 | — | — | 0.618 ± 0.065 | 0.668 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 8-30 | 50 | — | — | — | 0.725 | — | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 3-30 | 95 | 9 ± 5 | 0.395 ± 0.436 | 0.608 ± 0.023 | 0.635 | 0.649 ± 0.719 | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 3-30 | 80 | 9 ± 4 | 0.539 ± 0.388 | 0.675 ± 0.027 | 0.694 | 0.798 ± 0.575 | niemierzalne (δ = 0,5 a priori) | — |
| 427929:MB | 3-30 | 50 | — | — | 0.772 ± 0.037 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 3-5 | 95 | -5 ± 12 | — | 0.646 ± 0.034 | 0.647 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 3-5 | 80 | 11 ± 7 | 0.404 ± 0.549 | 0.701 ± 0.037 | 0.702 | 0.577 ± 0.784 | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 3-5 | 50 | -1 ± 2 | — | 0.818 ± 0.052 | 0.798 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 5-8 | 95 | 1 ± 6 | 1.420 ± 12.731 | 0.565 ± 0.039 | 0.619 | 2.512 ± 22.529 | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 5-8 | 80 | 7 ± 4 | 0.220 ± 0.455 | 0.656 ± 0.047 | 0.690 | 0.335 ± 0.694 | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 5-8 | 50 | — | — | 0.763 ± 0.067 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 8-30 | 95 | 16 ± 5 | — | 0.557 ± 0.056 | 0.618 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 8-30 | 80 | — | — | 0.618 ± 0.065 | 0.668 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 8-30 | 50 | — | — | — | 0.725 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 3-30 | 95 | 12 ± 14 | 0.390 ± 0.859 | 0.608 ± 0.023 | 0.635 | 0.641 ± 1.414 | niemierzalne (δ = 0,5 a priori) | — |
| 428770:all | 3-30 | 80 | 35 ± 9 | 0.409 ± 0.219 | 0.675 ± 0.027 | 0.694 | 0.606 ± 0.325 | zgodne | — |
| 428770:all | 3-30 | 50 | 17 ± 5 | 0.674 ± 0.301 | 0.772 ± 0.037 | 0.784 | 0.873 ± 0.392 | zgodne | — |
| 428770:J12only | 3-5 | 95 | -15 ± 11 | — | 0.646 ± 0.034 | 0.647 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 3-5 | 80 | 5 ± 6 | 0.308 ± 1.030 | 0.701 ± 0.037 | 0.702 | 0.440 ± 1.469 | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 3-5 | 50 | -2 ± 1 | — | 0.818 ± 0.052 | 0.798 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 5-8 | 95 | -3 ± 5 | — | 0.565 ± 0.039 | 0.619 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 5-8 | 80 | 5 ± 4 | — | 0.656 ± 0.047 | 0.690 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 5-8 | 50 | — | — | 0.763 ± 0.067 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 8-30 | 95 | 14 ± 4 | — | 0.557 ± 0.056 | 0.618 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 8-30 | 80 | — | — | 0.618 ± 0.065 | 0.668 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 8-30 | 50 | — | — | — | 0.725 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 3-30 | 95 | -4 ± 12 | — | 0.608 ± 0.023 | 0.635 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:J12only | 3-30 | 80 | 23 ± 8 | 0.269 ± 0.263 | 0.675 ± 0.027 | 0.694 | 0.399 ± 0.390 | zgodne | — |
| 428770:J12only | 3-30 | 50 | 10 ± 4 | 0.441 ± 0.346 | 0.772 ± 0.037 | 0.784 | 0.572 ± 0.449 | zgodne | — |
| 428770:MBTSbased | 3-5 | 95 | 0 ± 4 | — | 0.646 ± 0.034 | 0.647 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 3-5 | 80 | 1 ± 2 | — | 0.701 ± 0.037 | 0.702 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 3-5 | 50 | — | — | 0.818 ± 0.052 | 0.798 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 5-8 | 95 | 4 ± 3 | — | 0.565 ± 0.039 | 0.619 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 5-8 | 80 | — | — | 0.656 ± 0.047 | 0.690 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 5-8 | 50 | — | — | 0.763 ± 0.067 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 8-30 | 95 | — | — | 0.557 ± 0.056 | 0.618 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 8-30 | 80 | — | — | 0.618 ± 0.065 | 0.668 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 8-30 | 50 | — | — | — | 0.725 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 3-30 | 95 | 7 ± 6 | 0.235 ± 0.496 | 0.608 ± 0.023 | 0.635 | 0.386 ± 0.817 | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 3-30 | 80 | 6 ± 4 | 0.710 ± 0.581 | 0.675 ± 0.027 | 0.694 | 1.051 ± 0.861 | niemierzalne (δ = 0,5 a priori) | — |
| 428770:MBTSbased | 3-30 | 50 | — | — | 0.772 ± 0.037 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 3-5 | 95 | — | — | 0.646 ± 0.034 | 0.647 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 3-5 | 80 | — | — | 0.701 ± 0.037 | 0.702 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 3-5 | 50 | — | — | 0.818 ± 0.052 | 0.798 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 5-8 | 95 | — | — | 0.565 ± 0.039 | 0.619 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 5-8 | 80 | — | — | 0.656 ± 0.047 | 0.690 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 5-8 | 50 | — | — | 0.763 ± 0.067 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 8-30 | 95 | — | — | 0.557 ± 0.056 | 0.618 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 8-30 | 80 | — | — | 0.618 ± 0.065 | 0.668 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 8-30 | 50 | — | — | — | 0.725 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 3-30 | 95 | — | — | 0.608 ± 0.023 | 0.635 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 3-30 | 80 | — | — | 0.675 ± 0.027 | 0.694 | — | niemierzalne (δ = 0,5 a priori) | — |
| 428770:AFPstream | 3-30 | 50 | — | — | 0.772 ± 0.037 | 0.784 | — | niemierzalne (δ = 0,5 a priori) | — |

## 3. Fity masy: konfiguracja „jak w głównym łańcuchu” vs warianty (K3, K4)

Model binowany (10 MeV, okna jak w ficie TRExFitter głównego łańcucha): λ = N_S·(s + Σ_k μ_k κ_k t_k) + N_C·c. **C0** = „jak w głównym łańcuchu”: N_S i Comb swobodne, kształty wszystkich składowych i κ_k z MC, μ_k = 1 (w głównym łańcuchu normalizacja absolutna σ·L; tu względem sygnału, bo nie mamy L_eff per tryger — równoważne tylko przy μ_sygnału = 1). C0w = C0 z MC przeważonym do aktywności klasy; C1 = wolna μ refleksji (D⁰) / CharmMismatch (D⁺); C2 = wolna μ rozpadów częściowo zrekonstruowanych (D⁰: CharmMismatch + DecayModeMismatch wspólnie; D⁺: DecayModeMismatch); C3 = wolna μ PartComb; C4 = refleksja (kształt i R/S) z danych z tagu D*; C5 = wszystkie μ wolne; C6 = Comb jako exp(pol2) zamiast kształtu MC; C7 = C0 + swobodne przesunięcie i szerokość sygnału. Szablony = surowe histogramy MC (jak w głównym łańcuchu), bez niepewności statystycznej MC w likelihoodzie (uproszczenie; w głównym łańcuchu gammy); składowa z < 20 zdarzeń MC jest dołączana do Comb, a Comb z < 300 zdarzeń MC wygładzany exp(pol2). PartComb nie jest wygładzany (ma część pikującą).

**K4 (domknięcie na mb_nd):** 3/3 komórek w 1σ — OK. Szczegóły w results.json (`closure`).

### D0: ΔN_S/N_S [%] względem C0

| próbka | pT | WP | N_S ± σ (C0) | χ²/ndf C0 | C0w | C1 | C2 | C3 | C4 | C5 | C6 | C7 | μ̂ C1 | μ̂ C2 | μ̂ C3 | σ(C1)/σ(C0) | σ(C2)/σ(C0) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 427929:all | 3-5 | 95 | 0 ± 39 | 41/33 | nieokr. | | | | | | | | | | | | | |
| 427929:all | 3-5 | 80 | 26 ± 22 | 53/33 | nieokr. | | | | | | | | | | | | | |
| 427929:all | 3-5 | 50 | 18 ± 5 | 21/23 | -13.3 | -39.8 | +13.6 | +27.4 | — | -0.3 | -86.6 | +9.5 | 8.93 ± 6.70 | 0.00 ± 0.40 | 0.54 ± 0.45 | 1.23 | 1.04 |
| 427929:all | 5-8 | 95 | 0 ± 26 | 35/33 | nieokr. | | | | | | | | | | | | | |
| 427929:all | 5-8 | 80 | 4 ± 7 | 35/23 | nieokr. | | | | | | | | | | | | | |
| 427929:MB | 3-5 | 95 | 0 ± 39 | 41/33 | nieokr. | | | | | | | | | | | | | |
| 427929:MB | 3-5 | 80 | 26 ± 22 | 53/33 | nieokr. | | | | | | | | | | | | | |
| 427929:MB | 3-5 | 50 | 18 ± 5 | 21/23 | -13.3 | -39.8 | +13.6 | +27.4 | — | -0.3 | -86.6 | +9.5 | 8.93 ± 6.70 | 0.00 ± 0.40 | 0.54 ± 0.45 | 1.23 | 1.04 |
| 427929:MB | 5-8 | 95 | 0 ± 26 | 35/33 | nieokr. | | | | | | | | | | | | | |
| 427929:MB | 5-8 | 80 | 4 ± 7 | 35/23 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 3-5 | 95 | 35 ± 106 | 129/33 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 3-5 | 80 | 130 ± 46 | 118/33 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 3-5 | 50 | 91 ± 12 | 97/33 | -44.7 | -29.4 | +5.4 | +12.6 | — | -52.5 | -74.4 | +3.6 | 6.57 ± 2.17 | 0.00 ± 0.96 | 0.75 ± 0.25 | 1.18 | 1.03 |
| 428770:all | 5-8 | 95 | 111 ± 44 | 84/33 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 5-8 | 80 | 57 ± 23 | 157/33 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 5-8 | 50 | 33 ± 8 | 30/24 | -5.6 | -66.4 | +10.5 | +6.3 | — | -56.0 | -47.1 | -4.3 | 17.02 ± 18.51 | 0.00 ± 0.67 | 0.74 ± 0.55 | 1.15 | 0.96 |
| 428770:all | 8-30 | 95 | 225 ± 14 | 432/33 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 8-30 | 80 | 41 ± 9 | 129/33 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 8-30 | 50 | 77 ± 12 | 12/23 | +34.8 | -30.0 | — | — | — | -30.0 | -22.5 | +8.6? | 3.01 ± 1.13 | — | — | 1.14 | — |
| 428770:J12only | 3-5 | 95 | 33 ± 85 | 116/33 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 3-5 | 80 | 88 ± 40 | 109/33 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 3-5 | 50 | 69 ± 10 | 96/33 | -46.9 | -28.4 | +4.5 | +12.2 | — | -57.4 | -75.9 | -4.6 | 6.36 ± 2.43 | 0.03 ± 1.73 | 0.76 ± 0.30 | 1.18 | 1.04 |
| 428770:J12only | 5-8 | 95 | 76 ± 40 | 71/33 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 5-8 | 80 | 47 ± 21 | 138/33 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 5-8 | 50 | 27 ± 7 | 27/20 | -8.4 | -65.4 | +9.3 | +10.1 | — | -58.0 | -49.8 | -2.4? | 17.37 ± 18.57 | 0.00 ± 3.43 | 0.59 ± 0.55 | 1.16 | 0.98 |
| 428770:J12only | 8-30 | 95 | 203 ± 14 | 377/33 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 8-30 | 80 | 39 ± 9 | 115/33 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 8-30 | 50 | 69 ± 11 | 13/22 | +36.2 | -28.8 | — | — | — | -28.8 | -18.0 | +12.9 | 2.88 ± 1.15 | — | — | 1.15 | — |
| 428770:MBTSbased | 3-5 | 95 | 0 ± 29 | 36/33 | nieokr. | | | | | | | | | | | | | |
| 428770:MBTSbased | 3-5 | 80 | 35 ± 18 | 55/33 | nieokr. | | | | | | | | | | | | | |
| 428770:MBTSbased | 5-8 | 95 | 24 ± 15 | 42/33 | nieokr. | | | | | | | | | | | | | |
| 428770:MBTSbased | 5-8 | 80 | 9 ± 9 | 48/27 | nieokr. | | | | | | | | | | | | | |
| 428770:MBTSbased | 8-30 | 95 | 17 ± 3 | 48/24 | -21.1 | -40.0 | +15.1 | +24.3 | — | -15.8 | -100.0 | +0.0? | 20.00 ± 4.44 | 0.39 ± 0.43 | 0.67 ± 0.36 | 0.49 | 1.22 |
| 428770:AFPstream | 3-5 | 95 | 1 ± 41 | 46/33 | nieokr. | | | | | | | | | | | | | |

### Dplus: ΔN_S/N_S [%] względem C0

| próbka | pT | WP | N_S ± σ (C0) | χ²/ndf C0 | C0w | C1 | C2 | C3 | C4 | C5 | C6 | C7 | μ̂ C1 | μ̂ C2 | μ̂ C3 | σ(C1)/σ(C0) | σ(C2)/σ(C0) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 427929:all | 3-5 | 95 | 35 ± 18 | 50/48 | nieokr. | | | | | | | | | | | | | |
| 427929:all | 3-5 | 80 | 16 ± 11 | 61/48 | nieokr. | | | | | | | | | | | | | |
| 427929:all | 5-8 | 95 | 12 ± 6 | 81/48 | nieokr. | | | | | | | | | | | | | |
| 427929:MB | 3-5 | 95 | 35 ± 18 | 50/48 | nieokr. | | | | | | | | | | | | | |
| 427929:MB | 3-5 | 80 | 16 ± 11 | 61/48 | nieokr. | | | | | | | | | | | | | |
| 427929:MB | 5-8 | 95 | 12 ± 6 | 81/48 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 3-5 | 95 | 299 ± 39 | 125/48 | -74.6 | +0.8 | -0.4 | -100.0 | — | -100.0? | -100.0 | +3.4? | 0.00 ± 0.37 | 0.00 ± 0.55 | 0.00 ± 19.79 | 0.99 | 0.98 |
| 428770:all | 3-5 | 80 | 143 ± 26 | 121/48 | -46.7 | -0.1 | -1.4 | -4.8 | — | -5.2 | -84.0 | -4.2 | 1.12 ± 1.60 | 0.00 ± 1.97 | 1.07 ± 0.34 | 1.00 | 0.97 |
| 428770:all | 3-5 | 50 | 23 ± 9 | 42/46 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 5-8 | 95 | 118 ± 20 | 188/48 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 5-8 | 80 | 27 ± 10 | 56/48 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 5-8 | 50 | 23 ± 5 | 15/15 | -67.5 | — | — | +20.0 | — | +20.0 | -24.1 | +12.2 | — | — | 0.69 ± 0.37 | — | — |
| 428770:all | 8-30 | 95 | 246 ± 15 | 582/48 | nieokr. | | | | | | | | | | | | | |
| 428770:all | 8-30 | 80 | 37 ± 11 | 92/46 | +21.5 | — | — | -17.6 | — | -17.6 | +9.7 | +5.7 | — | — | 1.31 ± 0.84 | — | — |
| 428770:all | 8-30 | 50 | 66 ± 10 | 11/7 | -21.5 | — | — | — | — | +0.0 | -8.9 | +2.5 | — | — | — | — | — |
| 428770:J12only | 3-5 | 95 | 215 ± 35 | 113/48 | -76.0 | +1.0 | -0.4 | -100.0 | — | -76.3? | -100.0 | +4.7 | 0.00 ± 0.47 | 0.00 ± 0.72 | 0.00 ± 18.84 | 0.98 | 0.98 |
| 428770:J12only | 3-5 | 80 | 96 ± 23 | 104/48 | -45.8 | -0.1 | -0.2 | -0.6 | — | -1.7 | -95.6 | -0.8 | 1.07 ± 2.07 | 0.90 ± 2.47 | 1.01 ± 0.41 | 1.00 | 1.02 |
| 428770:J12only | 3-5 | 50 | 20 ± 8 | 30/38 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 5-8 | 95 | 80 ± 19 | 180/48 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 5-8 | 80 | 17 ± 9 | 56/48 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 8-30 | 95 | 222 ± 14 | 544/48 | nieokr. | | | | | | | | | | | | | |
| 428770:J12only | 8-30 | 80 | 31 ± 10 | 89/46 | +26.3 | — | — | -35.5 | — | -35.5 | +9.1 | +5.5 | — | — | 1.78 ± 1.53 | — | — |
| 428770:MBTSbased | 3-5 | 95 | 68 ± 15 | 84/48 | -73.9 | -1.0 | -1.2 | -77.9? | — | -78.5? | -49.0 | +2.4 | 0.00 ± 3.53 | 0.00 ± 3.19 | 4.69 ± 2.03 | 0.98 | 0.98 |
| 428770:MBTSbased | 3-5 | 80 | 32 ± 10 | 74/48 | -40.9 | +0.5 | +4.6 | -27.3 | — | -7.3 | -48.7 | +5.8 | 0.00 ± 12.85 | 0.00 ± 0.68 | 1.54 ± 1.08 | 1.03 | 0.98 |
| 428770:MBTSbased | 5-8 | 95 | 24 ± 7 | 47/48 | -14.0 | +0.0 | — | -55.3 | — | -55.6 | -34.6 | +11.3 | 1.08 ± 1.78 | — | 2.52 ± 3.74 | 1.00 | — |
| 428770:MBTSbased | 8-30 | 95 | 16 ± 4 | 48/27 | +16.3 | — | — | -83.5 | — | -83.5 | -65.7 | +8.3? | — | — | 7.26 ± 11.08 | — | — |
| 428770:AFPstream | 3-5 | 95 | 1 ± 4 | 48/48 | nieokr. | | | | | | | | | | | | | |

### 2b. Rozpady częściowo zrekonstruowane / inne mody (dolne pasmo masy) — μ̂ z fitu C2

μ̂ = normalizacja danych względem MC (przy κ z MC). D⁰: CharmMismatch + DecayModeMismatch (Kππ⁰, KK, ππ …); D⁺: DecayModeMismatch. To ten sam fit co w §3 — nie jest to obszar niezależny.

| próbka | gatunek | pT | WP | μ̂ ± σ | K2 |
|---|---|---|---|---|---|
| 427929:all | D0 | 3-5 | 50 | 0.00 ± 0.40 | NIEZGODNE — szablon wymaga korekty |
| 427929:MB | D0 | 3-5 | 50 | 0.00 ± 0.40 | NIEZGODNE — szablon wymaga korekty |
| 428770:all | D0 | 3-5 | 50 | 0.00 ± 0.96 | niemierzalne (δ = 0,5 a priori) |
| 428770:all | D0 | 5-8 | 50 | 0.00 ± 0.67 | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | D0 | 3-5 | 50 | 0.03 ± 1.73 | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | D0 | 5-8 | 50 | 0.00 ± 3.43 | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | D0 | 8-30 | 95 | 0.39 ± 0.43 | zgodne |
| 428770:all | Dplus | 3-5 | 95 | 0.00 ± 0.55 | niemierzalne (δ = 0,5 a priori) |
| 428770:all | Dplus | 3-5 | 80 | 0.00 ± 1.97 | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | Dplus | 3-5 | 95 | 0.00 ± 0.72 | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | Dplus | 3-5 | 80 | 0.90 ± 2.47 | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | Dplus | 3-5 | 95 | 0.00 ± 3.19 | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | Dplus | 3-5 | 80 | 0.00 ± 0.68 | niemierzalne (δ = 0,5 a priori) |

### 2c. D⁺: piki refleksji w danych vs MC (D*⁺→D⁰π⁺, φ→KK z D_s⁺/D⁺→KKπ, Λ_c→pKπ)

Liczby nad liniowym tłem, znormalizowane do N_S(D⁺) (dane: fit C0 po wecie D*; MC: prawda). D*: pik w Δm_min = m(Kππ) − m(Kπ) (przed wetem); φ: pik w m(KK) po przypisaniu masy K pionowi; D_s i D⁺→KKπ: piki w m(KKπ) w oknie φ; Λ_c: pik w m(pKπ). pT 3–30 łącznie. r = dane/MC.

| próbka | WP | składowa | N/N_S dane | N/N_S MC | r = dane/MC | pozycja dane [MeV] | pozycja MC [MeV] | K2 |
|---|---|---|---|---|---|---|---|---|
| 428770:all | 95 | D* (Δm) | 0.0290 ± 0.0430 | 0.0615 ± 0.0660 | 0.47 ± 0.86 | 146.5 | 146.9 | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 95 | φ (m_KK) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 95 | D_s→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 95 | D⁺→φπ (m_KKπ) | 0.0033 ± 0.0034 | — | — | 1872.5 | — | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 95 | Λ_c (m_pKπ) | 0.1442 ± 0.2482 | 0.1350 ± 0.4263 | 1.07 ± 3.84 | 2266.9 | 2235.1 | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 80 | D* (Δm) | 0.0346 ± 0.0350 | 0.0047 ± 0.0280 | 7.42 ± 45.13 | 145.4 | 152.6 | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 80 | φ (m_KK) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 80 | D_s→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 80 | D⁺→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:all | 80 | Λ_c (m_pKπ) | 0.0217 ± 0.1752 | — | — | 2238.2 | — | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 95 | D* (Δm) | 0.0138 ± 0.0528 | 0.0615 ± 0.0660 | 0.22 ± 0.89 | 148.2 | 146.9 | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 95 | φ (m_KK) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 95 | D_s→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 95 | D⁺→φπ (m_KKπ) | 0.0047 ± 0.0047 | — | — | 1872.5 | — | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 95 | Λ_c (m_pKπ) | 0.1106 ± 0.3045 | 0.1350 ± 0.4263 | 0.82 ± 3.43 | 2269.0 | 2235.1 | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 80 | D* (Δm) | 0.0298 ± 0.0451 | 0.0047 ± 0.0280 | 6.39 ± 39.53 | 145.5 | 152.6 | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 80 | φ (m_KK) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 80 | D_s→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 80 | D⁺→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:J12only | 80 | Λ_c (m_pKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 95 | D* (Δm) | 0.0850 ± 0.0542 | 0.0622 ± 0.0414 | 1.37 ± 1.26 | 145.7 | 146.0 | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 95 | φ (m_KK) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 95 | D_s→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 95 | D⁺→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 95 | Λ_c (m_pKπ) | 0.2284 ± 0.2773 | — | — | 2269.2 | — | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 80 | D* (Δm) | 0.0032 ± 0.0530 | 0.0176 ± 0.0328 | 0.18 ± 3.03 | 144.0 | 147.6 | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 80 | φ (m_KK) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 80 | D_s→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 80 | D⁺→φπ (m_KKπ) | — | — | — | — | — | niemierzalne (δ = 0,5 a priori) |
| 428770:MBTSbased | 80 | Λ_c (m_pKπ) | 0.0336 ± 0.3458 | — | — | 2183.6 | — | niemierzalne (δ = 0,5 a priori) |

### 2d. Składowa z b: |parametr zderzenia D| w oknie sygnału (po odjęciu pasm bocznych)

f_B = udział szablonu MC SigB (MatchedB) w dopasowaniu χ² rozkładu |impact| (0–0,2 mm, 40 binów; ostatni bin = przepełnienie) danych SR − (2/3)·(SB dolne + górne) szablonami SigC(+Reflection) i SigB z mb_nd + fd_nd. f_B MC = prawda mb_nd. WP80. Uwaga: H używa cech wierzchołkowych, więc cięcie rzeźbi ten rozkład jednakowo w danych i MC tylko o tyle, o ile MC je modeluje.

| próbka | gatunek | pT | f_B dane | f_B MC (prawda mb_nd) | χ²/ndf | K2 (różnica w σ) |
|---|---|---|---|---|---|---|

## 4. Proponowana systematyka na składowe z MC (w % N_S; K3)

syst_k = |S_k|·δ_k, S_k = [N_S(μ_k = 1,5) − N_S(μ_k = 0,5)]/N_S(C0). δ: Reflection — tag D* (§2a, ta sama próbka; przy „niemierzalne” pT 3–30); CharmMismatch D⁺ — największe δ z pików §2c; PR (D⁰) i DecayModeMismatch (D⁺) — √((μ̂−1)² + σ²) z C2 (**ten sam fit**); PartComb — założenie 25 %. „kształt R” = |N_S(C4) − N_S(C0)|/N_S. Pogrubione: > σ_stat/3 („istotna”).

| próbka | gatunek | pT | WP | σ_stat [%] | Reflection | kształt R | CharmMismatch | DMM / PR | PartComb | **suma** | C0w−C0 (tryger) | C6−C0 (kształt Comb) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 427929:all | D0 | 3-5 | 50 | 29.2 | 0.85 | — | — | **20.82** | **12.99** | **24.55** | -13.3 | -86.6 |
| 427929:MB | D0 | 3-5 | 50 | 29.2 | 0.85 | — | — | **20.82** | **12.99** | **24.55** | -13.3 | -86.6 |
| 428770:all | D0 | 3-5 | 50 | 13.1 | 0.63 | — | — | **6.42** | **11.53** | **13.21** | -44.7 | -74.4 |
| 428770:all | D0 | 5-8 | 50 | 25.4 | 3.85 | — | — | **16.95** | 6.03 | **18.40** | -5.6 | -47.1 |
| 428770:all | D0 | 8-30 | 50 | 15.1 | 3.65 | — | — | — | — | **3.65** | +34.8 | -22.5 |
| 428770:J12only | D0 | 3-5 | 50 | 15.1 | 0.60 | — | — | **5.89** | **11.43** | **12.87** | -46.9 | -75.9 |
| 428770:J12only | D0 | 5-8 | 50 | 27.2 | 5.12 | — | — | **16.25** | 6.67 | **18.30** | -8.4 | -49.8 |
| 428770:J12only | D0 | 8-30 | 50 | 15.9 | 4.62 | — | — | — | — | **4.62** | +36.2 | -18.0 |
| 428770:MBTSbased | D0 | 8-30 | 95 | 19.1 | 0.05 | — | — | **11.25** | **15.93** | **19.50** | -21.1 | -100.0 |
| 428770:all | Dplus | 3-5 | 95 | 13.2 | — | — | 0.65 | 0.06 | **19.62** | **19.63** | -74.6 | -100.0 |
| 428770:all | Dplus | 3-5 | 80 | 18.1 | — | — | 0.15 | 0.28 | **13.99** | **13.99** | -46.7 | -84.0 |
| 428770:all | Dplus | 5-8 | 50 | 23.4 | — | — | — | — | **14.11** | **14.11** | -67.5 | -24.1 |
| 428770:all | Dplus | 8-30 | 80 | 29.0 | — | — | — | — | **14.21** | **14.21** | +21.5 | +9.7 |
| 428770:all | Dplus | 8-30 | 50 | 15.7 | — | — | — | — | — | **0.00** | -21.5 | -8.9 |
| 428770:J12only | Dplus | 3-5 | 95 | 16.3 | — | — | 0.90 | 0.03 | **19.35** | **19.37** | -76.0 | -100.0 |
| 428770:J12only | Dplus | 3-5 | 80 | 23.6 | — | — | 0.15 | 1.05 | **14.62** | **14.66** | -45.8 | -95.6 |
| 428770:J12only | Dplus | 8-30 | 80 | 31.8 | — | — | — | — | **12.37** | **12.37** | +26.3 | +9.1 |
| 428770:MBTSbased | Dplus | 3-5 | 95 | 22.4 | — | — | 0.09 | 0.46 | **20.66** | **20.67** | -73.9 | -49.0 |
| 428770:MBTSbased | Dplus | 3-5 | 80 | 32.3 | — | — | 1.03 | 2.86 | 9.97 | **10.42** | -40.9 | -48.7 |
| 428770:MBTSbased | Dplus | 5-8 | 95 | 28.0 | — | — | 0.33 | — | **16.84** | **16.85** | -14.0 | -34.6 |
| 428770:MBTSbased | Dplus | 8-30 | 95 | 24.5 | — | — | — | — | **17.60** | **17.60** | +16.3 | -65.7 |

## 5. Co z tego wynika dla głównego łańcucha analizy

- Przenosi się (zależy od fizyki i rekonstrukcji, nie od klasyfikatora): definicje kategorii, R/S refleksji i jej kształt przy danym L_xy, pozycje i względne wielkości pików D*/φ/Λ_c przed cięciem klasyfikatora, f_B, wrażliwość S_k przy podobnym S/B.
- Nie przenosi się 1:1: wartości κ po cięciu (CharmNet ma inne wejścia, m.in. ptcone i mKπ/m_D dla D⁺, ciaśniejsze WP), L_xy ≥ 0,25 mm dla D⁰ (w głównym łańcuchu 0,10), nasze weto D* dla D⁺ (w głównym łańcuchu brak weta w preselekcji), bin 2–3 GeV (nie liczony), region z tagiem AFP (nie oglądany — skład tła może być tam inny), normalizacja absolutna σ·L (tu względna do sygnału).
- Hipotezy, nie wyniki: przyczyny ewentualnych niezgodności (modelowanie zmiennych wierzchołkowych, tryger).

