# SD_var: niepewność modelowa R = σ(D + p)/σ_incl(D) — POZIOM CZĄSTEK (kryterium: PLAN.md §4)

Wejście: `/home/k/Documents/charm/charm_work/SD_var_2026-11-19/runs/lxplus/merged`, `/home/k/Documents/charm/charm_work/SD_pred/runs/lxplus/merged`; scenariusze `/home/k/Documents/charm/charm_work/Q_lxplus/sd_var/scenarios.txt`. Brak danych dla: —.

Nadpisane z `/home/k/Documents/charm/charm_work/SD_var_2026-11-19/runs/lxplus/merged_close` (domknięcie): dv_CRqcd (10 części, 10000000 zd.), nv_CRqcd (22 części, 22000000 zd.).

## 1. Scenariusze

| scenariusz | zbiór | grupa | k_ND | R(D⁰) [%] | R/odniesienie | udział nie-SD [%] | R(D⁺)/R(D⁰) | R(D⁰) w binach pT 2–3 / 3–5 / 5–8 / 8–30 [%] | N(D⁰, 427929) | opis |
|---|---|---|---|---|---|---|---|---|---|---|
| A3off | N | M | 1.061 | 0.384 ± 0.005 | – | 11.8 | 1.011 | 0.433 / 0.375 / 0.324 / 0.268 | 18.6 | oficjalna konfiguracja ATLAS A3 (DL ε=0,07, mode 0, σ_SD 16 mb, bez CD) |
| A3alt | A | M | 1.061 | 0.639 ± 0.007 | 1.662 | 7.1 | 1.017 | 0.722 / 0.627 / 0.524 / 0.449 | 30.8 | ATLAS A3 ALT (ε=0,037) |
| SaS | A | M | 1.000 | 0.750 ± 0.006 | 1.952 | 10.2 | 1.001 | 0.868 / 0.719 / 0.603 / 0.540 | 35.5 | SaS (domyślny model Pythii; tune'y CMS, Monash) |
| MBR | A | M | 1.000 | 0.431 ± 0.005 | 1.121 | 13.9 | 1.006 | 0.501 / 0.417 / 0.343 / 0.275 | 20.3 | MBR domyślny (ε=0,104) |
| MBRe08 | A | M | 1.051 | 0.442 ± 0.005 | 1.150 | 13.6 | 1.025 | 0.514 / 0.430 / 0.336 / 0.309 | 20.9 | MBR ε=0,08 (CMS 7 TeV) |
| ABMST | A | M | 1.036 | 1.217 ± 0.009 | 3.165 | 4.2 | 0.987 | 1.387 / 1.180 / 0.977 / 0.927 | 58.3 | ABMST z rozszerzeniami Pythii (domyślne) |
| ABMSTnogap | B | M | 0.728 | 1.847 ± 0.019 | 4.805 | 5.9 | 1.013 | 2.101 / 1.794 / 1.476 / 1.405 | 88.5 | ABMST bez wymuszonej przerwy |
| MBRa0125 | B | M | 1.006 | 0.368 ± 0.005 | 0.956 | 16.1 | 0.996 | 0.426 / 0.357 / 0.291 / 0.233 | 17.4 | MBR α′=0,125 |
| SaSe07 | B | M | 1.000 | 0.460 ± 0.006 | 1.197 | 14.1 | 0.997 | 0.532 / 0.447 / 0.356 / 0.329 | 21.8 | SaS z ε=0,07 w widmie mas (eksploracyjne) |
| DL07own | Z | P | 1.000 | 0.356 ± 0.005 | – | 16.6 | 1.038 | 0.408 / 0.346 / 0.290 / 0.236 | 17.0 | odniesienie grup P i X: DL ε=0,07, σ jak SaS |
| DL04 | A | P | 1.000 | 0.548 ± 0.007 | 1.540 | 12.5 | 0.999 | 0.639 / 0.539 / 0.411 / 0.339 | 25.9 | ε=0,04 |
| DL058 | A | P | 1.000 | 0.419 ± 0.006 | 1.176 | 15.5 | 1.022 | 0.483 / 0.396 / 0.344 / 0.329 | 19.8 | ε=0,058 (ATLAS, przerwy w rapidity) |
| DL085 | A | P | 1.000 | 0.287 ± 0.006 | 0.805 | 19.7 | 1.001 | 0.328 / 0.280 / 0.232 / 0.189 | 13.7 | ε=0,085 (domyślne Pythii, DL) |
| ap010 | A | P | 1.000 | 0.298 ± 0.005 | 0.837 | 19.1 | 1.000 | 0.335 / 0.293 / 0.252 / 0.190 | 14.4 | α′=0,10 |
| ap040 | A | P | 1.000 | 0.389 ± 0.005 | 1.093 | 16.1 | 1.034 | 0.457 / 0.374 / 0.296 / 0.279 | 18.2 | α′=0,40 |
| BS | A | P | 1.000 | 0.323 ± 0.005 | 0.907 | 18.2 | 0.992 | 0.376 / 0.315 / 0.239 / 0.232 | 15.3 | Berger–Streng, ta sama trajektoria |
| DL02 | B | P | 1.000 | 0.704 ± 0.008 | 1.978 | 10.3 | 1.028 | 0.804 / 0.683 / 0.581 / 0.474 | 33.7 | ε=0,02 (granica Pythii) |
| DL10 | B | P | 1.000 | 0.216 ± 0.004 | 0.607 | 25.5 | 1.031 | 0.249 / 0.214 / 0.162 / 0.143 | 10.3 | ε=0,10 (poza zakresem pomiaru ATLAS) |
| BI | B | P | 1.000 | 0.584 ± 0.006 | 1.639 | 11.1 | 1.008 | 0.659 / 0.565 / 0.485 / 0.448 | 28.1 | Bruni–Ingelman |
| gap | B | P | 1.000 | 0.362 ± 0.006 | 1.016 | 16.5 | 0.987 | 0.417 / 0.355 / 0.281 / 0.240 | 17.2 | tłumienie małych przerw |
| H1A | B | P | 1.000 | 0.111 ± 0.003 | 0.312 | 43.3 | 1.025 | 0.133 / 0.106 / 0.080 / 0.071 | 5.1 | strumień H1 Fit A (dopasowany do twardej dyfrakcji ep) |
| H1B | B | P | 1.000 | 0.122 ± 0.003 | 0.344 | 40.5 | 1.004 | 0.135 / 0.122 / 0.105 / 0.089 | 6.0 | strumień H1 Fit B |
| SaSown | Z | P | 1.000 | 0.816 ± 0.008 | 2.290 | 9.3 | 1.028 | 0.935 / 0.792 / 0.661 / 0.548 | 38.9 | zamknięcie mode 0 wobec mode 1 (SaS) |
| sigPp5 | A | X | 1.000 | 0.651 ± 0.007 | 1.828 | 10.5 | 0.984 | 0.705 / 0.661 / 0.554 / 0.461 | 32.5 | σ_ref(Pom–p)=5 mb |
| sigPp25 | A | X | 1.000 | 0.167 ± 0.004 | 0.469 | 34.0 | 0.990 | 0.207 / 0.149 / 0.129 / 0.114 | 7.4 | σ_ref(Pom–p)=25 mb |
| pomH1A | A | X | 1.000 | 0.314 ± 0.005 | 0.881 | 19.2 | 0.992 | 0.352 / 0.292 / 0.284 / 0.271 | 15.1 | PDF Pomeronu H1 2006 Fit A |
| pomH1J | A | X | 1.000 | 0.418 ± 0.006 | 1.173 | 16.9 | 1.049 | 0.486 / 0.388 / 0.364 / 0.281 | 19.6 | PDF Pomeronu H1 2007 Jets |
| pomGKG | A | X | 1.000 | 0.518 ± 0.006 | 1.456 | 14.7 | 1.039 | 0.598 / 0.510 / 0.408 / 0.305 | 24.6 | PDF Pomeronu GKG18 LO B |
| pompi | B | X | 1.000 | 0.888 ± 0.010 | 2.492 | 13.0 | 1.006 | 0.945 / 0.914 / 0.769 / 0.604 | 44.8 | PDF Pomeronu jak π⁰ |
| pT0lo | A | X | 1.000 | 0.336 ± 0.005 | 0.944 | 15.6 | 1.053 | 0.398 / 0.312 / 0.271 / 0.235 | 15.5 | pT0Ref=2,30 (dyfrakcja i ND spójnie) |
| pT0hi | A | X | 1.000 | 0.376 ± 0.006 | 1.057 | 20.2 | 0.992 | 0.423 / 0.371 / 0.309 / 0.273 | 18.3 | pT0Ref=2,60 (dyfrakcja i ND spójnie) |
| CRqcd | A | X | 1.000 | 0.512 ± 0.007 | 1.437 | 39.5 | 0.988 | 0.580 / 0.497 / 0.444 / 0.331 | 24.6 | rekoneksja QCD + nowy model resztek (dyfrakcja i ND spójnie) |
| NDMonash | A | L | 1.000 | 0.326 ± 0.005 | 0.916 | 15.6 | 1.044 | 0.367 / 0.317 / 0.274 / 0.236 | 15.7 | proton wiodący i mianownik z ND Monash |
| NDCRqcd | A | L | 1.000 | 0.563 ± 0.008 | 1.581 | 35.8 | 1.025 | 0.635 / 0.558 / 0.468 / 0.358 | 27.3 | proton wiodący i mianownik z ND CRqcd |
| hard_t4 | N | H | – | 1.061 ± 0.007 | – | 14.1 | 1.041 | 1.414 / 0.803 / 0.403 / 0.289 | 35.3 | twarda dyfrakcja, strumień DL-A3 (nominał scenariusza twardego) |
| hard_SaS | A | H | – | 1.491 ± 0.008 | 1.405 | 10.0 | 1.026 | 1.978 / 1.144 / 0.556 / 0.398 | 50.0 | strumień SaS |
| hard_pT0lo | A | H | – | 0.945 ± 0.007 | 0.891 | 14.8 | 1.039 | 1.275 / 0.706 / 0.334 / 0.243 | 30.8 | pT0Ref=2,30 |
| hard_pT0hi | A | H | – | 1.183 ± 0.010 | 1.115 | 13.7 | 1.039 | 1.545 / 0.913 / 0.514 / 0.351 | 40.8 | pT0Ref=2,60 |
| hard_bSel2 | A | H | – | 1.084 ± 0.007 | 1.021 | 13.8 | 1.038 | 1.430 / 0.841 / 0.405 / 0.302 | 36.8 | bSelHard=2 |
| hard_pomGKG | A | H | – | 1.147 ± 0.008 | 1.081 | 13.0 | 1.044 | 1.545 / 0.856 / 0.409 / 0.272 | 37.3 | PDF Pomeronu GKG18 LO B |
| hard_pomH1A | A | H | – | 1.542 ± 0.010 | 1.453 | 9.7 | 1.034 | 1.983 / 1.214 / 0.728 / 0.588 | 55.4 | strumień i PDF H1 2006 Fit A |
| hard_H1B | A | H | – | 2.274 ± 0.012 | 2.143 | 6.6 | 1.018 | 3.040 / 1.703 / 0.869 / 0.640 | 75.1 | strumień H1 Fit B (PDF domyślne = Fit B LO) |
| hard_pT0x | B | H | – | 1.568 ± 0.009 | 1.477 | 12.4 | 1.028 | 1.985 / 1.251 / 0.752 / 0.632 | 57.2 | pT0Ref=3,00 (kraniec skanu 1901.05261) |
| hard_bSel3 | B | H | – | 1.198 ± 0.007 | 1.129 | 12.4 | 1.038 | 1.533 / 0.983 / 0.504 / 0.302 | 43.1 | bSelHard=3 (skrajny wg dokumentacji) |

## 2. Niepewność (obwiednia zbioru A w grupach; czynniki mnożące R nominalne)

| gatunek | R nom [%] | M (pełne modele) | P (parametry kształtu) | X (układ X) | L (proton ND) | normalizacja σ | **razem** | R zakres [%] | obwiednia A+B (M / P / X) |
|---|---|---|---|---|---|---|---|---|---|
| D⁰ | 0.384 ± 0.005 | ×1.00 … ×3.17 | ×0.81 … ×1.54 | ×0.47 … ×1.83 | ×0.92 … ×1.58 | ×0.63 … ×1.07 | **×0.40 … ×3.98** | 0.154 – 1.529 | ×0.96 … ×4.80 / ×0.31 … ×1.98 / ×0.47 … ×2.49 |
| D⁺ | 0.389 ± 0.007 | ×1.00 … ×3.09 | ×0.78 … ×1.48 | ×0.45 … ×1.73 | ×0.92 … ×1.56 | ×0.63 … ×1.07 | **×0.38 … ×3.79** | 0.149 – 1.474 | ×0.94 … ×4.81 / ×0.31 … ×1.96 / ×0.45 … ×2.42 |
| D*⁺ | 0.390 ± 0.009 | ×1.00 … ×3.20 | ×0.81 … ×1.54 | ×0.45 … ×1.79 | ×0.92 … ×1.59 | ×0.63 … ×1.07 | **×0.39 … ×3.98** | 0.151 – 1.555 | ×0.95 … ×4.76 / ×0.31 … ×1.93 / ×0.45 … ×2.43 |

Scenariusz twardy (osobno, nie wchodzi do sumy): D⁰ R = 1.061 % ×0.89 … ×2.14; D⁺ R = 1.105 % ×0.89 … ×2.10; D*⁺ R = 1.103 % ×0.89 … ×2.14

Uzysk D⁰ w 427929 (pT 3–30, scenariusz M z TAGREG, tylko niepewność modelu R): 18.6 (7.4 – 73.8).
