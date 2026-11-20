# BKGCAT — rozpady D w prawdzie generatora mb_nd vs PDG 2026 (POST HOC)

mb_nd: 78 plików, 19220986 zdarzeń, 8320523 wpisów TruthDmesons; hadrony inne niż D⁰, D⁺, D*⁺: 0 (**D_s i Λ_c nie są zapisywane w TruthDmesons — tablicy dla nich nie da się zbudować**).
Region: |η| < 2,5 (klasa brakującej masy wymaga precyzji float; ułamki bez cięcia η i dla pT ≥ 3 GeV, |η| < 2,2 — kolumny kontrolne). Prompt + z b razem.
Lista córek zawiera tylko π±, K±, π⁰ (rezonanse rozwinięte; K⁰, η→γγ, leptony, ν, γ niewidoczne) → tryb = topologia widzialna + klasa braku.

## |pdg| = 421: 1992686 wpisów w |η| < 2,5 (udział z b: 9.5 %)

### Porównanie z PDG

| tryb (topologia, klasy) | MC [%] | ± stat | PDG [%] | ± | MC/PDG | źródło PDG | MC bez cięcia η | MC pT ≥ 3, \|η\| < 2,2 |
|---|---|---|---|---|---|---|---|---|
| K- pi+ | 3.974 | 0.014 | 3.936 | 0.03 | 1.010 | Γ47 K⁻π⁺ (PDG zawiera FSR) | 3.985 | 3.955 |
| K+ pi- | 0.014 | 0.001 | 0.0148 | 0.0007 | 0.946 | K⁺π⁻ DCS + mieszanie (BR_frag §1, PDG 2026) | 0.014 | 0.015 |
| K- pi+ pi0 | 14.660 | 0.027 | 14.4 | 0.5 | 1.018 | Γ68 K⁻π⁺π⁰ | 14.633 | 14.589 |
| K- pi+ pi+ pi- | 7.862 | 0.020 | 8.2 | 0.14 | 0.959 | Γ87 K⁻2π⁺π⁻ | 7.873 | 7.859 |
| K- pi+ pi0 pi0 | 11.130 | 0.024 | 8.86 | 0.23 | 1.256 | Γ101 K⁻π⁺2π⁰ | 11.106 | 11.146 |
| K- pi+ pi+ pi- pi0 | 4.042 | 0.014 | 4.3 | 0.4 | 0.940 | Γ104 K⁻2π⁺π⁻π⁰ | 4.046 | 4.052 |
| K- K+ | 0.415 | 0.005 | 0.407 | 0.006 | 1.020 | Γ266 K⁺K⁻ | 0.413 | 0.410 |
| pi+ pi- | 0.148 | 0.003 | 0.1451 | 0.0024 | 1.019 | Γ151 π⁺π⁻ | 0.147 | 0.149 |
| pi+ pi- pi0 | 1.545 | 0.009 | 1.49 | 0.06 | 1.037 | Γ153 π⁺π⁻π⁰ | 1.552 | 1.560 |
| pi+ pi- pi0 pi0 | 0.976 | 0.007 | 1.002 | 0.031 | 0.974 | Γ212 π⁺π⁻2π⁰ | 0.976 | 0.976 |
| K- K+ pi0 | 0.481 | 0.005 | 0.341 | 0.014 | 1.412 | Γ289 K⁺K⁻π⁰ | 0.479 | 0.476 |
| K- K+ pi+ pi- | 0.354 | 0.004 | 0.278 | 0.017 | 1.275 | Γ303 K⁺K⁻π⁺π⁻ | 0.350 | 0.358 |
| pi0 pi0 | 0.084 | 0.002 | 0.0825 | 0.0025 | 1.019 | Γ152 2π⁰ | 0.084 | 0.086 |
| K- [wszystkie klasy] | 6.949 | 0.019 | 6.955 | 0.025 | 0.999 | Γ20 + Γ21 K⁻e⁺ν + K⁻μ⁺ν (topologia „samo K⁻”, wszystkie klasy) | 6.945 | 6.951 |
| pi- [wszystkie klasy] | 3.684 | 0.014 | 3.405 | 0.07 | 1.082 | Γ39 + Γ40 + Γ25 + Γ29: π⁻ℓ⁺ν + K̄⁰π⁻ℓ⁺ν (topologia „samo π⁻”) | 3.700 | 3.706 |
| pi+ pi-+K0 | 5.397 | 0.016 | 5.68 | 0.32 | 0.950 | 2 × Γ53 K⁰_Sπ⁺π⁻ (założenie K⁰_L = K⁰_S) | 5.221 | 5.382 |
| pi0+K0 | 2.247 | 0.011 | 2.215 | 0.04 | 1.014 | Γ48 + Γ49 K⁰_Sπ⁰ + K⁰_Lπ⁰ | 2.209 | 2.248 |
| K- pi++K0 | 0.638 | 0.006 | 0.68 | 0.1 | 0.938 | 2 × Γ268 K⁰_SK⁻π⁺ (założenie K⁰_L = K⁰_S) | 0.621 | 0.645 |
| K- pi++eta | 0.402 | 0.004 | 0.74 | 0.02 | 0.543 | Γ113 K⁻π⁺η × B(η→2γ) = 1,88 % × 39,36 % | 0.389 | 0.398 |

### Wszystkie tryby ≥ 0,2 % (topologia × klasa braku)

| tryb | MC [%] | N |
|---|---|---|
| pi+ pi- pi0+K0 | 16.401 | 326818 |
| K- pi+ pi0 | 14.230 | 283552 |
| K- pi+ pi0 pi0 | 11.061 | 220410 |
| K- pi+ pi+ pi- | 7.457 | 148593 |
| K-+X | 6.829 | 136085 |
| pi+ pi-+K0 | 5.397 | 107545 |
| K- pi+ pi+ pi- pi0 | 3.985 | 79401 |
| K- pi+ | 3.841 | 76547 |
| pi-+X | 3.679 | 73304 |
| pi0 pi0+K0 | 2.312 | 46073 |
| pi0+K0 | 2.247 | 44770 |
| K- pi0+X | 1.545 | 30796 |
| pi+ pi-+X | 1.536 | 30601 |
| pi+ pi- pi0 | 1.443 | 28750 |
| pi+ pi- pi0 pi0+K0 | 1.279 | 25487 |
| pi+ pi- pi0 pi0 | 0.957 | 19065 |
| K- K++K0 | 0.887 | 17672 |
| (nic)+X | 0.876 | 17462 |
| pi0+X | 0.755 | 15050 |
| pi+ pi+ pi- pi- | 0.741 | 14756 |
| pi0 pi0 pi0+K0 | 0.727 | 14494 |
| K- pi+ pi0 pi0 pi0 | 0.681 | 13567 |
| K- pi++K0 | 0.638 | 12711 |
| pi+ pi- pi0 pi0 pi0 | 0.589 | 11740 |
| pi+ pi- pi0+X | 0.588 | 11726 |
| pi+ pi+ pi- pi-+K0 | 0.540 | 10757 |
| K- K+ pi0 | 0.480 | 9567 |
| K- pi+ pi0+X0 | 0.430 | 8577 |
| pi- pi0+X | 0.422 | 8406 |
| K+ pi-+K0 | 0.415 | 8273 |
| K- K+ | 0.408 | 8126 |
| K- pi+ pi+ pi-+X0 | 0.405 | 8075 |
| K- pi++eta | 0.402 | 8005 |
| pi0 pi0 pi0 pi0+K0 | 0.387 | 7721 |
| K- pi++X | 0.384 | 7651 |
| pi+ pi- pi0 pi0 pi0+K0 | 0.357 | 7119 |
| K- K+ pi+ pi- | 0.351 | 6998 |
| pi0 pi0+X | 0.351 | 6989 |
| pi+ pi+ pi- pi- pi0 | 0.334 | 6649 |
| K- K+ pi+ pi- pi0 | 0.308 | 6139 |
| K- K+ pi0+K0 | 0.284 | 5663 |
| K- K+ pi0 pi0 | 0.276 | 5493 |

## |pdg| = 411: 1028169 wpisów w |η| < 2,5 (udział z b: 8.4 %)

### Porównanie z PDG

| tryb (topologia, klasy) | MC [%] | ± stat | PDG [%] | ± | MC/PDG | źródło PDG | MC bez cięcia η | MC pT ≥ 3, \|η\| < 2,2 |
|---|---|---|---|---|---|---|---|---|
| K- pi+ pi+ | 9.404 | 0.030 | 9.31 | 0.14 | 1.010 | Γ68 K⁻2π⁺ | 9.401 | 9.336 |
| K- pi+ pi+ pi0 | 7.064 | 0.026 | 6.35 | 0.08 | 1.112 | Γ93 K⁻2π⁺π⁰ | 7.089 | 7.091 |
| K- K+ pi+ | 1.013 | 0.010 | 0.961 | 0.017 | 1.054 | Γ177 K⁺K⁻π⁺ | 1.004 | 1.013 |
| K- K+ pi+ pi0 | 3.000 | 0.017 | 0.662 | 0.032 | 4.531 | Γ186 K⁺K⁻π⁺π⁰ | 2.985 | 2.999 |
| pi+ pi+ pi- | 0.480 | 0.007 | 0.327 | 0.009 | 1.467 | Γ131 2π⁺π⁻ | 0.479 | 0.471 |
| pi+ pi0 | 0.128 | 0.004 | 0.1243 | 0.0033 | 1.032 | Γ130 π⁺π⁰ | 0.128 | 0.138 |
| pi+ pi+ pi- pi0 | 1.444 | 0.012 | 1.165 | 0.03 | 1.240 | Γ147 2π⁺π⁻π⁰ | 1.433 | 1.453 |
| K- pi+ pi+ pi+ pi- | 0.515 | 0.007 | 0.57 | 0.05 | 0.903 | Γ117 K⁻3π⁺π⁻ | 0.523 | 0.516 |
| K- pi+ [wszystkie klasy] | 7.979 | 0.028 | 7.62 | 0.37 | 1.047 | Γ22 + Γ33 K⁻π⁺e⁺ν + K⁻π⁺μ⁺ν (topologia „K⁻π⁺”) | 7.991 | 7.983 |
| pi++K0 | 2.893 | 0.017 | 3.014 | 0.06 | 0.960 | Γ66 + Γ67 K⁰_Sπ⁺ + K⁰_Lπ⁺ | 2.847 | 2.888 |
| pi+ pi0+K0 | 14.505 | 0.038 | 14.62 | 0.4 | 0.992 | 2 × Γ78 K⁰_Sπ⁺π⁰ (założenie K⁰_L = K⁰_S) | 14.000 | 14.517 |
| pi+ pi+ pi-+K0 | 4.310 | 0.020 | 6.16 | 0.16 | 0.700 | 2 × Γ107 K⁰_S2π⁺π⁻ (założenie K⁰_L = K⁰_S) | 4.150 | 4.284 |
| K++K0 | 0.598 | 0.008 | 0.623 | 0.018 | 0.960 | Γ170 + Γ171 K⁰_SK⁺ + K⁰_LK⁺ | 0.586 | 0.608 |
| (nic) [wszystkie klasy] | 17.680 | 0.041 | 17.49 | 0.12 | 1.011 | Γ17 + Γ18 K̄⁰e⁺ν + K̄⁰μ⁺ν (dolne ograniczenie topologii „nic widzialnego”) | 17.696 | 17.696 |

### Wszystkie tryby ≥ 0,2 % (topologia × klasa braku)

| tryb | MC [%] | N |
|---|---|---|
| (nic)+X | 17.680 | 181784 |
| pi+ pi0+K0 | 14.505 | 149141 |
| K- pi+ pi+ | 9.220 | 94798 |
| K- pi++X | 7.804 | 80241 |
| K- pi+ pi+ pi0 | 6.967 | 71637 |
| pi+ pi0 pi0+K0 | 5.126 | 52707 |
| pi0+X | 4.707 | 48398 |
| pi+ pi+ pi-+K0 | 4.310 | 44319 |
| K- K+ pi+ pi0 | 2.987 | 30712 |
| pi++K0 | 2.893 | 29740 |
| pi++X | 1.706 | 17545 |
| pi+ pi+ pi- pi0 | 1.347 | 13846 |
| pi+ pi0+X | 1.232 | 12669 |
| pi+ pi+ pi- pi0+K0 | 1.178 | 12111 |
| K- K+ pi+ | 1.002 | 10306 |
| pi+ pi+ pi- pi0 pi0 pi0 | 0.964 | 9912 |
| pi+ pi+ pi- pi0 pi0 pi0 pi0 | 0.898 | 9230 |
| pi+ pi+ pi-+eta | 0.865 | 8889 |
| K+ pi0+K0 | 0.810 | 8325 |
| pi+ pi+ pi- pi0+eta | 0.802 | 8248 |
| K- pi+ pi+ pi0 pi0 | 0.725 | 7456 |
| K++K0 | 0.598 | 6151 |
| K- pi+ pi++K0 | 0.563 | 5787 |
| K++X | 0.525 | 5400 |
| K- pi+ pi+ pi+ pi- | 0.506 | 5203 |
| pi+ pi+ pi+ pi- pi- pi0 | 0.500 | 5136 |
| pi+ pi0 pi0 pi0 | 0.500 | 5136 |
| pi+ pi-+X | 0.492 | 5058 |
| pi+ pi+ pi+ pi- pi- pi0 pi0 | 0.472 | 4851 |
| pi+ pi0 pi0 | 0.458 | 4711 |
| pi+ pi0 pi0+eta | 0.439 | 4510 |
| K+ pi+ pi-+K0 | 0.435 | 4470 |
| pi+ pi0 pi0 pi0+K0 | 0.434 | 4461 |
| pi+ pi+ pi-+X | 0.424 | 4363 |
| pi+ pi0 pi0 pi0+eta | 0.416 | 4277 |
| pi+ pi+ pi- pi0 pi0 | 0.378 | 3887 |
| pi+ pi0 pi0 pi0 pi0 pi0 | 0.358 | 3681 |
| pi+ pi- pi0+X | 0.354 | 3637 |
| pi+ pi0 pi0 pi0 pi0 pi0 pi0 | 0.327 | 3360 |
| pi+ pi+ pi- | 0.308 | 3165 |
| K- pi+ pi+ pi+ pi- pi0 | 0.228 | 2349 |

## |pdg| = 413: 714719 wpisów w |η| < 2,5 (udział z b: 10.9 %)

D*⁺ → D⁰π⁺ (lista zawiera D⁰): 67.72 ± 0.10 %; PDG Γ1 67,7 ± 0,5 % (MC/PDG = 1.000); bez D⁰ (D⁺π⁰, D⁺γ): 32.28 %; PDG Γ2 + Γ3 = 32,3 %.

### Wszystkie tryby ≥ 0,2 % (topologia × klasa braku)

| tryb | MC [%] | N |
|---|---|---|
| D0 pi+ pi+ pi- pi0+K0 | 11.080 | 79193 |
| D0 K- pi+ pi+ pi0 | 9.666 | 69087 |
| D0 K- pi+ pi+ pi0 pi0 | 7.460 | 53320 |
| pi0+X | 5.479 | 39160 |
| D0 K- pi+ pi+ pi+ pi- | 5.019 | 35871 |
| D0 K- pi++X | 4.635 | 33128 |
| pi+ pi0 pi0+K0 | 4.477 | 31998 |
| D0 pi+ pi+ pi-+K0 | 3.659 | 26149 |
| K- pi+ pi+ pi0 | 2.836 | 20271 |
| D0 K- pi+ pi+ pi+ pi- pi0 | 2.696 | 19269 |
| D0 K- pi+ pi+ | 2.610 | 18652 |
| D0 pi+ pi-+X | 2.494 | 17826 |
| K- pi+ pi0+X | 2.401 | 17159 |
| K- pi+ pi+ pi0 pi0 | 2.130 | 15226 |
| pi+ pi0 pi0 pi0+K0 | 1.584 | 11324 |
| D0 pi+ pi0 pi0+K0 | 1.573 | 11240 |
| D0 pi+ pi0+K0 | 1.498 | 10707 |
| pi0 pi0+X | 1.443 | 10315 |
| pi+ pi+ pi- pi0+K0 | 1.330 | 9504 |
| D0 K- pi+ pi0+X | 1.063 | 7595 |
| D0 pi+ pi+ pi-+X | 1.056 | 7546 |
| D0 pi+ pi+ pi- pi0 | 1.002 | 7159 |
| K- K+ pi+ pi0 pi0 | 0.905 | 6465 |
| pi+ pi0+K0 | 0.879 | 6285 |
| D0 pi+ pi+ pi- pi0 pi0+K0 | 0.873 | 6240 |
| pi+ pi0+X | 0.776 | 5545 |
| D0 pi+ pi+ pi- pi0 pi0 | 0.650 | 4648 |
| D0 K- K+ pi++K0 | 0.608 | 4342 |
| D0 pi++X | 0.592 | 4230 |
| D0 pi+ pi0+X | 0.514 | 3671 |
| D0 pi+ pi+ pi+ pi- pi- | 0.497 | 3551 |
| D0 pi+ pi0 pi0 pi0+K0 | 0.491 | 3509 |
| pi+ pi0 pi0+X | 0.466 | 3328 |
| D0 K- pi+ pi+ pi0 pi0 pi0 | 0.464 | 3319 |
| D0 K- pi+ pi++K0 | 0.426 | 3042 |
| pi+ pi+ pi- pi0 pi0 | 0.415 | 2965 |
| D0 pi+ pi+ pi- pi0 pi0 pi0 | 0.406 | 2901 |
| D0 pi+ pi+ pi- pi0+X | 0.392 | 2799 |
| pi+ pi+ pi- pi0 pi0+K0 | 0.369 | 2635 |
| D0 pi+ pi+ pi+ pi- pi-+K0 | 0.366 | 2615 |
| D0 K- K+ pi+ pi0 | 0.325 | 2323 |
| K- K+ pi+ pi0 | 0.313 | 2234 |
| pi+ pi+ pi- pi0 pi0 pi0 pi0 | 0.302 | 2157 |
| D0 K- pi+ pi+ pi0+X0 | 0.284 | 2029 |
| (nic)+X | 0.283 | 2025 |
| D0 K+ pi+ pi-+K0 | 0.280 | 2004 |
| D0 K- K+ pi+ | 0.276 | 1972 |
| D0 pi+ pi- pi0+X | 0.276 | 1971 |
| D0 K- pi+ pi++eta | 0.275 | 1962 |
| D0 K- pi+ pi+ pi+ pi-+X0 | 0.268 | 1918 |
| pi+ pi+ pi- pi0+eta | 0.266 | 1900 |
| D0 K- pi+ pi++X | 0.265 | 1894 |
| pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0 | 0.265 | 1892 |
| D0 pi+ pi0 pi0 pi0 pi0+K0 | 0.260 | 1857 |
| pi+ pi+ pi- pi0 pi0+eta | 0.251 | 1794 |
| K+ pi0 pi0+K0 | 0.249 | 1782 |
| D0 pi+ pi+ pi- pi0 pi0 pi0+K0 | 0.245 | 1750 |
| D0 K- K+ pi+ pi+ pi- | 0.229 | 1635 |
| D0 pi+ pi0 pi0+X | 0.227 | 1623 |
| D0 pi+ pi+ pi+ pi- pi- pi0 | 0.227 | 1622 |
| K- pi+ pi+ pi0 pi0 pi0 | 0.223 | 1591 |

## D⁰, decayType 10 („K π”): skład klas braku (|η| < 2,5)

Wpisów: 107549 = 5.397 % wszystkich; topologie w tej klasie: ['K- pi+'].

| klasa | N | udział | × ułamek decayType = BR_MC [%] |
|---|---|---|---|
| ekskluzywny (m_D − m_wid < 30 MeV) | 76547 | 71.17 % | 3.841 |
| +X0 (radiacyjny, m_brak < 100 MeV) | 2635 | 2.45 % | 0.132 |
| +K0 | 12711 | 11.82 % | 0.638 |
| +eta | 8005 | 7.44 % | 0.402 |
| +X (inne) | 7651 | 7.11 % | 0.384 |

**Ułamek prawdziwego trybu sygnału (ekskl. + radiacyjny) w decayType: 73.62 ± 0.13 %; BR_MC(tryb sygnału) = 3.974 % wobec PDG 3.936 % (MC/PDG = 1.010); „cały decayType”/PDG = 1.371.**

## D⁺, decayType 4 („K π π”): skład klas braku (|η| < 2,5)

Wpisów: 102854 = 10.004 % wszystkich; topologie w tej klasie: ['K- pi+ pi+'].

| klasa | N | udział | × ułamek decayType = BR_MC [%] |
|---|---|---|---|
| ekskluzywny (m_D − m_wid < 30 MeV) | 94798 | 92.17 % | 9.220 |
| +X0 (radiacyjny, m_brak < 100 MeV) | 1890 | 1.84 % | 0.184 |
| +K0 | 5787 | 5.63 % | 0.563 |
| +eta | 3 | 0.00 % | 0.000 |
| +X (inne) | 376 | 0.37 % | 0.037 |

**Ułamek prawdziwego trybu sygnału (ekskl. + radiacyjny) w decayType: 94.01 ± 0.07 %; BR_MC(tryb sygnału) = 9.404 % wobec PDG 9.31 % (MC/PDG = 1.010); „cały decayType”/PDG = 1.075.**

## Skutek uboczny dla XSEC (C = N_reco^prawda / N_gen(decayType sygnału)) — do weryfikacji w studium XSEC; XSEC nie zmieniany

f_gen = ułamek prawdziwego trybu sygnału (ekskl. + radiacyjny) w mianowniku (pT ≥ 3 GeV, |η| < 2,2); f_reco = to samo w liczniku (Sig w zakresie mas, WP70).
Jeżeli σ ma dotyczyć trybu z BR PDG, to σ_popr/σ_XSEC = f_gen/f_reco.

| gat. | bin | f_gen | f_reco (WP70) | N_reco Sig: ekskl. + rad. / z brakiem (K0, η, X) | σ_popr/σ_XSEC |
|---|---|---|---|---|---|
| D0 | 3-5 | 0.7342 | 1.0000 | 3653 / 0 | 0.734 |
| D0 | 5-8 | 0.7342 | 1.0000 | 1919 / 0 | 0.734 |
| D0 | 8-30 | 0.7342 | 1.0000 | 848 / 0 | 0.734 |
| Dplus | 3-5 | 0.9387 | 1.0000 | 3245 / 0 | 0.939 |
| Dplus | 5-8 | 0.9387 | 1.0000 | 1827 / 0 | 0.939 |
| Dplus | 8-30 | 0.9387 | 1.0000 | 843 / 0 | 0.939 |

m_brak dla topologii „K- pi+” (421), klasy ≠ ekskluzywna: największe biny 5 MeV: 0–5 MeV: 2461, 495–500 MeV: 12690, 545–550 MeV: 7992, 955–960 MeV: 356, 995–1000 MeV: 384, 1000–1005 MeV: 410, 1010–1015 MeV: 367, 1015–1020 MeV: 361; razem 31002.

m_brak dla topologii „K- pi+ pi+” (411), klasy ≠ ekskluzywna: największe biny 5 MeV: 0–5 MeV: 1760, 5–10 MeV: 39, 495–500 MeV: 5779, 500–505 MeV: 38, 995–1000 MeV: 86, 1000–1005 MeV: 63, 1005–1010 MeV: 33, 1010–1015 MeV: 27; razem 8056.

