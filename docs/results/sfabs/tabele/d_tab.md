### D1. Połączenie per gatunek/bin (SF drzewa H przy WP70)

| gatunek | bin | A drabina | B tag D* 2D | C przeniesienie | zgodność (max odchylenie par, χ²/ndf, p) | **SF łączny** | tylko metody bezpośrednie (A, B) | waga C |
|---|---|---|---|---|---|---|---|---|
| Dplus | 3-5 | 0.741 ± 0.165 | n.d. | 0.968 ± 0.050 | 1.31σ, 1.72/1, 0.19 | **0.949 ± 0.048**  | 0.741 ± 0.165 | 0.92 |
| Dplus | 5-8 | 1.173 ± 0.147 | n.d. | 0.968 ± 0.055 | 1.31σ, 1.71/1, 0.19 | **0.993 ± 0.052**  | 1.173 ± 0.147 | 0.88 |
| Dplus | 8-30 | 1.085 ± 0.117 | n.d. | 0.968 ± 0.057 | 0.91σ, 0.82/1, 0.36 | **0.990 ± 0.051**  | 1.085 ± 0.117 | 0.81 |
| D0 | 3-5 | — (bramki) | 0.962 ± 0.179 | 0.968 ± 0.034 | 0.03σ, 0.00/1, 0.98 | **0.967 ± 0.033**  | 0.962 ± 0.179 | 0.97 |

Test założenia C (opisowo): drabina D⁺ uśredniona po 3 binach = **1.030 ± 0.080** (χ² zgodności binów 4.23/2) wobec SF_ref = 0.968 ± 0.029 → różnica +0.74σ; bezpośredni test ogranicza efekt specyficzny dla D⁺ (niewidoczny dla C) do ±0.08 (1σ).

### D2. Przekrój 427929 z nowym SF (własna kopia rachunku; wejście: `XSEC_2026-11-19/out/xsec_results.json` — tylko odczyt, stan z 19.11 23:49, D⁰ z modelem M8)

| gatunek | bin | σ przy SF = 1 [µb] | SF | **σ [µb] ± stat ± syst ± L** | syst wzgl. (z SF) | SF wzgl. | razem | było (scen. 1 ± 10 %): razem | Pythia A3 |
|---|---|---|---|---|---|---|---|---|---|
| D0 | 3-5 | 1365 | 0.967 ± 0.033 | **1411 ± 159 ± 378 ± 32** | 26.8 % | 3.4 % | 29.2 % | 30.7 % | 1230 (dane/Pythia 1.15 ± 0.33) |
| Dplus | 3-5 | 542 | 0.949 ± 0.048 | **571 ± 45 ± 60 ± 13** | 10.4 % | 5.0 % | 13.2 % | 15.8 % | 641 (dane/Pythia 0.89 ± 0.12) |
| Dplus | 5-8 | 176 | 0.993 ± 0.052 | **178 ± 8 ± 25 ± 4** | 14.0 % | 5.2 % | 14.9 % | 17.2 % | 254 (dane/Pythia 0.70 ± 0.10) |
| Dplus | 8-30 | 55 | 0.990 ± 0.051 | **55 ± 5 ± 6 ± 1** | 11.3 % | 5.2 % | 14.4 % | 16.8 % | 85 (dane/Pythia 0.65 ± 0.09) |

**D⁺/D⁰ (3–5 GeV, 427929): R = 0.397 (SF⁺ = SF⁰) → 0.405 ± 0.056 (stat) ± 0.115 (syst)**; SF⁺/SF⁰ = 0.981 ± 0.042 (wspólny SF_ref skraca się: kowariancja 0.00088); błąd względny razem 31.6 % (było 32.9 % przy scenariuszu SF⁺/SF⁰ = 1 ± 10 %; bez SF 31.3 %). Pythia A3: 0.521.

