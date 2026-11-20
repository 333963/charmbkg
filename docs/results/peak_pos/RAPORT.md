# PEAK_pos — pozycja i szerokość piku masy D vs D̄ z danych (16.11.2026)

**Cel.** (a) rozmiar systematyk PEAK_POSITION (placeholder w głównym łańcuchu analizy ±0,2 MeV, też dla Δm D*) i PEAK_WIDTH (±10 %);
(b) oszacowanie z danych „weak-mode bias” torów (sagitta, brak map Run 3; punkt M17 głównego łańcucha) z m(D⁺)−m(D⁻) i m(D⁰)−m(D̄⁰).

## Krótko
| wielkość (WP70, średnia ważona binów pT 3–5/5–8/8–30) | Δμ = dane − MC | σ_dane/σ_MC | (D−D̄)_dane − (D−D̄)_MC |
|---|---|---|---|
| D⁺ m(Kππ) (30 % lumi) | −0,79 ± 0,74 MeV | 0,966 ± 0,048 | +0,90 ± 1,49 MeV |
| D⁰ m(Kπ) z tagu D* (pełna lumi, bez refleksji) | −0,10 ± 0,71 MeV | 0,947 ± 0,047 | −1,39 ± 1,38 MeV |
| **D⁺ ⊕ D⁰(tag)** | **−0,43 ± 0,51 MeV** | **0,956 ± 0,034** | — |
| Δm(D*) (pełna lumi) | −55 ± 24 keV (χ²=11,7/2; 3–5 GeV: −211 ± 63 keV) | 1,038 ± 0,030 | −47 ± 48 keV |
| D⁰ inkluzywnie (weto refl., 30 % lumi) — tylko kontrola | −1,7 ± 1,3 MeV (2 biny; 3–5 niestabilny) | 0,83 ± 0,09 | −0,5 ± 2,8 MeV |

WP50 (stabilność): D⁺ Δμ +0,09 ± 0,76 MeV, σ 1,064 ± 0,051; D⁰(tag) −0,23 ± 1,04 MeV, σ 1,03 ± 0,08; Δm −18 ± 26 keV, σ 1,050 ± 0,032.
Rozrzut od modelu tła/okna/binowania (cheb2→cheb3/exp·pol1, węższe okno, bin 1 MeV; Δm: pol2 tła, zakres 165 MeV,
okno m(Kπ) ±60 MeV; m(Kπ)-tag: okno Δm ±0,8 MeV): Δμ ≤ 0,3–0,8 MeV (m(D)), ≤ 34 keV (Δm 3–5), 11–18 keV (Δm > 5 GeV); σ: 3–15 %/bin.
Bezwzględnie: D⁺ dane 1869,7 MeV (PDG 1869,66), D⁰ z tagu 1865,6 MeV (MC też 1865,6; PDG 1864,84), Δm dane ≈ 145,44 MeV (PDG 145,43).

## Rekomendacja
- **PEAK_POSITION m(Kπ)/m(Kππ): ±0,5 MeV** (= 1σ dokładność danych: −0,43 ± 0,51 MeV, zgodne z 0; 95 % CL |Δμ| ≲ 1,4 MeV).
  Placeholder ±0,2 MeV nie jest wykluczony, ale też nie jest wykazany danymi — dane go nie uzasadniają przy tej statystyce.
  Wpływ na liczbę sygnału przy σ ≈ 26–34 MeV jest pomijalny (przy ustalonym μ: ΔN/N ≈ (Δμ/σ)²/2 ≈ 0,01–0,02 %; przy wolnym μ — zero).
- **PEAK_POSITION Δm(D*): zależna od pT** — 3–5 GeV **±0,2 MeV** (dane faktycznie przesunięte o −0,21 ± 0,06(stat) ± 0,03(model) MeV przy WP70,
  −0,08 ± 0,05 przy WP50), 5–8 GeV ±0,1 MeV (−0,11 ± 0,05), 8–30 GeV ±0,05 MeV (+0,01 ± 0,03). Lepiej: w fitach D* puścić μ wolno
  albo przesunąć szablon MC o zmierzone Δμ i jako syst. wziąć jego błąd. Przy σ_Δm ≈ 0,9–1,0 MeV przesunięcie 0,2 MeV to ~0,2σ → przy ustalonym μ ΔN/N ≈ 2 %.
- **PEAK_WIDTH: ±10 % jest bezpieczne (zachowawcze)**; dane: σ_dane/σ_MC = 0,96 ± 0,03 (m(D), WP70), 1,06 ± 0,05 (D⁺ WP50), Δm 1,04 ± 0,03.
  Wersja „z danych”: m(D) ±6 % (|1−σ| ⊕ różnica WP70/WP50), Δm ±5 %. Brak różnicy szerokości D vs D̄ (poza fluktuacjami pojedynczych binów, ≤ 2σ).
- **Weak-mode (M17): brak efektu na poziomie czułości; |δ_sagitta| ≲ 1,4 TeV⁻¹ (95 %, uśrednione po η)** — patrz §3. Dla analizy pomijalne:
  m(D⁺) i m(D⁻) przesuwają się przeciwnie (±1,3 MeV na 1 TeV⁻¹), więc w piku D⁺+D⁻ to tylko poszerzenie (≤ 0,3 % σ przy δ = 1,5 TeV⁻¹);
  D⁰ jest praktycznie nieczuły (0,06–0,1 MeV/TeV⁻¹, córki o przeciwnych ładunkach). Przy fitach osobno D⁺ i D⁻ (LS/OS) z ustalonym μ:
  ΔN/N ≲ (2/26)²/2 ≈ 0,3 % — pokryte przez PEAK_POSITION ±0,5 MeV ⊕ PEAK_WIDTH.

## 1. Dane, selekcja, metoda
- **D⁰**: `HYP_d0/runs/lxplus/data_428770.parquet` (30 % lumi, model analizy D0_H, progi `models/D0_H/meta.json`), MC `HYP_d0/runs/lxplus/fd_nd.parquet`
  (prawda truth 0/1, `hypo_swap`). Zapach hipotezy q = −ładunek toru K (`d_charge0`).
- **D⁺**: `SWAP_d0/runs/lxplus/data_428770_dp.parquet` (30 % lumi, model Dplus_H, wstępnie score ≥ 0,59), weto D* (`dstar_veto==0`).
  MC: `FEAT_charmnet/frames/fd_nd_Dplus.pkl` (57k prawdziwych D⁺, weto D*) z wynikiem XGB H z FEAT (te same 18 cech, OOF) i **własnymi** progami WP70;
  kontrola: prawdziwe D⁺ z `SWAP_d0/runs/lxplus/mb_nd_dp.parquet` z modelem analizy (tabela na końcu: μ zgodne w 0,5 MeV, σ_eff w ~1–2 MeV przy małej statystyce).
- **D*-tag**: `HYP_d0/runs/lxplus_dstar` (pełna lumi), wiersze `hypo_from_tag==0`, bez konfliktu, WP D0_H; zapach z ładunku π_s.
  (1) Δm w 139,6–160 MeV dla |m(Kπ)−m₀| < 40 MeV; (2) m(Kπ) dla |Δm − 145,43| < 1,2 MeV.
- **Uwaga (FEAT frames)**: ramki `FEAT_charmnet/frames/data_428770_*` to tylko 2/39 plików (~5 % lumi) — za mało (sygnał D⁰ WP70 ~300); użyto większych wyciągów jw.
- **Fit**: binowany rozszerzony Poisson NLL (iminuit). Sygnał = 2 Gaussy o wspólnym μ; kształt (σ₁, σ₂, f) z MC w danym binie (oba ładunki);
  w danych wolne μ i skala k (σᵢ→kσᵢ). Tło: Czebyszew 2. st. (warianty: 3. st., exp·(1+bx), węższe okno, bin 1 MeV). Okna: D⁰ 1,76–1,98, D⁺ 1,74–2,04 GeV.
  MC: ten sam model (bez tła) osobno dla q=±1 → μ_MC(q), σ_MC(q). Δm: 2 Gaussy + u^a(1+c₁u), u = Δm − m_π; MC z wolnym kształtem i tłem (miękkie π kombinatoryczne).
  Szerokość = σ_eff = √(fσ₁² + (1−f)σ₂²). Błąd σ_dane/σ_MC — tylko z danych.
- **Refleksja D⁰ K↔π**: (i) nominalnie inkluzywnie weto |m_swap − m₀| > 25 MeV w danych i MC (jak w FEAT/charmbkg; usuwa prawie całe odbicia, bo dla
  odbicia m_swap = poprawna masa ≈ m₀); (ii) wariant keep_both: składnik refleksji o kształcie z MC (`hypo_swap==1`, ten sam q-label, bo odbicie D⁰
  ląduje jako „D̄⁰”), normalizacja R = N_refl/N_sig w oknie z MC (≈ 0,53–0,58); (iii) **D*-tag: hipoteza z ładunku π_s → refleksji brak** — to jest
  referencja D⁰. Inkluzywne D⁰ przy 30 % lumi ma S/B ≈ 2–3 % → fity 3–5 GeV i |η|>1 zjeżdżają na granicę k = 0,4 (oznaczone ⚠, niewliczane).

## 2. Tabele (WP70). Δμ = μ_dane − μ_MC dla danego q; (D−D̄) = μ(q=+1) − μ(q=−1); błędy stat. (dane ⊕ MC)
### D⁺→Kππ, m(Kππ) — WP70  (Δμ = μ_dane − μ_MC [MeV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [MeV] | σ_MC [MeV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 | 1341±216 | 1869.95 | 25.75 | -0.09±2.57 | -1.29±3.40 | -0.59±2.07 | +1.65±4.25 | +0.46±0.30 | 0.83±0.14 | 0.89±0.22 | 0.85±0.12 | 0.63 / 9 % |
| pT5-8 | 1020±114 | 1870.07 | 25.90 | +0.64±2.34 | -0.36±2.09 | +0.09±1.58 | +1.07±3.11 | +0.07±0.38 | 1.05±0.14 | 0.82±0.11 | 0.93±0.09 | 0.84 / 15 % |
| pT8-30 | 1925±125 | 1870.60 | 26.15 | -0.70±1.37 | -1.52±1.25 | -1.13±0.93 | +0.06±1.76 | -0.75±0.58 | 1.08±0.10 | 0.95±0.08 | 1.01±0.06 | 0.37 / 6 % |
| pT3-30 | 4204±289 | 1870.08 | 25.86 | +0.05±1.27 | -0.75±1.30 | -0.35±0.91 | +0.97±1.81 | +0.17±0.22 | 0.95±0.08 | 0.89±0.08 | 0.92±0.06 | 0.31 / 6 % |
| abseta<1 | 2406±171 | 1870.20 | 20.13 | -1.21±1.29 | +1.31±1.33 | +0.03±0.93 | -2.39±1.83 | +0.13±0.23 | 0.95±0.09 | 1.01±0.09 | 0.98±0.06 | nan / nan % |
| abseta>1 | 1664±241 | 1869.91 | 31.88 | +2.66±2.65 | -5.52±2.75 | -1.21±2.07 | +8.42±3.80 | +0.24±0.43 | 0.85±0.12 | 0.64±0.14 | 0.81±0.10 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = -0.79 ± 0.74 MeV (χ²=0.5/2); (D−D̄)_dane − (D−D̄)_MC = +0.90 ± 1.49 MeV; σ_d/σ_MC = 0.966 ± 0.048 (χ²=1.7/2).

### D⁰→Kπ inkluzywnie, weto refl. — WP70  (Δμ = μ_dane − μ_MC [MeV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [MeV] | σ_MC [MeV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 ⚠ | 730±157 | 1865.42 | 31.66 | -5.21±2.49 | +4.73±5.70 | -3.87±1.92 | -9.95±6.21 | -0.00±0.20 | 0.40±0.12 | 0.94±0.32 | 0.40±0.06 | 0.15 / 0 % |
| pT5-8 | 668±173 | 1866.02 | 32.03 | +3.68±4.89 | -3.13±3.42 | -0.03±2.53 | +6.60±5.96 | -0.20±0.28 | 0.75±0.28 | 0.79±0.20 | 0.70±0.14 | 0.34 / 15 % |
| pT8-30 | 1614±271 | 1866.41 | 33.77 | -3.32±2.06 | -0.80±2.35 | -2.22±1.50 | -1.80±3.09 | +0.72±0.45 | 0.88±0.14 | 1.01±0.21 | 0.92±0.11 | 0.19 / 8 % |
| pT3-30 | 3068±595 | 1865.72 | 32.01 | -2.07±2.20 | -0.05±2.31 | -1.13±1.59 | -1.99±3.19 | +0.03±0.15 | 0.70±0.17 | 0.80±0.19 | 0.74±0.12 | 0.38 / 8 % |
| abseta<1 | 2344±349 | 1865.69 | 23.39 | +0.40±2.04 | +0.55±2.17 | +0.48±1.49 | -0.17±2.97 | -0.02±0.16 | 0.92±0.14 | 0.92±0.17 | 0.92±0.11 | nan / nan % |
| abseta>1 ⚠ | 491±151 | 1865.98 | 39.74 | -7.61±4.86 | -3.54±4.65 | -5.74±3.38 | -3.94±6.72 | +0.12±0.31 | 0.40±0.11 | 0.40±0.18 | 0.40±0.06 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = -1.65 ± 1.29 MeV (χ²=0.6/1); (D−D̄)_dane − (D−D̄)_MC = -0.51 ± 2.77 MeV; σ_d/σ_MC = 0.829 ± 0.088 (χ²=1.5/1).

### D⁰ keep_both + szablon refl. z MC — WP70  (Δμ = μ_dane − μ_MC [MeV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [MeV] | σ_MC [MeV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 ⚠ | 996±167 | 1865.28 | 31.17 | -4.51±2.10 | +2.09±3.87 | -3.09±1.61 | -6.55±4.40 | +0.05±0.18 | 0.40±0.15 | 0.84±0.28 | 0.40±0.12 | nan / nan % |
| pT5-8 | 900±197 | 1865.89 | 31.67 | +4.84±3.83 | -4.12±2.91 | +0.01±2.17 | +8.61±4.81 | -0.35±0.26 | 0.74±0.21 | 0.83±0.17 | 0.73±0.13 | nan / nan % |
| pT8-30 | 1957±278 | 1866.25 | 33.55 | -4.14±1.98 | -2.29±1.97 | -3.07±1.37 | -1.15±2.77 | +0.70±0.42 | 0.93±0.12 | 0.88±0.15 | 0.92±0.10 | nan / nan % |
| pT3-30 | 4320±725 | 1865.58 | 31.60 | -1.28±1.90 | -1.05±1.85 | -1.13±1.33 | -0.22±2.65 | +0.00±0.14 | 0.75±0.13 | 0.85±0.14 | 0.79±0.11 | nan / nan % |
| abseta<1 | 2931±348 | 1865.52 | 23.18 | +0.44±1.71 | -0.16±1.81 | +0.17±1.24 | +0.65±2.48 | +0.05±0.15 | 0.88±0.11 | 0.95±0.14 | 0.90±0.09 | nan / nan % |
| abseta>1 ⚠ | 672±162 | 1865.86 | 39.40 | -6.39±4.62 | -3.65±5.06 | -4.87±3.02 | -2.78±6.85 | -0.04±0.28 | 0.40±0.18 | 0.63±0.40 | 0.40±0.18 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = -2.19 ± 1.16 MeV (χ²=1.4/1); (D−D̄)_dane − (D−D̄)_MC = +0.88 ± 2.42 MeV; σ_d/σ_MC = 0.849 ± 0.077 (χ²=1.5/1).

### D⁰ z tagu D* (π_s → zapach, bez refl.) — WP70  (Δμ = μ_dane − μ_MC [MeV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [MeV] | σ_MC [MeV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 | 842±176 | 1865.41 | 28.29 | +0.99±2.68 | +2.18±2.30 | +2.04±1.89 | -0.95±3.52 | +0.24±0.36 | 1.42±0.23 | 0.77±0.15 | 1.08±0.16 | 0.46 / 7 % |
| pT5-8 | 622±80 | 1865.72 | 28.34 | -0.77±1.66 | +0.69±1.90 | -0.34±1.33 | -2.02±2.48 | -0.56±0.46 | 0.70±0.12 | 1.15±0.14 | 0.97±0.11 | 0.23 / 8 % |
| pT8-30 | 1134±75 | 1865.96 | 29.59 | -1.18±1.30 | +0.23±1.35 | -0.51±0.94 | -0.63±1.73 | +0.78±0.74 | 0.87±0.07 | 0.98±0.09 | 0.92±0.06 | 0.54 / 7 % |
| pT3-30 | 2626±207 | 1865.59 | 28.47 | -0.12±1.09 | +0.99±1.08 | +0.51±0.77 | -1.05±1.52 | +0.07±0.27 | 1.08±0.09 | 0.95±0.08 | 1.01±0.06 | 0.23 / 3 % |
| abseta<1 | 1644±108 | 1865.52 | 22.43 | +1.78±1.09 | +0.61±1.06 | +1.19±0.76 | +1.21±1.49 | +0.05±0.28 | 1.06±0.09 | 1.01±0.08 | 1.03±0.06 | nan / nan % |
| abseta>1 | 919±253 | 1865.90 | 36.04 | -5.79±2.87 | +2.93±3.22 | -1.70±2.18 | -8.63±4.27 | +0.08±0.57 | 0.96±0.21 | 0.92±0.24 | 0.97±0.17 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = -0.10 ± 0.71 MeV (χ²=1.5/2); (D−D̄)_dane − (D−D̄)_MC = -1.39 ± 1.38 MeV; σ_d/σ_MC = 0.947 ± 0.047 (χ²=0.9/2).

### Δm(D*) = m(Kππ_s) − m(Kπ) — WP70  (Δμ = μ_dane − μ_MC [keV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [keV] | σ_MC [keV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 | 806±75 | 145544.81 | 947.72 | -325.88±91.87 | -95.58±79.41 | -211.05±63.20 | -245.38±120.67 | -15.07±13.61 | 1.02±0.14 | 0.79±0.12 | 0.95±0.09 | 34.24 / 6 % |
| pT5-8 | 672±40 | 145536.69 | 904.34 | -94.86±69.68 | -119.72±62.64 | -109.59±46.29 | +20.97±92.22 | -3.88±16.58 | 1.04±0.10 | 1.03±0.10 | 1.03±0.07 | 18.24 / 2 % |
| pT8-30 | 1349±44 | 145527.62 | 1010.90 | -5.09±44.46 | +24.21±45.12 | +10.08±31.74 | -21.56±58.03 | +7.74±25.39 | 1.00±0.05 | 1.11±0.05 | 1.05±0.04 | 10.72 / 5 % |
| pT3-30 | 2803±99 | 145539.67 | 936.22 | -130.82±39.00 | -73.07±37.57 | -101.39±27.09 | -64.84±53.27 | -7.10±9.73 | 1.04±0.06 | 1.00±0.05 | 1.02±0.04 | 4.19 / 4 % |
| abseta<1 | 1684±68 | 145533.10 | 806.93 | -123.75±41.48 | -32.17±39.34 | -75.82±28.52 | -92.15±56.15 | -0.58±10.75 | 1.03±0.06 | 0.97±0.06 | 1.00±0.04 | nan / nan % |
| abseta>1 | 1142±74 | 145561.36 | 1291.01 | -145.11±82.30 | -181.52±84.44 | -163.54±59.00 | +22.84±116.33 | -13.56±19.27 | 0.99±0.09 | 1.03±0.09 | 1.01±0.07 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = -54.97 ± 24.19 keV (χ²=11.7/2); (D−D̄)_dane − (D−D̄)_MC = -46.62 ± 48.17 keV; σ_d/σ_MC = 1.038 ± 0.030 (χ²=1.1/2).

## 3. Weak-mode (sagitta) z m(D) − m(D̄)
Czułość z MC (prawdziwe kandydaty WP70, `sagitta.py`): q/pT → q/pT + δ dla każdego toru (K: −q, π: +q):
D⁺: dm/dδ = −0,93 / −1,37 / −2,50 MeV na TeV⁻¹ (pT 3–5/5–8/8–30; −1,29 inkluzywnie), m(D⁺)−m(D⁻) = 2·dm/dδ·δ.
D⁰: −0,09 / −0,08 / +0,11 MeV/TeV⁻¹ — praktycznie nieczuły (kasowanie przeciwnych ładunków). Skala pędu (symetryczna) 1 ‰ → +1,5 MeV (D⁰), +1,35 MeV (D⁺);
zmierzone Δμ = −0,43 ± 0,51 MeV ⇒ skala pędu dane/MC = (−0,3 ± 0,4) ‰.

| D⁺ (WP70) | [(D⁺−D⁻)_dane − MC] [MeV] | δ_sagitta [TeV⁻¹] |
|---|---|---|
| pT 3–5 | +1,2 ± 4,3 | −0,6 ± 2,3 |
| pT 5–8 | +1,0 ± 3,1 | −0,4 ± 1,1 |
| pT 8–30 | +0,8 ± 1,9 | −0,16 ± 0,37 |
| pT 3–30 | +0,8 ± 1,8 | **−0,31 ± 0,70** |
| \|η\| < 1 | −2,5 ± 1,8 | +0,99 ± 0,72 |
| \|η\| > 1 | +8,2 ± 3,8 (k=1: +7,5 ± 4,5) | −3,1 ± 1,5 |

- Inkluzywnie δ = −0,3 ± 0,7 TeV⁻¹ (zgodne z 0; 8–30 GeV najczulszy: |δ| < 0,9 TeV⁻¹ 95 %).
- Podział |η|: centralny +1,4σ, przedni −2,1σ, różnica ~2,1σ; przedni wynik stabilny wzgl. tła/okna (±1 MeV) i przy ustalonym k, obie strony η
  ten sam znak (η<−1: +4,9 ± 6,0; η>1: +10,9 ± 6,6 MeV). Przy ~8 testowanych podziałach η to nie jest istotne, ale **jedyny ślad** — warto powtórzyć
  na pełnej lumi D⁺ (obecnie 30 %; błąd ×0,55).
- D⁰ (tag D*): (D⁰−D̄⁰)_dane − MC = −1,1 ± 1,5 MeV (pT 3–30), |η|>1: −8,7 ± 4,3 MeV — przy czułości 0,1 MeV/TeV⁻¹ wymagałoby δ ~ 80 TeV⁻¹,
  więc to fluktuacja / efekt nie-sagittowy (np. dopasowanie tła w małej próbce), nie weak-mode. Δm(D*⁺)−Δm(D*⁻) − MC = −47 ± 48 keV — zgodne z 0.

## 4. Stabilność (WP50) i kontrola MC
### D⁺→Kππ, m(Kππ) — WP50  (Δμ = μ_dane − μ_MC [MeV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [MeV] | σ_MC [MeV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 | 955±140 | 1869.69 | 25.29 | -1.14±2.76 | +0.50±2.45 | -0.42±1.84 | -1.33±3.67 | +0.30±0.35 | 1.09±0.21 | 0.90±0.17 | 0.99±0.14 | nan / nan % |
| pT5-8 | 830±87 | 1869.97 | 25.64 | +1.71±2.15 | +0.16±2.38 | +1.00±1.59 | +1.90±3.18 | +0.35±0.44 | 1.13±0.14 | 1.15±0.16 | 1.15±0.11 | nan / nan % |
| pT8-30 | 1254±78 | 1870.42 | 25.73 | +0.22±1.46 | -0.42±1.30 | -0.10±0.97 | +0.11±1.83 | -0.52±0.67 | 1.17±0.10 | 0.96±0.08 | 1.05±0.06 | nan / nan % |
| pT3-30 | 3030±185 | 1869.89 | 25.47 | +0.35±1.22 | +0.16±1.16 | +0.24±0.84 | +0.39±1.67 | +0.21±0.25 | 1.13±0.09 | 0.99±0.08 | 1.06±0.06 | nan / nan % |
| abseta<1 | 1747±107 | 1870.12 | 20.11 | -0.01±1.33 | +0.49±1.19 | +0.23±0.89 | -0.31±1.76 | +0.18±0.27 | 1.17±0.10 | 1.09±0.08 | 1.13±0.06 | nan / nan % |
| abseta>1 | 1267±169 | 1869.45 | 31.43 | +0.32±2.55 | -1.15±2.97 | -0.35±1.94 | +1.78±3.88 | +0.31±0.51 | 0.98±0.13 | 0.97±0.20 | 0.98±0.11 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = +0.09 ± 0.76 MeV (χ²=0.4/2); (D−D̄)_dane − (D−D̄)_MC = +0.46 ± 1.52 MeV; σ_d/σ_MC = 1.064 ± 0.051 (χ²=1.0/2).

### D⁰→Kπ inkluzywnie, weto refl. — WP50  (Δμ = μ_dane − μ_MC [MeV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [MeV] | σ_MC [MeV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 | 640±227 | 1865.14 | 31.36 | +7.37±7.29 | +6.21±3.33 | +7.02±3.26 | +1.15±8.01 | +0.00±0.24 | 0.96±0.39 | 0.61±0.25 | 0.82±0.21 | nan / nan % |
| pT5-8 | 469±93 | 1865.79 | 31.57 | +0.28±2.64 | +3.51±2.80 | +0.39±1.86 | -3.26±3.84 | -0.03±0.33 | 0.80±0.18 | 0.44±0.18 | 0.67±0.11 | nan / nan % |
| pT8-30 | 1075±159 | 1866.24 | 33.37 | -0.79±2.11 | -0.12±1.89 | -0.51±1.41 | -0.12±2.78 | +0.55±0.53 | 0.88±0.14 | 0.92±0.14 | 0.90±0.10 | nan / nan % |
| pT3-30 | 2080±335 | 1865.47 | 31.67 | +1.87±2.11 | +3.17±1.46 | +2.35±1.30 | -1.24±2.56 | +0.06±0.18 | 0.91±0.14 | 0.61±0.14 | 0.81±0.10 | nan / nan % |
| abseta<1 | 1307±160 | 1865.49 | 23.39 | +1.77±1.81 | +2.28±1.37 | +1.94±1.17 | -0.57±2.26 | -0.05±0.19 | 0.98±0.12 | 0.68±0.10 | 0.86±0.09 | nan / nan % |
| abseta>1 ⚠ | 674±469 | 1865.62 | 39.51 | -4.15±9.05 | +2.90±5.87 | +1.85±5.52 | -6.75±10.78 | +0.30±0.37 | 2.49±1.86 | 0.72±0.31 | 0.80±0.35 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = +0.58 ± 1.06 MeV (χ²=4.5/2); (D−D̄)_dane − (D−D̄)_MC = -1.37 ± 2.19 MeV; σ_d/σ_MC = 0.799 ± 0.072 (χ²=2.3/2).

### D⁰ z tagu D* (π_s → zapach, bez refl.) — WP50  (Δμ = μ_dane − μ_MC [MeV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [MeV] | σ_MC [MeV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 | 524±79 | 1865.23 | 28.13 | -1.71±2.42 | +1.60±2.29 | +0.38±1.61 | -3.10±3.31 | +0.21±0.42 | 1.10±0.20 | 0.88±0.17 | 0.97±0.12 | 0.86 / 10 % |
| pT5-8 | 535±63 | 1865.62 | 28.33 | -1.17±1.96 | -0.66±1.86 | -0.66±1.35 | -1.46±2.65 | -0.94±0.54 | 0.90±0.18 | 1.16±0.15 | 1.07±0.11 | 0.58 / 8 % |
| pT8-30 ⚠ | 805±55 | 1865.70 | 29.11 | -0.70±1.44 | +1.90±1.66 | +0.46±1.07 | -1.97±2.03 | +0.62±0.85 | 0.86±0.07 | 1.06±0.13 | 0.95±0.06 | 0.67 / 4 % |
| pT3-30 | 1877±120 | 1865.43 | 28.33 | -0.85±1.03 | +1.10±1.09 | +0.11±0.75 | -2.05±1.47 | -0.10±0.31 | 1.00±0.07 | 1.02±0.08 | 1.01±0.06 | 0.45 / 3 % |
| abseta<1 | 1254±68 | 1865.44 | 22.37 | -0.00±1.02 | -0.47±1.08 | -0.23±0.74 | +0.50±1.45 | +0.03±0.33 | 1.06±0.07 | 1.11±0.08 | 1.08±0.05 | nan / nan % |
| abseta>1 | 564±128 | 1865.49 | 36.12 | -4.00±2.88 | +6.25±2.91 | +1.00±2.14 | -10.68±4.04 | -0.43±0.68 | 0.79±0.16 | 1.01±0.29 | 0.91±0.15 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = -0.23 ± 1.04 MeV (χ²=0.2/1); (D−D̄)_dane − (D−D̄)_MC = -1.62 ± 2.10 MeV; σ_d/σ_MC = 1.027 ± 0.082 (χ²=0.3/1).

### Δm(D*) = m(Kππ_s) − m(Kπ) — WP50  (Δμ = μ_dane − μ_MC [keV]; σ = σ_eff 2 Gaussów)

| bin | N_S dane | μ_MC [keV] | σ_MC [keV] | Δμ D | Δμ D̄ | Δμ D+D̄ | (D−D̄) dane | (D−D̄) MC | σ_d/σ_MC D | σ_d/σ_MC D̄ | σ_d/σ_MC D+D̄ | rozrzut modelu Δμ / σ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pT3-5 | 611±42 | 145539.67 | 933.63 | -130.27±81.10 | -26.41±71.36 | -77.64±54.47 | -120.14±106.85 | -16.28±15.90 | 1.10±0.11 | 0.90±0.10 | 1.01±0.08 | 13.11 / 3 % |
| pT5-8 | 513±29 | 145532.86 | 896.60 | -1.92±72.11 | -70.66±60.53 | -37.31±46.60 | +58.45±92.11 | -10.29±19.49 | 1.11±0.11 | 0.96±0.09 | 1.03±0.07 | 10.65 / 3 % |
| pT8-30 | 931±35 | 145520.55 | 991.07 | +16.64±51.69 | +23.26±53.83 | +21.95±37.39 | -14.05±68.47 | -7.42±29.69 | 1.01±0.06 | 1.12±0.06 | 1.07±0.04 | 6.99 / 6 % |
| pT3-30 | 2044±63 | 145534.61 | 927.45 | -38.50±37.39 | -18.21±36.21 | -28.38±26.03 | -31.81±50.79 | -11.52±11.40 | 1.07±0.05 | 1.02±0.05 | 1.05±0.04 | 6.46 / 4 % |
| abseta<1 | 1248±45 | 145530.50 | 800.14 | -65.47±40.21 | +8.38±39.35 | -27.71±28.15 | -70.78±54.85 | +3.07±12.54 | 1.07±0.06 | 1.01±0.05 | 1.05±0.04 | nan / nan % |
| abseta>1 | 831±46 | 145552.21 | 1303.44 | +27.09±78.91 | -68.03±78.35 | -20.56±55.57 | +58.79±108.83 | -36.33±22.87 | 1.02±0.08 | 1.05±0.09 | 1.03±0.06 | nan / nan % |

Średnia ważona binów pT (tylko fity bez ⚠): Δμ = -18.28 ± 25.71 keV (χ²=2.5/2); (D−D̄)_dane − (D−D̄)_MC = -6.18 ± 51.43 keV; σ_d/σ_MC = 1.050 ± 0.032 (χ²=0.5/2).
### Kontrola MC D⁺: fd (model FEAT H, własne progi) vs mb_nd (model analizy Dplus_H), WP70

| próbka | pT | N | μ [MeV] | σ_eff [MeV] |
|---|---|---|---|---|
| mc | 3-5 | 21502 | 1869.95±0.15 | 25.75 |
| mc | 5-8 | 12936 | 1870.07±0.19 | 25.90 |
| mc | 8-30 | 5758 | 1870.60±0.29 | 26.15 |
| mbmc | 3-5 | 2761 | 1870.44±0.41 | 26.00 |
| mbmc | 5-8 | 1590 | 1869.63±0.55 | 25.67 |
| mbmc | 8-30 | 732 | 1870.24±0.90 | 28.34 |

## 5. Ograniczenia
- Dane D⁰/D⁺ inkluzywnie = 30 % lumi run 428770; D* = pełna lumi. Dokładność Δμ ~0,7 MeV na kanał — placeholderu 0,2 MeV dla m(D) nie da się zweryfikować bez pełnej lumi
  (pełna lumi D⁺: ~0,4 MeV; łącznie z D⁰-tag ~0,35 MeV).
- MC D⁺ z innym (FEAT) modelem XGB i własnymi progami WP70; kontrola na mb_nd z modelem analizy zgodna, ale mała statystyka.
- Błędy statystyczne z Hesse; błąd σ_dane/σ_MC bez błędu MC (MC ~10× większe). Rozrzut modelu tła podany osobno (nie dodany).
- Inkluzywne D⁰ (S/B ~2–3 %) daje niestabilne fity w 3–5 GeV i |η|>1 — nie używać go do tych systematyk.

## Pliki
`prep.py`, `prep2.py`, `prep_dstar.py` (wyciągi → `out/*.parquet`), `common.py`, `fitlib.py`, `fit_peaks.py` (→ `out/peaks_rows.csv`),
`fit_dstar.py` (→ `out/dstar_rows.csv`), `sagitta.py` (→ `out/sagitta.json`), `check_eta_kfix.py` (→ `out/check_eta_kfix.txt`),
`summarize.py` (→ `out/tables.md`, `out/summary.json`), `look.py` (→ `out/look_data.png`). Python: `.venv_h5/bin/python`. Nic poza `PEAK_pos/` nie zmieniono.
