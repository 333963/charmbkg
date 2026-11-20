# D23/ANALIZA — pierwszy pomiar D⁺→Kππ w binie pT 2–3 GeV (projekt B), 18.11.2026

Dane inkluzywne, bez AFP; runy osobno (427929: GRL ∧ HLT_mb_mbts, 82 % kandydatów; 428770: wszystko, 74 % z L1J12). Klasyfikator H (`H_s12345`, trening pT ≥ 3 GeV) użyty poza zakresem treningu.
Kryteria: `KRYTERIA.md` (zapisane przed otwarciem plików). Tabele: `tables.md` (T1–T5). Liczby: `out/{wp,fit,verdict,d0_dstar,retrain}.json`. Liczenie lokalne ~2 min.

## 0. Werdykt wg K1–K3
| run | WP (z MC + pasm) | N_S | Z = N/σ | S/B (±45 MeV) | K1 | Δμ dane − MC [MeV] | σ_dane/σ_MC | K3: y = N/ε vs sąsiedzi |
|---|---|---|---|---|---|---|---|---|
| 427929 | WP20 (score ≥ 0,992; ε_S^MC 0,205) | 705 ± 53 | **13,2** | 0,48 | TAK | −0,9 ± 1,4 | 0,99 ± 0,09 | stabilny (WP30 −1,0σ; WP10 +1,7σ); y = 3441 ± 260 |
| 428770 | WP20 (ten sam próg) | 626 ± 74 | **8,4** | 0,19 | TAK | +0,9 ± 2,6 | 1,01 ± 0,21 | stabilny (WP30 +0,9σ; WP10 0,0σ); y = 3055 ± 363 |

Pik D⁺ w 2–3 GeV jest zmierzony w obu runach. σ(N) = max(HESSE, rms bootstrapu po zdarzeniach, 300 replik) — obie zgodne w 3 %. Rozrzut tła pol2/pol4: 427929 686/706 (−3 %/0 %), 428770 667/606 (+7 %/−3 %).
Błędy Δμ, σ: bootstrap ⊕ błąd MC (0,33 MeV; 1,8 %). Dla 428770 błąd z wszystkich replik (post hoc, §7), bo 16 % replik trafia w granice fitu.

## 1. Kontrola wejść
- Wiersze = `summary.json` (2 680 368 / 1 547 360 / 243 118), pliki 96 / 450 / 15, pT ∈ [2,000; 3,000) we wszystkich próbkach, score bez NaN.
- Σw_ps vs histogram pełny: pasma (score < 0,93) 0,9934 (427929) i 0,9991 (428770); okno 1,0027 i 1,0023. Odchylenie 427929 (−0,66 %) to ~4σ przy liczeniu wierszy jako niezależnych —
  preskala jest po zdarzeniach (grona kandydatów), więc nie traktuję tego jako błędu; fity nie używają wierszy preskalowanych. Score ≥ 0,93: wiersze == histogram co do binu (118 734 i 253 175; w_ps = 1).
- Sygnał MC fd_nd: 27 709 po wecie D* (weto usuwa 180 = 0,65 %), wagi = 1, wiersze == histogram.
- **Strata klasyfikatora poza zakresem treningu** — ε_S^MC w 2–3 GeV przy progach WP binu 3–5:

| WP binu 3–5 (próg) | 95 (0,291) | 90 (0,572) | 80 (0,840) | 70 (0,933) | 50 (0,983) |
|---|---|---|---|---|---|
| ε_S 2–3 GeV | 0,867 | 0,769 | 0,619 | 0,495 | 0,298 |
| 2,0–2,5 / 2,5–3,0 | 0,82 / 0,91 | 0,71 / 0,83 | 0,54 / 0,69 | 0,41 / 0,57 | 0,22 / 0,37 |

  Próg „WP70” z 3–5 GeV to w 2–3 GeV WP50. `score_avg` daje to samo w 0,6 %. Kształt sygnału MC: μ 1869,4–1870,5 MeV, σ_eff 23,1–24,6 MeV (3–5 GeV: 23,6).

## 2. Punkt pracy (MC + pasma boczne, okno ±60 MeV zaślepione) — T1
Z_oczek rośnie monotonicznie od WP90 (2,4 / 2,2) do WP20 (11,5 / 11,1) i spada przy „WP10” (9,1 / 9,8) → **WP20 w obu runach**. Przy WP20: B(±45) = 1374 / 3110, S_oczek = 499 / 682.
Uwaga: punkt „WP10” po zaokrągleniu progu do krawędzi histogramu (0,998) ma ε = 0,067, nie 0,10 (reguła z KRYTERIA §1, bin score 0,002).

## 3. Fity — T2, T3
- Zaobserwowane vs oczekiwane przy WP20: N = 705 vs 533 (427929), 626 vs 728 (428770) — w granicach czynnika 1,5 ekstrapolacji.
- **Stosunek uzysków do 3–5 GeV**: y(2–3)/y(3–5, WP70) = **0,79 ± 0,08** (427929), **0,52 ± 0,07** (428770); oczekiwane 0,6 (×/÷ 1,5). Bez przekroju: ε_S^MC z fd_nd, akceptancja/rekonstrukcja nieuwzględnione;
  zestawy plików 2–3 i 3–5 GeV nie są identyczne (427929: 450 vs 458 plików; 428770: 96/97).
- Przy luźniejszych cięciach (błąd tylko HESSE, bez bootstrapu): 427929 Z = 3,2 / 4,5 / 6,0 / 8,0 / 10,0 dla WP90…50 (S/B 0,004 → 0,077); 428770 Z = 1,0 / 2,5 / 3,4 / 3,9 / 3,4 (S/B 0,001 → 0,017).
  Prognoza z planu („Z ≈ 4–5, S/B 0,01–0,025”) odpowiada cięciu ε ≈ 0,7; przy optymalnym (ciaśniejszym) cięciu jest ~2–3× lepiej.
- Pozycja/szerokość przy innych WP (427929): Δμ = +1,2 ± 1,5 (WP30), +1,8 ± 1,5 (WP40) MeV; σ_d/σ_MC = 0,97 ± 0,11, 0,92 ± 0,13. WP50 i WP10 — fit z wolnym kształtem niestabilny (55 % i 74 % udanych replik).
  428770: fit z wolnym kształtem stabilny tylko przy WP20 (84 %); przy WP50 na granicy.
- Domknięcie fitu na fd_nd (sygnał + tło filtrowanego MC; S/B dużo większe niż w danych): N_fit/N_prawda = 0,992 ± 0,017 (WP20), 1,00–1,01 (WP50–30), 1,035–1,044 (WP70–90).
  Klasa „inne” (truth 2–7) pod pikiem: 17 kandydatów na 5678 sygnału przy WP20, bez piku. Obciążenia fitu przy S/B danych nie da się sprawdzić — MC minbias 2–3 GeV nie istnieje.

## 4. Skan cięcia — T4
| run | ρ(20/30) | ρ(20/40) | ρ(20/50) | ρ(20/10) | czułość 2σ |
|---|---|---|---|---|---|
| 427929 | 0,943 ± 0,056 | 0,907 ± 0,074 | 0,894 ± 0,095 | 1,168 ± 0,118 | 12 / 16 / 21 / 20 % |
| 428770 | 1,100 ± 0,126 | 1,142 ± 0,182 | 1,67 ± 0,64 | 1,001 ± 0,142 | 23 / 32 / 77 / 28 % |

- 427929: y spada monotonicznie 3849 → 3794 → 3647 → 3441 → 2945 (WP50 → WP10); każda para w 1,7σ, ale kierunek stały (sugestia ε_S^dane < ε_S^MC przy bardzo ciasnym cięciu, nieistotna).
- 428770: y(WP50) = 1831 ± 543 odstaje (2,2σ od WP20; S/B 0,017) — luźny koniec w runie dżetowym jest nieużyteczny.
- **Czy da się coś powiedzieć o ε_S dane/MC: NIE jako pomiar.** Skan widzi tylko zmianę względną między WP; czułość 2σ = 12–21 % (427929) i 23–77 % (428770), a punkt odniesienia (WP50) ma S/B < 0,08.
  Wartość bezwzględna SF niewidoczna; brak tagu D* (§5). Wynik: „brak zależności SF od WP większej niż ~20 % (2σ) w 427929”, nic więcej.

## 5. D⁰ inkluzywne i D* (histogramy pełnej statystyki, bez klasyfikatora) — `out/d0_dstar.json`
- **D⁰→Kπ inkluzywne**: 34,7 M / 63,7 M kandydatów; Z = 1,9–3,1 (427929) i −0,3…1,9 (428770) zależnie od rzędu tła (pol3–5) → **piku nie widać** (Z < 5; wartości nie podaję).
- **D\***: Δm przy |m(D⁰) − 1,865| < 40 MeV, tylko RS (WS w produkcji brak): Z = 1,5–3,9 (427929), 1,9–5,4 (428770) zależnie od modelu tła progowego;
  jedyne Z > 5 pochodzi z najsztywniejszego tła o złym χ² (130/96), przy elastyczniejszym Z ≈ 2 → **nie spełnia Z ≥ 5 w sposób odporny; wartości nie podaję**. S/B ~ 0,5–1,7 %.

## 6. Wnioski dla fitu pojedynczego D w 2–3 GeV
1. **PEAK_POSITION / PEAK_WIDTH dla D⁺ — placeholdery da się zastąpić pomiarem**: Δμ = −0,9 ± 1,4 MeV (427929), +0,9 ± 2,6 MeV (428770); średnia ważona **−0,5 ± 1,3 MeV**; σ_d/σ_MC = 0,99 ± 0,09 i 1,01 ± 0,21, średnia **1,00 ± 0,08**.
   Zgodne z zerem/jedynką i z bin 3–5 GeV (PEAK_pos: −0,6 ± 2,1 MeV; 0,85 ± 0,12). Dokładność danych nie uzasadnia ±0,2 MeV; uczciwe liczby z danych: **±1,3 MeV i ±8–10 %** (placeholder ±10 % szerokości w porządku).
   Dla D⁰ i Δm(D*) w 2–3 GeV pomiaru nie ma.
2. **S/B**, jakie go czeka (D⁺, ±45 MeV): bez klasyfikatora ≲ 0,002; ε_S ≈ 0,7: 0,019 (427929) / 0,008 (428770); ε_S ≈ 0,5: 0,077 / 0,017; ε_S ≈ 0,2: 0,48 / 0,19. Run dżetowy ma S/B ~2,5–4× gorsze.
   Istotność maksymalna przy ε_S ≈ 0,2–0,3 (płaskie maksimum: Z = 12,6 przy WP30, 13,2 przy WP20 w 427929). Tło exp(pol3); pol2 przy luźnych cięciach przesuwa N o 12–40 % (WP80–90).
3. **Czego nie da się zmierzyć i dlaczego**: (a) ε_S dane/MC bezwzględnie — brak próbki tagowanej (D* niewidoczne), skan daje tylko zmianę względną z czułością ≥ 12 %; (b) odbicie K↔π D⁰ — D⁰ inkluzywne niewidoczne bez klasyfikatora,
   tag D* (π_s pod progiem pT) nie daje piku; (c) wydajność weta D* na sygnale D⁺ — tylko z MC (0,65 % sygnału usunięte), w danych brak czułości; (d) obciążenie fitu przy S/B danych — brak MC minbias 2–3 GeV.
4. **Klasyfikator trenowany poniżej 3 GeV** (`s05_retrain.py`, szybki XGB, 18 cech, CV 5-krotne po zdarzeniach, ocena na pasmach bocznych 427929): przy stałym ε_S tło spada o 10 / 7 / 6 / 3 % (ε = 0,5 / 0,4 / 0,3 / 0,2)
   dla treningu sygnał MC vs pasma danych (wariant A, może korzystać z różnic dane/MC ⇒ górne ograniczenie) a dla sygnał MC vs tło fd_nd (wariant B) jest o 1–9 % większe niż dla H; Z_max 11,6 vs 11,5; AUC 0,976 / 0,975 vs 0,974.
   **Dotrenowanie nic istotnego nie daje** — H ekstrapoluje dobrze; spadek ε_S przy stałym progu (0,70 → 0,50) to przesunięcie rozkładu score (miększa kinematyka), nie utrata separacji. Zastrzeżenie: 27,7 tys. sygnału, 400 drzew bez strojenia.
   Poboczne: cięcie H zmienia nachylenie tła (stosunek pasmo dolne/górne 0,41 → 1,26 przy WP20; wariant A 1,07) — exp(pol3) to opisuje (χ² 100/95, 86/95).

## 7. Post hoc (po obejrzeniu wyników, poza kryteriami)
- 428770: błąd Δμ i σ z wszystkich 300 replik (2,6 MeV; 0,21) zamiast tylko udanych (2,25; 0,17), bo odrzucanie replik na granicy obcina rozkład. 427929 bez różnicy (1,40 → 1,42).
- Średnie ważone obu runów w §6.1 nie były zapowiedziane w kryteriach.
- Przedziały pT 2–2,5 / 2,5–3 policzone tylko dla ε_S MC i D⁰/D*, nie dla fitu D⁺.

## 8. Zastrzeżenia
Z = N/σ(N), nie test ilorazu wiarygodności. ε_S^MC z fd_nd (MC filtrowane, 15 plików), score poza zakresem treningu; fold MC ≠ fold treningowy (bez znaczenia, pT < 3 nie było w treningu). Brak testu wstrzyknięcia z pełnym tłem.
Bootstrap po zdarzeniach tylko dla WP ≤ 50 (wiersze kompletne); luźniejsze WP: błędy HESSE. Stosunek do 3–5 GeV: rząd wielkości. D⁰/D*: szablon MC + gładkie tło, bez odbicia w modelu. 428770: 96/97 plików. Wykresów nie robiłem.

## 9. Pliki
`KRYTERIA.md`, `RAPORT.md`, `tables.md`, `scripts/{s01_wp.py, s02_fit.py, s03_tables.py, s04_d0_dstar.py, s05_retrain.py}`, `out/{wp.json, fit.json, verdict.json, d0_dstar.json, retrain.json}`, `logs/`.
