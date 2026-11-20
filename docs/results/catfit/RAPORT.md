# CATFIT — fit jednoczesny masy w kategoriach score drzewa H zamiast cięcia WP70 (propozycje B1 / B2), run 427929 (projekt B)

20.11.2026, ~03:25 CET. **Status: PRELIMINARY, w całości POST HOC względem XSEC_final.** Kryteria: `KRYTERIA.md` (zapisane ~02:05, przed wszystkimi fitami tego zadania).
Tylko piki inkluzywne D⁰→K⁻π⁺ i D⁺→K⁻π⁺π⁺, bez tagu protonu AFP. Praca lokalna (bez lxplus, bez condora); nowa ekstrakcja nie była potrzebna (`wps_hist.npz` ma 19 przedziałów score, w tym „poniżej WP95”).
Inne katalogi tylko czytane / importowane; `final.json`, nota — bez zmian.
Tabele: `tabele/catfit.md` (T1–T8), liczby: `out/catfit.json`; surowe: `out/{mc,toys,data}.json`.

## 0. Odpowiedź

**Nie. Fit w kategoriach score nie zmniejsza błędu całkowitego przekroju w żadnym z 6 binów; ze swobodnymi ułamkami zwiększa go 1,5–10 ×.** Werdykt wg kryteriów: 0 / 6 komórek „ZYSK” (F0 i F70), dla obu modeli D⁺ (M0 i A2).

| gat. | bin | błąd obecny: stat / syst / razem [%] | **F0 (ułamki swobodne): stat / syst / razem [%]** | g = F0 / obecny | F70 (jeden więz SF₇₀): razem [%] | g(F70) | obciążenie N_tot na mb_nd [%] | pokrycie F0 (zabawki) | σ: final → F0 [µb] |
|---|---|---|---|---|---|---|---|---|---|
| D⁰ | 3–5 | 11,5 / 26,7 / 29,2 | **30,4 / 53,2 / 61,3** | 2,10 | 28,4–45,7 | 0,97–1,57 | +12,7 ± 41,1 | 0,77 | 1028 ± 119 ± 274 → 1403 ± 426 ± 746 |
| D⁰ | 5–8 | 10,5 / 16,8 / 20,0 | **59,8 / 51,1 / 78,7** | 3,94 | 19,3–34,2 | 0,97–1,71 | +14,9 ± 25,3 | 0,79 | 278 ± 29 ± 47 → 157 ± 94 ± 80 |
| D⁰ | 8–30 | 10,1 / 20,5 / 22,9 | **25,7 / 23,8 / 35,1** | 1,53 | 20,0–27,6 | 0,87–1,20 | +4,1 ± 21,4 | 0,77 | 100 ± 10 ± 21 → 96 ± 25 ± 23 |
| D⁺ (M0) | 3–5 | 7,8 / 13,7 / 15,9 | **46,0 / 154,7 / 161** | 10,2 | 15,9–48,6 | 1,00–3,06 | −15,5 ± 43,9 | 0,68 | 536 ± 42 ± 73 → 522 ± 240 ± 807 |
| D⁺ (M0) | 5–8 | 4,6 / 16,5 / 17,3 | **24,9 / 44,5 / 51,0** | 2,95 | 17,7–30,6 | 1,02–1,77 | +12,9 ± 24,7 | 0,73 | 167 ± 8 ± 28 → 151 ± 38 ± 67 |
| D⁺ (M0) | 8–30 | 8,7 / 14,3 / 16,9 | **25,4 / 30,6 / 39,9** | 2,36 | 17,0–26,7 | 1,01–1,58 | +0,5 ± 21,3 | 0,68 | 52,0 ± 4,5 ± 7,4 → 38,8 ± 9,9 ± 11,9 |

(F70: zakres = obciążenie fitu wzięte z WP70 [dolna granica, korzystna] albo z F0 [górna]. D⁺ wariant A2: g(F0) = 8,4 / 2,75 / 1,95; g(F70) = 1,00–2,88 / 1,00–1,51 / 0,86–1,38 — ten sam werdykt.)

Przekroje z F0 są zgodne z `final.json` w granicach własnych (dużych) błędów: F0/final = 1,37 / 0,56 / 0,96 (D⁰), 0,97 / 0,90 / 0,75 (D⁺). Nie niosą nowej informacji o napięciach z XSEC_final.

**Dlaczego:** kategoria L = [100, 70) zawiera 30 % sygnału i 94–99,3 % wszystkich kandydatów (D⁺ 3–5 GeV: N_L = 1922 ± 2094 na 1,32·10⁷ kandydatów; D⁰ 3–5: 4424 ± 1918 na 4,8·10⁶). Uzysk w L ma błąd ≥ 40 % (w trzech komórkach > 100 %) i to on jest błędem N_tot;
kategorie czyste (A, B, C) „dominują statystycznie” tylko wtedy, gdy ułamki f_k są znane — a wtedy błąd SF wraca tylnymi drzwiami. Usunięcie wiersza SF (3–12 %) kosztuje zastąpienie go błędem uzysku w L (20–60 % stat., 20–150 % model tła).

## 1. Co zrobiono

- Kategorie ustalone z góry: **L = [100, 70), C = [70, 50), B = [50, 30), A = [30, 0)** % wydajności sygnału MC (fd_nd), per gatunek i bin. Fit łączny: wspólne dmu i skala szerokości, szablon sygnału per kategoria (fd_nd), swobodne N_k i tło per kategoria.
  D⁰: model M8 jak w XSEC_final (szablon DMM + CM i κ per kategoria wprost z `bkgcat.npz`, r = 1). D⁺: M0 (jak `final.json`) i obok A2 (wybór `FITBIAS_Dplus/KRYTERIA.md` §8; raportu tamtego studium w chwili startu nie było).
- Kod: `skrypty/cf_fit.py` (kopia `fb_fit.fit` + więzy Gaussa na ułamki skumulowane + tryb ułamków stałych), `cf_lib.py`, `cf_run.py` (mc → toys → data, w tej kolejności), `cf_report.py` (własna kopia rachunku przekroju).
- Mechanika: fit jednej kategorii WP70 odtwarza `final.json` co do cyfry (4059 ± 357 / 1536 ± 124 / 700 ± 54; 3547 ± 257 / 1651 ± 75 / 687 ± 40); liczniki prawdy WP70 w mb_nd = `N_reco_truth` (asercja).
  Tożsamość WPTIGHT: K = 3 (A + B + C) = ΣCBA z K = 4 = N(WP70) w granicach błędu (T7: 3981 ± 363 wobec 4059; 664 ± 53 wobec 700; D⁺ 3263 ± 194 wobec 3547, 1671 ± 74 wobec 1651, 662 ± 38 wobec 687).
  Wyjątek D⁰ 5–8 GeV: 1320 ± 115 wobec 1536 — fit łączny ląduje w znanym drugim minimum skali szerokości (sc = 0,94; XSEC_final §1), nie jest to informacja z kategorii.
- Metryka i progi wg `KRYTERIA.md` §4–5: σ_tot = stat ⊕ syst ⊕ L; reguły stat. / model / obciążenie identyczne jak w `final.json` per gatunek; wiersze wspólne przepisane; C_tot = 32,8 / 43,0 / 55,6 % (D⁰), 23,0 / 32,8 / 45,0 % (D⁺).

## 2. MC mb_nd z pełnym tłem (T3)

| komórka | b(N_tot) ± σ_b [%] | b(N₇₀) ± σ_b [%] | f_L fit / prawda | E₇₀ fit / prawda |
|---|---|---|---|---|
| D⁰ 3–5 | +12,7 ± 41,1 | −22,0 ± 12,7 | 0,54 ± 0,19 / 0,32 | 0,46 / 0,68 |
| D⁰ 5–8 | +14,9 ± 25,3 | −4,6 ± 8,5 | 0,39 ± 0,15 / 0,30 | 0,61 / 0,70 |
| D⁰ 8–30 | +4,1 ± 21,4 | +7,2 ± 10,1 | 0,28 ± 0,11 / 0,29 | 0,72 / 0,71 |
| D⁺ 3–5 (M0) | −15,5 ± 43,9 | +4,8 ± 8,0 | 0,21 ± 0,93 / 0,30 | 0,79 / 0,70 |
| D⁺ 5–8 (M0) | +12,9 ± 24,7 | +11,4 ± 6,4 | 0,34 ± 0,12 / 0,30 | 0,66 / 0,70 |
| D⁺ 8–30 (M0) | +0,5 ± 21,3 | −0,8 ± 7,2 | 0,32 ± 0,11 / 0,32 | 0,68 / 0,69 |

- Kryterium obciążenia (|b| ≤ 5 % lub ≤ 1σ_b) jest formalnie spełnione wszędzie — ale tylko dlatego, że σ_b = 19–44 %: domknięcie N_tot na jednej realizacji mb_nd **nie ma mocy** (błąd 2–5 × większy niż dla WP70). Do budżetu wchodzi |b| ⊕ σ_b = 22–47 % (WP70: 6–24 %).
- f_k na MC odtwarzane z błędem ±0,04–0,05 (C, B, A) i ±0,10–0,19 (L) przy ≥ 5 GeV; przy 3–5 GeV dla D⁰ fit przenosi sygnał do L (0,54 wobec 0,32; 1,2σ), dla D⁺ błąd f_L ≈ 0,9. p(χ²) ≥ 0,075 wszędzie.

## 3. Zabawki przy czystości danych (T4; 150 / komórkę, tło z pasm danych, okno zaślepione)

- F0: obciążenie −6,9 … +2,7 %, rozrzut **19–41 %** (WP70: 5–13 %), ważne ≥ 0,98; pokrycie 0,66–0,79 (D⁰ 0,77–0,79 — 1,3 × HESSE lekko zawyża; D⁺ 0,66–0,73) — warunek 0,60–0,76 spełniają tylko D⁺.
- Błąd ścieżki (fit ⊕ SF) w zabawkach, cięcie ⊕ SF → F0 → F70: D⁰ 12,5 → 35,7 → 16,7 / 17,7 → 27,9 → 15,5 / 15,9 → 23,4 → 14,5 %; D⁺ (M0) 9,1 → 40,9 → 14,2 / 8,2 → 20,4 → 6,8 / 9,1 → 20,9 → 8,5 %.
  F70 jest gorsze od cięcia przy 3–5 GeV (+34 … +56 %), lepsze o 7–17 % przy ≥ 5 GeV — w samej części „fit ⊕ SF”, przy tle z rodziny modelu (zabawki nie testują tła pod pikiem).
- Opisowo, prawda tła rzędu + 1: obciążenie F0 do +9 % (D⁰ 5–8), WP70 do +6,6 %; w granicach błędów zabawek.

## 4. Dane 427929 (T1, T2, T5–T7)

- **N_tot (F0):** 8412 ± 2557 / 1235 ± 738 / 980 ± 252 (D⁰); 5188 ± 2388 / 2140 ± 533 / 755 ± 192 (D⁺ M0). Fity ważne, bez granic; p(χ²) łączne 0,02–0,65, najniższe per kategoria 0,013 (D⁺ 3–5, kat. B) — bramka 10⁻³ spełniona.
- **Model tła dominuje tam, gdzie tło jest największe:** D⁺ 3–5 GeV rząd − 1 daje N_tot = 12 827 wobec 5188 (+147 %); D⁺ 5–8: 1424 wobec 2140 (−33 %); D⁰ 3–5: 5852–6069 wobec 8412 (−30 %); D⁰ 5–8: 1679–1740 wobec 1235 (+41 %).
  Opisowo rząd + 2 tylko w L: D⁰ 3–5 8412 → 4373 ± 2591; K = 5 zamiast K = 4: D⁰ 5–8 1235 → 2050 ± 482. Uzysk bez cięcia nie jest określony lepiej niż ±30–150 %.
- **f_k dane vs MC (= „SF_k” z fitu; T5):** przy ≥ 5 GeV zgodne z 1 w granicach ±0,2–0,5 na kategorię (np. D⁰ 8–30: 1,06 ± 0,47 / 1,02 ± 0,27 / 1,11 ± 0,27 / 0,85 ± 0,18 dla L / C / B / A; D⁺ 5–8 M0: 0,72 ± 0,75 / 1,29 ± 0,38 / 1,03 ± 0,34 / 1,07 ± 0,31).
  D⁰ 3–5: SF_L = 1,72 ± 0,39, SF_A = 0,64 ± 0,16 — **to samo przesunięcie co fit na mb_nd względem prawdy (1,69 ± 0,58; 0,87 ± 0,37)**, czyli obciążenie fitu, nie różnica dane/MC. D⁰ 5–8: N_L = −152 ± 555, ułamki nieokreślone.
- **SF₇₀ z fitu (E₇₀^dane/E₇₀^MC) wobec niezależnych:** D⁰ 0,68 ± 0,17 (MC fit/prawda 0,68) / nieokreślone / 0,97 ± 0,20 wobec `final.json` 0,967 ± 0,033 / 1,005 ± 0,126 / 0,965 ± 0,097 i tagu D* 427929 (WP70) 1,31 ± 0,43 / 1,00 ± 0,13 / 0,96 ± 0,10;
  D⁺ (M0) 0,91 ± 0,74 / 1,12 ± 0,33 / 1,25 ± 0,50 wobec SFABS 0,95 / 0,99 / 0,99 ± 0,05. Zgodne, ale 2–15 × mniej precyzyjne niż SF już używane — fit K = 4 nie jest konkurencyjnym pomiarem SF (potwierdza SFABS „mianownik zawsze za brudny”).
- **Przekrój (T7):** tabela w §0. F70: 1004 / 236 / 97 µb (D⁰), 499 / 169 / 48 µb (D⁺ M0) — F70/final = 0,98 / 0,85 / 0,97 i 0,93 / 1,01 / 0,92; odchyłka D⁰ 5–8 (0,85) to drugie minimum sc, jak wyżej.

## 5. Pytanie (ii): ile zysku z danych, ile z więzów (T6; σ(N)/N z samego fitu, dane)

| komórka | σ_SF | cięcie WP70 ⊕ SF | F0 swobodne | **F70** (jeden więz = jedyny SF, jaki istnieje) | F3 (trzy niezależne więzy — założenie toya) | Ffix (ułamki z MC, SF bez błędu) |
|---|---|---|---|---|---|---|
| D⁰ 3–5 | 3,4 | 11,9 | 30,4 | 11,2 (−6 %) | 7,2 (−39 %) | 6,8 |
| D⁰ 5–8 | 12,5 | 16,3 | 59,8 | 14,6 (−10 %) | 11,0 (−33 %) | 8,9 |
| D⁰ 8–30 | 10,1 | 14,2 | 25,7 | 13,2 (−7 %) | 11,2 (−21 %) | 9,9 |
| D⁺ 3–5 | 5,0 | 8,8 | 40,6 | 8,0 (−9 %) | 5,2 (−41 %) | 4,2 |
| D⁺ 5–8 | 5,2 | 6,9 | 24,7 | 6,6 (−4 %) | 5,0 (−28 %) | 4,0 |
| D⁺ 8–30 | 5,2 | 7,8 | 23,5 | 7,4 (−5 %) | 6,2 (−21 %) | 5,4 |

- **Z danych (F70 wobec cięcia ⊕ SF): −4 … −10 % części „fit ⊕ SF”** — informacja z kategorii L i ze wspólnego kształtu. W budżecie całkowitym (obciążenie fitu 6–24 %, śledzenie 5–8 %) to g = 0,97–1,02 (dolna granica), czyli nic.
- **Z więzów (F70 → F3 → Ffix): −21 … −41 %** — dokładnie zakres toya (−23 … −54 %). Pochodzi wyłącznie z założenia, że wydajności przy WP50 i WP30 są znane niezależnie z tym samym błędem co przy WP70. Takich pomiarów nie ma
  (D⁺: brak tagu; D⁰ tag D* w 427929 przy WP50: 0,95 ± 0,13 / 0,97 ± 0,12, przy 3–5 GeV ± 0,43); zmierzone SF na sąsiednich WP są ponadto silnie skorelowane, nie niezależne. F3 z SF wziętym z tego samego piku byłby tożsamością WPTIGHT.

## 6. Skąd różnica względem toya z raportu (`FIT_scorebins`, 17.11; „−23 … −54 % błędu całkowitego”)

1. „Błąd całkowity” toya = stat ⊕ SF. W `final.json` SF to 3–12 % z 16–29 %; resztę (obciążenie fitu, model, śledzenie) toy pomija, a te wiersze się nie zmniejszają (dla F0 rosną).
2. Zysk toya dotyczy ułamków **z więzami Gaussa na każdej granicy kategorii, niezależnymi, σ = 5–10 %** (tu F3: odtworzony, −21 … −41 %). Sam toy podawał dla ułamków swobodnych σ = 8–84 % — to jest wynik F0 (21–60 %).
3. Statystyka i czystość: toy na 30 % runu 428770 z N_S z S/B MC (6–9 tys. D⁰ na bin; 5–8 GeV zawyżone × 1,7), kształt sygnału stały, tło z rodziny modelu. W 427929 N₇₀ = 0,7–4 tys., a szerokość jest wolna / ograniczona.
4. Model tła w kategorii brudnej: w toyu exp(pol2) dopasowany do pasm i nim generowany; na danych zmiana rzędu o 1 przesuwa N_tot o 7–147 %.

## 7. B2 — pojedynczy punkt pracy wg metryki z fitu piku (T8)

σ_reszta(WP) = stat ⊕ model ⊕ obciążenie(mb_nd) ⊕ wspólne ⊕ L (bez SF); σ_SF^max = największy błąd niezależnego SF, przy którym błąd całkowity spada o ≥ 15 % względem `final.json`.

| komórka | σ_reszta [%]: WP70 / 60 / 50 / 40 / 30 | σ_SF^max [%]: WP50 / 40 / 30 | niezależny SF przy tym WP | werdykt |
|---|---|---|---|---|
| D⁰ 3–5 | 30,1 / 28,3 / **16,5 / 15,0 / 16,3** | 18,5 / 19,7 / 18,7 | brak (tag D* 427929 ± 0,43; SFABS dotyczy WP70) | **wskazanie warunkowe** (niżej) |
| D⁰ 5–8 | 16,5 / 15,3 / 15,5 / 15,5 / 17,1 | 6,9 / 6,9 / — | tag D* WP50: ± 13,6 % → razem 20,6 % (g = 1,03) | brak zysku |
| D⁰ 8–30 | 21,6 / 23,7 / 20,6 / 24,3 / 20,5 | niemożliwe | tag D* WP50: ± 12,3 % → 24,0 % (g = 1,05) | brak zysku |
| D⁺ 3–5 (M0) | 15,9 / 18,9 / 18,6 / 16,2 / 18,0 | niemożliwe | brak | brak zysku |
| D⁺ 5–8 (M0) | 17,0 / 14,8 / **12,7 / 12,9 / 13,0** | 7,4 / 7,2 / 6,8 | brak | wskazanie warunkowe |
| D⁺ 8–30 (M0) | 16,7 / 17,9 / 16,6 / 16,3 / 16,7 | niemożliwe | brak | brak zysku |

- **Wg kryteriów: 0 komórek „ZYSK”** (nigdzie nie istnieje niezależny SF dla WP ≠ 70 o wystarczającej precyzji).
- D⁰ 3–5 GeV: σ_reszta spada z 30 % do 15–17 % dla WP ≤ 50 w trzech sąsiednich punktach. **Cały spadek pochodzi z wiersza „obciążenie fitu”**: domknięcie mb_nd −22,0 ± 12,7 % (WP70), −21,7 ± 9,9 % (WP60) wobec −5,4 ± 8,6 / +4,2 ± 8,0 / −0,6 ± 8,1 % (WP50 / 40 / 30); część statystyczna maleje tylko 11,5 → 7–9 %.
  BKGCAT oceniał dołek przy WP70 jako fluktuację 2,4σ jednej realizacji mb_nd zmieniającą znak ze score — wybór WP tam, gdzie pojedyncza realizacja akurat domyka, jest wyborem na szumie. Zysk byłby realny dopiero z (a) niezależnym MC potwierdzającym brak obciążenia przy WP ≤ 50 i (b) SF(WP50) przeniesionym z tagu D* ≥ 5 GeV z błędem ≤ 18 % (dla WP70 przeniesienie dało ± 3 %; dla WP50 niezrobione).
  N/ε_MC względem WP70: 0,99 / 0,92 / 0,91 (WP50 / 40 / 30) — dryf −8 % przy ciasnych WP, w granicach błędu.
- D⁺ 5–8 GeV (M0): 17,0 → 12,7 % przy WP50, znów przez domknięcie (+11,4 ± 6,4 → +3,0 ± 6,0 %); potrzebny niezależny SF(WP50) z błędem ≤ 7 % — dla D⁺ nie istnieje (jedyna droga: przeniesienie z D⁰, przy WP70 ± 5 %). W modelu A2 ten sam obraz dopiero przy WP ≤ 40.
- D⁺ 3–5: N/ε spada do 0,79–0,85 przy ciasnych WP przy domknięciu MC −10 … −12 % — znany „dryf” (WPTIGHT / FITBIAS), zaostrzanie szkodzi.

## 8. Zastrzeżenia

- Post hoc; 6 komórek × 2 warianty ułamków + 24 punkty B2, bez poprawki na wielokrotność. Wynik negatywny B1 jest odporny (g ≥ 1,5 dla F0; dla F70 dolna granica g ≥ 0,97 poza dwiema komórkami na 0,86–0,87 — D⁰ 8–30 i D⁺ 8–30 A2 — gdzie różnicę robi mniejszy rozrzut wariantów tła w ficie łącznym, a górna granica to 1,2–1,4; D⁺ 8–30 A2 nie spełnia pokrycia: 0,53).
- Kategorii nie optymalizowano (K = 4 ustalone z góry); inny podział L (K = 5) nie pomaga (D⁰ 3–5: 7979 ± 1863).
- Domknięcie na mb_nd to jedna realizacja o statystyce ≈ danych; szablony pikujące są „w próbce”.
- Dla F70 obciążenie fitu podano jako zakres (WP70 … F0), nie zmierzono go osobno. Część „ułamki” błędu F70 zawiera też swobodę podziału C / B / A.
- Reguła stat. D⁺ (max(HESSE, bootstrap)) zaniża błąd także dla WP70 (pokrycie 0,57–0,63) — dotyczy obu stron porównania jednakowo.
- Odstępstwa od `KRYTERIA.md`: (i) po etapie MC poprawiono pamięć podręczną szablonów w `cf_fit.py` (przyspieszenie; wartości funkcji bez zmian, MC nie powtarzano); (ii) pierwsze uruchomienie zabawek przerwano po ~11 min bez wyników (za wolne) i uruchomiono ponownie z tymi samymi ziarnami.
- Model D⁺ „A2” wzięty z `FITBIAS_Dplus/KRYTERIA.md` §8, nie z raportu (nie istniał w chwili startu) — wariant, nie nominał.

## 9. Pliki

`KRYTERIA.md`, `RAPORT.md`, `NOTA_tekst_EN.md`, `skrypty/{cf_fit,cf_lib,cf_run,cf_report}.py`, `tabele/catfit.md`, `out/{catfit.json, mc.json, toys.json, data.json, *.log}`.
Baza wiedzy: `charm-bkg-notes/podejscia/Fit jednoczesny w kategoriach score na danych 427929 (CATFIT) - brak zysku, uzysk bez cięcia ma błąd 25-60 procent.md` + wiersz w `00_INDEX.md` + wpis w dzienniku.
Pozostałości poza katalogiem: brak (żadnych plików w `/tmp/klyczak_catfit`, nic na lxplus / EOS).
