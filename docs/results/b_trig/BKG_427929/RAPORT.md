# BKG_427929 — ile niezgodności tła dane/MC to tryger, a ile zostaje w czystym minbias (17.11.2026)

Projekt B (tło podobne do sygnału), D⁰→Kπ, klasyfikator H. Praca lokalna na `B_trig/runs/lxplus/bt_*.npz` (tylko odczyt).

## 0. Odpowiedź

- **Tryger odpowiada za około 1/3 niezgodności ε_B dane/MC przy 3–5 GeV i około 1/2 przy 5–8 GeV. Reszta zostaje w czystym minbias.**
  W runie 427929 (MBTS) R = ε_B(dane)/ε_B(MC) = 0,86–0,88 przy 3–5 GeV (WP95–80), wobec 0,79–0,82 w 428770. Odchylenie od 1 w 427929 ma ponad 8σ (WP95, 3–5 GeV, stat ⊕ syst).
- Przy 8–30 GeV pomiar pod pikiem nie rozstrzyga (błąd R w 427929: 0,08–0,38). Zliczenie w samych pasmach bocznych pokazuje różnicę runów także tam (3–5σ).
- **Część trygerowa daje się usunąć przeważeniem MC w zmiennych aktywności (bez score).** Po takim przeważeniu 428770 ląduje na poziomie surowego 427929: 0,86–0,90 (3–5 GeV), 0,89–0,94 (5–8 GeV).
- **Reszta w 427929 nie siedzi w aktywności zdarzenia, izolacji ani pT.** To samo przeważenie nie zmienia tam R (0,85–0,88 przy 3–5 GeV). Domyka dopiero przeważenie zawierające score: R = 0,98–1,02 przy 3–5 GeV (WP95–70), χ² = 13,2/15.
- Rekomendacja: systematykę modelu tła mierzyć na 427929, nie na 428770. Liczby w §5.

## 1. Dane, definicje, czego brakuje

| próbka | wybór | kandydaci D⁰ (1,75–2,00 GeV) | w pasmach bocznych | pod pikiem |
|---|---|---|---|---|
| dane 427929 | GRL i `HLT_mb_mbts_L1MBTS_1` (klasa MB427929) | 1 197 340 | 472 974 | 532 145 |
| dane 428770 | wszystko (82 % L1_J12) | 4 414 546 | 1 751 934 | 1 953 452 |
| MC mb_nd | inkluzywne, wagi 1 | 932 616 | 368 238 | 414 406 |
| MC mb_nd, truth ≥ 8 | definicja btrig/charmbkg | 920 112 | 364 310 | 407 607 |

- Źródło: tablice `all_*` (inkluzywne D⁰, 30 % zdarzeń wg hasza). Z 427929 odpada 22 % kandydatów bez GRL (część z μ do 32).
- Regiony jak w `val_epsb.default_regions('D0')`: SR |m−m₀| < 55 MeV, SB |m−m₀| ≥ 75 MeV w oknie 1,75–2,00 (lewe pasmo ma tylko 40 MeV).
- Progi WP: kwantyle score prawdziwych D⁰ w fd_nd, per bin pT — odtwarzają tabelę z `RAPORT_DANE.md`.
- **MC tła w kodzie** (`btrig_analyze.cmd_epsb`): `all_truth ≥ 8` w SR, zliczenie wprost, wagi 1 dla klas `all`, `all427929`, `MB427929` (emulacja tylko dla klas trygera 428770). Dwie rzeczy, których ta definicja nie mówi wprost:
  - truth 100/101 to odbicia K↔π prawdziwego D⁰ (`mark_swaps`: truth += 100) i wchodzą do „≥ 8”; średni score 0,91.
  - truth 2–7 (częściowo zrekonstruowany powab, gładki w masie) jest wykluczony, a w danych interpolacja z SB go zawiera. W SB dodaje 13 % (3–5), 47 % (5–8) i 54 % (8–30 GeV) kandydatów przechodzących WP50 ponad truth ≥ 8.
  Dlatego podaję dwa mianowniki: truth ≥ 8 (ciągłość z btrig) i truth ≥ 2 (całe tło niepikujące; spójny z tym, co mierzy interpolacja).
- **Czego brakuje:** w `all_*` są tylko score, pT, masa, ptcone40/100, nSelBL, lead pT, ntrk500. Nie ma η ani 18 cech H. Pełne cechy są wyłącznie w parach D⁰←D*, gdzie masa sięga tylko ±75 MeV — pasm bocznych tam nie ma. Kroki 1 i 3 są więc na 6 zmiennych + score; 18 cech tylko w C2ST na paśmie 55–75 MeV (§3).

## 2. ε_B pod pikiem, dane/MC (`s1_epsb.py`, `out/s1.md`)

Estymator: ε_B(SR) = N_pass/N_all z dopasowań ML exp(pol2) do SB, kształt wolny (jak `eps_b_pol2_ml`). Odtwarza `btrig_epsb.md` co do 3 cyfr.
Błędy: bootstrap Poissona histogramów (kategoria score × masa), 400 replik. Systematyka interpolacji: największe odchylenie wariantów {sym, liniowa interpolacja lewe/prawe, exp(pol1)} od pol2.

R = ε_B(dane)/ε_B(MC truth ≥ 8), ± stat danych ± stat MC ± syst:

| bin | WP | 427929 | 428770 | D = ε(427929)/ε(428770) | różnica runów |
|---|---|---|---|---|---|
| 3–5 | 0,95 | 0,883 ± 0,007 ± 0,002 ± 0,012 | 0,822 ± 0,004 ± 0,002 ± 0,014 | 1,074 ± 0,010 ± 0,004 | 7,1σ |
| 3–5 | 0,90 | 0,863 ± 0,010 ± 0,003 ± 0,017 | 0,800 ± 0,006 ± 0,003 ± 0,008 | 1,079 ± 0,015 ± 0,010 | 4,4σ |
| 3–5 | 0,80 | 0,867 ± 0,016 ± 0,005 ± 0,033 | 0,795 ± 0,009 ± 0,004 ± 0,017 | 1,090 ± 0,023 ± 0,018 | 3,1σ |
| 3–5 | 0,70 | 0,899 ± 0,026 ± 0,008 ± 0,034 | 0,816 ± 0,014 ± 0,007 ± 0,018 | 1,101 ± 0,036 ± 0,017 | 2,6σ |
| 3–5 | 0,50 | 1,015 ± 0,066 ± 0,020 ± 0,063 | 0,813 ± 0,031 ± 0,016 ± 0,023 | 1,248 ± 0,093 ± 0,062 | 2,2σ |
| 5–8 | 0,95 | 0,922 ± 0,025 ± 0,007 ± 0,049 | 0,796 ± 0,009 ± 0,006 ± 0,020 | 1,159 ± 0,034 ± 0,035 | 3,3σ |
| 5–8 | 0,90 | 0,879 ± 0,033 ± 0,009 ± 0,033 | 0,773 ± 0,013 ± 0,008 ± 0,036 | 1,137 ± 0,047 ± 0,042 | 2,2σ |
| 5–8 | 0,80 | 0,917 ± 0,058 ± 0,016 ± 0,073 | 0,722 ± 0,019 ± 0,013 ± 0,040 | 1,271 ± 0,087 ± 0,036 | 2,9σ |
| 5–8 | 0,70 | 0,909 ± 0,083 ± 0,024 ± 0,068 | 0,686 ± 0,026 ± 0,018 ± 0,027 | 1,325 ± 0,133 ± 0,048 | 2,3σ |
| 5–8 | 0,50 | 0,913 ± 0,170 ± 0,047 ± 0,079 | 0,707 ± 0,055 ± 0,036 ± 0,088 | 1,291 ± 0,271 ± 0,169 | 0,9σ |
| 8–30 | 0,95 | 0,818 ± 0,058 ± 0,017 ± 0,059 | 0,827 ± 0,013 ± 0,017 ± 0,034 | 0,989 ± 0,071 ± 0,103 | −0,1σ |
| 8–30 | 0,90 | 0,832 ± 0,079 ± 0,024 ± 0,055 | 0,785 ± 0,018 ± 0,023 ± 0,028 | 1,060 ± 0,106 ± 0,091 | 0,4σ |
| 8–30 | 0,80 | 0,927 ± 0,133 ± 0,038 ± 0,019 | 0,764 ± 0,027 ± 0,031 ± 0,070 | 1,213 ± 0,180 ± 0,125 | 1,0σ |
| 8–30 | 0,70 | 0,972 ± 0,186 ± 0,055 ± 0,064 | 0,724 ± 0,035 ± 0,041 ± 0,059 | 1,343 ± 0,272 ± 0,162 | 1,1σ |
| 8–30 | 0,50 | 0,980 ± 0,310 ± 0,096 ± 0,191 | 0,700 ± 0,064 ± 0,068 ± 0,149 | 1,399 ± 0,473 ± 0,616 | 0,5σ |

D to stosunek dane/dane — MC się skraca. Kolumna „różnica runów” = (D−1)/(stat ⊕ syst).

χ² zgodności z 1 (15 punktów, pełna kowariancja między WP z replik; syst dodana diagonalnie — przybliżenie):

| | stat | stat ⊕ syst |
|---|---|---|
| 427929, MC truth ≥ 8 | 323/15 | 126/15 (p = 1,5·10⁻¹⁹) |
| 428770, MC truth ≥ 8 | 1995/15 | 713/15 |
| 427929, MC truth ≥ 2 | 341/15 | 133/15 |
| 428770, MC truth ≥ 2 | 2089/15 | 833/15 |
| D = 427929/428770 względem 1 | 94,6/15 | 72,2/15 (p = 1,8·10⁻⁹) |
| D = 427929/(428770, tylko łańcuchy MBTS) | 27,6/15 | 18,4/15 (p = 0,24) |

- Żaden run nie jest zgodny z MC; 427929 jest bliżej, ale χ² zdominowane przez bin 3–5 GeV nadal odrzuca R = 1.
- Z mianownikiem truth ≥ 2 R w 427929 nie rośnie do 1 przy ciasnych WP: 3–5 GeV 0,88 / 0,86 / 0,86 / 0,88 / 0,94; 5–8 GeV 0,92 / 0,86 / 0,87 / 0,82 / 0,73. „0,87–1,02” z pierwszych liczb to częściowo skutek braku truth 2–7 w mianowniku.
- Stosunek dane/dane liczony w samych SB (zliczenie, bez ekstrapolacji; `out/s4.md`) jest dużo dokładniejszy i rośnie z ciasnotą WP: 3–5 GeV 1,077 ± 0,003 → 1,201 ± 0,029; 5–8 GeV 1,126 ± 0,010 → 1,393 ± 0,078; 8–30 GeV 1,080 ± 0,025 → 1,92 ± 0,21. To miara efektu trygera w SB, nie pod pikiem.
- Podzbiór MBTS w 428770 (głównie j20 i sptrk_pt4) nadal różni się od 427929 o 5–9 % przy 3–5 GeV (2–3,5σ na punkt): to nie jest próbka minbias.
- Domknięcie estymatora na MC (interpolacja z SB truth ≥ 8 / prawda w SR): pol2 0,965–0,986 ± 0,01–0,03 przy 3–5 GeV (WP95–70), sym 0,98–1,03. Przy 8–30 GeV pol2 na MC ma błąd 7–60 % — dlatego mianownik MC to zliczenie w SR, a nie ta sama interpolacja.

Udział trygera f = (R₄₂₇ − R₄₂₈)/(1 − R₄₂₈), mianownik truth ≥ 2 (z truth ≥ 8 podobnie do WP70):

| bin | WP95 | WP90 | WP80 | WP70 | WP50 |
|---|---|---|---|---|---|
| 3–5 | 0,34 ± 0,05 | 0,31 ± 0,06 | 0,34 ± 0,09 | 0,40 ± 0,15 | 0,74 ± 0,29 |
| 5–8 | 0,60 ± 0,13 | 0,43 ± 0,15 | 0,59 ± 0,19 | 0,53 ± 0,21 | 0,37 ± 0,33 |
| 8–30 | −0,04 ± 0,28 | 0,17 ± 0,30 | 0,46 ± 0,38 | 0,54 ± 0,42 | 0,50 ± 0,58 |

Błędy f są statystyczne (bez systematyki interpolacji).

## 3. Kształty i C2ST (`s2_shapes.py`, `out/s2.json`)

Pasma boczne, rozkłady unormowane, MC inkluzywne. Odległość KS (wartości p są zerowe przy tych licznościach — KS służy jako miara wielkości efektu):

| zmienna (3–30 GeV) | 427929 vs MC | 428770 vs MC | 427929 vs 428770 |
|---|---|---|---|
| score H | 0,042 | 0,077 | 0,035 |
| pT | 0,012 | 0,168 | 0,179 |
| ptcone40 | 0,069 | 0,356 | 0,399 |
| ptcone100 | 0,038 | 0,440 | 0,446 |
| nSelBL | 0,052 | 0,300 | 0,264 |
| lead pT | 0,0015 | 0,533 | 0,533 |
| ntrk500 | 0,062 | 0,325 | 0,277 |
| masa | 0,005 | 0,012 | 0,008 |
| N kandydatów | 473 001 / 368 181 | 1 751 720 / 368 181 | 473 001 / 1 751 720 |

- Score per bin pT (KS): 427929 vs MC 0,042 / 0,047 / 0,062; 428770 vs MC 0,067 / 0,077 / 0,096; dane–dane 0,026 / 0,033 / 0,036. χ²/ndf dla score przy 3–5 GeV: 1563/19, 5559/19, 995/19.
- W 427929 lead pT zgadza się z MC (KS 0,0015, χ² 21/19), a krotność jest w danych wyższa (średnie nSelBL 39,5 vs 37,5). W 428770 lead pT jest dwa razy większe (9,6 vs 4,6 GeV).
- Stosunek gęstości score (fig. 2): w 427929 spada do ≈ 0,87 dla score > 0,4 i jest tam płaski (3–5 GeV); w 428770 ≈ 0,79 i dalej maleje.
- η (tylko pary D*, pasmo 55–75 MeV, Δm > 150 MeV): KS 0,012 (427929 vs MC), 0,018 (428770 vs MC), 0,010 (dane–dane) — małe.

C2ST: HGB, 5-fold po zdarzeniach, 100 tys. na klasę (mniej w wyższych binach); błąd z rozrzutu foldów:

| cechy (SB, 3–30 GeV) | 427929 vs MC | 428770 vs MC | 427929 vs 428770 |
|---|---|---|---|
| score | 0,525 ± 0,001 | 0,544 ± 0,002 | 0,522 ± 0,001 |
| score + ptcone40/100 | 0,563 ± 0,001 | 0,787 ± 0,001 | 0,800 ± 0,001 |
| nSelBL, lead pT, ntrk500 | 0,560 ± 0,002 | 0,865 ± 0,001 | 0,850 ± 0,001 |
| 6 zmiennych bez score | 0,592 ± 0,002 | 0,876 ± 0,001 | 0,870 ± 0,001 |
| 6 zmiennych + score | 0,596 ± 0,002 | 0,876 ± 0,001 | 0,870 ± 0,001 |
| pary D*: 18 cech H (N = 38 801 / 38 801 / 53 938 na klasę) | 0,563 ± 0,002 | 0,624 ± 0,002 | 0,611 ± 0,002 |
| pary D*: 18 cech + ptcone + zdarzenie | 0,608 ± 0,002 | 0,872 ± 0,001 | 0,862 ± 0,001 |

- Dane–dane (0,870) ≈ 428770 vs MC (0,876): prawie cała rozróżnialność 428770 od MC to tryger, i siedzi w aktywności zdarzenia oraz izolacji, nie w score (0,522).
- 427929 vs MC nie jest 0,50: 0,596 na 6 zmiennych, 0,563 na 18 cechach H. To jest niezgodność modelowania w czystym minbias.
- Per bin pT wyniki są takie same jakościowo (np. 3–5 GeV, 6 zmiennych + score: 0,595 / 0,873 / 0,864; 8–30 GeV: 0,563 / 0,741 / 0,748).
- Zastrzeżenie do wierszy „pary D*”: to kandydaci z miękkim pionem w paśmie bliskim 55–75 MeV, nie pełne tło pod pikiem.

## 4. Przeważenie MC, test domknięcia (`s3_reweight.py`, `out/s3.json`)

Klasyfikator dane-SB vs MC-SB (inkluzywne), 5-fold po zdarzeniach, w = p/(1−p), obcięcie 0,05–20. Trening tylko na SB; MC pod pikiem nie bierze udziału w treningu (wagi = średnia z 5 modeli). Test: ε_B(SR) z ważonego MC wobec interpolacji pol2 danych.

Warianty cech: V0 = pT + nSelBL + lead pT + ntrk500; V1 = V0 + ptcone40/100 (bez score, bez masy); V2 = V1 + score; V3 = V2 + masa.

R po przeważeniu, 427929, mianownik truth ≥ 2 (błąd: stat danych ⊕ stat ważonego MC):

| bin | WP | przed | V1 (bez score) | V2 (+ score) |
|---|---|---|---|---|
| 3–5 | 0,95 / 0,90 / 0,80 / 0,70 / 0,50 | 0,884 / 0,862 / 0,859 / 0,880 / 0,936 | 0,876 / 0,854 / 0,850 / 0,867 / 0,913 | 0,987 ± 0,008 / 0,984 ± 0,012 / 0,986 ± 0,019 / 1,015 ± 0,031 / 1,097 ± 0,075 |
| 5–8 | jw. | 0,916 / 0,863 / 0,872 / 0,823 / 0,725 | 0,927 / 0,892 / 0,921 / 0,886 / 0,776 | 1,041 ± 0,029 / 1,019 ± 0,039 / 1,064 ± 0,069 / 1,034 ± 0,100 / 0,938 ± 0,180 |
| 8–30 | jw. | 0,789 / 0,784 / 0,826 / 0,820 / 0,775 | 0,792 / 0,801 / 0,864 / 0,875 / 0,829 | 0,904 ± 0,067 / 0,935 ± 0,090 / 1,023 ± 0,152 / 1,037 ± 0,208 / 1,018 ± 0,333 |

| 427929 | AUC w SB | n_eff/N w SR | χ² vs 1 przed → po (truth ≥ 2) | (truth ≥ 8) |
|---|---|---|---|---|
| V0 | 0,566 | 0,885 | 338 → 303 | 319 → 305 |
| V1 | 0,595 | 0,847 | 320 → 349 | 325 → 354 |
| V2 | 0,599 | 0,838 | 327 → 13,2 (p = 0,59) | 325 → 19,6 (p = 0,19) |
| V3 | 0,599 | 0,838 | 329 → 13,0 (p = 0,60) | 334 → 18,9 (p = 0,22) |

χ² „przed” różni się między wierszami o kilka procent — to szum losowania replik, nie efekt wariantu.

- **V0 i V1 nie domykają w 427929** (χ² bez zmian, R przesuwa się o < 0,02 przy 3–5 GeV). Rozkład score w SB też się nie poprawia (KS 0,046 → 0,050).
- **V2 domyka**: wszystkie 15 punktów zgodne z 1. Przy 3–5 GeV zostaje −1,3…−1,6 % przy WP95–80 (1–1,6σ), tyle co obciążenie estymatora pol2 zmierzone na MC (0,965–0,986).
- Z mianownikiem truth ≥ 8 V2 przestrzeliwuje przy ciasnych WP (3–5 GeV WP50: 1,19 ± 0,08; 5–8 GeV WP80: 1,12 ± 0,07), bo wagi uczone są na inkluzywnym MC. Spójna para to trening inkluzywny + mianownik truth ≥ 2.
- V3 (z masą) = V2 w granicach 0,03.
- n_eff/N per bin w SR (V2): 0,84 / 0,85 / 0,85. Udział wag obciętych: 0.
- Wyciek masy: regresja m z cech na tle MC (truth ≥ 8, całe okno, OOF) daje R² = 0,006 (V1) i 0,013 (V2); korelacja wagi z masą w SB 0,004 (V1) i 0,020 (V2); średnia waga SR/SB = 0,984–0,986. Wycieku nie widać.

To samo dla 428770 (mianownik truth ≥ 2):

| bin | przed | V1 (bez score) | V2 (+ score) | n_eff/N |
|---|---|---|---|---|
| 3–5 | 0,823 / 0,799 / 0,789 / 0,799 / 0,750 | 0,880 ± 0,007 / 0,858 ± 0,010 / 0,874 ± 0,017 / 0,897 ± 0,026 / 0,841 ± 0,054 | 0,993 / 0,988 / 1,013 / 1,044 / 0,986 | 0,13 |
| 5–8 | 0,790 / 0,759 / 0,687 / 0,621 / 0,562 | 0,891 ± 0,019 / 0,914 ± 0,027 / 0,909 ± 0,041 / 0,898 ± 0,062 / 0,935 ± 0,125 | 1,001 / 1,035 / 1,033 / 1,016 / 1,063 | 0,23 |
| 8–30 | 0,798 / 0,739 / 0,681 / 0,611 / 0,553 | 0,943 ± 0,031 / 0,937 ± 0,043 / 1,009 ± 0,072 / 1,095 ± 0,109 / 1,301 ± 0,232 | 1,059 / 1,071 / 1,161 / 1,276 / 1,513 | 0,50 |

- V1 przenosi 428770 z 0,79–0,82 na 0,86–0,90 (3–5 GeV), czyli na poziom surowego 427929 (0,86–0,88). To niezależne potwierdzenie rozkładu „tryger + reszta”.
- Cena: AUC 0,876, n_eff/N = 0,13–0,23, 7–8 % wag na obcięciu. Przy 8–30 GeV V2 przestrzeliwuje (1,16–1,51).

## 5. Wniosek ilościowy

R = ε_B(dane)/ε_B(MC), mianownik truth ≥ 2, błąd stat ⊕ syst interpolacji (`out/s4.md`; tam też wersja truth ≥ 8):

| bin | WP | 428770 | 427929 | 427929 po przeważeniu V2 |
|---|---|---|---|---|
| 3–5 | 0,95 | 0,823 ± 0,014 | 0,884 ± 0,014 | 0,987 ± 0,015 |
| 3–5 | 0,90 | 0,799 ± 0,010 | 0,862 ± 0,019 | 0,984 ± 0,022 |
| 3–5 | 0,80 | 0,789 ± 0,019 | 0,860 ± 0,037 | 0,986 ± 0,042 |
| 3–5 | 0,70 | 0,799 ± 0,023 | 0,880 ± 0,043 | 1,015 ± 0,049 |
| 3–5 | 0,50 | 0,750 ± 0,038 | 0,936 ± 0,086 | 1,097 ± 0,101 |
| 5–8 | 0,95 | 0,790 ± 0,023 | 0,916 ± 0,055 | 1,041 ± 0,063 |
| 5–8 | 0,90 | 0,759 ± 0,038 | 0,863 ± 0,046 | 1,019 ± 0,055 |
| 5–8 | 0,80 | 0,687 ± 0,043 | 0,872 ± 0,090 | 1,064 ± 0,109 |
| 5–8 | 0,70 | 0,621 ± 0,037 | 0,823 ± 0,099 | 1,034 ± 0,126 |
| 5–8 | 0,50 | 0,562 ± 0,086 | 0,725 ± 0,153 | 0,938 ± 0,197 |
| 8–30 | 0,95 | 0,797 ± 0,039 | 0,789 ± 0,081 | 0,904 ± 0,093 |
| 8–30 | 0,90 | 0,740 ± 0,038 | 0,784 ± 0,093 | 0,935 ± 0,109 |
| 8–30 | 0,80 | 0,683 ± 0,072 | 0,828 ± 0,124 | 1,023 ± 0,154 |
| 8–30 | 0,70 | 0,612 ± 0,066 | 0,823 ± 0,172 | 1,037 ± 0,219 |
| 8–30 | 0,50 | 0,558 ± 0,137 | 0,780 ± 0,298 | 1,018 ± 0,388 |

Średnie po binach pT (wagi 1/σ², zdominowane przez 3–5 GeV):

| WP | 428770 | 427929 | 427929 po V2 |
|---|---|---|---|
| 0,95 | 0,812 ± 0,012 | 0,883 ± 0,013 | 0,988 ± 0,015 |
| 0,90 | 0,793 ± 0,009 | 0,859 ± 0,018 | 0,987 ± 0,020 |
| 0,80 | 0,767 ± 0,017 | 0,859 ± 0,033 | 0,998 ± 0,038 |
| 0,70 | 0,738 ± 0,019 | 0,869 ± 0,038 | 1,019 ± 0,045 |
| 0,50 | 0,709 ± 0,034 | 0,879 ± 0,073 | 1,061 ± 0,088 |

**Wyniki (zmierzone):**
1. Niezgodność w 428770: 19–29 % (średnio po pT). W 427929: 12–14 %. Różnica runów 7–10 % przy 3–5 GeV (do 7σ na punkt, χ² = 72/15 łącznie).
2. W czystym minbias niezgodność zostaje i jest istotna: 0,884 ± 0,014 przy WP95, 3–5 GeV.
3. Przeważenie w aktywności/izolacji/pT usuwa część trygerową (428770 → poziom 427929), ale nie rusza reszty w 427929.
4. Przeważenie ze score, uczone w SB, przenosi się pod pik w 427929 z dokładnością pomiaru: ±1,5 % (WP95), ±2 % (WP90), ±4 % (WP80), ±5 % (WP70), ±10 % (WP50) przy 3–5 GeV; ±6–20 % przy 5–8 GeV; ±9–39 % przy 8–30 GeV.

**Hipotezy (niesprawdzone tutaj):**
- Reszta 12–14 % w minbias pochodzi z modelowania cech kandydata (wierzchołek, parametry zderzenia), a nie ze składu tła. Poszlaka: w parach D* największe KS dla 427929 vs MC mają `pointing` (0,034), `impact_sig` i `abs_impact_z0` (0,024). Skład (udział truth 2–7, odbić) nie był tu zmieniany.
- Pełny SALAD na 18 cechach domknie 427929 tak samo jak V2. V2 używa score jako cechy, więc w SB domyka się z konstrukcji; test dotyczy tylko przeniesienia stosunku dane/MC z SB pod pik.
- Część trygerowa to zależność ε_B od aktywności/dżetu przy D⁰ (jak w `ISO_origin`); zmiennych dżetowych w npz nie ma.

**Rekomendacja:**
- Systematykę modelu tła (SALAD vs SB, domknięcie pod pikiem) mierzyć na 427929. Na 428770 ta sama liczba miesza modelowanie z trygerem, którego MC nie ma, a domknięcie wymaga wag o n_eff/N = 0,13–0,23.
- Obecny budżet SALAD ±5 / ±10 / ±25 % (WP80/70/50) jest zgodny z domknięciem w 427929 przy 3–5 GeV (±4 / ±5 / ±10 %) i 5–8 GeV w granicach błędów. Przy 8–30 GeV 427929 ma za mało statystyki, żeby go potwierdzić.
- MC tła dla D⁰ definiować jako truth ≥ 2 (`mc_truth_min: 2`), zgodnie z B1-D⁰; z truth ≥ 8 liczby przy WP70/50 i pT > 5 GeV są zawyżone.
- Dla 428770 stosować przeważenie w aktywności (lub per klasa trygera) jako osobny krok przed SALAD; sam w sobie daje ten sam poziom co minbias.

## 6. Zastrzeżenia

- Brak 18 cech i η w próbce inkluzywnej; SALAD v1.0 w konfiguracji produkcyjnej nie był uruchamiany.
- `all_*` to 30 % zdarzeń; 427929 po GRL to 78 % z tego. MC przy 8–30 GeV: 5,8 tys. kandydatów w SR (truth ≥ 2); w SB WP50 przechodzi 122.
- Błąd R po przeważeniu nie zawiera wariancji treningu wag (jedno ziarno, jeden podział).
- Systematyka interpolacji dodana diagonalnie w χ²; między WP jest skorelowana, więc p-wartości „stat ⊕ syst” są przybliżone.
- f zakłada addytywny podział (1 − R); przy 8–30 GeV jest nieokreślone.
- Lewe pasmo boczne ma 40 MeV, prawe 60 MeV — ekstrapolacja pol2 jest asymetryczna; warianty sym/liniowy/pol1 różnią się od pol2 o 1–8 % (do 19 % przy 8–30 GeV, WP50).
- Podzbiór MBTS w 428770 nie jest próbką minbias (j20, sptrk_pt4); zgodność z 427929 po dodaniu systematyki (p = 0,24) nie znaczy równoważności.

## 7. Pliki

| plik | rola |
|---|---|
| `bkg_common.py` | wczytanie, regiony, progi, estymatory ε_B (ML exp(pol-k), sym, liniowy) |
| `s1_epsb.py` → `out/s1.md`, `out/s1.json`, `out/s1_rep.npz` | ε_B pod pikiem, R, χ², różnica runów, domknięcie na MC, liczności |
| `s2_shapes.py` → `out/s2.json` | KS/χ² kształtów, stosunki score, C2ST (inkluzywne i pary D*) |
| `s3_reweight.py` → `out/s3.json`, `out/s3_w_d427.npz`, `out/s3_w_d428.npz` | przeważenie V0–V3, test domknięcia, wyciek masy, wagi |
| `s4_summary.py` → `out/s4.md`, `out/s4.json` | tabele zbiorcze i rysunki |
| `fig1_epsb_ratio.png` | R per WP i bin pT: 428770, 427929, 427929 po V2 (oba mianowniki) |
| `fig2_score_ratio_sb.png` | stosunek gęstości score w SB: dane/MC w obu runach i dane/dane |
| `fig3_c2st.png` | C2ST AUC per zestaw cech, trzy porównania |

Uruchomienie: `cd B_trig/BKG_427929 && P=~/Documents/charm/.venv_h5/bin/python; $P s1_epsb.py; $P s2_shapes.py; $P s3_reweight.py; $P s4_summary.py` (ok. 1 + 4 + 6 + 0,2 min, 4 wątki).
