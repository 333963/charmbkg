# FITBIAS — obciążenie fitu masy inkluzywnego piku przy małym S/B i próba naprawy (projekt B)

19.11.2026, ~22:35 CET. Kryteria: `KRYTERIA.md` (cz. 1 przed wszystkimi fitami; cz. 2 po diagnozie na MC, przed fitami danych; cz. 3 = poprawka BR przed fitami okna danych).
**Całość post hoc** względem etapów 4 / 5a, WPSCAN, WPTIGHT. Praca lokalna, bez jobów; bez AFP; pojedyncze D; drzewo H.
Wejścia: `WPSCAN_2026-11-19/runs/lxplus/wps_hist.npz` (W: dane + MC, prawda tylko sygnał/odbicie/reszta) i `MC_templates/runs/lxplus/mctmpl_{mb,fd}_nd.npz` (T: 24 kategorie prawdy, klasyfikator pokrewny, WP95/80/50).
Kontrola mechaniki: moduł `skrypty/fb_fit.py` odtwarza `wps_fit.json` (mb_nd, WP70) co do cyfry; M0 odtwarza WPTIGHT (SF_rel^MC D⁺ 0,917 / 0,966 / 0,959) i XSEC (b D⁰ −27 / −16 / −14 %).
Liczba: diagnoza 30 komórek (T) + 54 (W); wybór: 10 ustawień × 2 tryby (S, J) = 20 wariantów, 2 decyzje (D⁰, D⁺); zabawki 24 komórki × 150; dane 12 komórek. Bez poprawki na wielokrotność.

## Odpowiedź

1. **Nie ma jednego „zawyżenia”.** D⁺: sygnał zawyżony w kategorii luźnej (stąd SF_rel^MC 0,93–0,95). D⁰: uzysk bezwzględny **zaniżony** (−14…−27 % przy WP70) — to inny efekt i inny znak; stosunek pass/fail go nie widzi, bo jest wspólny dla kategorii.
2. **Przyczyny (liczby niżej):** D⁰ — brak w modelu składowej DecayModeMismatch (−5 % przy pT < 8 GeV) wzmacniany przez wolną szerokość (razem −10…−17 %), reszta = fluktuacja jednej realizacji tła mb_nd.
   D⁺ — CharmMismatch (głównie D_s) +5,5 ± 1,3 % w kategorii luźnej, PartComb + Comb +14 ± 6 % (niepikujące w sposób wykazany; 2σ), wolna szerokość bez wpływu na SF_rel. Odbicie K↔π i szablon sygnału fd vs mb — nie są przyczyną.
3. **Wybór wg kryteriów:** żadne z 10 ustawień nie spełnia D1 (pokrycie: rms pull 1,2–1,5 — HESSE zaniża błąd o ~25 %). Wg Q: **D⁰ → M8** (ograniczenie Gaussa + szablon DMM + CM + zakres 1,70–2,05), **D⁺ → M0 (nominał zostaje)**. Tryb łączny (d) nie lepszy od S.
4. **Kalibracja D⁰ na inkluzywnym piku się NIE odblokowuje:** σ(SF_rel^MC) = 0,127 / 0,057 / 0,053 (ε_t = 0,50; M0: 0,150 / 0,070 / 0,075) — bramka 0,05 niezaliczona w żadnym binie i żadnym ε. Błąd f ogranicza statystyka tła w kategorii luźnej, nie kształt sygnału.
5. **„Dryf” nie znika:** D⁺ bez zmian (model ten sam; SF_rel^MC 0,948 ± 0,023, dane/MC 1,00 ± 0,03); D⁰: SF_rel^MC = 1,053 ± 0,037 (prawda 0,99), dane 0,975 ± 0,039 / 0,976 ± 0,022, stosunek podwójny 0,93 ± 0,05 / 0,93 ± 0,04.

## 1. Diagnoza na MC z prawdą (tabele `tabele/diag_T.md`, `diag_T2.md`, `diag_W.md`, `diag_sig.md`, `decomp.md`)

Prawda w plikach: W — tylko sygnał / odbicie / reszta; pełny rozkład tła na kategorie odtworzony z T (inny, pokrewny klasyfikator; liczności zgodne z W na ~5 %).
b = N_fit/N_prawda − 1. Przyrosty przy dokładaniu składowych prawdy po kolei, model nominalny z kształtem stałym (reżim addytywny), średnie po binach lub per bin (3–5 / 5–8 / 8–30 GeV):

| przyczyna | D⁰ | D⁺ |
|---|---|---|
| szablon sygnału fd vs sygnał prawdy mb (bez tła) | zgodny: χ² p ≥ 0,05, sc 0,93–1,03, dmu ≤ 1 MeV poza 8–30 (do +3,5 MeV) | zgodny w 3–8 GeV; 8–30: sygnał mb szerszy (sc 1,04–1,10) → −7…−15 % w kategorii ciasnej |
| odbicie K↔π: f_refl fd vs mb | 0,744 vs 0,743 (WP70) — bez wpływu (Δb ≤ 0,001); Sig + Refl: b = −0,00…−0,01 (3–8 GeV) | — |
| **DecayModeMismatch** (prawdziwe D w innej modzie; D⁰: κ ≈ 0,9 N_S, ~20 % sygnału w oknie) | **−5,4 / −6,0 / +3,2 %** (WP80; to samo w WP50 i luźnej) przy kształcie stałym; **−9 / −16 / −6 %** przy wolnym | −1,2 ± 0,5 % |
| CharmMismatch | +1 % | **+2,4 ± 0,6 %** (WP80), **+5,5 ± 1,3 %** (luźna [80, 50)), +0,4 % (WP50); głównie CM_431 (D_s): frakcja sygnałopodobna 12 ± 5 % |
| PartComb | +1 ± 2 % (WP80) | +3,6 ± 1,9 % (WP80), **+9 ± 4,5 %** (luźna); 3–5: 0 ± 9, 5–8: +12 ± 7, 8–30: +15 ± 8 % |
| Comb (kombinatoryka) | −1 ± 5 % (WP80); D⁰ 3–5 WP70 (W): „reszta” −22 ± 10 % | +2,3 ± 1,5 % (WP80), +5 ± 4 % (luźna) |
| wolna szerokość/przesunięcie na tej samej realizacji (wolny − stały) | WP70 (W): −2 / −18 / −7 %; sc = 0,94 / 0,82 / 0,92 | WP70: 0 / +6,5 / +2 %; na SF_rel bez wpływu (0,917 → 0,909; 0,966 → 0,974; 0,959 → 0,961) |
| obciążenie własne estymatora z wolnym kształtem (zabawki z modelu, prawda = model) | +1,6 % (S/B 0,19) → +6,6 % (0,04) → +13 % (0,016); pull rms 1,2–1,5; rozrzut N 1,4–2,3 × większy niż przy kształcie stałym | +0,2…+1,4 % dla S/B ≥ 0,09; +6 % przy 0,036; pull rms 1,1–1,3 |
| rząd tła (nom + 1, + 2), zakres | nie usuwa (N_spur bez zmian w błędzie); exp(Czebyszew 2) na 250 MeV jest zdegenerowany z ogonami sygnału (fit samego sygnału: −4…−22 %) | nom + 1, + 2 bez zmian; zakres 1,75–2,10 przesuwa SF_rel o +2 pkt |

- **D⁺, pass/fail (to jest „zawyżenie” z WPTIGHT):** N_B (kategoria [70, 50)) = +30 ± 16 / +8 ± 11 / +10 ± 13 %, N_A = −5 ± 4 / +3 ± 4 / −7 ± 6 % (kształt stały). Z −5,3 ± 2,3 % (średnia SF_rel^MC − prawda, ε = 0,50):
  ~1 pkt składowa pikująca DMM + CM (M2 → M5), ~2 pkt opis tła w szerokim zakresie (M0 → M7), reszta (−2 ± 2,7 %; w 3–5 GeV −6 ± 5 %) nieodróżnialna od fluktuacji jednej realizacji mb_nd. Wolna szerokość: 0.
- **D⁰, uzysk WP70:** M0 −26,5 ± 11,6 / −16,2 ± 7,3 / −13,6 ± 8,8 % (identycznie jak XSEC). Z tego deterministycznie: DMM −5 % (stały) → −10…−17 % z wolną szerokością (5–8, 8–30 GeV: po ograniczeniu kształtu i dodaniu szablonu −1,5 ± 8 / −1 ± 10 %);
  3–5 GeV: −23 ± 10…14 % zostaje w każdym wariancie — „reszta” tła ma tam dołek −22 ± 10 % sygnału względem exp(pol2) (fluktuacja 2σ albo struktura PartComb; nierozstrzygalne tą statystyką).
- **BR w MC ≠ PDG (XSEC):** D⁰→Kπ ×1,364, D⁺→Kππ ×1,068 → w danych κ(DMM) większe o ten czynnik: efekt DMM dla D⁰ w danych ≈ −7 % zamiast −5 % (pT < 8 GeV), jeśli pozostałe mody w MC są poprawne (niesprawdzone).
- Frakcje sygnałopodobne PartComb z fd_nd nie przenoszą się na mb_nd (skład wymuszonych rozpadów: D⁺ 3–5 WP50 fd +7 ± 1,4 % vs mb −3 ± 2,5 %) — fd nie zwiększa statystyki tej diagnozy.

## 2. Warianty i wybór (tabela `tabele/select.md`; kryteria cz. 2)

| ustawienie | D⁰: b_A (Asimov, WP80) | b mb_nd WP70 (średnia) | zabawki: obciążenie / rms pull / σ_toy | D1 | Q | D⁺: b_A | b mb_nd WP70 | zabawki | D1 | Q |
|---|---|---|---|---|---|---|---|---|---|---|
| M0 nominał | 0,115 | −0,173 ± 0,051 | −0,039 / 1,46 / 0,143 | 1/12 | 0,184 | 0,015 | +0,057 ± 0,032 | +0,001 / 1,26 / 0,068 | 5/12 | 0,070 |
| M1 (a) Gauss | 0,035 | −0,082 ± 0,057 | −0,036 / 1,25 / 0,136 | 4/12 | 0,140 | 0,018 | +0,059 ± 0,033 | −0,002 / 1,21 / 0,067 | 8/12 | 0,069 |
| M2 (a′) stały | 0,034 | −0,071 ± 0,051 | −0,011 / 1,63 / 0,149 | 3/12 | 0,153 | 0,026 | +0,029 ± 0,028 | +0,006 / 1,79 / 0,082 | 3/12 | 0,086 |
| M4 (a + b) | 0,021 | −0,048 ± 0,057 | +0,003 / 1,26 / 0,147 | 6/12 | 0,148 | 0,016 | +0,041 ± 0,033 | −0,001 / 1,22 / 0,068 | 9/12 | 0,070 |
| M6 (c1) rząd + 1 | 0,112 | −0,198 ± 0,048 | −0,037 / 1,43 / 0,149 | 1/12 | 0,187 | 0,023 | +0,057 ± 0,032 | +0,003 / 1,28 / 0,071 | 3/12 | 0,074 |
| M7 (c2) zakres alt. | 0,116 | −0,150 ± 0,040 | −0,059 / 1,39 / 0,111 | 0/12 | 0,160 | 0,011 | +0,064 ± 0,037 | +0,012 / 1,31 / 0,089 | 2/12 | 0,090 |
| **M8 (a + b + c2)** | 0,049 | −0,049 ± 0,046 | −0,003 / 1,28 / 0,116 | 5/12 | **0,126** | 0,015 | +0,057 ± 0,037 | −0,002 / 1,27 / 0,083 | 4/12 | 0,084 |
(M3, M5, M9 w `select.md`.)
- **D1 niespełnione przez żadne ustawienie** (wymagało ≥ 10/12 komórek z |⟨pull⟩| ≤ 0,25 i rms pull ≤ 1,25): w najlepszych (M4, M9) obciążenie w zabawkach ≤ 0,3 %, ale HESSE zaniża błąd o 20–30 % (rms pull 1,2–1,3; pokrycie 68 % → 58–62 %).
  Kształt stały (M2, M5) ma rms pull do 4 — płaci za rzeczywistą niepewność szerokości. Przypadek „nic nie spełnia D1” nie był przewidziany w kryteriach — **odstępstwo: wybór wg samego Q**.
- **D⁰ → M8:** Q 0,126 (M0: 0,184); rozrzut N w zabawkach 0,116 vs 0,143; obciążenie w zabawkach −0,3 % vs −3,9 % (M0 nie ma składowej DMM, którą zabawki wstrzykują). D2 spełnione (3–5 GeV: −23 ± 10 %, 2,4σ — na granicy).
- **D⁺ → M0:** najmniejsze Q ma M1 (0,069), różnica do M0 < 5 % → reguła prostoty zostawia nominał. Dla D⁺ żaden wariant nie poprawia fitu jednej kategorii w mierzalny sposób.
- (d) fit łączny A + B: rms(E) ≤ 0,95 rms(N_70) w 3/12 (D⁰ M8) i 7/12 (D⁺ M0) komórek — kryterium niespełnione; zgodne z WPTIGHT (algebra).

## 3. Skutek dla danych (`tabele/data.md`)

**Uzyski WP70** (σ = max(HESSE, bootstrap po zdarzeniach); D⁰ M8 z κ × 1,364):

| gat. | run | pT | N stary (M0) | N nowy | ΔN/N | σ_N/N stary → nowy | p(χ²) nowy |
|---|---|---|---|---|---|---|---|
| D⁰ | 427929 | 3–5 | 4527 ± 970 | 4105 ± 463 | −9,3 % | 0,214 → 0,113 | 0,25 |
| D⁰ | 428770 | 3–5 | 6425 ± 869 | 6403 ± 639 | −0,3 % | 0,135 → 0,100 | 0,24 |
| D⁰ | 427929 | 5–8 | 1579 ± 178 | 1530 ± 126 | −3,1 % | 0,113 → 0,082 | 0,14 |
| D⁰ | 428770 | 5–8 | 4548 ± 380 | 4846 ± 302 | +6,5 % | 0,083 → 0,062 | 0,46 |
| D⁰ | 427929 | 8–30 | 764 ± 79 | 681 ± 55 | −10,9 % | 0,103 → 0,080 | 0,038 |
| D⁰ | 428770 | 8–30 | 7716 ± 353 | (6829 ± 261) | (−11,5 %) | (0,046 → 0,038) | **1·10⁻⁷ — model odrzucony, zostaje M0** |
| D⁺ | oba | wszystkie | bez zmian (M0) | — | 0 | bez zmian | — |
- Zmiany N mieszczą się w starym błędzie (≤ 1,1σ_stary); błąd D⁰ maleje ×0,53–0,84 (głównie ograniczenie kształtu + szerszy zakres). Uwaga: zabawki mówią, że rzeczywisty rozrzut jest ~1,3 × HESSE także dla M8.
- **428770, 8–30 GeV: dane odrzucają szablon DMM + CM** (p = 1·10⁻⁷ z κ × 1,364; 3·10⁻⁴ z κ_MC; M0: 0,43). Skala BR pogarsza opis; skład/kształt z mb_nd (≈ 900 zdarzeń szablonu, minimum-bias) nie przenosi się na próbkę dżetową przy wysokim pT.
  κ_MC zamiast κ × BR zmienia N o ≤ 1,2 % przy pT < 8 GeV, +3,4 % / +8,4 % przy 8–30.
- D⁺ opisowo (M4, niewybrany): N zmienia się o −6…+0 %, błąd bez zmian (`data.md` R1b).

**(i) Odblokowanie kalibracji D⁰ — NIE.** Bramka (σ ≤ 0,05 ∧ |Δ| ≤ 2,5σ) na mb_nd, M8: ε_t = 0,50: 1,209 ± 0,127 / 1,090 ± 0,057 / 0,994 ± 0,053; 0,35: 1,324 ± 0,168 / 1,061 ± 0,066 / 1,053 ± 0,076 → 0/9 punktów (M0: 0/9; M4: 0/9).
Poprawa σ o 15–30 % pochodzi z szerszego zakresu, nie z ograniczenia kształtu (M1, M4: σ bez zmian). 5–8 i 8–30 GeV są blisko progu (0,053–0,057), 3–5 GeV daleko (0,13; S/B w kategorii luźnej 0,02).
Oczekiwana σ(SF_rel) przy statystyce danych (zabawki, M8): 427929 0,083 / 0,065 / 0,068; 428770 0,050 / 0,027 / 0,022 — w runie dżetowym precyzja by wystarczyła, ale mb_nd (statystyka ≈ 427929, skład minimum-bias) nie może tego domknąć, a przy 8–30 GeV dane odrzucają model.
Metody R (sieć × drzewo) lokalnie nie powtórzyłem — kategorii łącznych nie ma w histogramach; pass/fail zagnieżdżony jest jej odpowiednikiem (te same σ w WPTIGHT). D⁺: bramka jak dotąd (ε = 0,50: 3/3).

**(ii) Dryf — nie znika.**
- D⁺ (model bez zmian): SF_rel^MC = 0,948 ± 0,023 / 0,931 ± 0,025 / 0,933 ± 0,028, dane 0,948 ± 0,023 (427929), 0,950 ± 0,015 (428770); podwójny stosunek 1,00 ± 0,03. Kryterium |SF_rel^MC − 1| ≤ 0,03 niespełnione.
  Opisowo M8 dla D⁺: SF_rel^MC = 0,942 / 0,995 / 1,002 (średnio −1,8 ± 2,7 %) — obciążenie maleje, kosztem σ_toy 0,083 vs 0,068; niewybrane.
- D⁰ (M8): SF_rel dane 0,975 ± 0,039 / 0,941 ± 0,045 / 0,945 ± 0,047 (427929), 0,976 ± 0,022 / 0,942 ± 0,026 / 0,937 ± 0,028 (428770); MC-fit 1,053 ± 0,037 / 1,073 ± 0,045 / 1,079 ± 0,048 przy prawdzie 0,99–1,00;
  podwójny stosunek 0,93 ± 0,05 / 0,88 ± 0,06 / 0,88 ± 0,06 (427929), 0,93 ± 0,04 / 0,88 ± 0,04 / 0,87 ± 0,05 (428770) — MC wspólne dla obu runów, więc to jedno odchylenie ~2–3σ, nie dwa; bez bramki bez wniosku.
  Dryf z fitów jednej kategorii d(0,50): 427929 −0,051 ± 0,062 → +0,005 ± 0,047; 428770 −0,111 ± 0,043 → +0,070 ± 0,026 (zmiana znaku pochodzi z 8–30 GeV, gdzie model jest odrzucony).

**XSEC:** `XSEC_2026-11-19/fitbias_override.json` zapisany w formacie haka `xsec.py` (klucze `D0|run|bin`: N, stat, fit_model_rel, fit_bias_rel) dla 5 komórek D⁰; **bez** `D0|428770|8-30` i bez D⁺. fit_bias_rel = |b| ⊕ σ_b (mb_nd, M8): 0,253 / 0,082 / 0,094
(M0: 0,29 / 0,18 / 0,16); fit_model_rel = max |ΔN/N| wobec M4, M9 i κ bez skali BR: 0,8–11 % (odstępstwo od cz. 3: zbiór „dopuszczalnych wg D1” jest pusty). `xsec.py` nie zmieniany i nie uruchamiany.

## Zastrzeżenia
- Post hoc; 20 wariantów oceniane na jednej realizacji mb_nd (statystyka ≈ jeden run 427929) — domknięcie b ma błąd 5–12 % i nie rozróżnia wariantów; ranking oparty na Asimovie i zabawkach.
- D1 niespełnione przez nic; wybór wg Q to rozszerzenie kryteriów po fakcie. Wybrany M8 nie ma poprawnego pokrycia (rms pull 1,28): błąd HESSE trzeba skalować ~1,3 albo brać z bootstrapu.
- „Asimov” nie jest całkiem bez szumu: sygnał, odbicie, DMM, CM to rzeczywiste zdarzenia mb_nd, a wygładzone tło (exp nom + 2) zachowuje część fluktuacji pod pikiem.
- Szablon (b): kształt i κ z tego samego mb_nd, na którym liczone jest domknięcie (częściowo w próbce); klasyfikator w T ≠ drzewo H z W, κ interpolowane w ε (brak WP70 w T); D⁰ 8–30: ~900 zdarzeń szablonu.
- Zabawki: tło z tej samej rodziny co model (nie testują tła pod pikiem); prawdziwy kształt losowany z priorytetu (a) — sprzyja wariantom Gaussa z założenia; S/B ustawione jak w MC (i 0,5×), nie z okna danych.
- Ograniczenie (a) dla D⁺ (0,5 MeV, 10 %) pochodzi z PEAK_pos (ten sam pik inkluzywny), nie z D*.
- Skala BR ×1,364 zakłada poprawne pozostałe mody w MC; dane przy 8–30 GeV jej nie popierają.
- Błędy d dla mb_nd w `data.md` R3 są przybliżone (wzór zagnieżdżony), nie do cytowania. 30 replik bootstrapu (błąd rms ±13 %).

## Czego nie dało się ustalić lokalnie
- Czy nadwyżka PartComb + Comb w kategorii luźnej D⁺ (+14 ± 6 %) i dołek tła D⁰ 3–5 GeV (−22 ± 10 %) są strukturą, czy fluktuacją — potrzeba mb_nd ≥ 4× większego albo kategorii prawdy PartComb w histogramach drzewa H przy WP70.
- Rozkładu tła na kategorie prawdy dla właściwego drzewa H i WP70 (W ma tylko sumę) — wymaga joba (histogramy masa × score × kategoria prawdy z `app/d150`).
- Metody R sieć × drzewo z nowym modelem (brak kategorii łącznych lokalnie).
- BR pozostałych mod D⁰/D⁺ w mb_nd (czy κ w danych to rzeczywiście κ_MC × 1,364).

## Rekomendacja
1. D⁺: nominał bez zmian. D⁰ przy pT < 8 GeV: M8 jako wariant/kandydat na nominał (mniejszy rozrzut, brak deficytu DMM), z błędem z bootstrapu; przy 8–30 GeV zostaje M0 (w 428770 model odrzucony).
2. Kalibracja D⁰ na inkluzywnym piku pozostaje zablokowana; zostaje tag D*. Do odblokowania 5–8 / 8–30 GeV wystarczyłoby ~1,3× więcej mb_nd, do 3–5 GeV ~6×.
3. Zgłosić autorom próbek BR ≠ PDG; składową DMM w fitach D⁰ normalizować z PDG, nie z tablicy generatora.

## Pliki
`KRYTERIA.md`, `RAPORT.md`; `skrypty/`: `fb_fit.py` (moduł), `fb_diag_T.py`, `fb_diag_T2.py`, `fb_diag_W.py`, `fb_diag_sig.py`, `fb_rep_diag.py`, `fb_rep_decomp.py` (diagnoza), `fb_var.py` (asimov | mc | toys | data), `fb_select.py`, `fb_rep_data.py`;
`tabele/`: `diag_T.md`, `diag_T2.md`, `diag_W.md`, `diag_sig.md`, `decomp.md`, `select.md`, `data.md`; `out/`: json + logi (`var_data_k1.json` = κ bez skali BR; `fitbias_override_copy.json`).
