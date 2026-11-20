# FIT_scorebins — jednoczesny fit masy w kategoriach score H vs pojedyncze cięcie (projekt B, 17.11.2026)

## 0. Kryteria „warto” — zapisane PRZED obejrzeniem wyników (17.11.2026, przed pierwszym uruchomieniem pseudo-eksperymentów)

Miara: σ(N_S)/N_S z pseudo-eksperymentów (odchylenie standardowe N̂_S), zysk g = 1 − σ_kat/σ_cięcie(najlepszy WP), sparowane toye.

1. **WARTO (rekomendacja wdrożenia)**, jeśli jednocześnie:
   - g ≥ 10 % (stat.) w co najmniej dwóch z trzech binów pT (3–5, 5–8, 8–30) dla danego gatunku, istotnie > 0 (g/σ_g ≥ 3);
   - |obciążenie N_S| < 0,2 σ_stat w nominale (średni pull), pokrycie 1σ w 68 ± 4 %;
   - przy ε_S per kategoria przesuniętych jak zmierzone SF (0,83–1,05) obciążenie konfiguracji rekomendowanej (z ułamkami
     ograniczonymi Gaussem lub po poprawce SF) < 0,5 σ_całk i **nie większe** niż dla pojedynczego cięcia z tym samym SF;
   - przy rzeźbieniu tła (R_sculpt 0,89–0,97) obciążenie nie większe niż dla pojedynczego cięcia (różnica < 0,2 σ).
2. **WARTO WARUNKOWO** (tylko jako kontrola krzyżowa / wariant), jeśli g ∈ [5, 10) % albo g ≥ 10 % tylko w jednym binie,
   albo zysk istnieje tylko przy ułamkach ustalonych z MC (bez odporności na SF).
3. **NIE WARTO**, jeśli g < 5 % wszędzie albo zysk statystyczny jest mniejszy niż dodatkowa systematyka (ε_S per kategoria).
4. Osobno oceniamy błąd **całkowity** (stat ⊕ niepewność SF 5 % i 10 %): cięcie z SF vs kategorie z ułamkami ograniczonymi/swobodnymi.
   Tu próg „warto” ten sam: ≥ 10 % mniejszy błąd całkowity przy obciążeniu < 0,2σ.
5. Ślepota: w danych używane wyłącznie pasma boczne (|m − m0| > 75 MeV D⁰, > 70 MeV D⁺) — kształt i poziom tła; piku w danych nie fitujemy.

(Wyniki — poniżej, dopisane po uruchomieniu.)

---

## 1. Czy to już było robione?

- **Pojedyncze D (D⁰, D⁺) — nie.** W `charm-bkg-notes` (podejscia/, wyniki/, 00_INDEX.md), `charmbkg/docs` i `charm_work/*` nie ma fitu
  jednoczesnego masy w kategoriach score dla kandydatów inkluzywnych. Wszystkie fity pojedynczych D są „cięcie WP → fit 1D masy”.
- **Pary D–D (projekt A) — tak, w innej formie**: `fit2d` par (m7: SS/SB/BS/BB) używa kategorii score nóg ([0,30), [30,45), [45,65) % ε_S z MC)
  i tagu D* jako osobnej kategorii (`DOCS/double_charm.md` p. 4, `DOCS/history.md` 26.09, `S9_pairs_ad/README_S9.md`: „fit2d z kategoriami score
  wyczerpuje informację”). Tam nie porównywano zysku względem pojedynczego cięcia dla pojedynczych D ani odporności na SF per kategoria.
- Pokrewne, ale inne: fit pass/fail Δm z tagu D* (pomiar SF), KF (sufit zysku z tła znanego: σ_N ≥ √(S+B)).
- Brakowało więc wszystkiego z p. 1–4 planu; zrobione poniżej.

## 2. Metoda (skrypty w tym katalogu)

| plik | rola |
|---|---|
| `fsb_prep.py` | wejścia per gatunek × bin pT → `out/inputs_*.npz/json` |
| `fsb_common.py` | silnik: binowany fit Poissona (Fisher scoring + LM), tryby `fixed` / `free` / `constr` |
| `fsb_run.py` | Asimov + optymalizacja granic (programowanie dynamiczne), 1000 pseudo-eksperymentów, scenariusze SF i rzeźbienia → `out/results_*.json`, `out/toys_*.npz` |
| `fsb_toys_sf.py` | 500 toyów z **losowaną prawdziwą wydajnością** (błąd całkowity stat ⊕ SF) → `out/toys_sf.json` |
| `fsb_report.py` | `out/tables.md` (pełne tabele), `fig1_sigma_vs_wp.png`, `fig2_robustness.png` |

**Próbki** (największa lokalna statystyka z prawdą MC i danymi):
- D⁰: `B_trig/runs/lxplus/bt_*.npz`, klucze `all_*` = wszyscy kandydaci D⁰ w oknie 1,75–2,00 GeV, pT ≥ 3 GeV, score H nominalny;
  dane run 428770 (30 % zdarzeń, 4,42 M kandydatów), sygnał fd_nd (prawda 0/1: 216 k; refleksje K↔π prawda 100/101: 105 k), mb_nd (30 %).
- D⁺: `DLZ_dplus/frames/*_Dplus.pkl` + `DLZ_dplus/out/Dplus_mc/scores_H.npz` (H, OOF, po wecie D*), okno 1,72–2,02 GeV;
  dane 428770 2,79 M kandydatów, sygnał fd 57 k.
- **Ślepota**: z danych wzięto wyłącznie pasma boczne (|m − m0| > 75 MeV D⁰, > 70 MeV D⁺); piku w danych nie fitowano ani nie oglądano.

**Model prawdy toyów.** 17 komórek score o granicach z kwantyli sygnału MC (ε_S skumulowana 100, 95, 90, …, 30, 20, 10, 0 %), per bin pT.
Sygnał: histogram masy MC per komórka (D⁰: + refleksja o stałym stosunku do sygnału z MC). Tło: exp(pol2) dopasowane do pasm bocznych
danych **osobno w każdej komórce** i ekstrapolowane pod pik (χ²/ndf SB 0,5–1,7; najbrudniejsza komórka 3–5 GeV: 4,2 D⁰ / 3,3 D⁺ — exp(pol2)
nie wystarcza przy ~10⁶ wpisów). Normalizacja: B z SB danych; **N_S = S_MC(mb, prawda 0/1) × N_SB(dane)/N_SB(mb)** (źródło: `out/inputs_*.json`),
czyli S/B jak w MC minbias: N_S = 6365 / 7244 / 8843 (D⁰) i 1068 / 1199 / 1724 (D⁺) w próbce lokalnej. Kontrola względem istniejącego fitu
danych (`HYP_d0/nominal_sigma_check.md`: σ_rel WP70 keep_both 7,6 / 7,0 / 3,1 % na tych samych 30 % runu): u mnie WP70 daje 8,6 / 4,1 / 3,3 % —
zgodne w 3–5 i 8–30, w 5–8 moja normalizacja jest optymistyczna ×1,7 (zysk względny jest na to mało czuły — skan S×0,5…×2 niżej).

**Fity.** Sygnał = szablon MC kategorii (kształt ustalony — tak samo w cięciu i w kategoriach), tło exp(pol2) swobodne (3 parametry) w każdej kategorii.
(a) cięcie: jedna kategoria „score > próg WP”, N_S = N_fit/ε_MC, WP 95…30 co 5; (b) K kategorii pokrywających **cały** zakres score, wspólne N_S,
ułamki z MC; granice optymalizowane (DP po 17 komórkach, Asimov); (c) ułamki swobodne (N_S = suma) oraz ograniczone Gaussem: skumulowana
wydajność na każdej wewnętrznej granicy E_k → E_k(1+δ_k), δ_k ~ G(0, σ_SF), pokrycie całości = 1; (d) fit 2D nie-binowany nie był robiony —
jego granicą jest K = 17 (wszystkie komórki), co pokazuje nasycenie. Suma tła exp(pol2) z komórek nie jest dokładnie exp(pol2) → realistyczne,
małe niedopasowanie modelu tła w szerokich kategoriach i luźnych cięciach jest w toyach obecne.

## 3. Wyniki

### 3.1 Zysk statystyczny (1000 toyów sparowanych; błąd zysku z bootstrapu toyów) — rys. `fig1_sigma_vs_wp.png`

| gatunek, pT | najlepsze cięcie: σ/N_S | K=2 | K=3 | K=5 | K=8 | K=17 (≈2D) | cięcie WP70 → K=5 |
|---|---|---|---|---|---|---|---|
| D⁰ 3–5 | WP40: 5,35 % | +0,6 ± 0,3 | +11,0 ± 1,6 | **+14,8 ± 1,3** | +16,9 ± 1,4 | +17,2 ± 1,5 | 8,60 → 4,55 % |
| D⁰ 5–8 | WP55: 3,90 % | +0,3 ± 0,3 | +5,8 ± 1,3 | **+6,5 ± 1,3** | +7,1 ± 1,2 | +5,8 ± 1,5 | 4,12 → 3,64 % |
| D⁰ 8–30 | WP60: 3,18 % | +0,9 ± 0,5 | +5,6 ± 1,0 | **+5,5 ± 1,1** | +5,3 ± 1,3 | +6,5 ± 1,5 | 3,34 → 3,00 % |
| D⁺ 3–5 | WP35: 8,02 % | +0,2 ± 0,2 | +7,2 ± 1,9 | **+10,7 ± 1,5** | +11,8 ± 1,4 | +12,3 ± 1,6 | 11,73 → 7,16 % |
| D⁺ 5–8 | WP55: 6,12 % | −0,1 ± 0,3 | +7,7 ± 1,5 | **+9,5 ± 1,1** | +10,6 ± 1,1 | +10,7 ± 1,2 | 6,32 → 5,54 % |
| D⁺ 8–30 | WP70: 4,69 % | +0,7 ± 0,3 | +6,5 ± 1,4 | **+9,4 ± 1,4** | +10,3 ± 1,1 | +11,5 ± 1,3 | brak zmiany WP |

Zysk = 1 − σ_K/σ_cięcie [%], ułamki ustalone z MC. Uwagi:
- K = 2 nic nie daje (optimum to „najlepsze cięcie + reszta”); zysk zaczyna się od K = 3 i nasyca przy K = 5–8 (K = 5 daje ≥ 80 % tego, co K = 17).
- Względem **optymalnego** cięcia zysk jest umiarkowany (5–15 %); względem używanego dziś WP70 duży tam, gdzie WP70 jest dalekie od optimum
  (3–5 GeV: −47 % D⁰, −39 % D⁺) — ale większość tej różnicy daje samo przejście na WP35–40.
- Obciążenie nominalne (średni pull): −0,03…+0,25 σ; +0,16…+0,25 σ w 5–8 i 8–30 GeV pochodzi z niedopasowania sumy kształtów tła i jest
  **takie samo w cięciu** (D⁰ 8–30: cięcie +0,25, K=5 +0,25, K=8 +0,03, K=17 −0,07). Pokrycie 1σ 65–71 %, szerokość pulla 0,94–1,03, zbieżność 100 %.
- Skan normalizacji sygnału (Asimov, K=5): S×0,5 / ×1 / ×2 → D⁰ +15,5/+13,9/+11,6 %, +9,3/+7,5/+5,6 %, +6,4/+6,9/+5,9 %; D⁺ +13,6/+12,7/+12,0 %,
  +10,9/+9,6/+10,0 %, +11,4/+8,8/+7,4 % — wniosek nie zależy od założonego S.
- Granice K=5 (dolne krawędzie ε_S): D⁰ 3–5: 100/60/40/20/10; 5–8: 100/85/65/50/20; 8–30: 100/95/85/75/60;
  D⁺ 3–5: 100/70/50/35/20; 5–8: 100/85/75/55/30; 8–30: 100/90/75/45/20. Wariant K=3 bez najbrudniejszej komórki (ε ≤ 95 %) daje to samo co K=3 (±1 %).

### 3.2 Ułamki swobodne i ograniczone

- **Swobodne (bez SF)**: σ/N_S = 36 / 14 / 8 % (D⁰), 84 / 30 / 14 % (D⁺) — tyle, co fit bez cięcia. Informacja ze score wchodzi wyłącznie przez
  znajomość ułamków; „bez SF” = „bez klasyfikatora”. Wariant nieużyteczny jako pomiar (tylko kontrola).
- **Ograniczone Gaussem — błąd całkowity** (500 toyów z losowaną prawdziwą wydajnością, δ_k ~ G(0, σ_SF) na każdej granicy; cięcie: σ_stat ⊕ σ_SF):

| σ_SF | gatunek | pT 3–5: cięcie → K=5 (K=8) | pT 5–8 | pT 8–30 |
|---|---|---|---|---|
| 5 % | D⁰ | 7,32 → 5,39 % (−26 %) | 6,34 → 4,43 % (−30 %) | 5,93 → 3,72 % (−37 %) |
| 10 % | D⁰ | 11,34 → 7,16 % (−37 %); K=8 6,18 % | 10,73 → 5,87 % (−45 %); K=8 5,03 % | 10,49 → 4,88 % (−54 %); K=8 4,56 % |
| 30 % (brak pomiaru SF) | D⁰ | 30,5 → 13,6 %; K=8 11,1 % | 30,3 → 10,4 %; K=8 9,1 % | 30,2 → 7,5 %; K=8 6,4 % |
| 5 % | D⁺ | 9,45 → 7,32 % (−23 %) | 7,90 → 6,08 % (−23 %) | 6,86 → 4,74 % (−31 %) |
| 10 % | D⁺ | 12,82 → 8,71 % (−32 %); K=8 8,17 % | 11,72 → 7,11 % (−39 %); K=8 6,74 % | 11,05 → 6,04 % (−45 %); K=8 5,48 % |
| 30 % | D⁺ | 31,1 → 15,0 %; K=8 12,8 % | 30,6 → 13,7 %; K=8 10,7 % | 30,4 → 10,2 %; K=8 8,9 % |

  W tych toyach: obciążenie +0,03…+0,23 σ (dla σ_SF = 30 % do +0,39 σ), pokrycie 65–72 % (30 %: 68–75 %), σ_HESSE ≈ odch. std. (pull 0,93–1,03) —
  błąd z fitu jest poprawnym błędem całkowitym. Fit z ułamkami **ustalonymi** na tych samych toyach ma rozrzut tylko trochę większy przy 5 %
  (D⁰ K=5: 5,45 / 4,52 / 4,01 %), wyraźnie większy przy 10 % (7,70 / 6,72 / 6,12 %) i jest obciążony (do −1,4 σ przy K=8).
- Skąd tak duży zysk: kategorie pokrywają 100 % sygnału, więc błąd wydajności na granicy tylko **przesuwa** sygnał między kategoriami
  (częściowo się znosi i jest mierzony na miejscu), a w cięciu wchodzi w N_S w całości. Przy zmierzonych SF ± 0,10–0,13 (5–8 GeV), ± 0,04–0,05
  (8–30) i braku pomiaru w 3–5 GeV błąd SF dominuje nad statystycznym — i to on jest główną stawką, nie 5–15 % statystyki.
- **Zastrzeżenie (zależność od modelu):** liczby zakładają niezależne niepewności względne σ_SF na każdej granicy i dokładnie znane pokrycie
  (ε = 1 przed cięciem na score). Korelacje między granicami ani niepewność kształtu sygnału nie były badane.

### 3.3 Odporność — rys. `fig2_robustness.png` (Asimov; przesunięcie N̂_S/N_S względem nominału)

**(i) Prawdziwe ε_S ≠ MC** — scenariusze z `B_trig/runs/lxplus/out/btrig_sf.md` (wiersze „all”): A = 428770 5–8 (SF 0,96/0,87/0,84 przy WP50/80/95),
B = 428770 8–30 (0,94/0,97/1,01), C = 427929 5–8 (0,99/1,03/1,05); SF liniowo w ε, → 1 przy ε = 1.

| scenariusz | metoda | D⁰ 3–5 | D⁰ 5–8 | D⁰ 8–30 | D⁺ 3–5 | D⁺ 5–8 | D⁺ 8–30 |
|---|---|---|---|---|---|---|---|
| A | cięcie (najl. WP, ε z MC) | −4,0 % (−0,75σ) | −5,5 % (−1,40σ) | −6,8 % (−2,14σ) | −4,0 % (−0,50σ) | −5,4 % (−0,89σ) | −9,8 % (−2,09σ) |
| A | K=5, ułamki MC | −4,6 % (−1,00σ) | −7,0 % (−1,92σ) | −9,9 % (−3,29σ) | −4,8 % (−0,67σ) | −7,7 % (−1,39σ) | −9,1 % (−2,15σ) |
| A | K=5, ogr. 10 % | −4,7 % | −6,7 % | −8,1 % | −5,1 % | −8,1 % | −7,4 % |
| A | K=5, swobodne | −4,0 %* | −1,9 % | +0,1 % | −1,9 % | −0,6 % | −1,3 % |
| B | cięcie | −6,0 % | −5,5 % | −5,1 % | −6,0 % | −5,5 % | −4,1 % |
| B | K=5, ułamki MC | −5,8 % | −4,9 % | −3,7 % | −5,7 % | −4,7 % | −4,0 % |
| B | K=5, ogr. 10 % | −5,6 % | −4,1 % | −1,8 % | −5,6 % | −4,2 % | −3,4 % |
| C | cięcie / K=5 MC | −1,0 / −0,7 % | −0,4 / +0,3 % | +0,2 / +1,6 % | −1,0 / −0,6 % | −0,4 / +0,6 % | +1,6 / +1,3 % |

  (* przy σ = 36 %, czyli −0,11σ.) Wnioski: z ułamkami MC kategorie są obciążone **tak samo jak nieskorygowane cięcie** (−4…−10 %), a w jednostkach
  własnego (mniejszego) σ_stat gorzej; ograniczenie Gaussa 10 % wokół wartości MC pochłania przesunięcie tylko częściowo (dane w brudnych kategoriach
  słabo mierzą ułamki). **Kategorie nie zwalniają z pomiaru SF** — ułamki centralne muszą być MC × SF, a ograniczenie ma szerokość błędu SF.
  Wtedy obciążenie znika z konstrukcji, a zysk to mniejsza propagacja błędu SF (§3.2).

**(ii) Rzeźbienie tła zależne od score** — tło prawdziwe w komórce × [1 − a_j·G(m; m0, 55 MeV)], a_j rośnie od 0 (score → 0) do 1 − R powyżej WP70:

| R | metoda | D⁰ 3–5 | D⁰ 5–8 | D⁰ 8–30 | D⁺ 3–5 | D⁺ 5–8 | D⁺ 8–30 |
|---|---|---|---|---|---|---|---|
| 0,89 | cięcie (najl. WP) | −3,7 % (−0,69σ) | −2,1 % (−0,54σ) | −1,6 % (−0,49σ) | −2,0 % (−0,25σ) | −1,7 % (−0,27σ) | −1,7 % (−0,37σ) |
| 0,89 | cięcie WP70 | −22,5 % (−2,61σ) | −4,2 % (−1,02σ) | −2,3 % (−0,68σ) | −13,4 % (−1,14σ) | −3,1 % (−0,50σ) | = |
| 0,89 | K=3 | −4,5 % (−0,94σ) | −3,1 % (−0,85σ) | −2,5 % (−0,82σ) | −3,3 % (−0,44σ) | −2,6 % (−0,46σ) | −1,7 % (−0,39σ) |
| 0,89 | K=5 | −5,2 % (−1,15σ) | −3,3 % (−0,90σ) | −2,5 % (−0,85σ) | −4,3 % (−0,59σ) | −2,6 % (−0,46σ) | −1,8 % (−0,42σ) |
| 0,89 | K=5 ogr. 10 % | −7,0 % | −5,1 % | −4,3 % | −5,0 % | −3,3 % | −2,4 % |
| 0,97 | cięcie / K=5 | −0,19 / −0,30σ | −0,14 / −0,24σ | −0,13 / −0,23σ | −0,07 / −0,16σ | −0,07 / −0,12σ | −0,10 / −0,11σ |

  Swobodne tło per kategoria **nie pochłania** rzeźbienia o szerokości porównywalnej z pikiem (żaden model gładki tego nie zrobi — ani w cięciu,
  ani w kategoriach). Kategorie są na nie nieco bardziej czułe niż najlepsze (ciasne) cięcie: o 0,2–0,45 σ dla D⁰ i 0,05–0,35 σ dla D⁺ przy R = 0,89,
  o ≤ 0,1 σ przy R = 0,97; dużo mniej czułe niż cięcie WP70 przy niskim pT. Wariant z ograniczeniem jest czulszy (rzeźbienie udaje migrację sygnału).

## 4. Ocena względem kryteriów z §0 i wniosek

| kryterium | wynik |
|---|---|
| zysk stat. ≥ 10 % w ≥ 2 z 3 binów | D⁰: tylko 3–5 GeV (K=5: 14,8 / 6,5 / 5,5 %) → **niespełnione**; D⁺: K=5 10,7 / 9,5 / 9,4 %, K=8 11,8 / 10,6 / 10,3 % → **spełnione dla K=8, na granicy dla K=5** |
| obciążenie nominalne < 0,2σ, pokrycie 68 ± 4 % | spełnione poza 8–30 GeV (+0,16…+0,25σ, identyczne w cięciu — to model tła, nie kategorie) |
| odporność na SF (konfiguracja rekomendowana < 0,5σ_całk, nie gorzej niż cięcie) | ułamki z MC: **niespełnione** (jak nieskorygowane cięcie, −4…−10 %); z ułamkami MC×SF i ograniczeniem — spełnione z konstrukcji |
| rzeźbienie: nie gorzej niż cięcie (< 0,2σ różnicy) | R = 0,97: spełnione; R = 0,89: **niespełnione dla D⁰** (różnica 0,3–0,45σ), D⁺ na granicy |
| błąd całkowity (stat ⊕ SF) ≥ 10 % mniejszy, obciążenie < 0,2σ | **spełnione wyraźnie**: −23…−37 % przy σ_SF = 5 %, −32…−54 % przy 10 % (obciążenie +0,03…+0,23σ) |

**Werdykt: WARTO WARUNKOWO jako zysk statystyczny (5–15 % względem optymalnego cięcia), WARTO jako sposób na błąd SF.**
Sam zysk statystyczny nie uzasadnia złożoności dla D⁰ powyżej 5 GeV. Prawdziwa wartość to to, że przy kategoriach pokrywających cały score
niepewność wydajności klasyfikatora (dziś 5–13 %, przy 3–5 GeV niezmierzona) przestaje wchodzić w N_S w całości: błąd całkowity spada o 1/4–1/2.
Najtańszy krok, niezależny od wszystkiego: **przesunąć punkt pracy przy 3–5 GeV z WP70 na WP35–40** (σ_stat 8,6 → 5,4 % D⁰, 11,7 → 8,0 % D⁺;
ale wymaga SF przy WP40, którego dziś nie ma).

**Rekomendowana konfiguracja:** K = 5 kategorii score per bin pT, pokrywających cały zakres (granice z §3.1), wspólne N_S, ułamki = MC × SF
z ograniczeniem Gaussa o szerokości błędu SF na każdej granicy, tło exp(pol2) swobodne per kategoria, kształt sygnału wspólny. K = 3 jako wersja
minimalna (daje ~2/3 zysku), K = 8 tylko dla D⁺. Wariant ze swobodnymi ułamkami — wyłącznie jako kontrola zgodności (χ² ułamków dane/MC).

**Koszt / systematyka do dodania:** 3K parametrów tła zamiast 3; SF potrzebne na K−1 granicach zamiast jednej (z tagu D*, fit pass/fail — istnieje
dla WP50/80/95); systematyka kształtu tła per kategoria (alternatywa pol3 / exp(pol3)), szczególnie w najbrudniejszej kategorii (χ²/ndf SB 3–4
przy 3–5 GeV); test rzeźbienia per kategoria; statystyka MC szablonów; korelacja SF między granicami.

**Przeniesienie do TRExFitter:**
- `Region` = kategoria score × bin pT (× gatunek), zmienna = masa (te same biny 5 MeV), `Type: SIGNAL` dla wszystkich; fit jednoczesny.
- `Sample` sygnał: histogram MC regionu; jeden `NormFactor` μ_S na bin pT wspólny dla wszystkich kategorii (D⁰: refleksja jako osobny `Sample`
  z tym samym `NormFactor`).
- Niepewność wydajności: K−1 parametrów uciążliwych „eff_granica_k”, każdy jako `OverallSys` na sygnale **antyskorelowany** w dwóch sąsiednich
  regionach (+E_k σ/φ_k w kategorii k, −E_k σ/φ_{k−1} w k−1) — migracja zachowująca sumę, dokładnie parametryzacja `constr` z `fsb_common.py`;
  wartości centralne ułamków już po SF.
- Tło: TRExFitter nie ma funkcji analitycznych — `Sample` tła z szablonu (fit exp(pol2) do SB) + `NormFactor` per region + dwa `HistoSys`
  (nachylenie, krzywizna) zdekorelowane między regionami, z szerokim/płaskim ograniczeniem; nie `ShapeFactor` per bin (zabiera całą czułość).
  Alternatywnie zostawić fit w charmtag (`fit2d` par ma już kategorie score) i do TRExFittera przekazać tylko N_S.
- `MCstatThreshold` (gammy) dla szablonów sygnału w wąskich kategoriach.

## 5. Ograniczenia

- Kształt sygnału ustalony z MC we wszystkich wariantach (bez swobodnej pozycji/szerokości piku, bez statystyki MC) — w realnym ficie wspólne
  parametry kształtu między kategoriami powinny dodatkowo sprzyjać kategoriom, ale tego nie liczyłem.
- N_S z S/B MC, nie z danych; w 5–8 GeV D⁰ prawdopodobnie zawyżone (kontrola w §2). D⁰: 30 % runu 428770; D⁺: próbka z ramek DLZ (ułamek runu nieustalony).
- Fit 2D nie-binowany (d) niezrobiony; zastępuje go granica K = 17.
- Odporność (§3.3) liczona na Asimovie, nie toyami; scenariusz SF dla 3–5 GeV przeniesiony z wyższych binów (brak pomiaru).
- Model niepewności SF (niezależne δ na granicach, pokrycie = 1) jest założeniem — od niego zależą liczby §3.2.
- Rzeźbienie: jeden kształt (dołek Gaussa 55 MeV); R_sculpt 0,89 to wartość H dla D⁺ z FEAT/DLZ.
- Refleksja D⁰ jest częściowo zawarta także w tle z SB (lekkie podwójne liczenie, wspólne dla wszystkich wariantów).

## Uzupełnienie 18.11: kalibracja kształtem score na danych D⁰←D* (R3; `R3_2026-11-18/`) — nie daje lepszej precyzji niż pass/fail
Kryteria zapisane z góry (`R3_2026-11-18/KRYTERIA.md`). Fit jednoczesny Δm w 6 kategoriach score (granice WP95/80/70/50/25 + odrzucone), Σf = 1; błędy ze sparowanego bootstrapu po zdarzeniach (200 replik); silnik odtwarza oficjalny `fit_passfail` (≤ 0,08σ).
| run | pT | SF70 kształt | SF70 pass/fail | σ_k/σ_pf | SF95 kształt |
|---|---|---|---|---|---|
| 428770 | 5–8 | 0,909 ± 0,133 | 0,913 ± 0,134 | 1,00 | 0,836 ± 0,099 |
| 428770 | 8–30 | 0,991 ± 0,041 | 0,986 ± 0,043 | 0,97 | 1,008 ± 0,035 |
| 427929 | 5–8 | 0,995 ± 0,123 | 0,995 ± 0,131 | 0,94 | 1,045 ± 0,111 |
| 427929 | 8–30 | 0,983 ± 0,097 | 1,002 ± 0,106 | 0,91 | 0,913 ± 0,069 |
- (b) precyzja: niespełnione (próg 0,8); obie metody skorelowane 0,95–0,99 — informacja jest ta sama, ogranicza statystyka D* i tło w kategorii odrzuconej. (c) zgodność: spełnione. (a) wstrzyknięcie: obciążenie ≤ 3 % poza 428770 5–8 (+3,4 ± 0,8 %, to samo w pass/fail → model fitu). 3–5 GeV: bez wartości pomiarowej w obu runach.
- Macierz korelacji SF między WP (pełna w `R3_2026-11-18/verdict.json`): ρ(95,80) 0,78–0,87; ρ(95,70) 0,75–0,84; ρ(80,70) 0,94–0,97; ρ(70,50) 0,86–0,97 → zastępuje założenie niezależnych δ z §3.2.
- Post hoc, deficyt 5–8 GeV w 428770: w kształcie nieistotny (χ² 7,2/5, p = 0,21); dane/MC ułamków sygnału: 3,9 ± 1,8 (< WP95), 0,70 ± 0,40, 0,56 ± 0,24, 0,80 ± 0,16, 0,91 ± 0,14, 0,99 ± 0,16 (górne 50 %). „0,83” to sygnał przypisany kategorii odrzuconej (S/B ≲ 1,5 %), gdzie decyduje model tła.
- Znalezisko uboczne, 8–30 GeV: ułamek sygnału przy najwyższym score (> WP25) dane/MC 0,72 ± 0,05 (428770, −5,5σ), 0,74 ± 0,12 (427929); potwierdzone zliczeniem w oknie. Przyczyna niezbadana; dotyczy WP ciaśniejszych niż 50.
- Zastrzeżenia: toye mają tło z tej samej rodziny funkcji co fit; bez klas triggera; błąd HESSE zaniżony 0–15 % (brać bootstrap).

### Uzupełnienie 18.11 (15:50): źródło deficytu przy najwyższym score w 8–30 GeV (`HIGHSCORE_2026-11-18/`, eksploracyjne, post hoc)
- Ok. 2/3 efektu w 428770 to widmo pT wewnątrz szerokiego binu: sygnał D* w danych J12 ma 64 % powyżej 12 GeV (MC 28 %), a f(>WP25) w MC przy jednym progu spada z 0,283 (8–10 GeV) do 0,082 (20–30 GeV). ρ50 = [N(>WP25)/N(>WP50)] dane/MC: 8–12 GeV 1,01 ± 0,06; 12–30 GeV 0,655 ± 0,053. Po przeważeniu MC do widma danych kształt score zgodny (χ² 97/13 → 16,5/13).
- Reszta: 0,905 ± 0,048 (428770), 0,76 ± 0,10 (427929; tam widma pT zgodne), razem ≈ 0,88 ± 0,04 (−2,8σ lokalnie) — nierozstrzygnięte.
- Wykluczone jako główna przyczyna: okno/pasma Δm, model tła, klasa triggera, pile-up, udział nieprompt (wymagałby ×35), swap K↔π, przetrenowanie, pojedyncza cecha H.
- Wniosek praktyczny: SF dla WP ciaśniejszych niż 50 w 428770 liczyć w wąskich binach pT albo z MC przeważonym do widma danych.
