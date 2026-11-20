# FORCED — f_acc (AFP-MIX v2) w 428770 po wyłączeniu „tagów wymuszonych”, 18.11.2026

Projekt B. Kryteria zapisane przed liczeniem: `KRYTERIA.md`. Praca lokalna (ok. 15 min liczenia),
istniejące pliki bez zmian. W regionie tagowanym tylko frakcje tagu przypadkowego, liczniki i bity (żadnych fitów mas, uzysków, przekrojów).

## 0. Krótko

1. **Odtworzenie:** `run_data.py` bez zmian (`--samples data_428770`) daje sekcję 428770 raportu v2 identyczną co do znaku (`diff` pusty):
   N_tag 16065, N_acc(A) 766,1 ± 25,3, B(i) 1160,7 ± 218,2, B(ii) 857,7 ± 50,0, f_acc 0,0477 ± 0,0135, bramki G0–G8 PASS, G9 N/A.
2. **Nominalny f_acc (estymator A, klasa „z D”) się nie zmienia:** 0,0477 → 0,0473 (W1) / 0,0476 (W2); Δ = −0,0004 / −0,0001 = −0,03 σ / −0,005 σ
   pełnego błędu nominalnego (−0,22 / −0,04 σ_stat). Kryterium |Δ| < 0,5 σ spełnione w obu wariantach z zapasem ×19. **Werdykt: bez zmian.**
3. **Systematyki nie trzeba powiększać** (dodanie |Δ| w kwadracie: +0,03 % błędu). Odwrotnie: tagi wymuszone zawyżały estymator B(i) —
   9 z 43 otagowanych dawców puli B(i) to tagi wymuszone. Po wyłączeniu B(i)/A spada z 1,52 ± 0,29 do 1,13 ± 0,23, składnik „½Δ metod” z 25,8 % do 6,3–7,0 %,
   a błąd całkowity z 28,4 % do 13,2 % (0,0135 → 0,0063). To zmiana błędu, nie wartości — decyzja, czy ją przyjąć, pozostaje otwarta (zastrzeżenie 1).
4. Wzrost częstości tagu przy nSelBL < 5 względem 40–60: ×30,8 ± 0,9 → ×21,9 ± 0,7 (W1), ×22,1 ± 0,8 (W2) — potwierdza liczbę przeglądu (×30 → ×22).
5. 427929: 0 zdarzeń próbki mieszania zmienia status (87 zdarzeń tagowanych bez łańcucha nie-AFP z S w całym drzewie, żadne w próbce „trygery nie-AFP”).

## 1. Definicje wykluczenia (`f1_define.py`, `out/f1_428770.log`, `out/t1_liczniki_428770.csv`)

Klucz: `event_number` — jednoznaczny w `data_428770_events.parquet` (7 989 496) i w `bits_428770.parquet` (1 338 608); wszystkie zdarzenia z bitów mają wiersz
w tabeli zdarzeń; LB, nSelBL, npu i 4 bity kontrolne zgodne dla wszystkich. Zbiór zdarzeń z bitami = zdarzenia z `afp_tag_A ∨ afp_tag_C` (dokładnie).
Zbiór S z `TAGREG_2026-11-18/v2/out/strumien_S_428770.csv`. Wykluczenie działa u biorców i dawców (zdarzenie traci status „tryger nie-AFP”).

| wariant | zdarzenia tagowane | zdarzenia nietagowane (brak pełnych bitów) | wykluczone z próbki (5 449 549) | z nich tagowane | z D | z D ∧ tag |
|---|---|---|---|---|---|---|
| W1a (definicja przeglądu) | pełne 75 bitów: brak łańcucha nie-AFP z S ∧ jest łańcuch AFP z S | bez zmian | 1 148 | 1 148 | 129 | 129 |
| **W1** | jak W1a | wykluczone, jeśli brak bitu J12 / j20_L1MBTS_2 ∧ L1_AFP_A_OR_C przed preskalą = 1 (mogły wejść łańcuchem AFP) | 1 393 | 1 148 | 208 | 129 |
| **W2** (symetryczny) | wymóg bitu łańcucha nie-AFP z S wśród 11 zapisanych: `noalg_L1J12` ∨ `j20_L1MBTS_2` (RD0 = 0 zdarzeń) | ten sam wymóg | 99 257 (1,8 %) | 1 572 | 71 567 | 338 |

- W1 jest dokładny dla tagowanych, przybliżony dla nietagowanych: 245 zdarzeń nietagowanych (79 z D; 158 z ≥ 1 protonem AFP) to górna granica zdarzeń
  zapisanych łańcuchem AFP bez tagu offline. W1a i W1 dają te same liczby do 4. cyfry.
- W2 traktuje tagowane i nietagowane jednakowo (bez asymetrii licznik/mianownik), kosztem 1,8 % próbki; traci też 424 zdarzenia tagowane (209 z D)
  z prawdziwym łańcuchem nie-AFP z S spoza 11 bitów (j10f_L1MBTS_2, EM, MU) i 97 685 nietagowanych bez bitu S13. Zawiera wszystkie 1 148 wymuszonych.
- Tagowanych bez żadnego łańcucha S w próbce mieszania: 0.

## 2. Wyniki (`f2_run.py` → `W1a/`, `W1/`, `W2/`; `f3_table.py` → `out/t3_facc_warianty.csv`, `out/t4_gatunki_pik.csv`)

`f2_run.py` importuje `run_data.py` bez zmian i opakowuje tylko `trigger_classes` (maska wariantu). Błędy jak w v2 (bootstrap 200, to samo ziarno).

| wielkość (A∨C, po cięciu ξ_cal) | nominał (= raport v2) | W1 (W1a identyczny) | W2 | Δ W1 / W2 |
|---|---|---|---|---|
| dawcy / biorcy | 5 449 549 / 4 475 334 | 5 448 156 / 4 475 126 | 5 350 292 / 4 403 767 | |
| N_tag (biorcy „z D”) | 16 065 | 15 970 | 15 833 | −95 / −232 |
| N_acc (A) | 766,1 ± 25,3 | 755,9 ± 27,7 | 754,0 ± 24,6 | |
| **f_acc (A) ± razem** | **0,0477 ± 0,0135** | **0,0473 ± 0,0063** | **0,0476 ± 0,0061** | **−0,03 σ / −0,005 σ** (σ = 0,0135); −0,22 / −0,04 σ_stat |
| f_acc (A) strona A / C | 0,0455 / 0,0500 | 0,0449 / 0,0499 | 0,0452 / 0,0503 | |
| B(i): N_acc; f | 1160,7 ± 218,2; 0,0722 ± 0,0136 | 850,7 ± 177,9; 0,0533 ± 0,0111 | 859,9 ± 205,1; 0,0543 ± 0,0130 | −1,40 σ / −1,32 σ (σ B(i) nom.) |
| B(i): otagowani dawcy / pula | 43 / 23 239 | 34 / 23 225 | 25 / 17 966 | |
| B(ii): N_acc; f | 857,7 ± 50,0; 0,0534 ± 0,0031 | 846,2 ± 53,3; 0,0530 ± 0,0033 | 838,9 ± 52,1; 0,0530 ± 0,0033 | −0,13 σ / −0,13 σ |
| błąd wzgl.: stat / pula / ekstrapolacja / ½Δ metod | 3,3 / 4,8 / 10,4 / 25,8 % | 3,7 / 4,7 / 10,0 / 6,3 % | 3,3 / 4,5 / 9,3 / 7,0 % | |
| błąd wzgl. razem | 28,4 % | 13,2 % | 12,9 % | |
| G1 p / G5 p | 0,996 / 0,756 | 0,998 / 0,823 | 1,000 / 0,844 | |
| G2 (rozłączne klasy) obs/pred, z | 0,931, −1,67 | 0,974, −0,56 | 0,963, −0,76 | |
| G7 = B(ii)/A, z | 1,120, +2,04 | 1,119, +1,95 | 1,113, +1,85 | |
| G8 = B(i)/A, z | 1,515, +1,79 | 1,125, +0,54 | 1,140, +0,52 | |
| G9 obs / pred (N/A) | 10 / 4,08 | 10 / 4,02 | 10 / 3,97 | |
| werdykt bramek | PASS (G9 N/A) | PASS (G9 N/A) | PASS (G9 N/A) | status żadnej bramki się nie zmienia |

- Gatunki × pT (pik, bez fitu): największa zmiana f_acc w którymkolwiek z 21 przedziałów = 0,29 σ_stat (D⁺ 3–150 GeV: 0,0507 → 0,0503); `out/t4_gatunki_pik.csv`.
- Poprawiają się domknięcia G2 (−1,7 → −0,6) i G8 (+1,8 → +0,5): obie bramki używają zdarzeń o małej aktywności / klas MBTS, gdzie siedzą tagi wymuszone.
  G7 (obszar kontrolny, biorcy „z D”) praktycznie bez zmian — nadwyżka +12 % zostaje.

### Częstość tagu vs nSelBL (`out/t2_ptag_nselbl_428770.csv`; tag = afp_tag_A ∨ C, próbka mieszania)

| wariant | N_tag nSelBL < 5 | P(< 5) | P(5–10) | P(40–60) | wzrost < 5 / 40–60 | wzrost 5–10 / 40–60 |
|---|---|---|---|---|---|---|
| nominał | 1 772 | 0,0578 | 0,0302 | 0,00187 | **×30,8 ± 0,9** | ×16,1 |
| W1a / W1 | 1 232 | 0,0409 / 0,0410 | 0,0275 | 0,00187 | **×21,9 / ×22,0 ± 0,7** | ×14,7 |
| W2 | 1 219 | 0,0415 | 0,0276 | 0,00188 | **×22,1 ± 0,8** | ×14,7 |

Tagi wymuszone tłumaczą ok. 30 % wzrostu przy nSelBL < 5; pozostaje ×22 (w 427929, czysty MBTS: ×68 — tam wzrost jest fizyczny, zdarzenia dyfrakcyjne mają małą krotność).

## 3. 427929 (`out/f1_427929.log`)

Sprawdzenie niepotrzebne: z 703 247 zdarzeń tagowanych próbki mieszania **0** zmienia status (w całym drzewie 87 zdarzeń tagowanych ma tylko łańcuch AFP — już wykluczone w v2 jako „wyłącznie L1AFP”).

## 4. Wniosek

- Wg kryterium (|Δ| < 0,5 σ, oba warianty): **nominalny f_acc w 428770 bez zmian** — 0,048 ± 0,014 pozostaje poprawny.
- Systematyki nie powiększać. Jeśli wykluczenie zostanie przyjęte jako część definicji próbki (rekomendowane W1 — dokładne dla tagowanych, W2 potwierdza),
  wynik to 0,047 ± 0,006; zmienia się wtedy tylko składnik „½Δ metod”.
- Zdanie do noty (§6.3.3, zamiast końcówki proponowanej w SKIM_review §3.1, l. 808–812): „Removing these events from targets and donors changes the nominal estimate
  from 0.0477 to 0.0473 (0.0476 if a stored non-AFP chain of the skim set is required for all events), i.e. by less than 0.03 of its uncertainty.
  The forced tags did bias the cross-check B(i), whose pool is dominated by low-activity events: 9 of its 43 tagged donors are forced, and without them
  B(i)/A goes from 1.52 ± 0.29 to 1.13 ± 0.23. The quoted uncertainty, which is dominated by half of this difference, is therefore conservative.”
  Dla l. 819–821: „without them the rise relative to n_sel^BL 40–60 is a factor of 22 instead of 31” (dokładnie 30,8 → 21,9).

## 5. Zastrzeżenia

1. Zmniejszenie błędu (0,0135 → 0,0063) opiera się na B(i) z 34 (W1) / 25 (W2) otagowanymi dawcami: B(i)/A = 1,13 ± 0,23, czyli ½Δ = 6 % jest poniżej
   własnej precyzji (~12 %). Nominalne 25,8 % też było efektem 1,8 σ. Ostrożnie: zostawić 0,048 ± 0,014 i napisać, że błąd jest konserwatywny.
2. Dla zdarzeń nietagowanych pełnych bitów nie ma. W1 używa bitu L1_AFP przed preskalą (górna granica 245 zdarzeń); W2 omija problem wymogiem symetrycznym,
   ale usuwa też zdarzenia zapisane przez j10f/EM/MU. Zgodność W1 i W2 (różnica f_acc 0,0003) pokazuje, że wybór nie ma znaczenia.
3. Zbiór S pochodzi z TAGREG v2 (kryterium względne, wyznaczone w regionie tagowanym); nazwa strumienia/skimu nadal nieznana.
4. Zakres metody bez zmian: tylko zdarzenia z trygerów nie-AFP (4,4 % tagów po GRL). Nadwyżka G7 (+12 %, 2 σ) i G9 (10 vs 4) nie są związane z tagami wymuszonymi.
5. Nie liczono żadnych fitów ani uzysków; „z D” = zdarzenie z ≥ 1 kandydatem slim.

## 6. Pliki (`~/Documents/charm/charm_work/AFP_mix/FORCED_2026-11-18/`)

- `KRYTERIA.md`, `RAPORT.md`; skrypty: `f1_define.py [run]`, `f2_run.py <W1a|W1|W2>`, `f3_table.py` (`~/Documents/charm/.venv_h5/bin/python`).
- `nominal/`, `W1a/`, `W1/`, `W2/`: `RAPORT_DANE_v2.md`, `result_v2.json`, `stdout.log`.
- `out/`: `excl_428770.parquet` (maski), `f1_428770.log`, `f1_427929.log`, `f3.log`, `t1_liczniki_{428770,427929}.csv`, `t2_ptag_nselbl_{428770,427929}.csv`,
  `t3_facc_warianty.csv`, `t4_gatunki_pik.csv`.
