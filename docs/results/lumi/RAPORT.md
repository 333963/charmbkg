# Świetlność scałkowana runów 427929 i 428770 (projekt B), 19.11.2026

Zakres: tylko odczyt (COOL przez Frontier, cvmfs, lokalne tabele TAGREG); bez condora, bez zmian w ORCH/etapach/nocie. W regionie tagowanym użyto wyłącznie liczników zdarzeń.
Oznaczenia: **[Z]** = założenie, **[R]** = mój rachunek.

## 0. Wynik

| run | łańcuch / region | preskala | L [nb⁻¹] | błąd | LB |
|---|---|---|---|---|---|
| 427929 | `HLT_mb_mbts_L1MBTS_1` ≡ `HLT_noalg_mb_L1MBTS_1` | 1 | **0,4737** | ± 0,0105 (2,2 % [Z]) | 119–324, 189 LB |
| 427929 | to samo, nasza kopia plików (ubytki plików w 16 LB) | 1 | 0,4648 | jw. + kompletność [Z] | jw. |
| 428770 | R3 `HLT_noalg_L1J12` (cały run w GRL, bez preskali) | 1 | **29,41** | ± 0,68 (2,3 %) | 4–25, 210–323, 136 LB |
| 428770 | R1 `HLT_noalg_L1AFP_A_OR_C` | 60 (120 w LB 4–16) | **0,4687** | ± 0,0104 (2,2 %) | jw. |
| 428770 | `HLT_j20_L1AFP_A_OR_C` (sam łańcuch) | 1 (2 w LB 4–16) | 28,12 | ± 0,62 | jw. |
| 428770 | R2 `j20_L1AFP_A_OR_C ∧ ¬R1` | jw. × (1 − 1/60) | **27,65** | ± 0,61 | jw. |

Tag świetlności: `OflLumi-Run3-003`. Jednostka nb⁻¹. Pełny zapis: `lumi.json`; per LB: `tabele/lumi_final_perLB_<run>.csv`.

L = Σ_LB L_del × live / preskala. Składniki:
- **L_del per LB**: oficjalne pliki `ilumicalc_histograms_None_427929-428770_OflLumi-Run3-003*.root` z cvmfs (`GroupData/GoodRunsLists/data22_13p6TeV/20231005` i `20231024`, GRL `PHYS_StandardGRL_All_Good_25ns_AFP_SIT`), `Inst_m_Lumi × DeltaT`. Własny przebieg `iLumiCalc` (Athena 24.0.151) na GRL analizy daje tę samą sumę (31 135,1 μb⁻¹ = 0,4755 + 30,660 z LB 324).
- **Czas martwy**: `iLumiCalc` w tym wydaniu nie znajduje pozycji L1 ani łańcuchów HLT Run 3 w COOL (`/TRIGGER/LVL1/Menu` pusty dla tych runów), więc ułamek żywy policzyłem wprost z COOL `/TRIGGER/LUMI/LVL1COUNTERS` (TAV/TAP per LB; CTP id z TriggerDB, SMK 3094). Nominalnie `L1_RD0_FILLED` (losowy, 5·10⁷ zliczeń): 427929 **0,99625**; 428770 **0,9665** (LB 4: 0,13; LB 17–20: 0,87–0,89).
  Inne pozycje: 427929 MBTS_1 0,99626, AFP 0,99639, J12 0,9934; 428770 AFP 0,9692, J12 0,9591, EM22VHI 0,9569. Rozrzut wchodzi do błędu (0,3 % / do 0,8 %).
  Oficjalne pliki używają `L1_EM22VHI`; w 427929 to 1374 zliczenia na 189 LB, więc ich „recorded” (0,4619 nb⁻¹, LB 265 z zerem) jest nieużywalne per LB.
- **Preskale**: z drzewa = L1 × HLT z COOL we wszystkich LB (TAGREG v2 §4.1).
- **LB**: `GRL_MinBiasAFP_lowmu2022.xml` (gałąź minimum-bias CharmAnalysis) ∩ LB obecne w naszych plikach. 427929: 119–140, 150–154, 157–163, 165–167, 170–223, 226–279, 281–324 (189, wszystkie obecne). 428770: 4–25, 210–323 (136); LB 324 jest w GRL, ale nie ma go w plikach (L_del 0,230 nb⁻¹, live = 0 w pliku oficjalnym).

## 1. Skąd 29,41 nb⁻¹

Źródło liczby: `charmpp`, `config/trex/minbias_trex_example.yaml:15`, `lumi_nb: 29.41`, opis „PLACEHOLDER: recorded GRL luminosity without trigger/prescale” (GITLAB_2026-11-17/RAPORT.md).
Odtworzona: L_del(428770, 136 LB, tag 003) × live(RD0_FILLED) = 30,4297 × 0,96650 = **29,4105 nb⁻¹**. To świetlność zapisana runu 428770 w GRL bez preskali; dotyczy tylko łańcuchów o preskali 1 (`noalg_L1J12`). Dla 427929 nie ma zastosowania.
Zgodność do 4 cyfr traktuję jako niezależne potwierdzenie naszej liczby dla 428770 (tam ktoś policzył ją narzędziem z poprawnym czasem żywym).

## 2. Sprawdzian krzyżowy MBTS (427929)

- N(L1_MBTS_1 przed preskalą, 189 LB) = 28 914 481; L_del = 0,4755 nb⁻¹ → σ_wid = **60,8 mb** (per LB 59,0–65,5, rozrzut 1,3 %).
- Niezależnie od świetlności: ułamek przecięć z MBTS_1 = 28,91·10⁶ / 8,143·10⁹ (RD0_FILLED) = 0,355 %; przy ⟨μ⟩ = 0,00469 (σ_inel = 80 mb) daje wydajność 0,757.
- Oczekiwanie: σ_fid(ξ > 10⁻⁶) = 68,1 ± 1,4 mb przy 13 TeV, σ_inel = 78,1 ± 2,9 mb (ATLAS, arXiv:1606.02625; wydajność trygera MBTS 99,7 %); przy 13,6 TeV ok. +1 % → 69–72 mb **[Z]** (liczby 68,1 mb nie potwierdziłem dziś u źródła, 78,1 mb tak).
- Stosunek 0,85–0,88: świetlność z licznika MBTS wyszłaby 0,40–0,42 nb⁻¹, czyli **12–16 % mniej** niż z narzędzia. Mieści się w oknie 10–20 %, nie w 10 %. Kierunek: albo wydajność MBTS_1 w Run 3 niższa niż w 2015, albo kalibracja tagu 003 dla tego runu za wysoka. Z tagiem 002 różnica byłaby większa (0,507), z 001 mniejsza (0,461).
- W 428770 sprawdzian MBTS nie działa: MBTS_1 odpala w 12,5 % przecięć przy μ = 0,045 (σ „widzialny” 222 mb), czyli sygnały spoza przecięcia w pociągach paczek (590 paczek wobec 64 w 427929).

## 3. Niepewność

- Kalibracja: przyjąłem **2,2 %** (wstępna wartość dla danych pp 2022, ATL-DAPR-PUB-2023-001; nowsza nota 2024 podaje dla 2022 mniej) **[Z]** — dokument dotyczy danych fizycznych przy dużym μ; dla runów niskiego μ nie znalazłem osobnej liczby.
- Rozrzut tagów dla tych samych LB: 427929 0,461 / 0,507 / 0,4755 (001 / 002 / 003); 428770 32,27 / 32,71 / 30,66. Tag 004 = 003; 005–007 dają 0. Konserwatywnie ±7 % do czasu potwierdzenia przez grupę świetlności.
- Czas żywy: ±0,3 % (427929), ±0,8 % (428770, J12).

## 4. Kompletność naszych plików

- 427929: N_drzewo(`noalg_mb_L1MBTS_1`) / TAV(L1_MBTS_1) = 0,828; w 173 z 189 LB płaskie 0,84 (maks. 0,855, żaden LB nie dochodzi do 1), w 16 LB spadki do ~0,70 (= 5/6 plateau), LB 271: 0,57, LB 272: 0,28.
  Interpretacja **[Z]**: plateau 0,84 to wymóg na poziomie zdarzenia przy produkcji ntupli (w `cp_config_minbias_base.yaml` selekcja wierzchołka pierwotnego nie jest wyłączona, w `cp_config_spg.yaml` jest), czyli wydajność, nie świetlność; spadki to brakujące pliki. Wtedy kompletność plików w GRL = 0,981 i L(kopia) = 0,4648 nb⁻¹. Gdyby plateau było ubytkiem plików, L(kopia) = 0,392 nb⁻¹.
  Do rozstrzygnięcia jednym zliczeniem: liczba zdarzeń drzewa bez wierzchołka pierwotnego (powinna być 0).
  Brakujące 125 z 583 plików leżą więc głównie poza LB z GRL.
- 428770: N_drzewo(`noalg_L1J12`) / TAV(L1_J12) = 0,966, płaskie (0,933–0,985); kompletność plików 0,998.
- Uwaga dla R1: N_drzewo(`noalg_L1AFP_A_OR_C`) / (TAV(L1_AFP)/60) = 0,44 (0,42–0,47 per LB). Preskala 60 jest potwierdzona, więc 56 % zdarzeń tego łańcucha nie ma w drzewie — ten sam wymóg na poziomie zdarzenia **[Z]**. Dla przekroju w R1 to wydajność selekcji zdarzeń, którą trzeba mieć pod kontrolą (dla zdarzeń z D bliska 1, jeśli to wymóg wierzchołka).

## 5. 428770: czy jedna świetlność na region

- R1: tak, 0,4687 nb⁻¹ (preskala zmienia się 120 → 60 po LB 16, ale to tylko waga LB; 4,6 % świetlności R1 pochodzi z LB 4–16).
- R2: tak, 27,65 nb⁻¹, pod warunkiem definicji bitowej `j20_L1AFP_A_OR_C ∧ ¬noalg_L1AFP_A_OR_C` (czynnik 1 − 1/60 jest czystą preskalą).
- R3: 29,41 nb⁻¹ dla samego łańcucha; weto ¬R2 jest cięciem fizycznym zależnym od preskali L1_AFP (2 w LB 4–16), więc definicja regionu nie jest jednorodna w LB. Najprościej liczyć R3 bez LB 4–16 (L = 26,83 nb⁻¹ **[R]**: 1,90 + 24,92) albo z wagą per LB.
- Region „bez wymogu trygera” nadal nie ma jednej świetlności (rozpiętość 63×).

## 6. Przy okazji

- ⟨μ⟩ oficjalne: 0,00469 (427929) i 0,0452 (428770) = `average_mu` w drzewie. Wartości 0,0085 i 0,054 używane w nocie to średnia `actual_mu` po zdarzeniach (per paczka, ważona zdarzeniami).
- 427929 ma w oficjalnym GRL AFP LB 88–485 i 2085 nb⁻¹ — to część przy dużym μ; nasz zakres 119–324 to 0,023 % tej świetlności.

## 7. Zablokowane / niezrobione

- Lista wszystkich tagów w COOL (`listTags` na `/TRIGGER/OFLLUMI/OflPrefLumi`) zawiesiła się na > 8 min — przerwałem własny proces; tagi sprawdzone próbami 001–007.
- Strona LuminosityForPhysics (twiki, SSO) i atlas-lumicalc (SSO) — pominięte; aktualność tagu 003 dla runów niskiego μ wynika tylko z plików cvmfs (2023) i z tego, że 004 daje to samo.
- Nic nie wymagało certyfikatu grid.

## 8. Pliki

- `lumi.json`; `tabele/lumi_final_perLB_{427929,428770}.csv`, `tabele/lumi_perLB_*.csv`, `tabele/cvmfs_*__LumiMetaData.csv`, `tabele/a3.json`.
- `skrypty/`: lxplus `l0_help.sh`, `l1_tags.sh`, `l2_ilumicalc.sh`, `l3_cool_counters.py` + `l3_run.sh`, `l4_tags_scan.sh`; lokalnie `a1_cvmfs_ilumicalc.py`, `a2_lumi_ourLB.py`, `a3_lumi_perLB.py`, `a4_final.py` (kolejność a1 → a3 → a4; `.venv_h5/bin/python -I`).
- `surowe/`: `cvmfs_20231005/`, `cvmfs_20231024/` (GRL, notes, ilumicalc), `lxplus/` (kopia `/tmp/klyczak_lumi`: logi iLumiCalc, `cool/l1counters_<run>.csv`, `tagscan/`), `GRL_MinBiasAFP_lowmu2022.xml`, logi a1–a4.
- lxplus: `/tmp/klyczak_lumi` na lxplus9103 (ok. 1 MB) zostawiony.
