# TAGREG v2 — pełne bity trygera w regionie tagowanym AFP (wynik etapu `tagtrig`), 18.11.2026

Projekt B, wkład do pomiaru pojedynczego D z tagiem protonu AFP. Wejście: `../runs/lxplus/` (tylko odczyt). Kryteria zapisane przed liczeniem: `KRYTERIA.md`.
Praca lokalna.

**Granice dotrzymane:** w regionie tagowanym tylko liczby zdarzeń, bity, preskale i liczniki kandydatów (`n_D0`, `n_Dp`, `n_Dst` > 0). Żadnych mas, fitów, uzysków,
przekrojów, stosunków σ_tag/σ_incl. „Zdarzenie z D⁰” = zdarzenie z ≥ 1 kandydatem slim, nie sygnał.

## 0. Wynik w skrócie

1. **Kontrola przeszła.** 97 / 458 plików ROOT, 0 duplikatów; T0/T1/T2, T2∧D⁰, podział A/C/AC i tabela 8 klas z 11 bitów identyczne co do zdarzenia z `../RAPORT.md`.
   Uwaga formatowa: „106 łańcuchów” w drzewie to **75 łańcuchów HLT + 31 pozycji L1 po preskali**; osobno 31 bitów L1 przed preskalą.
2. **Klasa `none` w 428770 jest rozpoznana w 99,9 %.** Zdarzenia (459 662): 69,7 % `HLT_noalg_L1AFP_A_AND_C`, 19,8 % `HLT_j20_L1AFP_A_OR_C`, 2,4 % `j20f_L1AFP`,
   7,0 % `noalg_L1AFP_{FSA,NSA,FSC,NSC}_BGRP12`. Zdarzenia z D⁰ (16 331): **80,1 % `HLT_j20_L1AFP_A_OR_C`**, 7,1 % `noalg_L1AFP_A_AND_C`, 6,8 % `j20f_L1AFP`.
   Bez żadnego bitu HLT zostaje 571 zdarzeń T2 (0,05 %; 124 z D⁰).
3. **Niespójność preskal z poprzedniego raportu jest wyjaśniona: to nie preskala, tylko skład drzewa.** Drzewo 428770 zawiera zdarzenie tylko wtedy, gdy przeszło
   ono łańcuch z określonego zbioru S (dżetowe, noalg AFP, EM, MU, J — nie ma w nim `mb_sptrk_*` ani łańcuchów MBTS). Zdarzenia `sptrk`/MBTS są w drzewie tylko
   „przy okazji”, więc nie są próbką niezależną od `noalg_L1AFP`. Po dodaniu warunku „przeszedł też inny łańcuch z S” wszystkie odniesienia dają 50–59, zgodnie z 60,
   a zależność od krotności znika. **Preskala `noalg_L1AFP_A_OR_C` = L1 1 (2 w LB 4–16) × HLT 60 — potwierdzona trzema drogami: drzewo, COOL, bity.**
4. **Wniosek dla pomiaru: w tym drzewie `HLT_mb_sptrk_pt4_L1AFP_A_OR_C` nie może być regionem trygera.** Drzewo ma tylko ~24 % zdarzeń tego łańcucha
   (dla zdarzeń z kandydatem D ~34 %), i to tych, które przeszły też dżet. Rekomendacja „proton + tor, 28 nb⁻¹” z `../RAPORT.md` §2.5 jest **nieaktualna**.
   Zdarzenia z D⁰ przypisane tam do „tylko sptrk_pt4” (52 %) są w rzeczywistości zdarzeniami `j20_L1AFP_A_OR_C`.
5. **Regiony możliwe w tym drzewie (428770, T2):** R1 `noalg_L1AFP_A_OR_C` (preskala 60/120, F = 0,0159; 48,8 % zdarzeń, 10,9 % zdarzeń z D⁰);
   R2 `j20_L1AFP_A_OR_C ∧ ¬R1` (1/2, F = 0,956; 17,3 %, 77,0 %); R3 `noalg_L1J12 ∧ ¬R1 ∧ ¬R2` (1, F = 1; 0,6 %, 4,5 %). Razem 66,8 % zdarzeń i 92,4 % zdarzeń z D⁰.
   R2 jest silnie obciążony: w próbce nieobciążonej (R1) bit `j20` ma tylko 13 % zdarzeń z D⁰, a bit `sptrk_pt4` 32 %.
6. **427929:** 99,99 % zdarzeń T2 ma `HLT_noalg_mb_L1MBTS_1`, wszystkie preskale 1 (drzewo = COOL). Jeden region MBTS wystarcza.
7. **Nadal brakuje** świetlności per LB (węzeł COOL jej nie liczył — brak tagu świetlności i GRL xml) oraz informacji, z jakiego strumienia / skimu pochodzi drzewo 428770.

## 1. Kontrola (`a1_control_none_fam.py`, `out/a1.log`)

| | 428770 | 427929 | zgodność z `../RAPORT.md` |
|---|---|---|---|
| pliki ROOT / duplikaty numeru zdarzenia | 97 / 0 | 458 / 0 | – |
| zdarzenia drzewa / w GRL | 7 989 496 / 7 989 496 | 28 962 571 / 23 874 370 | tak |
| T0 / T1 / T2 | 1 338 608 / 1 233 602 / 1 134 024 | 707 688 / 703 350 / 701 183 | identyczne |
| T2 z D⁰ / D⁺ / D* / dowolnym D | 76 030 / 75 673 / 14 511 / 114 120 | 10 438 / 11 533 / 1 479 / 18 872 | D⁰ i „D” identyczne |
| T2: tylko A / tylko C / A∧C | 570 220 / 321 855 / 241 949 | 354 983 / 338 503 / 7 697 | identyczne |
| 8 klas z 11 bitów (zdarzenia, z D, z D⁰) | odtworzone co do zdarzenia | odtworzone | identyczne |
| lista łańcuchów w częściach / brakujące preskale | identyczna / brak | identyczna / brak | – |
| preskala stała w LB (`ps_nvals` = 1) | wszystkie 136 LB × 106 | wszystkie 189 LB × 106 | – |
| bit L1 po preskali ⊂ bit przed preskalą | 0 naruszeń | 0 naruszeń | – |

Klasę `none` poprzedniego raportu odtwarzam dokładnie: 459 662 zdarzeń, 16 331 z D⁰, 231 333 dwustronnych.

## 2. Klasa `none` (żaden z 11 bitów) w 428770, region T2 (`out/none11_ranking_428770.csv`)

„Łącznie” = bit ustawiony; „tylko” = jedyny bit HLT zdarzenia; „zachłannie” = kolejno łańcuch z największą liczbą jeszcze nieprzypisanych zdarzeń.

| łańcuch | zdarzenia łącznie | tylko | zachłannie [%] | z D⁰ łącznie | tylko | zachłannie [%] | z D⁺ łącznie | z D* łącznie |
|---|---|---|---|---|---|---|---|---|
| HLT_noalg_L1AFP_A_AND_C | 320 307 | 293 374 | 69,68 | 1 184 | 13 | 7,05 | 1 325 | 138 |
| HLT_j20_L1AFP_A_OR_C | 91 233 | 78 245 | 19,76 | 13 086 | 10 588 | **80,13** | 17 010 | 2 600 |
| HLT_mb_sptrk_pt2_L1AFP_A_AND_C | 26 132 | 0 | – | 1 171 | 0 | – | 1 204 | 136 |
| HLT_j20f_L1AFP_A_OR_C | 13 872 | 10 255 | 2,44 | 1 563 | 961 | 6,82 | 2 221 | 306 |
| HLT_noalg_L1AFP_FSA / NSA / FSC / NSC_BGRP12 | 9 940 / 9 278 / 7 304 / 6 044 | 9 646 / 9 008 / 7 089 / 5 837 | 2,13 / 1,99 / 1,56 / 1,28 | 89 / 102 / 91 / 69 | 81 / 95 / 85 / 58 | ≤ 0,6 każdy | 128 / 127 / 103 / 80 | 16 / 16 / 10 / 13 |
| HLT_mb_sptrk_pt2_L1AFP_A_OR_C | 9 753 | 8 | – | 1 512 | 4 | – | 1 946 | 283 |
| HLT_2j10_20detaAA_L1AFP_A_OR_C | 4 379 | 259 | – | 949 | 56 | – | 1 156 | 205 |
| HLT_j10f_L1MBTS_2 | 1 475 | 1 175 | 0,26 | 126 | 96 | 0,64 | 187 | 18 |
| HLT_noalg_L1EM3 / EM7 / MU3V / MU5VF_AFP | 879 / 438 / 835 / 786 | 750 / 233 / 748 / 705 | – | 130 / 100 / 89 / 38 | 110 / 69 / 77 / 33 | 0,7 / – / 0,5 / – | 156 / 108 / 87 / 44 | 29 / 12 / 8 / 8 |
| **bez żadnego bitu HLT** | **571** (0,12 % klasy) | | | **124** (0,76 %) | | | 155 | 30 |

- Liczba bitów HLT w zdarzeniu klasy `none`: 0 → 571; 1 → 417 693; 2 → 37 689; ≥ 3 → 3 709. Klasa jest więc w 91 % jednołańcuchowa.
- Dla D⁺ i D* obraz jak dla D⁰: `j20_L1AFP_A_OR_C` 80,7 % i 83,4 %, `j20f` 7,7 % i 6,9 %, `noalg_L1AFP_A_AND_C` 6,1 % i 4,1 %.
- Zdarzenia dwustronne (231 tys. w klasie) to `noalg_L1AFP_A_AND_C` (preskala 1) — stąd ich przewaga i brak kandydatów.
- Reszta bez bitu HLT (571): 563 mają bity L1 przed preskalą AFP, EM3 i EM7, 159 ma L1_AFP po preskali, 391 nie ma żadnego bitu L1 po preskali; skupione w LB 9–14
  (32–42 na LB). To zdarzenia zapisane przez łańcuchy spoza 75 przechowywanych (menu ma 380) — nie do rozpoznania lokalnie; wpływ pomijalny.

## 3. Skład regionu tagowanego wg rodzin (`out/rodziny_T2_<run>.csv`, `out/lancuchy_<run>.csv`)

Hierarchia rozłączna (zapisana z góry): zdarzenie trafia do pierwszej rodziny, której bit ma. Porządek = rosnąca siła warunku centralnego przy protonie w L1,
potem trygery bez protonu. Rodzina mówi więc „najmniej obciążony tryger, którym to zdarzenie da się opisać”.
1 `AFPnoalg` = `noalg_L1AFP_A_OR_C`; 2 `AFPnoalg-inne` = pozostałe `noalg_L1AFP_*` (A_AND_C, BGRP12, EMPTY…); 3 `AFP+tor` = `mb_sptrk_pt{2,4,6,8}_L1AFP_*`;
4 `AFP+dżet` = `j20`, `j20f`, `2j10_20detaAA` na L1AFP, `*_L1AFP_A_OR_C_J12`, `noalg_L1MU5VF_AFP`; 5 `dżet` = `noalg_L1J*`, `j20_L1J12`, dżetowe na MBTS_2 i RD0;
6 `MBTS`; 7 `pozostałe` (RD0, EM, MU, TAU); 8 bez bitu HLT.

### 3.1 428770, T2 — rozłącznie [% zbioru] (w nawiasie: łącznie, czyli zdarzenia z bitem rodziny)

| rodzina | zdarzenia 1 134 024 | z D 114 120 | z D⁰ 76 030 | z D⁺ 75 673 | z D* 14 511 |
|---|---|---|---|---|---|
| AFPnoalg | 553 760 = 48,83 (48,83) | 12,88 (12,88) | 10,87 (10,87) | 11,89 | 8,59 |
| AFPnoalg-inne | 355 990 = 31,39 (31,89) | 3,41 (3,46) | 2,92 (2,97) | 2,82 | 1,87 |
| AFP+tor | 111 555 = 9,84 (13,93) | 56,37 (63,20) | 62,33 (69,03) | 55,82 | 63,73 |
| AFP+dżet | 100 774 = 8,89 (19,09) | 22,19 (80,07) | 18,21 (82,13) | 23,86 | 19,09 |
| dżet | 9 634 = 0,85 (2,47) | 4,62 (14,36) | 5,16 (16,23) | 5,05 | 6,25 |
| MBTS | 6 = 0,00 (0,49) | 0,00 (2,09) | 0,00 (2,33) | 0,00 | 0,00 |
| pozostałe | 1 734 = 0,15 (0,38) | 0,35 (1,55) | 0,34 (1,68) | 0,36 | 0,27 |
| bez bitu HLT | 571 = 0,05 | 0,18 | 0,16 | 0,20 | 0,21 |

- W T1 (przed cięciem ξ_cal) udziały różnią się o ≤ 2,5 pp (AFPnoalg 46,4 %, AFP+dżet 11,3 %; `out/rodziny_T1_428770.csv`).
- **Hierarchia wprowadza tu w błąd w jednym punkcie:** rodzina `AFP+tor` (62 % zdarzeń z D⁰) nie jest samodzielnym źródłem zdarzeń. Z 52 487 zdarzeń z D⁰ z bitem
  `AFP+tor` 48 314 ma też bit `AFP+dżet`; łańcuchy sptrk nie występują w drzewie same (§4.3). Faktycznym źródłem jest `AFP+dżet` (łącznie 82 % zdarzeń z D⁰).
- Łańcuchy wiodące (T2, zdarzenia / z D⁰): `noalg_L1AFP_A_OR_C` 553 760 / 8 267; `noalg_L1AFP_A_AND_C` 328 152 / 1 753; `j20_L1AFP_A_OR_C` 199 976 / 59 641;
  `mb_sptrk_pt4_L1AFP_A_OR_C` 116 296 / 49 274; `mb_sptrk_pt6_L1AFP` 54 736 / 23 955; `j20f_L1AFP` 22 294 / 5 093; `noalg_L1J12` 19 779 / 10 273
  (`j20_L1J12` ma ten sam zbiór: 19 772); `2j10_20detaAA_L1AFP` 18 739 / 8 711; `j20_L1MBTS_2` 6 463 / 1 979.

### 3.2 427929, T2

| rodzina | zdarzenia 701 183: rozłącznie [%] (łącznie [%]) | z D⁰ 10 438: rozłącznie (łącznie) |
|---|---|---|
| AFPnoalg | 476 064 = 67,89 (67,89) | 68,39 (68,39) |
| AFP+tor / AFP+dżet | 0 (10,44 / 0,41) | 0 (68,03 / 9,06) |
| dżet | 1 699 = 0,24 (0,73) | 4,51 (14,26) |
| MBTS | 223 409 = 31,86 (**99,99**) | 27,09 (**99,99**) |
| AFPnoalg-inne / pozostałe / bez bitu | 5 / 4 / 2 | 0 / 0 / 0 |

Tu właściwy opis to „łącznie”: MBTS ma 99,99 % zdarzeń (701 080 `noalg_mb_L1MBTS_1`; 103 bez), a 67,9 % `AFPnoalg` to wydajność L1 AFP, nie osobny strumień.

### 3.3 Zależność od LB (8 bloków o równej liczbie T1; `out/rodziny_LB_<run>.csv`)

| run | LB | zdarzenia: AFPnoalg / AFPnoalg-inne / AFP+tor / AFP+dżet / dżet [%] | z D⁰: AFPnoalg / AFP+tor / AFP+dżet / dżet [%] |
|---|---|---|---|
| 428770 | 4–24 | 43,1 / 37,9 / 8,8 / 8,3 / 1,5 | 9,7 / 57,2 / 19,5 / 8,5 |
| 428770 | 25–323 (7 bloków) | 49,4–50,0 / 30,0–30,9 / 9,9–10,0 / 8,8–9,1 / 0,7–0,8 | 10,8–11,4 / 62,6–63,8 / 17,2–18,5 / 4,1–5,4 |
| 427929 | 119–324 (8 bloków) | AFPnoalg 67,5–68,3; MBTS 31,5–32,3 | AFPnoalg 65,9–70,2 |

Skład stały; jedyna zmiana to blok LB 4–24 w 428770 (w LB 4–16 preskala L1_AFP = 2). W bloku LB 4–24 jest też 0,3 % zdarzeń bez bitu HLT (skupione w LB 9–14).

## 4. Preskale (`a2_prescale_stream.py`, `out/a2.log`)

### 4.1 Status COOL i porównanie z drzewem

Węzeł COOL zakończył się kodem 0, 0 błędów; zwrócił klucze i preskale L1 i HLT osobno dla 11 łańcuchów HLT i 5 pozycji L1 (131 / 154 wierszy). Świetlności nie liczył.
Zakresy kluczy pokrywają wszystkie LB GRL (428770: HLT 1586, 1587; L1 1631, 1630, 1639; 427929: HLT 1390, L1 1384); pozycje poza GRL są wyłączone.

| łańcuch (428770) | drzewo `trigPS` | COOL L1 | COOL HLT | L1 × HLT = drzewo |
|---|---|---|---|---|
| HLT_noalg_L1AFP_A_OR_C | 60 (LB 17–323), 120 (LB 4–16) | 1; 2 w LB 4–16 | 60 | 136 / 136 LB |
| HLT_mb_sptrk_pt4_L1AFP_A_OR_C | 1, 2 | 1; 2 | 1 | 136 / 136 |
| HLT_noalg_L1J12 | 1 | 1 | 1 | 136 / 136 |
| HLT_j20_L1MBTS_2, HLT_mb_sptrk_pt4_L1MBTS_2 | 50 (LB 21–323), 25 (LB 4–20) | 50; 25 | 1 | 136 / 136 |
| HLT_mb_sptrk_pt2_L1MBTS_2, HLT_noalg_mb_L1MBTS_2 | 1250, 625 | 50; 25 | 25 | 136 / 136 |
| HLT_noalg_mb_L1MBTS_1 / HLT_mb_mbts_L1MBTS_1 | 5000,66 / 50 006,6 | 5000,66 | 1 / 10 | 136 / 136 |
| HLT_mb_sptrk_L1MBTS_1, HLT_noalg_L1RD0_FILLED | 0 | – | wyłączone | zgodne |

Największa różnica względna 2·10⁻⁸. 427929: wszystkie 1 × 1 w 189 LB GRL (RD0 wyłączony). **`trigPS` w drzewie = preskala L1 × preskala HLT.**
Ponieważ drzewo ma też `trigPS` pozycji L1, czynnik HLT pozostałych 64 łańcuchów wynika z samego drzewa (`out/ps_drzewo_rozklad_<run>.csv`), np. 428770:
`j20_L1AFP_A_OR_C`, `j20f`, `2j10`, `sptrk_pt6/pt8_L1AFP` 1 (całkowita 1/2); `sptrk_pt2_L1AFP_A_OR_C` HLT 10; `noalg_L1AFP_A_AND_C` 1 × 1;
`noalg_L1AFP_*_BGRP12` L1 50 × HLT 50; `noalg_L1EM7` 20; `EM3` 250; `MU3V` 15. Dla tych łańcuchów COOL nie był pytany — status „tylko drzewo”.

### 4.2 Efektywne preskale z bitów (kontrola niezależna od drzewa i od COOL; zdarzenia T0, 428770)

L1: p = N(bit przed preskalą ∧ odniesienie) / N(bit po preskali ∧ odniesienie), odniesienie o innym L1.

| pozycja L1 | odniesienie | blok LB | N przed / po | p_eff | drzewo / COOL |
|---|---|---|---|---|---|
| L1_AFP_A_OR_C | noalg_L1J12 | 4–16 | 1 791 / 910 | **1,97 ± 0,05** | 2 |
| L1_AFP_A_OR_C | noalg_L1J12 | 17–323 | 19 092 / 19 092 | **1,000** | 1 |
| L1_MBTS_2 | noalg_L1AFP (OR, AND) | 4–20 | 67 622 / 2 730 | 24,8 ± 0,5 | 25 |
| L1_MBTS_2 | noalg_L1AFP (OR, AND) | 21–323 | 700 565 / 14 112 | 49,6 ± 0,4 | 50 |
| L1_MBTS_1 | noalg AFP ∨ J12 | wszystkie | 895 526 / 194 | 4 620 ± 330 | 5000,66 |
| L1_J12 / L1_AFP_A_AND_C | noalg AFP / J12 | wszystkie | 459 / 459; 96 / 96 | 1,000 | 1 |
| L1_AFP_A / L1_AFP_FSA_BGRP12 | noalg_L1J12 | wszystkie | 14 174 / 273; 15 389 / 334 | 52 ± 3 / 46 ± 3 | 50 |
| L1_EM7 / L1_EM3 / L1_MU3V / L1_RD0_FILLED | noalg AFP (∨ J12) | wszystkie | – | 24 ± 5 / 330 ± 70 / 16 ± 2 / 995 ± 32 | 20 / 250 / 15 / 1000 |

HLT, łańcuch `noalg_L1AFP_A_OR_C`: p_HLT = N(odniesienie ∧ L1_AFP po preskali) / N(… ∧ noalg), LB ≥ 17 (L1 = 1, więc to także preskala całkowita):

| odniesienie | w zbiorze S | n odniesienia | p_eff | p_eff z warunkiem „inny łańcuch z S” | udział odniesienia z innym łańcuchem S |
|---|---|---|---|---|---|
| noalg_L1AFP_A_AND_C | tak | 331 570 | **60,6 ± 0,8** | – | – |
| j20_L1AFP_A_OR_C | tak | 300 929 | **58,0 ± 0,8** | 56,6 ± 1,6 | 0,22 |
| j20f_L1AFP_A_OR_C | tak | 33 510 | 61,5 ± 2,6 | 59,6 ± 3,8 | 0,42 |
| 2j10_20detaAA_L1AFP_A_OR_C | tak | 40 073 | 57,0 ± 2,1 | 57,4 ± 2,2 | 0,92 |
| noalg_L1J12 | tak | 19 092 | 55,5 ± 3,0 | – | – |
| j20_L1MBTS_2 | tak | 6 031 | 59,1 ± 5,8 | – | – |
| noalg_L1AFP_FSA_BGRP12 | tak | 9 234 | 60,4 ± 4,8 | – | – |
| mb_sptrk_pt4_L1AFP_A_OR_C | **nie** | 169 721 | **14,4 ± 0,1** | **57,6 ± 1,1** | 0,947 |
| mb_sptrk_pt6_L1AFP_A_OR_C | nie | 78 497 | 28,6 ± 0,5 | 58,1 ± 1,6 | 0,981 |
| mb_sptrk_pt8_L1AFP_A_OR_C | nie | 37 454 | 37,1 ± 1,2 | 57,9 ± 2,3 | 0,990 |
| mb_sptrk_pt2_L1AFP_A_OR_C | nie | 39 840 | 4,11 ± 0,04 | 58,6 ± 2,5 | 0,769 |
| mb_sptrk_pt4_L1MBTS_2 | nie | 3 615 | 15,0 ± 0,9 | 52,9 ± 6,5 | 0,951 |
| noalg_mb_L1MBTS_2 | nie | 920 | 2,41 ± 0,09 | 50 ± 15 | 0,597 |

LB 4–16 (L1 = 2): p_HLT wobec `j20_L1AFP` 61,6 ± 4,0, wobec `A_AND_C` 56,7 ± 3,4 — czynnik HLT 60 niezależny od bloku, czynnik 2 siedzi w L1.
Inne łańcuchy noalg wobec odniesień ze zbioru S: `noalg_L1EM7` 20,4 ± 0,6 (20), `EM3` 249 ± 11 (250), `MU3V` 16,0 ± 1,5 (15), `noalg_L1AFP_FSA/NSA_BGRP12` HLT 52 ± 3 / 48 ± 3 (50),
`noalg_mb_L1MBTS_2` całkowita 1 184 ± 46 (średnia z 625/1250), HLT 25,8 ± 1,0 (25), `noalg_mb_L1MBTS_1` 5 020 ± 390 (5000,66).

### 4.3 Wyjaśnienie niespójności z `../RAPORT.md` §2.4 (hipoteza H-strumień, postawiona przed liczeniem — potwierdzona)

- **(a) Zbiór S.** Jako jedyny bit HLT zdarzenia występują tylko: `noalg_L1AFP_A_OR_C` (562 076 z 593 058 zdarzeń T0), `noalg_L1AFP_A_AND_C` (315 936), `j20_L1AFP_A_OR_C`
  (123 789 z 315 462), `j20f_L1AFP` (14 853), `2j10_20detaAA_L1AFP` (1 701), cztery `noalg_L1AFP_*_BGRP12`, `noalg_L1EM3/EM7/MU3V/MU5VF_AFP`, `j20/j10f_L1MBTS_2`,
  `j20/j10f_L1RD0_FILLED`, warianty EMPTY/FIRSTEMPTY/UNPAIRED; po ich wyłączeniu jeszcze pary `noalg_L1J12` + `j20_L1J12`, `*_L1AFP_A_OR_C_J12`, `noalg_L1J75/J85/J100/TAU60/EM22VHI`.
  Łańcuchy `mb_sptrk_*` i `mb_mbts_*`/`noalg_mb_*` na MBTS **nigdy lub prawie nigdy** nie występują bez łańcucha z S: `sptrk_pt4_L1AFP_A_OR_C` 82 z 177 876 (0,05 %),
  `pt6` 48 z 82 181, `pt2` 24 z 41 828, `sptrk_pt4_L1MBTS_2` 0 z 6 018, `noalg_mb_L1MBTS_2` 0 z 1 203 (`out/strumien_S_428770.csv`).
  Poziom 0,05 % odpowiada częstości zdarzeń bez żadnego bitu HLT (0,06 %), czyli łańcuchom nieprzechowywanym.
- **(b)** 99,93 % zdarzeń T0 ma ≥ 1 bit z S (874 bez; 779 z nich bez żadnego bitu HLT).
- **(c)** Z warunkiem „przeszedł też łańcuch z S inny niż `noalg_L1AFP`” p_eff wobec odniesień sptrk/MBTS = 57,6 ± 1,1; 58,1 ± 1,6; 57,9 ± 2,3; 58,6 ± 2,5; 52,9 ± 6,5;
  50 ± 15 — wszystkie w 10 % od 60 (kryterium spełnione). Bez warunku: 14,4; 28,6; 37,1; 4,1; 15,0; 2,4, czyli dokładnie 60 × udział zdarzeń odniesienia, które i tak
  byłyby w drzewie (np. 60 × 0,24 = 14,4).
- **(d)** Zależność od krotności znika: odniesienie `sptrk_pt4_L1AFP`, nSelBL < 20 / 20–40 / ≥ 40: bez warunku 10,0 / 18,9 / 27,9; z warunkiem 58,8 ± 1,7 / 57,1 ± 1,7 / 56,2 ± 2,4.
  Z kandydatem D / bez: 20,3 / 9,4 → 57,8 ± 1,3 / 57,1 ± 1,8. Rosła z krotnością, bo z krotnością rośnie szansa, że zdarzenie sptrk przeszło też `j20`.
- **Werdykt: PASS (a)–(d).** Niespójność nie dotyczy preskali. Drzewo 428770 jest wyborem zdarzeń według zbioru S (strumień albo skim — lokalnie nierozróżnialne;
  pamięć projektu: „428770 = jetowy physics_Main”; 56 % zdarzeń drzewa ma `noalg_L1J12`). Liczba „≈ 23” z `AFPTRIG` i 14,6 / 16,8 / 2,5 z poprzedniego raportu to artefakty
  tego wyboru.
- Zastrzeżenie: odniesienia `j20_L1AFP` i sptrk z warunkiem dają 57–58, czyli 3–4 % poniżej 60 (2,5σ dla `j20`), a `A_AND_C` daje 60,6 ± 0,8. Mieści się w kryterium 10 %,
  ale nie wykluczam małej niekompletności także dla zdarzeń dżetowych (łańcuchy nieprzechowywane). Rozstrzygnie to licznik zdarzeń strumienia z AMI/COOL.

### 4.4 Skutek dla łańcuchów spoza S

W próbce nieobciążonej względem HLT (`noalg_L1AFP_A_OR_C` ∧ T2; `out/a4_noalg_fractions.log`) bit ma:

| bit | wszystkie zdarzenia | z D | z D⁰ | z D⁺ | z D* |
|---|---|---|---|---|---|
| mb_sptrk_pt4_L1AFP_A_OR_C | 1,66 % | 24,4 % | **32,4 %** | 21,9 % | 32,7 % |
| j20_L1AFP_A_OR_C | 0,63 % | 10,6 % | **13,3 %** | 11,7 % | 16,9 % |
| mb_sptrk_pt6_L1AFP_A_OR_C | 0,37 % | 6,1 % | 8,0 % | 6,3 % | 10,7 % |
| noalg_L1J12 | 0,04 % | 1,2 % | 1,7 % | 1,4 % | 3,1 % |

`sptrk_pt4_L1AFP` wybrałby 2,4× więcej zdarzeń z D⁰ niż `j20_L1AFP`, ale drzewo ma z niego tylko część, która przeszła też łańcuch z S (24 % wszystkich, ~34 % z kandydatem D: 14,4/60 i 20,3/60).
Wagi 1/p ani L_eff łańcuchów spoza S nie opisują ich zdarzeń w tym drzewie (N·p zaniżone 1,6–25×).

## 5. Świetlność względna i regiony (`a3_lumi_regions.py`, `out/lancuchy_F_pokrycie_<run>.csv`, `out/regiony_<run>.csv`)

F = Σ_LB w_LB / p_LB; w_LB: (a) liczba zdarzeń `noalg_L1J12` w LB (wszystkie zdarzenia GRL), (b) ⟨μ⟩ przy stałej długości LB. Świetlności bezwzględnej nie ma.
„Znana preskala” = stała w LB i > 0 w całym GRL, łańcuch w S, zgodność z COOL (gdzie pytany) i z bitami.

### 5.1 Główne łańcuchy, 428770, T2

| łańcuch | preskala (L1 × HLT) | F (a) / (b) | w S | potwierdzenie | zdarzenia [%] | z D⁰ [%] | z D⁺ [%] | z D* [%] |
|---|---|---|---|---|---|---|---|---|
| HLT_noalg_L1AFP_A_OR_C | 60 (1 × 60); 120 w LB 4–16 | 0,01594 / 0,01586 | tak | drzewo, COOL, bity | 48,83 | 10,87 | 11,89 | 8,59 |
| HLT_j20_L1AFP_A_OR_C | 1; 2 w LB 4–16 | 0,9563 / 0,9518 | tak | drzewo, L1 z COOL i bitów | 17,63 | 78,44 | 77,18 | 80,77 |
| HLT_noalg_L1J12 (≡ HLT_j20_L1J12) | 1 | 1,000 | tak | drzewo, COOL, bity | 1,74 | 13,51 | 12,84 | 17,69 |
| HLT_2j10_20detaAA_L1AFP_A_OR_C | 1; 2 | 0,9563 | tak | drzewo | 1,65 | 11,46 | 11,02 | 13,86 |
| HLT_j20f_L1AFP_A_OR_C | 1; 2 | 0,9563 | tak | drzewo | 1,97 | 6,70 | 7,03 | 7,06 |
| HLT_noalg_L1AFP_A_AND_C | 1 | 1,000 | tak | drzewo, bity | 28,94 | 2,31 | 2,14 | 1,36 |
| HLT_j20_L1MBTS_2 | 50; 25 w LB 4–20 | 0,02229 / 0,02254 | tak | drzewo, COOL, bity (L1) | 0,57 | 2,60 | 2,53 | 2,43 |
| HLT_noalg_L1AFP_{FSA,NSA,FSC,NSC}_BGRP12 | 2500 (50 × 50) | 0,0004 | tak | drzewo, bity | 0,55–0,90 | 0,12–0,21 | – | – |
| HLT_mb_sptrk_pt4_L1AFP_A_OR_C | 1; 2 | 0,9563 | **nie** | preskala zgodna, próbka niekompletna | 10,26 | 64,81 | 56,38 | 65,07 |
| HLT_mb_sptrk_pt6 / pt8 / pt2_L1AFP_A_OR_C | 1; 1; 10 (×2 w LB 4–16) | 0,956; 0,956; 0,0956 | nie | jw. | 4,83 / 2,32 / 2,57 | 31,5 / 15,4 / 9,2 | – | – |

F zgadza się z `../RAPORT.md` §2.5 dla wspólnych łańcuchów (1,000; 0,9563; 0,01594; 0,02229). Różnica (a) − (b) ≤ 0,5 % (1 % dla MBTS_2).

### 5.2 Ile zdarzeń da się opisać jednym łańcuchem o znanej preskali (428770, T2)

- ≥ 1 bit z 30 łańcuchów o znanej preskali: 99,94 % zdarzeń, 99,80 % z D, 99,82 % z D⁰, 99,78 % z D⁺, 99,77 % z D*.
- Jednym łańcuchem: najwyżej 48,8 % zdarzeń (`noalg_L1AFP_A_OR_C`) albo 78,4 % zdarzeń z D⁰ (`j20_L1AFP_A_OR_C`).
- Trzema rozłącznymi: 66,8 % zdarzeń, **92,4 % zdarzeń z D⁰**, 92,0 % z D⁺, 93,8 % z D*; pięcioma: 96,4 % i 97,1 %.

### 5.3 Proponowane regiony rozłączne

**428770** (R1 ustalony z góry jako jedyny bez warunku centralnego; R2, R3 zachłannie wg zdarzeń z D⁰):

| region | definicja | preskala | F | zdarzenia T2 | udział | z D⁰ | udział | z D⁺ / D* [%] |
|---|---|---|---|---|---|---|---|---|
| R1 „proton” | HLT_noalg_L1AFP_A_OR_C | 60; 120 w LB 4–16 | 0,0159 | 553 760 | 48,83 % | 8 267 | 10,87 % | 11,89 / 8,59 |
| R2 „proton + dżet” | HLT_j20_L1AFP_A_OR_C ∧ ¬R1 | 1; 2 w LB 4–16 | 0,956 | 196 501 | 17,33 % | 58 542 | 77,00 % | 75,79 / 79,32 |
| R3 „dżet L1” | HLT_noalg_L1J12 ∧ ¬R1 ∧ ¬R2 | 1 | 1,000 | 6 635 | 0,59 % | 3 415 | 4,49 % | 4,33 / 5,84 |
| (R4) | HLT_j20f_L1AFP_A_OR_C ∧ ¬R1–R3 | 1; 2 | 0,956 | 13 853 | 1,22 % | 2 072 | 2,73 % | 3,13 / 2,61 |
| (R5) | HLT_noalg_L1AFP_A_AND_C ∧ ¬R1–R4 | 1 | 1,000 | 322 160 | 28,41 % | 1 550 | 2,04 % | 1,94 / 1,10 |
| poza | – | – | – | 41 115 | 3,63 % | 2 184 | 2,87 % | 2,93 / 2,54 |

- Udziały stałe w LB 25–323 (R1 49,4–50,0 % zdarzeń, 10,8–11,4 % z D⁰; R2 17,4–17,8 %, 77,4–78,6 %); w LB 4–24 R2 ma 70,2 % zdarzeń z D⁰, R3 8,6 %.
- Weto ¬R1 w R2 zmniejsza świetlność R2 o czynnik (1 − 1/60) = 0,983 (czysta preskala, niezależna). Weto ¬R2 w R3 jest cięciem fizycznym (brak L1_AFP albo brak dżetu HLT)
  i musi wejść do wydajności R3.
- R2 wymaga wydajności `j20` względem zdarzenia z D i protonem: tylko 13 % zdarzeń z D⁰ próbki nieobciążonej ma ten bit, więc poprawka jest rzędu 7× i zależy od pT D.
  R1 może ją zmierzyć na danych (ok. 1 100 zdarzeń z D⁰ z oboma bitami w T2), ale to już krok pomiarowy — nie wykonany.
- R3 nie ma protonu w trygerze: 5 % zdarzeń z D⁰, obciążenie dżetowe L1; sensowny raczej jako kontrola niż region pomiaru.
- R5 to zdarzenia dwustronne z preskalą 1 — dużo zdarzeń, prawie bez kandydatów; region kontrolny dla tagu podwójnego.
- **Poprzednia propozycja „proton + tor” (`sptrk_pt4_L1AFP`) odpada w tym drzewie.** Byłaby lepsza (32 % zamiast 13 % zdarzeń z D⁰, słabsze obciążenie),
  ale wymaga danych ze strumienia, do którego ten łańcuch pisze.

**427929:** jeden region `HLT_noalg_mb_L1MBTS_1` (preskala 1, F = 1): 99,99 % zdarzeń T2 i 99,99 % zdarzeń z D⁰ (10 437 z 10 438). Alternatywnie `noalg_mb_L1MBTS_2`
(98,2 % / 99,97 %). Bit `noalg_L1AFP_A_OR_C` (67,9 %) nie jest potrzebny.

## 6. Czego nadal brakuje i pytania

| brak | stan |
|---|---|
| Świetlność per LB (dostarczona, zapisana, czas żywy) dla 428770 i 427929 | węzeł `lumi` jej nie liczył: brak tagu świetlności 2022 low-μ i pliku GRL xml w materiałach; `iLumiCalc` nieuruchomiony |
| Świetlność strumienia, z którego jest drzewo 428770, i liczba zdarzeń tego strumienia (kontrola kompletności; 3–4 % z §4.3) | brak; potrzebna nazwa strumienia / datasetu |
| Preskale z COOL dla 64 pozostałych łańcuchów (w tym `j20_L1AFP_A_OR_C`) | COOL pytany tylko o 11 + 5; czynnik HLT znany z drzewa i z bitów |
| Decyzje HLT przed preskalą | nie ma; przy regionach rozłącznych niepotrzebne |
| Tożsamość 571 zdarzeń T2 bez bitu HLT | łańcuchy spoza 75 przechowywanych; pomijalne |
| 427929: kompletność kopii (458 z 583 plików) | bez zmian |

**Pytania do autora pomiaru:**
1. Z jakiego strumienia (albo po jakim skimie trygerowym) zrobione są ntuple 428770? W drzewie nie ma samodzielnych zdarzeń `mb_sptrk_*` ani MBTS, a 56 % zdarzeń ma
   `noalg_L1J12`. Czy istnieje dataset strumienia minimum-bias tego runu, w którym `HLT_mb_sptrk_pt4_L1AFP_A_OR_C` (preskala 1) jest kompletny?
2. Jaką świetlność i jaki tryger zakłada dla regionu tagowanego w 428770? Liczba 29,41 nb⁻¹ pasuje tylko do łańcuchów o preskali 1 (`j20_L1AFP_A_OR_C` w 96 %, `noalg_L1J12`);
   dla 49 % zdarzeń regionu (`noalg_L1AFP_A_OR_C`) właściwa jest 63× mniejsza. Prośba o tag świetlności i plik GRL xml do `iLumiCalc`.

## 7. Odstępstwa od `KRYTERIA.md` (po fakcie)

- Zbiór S: próg bezwzględny „> 20 zdarzeń z jedynym bitem” zastąpiony względnym „> 1 % zdarzeń łańcucha (i ≥ 5)”. Próg bezwzględny wpuszczał `mb_sptrk_pt4_L1AFP_A_OR_C`
  (28 z 177 876 = 0,016 %, poniżej częstości zdarzeń bez żadnego bitu HLT). Obie listy są w `out/a2.log`.
- Warunek testu (c): „łańcuch z S inny niż `noalg_L1AFP`” doprecyzowany na „…i inny niż samo odniesienie”.
- Blok preskal L1_MBTS_2 w tabeli §4.2 podany jako LB 4–20 (suma bloków 4–16 i 17–20 z `out/peff_L1_428770.csv`).
- Dodane po obejrzeniu wyników: udziały bitów w próbce `noalg_L1AFP` (§4.4), regiony R4–R5, tabela 427929 w ujęciu „łącznie”.
- 427929: wyznaczenie zbioru S tą metodą nie działa (łańcuchy MBTS zawsze przechodzą razem), więc kolumna „w S” w `out/lancuchy_F_pokrycie_427929.csv` nic nie znaczy;
  kompletność MBTS wynika stąd, że 23 850 567 z 23 874 370 zdarzeń GRL ma `noalg_mb_L1MBTS_1`.
- Kontrola klas 11-bitowych: pierwsze uruchomienie miało inną kolejność priorytetów (J12 przed MBTS) niż `../common.py`; poprawione, wynik identyczny.

## 8. Pliki (`~/Documents/charm/charm_work/TAGREG_2026-11-18/v2/`)

- `KRYTERIA.md` (przed liczeniem), `RAPORT.md`.
- Skrypty (kolejność; `~/Documents/charm/.venv_h5/bin/python`, razem < 2 min): `v2common.py`, `a1_control_none_fam.py`, `a2_prescale_stream.py`, `a3_lumi_regions.py`.
- `out/`: logi `a1.log`, `a2.log`, `a3.log`, `a4_noalg_fractions.log`; `a1.json`, `a2.json`, `a3.json`; tabele `none11_ranking_<run>.csv`, `rodziny_T1|T2_<run>.csv`,
  `rodziny_LB_<run>.csv`, `lancuchy_<run>.csv`, `strumien_S_<run>.csv`, `ps_drzewo_vs_cool_<run>.csv`, `ps_drzewo_rozklad_<run>.csv`, `peff_L1_428770.csv`,
  `peff_noalgAFP_428770.csv`, `peff_noalgAFP_przekroje_428770.csv`, `peff_inne_428770.csv`, `lancuchy_F_pokrycie_<run>.csv`, `regiony_<run>.csv`.
- Katalog ma 0,3 MB; żadnych plików tymczasowych ani cache.
