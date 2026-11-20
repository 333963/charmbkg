# Ntuple runu 428770 pochodzą ze strumienia `physics_Main`; łańcuchy `mb_*` (AFP + tor, MBTS) piszą do innego strumienia, którego nie przetworzono

Przypisanie łańcuchów `mb_*` do strumienia `MinBias` potwierdzono w źródle menu trygera, nie w bazie trygera dla tego runu.

Źródła: analiza kierunku z 20.11.2026 (część wspólna §3.1 i przegląd kodu charmpp §1.3 — oba poza repo);
[`results/tagreg/v2/RAPORT.md`](results/tagreg/v2/RAPORT.md) §4.3.

**Liczba:** dataset `data22_13p6TeV.00428770.physics_Main.merge.AOD.r16487_p5949`, 193 pliki, 11 503 705 zdarzeń, CharmAnalysis zachowuje 7 989 496 (źródło: charmpp (gałąź minimum-bias): `docs/CHARMPP_v2.0_PRODUCTION.md` l. 60, przeczytane w lustrze 20.11; liczba zachowanych zdarzeń wg przeglądu kodu charmpp §1.3; błędy nie dotyczą). W drzewie `mb_sptrk_pt4_L1AFP_A_OR_C` występuje bez łańcucha z zbioru S w 82 zdarzeniach na 177 876 (0,05 %); drzewo ma ok. 24 % zdarzeń tego łańcucha (14,4 / 60), ok. 34 % wśród zdarzeń z kandydatem D (20,3 / 60) (`TAGREG_2026-11-18/v2/RAPORT.md` §4.3, l. 182–208).

## Kontekst
- Odpowiada na otwarte pytanie z TAGREG v2 pkt 7 i tabeli braków: „z jakiego strumienia / skimu pochodzi drzewo 428770”. Nazwa strumienia jest w dokumentacji produkcji v2.0 na GitLab (lustro `~/Projects/charm/gitlab_mirror/charmpp.git`).
- Run 427929 jest w produkcji brany z `data22_13p6TeV.00427929.physics_MinBias.merge.AOD.r15138_p5949` (583 pliki, 57 878 241 zdarzeń, status „pending”; `docs/CHARMPP_v2.0_PRODUCTION.md` l. 61, przeczytane w lustrze). Dla 428770 strumienia `physics_MinBias` w liście próbek nie ma.
- Wniosek (nie fakt): łańcuchy `HLT_mb_sptrk_*` i `HLT_*_L1MBTS_*` zapisują do innego strumienia, najpewniej `physics_MinBias`; zdarzenia trafiają do drzewa tylko wtedy, gdy przeszły też łańcuch strumienia `physics_Main`.
- Skutek dla regionu z tagiem: `mb_sptrk_pt4_L1AFP_A_OR_C` ma preskalę 1 (2 w LB 4–16), ułamek czasu żywego 0,9563, czyli ok. 28 nb⁻¹ wobec 0,469 nb⁻¹ regionu R1. W próbce nieobciążonej (R1) ten bit ma 32,4 % zdarzeń z D⁰ wobec 13 % dla `j20`; pełny strumień wybrałby 2,4 × więcej zdarzeń z D⁰ niż `j20_L1AFP` (`TAGREG v2` l. 203–208). R2 dał 278 ± 49 ± 48 tagowanych D⁰ i 185 ± 30 ± 24 D⁺ (`TAGGED_2026-11-19/RAPORT.md` §0.4), więc pełny strumień AFP + tor to rząd kilkuset tagowanych D na kanał (szacunek, nie pomiar).
- Naturalny mianownik w tym samym strumieniu: `mb_sptrk_pt4_L1MBTS_2`, preskala 50 (25 w LB 4–20), ok. 0,6 nb⁻¹; ten sam wymóg HLT, inne ziarno L1.

## Sprawdzenie w źródle menu (20.11, ~04:30)
- Plik `TriggerMenuMT/HLT/Menu/PhysicsP1_pp_lowMu_run3_v1.py` w wydaniach Athena 22.0.70, 22.0.73, 22.0.76, 22.0.79 i 22.0.85 na cvmfs (daty plików 12.06–16.08.2022; odczyt przez `bin/lx`, grep):
  - `HLT_mb_sptrk_pt4_L1AFP_A_OR_C` → `stream=['MinBias']` (l. 121; l. 126 w 22.0.85),
  - `HLT_mb_sptrk_pt4_L1MBTS_2` → `stream=['MinBias']` (l. 115; l. 120 w 22.0.85),
  - `HLT_noalg_L1J12` → `stream=['Main']` (22.0.85 l. 268),
  - `HLT_j20_L1AFP_A_OR_C` i `HLT_noalg_L1AFP_A_OR_C` → `stream=[PhysicsStream]` (22.0.85 l. 60, 298).
- Wniosek z sekcji „Kontekst” jest więc potwierdzony na poziomie definicji menu: łańcuch AFP + tor zapisuje do strumienia `MinBias`, a nasze ntuple runu 428770 powstały z `physics_Main`.
- Nie sprawdzono, którego wydania użyto online w runie 428770; przypisanie jest jednakowe we wszystkich pięciu sprawdzonych.

## Zastrzeżenia
- Przypisanie łańcuchów do strumieni sprawdzono 20.11 (~04:30) w źródle menu trygera, nie w bazie trygera dla tego runu (zob. „Sprawdzenie w źródle menu”). Istnienia AOD `physics_MinBias` runu 428770 nadal nie sprawdzono: wymaga proxy gridowego (`~/Documents/charm/START_grid.command`, wynik w `runs/grid/rucio_streams.txt`).
- Szacunek „kilkuset” skaluje uzysk R2 czynnikiem 2,4 liczonym na zdarzeniach z kandydatem D⁰, nie na sygnale; nie uwzględnia innego składu tła.
- TAGREG v2 §4.3 dopuszcza małą niekompletność (3–4 %) także dla zdarzeń dżetowych.

## Powiązane (w repo)
- Region tagowany — skład trygerów, L_eff, oczekiwany uzysk w 427929, wydajność ξ_cal: [`results/tagreg/RAPORT.md`](results/tagreg/RAPORT.md), [`results/tagreg/v2/RAPORT.md`](results/tagreg/v2/RAPORT.md).
- Skład trygerowy runu 428770: [`results/b_trig/`](results/b_trig/).
- Odczyt z tagiem (R2: 278 ± 49 ± 48 D⁰, 185 ± 30 ± 24 D⁺): [`results/tagged/RAPORT.md`](results/tagged/RAPORT.md) §0.4.
