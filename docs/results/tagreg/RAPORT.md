# TAGREG — region tagowany AFP: skład trygerów i świetlność (A), przewidywany uzysk w 427929 (B), wydajność cięcia ξ_cal (C) — 18.11.2026

Projekt B, wkład do pomiaru pojedynczego D z tagiem protonu AFP. Zakres A, B, C w regionie tagowanym przyjęty 18.11.2026.
Praca lokalna, tylko odczyt `AFP_mix/runs/lxplus/*.parquet`, `*.meta.json`, `B_trig/runs/lxplus/*.npz`.
Kryteria zapisane przed liczeniem: `KRYTERIA.md`.

**Czego tu NIE ma (granice dotrzymane):** żadnego fitu masy / Δm ani uzysku sygnału D w regionie tagowanym (także żadnego „pik minus pasma”), żadnego przekroju,
żadnego R = σ_tag/σ_incl z danych. W regionie tagowanym są tylko liczby zdarzeń i kandydatów (całe okno slim albo pasma boczne), bity trygera, preskale,
rozkłady ξ, ξ_cal, frakcje tagu przypadkowego. Część B jest przewidywaniem (uzyski inkluzywne × założony ułamek tagu), nie pomiarem.

## 0. Wynik w skrócie

1. **Lokalnie jest wszystko oprócz świetlności per LB i 95 ze 106 łańcuchów.** Bity 11 łańcuchów (po preskali), 5 bitów L1 „before prescale”, preskale per LB
   (zapisane `trigPS`), LB, μ, GRL, npu, protony (n, ξ1, ξ2, tag), ξ_cal — dla wszystkich zdarzeń drzewa obu runów. Score klasyfikatora tylko dla 30 % zdarzeń (D⁰).
2. **428770, region sygnału (T2): 1 134 024 zdarzeń, 76 030 z kandydatem D⁰.** Zdarzenia: 48 % tylko `noalg_L1AFP`, 40,5 % klasa `none` (żaden z 11 zapisanych bitów;
   98 % ma bit L1_AFP), 8 % tylko `sptrk_pt4_L1AFP`. Zdarzenia z D⁰: 52 % tylko `sptrk_pt4_L1AFP`, 21,5 % `none`, 8 % AFP ∧ J12, 7 % tylko noalg, 5 % tylko J12.
   Skład stały w LB (poza LB 4–24, gdzie preskale były 2× większe). „96 % tagów ze strumienia AFP” nie potwierdza się w tej postaci: z zapisanymi bitami AFP jest
   58 % zdarzeń i 72 % zdarzeń z D⁰; reszta to głównie `none`.
3. **427929, T2: 701 183 zdarzeń, 10 438 z D⁰; 100 % MBTS** (68 % ma też bit AFP — to wydajność L1 AFP, nie osobny strumień); wszystkie łańcuchy bez preskali w GRL.
4. **Świetlność efektywna 428770** (ułamek F = L_eff/L z preskal per LB × 29,41 nb⁻¹ z głównego łańcucha; samej liczby 29,41 nie da się lokalnie zweryfikować):
   J12 29,4; `sptrk_pt4_L1AFP` 28,1; `noalg_L1AFP` 0,47; `j20`/`sptrk_pt4_L1MBTS_2` 0,66; `noalg_mb_L1MBTS_2` 0,026; `noalg_mb_L1MBTS_1` 0,0059 nb⁻¹.
   Jedna liczba 29,41 nb⁻¹ jest poprawna tylko dla zdarzeń J12 i (w 96 %) `sptrk_pt4_L1AFP`; dla 48 % zdarzeń regionu (tylko noalg) właściwa jest 63× mniejsza.
   Z wagami preskali region tagowany jest w > 99 % „noalg-podobny”. **Brakuje:** L per LB (oba runy), osobnych preskal L1/HLT, tożsamości łańcuchów klasy `none`.
5. **Oczekiwany uzysk tagowanych D w 427929 (nasza kopia runu, selekcja H WP80, pT 3–30 GeV): D⁰ ≈ 35 (18–71), D⁺ ≈ 18 (9–36), D* ≈ 2 (1–4) prawdziwych tagów;
   przypadkowych 0,3 / 0,16 / 0,02.** Model (R_okno 0,75 % × A·ε 0,55) daje 39 / 31 / 4 z zakresem ×0,17…×3,9. Istotność rzędu 2–4σ na gatunek po zsumowaniu pT,
   ≤ 2σ w pojedynczym binie → pomiar per bin pT w 427929 nie jest wykonalny; możliwy tylko jeden punkt pT-inkluzywny na ~30–50 %. 428770 ma 6–7× więcej zdarzeń
   tagowanych z kandydatem (przy pT > 8 GeV 30–60×), ale z trygerem obciążającym i f_acc ≈ 5–6 %.
6. **Wydajność cięcia ξ_cal dla prawdziwych tagów:** 427929: **1,000 +0,000 −0,005** (pomiar; składowa przypadkowa tłumaczy cały ogon); 428770: **0,99 +0,01 −0,03**
   dla nSelBL < 30 (pomiar z modelem mieszania; rozstęp wariantów 0,96–1,00), przy nSelBL ≥ 30 metoda nie działa (ε > 1). Surowy ułamek zdarzeń z D przechodzących cięcie:
   96,6 % (427929) i 77,1 % (428770) — różnica do ε_true to tagi przypadkowe (4,4 % i 27 % tagów przed cięciem). Przeniesienie na sygnał D jest założeniem.

## 1. Dane lokalne

| co | gdzie | 428770 | 427929 |
|---|---|---|---|
| zdarzenia drzewa (wszystkie, także bez kandydata) | `AFP_mix/runs/lxplus/data_<run>_events.parquet` | 7 989 496 (wszystkie w GRL) | 28 962 571 (23 874 370 w GRL; poza GRL ⟨μ⟩ = 0,42) |
| protony: `afp_nProtons`, `afp_nInWindow`, `afp_xi1/2`, `afp_tag` (A, C); `xi_cal_A/C`; LB, μ, GRL, npu, nSelBL | tamże | tak | tak |
| `trigPassed` 11 łańcuchów (po preskali); `trigTBP` 5 pozycji L1 (MBTS_1, MBTS_2, J12, AFP_A_OR_C, RD0_FILLED) | tamże | tak | tak |
| preskale per LB (11 łańcuchów) | `data_<run>.meta.json::trigps`; także `ps_*` per zdarzenie w `B_trig/runs/lxplus/trg_data_<run>.npz` (13 łańcuchów) | 136 LB | 249 LB (189 w GRL) |
| kandydaci D⁺/D⁰/D* slim v2.0 (bez score; pT ≥ 3 GeV, D⁰ L_xy ≥ 0,25 mm) | `data_<run>_cands.parquet` | 50,5 M | 17,6 M |
| score H dla D⁰ | `B_trig/runs/lxplus/bt_data_<run>.npz` (`all_*`) | 30 % zdarzeń z D⁰ | 30 % |
| świetlność per LB, preskale L1 i HLT osobno, pozostałe 95 łańcuchów, decyzje HLT „before prescale” | — | **brak** | **brak** |
| MC z protonem (prawda / AFP) i użytecznym ξ_cal | — | brak lokalnie (§4.4) | — |

Kontrole (`out/s0.log`): `afp_tag_s` ≡ (nProtons_s = 1 ∧ nInWindow_s = 1) w 100 %; ξ1 zawsze w oknie; numer zdarzenia unikalny. Binu 2–3 GeV w slim nie ma (pT ≥ 3).
Definicje: T0 = GRL ∧ tag(A ∨ C); T1 = T0 ∧ npu = 0; **T2 = T1 ∧ ξ_cal,s ≤ ξ_AFP,s + 0,02 na każdej tagowanej stronie** (region sygnału głównego łańcucha, bez wymogu trygera).
GRL: w głównym łańcuchu dekoracja — tu stosowany (w 427929 bez GRL doszłoby 146 tys. zdarzeń z μ ~ 0,4–22).

| run | GRL | T0 | T1 | T2 | T2 z kand. D (slim) | T2 z D⁰ (pT ≥ 3) | T2 z D⁰ score > 0,6 (30 % próbki) | T2: tylko A / tylko C / A ∧ C |
|---|---|---|---|---|---|---|---|---|
| 428770 | 7 989 496 | 1 338 608 | 1 233 602 | 1 134 024 | 114 120 | 76 030 | 2 555 | 570 220 / 321 855 / 241 949 |
| 427929 | 23 874 370 | 707 688 | 703 350 | 701 183 | 18 872 | 10 438 | 353 | 354 983 / 338 503 / 7 697 |

Ułamek T2/GRL: 428770 — zdarzenia 14,2 %, z D 2,13 %, z D⁰ 1,76 %; 427929 — 2,94 %, 0,55 %, 0,51 %. Tag jest dużo częstszy w zdarzeniach bez kandydata (mała krotność).
W 428770 zdarzenia dwustronne (242 tys.) są prawie wyłącznie w klasie `none` (231 tys.) i bez kandydata (z D: 2 388).

## 2. A — skład trygerów i świetlność (`s1_trig.py`, `s1b_psdiag.py`; `out/s1.json`, `out/s1.log`, `out/s1b.log`)

### 2.1 Klasy rozłączne w regionie sygnału T2 (A = noalg_L1AFP ∨ sptrk_pt4_L1AFP, J = noalg_L1J12, M = 7 łańcuchów L1MBTS)

| klasa | 428770 zdarzenia | udział | 428770 z kand. D | udział | 428770 z D⁰ | udział | z D⁰ score > 0,6 | 427929 zdarzenia | udział | 427929 z D⁰ | udział |
|---|---|---|---|---|---|---|---|---|---|---|---|
| AFPnoalg (tylko noalg_L1AFP) | 544 031 | 47,97 % | 11 073 | 9,70 % | 5 563 | 7,32 % | 183 | 92 | 0,01 % | 1 | 0,01 % |
| AFPsptrk4 (tylko sptrk_pt4_L1AFP) | 94 359 | 8,32 % | 53 435 | 46,82 % | 39 454 | 51,89 % | 1 237 | 0 | – | 0 | – |
| AFPboth | 8 769 | 0,77 % | 3 355 | 2,94 % | 2 492 | 3,28 % | 77 | 0 | – | 0 | – |
| AFP ∧ J12 | 10 623 | 0,94 % | 7 851 | 6,88 % | 6 257 | 8,23 % | 305 | 0 | – | 0 | – |
| AFP ∧ MBTS | 3 095 | 0,27 % | 1 474 | 1,29 % | 1 098 | 1,44 % | 33 | 475 972 | 67,88 % | 7 138 | 68,38 % |
| J12only ¬A | 8 725 | 0,77 % | 5 131 | 4,50 % | 3 784 | 4,98 % | 137 | 0 | – | 0 | – |
| MBTS ¬A | 4 760 | 0,42 % | 1 609 | 1,41 % | 1 051 | 1,38 % | 33 | 225 108 | 32,10 % | 3 299 | 31,61 % |
| none (żaden z 11 bitów) | 459 662 | 40,53 % | 30 192 | 26,46 % | 16 331 | 21,48 % | 550 | 11 | 0,00 % | 0 | – |
| **suma** | 1 134 024 | | 114 120 | | 76 030 | | 2 555 | 701 183 | | 10 438 | |

Przed cięciem ξ_cal (T1) udziały są takie same w granicach 1–2 pp (pełne liczby w `out/s1.log`).

### 2.2 Udziały łączne (każdy bit osobno), T2

| łańcuch | 428770 zdarzenia | 428770 z D⁰ | preskala zapisana (LB w GRL) | N·p zdarzenia | N·p z D⁰ | 427929 zdarzenia | 427929 z D⁰ |
|---|---|---|---|---|---|---|---|
| HLT_noalg_L1AFP_A_OR_C | 553 760 | 8 267 | 60 (123 LB), 120 (13 LB) | 34,7 M | 519 540 | 476 064 | 7 139 |
| HLT_mb_sptrk_pt4_L1AFP_A_OR_C | 116 296 | 49 274 | 1 (123 LB), 2 (13 LB) | 121 617 | 51 504 | 7 853 | 2 337 |
| HLT_noalg_L1J12 | 19 779 | 10 273 | 1 | 19 779 | 10 273 | 307 | 176 |
| HLT_j20_L1MBTS_2 | 6 463 | 1 979 | 50 (119 LB), 25 (17 LB) | 290 450 | 88 525 | 3 886 | 1 329 |
| HLT_mb_sptrk_pt4_L1MBTS_2 | 3 883 | 1 675 | 50, 25 | 174 400 | 74 700 | 11 576 | 3 485 |
| HLT_mb_sptrk_pt2_L1MBTS_2 | 358 | 98 | 1250, 625 | 403 125 | 109 375 | 107 363 | 10 381 |
| HLT_noalg_mb_L1MBTS_2 | 958 | 98 | 1250, 625 | 1,10 M | 109 375 | 688 596 | 10 435 |
| HLT_noalg_mb_L1MBTS_1 | 224 | 13 | 5000,66 | 1,12 M | 65 009 | 701 080 | 10 437 |
| HLT_mb_mbts_L1MBTS_1 | 20 | 1 | 50 006,6 | 1,00 M | 50 007 | 700 765 | 10 437 |
| HLT_mb_sptrk_L1MBTS_1 | 0 | 0 | 0 (wyłączony) | – | – | 701 045 | 10 437 |
| HLT_noalg_L1RD0_FILLED | 0 | 0 | 0 | – | – | 0 | 0 |
| bit L1 „before prescale”: AFP_A_OR_C / MBTS_1 / MBTS_2 / J12 | 1 118 935 / 1 044 660 / 939 643 / 19 779 | 71 991 / 75 056 / 73 674 / 10 273 | | | | 476 064 / 701 080 / 688 596 / 307 | 7 139 / 10 437 / 10 435 / 176 |

427929: wszystkie preskale = 1 w 189 LB z GRL (AFP: 10 tylko w 26 LB poza GRL). N·p = liczba zdarzeń × zapisana preskala LB (równoważnik „bez preskali”).
- 428770: bez wag dominuje `sptrk_pt4_L1AFP` wśród zdarzeń z D⁰ (65 % łącznie); z wagami preskali noalg_L1AFP daje 10× więcej (520 tys. vs 52 tys.).
  N·p z D⁰ dla łańcuchów MBTS (50–110 tys.) są zgodne co do rzędu między sobą — to ta sama populacja minimum bias z tagiem, próbkowana z różnymi preskalami.
- W 427929 bit L1_AFP ma 67,9 % zdarzeń T2 (strona A ~0,9, C ~0,45; `AFPTRIG` §2) — tam tryger AFP nie jest potrzebny.

### 2.3 Zależność od LB (8 bloków o równej liczbie T1; udziały w T2)

| run | LB | ⟨μ⟩ | T2/T1 | zdarzenia: noalg / sptrk4 / none / J12only / MBTS¬A [%] | z D⁰: noalg / sptrk4 / AFP∧J12 / none / J12only [%] |
|---|---|---|---|---|---|
| 428770 | 4–24 | 0,0549 | 0,919 | 42,3 / 7,3 / 46,3 / 1,2 / 0,9 | 6,6 / 46,5 / 7,6 / 21,9 / 8,8 |
| 428770 | 25–225 | 0,0518 | 0,921 | 48,5 / 8,5 / 40,0 / 0,7 / 0,4 | 7,3 / 52,3 / 8,4 / 21,4 / 4,3 |
| 428770 | 226–241 | 0,0526 | 0,920 | 48,7 / 8,4 / 40,0 / 0,6 / 0,3 | 7,6 / 53,1 / 8,7 / 21,4 / 3,8 |
| 428770 | 242–257 | 0,0531 | 0,920 | 48,6 / 8,4 / 39,9 / 0,7 / 0,3 | 7,3 / 52,7 / 8,4 / 21,8 / 4,3 |
| 428770 | 258–274 | 0,0543 | 0,918 | 48,5 / 8,5 / 39,9 / 0,7 / 0,4 | 7,4 / 52,9 / 8,0 / 21,6 / 4,6 |
| 428770 | 275–290 | 0,0558 | 0,919 | 48,9 / 8,5 / 39,5 / 0,7 / 0,4 | 7,5 / 52,1 / 8,2 / 21,8 / 4,6 |
| 428770 | 291–306 | 0,0556 | 0,917 | 49,1 / 8,5 / 39,3 / 0,7 / 0,3 | 7,7 / 52,9 / 8,1 / 20,5 / 5,0 |
| 428770 | 307–323 | 0,0549 | 0,920 | 48,9 / 8,5 / 39,5 / 0,7 / 0,3 | 7,1 / 52,3 / 8,4 / 21,4 / 4,8 |
| 427929 | 119–324 (8 bloków) | 0,0073–0,0092 | 0,997 | AFP∧MBTS 67,5–68,3; MBTS¬A 31,7–32,5 | AFP∧MBTS 65,9–70,1; MBTS¬A 29,9–34,1 |

Skład stabilny; jedyna zmiana to pierwszy blok 428770 (LB 4–16: preskale AFP 2× większe → mniej AFP, więcej `none` i J12). Blok „25–225” obejmuje przerwę w GRL (wg `GRL_MinBiasAFP_lowmu2022.xml`: LB 4–25 i 210–324; samej przerwy nie sprawdzałem).

### 2.4 Preskale: zapisane a efektywne (kontrola K-A2)

P(X przeszedł | zdarzenie zapisane przez łańcuch Y ∧ bit L1 „before prescale” X) = 1/p_eff. Dla X = `HLT_noalg_L1AFP_A_OR_C` w 428770:

| odniesienie Y | n_ref | p_eff | zapisana | werdykt |
|---|---|---|---|---|
| noalg_L1J12 (preskala 1, inny L1) | 26 048 | **57,4 ± 2,7** (LB > 16: 55,4 ± 2,7; LB ≤ 16: 93 ± 19) | 60 (120) | zgodne (1,7σ; −8 %) |
| j20_L1MBTS_2 | 8 722 | 63,7 ± 5,4 | 60 | zgodne |
| sptrk_pt4_L1AFP (ten sam L1) | 221 191 | 14,59 ± 0,12 | 60 | **niezgodne** |
| sptrk_pt4_L1MBTS_2 | 5 268 | 16,8 ± 1,0 | 60 | **niezgodne** |
| noalg_mb_L1MBTS_2 | 1 223 | 2,47 ± 0,11 | 60 | **niezgodne** |

- Pozostałe łańcuchy noalg wobec J12: `noalg_mb_L1MBTS_2` 1112 ± 18 (zapisane 1250/625, średnia harmoniczna ≈ 1150), `noalg_mb_L1MBTS_1` 4811 ± 158 (5000,7),
  `mb_mbts_L1MBTS_1` 52 400 ± 5 700 (50 006,6), J12 1,000. 427929: wszystkie 1,000. Zapisane `trigPS` zachowują się więc jak pełna preskala łańcucha.
- Niezgodność dla odniesień `sptrk_pt4_*` zależy od zawartości zdarzenia: p_eff = 10,0 (nSelBL < 20) → 28,6 (≥ 40); 20,5 z kandydatem D, 9,4 bez; **53,9 ± 2,7 gdy
  zdarzenie ma też L1_J12**. Preskala nie może zależeć od krotności, więc to nie jest inna preskala, tylko zdarzenia odniesienia nie są próbką niezależną
  (hipotezy, nietestowalne lokalnie: skim drzewa zależny od trygera/zawartości — w drzewie jest 69 % zdarzeń AOD; albo skorelowane losowania preskal HLT).
  Liczba „≈ 23” z `AFPTRIG` §2 jest mieszaniną takich odniesień. **Wniosek: zapisane 60/120 potwierdza jedyne czyste odniesienie (J12), ale zanim ktokolwiek użyje wag 1/p,
  trzeba to sprawdzić w COOL i zrozumieć skim** — w przeciwnym razie N·p dla łańcuchów, których zdarzenia skim przerzedza, jest zaniżone.

### 2.5 Efektywna świetlność

Świetlności per LB lokalnie nie ma i nie była pobierana. Policzony jest ułamek F_i = Σ_LB w_LB/p_i,LB, w_LB = udział LB w świetlności przybliżony (a) liczbą zdarzeń
łańcucha bez preskali (noalg_L1J12) w LB, (b) ⟨μ⟩ LB przy stałej długości LB. L_eff = F × 29,41 nb⁻¹ (placeholder w głównym łańcuchu, nieweryfikowany).

| łańcuch (428770) | F (a) | F (b) | L_eff [nb⁻¹] | 29,41 / L_eff |
|---|---|---|---|---|
| HLT_noalg_L1J12 | 1,000 | 1,000 | 29,41 | 1 |
| HLT_mb_sptrk_pt4_L1AFP_A_OR_C | 0,9563 | 0,9518 | 28,1 ± 0,1 | 1,05 |
| HLT_noalg_L1AFP_A_OR_C | 0,01594 | 0,01586 | 0,468 ± 0,002 | 63 |
| HLT_j20_L1MBTS_2, HLT_mb_sptrk_pt4_L1MBTS_2 | 0,02229 | 0,02254 | 0,66 ± 0,01 | 45 |
| HLT_mb_sptrk_pt2_L1MBTS_2, HLT_noalg_mb_L1MBTS_2 | 8,92·10⁻⁴ | 9,01·10⁻⁴ | 0,0264 ± 0,0003 | 1 120 |
| HLT_noalg_mb_L1MBTS_1 | 2,00·10⁻⁴ | 2,00·10⁻⁴ | 0,00588 | 5 000 |
| HLT_mb_mbts_L1MBTS_1 | 2,00·10⁻⁵ | 2,00·10⁻⁵ | 0,000588 | 50 000 |

Błąd = |F(a) − F(b)| (samo przybliżenie wag LB); nie zawiera niepewności liczby 29,41, czasu martwego, ani kwestii z §2.4.
427929: F = 1 dla wszystkich łańcuchów w GRL, więc L_eff = L_GRL — której lokalnie nie ma (rząd: 23,9 M zdarzeń MBTS w naszej kopii / σ_widz ~ 70 mb ≈ 0,3 nb⁻¹ dla kopii;
kopia to ~50–72 % runu; to oszacowanie rzędu wielkości, nie liczba do użycia).

**Suma logiczna — wybór metody.** Waga zdarzenia w = 1/(1 − Π_i(1 − 1/p_i)) po łańcuchach, które zdarzenie przeszłoby przed preskalą, wymaga decyzji „before prescale”
każdego łańcucha. Lokalnie są tylko bity L1 TBP (5 pozycji) i decyzje po preskali; dla łańcuchów z warunkiem HLT (sptrk_pt4, j20, mb_mbts) decyzji bez preskali nie ma,
a klasa `none` (40 % zdarzeń) pochodzi z łańcuchów, których nie znamy. Wagi sumy logicznej nie da się więc policzyć poprawnie. **Rekomendacja: rozłączne regiony trygera,
jeden łańcuch na region, L_eff tego łańcucha** (mechanizm `minbias_trigger_regions` już to obsługuje):
- region „proton” = `HLT_noalg_L1AFP_A_OR_C` (bez obciążenia centralnego; 0,47 nb⁻¹; 8 267 zdarzeń T2 z D⁰);
- region „proton + tor” = `HLT_mb_sptrk_pt4_L1AFP_A_OR_C` ∧ ¬noalg (28 nb⁻¹; 49 tys. zdarzeń z D⁰; wymaga wydajności warunku „tor pT > 4 GeV” względem kandydata D — obciążenie
  jak w `B_trig/AFPSTREAM`: SF izolacji 0,65–0,78);
- J12 i klasa `none` poza pomiarem albo jako osobne regiony kontrolne.
Stan obecny (`alltrig`, waga 1, 29,41 nb⁻¹) miesza populacje o L_eff różnych 63× — przekrój z takiego regionu nie jest zdefiniowany, niezależnie od MC.

### 2.6 Czego brakuje i skąd to wziąć (nic nie uruchomiono)

| brak | źródło | narzędzie |
|---|---|---|
| L per LB (dostarczona, zapisana, live fraction) dla 428770 i 427929 z GRL `GRL_MinBiasAFP_lowmu2022.xml` | COOL (`/TRIGGER/OFLLUMI/…`, tag świetlności 2022 dla niskiego μ — do potwierdzenia na stronie LuminosityForPhysics) | `iLumiCalc --xml=<GRL> --lumitag=<tag> --livetrigger=L1_MBTS_1 --trigger=<łańcuch>` (Athena na lxplus) albo serwis atlas-lumicalc (logowanie SSO); wynik: tabela per LB z L po preskali per łańcuch |
| preskale L1 i HLT osobno per LB, klucze preskal | COOL `/TRIGGER/LVL1/Lvl1ConfigKey`, `/TRIGGER/HLT/PrescaleKey` + TriggerDB | `AtlCoolTrigger.py -r 428770 -v`, atlas-runquery; rozstrzyga §2.4 |
| tożsamość łańcuchów klasy `none` (106 łańcuchów jest w ntuplach v2.0, u nas 11) | ROOT v2.0 na EOS (`charm_raw_v20`) | etap ORCH — szkic `SZKIC_stage_tagtrig.sh.txt` (niewysłany, niedopisany do pipeline) |
| decyzje HLT przed preskalą (do wag sumy logicznej) | tylko nowa produkcja CharmAnalysis (`isPassedBits`) | do ustalenia: nowa produkcja; przy rozłącznych regionach niepotrzebne |
| kryterium skimu drzewa (69 % AOD w 428770, 50 % w 427929) | konfiguracja derywacji / CharmAnalysis | do ustalenia: kryterium skimu |

## 3. B — oczekiwany uzysk tagowanych D w 427929 (PRZEWIDYWANIE; `s2_yield.py`, `out/s2_tables.md`)

N_tag = N_incl × P. N_incl z istniejących fitów inkluzywnych (`MC_templates/runs/lxplus/ana/WYNIKI.md` §3, fit C0, „427929:all” = GRL; D* z tagu, §2a) — nasza kopia runu,
nasz klasyfikator H (nie CharmNet), pT 3–30 GeV, D⁰ L_xy ≥ 0,25 mm. Dwa scenariusze P, oba są założeniami:
- **S (kotwica z pasm bocznych):** P = udział kandydatów z pasm bocznych masy w T2 w tym samym gatunku i binie pT (liczby kandydatów, bez piku) × (1 − f_acc);
  założenie „sygnał D tagowany tak często jak kandydaci tła”, przypisany czynnik 0,5–2.
- **M (model):** P = R_okno × A·ε_AFP × ε_ξcal = 0,75 % × 0,55 × 0,99 = 0,41 % (zakres 0,07–1,6 %). R_okno: Pythia poziom cząstek z oknem 0,03–0,10, obie strony
  (`SD_pred` §8: 0,11–0,82 % miękka zależnie od strumienia, 1,1–2,3 % twarda typ 4); literaturowe 1,6–2 % dotyczy ξ < 0,1 bez dolnej granicy (udział okna 20–50 %).
  A·ε_AFP: tory NEAR × FAR 0,90–0,95, dopasowanie > 0,9 (`RECON_program` D2); akceptancji w t/pT protonu w źródłach nie ma — przyjęte 0,6–0,9. To jest główne założenie.
- Przypadkowe: N_incl × P_SB × f_acc, f_acc = 0,0087 ± 0,0044 (`AFP_mix` v2; tu z metody kształtu 0,0083, §4).

| gatunek | pT [GeV] | N_incl WP80 | P_tag w pasmach bocznych [%] | prawdziwe, S (×0,5–2) | prawdziwe, M (zakres) | przypadkowe | √S | Z przy S/B inkl. |
|---|---|---|---|---|---|---|---|---|
| D⁰ | 3–5 | 6773 ± 388 | 0,365 ± 0,005 | 24,5 (12–49) | 27,7 (4,7–107) | 0,21 ± 0,11 | 4,9 | – |
| D⁰ | 5–8 | 2010 ± 120 | 0,354 ± 0,014 | 7,0 (3,5–14) | 8,2 (1,4–32) | 0,06 ± 0,03 | 2,7 | – |
| D⁰ | 8–30 | 813 ± 51 | 0,470 ± 0,045 | 3,8 (1,9–7,6) | 3,3 (0,6–13) | 0,03 ± 0,02 | 1,9 | – |
| **D⁰** | **3–30** | 9596 | | **35 (18–71)** | **39 (7–152)** | 0,31 | 5,9 | – |
| D⁺ | 3–5 | 4999 ± 156 | 0,233 ± 0,003 | 11,5 (5,8–23) | 20,4 (3,5–79) | 0,10 ± 0,05 | 3,4 | 1,3 |
| D⁺ | 5–8 | 1833 ± 62 | 0,232 ± 0,010 | 4,2 (2,1–8,4) | 7,5 (1,3–29) | 0,04 ± 0,02 | 2,1 | 1,3 |
| D⁺ | 8–30 | 767 ± 36 | 0,324 ± 0,033 | 2,5 (1,2–4,9) | 3,1 (0,5–12) | 0,02 ± 0,01 | 1,6 | 1,2 |
| **D⁺** | **3–30** | 7599 | | **18 (9–36)** | **31 (5–120)** | 0,16 | 4,3 | 2,2 |
| D*⁺ (tag) | 3–5 / 5–8 / 8–30 | 513 / 265 / 107 | 0,23 / 0,26 / 0,21 | 1,2 / 0,7 / 0,2 | 2,1 / 1,1 / 0,4 | < 0,02 | ≤ 1,1 | – |
| **D*⁺** | **3–30** | 885 | | **2 (1–4)** | **4 (0,6–14)** | 0,02 | 1,4 | – |

Przy WP50 uzyski są ~2× mniejsze (D⁰ 17, D⁺ 9, D* 1; pełna tabela w `out/s2_tables.md`). √S = górna granica istotności (zero tła); „Z przy S/B inkl.” = S/√(S + B)
z S/B fitu inkluzywnego D⁺ przy WP70 (0,157 / 0,698 / 1,336; `DLZ_dplus/DATA_2026-11-17`), przy założeniu takiego samego S/B w regionie tagowanym.

**Wniosek (wg kryterium z `KRYTERIA.md`: Z ≥ 5 tak, 3–5 na granicy, < 3 nie):**
- per bin pT: **nie** dla każdego gatunku i binu (Z ≤ 2 przy realnym tle; nawet bez tła tylko D⁰ 3–5 GeV sięga 5);
- pT-inkluzywnie: D⁰ „na granicy” tylko przy zerowym tle (√S ≈ 6, realnie 2–4σ, bo S/B D⁰ w 3–5 GeV jest małe), D⁺ ≈ 2σ, D* nie;
- 427929 może dać jeden punkt pT-inkluzywny (D⁰ + D⁺) z błędem statystycznym ~30–50 % oraz kontrolę czystości (f_acc < 1 %, ε_ξcal = 1, tryger nieobciążony) — nie pomiar różniczkowy.
  Kompletny run (nasza kopia to 458 z 583 plików; rząd 50–72 % zdarzeń) dałby najwyżej ×1,4–2 w uzysku, czyli ×1,2–1,4 w istotności; bin 2–3 GeV i CharmNet nie były liczone.

**Wobec 428770** (liczby zdarzeń i kandydatów, nie uzyski): T2 z kandydatem D 114 120 vs 18 872 (×6,0); z D⁰ ×7,3; ze score > 0,6 ×7,2; kandydaci w pasmach bocznych ×8–10 (3–5 GeV),
×18–23 (5–8), ×33–58 (8–30 GeV). Ten sam rachunek kotwicy SB dla 428770 daje 4–7× więcej tagowanych D przy 3–8 GeV i 15–22× przy 8–30 GeV (D⁰ ≈ 220, D⁺ ≈ 135 łącznie) —
orientacyjnie, bo w 428770 skład trygerów sygnału i pasm bocznych nie musi być ten sam. Cena: f_acc 5–6 % (w klasie `none` 8,6 %), 27 % tagów przypadkowych przed cięciem ξ_cal,
tryger obciążający (52 % zdarzeń z D⁰ z `sptrk_pt4_L1AFP`), L_eff różne 63× między klasami. Próbka nieobciążona trygerem w 428770 (bit noalg_L1AFP) ma 8 267 zdarzeń T2 z D⁰ —
mniej niż 427929 (10 438).

## 4. C — wydajność cięcia ξ_cal ≤ ξ_AFP + 0,02 (`s3_xical.py`, `s4_tables.py`; `out/s3.json`, `out/s4_tables.md`)

Δ = ξ_cal − ξ_AFP po stronie tagu, zdarzenia T1. Składowa przypadkowa: szablon mieszania (pula centralna bez protonu po stronie × protony ze zdarzeń z ξ_cal > 0,14,
komórki blok LB × klasa nSelBL, jak `AFPTRIG`), normalizowany do danych w obszarze „proton niemożliwy”. ε_true = [N(Δ ≤ 0,02) − N_acc(Δ ≤ 0,02)] / [N − N_acc].
Nominał: normalizacja w Δ > 0,10; warianty: Δ > 0,045, Δ > 0,20, druga pula (428770: MBTS zamiast wszystkich nie-AFP), bez warstw nSelBL. Błąd stat.: Poisson danych (bootstrap);
statystyka puli (≤ 1,5 M zdarzeń, losowo) pominięta.

### 4.1 Rozkład Δ, zdarzenia tagowane z kandydatem D (A + C) — pomiar

| Δ | 427929 dane | 427929 przypadkowe | 428770 dane | 428770 przypadkowe |
|---|---|---|---|---|
| −0,10…−0,06 | 2 142 | 5,6 | 11 107 | 263 |
| −0,06…−0,02 | 9 902 | 42,5 | 56 714 | 1 975 |
| −0,02…0,00 | 6 426 | 46,8 | 42 734 | 2 026 |
| 0,00…0,02 | 466 | 62,0 | 6 356 | 2 576 |
| **0,02…0,045 (tuż nad cięciem)** | 85 | 92,8 | 4 399 | 3 753 |
| 0,045…0,10 | 205 | 220,0 | 9 474 | 9 015 |
| 0,10…0,20 | 268 | 275,8 | 13 119 | 12 976 |
| > 0,20 | 115 | 107,2 | 7 739 | 7 882 |

427929: nad cięciem dane = szablon (673 vs 696) — brak mierzalnego przecieku prawdziwych protonów. 428770: nadwyżka danych 646 ± 70 w 0,02–0,045 i 459 ± 100 w 0,045–0,10
(to „przeciek własnych protonów” z `AFPTRIG` §4.1), razem ~1 % prawdziwych tagów.

### 4.2 ε_true (A + C)

| run | strumień | próbka | N tagów | surowy ułamek Δ ≤ 0,02 | **ε_true (norm. > 0,10) ± stat** | norm. > 0,045 | norm. > 0,20 | f_acc w SR | udział przypadkowych przed cięciem |
|---|---|---|---|---|---|---|---|---|---|
| 427929 | wszystkie | wszystkie zdarzenia | 711 065 | 0,9969 | 0,9999 ± 0,0001 | 1,0000 | 0,9996 ± 0,0005 | 0,0058 | 0,9 % |
| 427929 | wszystkie | z kand. D (slim) | 19 609 | 0,9657 | **1,0012 ± 0,0014** | 1,0004 | 1,006 ± 0,005 | 0,0083 | 4,4 % |
| 427929 | wszystkie | z kand. w pasmach bocznych | 11 153 | 0,9580 | 1,0030 ± 0,0021 | 1,0003 | 1,016 ± 0,009 | 0,0093 | 5,4 % |
| 427929 | wszystkie | z D⁰ pT ≥ 3 | 10 848 | 0,9654 | 1,0018 ± 0,0020 | 1,0014 | 1,005 ± 0,006 | 0,0080 | 4,4 % |
| 428770 | wszystkie | wszystkie zdarzenia | 1 487 307 | 0,9306 | 0,9810 ± 0,0002 | 0,9931 | 0,9642 ± 0,0004 | 0,031 | 8,1 % |
| 428770 | wszystkie | z kand. D (slim) | 151 642 | 0,7710 | **0,9901 ± 0,0014** (pula MBTS: 0,9842) | 0,9966 | 0,969 ± 0,003 | 0,059 (0,053) | 26,7 % |
| 428770 | wszystkie | z kand. w pasmach bocznych | 103 826 | 0,7382 | 0,9947 ± 0,0018 | 0,9961 | 0,978 ± 0,004 | 0,061 | 30,3 % |
| 428770 | bity AFP | z kand. D | 99 959 | 0,7804 | 0,9835 ± 0,0017 (MBTS: 0,9773) | 0,9951 | 0,960 ± 0,004 | 0,047 | 24,4 % |
| 428770 | sptrk_pt4_L1AFP | z kand. D | 87 010 | 0,7682 | 0,9851 ± 0,0017 | 0,9945 | 0,966 ± 0,004 | 0,048 | 25,8 % |
| 428770 | noalg_L1AFP | z kand. D | 17 389 | 0,8503 | 0,9742 ± 0,0033 | 0,9971 | 0,924 ± 0,005 | 0,038 | 16,0 % |
| 428770 | J12/MBTS bez AFP | z kand. D | 8 776 | 0,7689 | 0,992 ± 0,006 | 0,994 | 0,993 ± 0,014 | 0,051 | 26,5 % |
| 428770 | none | z kand. D | 42 907 | 0,7493 | 1,006 ± 0,003 | 1,000 | 0,985 ± 0,007 | 0,086 | 31,9 % |

Strony A i C zgodne (427929: 1,005 ± 0,002 i 0,997 ± 0,002; 428770: 0,991 ± 0,002 i 0,988 ± 0,003). Pasma boczne per gatunek: `out/s4_tables.md`.

### 4.3 Zależności (ε_true przy normalizacji > 0,10; „z kand. D”)

| przekrój | 427929 | 428770 (wszystkie strumienie) | 428770 bity AFP |
|---|---|---|---|
| nSelBL < 15 | 0,9995 ± 0,0016 (surowy 0,995) | 0,979 ± 0,002 (0,941) | 0,977 ± 0,002 |
| 15–30 | 1,003 ± 0,002 (0,973) | 0,979 ± 0,002 (0,869) | 0,975 ± 0,002 |
| 30–50 | 0,998 ± 0,008 (0,856) | 1,024 ± 0,006 (0,596) | 1,005 ± 0,005 |
| ≥ 50 | 1,05 ± 0,10 (0,370; N = 154) | 1,35 ± 0,06 (0,148) | 1,23 ± 0,05 |
| ξ_AFP 0,03–0,05 / 0,05–0,075 / 0,075–0,10 | 0,997 / 1,002 / 1,005 (± 0,003) | 0,989 / 0,987 / 0,989 (± 0,003) | 0,978 / 0,982 / 0,986 |
| bloki LB (8) | 0,994–1,010 (± 0,004) | 0,983–0,995 (± 0,004) | 0,973–0,992 |
| bez warstw nSelBL | 1,0013 | 0,975 | 0,976 |

- Surowy ułamek spada z krotnością (427929: 0,995 → 0,37; 428770: 0,94 → 0,15), bo rośnie udział tagów przypadkowych i ich ξ_cal; ε_true w 427929 pozostaje 1 w każdej klasie.
- 428770, nSelBL ≥ 30: ε_true > 1 (1,02; 1,35) — szablon przeszacowuje składową przypadkową w SR tam, gdzie stanowi ona > 50 % tagów. Dla tych klas metoda nie jest
  pomiarem; łączna wartość 0,990 jest średnią z 0,979 (nSelBL < 30, 69 % tagów) i wartości niefizycznych. Dla nSelBL < 30 wynik: 0,979 ± 0,002 (stat), warianty 0,955–0,994.
- Brak zależności od ξ_AFP i LB/μ w granicach ±0,01.

### 4.4 Co jest pomiarem, a co założeniem

| stwierdzenie | status |
|---|---|
| Rozkłady Δ i surowe ułamki (4.1, 4.2) | **pomiar** (liczenie zdarzeń) |
| 427929: ogon Δ > 0,02 jest w całości opisany składową przypadkową; ε_true = 1,000 +0,000 −0,005 (stat 0,0014; warianty normalizacji 1,000–1,006, ucięte na 1) | **pomiar** z modelem mieszania; wg kryterium „zmierzona” |
| 428770: ε_true = 0,99 +0,01 −0,03 (nominał 0,990; warianty normalizacji 0,969–0,997, pula −0,006; nSelBL < 30: 0,979) | **pomiar z zastrzeżeniem** — metoda zawodzi przy nSelBL ≥ 30, kształt szablonu nie domyka się w 428770 (`AFPTRIG` K1/K5 FAIL); traktować jako przedział 0,96–1,00 |
| Szablon tła przypadkowego opisuje tagi przypadkowe przy małym ξ_cal (martwe pole ξ_cal ≲ 0,075) | **założenie** metody (nietestowalne; `AFP_mix` §6.3.3) |
| ε cięcia dla zdarzeń z prawdziwym D = ε_true zdarzeń tagowanych z kandydatem tła | **założenie**; wsparcie: brak zależności od nSelBL w 427929, zgodność próbek „z D” / „pasma boczne” / „z D⁰” w ±0,003 (427929) i ±0,005 (428770) |
| Na poziomie cząstek ξ_cal ≤ ξ_p zawsze, więc ε = 100 % i każda strata jest detektorowa (rozdzielczość, szum, pile-up) | **MC poziom cząstek** (`SD_pred` §8 pkt 6), nie detektor |
| Rozdzielczość ξ_cal z MC z pełną symulacją | **nie da się lokalnie**: plików MC z protonami nie mamy; wg dokumentacji produkcji v2.0 ξ_cal w MC jest nieużywalne (mediana 0,002 wobec 0,17 w danych), MC nie ma pile-upu |
| Normalizacja pomiaru: ε_ξcal jako czynnik przy sygnale | propozycja: 427929 — 1,00 z niepewnością 0,5 %; 428770 — 0,98 ± 0,02, per klasa trygera, razem z f_acc z tej samej metody (0,047–0,086 zależnie od klasy) |

## 5. Zastrzeżenia

1. Selekcja kandydatów jest nasza (slim v2.0, pT ≥ 3 GeV, bez score; score H tylko dla 30 % D⁰), nie CharmNet; bin 2–3 GeV nieobecny.
2. Tylko 11 ze 106 łańcuchów: klasa `none` (40,5 % zdarzeń, 21,5 % zdarzeń z D⁰ w T2 428770) jest niezidentyfikowana — 98 % ma bit L1_AFP, więc to inne łańcuchy na L1 AFP;
   nie są to `sptrk_pt6/pt8_L1MBTS_2` (0 zdarzeń).
3. Preskale: zapisane wartości potwierdza tylko odniesienie J12; trzy inne odniesienia dają 2,5–17 z zależnością od krotności (§2.4) — niewyjaśnione, do sprawdzenia w COOL.
4. L_eff to ułamki z preskal × nieweryfikowane 29,41 nb⁻¹; wagi LB przybliżone (różnica wariantów ≤ 0,5 % dla AFP, 1 % dla MBTS). Świetlności 427929 nie ma.
5. Część B: N_incl pochodzi z fitów inkluzywnych z naszym klasyfikatorem i naszą kopią 427929 (kompletność niepotwierdzona); fity D⁰ przy WP80 w 3–5 GeV mają złe χ² (223/33).
   Oba scenariusze ułamka tagu są założeniami; kotwica z pasm bocznych nie jest pomiarem R, a model ma rozrzut ×20. S/B w regionie tagowanym przyjęte jak inkluzywnie (tylko D⁺).
6. Część C: statystyka puli pominięta w błędzie; w 428770 pula centralna pochodzi z trygerów nie-AFP (dla klasy `none` i noalg niezwalidowana); ε > 1 w części klas.
7. Zdarzenia z tagiem po obu stronach liczone po każdej stronie w C; w T2 wymagane cięcie na każdej tagowanej stronie.
8. GRL zastosowany, choć w głównym łańcuchu jest tylko dekoracją; bez GRL w 427929 dochodzą zdarzenia z μ do ~22.

## 6. Propozycje wymagające osobnej decyzji (nie wykonane)

- Etap ORCH `tagtrig` (szkic `SZKIC_stage_tagtrig.sh.txt`): 106 bitów trygera w regionie tagowanym z ROOT v2.0 + L i preskale per LB z COOL (odczyt).
- Fit sygnału w regionie tagowanym per klasa trygera oraz R w 427929 — to byłby już pomiar; poza przyjętym zakresem.

## 7. Pliki (`~/Documents/charm/charm_work/TAGREG_2026-11-18/`)

- `KRYTERIA.md` (przed liczeniem), `RAPORT.md`, `SZKIC_stage_tagtrig.sh.txt`.
- Skrypty: `common.py`, `s0_load.py`, `s1_trig.py`, `s1b_psdiag.py`, `s2_yield.py`, `s3_xical.py`, `s4_tables.py` (kolejność uruchamiania; `~/Documents/charm/.venv_h5/bin/python`, razem ~25 min, RAM ≤ 12 GB).
- Wyniki: `out/s0.log`, `out/s1.json`, `out/s1.log`, `out/s1b.log`, `out/s2.json`, `out/s2_tables.md`, `out/s3.json`, `out/s3.log`, `out/s4_tables.md`; cache `cache/ev_<run>.npz` (kopia kolumn z parquet, 4,1 GB — do usunięcia tylko za zgodą).

## 8. Odstępstwa od `KRYTERIA.md` (po fakcie)

- C: krawędzie histogramu Δ zgrubione do 0,02 i pula centralna ograniczona losowo do 1,5 M zdarzeń (czas liczenia); pierwszy bieg z drobnymi binami przerwany, liczby 427929 identyczne do 4. cyfry.
- C: ε_true > 1 ucięte na 1 w podsumowaniu; dla 428770 dodany podział „nSelBL < 30” jako zakres ważności po obejrzeniu zależności od krotności.
- A: diagnostyka preskali (`s1b_psdiag.py`) i sprawdzenie `sptrk_pt6/pt8` dodane po obejrzeniu niezgodności; wariant „bez GRL” podany tylko jako liczba zdarzeń.
- B: istotność z S/B policzona tylko dla D⁺ (dla D⁰ i D* brak S/B inkluzywnego przy tym WP w raportach) — podana górna granica √S.
