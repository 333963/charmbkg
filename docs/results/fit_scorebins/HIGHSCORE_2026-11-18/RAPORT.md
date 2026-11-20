# HIGHSCORE — deficyt sygnału D⁰←D* przy score H > WP25, pT 8–30 GeV (projekt B, 18.11.2026)

**ANALIZA EKSPLORACYJNA PO FAKCIE.** Efekt znaleziono w danych (R3 §6), a wszystkie podziały, obserwable i testy poniżej wybrano po jego
zobaczeniu. Istotności są lokalne, bez poprawki na liczbę prób; nic tu nie jest kalibracją. Dane inkluzywne (bez tagu AFP). Nic nie pobierano,
nic na lxplus, istniejących plików nie zmieniano. Pełne tabele: `tables.md` (T1–T7); liczby: `out/*.json`.

## 0. Wejścia i metoda
- `B_trig/runs/lxplus/bt_{data_428770,data_427929,fd_nd,mb_nd}.npz` — **mają wszystkie 18 cech H** per kandydat (+ ptcone40/100, nSelBL, bity triggera, prawda MC);
  niczego nie trzeba było łączyć. **Brak ładunku** kandydata/torów → podział po znaku niewykonalny lokalnie (§5).
- Wydajki: zliczenie w oknie |Δm − μ| < 3 MeV minus tło z fitu pasm bocznych (kształt tła i reguła złożoności z `r3_common`), per kategoria score;
  MC w tym samym oknie. Zgodne z fitem jednoczesnym R3 (R_top 0,714 vs 0,724). Błędy: bootstrap po zdarzeniach (200–300 replik) ⊕ statystyka MC.
- Obserwable: R_top = f(>WP25)_dane/f_MC (jak w R3; zależy od N całkowitego, czyli od tła w brudnych kategoriach) oraz
  **ρ50 = [N(>WP25)/N(>WP50)]_dane/MC** i ρ70 (tylko czyste kategorie, niezależne od modelu tła poniżej WP50). ρ50 jest tu miarą główną.
- Prawda MC: `truth` = 1 potraktowano jako nieprompt (z b) — **wniosek z rozkładów** (impact_sig mediana 3,0 vs 1,25; L_xy 1,6 vs 0,77 mm; udział 15 %),
  nie z dokumentacji kodów.

## 1. Odporność (T1)
| wariant (428770, 8–30 GeV) | R_top | ρ50 |
|---|---|---|
| nominalny (okno ±3 MeV) | 0,71 ± 0,08 | **0,741 ± 0,036** (−7σ) |
| okno ±1,5 / ±2 / ±4 MeV | 0,73 / 0,72 / 0,68 | 0,762 / 0,746 / 0,731 (± 0,04) |
| pasma tylko < 155 MeV; tło c2 = 0; sam próg | 0,74 / 0,66 / 0,65 | 0,753 / 0,729 / 0,740 |
| pasma > 5 MeV (tło bez podparcia przy piku) | nieokreślone (N ×1,9) | 0,714 ± 0,038 |
| m(Kπ) ±30 / ±20 MeV | 0,71 / 0,71 | 0,758 / 0,752 |
| J12only / MBTSbased / other | 0,71 ± 0,09 / 0,75 ± 0,15 / 0,52 ± 0,16 | 0,734 ± 0,038 / 0,79 ± 0,10 / 0,79 ± 0,19 |
| nSelBL < 45 / 45–65 / > 65 | 0,63 / 0,72 / nieokr. | 0,71 ± 0,05 / 0,78 ± 0,06 / 0,82 ± 0,12 |
| \|η\| < 1 / ≥ 1; npu = 0; LB < 160 / ≥ 160 | 0,70 / 0,71; 0,71; 0,75 / 0,72 | 0,77 ± 0,05 / 0,69 ± 0,06; 0,74; 0,86 ± 0,09 / 0,72 ± 0,04 |
| **pT 8–12** | 1,01 ± 0,24 | **1,01 ± 0,06** |
| **pT 12–30** | 0,64 ± 0,07 | **0,655 ± 0,053** |
| pT 8–10 / 10–12 / 12–15 / 15–20 / 20–30 (MC z tego samego binu) | – | 1,06 ± 0,07 / 0,98 ± 0,09 / 0,68 ± 0,09 / 0,94 ± 0,12 / 0,72 ± 0,20 |
| J12only: 8–12 / 12–30; MBTSbased: 8–12 / 12–30 | – | 1,02 ± 0,06 / 0,67 ± 0,06; 0,95 ± 0,13 / 0,56 ± 0,17 |
| 5–8 GeV (428770 / 427929) | 0,99 ± 0,16 (fit R3) / 0,99 ± 0,15 | 1,06 ± 0,07 / 0,99 ± 0,07 |
| 427929, 8–30 | 0,94 ± 0,24 (zliczenie; fit R3: 0,74 ± 0,12) | 0,76 ± 0,10 (−2,3σ) |
| 427929: 8–12 / 12–30 | – | 0,84 ± 0,11 / 0,44 ± 0,20 |

- Efekt **nie zależy** od okna i pasm Δm, modelu tła (w ρ), okna masy, klasy triggera, pile-upu, η, okresu runu; zależność od krotności słaba (< 1,5σ).
- Efekt **zależy od pT**: brak go w 8–12 GeV, cały siedzi powyżej 12 GeV — we wszystkich trzech niezależnych próbkach (J12only, MBTSbased, 427929).
  W 5–8 GeV go nie ma (ρ50 = 1,06 ± 0,07 i 0,99 ± 0,07; samo R_top = 0,99 ± 0,16 nie wyklucza 0,72 mocniej niż na 1,7σ; ρ50 odbiega od 0,74 o 4,8σ i 3,6σ).
- R_top w 427929 jest kruche (zliczenie 0,94, fit 0,74 — różnica w N całkowitym 177 vs ~236); ρ50 = 0,76 ± 0,10 stabilne w wariantach (0,72–0,78; przy m ±20/30 MeV 0,85–0,88).

## 2. Główne źródło w 428770: widmo pT wewnątrz binu 8–30 (T2, T5, T6)
- W MC ułamek sygnału powyżej progu WP25 (jeden próg na cały bin 8–30) stromo maleje z pT: 0,283 (8–10), 0,260 (10–12), 0,227 (12–15), 0,153 (15–20), **0,082** (20–30 GeV).
- Widmo sygnału w 428770 (trigger L1_J12) jest dużo twardsze niż w MC: 64 % sygnału powyżej 12 GeV wobec 28 % w MC; wagi dane/MC w pT rosną od 0,37 ± 0,04 (8–9 GeV)
  do 7,4 ± 1,0 (23–30 GeV). W 427929 (MBTS) widmo zgadza się z MC (wagi 0,9–1,2 ± 0,1–0,4).
- MC przeważone do widma pT danych (wagi z wydajków przy score ≥ WP70; z ≥ WP50, ≥ WP80 i z grubszych binów to samo):

| 428770 | dane/MC | samo widmo pT | dane/MC(pT-waż.) |
|---|---|---|---|
| R_top, 8–30 | 0,714 ± 0,068 | **0,775 ± 0,020** | 0,92 ± 0,09 |
| ρ50, 8–30 | 0,741 ± 0,037 | **0,819 ± 0,018** | **0,905 ± 0,048** (−2,0σ) |
| ρ70, 8–30 | 0,713 ± 0,040 | 0,795 ± 0,019 | 0,898 ± 0,054 |
| ρ50, 8–12 | 1,012 ± 0,059 | 0,988 ± 0,006 | 1,025 ± 0,060 |
| ρ50, 12–30 | 0,655 ± 0,059 | 0,829 ± 0,029 | **0,790 ± 0,078** (−2,7σ) |
| χ² kształtu, 14 kategorii, 8–30 | 97,1/13 (p = 6·10⁻¹⁵) | – | 16,5/13 (p = 0,22) |
| χ² kształtu, 14 kategorii, 12–30 | 53,2/13 | – | 29,4/13 (p = 0,006) |

  Czyli z deficytu 0,26 w ρ50 około **0,18 (≈ 2/3) to niedopasowanie widma pT** między danymi J12 a MC przy jednym progu na szeroki bin; po jego
  usunięciu rozkład score w całym 8–30 jest zgodny z MC (p = 0,22), a w ρ50 zostaje 0,905 ± 0,048.
- Dodatkowe przeważenie MC w ptcone40 (otoczenie kandydata; dane J12 są dużo mniej izolowane, χ²/ndf rozkładu ≈ 16–24 nawet po wagach pT) przy zachowaniu widma pT:
  przewidywane ρ50 = 0,786, reszta **0,944 ± 0,053**, χ²(6 kat.) = 2,1/5. W MC przy pT > 12 i ptcone40 > 15 GeV ułamek > WP25 spada do 0,07–0,14 (także dla promptu).
  Samo ptcone40 bez pT nie wystarcza (χ² = 22/5). W 427929 żadne z tych przeważeń nic nie zmienia (ρ50 0,76 → 0,76).

## 3. Przesunięcie czy ogon (T3)
- Surowy wzór w 428770 (8–30): R_c rośnie do 1,2–1,3 w kategoriach ε = 0,25–0,5, potem monotonicznie spada: 0,84; 0,83; 0,87; 0,56 ± 0,08; 0,51 ± 0,08 (górne 5 %).
- Jedno przesunięcie w logit, a = −0,31, opisuje całość: χ² 97,1/13 → 21,2/12; dodanie skali (b = 0,99) lub kompresja samego ogona nie poprawiają istotnie
  (19,6/11; 17,8/11). Rozmycie logit: brak poprawy (95,3/12) → to **przesunięcie całego rozkładu, nie efekt samego ogona ani rozdzielczości score**.
- Po wagach pT przesunięcie resztkowe a = −0,05 (Δχ² = 3,9); w 12–30 GeV a = −0,21 (29,4/13 → 20,2/12, Δχ² = 9,2).
- 427929 (6 kategorii): tożsamość χ² = 8,0/5 (p = 0,16); nadmiar w WP50–25 (1,5 ± 0,3) większy niż deficyt na szczycie — kształt nierozstrzygnięty.

## 4. Cechy (T4) i skład MC (T3, T5)
- Bez wag pT najgorszą cechą H jest `log_pt` (χ²/ndf = 846/9 w regionie ≥ WP70; przeważenie w niej samej daje ρ50 = 0,835, R_top = 0,78), dalej L_xy (141/9)
  i pointing (83/9) — obie przestają odstawać po wagach pT (10/9, 10/9), czyli były wtórne do pT.
- Po wagach pT żadna z 18 cech H nie przekracza χ²/ndf ≈ 2,6 (log_p_K 23,7/9, log_p_pi_min 22,2/9, dedx_pi_max 20,6/9, |cos θ*| 20,0/9, impact_sig 19,6/9);
  przeważenie MC w pojedynczej cesze zmienia ρ50 o ≤ 0,02, raking dwóch najgorszych: 0,832 ± 0,009 (z samego pT: 0,824). **Żadna cecha H ani para cech nie zamyka reszty.**
  Jedyny spójny ślad: impact_sig przesunięte w danych ku większym wartościom (+2,5σ w 428770, +2,1σ w 427929; przeważenie → ρ70 0,763 zamiast 0,799).
  W 427929 odstaje pointing (33/9, −3,8σ; przeważenie → ρ70 0,90) — jedna z 18 cech × kilka regionów, 140 par sygnału.
- Zmienne spoza H (ptcone40/100, ntrk500, nSelBL, lead_pt) różnią się w 428770 bardzo (χ²/ndf 6–37 po wagach pT), w 427929 nie.
- **Nieprompt:** w MC udział 15 %, największy w kategoriach pośrednich (20–23 %), najmniejszy na szczycie (7–10 %); f(>WP25) = 0,265 (prompt) i 0,167 (nieprompt).
  Nawet 100 % nieprompt daje R_top ≥ 0,67; R_top = 0,72 wymagałoby przeskalowania składowej nieprompt ×35 (udział 87 %); ×2 daje 0,956.
  W kształcie: k = 34 bez wag pT (niefizyczne), po wagach pT k = 1,4 bez poprawy χ² (16,2 vs 16,5). **Zmiana udziału nieprompt nie tłumaczy efektu**;
  może dołożyć najwyżej kilka % do reszty.
- **Swap K↔π** (truth 100/101 z prawdziwym D*, 8,5 % sygnału w oknie) ma w MC rozkład score podobny do sygnału i dodany pogarsza zgodność (χ² 114 vs 97) — wykluczony.
- **Przetrenowanie:** prawdziwe pary z mb_nd (poza treningiem) mają f(>WP25) = 0,267 ± 0,032 (8–30, n = 187) i 0,264 ± 0,022 (5–8) wobec 0,250 w fd_nd —
  brak zawyżenia score w próbce treningowej (czułość ograniczona: 0,72 × 0,25 = 0,18 odrzucone na 2,7σ).

## 5. Wnioski i granice
**Najbardziej prawdopodobne źródło (428770):** nie błąd modelowania cech przy score → 1, tylko **niedopasowanie populacji**: sygnał w runie jetowym ma twardsze
widmo pT (i jest mniej izolowany) niż MC, a wydajność ciasnego cięcia silnie zależy od pT w obrębie szerokiego binu 8–30. Samo widmo pT daje R_top = 0,775 ± 0,020
(z obserwowanych 0,71–0,72), z ptcone40 0,75. Przy WP70–95 zależność wydajności od pT jest słaba (SF50 z samego widma: 0,947), stąd SF ≈ 1 tam.

**Wykluczone jako główna przyczyna:** definicja okna/pasm Δm i model tła; okno masy; klasa triggera sama w sobie (MBTSbased i 427929 pokazują ten sam wzór w pT);
pile-up, η, okres runu; efekt samego ogona/rozmycie score; udział nieprompt; swap K↔π; przetrenowanie; pojedyncza źle modelowana cecha H.

**Reszta nierozstrzygnięta:** po uwzględnieniu widma pT zostaje ρ50 = 0,905 ± 0,048 (428770) i 0,76 ± 0,10 (427929, gdzie widmo pT zgadza się z MC) —
po ok. 2σ każde, razem ≈ 0,88 ± 0,04 (−2,8σ lokalnie, po fakcie). Skupia się powyżej 12 GeV (428770: 0,79 ± 0,08; wąski bin 12–15: 0,68 ± 0,09, ale 15–20: 0,94 ± 0,12;
427929: 0,44 ± 0,20). Lokalnie **nie da się** rozdzielić: (i) fluktuacji, (ii) niedokładności MC w zależności score od pT powyżej 12 GeV (MC: 2129 par w 12–30, 282 w 20–30),
(iii) wpływu otoczenia jetowego (zamyka 428770 do 0,94 ± 0,05, ale nie dotyczy 427929), (iv) przesunięcia impact_sig (~2σ w obu runach), (v) kilku % od udziału nieprompt.

**Czego nie dało się zrobić lokalnie:**
- Znak ładunku: brak kolumny w `bt_*.npz` i `EV_dstar/lxplus/ev_*.npz`. Źródło: h5 v2.0 na EOS (`charm_h5_v20/parts/<próbka>`), złączenie po (`file_idx`, `row`),
  albo dopisanie ładunku do listy `keep` w `btrig_extract.py` i ponowna ekstrakcja jako etap ORCH. Nie pobierano.
- Test, czy MC poprawnie opisuje zależność score od pT powyżej 12 GeV, wymaga więcej MC przy wysokim pT lub więcej danych minbias (427929: 57 par w 12–30 GeV).
- Wagi pT zakładają SF(WP70) płaskie w pT (sprawdzone wariantami ≥ WP50/≥ WP80/grube biny: przewidywane ρ50 0,811–0,821); błąd wag pT wchodzi w bootstrap w T2,
  ale nie w χ² kształtu (T3, T6), gdzie wagi są ustalone z nominału.
- Rozkłady cech: odjęcie pasm zakłada rozkład cech tła niezależny od Δm; χ² cech z błędami diagonalnymi; kolumna KS/p w JSON orientacyjna.

**Konsekwencja praktyczna (do decyzji, nie wdrożone):** SF dla punktów pracy ciaśniejszych niż WP50 w 428770 nie należy brać z jednego binu 8–30 z MC o widmie minbias —
albo biny pT ≤ 2–3 GeV powyżej 8 GeV, albo MC przeważone do widma danych; wtedy R_top ≈ 0,92 ± 0,09, ρ50 ≈ 0,91 ± 0,05 zamiast 0,72.

## 6. Pliki
| plik | rola |
|---|---|
| `hs_common.py` | wczytanie, estymator zliczeniowy, wagi pT (importuje `R3_2026-11-18/r3_common.py`, bez zmian) |
| `hs_robust.py` → `out/robust.json`, `out/robust.log` | odporność, 73 warianty, 300 replik (T1) |
| `hs_ptw.py` → `out/ptw.json` | widmo pT: przewidywanie i reszta, bootstrap obejmuje wagi (T2) |
| `hs_shape.py` → `out/shape.json` | 14 kwantyli MC, modele transformacji logit, nieprompt, swap (T3) |
| `hs_feats.py` → `out/feats.json` | 18 cech H + 7 pozostałych, ranking, domknięcie przez przeważenie (T4) |
| `hs_env.py` → `out/env.json` | MC vs ptcone40, test pT ⊕ ptcone40, rozwiązanie dla udziału nieprompt (T5, T6) |
| `hs_tables.py` → `tables.md` | tabele zbiorcze + kontrola przetrenowania na mb_nd (T7) |

Liczenie lokalne ok. 30 min zegarowych (maszyna współdzielona z innym zadaniem; `hs_robust` 14 min).
