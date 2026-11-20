# ISO_ped — izolacja D⁰ z odjętym piedestałem zdarzenia podkładowego (R4), 18.11.2026

Projekt B. Liczenie lokalne (ok. 10 min). Kryteria: `KRYTERIA.md` (zapisane przed obejrzeniem SF wariantów).
Pełne tabele: `out/tables.md`; liczby: `out/iso_ped.json`; werdykt: `verdict.json`.

## Krótko
- **Dane były lokalnie**: `B_trig/TRKJET2/runs/lxplus/tj2_{data_427929,data_428770,fd_nd,mb_nd}.npz` (komplet plików; pary D⁰–π_s) mają stożek surowy
  (`ptcone40`, własny `v_s40_p2`), dwa stożki prostopadłe R = 0,4 przy φ_D ± π/2 (`p_s2_p/m`), stożek losowy (`r_s2`) i krotność (`ntrk05`).
  Torów z Δφ względem D nie ma, więc pełnego obszaru 60°–120° nie zbudowano — piedestał = średnia z dwóch stożków ±90° (podzbiór tego obszaru).
- **Odtworzenie odniesienia** (gałąź ptcone40 < 2 GeV, MC bez wag, pT 3–30): 9:MB 0,977 ± 0,043; 8:all 0,509 ± 0,025; J12only 0,455 ± 0,027 — zgodne z TRKJET2.
- **SF po odjęciu jest bliżej 1 w runie 428770 o 4–6σ (różnica sparowana), ale nie domyka**: 8:all 5–8 GeV 0,60 → 0,74–0,77; 8–30 GeV 0,51 → 0,61–0,68.
  W 427929 (MBTS) bez zmiany (SF zgodny z 1 przed i po).
- **Tło gorzej**: przy tej samej ε_S^MC ε_B na pasmach bocznych danych rośnie o 22–47 % w 8:all i o 16–25 % w 9:MB 5–8 → K3 (≤ 10 %) nie przechodzi.
- **Werdykt wg kryteriów: V1, V2, V3 = „POPRAWA SF, GORSZE TŁO”; R4 łącznie niepotwierdzona, nie domyka.**
- **Test bezpośredni**: hipoteza „pik bliski dobrze opisany, rozbieżność w piedestale” jest prawdziwa tylko w klasach MBTS runu 428770
  (piedestał dane/MC 1,64 ± 0,17, pik bliski 1,04 ± 0,22), nie w J12only (piedestał 1,71 ± 0,08, pik bliski 1,78 ± 0,12 — oba zawyżone).
  W 427929 oba zgodne z MC (1,01 ± 0,07 i 1,00 ± 0,12).

## 1. Definicje
Selekcja, klasy, fit Δm pass/fail (`fit_pf`, kształt z MC truth) i SBS jak w TRKJET2; MC fd_nd truth ∈ {0,1}, przeważone w pT do sygnału klasy.
cone = `v_s40_p2`; perp_avg = (p_s2_p + p_s2_m)/2. Warianty: **V1** cone − perp_avg; **V2** cone − P(ntrk05) (mapa per run z par Δm > 148,5 MeV, MC z fd_nd);
**V3** cone/(perp_avg + 1 GeV); V1r cone − stożek losowy (informacyjnie). Progi dobrane per bin pT tak, by ε_S^MC (bez wag) = ε_S^MC(cone < 2 GeV):

| bin | ε_S^MC | V1 < | V2 < | V3 < |
|---|---|---|---|---|
| 3–5 | 0,759 | 1,40 | 1,24 | 1,28 |
| 5–8 | 0,628 | 1,35 | 1,21 | 1,23 |
| 8–30 | 0,401 | 1,31 | 1,20 | 1,20 |

D = |SF_raw − 1| − |SF_var − 1|, σ_D z bootstrapu po zdarzeniach (300 replik, estymator SBS; MC losowane niezależnie; wagi pT MC stałe).

## 2. SF: surowa vs odjęta (ε_dane z fitu / ε_MC(pT))

| klasa | bin | N_sig | surowa | V1 | V2 | V3 | D/σ_D (V1 / V2 / V3) |
|---|---|---|---|---|---|---|---|
| 9:MB (427929) | 3–5 | 539 | 0,94 ± 0,06 | 0,98 ± 0,06 | 1,02 ± 0,06 | 1,01 ± 0,05 | 0,5 / 0,5 / 0,6 |
| 9:MB | 5–8 | 302 | 1,07 ± 0,06 | 1,04 ± 0,06 | 1,09 ± 0,06 | 1,03 ± 0,06 | 0,8 / −0,7 / 1,0 |
| 9:MB | 8–30 | 116 (niemiar.) | 0,90 ± 0,10 | 0,83 ± 0,09 | 0,93 ± 0,10 | 0,83 ± 0,09 | −0,9 / 0,4 / −0,8 |
| 8:all (428770) | 3–5 | 886 | 0,78 ± 0,07 | 0,96 ± 0,07 | 0,94 ± 0,07 | 0,98 ± 0,07 | 4,7 / 4,6 / 4,7 |
| 8:all | 5–8 | 664 | 0,60 ± 0,05 | 0,74 ± 0,05 | 0,74 ± 0,05 | 0,77 ± 0,05 | 4,8 / 5,9 / 5,1 |
| 8:all | 8–30 | 1246 | 0,51 ± 0,04 | 0,66 ± 0,05 | 0,61 ± 0,04 | 0,68 ± 0,05 | 5,1 / 4,2 / 5,4 |
| 8:J12only | 5–8 | 497 | 0,59 ± 0,06 | 0,75 ± 0,07 | 0,76 ± 0,07 | 0,77 ± 0,07 | 4,2 / 5,7 / 4,5 |
| 8:J12only | 8–30 | 1067 | 0,43 ± 0,04 | 0,60 ± 0,05 | 0,53 ± 0,04 | 0,61 ± 0,05 | 5,1 / 4,0 / 5,5 |
| 8:MBTSbased | 5–8 | 94 (niemiar.) | 0,61 ± 0,09 | 0,80 ± 0,10 | 0,69 ± 0,10 | 0,87 ± 0,10 | 2,2 / 1,5 / 2,6 |
| 8:MBTSbased | 8–30 | 123 (niemiar.) | 0,84 ± 0,11 | 0,89 ± 0,11 | 0,91 ± 0,11 | 0,89 ± 0,11 | 0,5 / 0,8 / 0,5 |
| 8:MBTSbased | 3–30 | 367 | 0,74 ± 0,08 | 0,92 ± 0,08 | 0,85 ± 0,08 | 0,96 ± 0,08 | 3,1 / 2,5 / 3,4 |
| 8:notJ12 | 3–30 | 560 | 0,78 ± 0,06 | 0,92 ± 0,06 | 0,91 ± 0,06 | 0,93 ± 0,06 | 3,5 / 4,0 / 3,5 |

- Estymator SBS daje te same znaki i istotności (kolumny „SF (SBS)”, „D (SBS)” w `out/tables.md`) — brak komórek niejednoznacznych w kryteriach.
- Jedyny fit z det = False: 8:other 8–30 surowa (ε_fit 0,318 ± 0,043; użyto SBS 0,281 ± 0,061). Pozostałe 122 fity det ∧ valid.
- Zależność od trygera maleje, ale nie znika: χ² stałości SF po klasach rozłącznych (ndf 3): 5–8: 38,5 → 12,6 (V1) / 9,7 (V3); 8–30: 30,4 → 13,3 / 12,6; 3–5: 3,4 → 0,5.
- Stożek losowy (V1r) poprawia o połowę mniej niż stożki prostopadłe (8:all 3–30: 0,67 → 0,73 vs 0,82).

## 3. Tło i zysk
ε_B = ułamek D⁰ z pasm bocznych masy (50–75 MeV od m_D⁰, score_H > 0,6, jeden wpis na D⁰) przechodzących cięcie; ta sama ε_S^MC.

| klasa | bin | ε_B surowa (dane) | ε_B var/surowa: V1 | V2 | V3 | to samo na MC (V1/V2/V3) | Z_var/Z_raw dane: V1 | V2 | V3 |
|---|---|---|---|---|---|---|---|---|---|
| 9:MB | 5–8 | 0,364 | 1,21 ± 0,04 | 1,16 ± 0,03 | 1,25 ± 0,04 | 1,05 / 1,10 / 1,06 | 0,88 ± 0,03 | 0,95 ± 0,03 | 0,86 ± 0,03 |
| 9:MB | 8–30 (n = 76) | 0,237 | 1,00 ± 0,12 | 1,11 ± 0,09 | 1,06 ± 0,14 | 1,17 / 1,33 / 1,33 | 0,92 ± 0,09 | 0,98 ± 0,09 | 0,90 ± 0,11 |
| 8:all | 5–8 | 0,113 | 1,38 ± 0,04 | 1,35 ± 0,03 | 1,47 ± 0,04 | jw. | 1,06 ± 0,05 | 1,07 ± 0,04 | 1,06 ± 0,06 |
| 8:all | 8–30 | 0,044 | 1,33 ± 0,12 | 1,22 ± 0,08 | 1,42 ± 0,14 | jw. | 1,12 ± 0,08 | 1,09 ± 0,06 | 1,12 ± 0,08 |

- Z = ε_S(dane)/√ε_B(dane). W 428770 odjęcie podnosi ε_S danych bardziej niż √ε_B, stąd Z_var/Z_raw ≈ 1,06–1,12 (1–1,5σ powyżej 1);
  w 427929 i na MC siła rozróżniania spada o 3–14 % (MC: 0,92–0,98).
- ε_B dane/MC (surowa): 9:MB 3–30 1,11 ± 0,02; 8:all 0,51 ± 0,01; po odjęciu 8:all 0,62–0,65 — tło ma ten sam problem co sygnał i poprawia się podobnie.

## 4. Test bezpośredni: piedestał i pik bliski sygnału D* (odjęcie pasm Δm; MC przeważone w pT), pT 3–30

| klasa | ⟨perp_avg⟩ dane | MC | dane/MC | ⟨cone − perp_avg⟩ dane | MC | dane/MC | F_ped | ⟨ntrk05⟩ dane / MC |
|---|---|---|---|---|---|---|---|---|
| 9:MB | 0,78 ± 0,06 | 0,77 | 1,01 ± 0,07 | 1,34 ± 0,16 | 1,34 | 1,00 ± 0,12 | – | 39,1 / 36,3 |
| 8:all | 1,35 ± 0,05 | 0,80 | 1,69 ± 0,07 | 4,62 ± 0,26 | 2,77 | 1,67 ± 0,10 | 0,23 ± 0,03 | 56,2 / 39,1 |
| 8:J12only | 1,37 ± 0,06 | 0,80 | 1,71 ± 0,08 | 5,35 ± 0,32 | 3,01 | 1,78 ± 0,12 | 0,19 ± 0,03 | 58,1 / 39,5 |
| 8:MBTSbased | 1,29 ± 0,13 | 0,78 | 1,64 ± 0,17 | 2,13 ± 0,45 | 2,05 | 1,04 ± 0,22 | – (mian. < 3σ) | 51,4 / 37,7 |

- Biny pT w 8:all: piedestał dane/MC 2,25 ± 0,19 / 1,79 ± 0,12 / 1,28 ± 0,06 (3–5 / 5–8 / 8–30); pik bliski 1,64 ± 0,83 / 2,19 ± 0,35 / 1,58 ± 0,07;
  F_ped 5–8: 0,24 ± 0,06, 8–30: 0,08 ± 0,02. Przy wysokim pT nadwyżka stożka to prawie w całości pik bliski (dżet wokół D⁰), nie piedestał.
- Mapa P(n) piedestału przy ustalonej krotności jest w danych i MC niemal taka sama (np. ntrk05 30–40: 0,67 / 0,67 / 0,70 GeV dla 427929 / 428770 / MC;
  60–80: 1,48 / 1,53 / 1,52). Zawyżony piedestał w 428770 to więc inny rozkład krotności (⟨ntrk05⟩ 56 vs 39), a nie inna gęstość przy danej krotności.
  (Mapy danych z par tła Δm, MC ze wszystkich par — porównanie przybliżone.)
- χ² rozkładu perp_avg w `out/tables.md` jest zawyżone przez ostatni bin (> 5 GeV), gdzie wagi SBS dają ułamek ≈ 0 z niedoszacowanym błędem
  (np. 9:MB 5–8: 95/7 przy zgodnych średnich) — nie używać; miarodajne są średnie i ułamki.

## 5. Werdykt (wg `KRYTERIA.md`)

| wariant | K1 poprawa ≥ 2σ w 8:all | K2 nie gorzej (8:all, 9:MB) | K3 ε_B ≤ 1,10 | domyka 8:all | werdykt |
|---|---|---|---|---|---|
| V1 | TAK (4,8σ; 5,1σ) | TAK | NIE (1,38; 1,33; 9:MB 5–8 1,21) | NIE | POPRAWA SF, GORSZE TŁO |
| V2 | TAK (5,9σ; 4,2σ) | TAK | NIE (1,35; 1,22; 9:MB 1,16; 1,11) | NIE | POPRAWA SF, GORSZE TŁO |
| V3 | TAK (5,1σ; 5,4σ) | TAK | NIE (1,47; 1,42; 9:MB 5–8 1,25) | NIE | POPRAWA SF, GORSZE TŁO |

**R4 łącznie: niepotwierdzona.** Odjęcie piedestału zdejmuje ok. 1/3 deficytu SF w runie dżetowym (reszta 0,6–0,77 w 5–30 GeV, nadal zależna od klasy),
kosztem 20–45 % większego ε_B. Hipoteza literaturowa pasuje do klas MBTS (tam SF 0,74 → 0,92–0,96, pT 3–30), nie do L1_J12, gdzie tryger zmienia sam pik bliski.
Zgodne z TRKJET2 T1 („piedestał podniesiony, lecz nie tłumaczy deficytu”); nowe: test przy równej ε_S^MC, biny pT, różnica sparowana, ε_B, rozdzielenie MBTS/J12.

## 6. Zastrzeżenia
- Region poprzeczny = dwa stożki R = 0,4 przy ±90° i η_D, nie pełne 60°–120°; szum per zdarzenie większy niż dla pełnego obszaru (V2 go nie ma i daje podobny wynik).
  Pełny obszar wymagałby nowej ekstrakcji z drzewa torów (rozszerzenie `Q_lxplus/trkjet2/tj2_kernel.py` o Σ pT i N w 60° < |Δφ| < 120°, |η| < 2,5; etap ORCH) — nie uruchamiano.
- 9:MB 8–30 (N_sig 116) i MBTSbased w binach (94, 123) poniżej progu 150: K2 w 9:MB 8–30 oparte na komórce niemiarodajnej (D/σ −0,9…0,4).
- ε_B to pary z pasm masy D⁰ z towarzyszącym π_s (zastępnik tła inkluzywnego D⁰); MC tła w 8–30 ma 55 wpisów.
- Wagi pT MC stałe w bootstrapie; błąd SF w tabeli = fit ⊕ MC, σ_D z SBS. V2 używa mapy wyznaczonej na danych tego samego runu (w MC własnej) — to kalibracja, nie czysta predykcja MC.
- Trzy warianty × dwa biny: przy 4–6σ efekt „look-elsewhere” bez znaczenia. Liczby ALICE z R4 nadal niezweryfikowane w PDF.
- Post hoc (poza kryteriami, informacyjnie): interpretacja MBTSbased vs J12only w §4; uwaga o χ². Technicznie: pierwszy bieg padł na bootstrapie
  (start procesów forkserver), fity policzone ponownie i zapisane w `out/fits_cache.pkl`; `out/run_fail1.log` zostawiony.

## Pliki
`KRYTERIA.md`, `iso_ped.py`, `make_tables.py`, `RAPORT.md`, `verdict.json`, `out/iso_ped.json`, `out/tables.md`, `out/fits_cache.pkl`, `out/run.log`, `out/run_fail1.log`, `out/make_tables.log`.
