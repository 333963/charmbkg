# PEAK_lxy010 — pik masy i odbicie K↔π przy selekcji głównego łańcucha (L_xy ≥ 0,10 mm, bin 2–3 GeV), 18.11.2026

Projekt B, zadanie #7 planu prac z 17.11 (`GITLAB_2026-11-17/RAPORT.md` §4 pkt 7). Dane inkluzywne, bez żadnej selekcji AFP.
Nic nie wysłano na lxplus, nic nie pobrano.

## 0. Najważniejsze
1. **Pełnego pomiaru przy selekcji głównego łańcucha lokalnie zrobić się nie da.** Brakuje MC z L_xy < 0,25 mm (każdy lokalny wyciąg MC ma L_xy ≥ 0,25)
   i brakuje jakichkolwiek kandydatów z pT < 3 GeV (lokalnie i na naszym EOS jest tylko produkcja `dpt3to150`). Przygotowana jest paczka ORCH
   (`Q_lxplus/peak_lxy010`, etap `peaklxy`, **niewysłana**), która daje MC i oba runy danych przy L_xy ≥ 0,10. Bin 2–3 GeV wymaga osobnego pobrania (§6).
2. **Co zmierzono lokalnie:** dane z L_xy 0,10–0,25 mm są — run 427929 z klasyfikatorem (h5 v2.0, 226 z 458 plików) i run 428770 bez klasyfikatora
   (kandydaci D* z `AFP_mix`, 97/97 plików). Zmierzono więc **zmianę pozycji i szerokości piku w danych** przy obniżeniu progu 0,25 → 0,10 mm oraz
   dane − MC z kształtem MC przy 0,25 mm. Kontrola silnika: Δm, 8–30 GeV, WP70 = **+10,08 ± 31,74 keV**, k = 1,053 ± 0,036 — identycznie jak w `PEAK_pos/RAPORT.md`.
3. **D⁺: nic się nie zmienia.** Produkcja v2.0 obniżyła próg tylko dla D⁰ i D⁰ w D*; D⁺ ma L_xy ≥ 0,5 mm w v1.1 i v2.0 (wiersze identyczne).
   Liczby PEAK_pos dla m(Kππ) obowiązują przy selekcji głównego łańcucha w 3–30 GeV.
4. **Obniżenie progu nie zmienia istotnie piku m(Kπ)**: δμ(≥0,10 − ≥0,25) = −1,0 ± 0,9 MeV (427929, WP70), −0,3 ± 0,7 MeV (428770, 8–30 GeV, bez klasyfikatora);
   stosunek szerokości 0,97 ± 0,11 i 1,05 ± 0,09.
5. **Δm reaguje na próg w 3–5 GeV:** δμ = −127 ± 53 keV (427929, WP70; 3–30 GeV: −75 ± 32 keV), przy 5–8 i 8–30 GeV zgodne z zerem (−8 ± 18, −28 ± 27 keV;
   428770 bez klasyfikatora: +62 ± 102, +17 ± 25 keV). Powyżej 0,25 mm widać ten sam kierunek: dane − MC = −211 ± 119 / −136 ± 53 / −48 ± 29 keV
   w plasterkach 0,25–0,40 / 0,40–0,70 / ≥ 0,70 mm. Placeholder ±0,2 MeV dla Δm **nie pokrywa** 3–5 GeV przy niskim L_xy.
6. **Odbicie:** przy ostrym cięciu klasyfikatora próg L_xy nie zmienia stosunku odbicie/sygnał (427929: 1,21 ± 0,11 → 1,23 ± 0,11). Przy luźnym cięciu plaster
   0,10–0,25 ma stosunek 2,03 ± 0,37 wobec 1,38 ± 0,14 powyżej 0,25 (1,6σ) — bez MC przy 0,10 nie da się tego zamienić na dane/MC.

## 1. Co jest lokalnie (stan 18.11)
| próbka | L_xy 0,10–0,25 mm | pT 2–3 GeV | uwagi |
|---|---|---|---|
| dane 427929, h5 v2.0 (`DLZ_dplus/DATA_2026-11-17/h5/data_427929`) | **tak** (51 % D⁰ i D* poniżej 0,25) | nie | 226 z 458 plików (rsync przerwany), pełne kolumny → score D0_H, tag D*, ładunki |
| dane 428770 i 427929, `AFP_mix/runs/lxplus/*_cands.parquet` | **tak, tylko wiersze D\*** (D⁰ przycięte do 0,25) | nie | komplet plików; tylko gatunek, masa, Δm, pT, η, L_xy (kwant 0,1 MeV); bez score i ładunków |
| dane 428770, `HYP_d0/runs/lxplus*`, `SWAP_d0/runs`, FEAT/DLZ/CN/B_trig/FMM/DSTWS/ISO | nie (≥ 0,25; D⁺ ≥ 0,5) | nie | wejścia PEAK_pos i SWAP_d0 |
| MC fd_nd, mb_nd — wszystkie wyciągi i h5 v1.2 (`FEAT_charmnet/data`) | **nie** | nie | |
| MC mb_nd ROOT v9 (`~/Documents/charm/mc23_minbiasND_Dpt3to150GeV_v9`, 40 plików) | pozornie tak (brak progu L_xy) | nie | **inna selekcja niż v2.0** — nie użyto, patrz niżej |

MC v9 sprawdzono na pliku `_000002` (te same 246 498 zdarzeń co v1.2 i v2.0): wiersze v1.2 są podzbiorem v9 z identycznymi wartościami, ale w plasterku
0,10–0,25 mm v9 ma 7 144 wierszy D⁰ i 21 prawdziwych, a v2.0 ma 44 624 i 58 (`Q_lxplus/h5_v20/README.md`). To nie jest ta sama selekcja, więc v9 nie zastępuje MC v2.0.

Bin 2–3 GeV: produkcja ntupli v2.0 `dpt2to3` istnieje na gridzie (`MULT_charm/V20_FORMAT.md`), ale etap `rucio_v20` pobrał tylko `dpt3to150`.
Minimalne pT kandydata we wszystkich lokalnych plikach to 3,000 GeV. Nie ma więc czego dopasować — nie mogę podać nawet istotności piku w 2–3 GeV.

## 2. Metoda
- Silnik fitu z PEAK_pos: `fitlib.fit_mc`, `fitlib.fit_data` importowane; `fit_dm` w `dmfit.py` to kopia 1:1 z `PEAK_pos/fit_dstar.py`.
  Sygnał = 2 Gaussy o wspólnym μ, kształt z MC fd_nd, w danych wolne μ i skala k = σ_dane/σ_MC. Δm: tło progowe, okno |m(Kπ) − m₀| < 40 MeV.
  m(Kπ): Czebyszew 2. st., okno |Δm − 145,43| < 1,2 MeV. Hipoteza K/π z ładunku π_s, więc bez odbić.
- **Kształt i pozycja MC zawsze z L_xy ≥ 0,25 mm** (fd_nd, `HYP_d0/runs/lxplus_dstar`). „Δμ przy ≥ 0,10” w tabelach to dane(≥ 0,10) − MC(≥ 0,25).
- 427929: `charmtag.hyp_extract --dstar-only --score-mode model` uruchomiony lokalnie na h5 v2.0 (ten sam moduł co etap `hyp_dstar`), progi WP70 D0_H.
  Wszystkie zdarzenia (bez GRL i bez wymogu trygera — tych kolumn w wyciągu nie ma).
- 428770 przy 0,10: wiersze D* z `AFP_mix` bez klasyfikatora, m(Kπ) = m(Kππ) − Δm, bin po pT(D*), kwantyzacja 0,1 MeV rozmyta losowo ±0,05 MeV.
- Trzy cięcia: ≥ 0,25, ≥ 0,10 i sam plaster 0,10–0,25 (rozłączny z ≥ 0,25). Różnica δ(≥ 0,10 − ≥ 0,25) z bootstrapu wierszy, 200 replik,
  błąd = połowa przedziału 16–84 %.
- „n.m.” (niemierzalne): fit niezbieżny, Z = N_S/σ(N_S) < 3, k na granicy albo σ(k) ≥ 0,5. Wtedy podana jest tylko istotność.
- Odbicie: metoda `SWAP_d0/swap_dstar.py` (okno Δm minus sidebandy), R = N_odbicie(1,75–2,00 GeV) / N_zgodne(1,84–1,89 GeV). Kontrola: 428770, score > 0,857:
  R = 1,022 ± 0,041, mb 1,014, fd 0,904, część niskomasowa 0,428 / 0,364 / 0,316 — jak w `swap_d0.md` §2 (tam 1,02 / 1,01 / 0,90 i 0,43 / 0,36 / 0,32).

## 3. Wyniki — pozycja i szerokość (pełne tabele: `out/tables.md`)
Δμ = dane − MC(fd, ≥ 0,25); Δm w keV, m(Kπ) w MeV; k = σ_dane/σ_MC.

| gatunek / obserwabla | run, selekcja | bin pT | Δμ przy ≥ 0,25 | Δμ przy ≥ 0,10 | k przy ≥ 0,25 | k przy ≥ 0,10 | δμ(0,10 − 0,25) | k(0,10)/k(0,25) | przyrost N_S |
|---|---|---|---|---|---|---|---|---|---|
| Δm(D*) | 428770, WP70 (PEAK_pos) | 3–5 | −211 ± 63 | brak danych | 0,95 ± 0,09 | — | — | — | — |
| | | 5–8 | −110 ± 46 | brak danych | 1,03 ± 0,07 | — | — | — | — |
| | | 8–30 | +10 ± 32 | brak danych | 1,05 ± 0,04 | — | — | — | — |
| Δm(D*) | 427929, WP70 | 3–5 | +60 ± 102 | −68 ± 95 | 1,07 ± 0,15 | 0,96 ± 0,14 | **−127 ± 53** | 0,90 ± 0,08 | 294 → 314 |
| | | 5–8 | −47 ± 82 | −55 ± 86 | 0,99 ± 0,11 | 1,02 ± 0,12 | −8 ± 18 | 1,04 ± 0,04 | 150 → 154 |
| | | 8–30 | +168 ± 126 | +141 ± 122 | 1,09 ± 0,17 | 1,12 ± 0,14 | −28 ± 27 | 1,03 ± 0,03 | 74 → 80 |
| | | 3–30 | +42 ± 64 | −33 ± 65 | 1,06 ± 0,09 | 1,01 ± 0,09 | **−75 ± 32** | 0,96 ± 0,05 | 519 → 549 |
| Δm(D*) | 428770, bez klasyfikatora | 3–5 | n.m. (Z 2,2) | n.m. (Z 2,3) | — | — | — | — | — |
| | | 5–8 | −42 ± 128 | +20 ± 138 | 1,13 ± 0,20 | 1,06 ± 0,20 | +62 ± 102 | 0,94 ± 0,13 | 1198 → 1358 |
| | | 8–30 | −39 ± 42 | −22 ± 43 | 0,95 ± 0,07 | 0,88 ± 0,07 | +17 ± 25 | 0,93 ± 0,05 | 1981 → 2103 |
| m(Kπ) D⁰ z tagu | 428770, WP70 (PEAK_pos) | 3–5 / 5–8 / 8–30 | +2,04 ± 1,89 / −0,34 ± 1,33 / −0,51 ± 0,94 | brak danych | 1,08 ± 0,16 / 0,97 ± 0,11 / 0,92 ± 0,06 | — | — | — | — |
| m(Kπ) D⁰ z tagu | 427929, WP70 | 3–5 | −1,35 ± 2,04 | −2,12 ± 2,07 | 0,59 ± 0,19 | 0,64 ± 0,18 | −0,77 ± 1,31 | 1,08 ± 0,20 | 169 → 224 |
| | | 5–8 | n.m. (Z 4,0, fit niezbieżny) | n.m. | — | — | — | — | — |
| | | 8–30 | n.m. (Z 0,5) | −2,29 ± 3,90 | — | 0,91 ± 0,18 | — | — | — |
| | | 3–30 | −1,52 ± 1,88 | −2,55 ± 1,81 | 0,96 ± 0,16 | 0,94 ± 0,15 | −1,03 ± 0,92 | 0,97 ± 0,11 | 436 → 489 |
| m(Kπ) D⁰ z tagu | 428770, bez klasyfikatora | 8–30 | −0,92 ± 1,18 | −1,19 ± 1,29 | 1,01 ± 0,11 | 1,05 ± 0,11 | −0,28 ± 0,69 | 1,05 ± 0,09 | 2028 → 2389 |
| m(Kππ) D⁺ | 428770, WP70 (PEAK_pos, 30 % lumi) | 3–5 / 5–8 / 8–30 | −0,59 ± 2,07 / +0,09 ± 1,58 / −1,13 ± 0,93 | **to samo** (L_xy ≥ 0,5 w obu) | 0,85 ± 0,12 / 0,93 ± 0,09 / 1,01 ± 0,06 | to samo | 0 | 1 | — |
| wszystko | oba runy | 2–3 | brak kandydatów | brak kandydatów | — | — | — | — | — |

Uwagi do tabeli:
- Sam plaster 0,10–0,25 mm jest mierzalny tylko w jednej komórce: 428770, bez klasyfikatora, 8–30 GeV, Δm: N_S = 186 ± 55 (Z 3,4), Δμ = +128 ± 160 keV,
  k = 0,63 ± 0,17. W 427929 z WP70 plaster ma 47 ± 17 zdarzeń sygnału (Z 2,7; Δm, 3–30 GeV) — nasz klasyfikator (trenowany przy L_xy ≥ 0,25,
  L_xy jest jego cechą) prawie całkowicie go usuwa, dlatego przyrost N_S przy WP70 to tylko +6 %.
- δμ(Δm) = −75 ± 32 keV przy przyroście sygnału o 6 % to więcej, niż może dać sama średnia ważona (plaster musiałby leżeć ~0,9 MeV niżej).
  Część efektu pochodzi więc z dopasowania tła po dodaniu plastra. Traktuję to jako wskazówkę kierunku (2,3σ), nie jako pomiar.
- 427929 przy ≥ 0,25 nie potwierdza przesunięcia Δm w 3–5 GeV z 428770: +60 ± 102 wobec −211 ± 63 keV (2,3σ różnicy; 3–30 GeV: +42 ± 64 wobec −101 ± 27, 2,1σ).
- m(Kπ) bez klasyfikatora poniżej 8 GeV i Δm bez klasyfikatora w 3–5 GeV są niestabilne (k na granicy, Z ≈ 2) — nie używać.

**Trend powyżej 0,25 mm** (428770, pełna lumi, WP70, 3–30 GeV; MC dopasowane w tym samym plasterku; `out/lxy_trend.csv`):

| L_xy [mm] | μ_MC(Δm) [keV] | σ_MC(Δm) [keV] | Δμ(Δm) dane − MC [keV] | k(Δm) | Δμ m(Kπ) [MeV] | k m(Kπ) |
|---|---|---|---|---|---|---|
| 0,25–0,40 | 145 558 ± 11 | 955 | −211 ± 119 | 0,97 ± 0,22 | n.m. (Z 2,9) | n.m. |
| 0,40–0,70 | 145 539 ± 8 | 954 | −136 ± 53 | 0,99 ± 0,08 | −0,57 ± 1,43 | 1,03 ± 0,12 |
| ≥ 0,70 | 145 531 ± 7 | 919 | −48 ± 29 | 1,03 ± 0,04 | +0,18 ± 0,87 | 0,95 ± 0,07 |

MC samo przesuwa się o +27 ± 13 keV między ≥ 0,70 a 0,25–0,40 mm, a dane idą w przeciwną stronę. Różnica skrajnych plasterków to 1,3σ, więc trend
nie jest istotny, ale jest zgodny co do znaku z δμ z 427929. W binach pT nie jest jednolity (8–30 GeV: +220 ± 176 w 0,25–0,40).

## 4. Wyniki — odbicie K↔π (pełna tabela: `out/tables.md`, `out/refl_rows.csv`)
R = odbicie / zgodne w oknie; MC zawsze L_xy ≥ 0,25.

| score > | bin pT | 428770 ≥ 0,25: R; R/R_mb | 427929 ≥ 0,25: R; R/R_mb | 427929 ≥ 0,10: R; R/R_mb | 427929 plaster 0,10–0,25: R | MC mb / fd |
|---|---|---|---|---|---|---|
| 0,857 | 3–30 | 1,02 ± 0,04; 1,01 ± 0,08 | 1,21 ± 0,11; 1,20 ± 0,13 | 1,23 ± 0,11; 1,21 ± 0,14 | 1,33 ± 0,46 | 1,01 ± 0,07 / 0,90 ± 0,01 |
| 0,857 | 3–5 | 1,25 ± 0,10; 1,09 ± 0,14 | 1,40 ± 0,17; 1,22 ± 0,19 | 1,39 ± 0,16; 1,21 ± 0,19 | 1,31 ± 0,51 | 1,15 ± 0,12 / 0,92 ± 0,01 |
| 0,857 | 5–8 | 0,99 ± 0,08; 1,12 ± 0,15 | 1,06 ± 0,16; 1,20 ± 0,22 | 1,07 ± 0,16; 1,21 ± 0,22 | n.m. (6 ± 3 zdarzeń) | 0,88 ± 0,09 / 0,89 ± 0,01 |
| 0,857 | 8–30 | 0,87 ± 0,04; 0,98 ± 0,16 | 0,79 ± 0,17; 0,89 ± 0,23 | 0,81 ± 0,17; 0,91 ± 0,24 | n.m. (1 zdarzenie) | 0,89 ± 0,13 / 0,87 ± 0,02 |
| 0,387 | 3–30 | 1,15 ± 0,05; 1,07 ± 0,10 | 1,38 ± 0,14; 1,28 ± 0,17 | 1,59 ± 0,15; 1,47 ± 0,19 | **2,03 ± 0,37** | 1,08 ± 0,09 / 0,79 ± 0,01 |
| 0,387 | 3–5 | 1,53 ± 0,12; 1,19 ± 0,17 | 1,68 ± 0,24; 1,31 ± 0,25 | 1,95 ± 0,24; 1,53 ± 0,27 | 2,42 ± 0,53 | 1,28 ± 0,16 / 0,81 ± 0,01 |

- Normalizacja: przy ostrym cięciu klasyfikatora obniżenie progu zmienia R o +0,01 (bo plaster prawie nie przechodzi cięcia). Przy luźnym cięciu
  R rośnie z 1,38 do 1,59, a sam plaster ma 2,03 ± 0,37, czyli 1,47 ± 0,31 razy więcej niż powyżej 0,25 mm (1,6σ).
- Stosunki dane/MC przy ≥ 0,25 w binach pT mieszczą się w ~1σ od 1 w obu runach; dokładniejszy pomiar per WP jest w nocie §6.6
  (0,97–0,99 ± 0,04–0,06 w 428770, 0,93–1,01 ± 0,07–0,09 w 427929). Wartości 3–30 GeV dla 427929 (1,20) zależą od widma pT i nie są miarą zgodności.
- MC fd daje R niższe niż MC mb (0,90 wobec 1,01; przy luźnym cięciu 0,79 wobec 1,08) — szablon z próbki forced-decay zaniżałby odbicie.
- Kształt (χ² znormalizowanego histogramu odbić, ndf 24): 428770 ≥ 0,25: 53 względem fd, 21 względem mb. 427929: 14 / 13 (≥ 0,25) i 11 / 14 (≥ 0,10)
  przy 480–540 zdarzeniach — za mało, żeby cokolwiek odróżnić. Plaster względem ≥ 0,25 w tych samych danych (score > 0,387): 19/24.

## 5. Porównanie z placeholderami i propozycja rozmiarów
| systematyka | placeholder | czy pokrywa to, co widać przy ≥ 0,10 | propozycja |
|---|---|---|---|
| PEAK_POSITION, m(Kππ) D⁺ | ±0,2 MeV | selekcja D⁺ bez zmian; dane −0,79 ± 0,74 MeV (PEAK_pos) — zgodne z 0, ale 0,2 MeV nie jest wykazane | **±0,5 MeV** (jak w PEAK_pos) |
| PEAK_POSITION, m(Kπ) D⁰ | ±0,2 MeV | nie wykazane: przy 0,25 −0,10 ± 0,71 MeV, a obniżenie progu dodaje δμ = −1,0 ± 0,9 (427929) i −0,3 ± 0,7 MeV (428770); średnia ważona −0,55 ± 0,55 MeV | **±0,8 MeV** (0,5 ⊕ 0,55) do czasu MC przy 0,10; wpływ na N_S < 0,1 % |
| PEAK_POSITION, Δm(D*) 3–5 GeV | ±0,2 MeV | **nie**: −211 ± 63 keV już przy 0,25 (428770), a próg 0,10 przesuwa dalej w dół (−127 ± 53 keV w 427929) | **±0,3 MeV** albo wolne μ |
| PEAK_POSITION, Δm(D*) 5–8 GeV | ±0,2 MeV | tak: −110 ± 46 (0,25); δμ −8 ± 18, +62 ± 102 keV | **±0,15 MeV** |
| PEAK_POSITION, Δm(D*) 8–30 GeV | ±0,2 MeV | tak, z zapasem: +10 ± 32 (0,25); δμ −28 ± 27, +17 ± 25 keV; plaster +128 ± 160 keV | **±0,10 MeV** |
| PEAK_WIDTH, m(Kππ), m(Kπ) | ±10 % | tak: D⁺ 0,97 ± 0,05; D⁰ 0,95 ± 0,05 (0,25); k(0,10)/k(0,25) = 0,97 ± 0,11 i 1,05 ± 0,09 | **±10 % zostawić** dla D⁰ (przy 0,10 dokładność danych to ~10 %); D⁺ można zwęzić do ±6 % |
| PEAK_WIDTH, Δm(D*) | ±10 % | tak: 1,04 ± 0,03 (0,25); k(0,10)/k(0,25) = 0,96 ± 0,05 (427929), 0,93 ± 0,05 (428770, 8–30), razem 0,945 ± 0,034; skrajna komórka przy 0,10: k = 0,88 ± 0,07 (428770, 8–30, bez klasyfikatora) | **±10 % zostawić** (przy 0,25 wystarczało ±5 %; przy 0,10 nie zwężać przed MC) |
| Reflection, normalizacja | brak (MC, sztywna) | przy ≥ 0,25: dane/MC = 1 w granicach 5 % (428770) i 7–9 % (427929) — nota §6.6. Plaster 0,10–0,25: niezweryfikowany, R wyższe o czynnik 1,47 ± 0,31 przy luźnym cięciu | **±10 %** dla ≥ 0,25; **±20 %** dopóki nie ma MC przy 0,10; lepiej związać `Reflection` z POI sygnału (R z MC minbias, nie forced-decay) |
| Reflection, kształt | brak | niepotwierdzony już przy 0,25 (χ²/ndf 53/24 względem fd w 428770); przy 0,10 brak statystyki | wariant kształtu: MC mb zamiast fd albo kształt z danych z tagu D* (nota §6.6: 4–17 % N_S) |
| wszystko, 2–3 GeV | jak wyżej | nieznane — brak danych i MC | bez liczby; patrz §6 |

Rozdzielenie Δm od masy D⁰ jest konieczne: skala niepewności pozycji to 0,1–0,3 MeV dla Δm i 0,5–0,8 MeV dla m(Kπ) i m(Kππ), a w głównym łańcuchu
jest jeden parametr ±0,2 MeV dla obu (`src/DSelection.cxx:27`). Dla Δm w 3–5 GeV przesunięcie 0,3 MeV przy σ ≈ 0,95 MeV i ustalonym μ to ΔN/N ≈ 5 %,
więc tam najlepiej puścić μ wolno albo przesunąć szablon o zmierzoną wartość.

## 6. Paczka ORCH (przygotowana, NIE wysłana) i co potem
- Paczka: `Q_lxplus/peak_lxy010/` (README, `pl_common.sh`, `pl_job.sh`, `pl.sub`, `make_dag_pl.sh`, `submit_pl.sh`, `src/` = kopia `hyp_extract/src`).
  Wzorzec: `hyp_extract` + etap `hyp_dstar`, czyli to, co zrobiło wejścia PEAK_pos i SWAP_d0. Zmiany: wejście h5 v2.0 (`charm_h5_v20/parts`),
  cztery próbki (fd_nd, mb_nd, data_428770, data_427929), score z modeli foldów, wyjście `charm_q/peak_lxy010`, batch `charm_peaklxy`,
  kontrola w `finish` (STOP bez DONE, gdy próbka nie ma L_xy < 0,25 albo score = NaN).
- Etap: `ORCH/stages/peaklxy.sh` → pobiera do `PEAK_lxy010_2026-11-18/runs/lxplus/`. Wpisu do `pipeline.yaml` nie dodano, `state/` nie utworzono.
- Sprawdzone lokalnie: `bash -n` (5 plików), `py_compile`, generacja `pl.dag` i komendy węzłów w piaskownicy (`DRYRUN=1`), węzeł `finish` na
  lokalnym parquecie 427929, a sama ekstrakcja `charmtag.hyp_extract --dstar-only --score-mode model` na 226 plikach h5 v2.0 (1,5 min, 645 tys. wierszy).
- **Po biegu:** w `fit_lxy.py` i `refl_lxy.py` podmienić źródła na `runs/lxplus/{fd_nd,mb_nd,data_428770,data_427929}.parquet`, liczyć kształt MC
  osobno dla każdego cięcia L_xy (teraz zawsze ≥ 0,25) i dodać dane − MC w samym plasterku. Wtedy 428770 ma WP70 i ładunki także przy 0,10.
- **Bin 2–3 GeV — osobny, nieprzygotowany etap:** pobranie kontenerów `…Dpt2to3GeV.v2.0.Oct28_EXT0` (dane obu runów i MC), konwersja h5 do osobnego
  katalogu (te same zdarzenia są w obu produkcjach), obniżenie twardego pT ≥ 3 w `charmtag.samples.in_range`. Model D0_H nie ma progów WP dla 2–3 GeV
  (trenowany przy pT ≥ 3), więc potrzebny byłby też wybór punktu pracy. Z tego, co widać w 3–5 GeV (Δm: największe przesunięcie; m(Kπ) z tagu: Z 3,4–4,8),
  w 2–3 GeV pik m(Kπ) z tagu D* będzie na granicy mierzalności.

## 7. Zastrzeżenia
1. Brak MC z L_xy < 0,25: każde „dane − MC przy 0,10” zakłada, że MC nie zmienia się poniżej 0,25. Powyżej 0,25 MC zmienia się o +27 keV (Δm) i do ~1,4 MeV
   w pojedynczych komórkach m(Kπ) (5–8 GeV), σ_MC o ±4 %. Liczby z §5 dla progu 0,10 są więc wstępne.
2. Nasz klasyfikator usuwa plaster 0,10–0,25, a CharmNet (główny łańcuch; trenowany na v2.0) może go zachować. Wtedy rzeczywisty efekt w głównym łańcuchu leży między naszymi
   liczbami WP70 i „bez klasyfikatora”, a te drugie są mierzalne tylko w 8–30 GeV.
3. 427929: 49 % plików, bez GRL i bez wymogu trygera; statystyka D* to 519 zdarzeń sygnału przy WP70 (428770: 2803), więc biny pT mają błędy 80–130 keV i 2–4 MeV.
4. 428770 przy 0,10: bez score i ładunków, bin po pT(D*) zamiast pT(D⁰), brak usuwania konfliktów tagu.
5. D⁰ mierzone tylko przez tag D*; inkluzywne D⁰ nie nadaje się do tych systematyk (PEAK_pos §5).
6. Bootstrap po wierszach, nie po zdarzeniach; w 8–30 GeV 40 z 200 replik odrzuconych (fit niezbieżny).
7. Błędy statystyczne z Hesse; rozrzut od modelu tła (PEAK_pos: 11–34 keV dla Δm, 0,3–0,8 MeV dla m) nie jest dodany.

## Pliki
`run_local.sh` (całość, ~12 min), `fit_lxy.py`, `dmfit.py`, `lxy_trend.py`, `refl_lxy.py`, `summarize.py`;
`out/tables.md`, `out/summary.json`, `out/lxy_rows.csv`, `out/lxy_boot.csv`, `out/lxy_trend.csv`, `out/refl_rows.csv`, `out/refl_shape.csv`,
`out/dstar_data_427929_v20.parquet` (+ `.meta.json`), `logs/`. Paczka: `../Q_lxplus/peak_lxy010/`, etap: `../ORCH/stages/peaklxy.sh`.

---

## §7 Wynik z MC przy 0,10 mm (18.11 wieczór)

Wejścia: `runs/lxplus/{fd_nd,mb_nd,data_428770,data_427929}.parquet` (etap `peaklxy`, h5 v2.0, L_xy ≥ 0,10, score D0_H z modeli foldów).
Skrypty i tabele: `v2/` (`run_v2.sh`, ~10 min). Sekcje 0–7 powyżej są wersją wstępną; tam, gdzie się różnią, obowiązuje §7.
Dane inkluzywne, bez selekcji AFP. Nic nie wysłano.

### 7.0 Kontrola wejść (`v2/check_inputs.py` → `v2/out/check_inputs.txt`)
| próbka | pliki (meta / oczekiwane) | wiersze (parquet = meta = summary) | L_xy min; udział < 0,25 | wiersze ≥ 0,25 wobec wyciągu v1.2 (`HYP_d0/runs/lxplus_dstar`) | score |
|---|---|---|---|---|---|
| fd_nd | 24 / 24, każdy z wierszami | 555 118 | 0,100; 44,8 % | 306 349 = 306 349, dopasowane 100 % | identyczny z OOF (max różnica 0) |
| mb_nd | 78 / 78 | 844 712 | 0,100; 54,8 % | 382 212 z 391 904 (v1.1 miało 80 plików; brak 2,5 %), 0 nowych | identyczny |
| data_428770 | 97 / 97 | 4 428 035 | 0,100; 50,1 % | 2 209 990 dopasowanych, 2 nowe, 0 brakujących | identyczny |
| data_427929 | 458 / 458 | 1 519 972 | 0,100; 51,4 % | lokalny wyciąg z 226 plików (644 976 wierszy) zawarty w 100 %, score identyczny | — (run poza treningiem) |

Brak NaN w score, brak zduplikowanych wierszy, pT 3–30 GeV, `frac = 1`.

**Klucz (run, event) w fd_nd.** Potwierdzone: `run = 0`, 10 000 wartości `event`; 51 % kluczy (plik, run, event) obejmuje więcej niż jedno zdarzenie,
a 72 % wierszy leży w takich kluczach (302 839 bloków zdarzeń wobec 172 100 kluczy).
- *Parowanie D⁰ ↔ π_soft nie używa tego klucza.* `charmtag.dstar_tag` łączy po `evblk` (ciągły blok wierszy h5 o tym samym (plik, run, event))
  i po zbiorze trackId obu córek. Dwa zdarzenia skleiłyby się tylko wtedy, gdy sąsiadują w pliku i mają ten sam numer. Test: w zdarzeniu D* stoi
  przed D⁰, więc para sklejona „do przodu” miałaby `dstar_row > row` — jest ich **0 z 555 118** (także 0 w mb i w danych); maks. odległość
  wierszy D⁰–D* to 107. Sklejenie „do tyłu” jest symetryczne, więc górna granica to ~6 wierszy (95 % CL), czyli < 1,1·10⁻⁵. Liczb to nie zmienia.
  (Dla porównania `bt_fd_nd.npz` łączone po (run, event): 0,062 % — `AUDYT_evjoin`.)
- *Fold wynika z (run, event)*: `fold = splitmix(run, event) mod 5`, w fd_nd więc tylko z numeru zdarzenia (5 grup po ~2000 numerów, udziały 19,8–20,2 %).
  Co to znaczy dla „out-of-fold”: wszystkie zdarzenia fizyczne o tym samym numerze trafiają do jednego foldu, tak samo w treningu
  (`models/D0_H/oof.h5`: sygnał `run = 0`, `event` 1…10000, ten sam hash — sprawdzone) i przy aplikacji. Kandydat jest więc oceniany modelem, który
  nie widział ani jego, ani żadnego zdarzenia o tym numerze — score **jest** poza foldem treningowym, tylko grupa foldu to „numer zdarzenia”,
  a nie „zdarzenie”. Dowód liczbowy: dla wszystkich 306 349 wierszy z L_xy ≥ 0,25 score z modeli foldów = score OOF z joinu (różnica 0).
  Wiersze z plastra 0,10–0,25 nie były w treningu w ogóle (trening na v1.1, L_xy ≥ 0,25). Dla danych 428770 (tło treningu) fold z (run, event) jest jednoznaczny.
- Czego nie wolno robić na tych plikach: grupować ani liczyć zdarzeń fd_nd po (run, event) lub (plik, event) — trzeba po `evblk`. Bootstrap w `v2/` idzie po `evblk`.

### 7.1 Kontrola i metoda
- Odtworzenie liczby opublikowanej z nowych plików: Δm, 428770, WP70, 8–30 GeV, ≥ 0,25: **+10,08 ± 31,74 keV, k = 1,053 ± 0,036** (identycznie jak
  `PEAK_pos/RAPORT.md`). Odbicie, 428770, score > 0,857, ≥ 0,25: R = 1,022 ± 0,041 (analityczny), fd 0,904, mb 1,016 (było 1,014 przy 80 plikach).
- Zmiany metody względem wersji wstępnej: kształt i pozycja MC fd_nd liczone **osobno dla każdego cięcia L_xy**; 428770 ma teraz score, ładunki
  i pT(D⁰) także przy 0,10; 427929 ma komplet 458 plików; bootstrap po zdarzeniach (Poisson na `evblk`, 200 replik) zamiast po wierszach;
  mb_nd jako wariant na trzy sposoby: sygnał z prawdy względem fd, wszystkie kandydaty dopasowane jak dane (pseudo-dane), R i kształt odbicia.
- Błąd R z bootstrapu jest o 15–35 % mniejszy od analitycznego (wiersz zgodny i odbity tej samej pary są w jednym zdarzeniu) — w tabelach bootstrapowy.

### 7.2 Pozycja i szerokość — dane − MC fd przy tym samym cięciu (WP70; pełne tabele `v2/out/tables_v2.md` §1–3)
Δm w keV, m(Kπ) w MeV. „oba” = średnia ważona runów.

| obserwabla | bin pT | 428770 ≥ 0,25 | 428770 ≥ 0,10 | 427929 ≥ 0,25 | 427929 ≥ 0,10 | oba ≥ 0,25 | **oba ≥ 0,10** | k oba ≥ 0,25 → ≥ 0,10 | δμ_dane − δμ_MC po dodaniu plastra (oba) |
|---|---|---|---|---|---|---|---|---|---|
| Δm(D*) | 3–5 | −211 ± 63 | −161 ± 73 | +21 ± 69 | −30 ± 67 | −105 ± 47 | **−89 ± 49** | 1,03 ± 0,07 → 1,01 ± 0,08 | −6 ± 26 |
| | 5–8 | −110 ± 46 | −108 ± 46 | −9 ± 58 | −10 ± 58 | −70 ± 36 | **−71 ± 36** | 1,04 ± 0,05 → 1,04 ± 0,05 | +1 ± 8 |
| | 8–30 | +10 ± 32 | +11 ± 32 | +92 ± 89 | +85 ± 88 | +19 ± 30 | **+19 ± 30** | 1,08 ± 0,03 → 1,08 ± 0,03 | +1 ± 2 |
| | 3–30 | −101 ± 27 | −86 ± 30 | +21 ± 44 | −10 ± 45 | −68 ± 23 | **−63 ± 25** | 1,05 ± 0,03 → 1,04 ± 0,04 | +2 ± 12 |
| m(Kπ) z tagu | 3–5 | +2,04 ± 1,89 | +2,65 ± 2,06 | −2,08 ± 1,54 | −2,30 ± 1,49 | −0,44 ± 1,19 | **−0,60 ± 1,20** | 0,96 ± 0,12 → 0,96 ± 0,11 | +0,05 ± 0,69 |
| | 5–8 | −0,34 ± 1,33 | +0,12 ± 1,34 | n.m. | n.m. | −0,34 ± 1,33 | **+0,12 ± 1,34** | 0,97 ± 0,11 → 0,98 ± 0,11 | +0,46 ± 0,35 |
| | 8–30 | −0,51 ± 0,94 | −0,52 ± 0,93 | −0,49 ± 2,62 | n.m. | −0,51 ± 0,88 | **−0,52 ± 0,93** | 0,93 ± 0,05 → 0,92 ± 0,06 | −0,01 ± 0,05 |
| | 3–30 | +0,51 ± 0,77 | +0,83 ± 0,82 | −1,57 ± 1,18 | −1,74 ± 1,20 | −0,11 ± 0,65 | **+0,00 ± 0,68** | 1,00 ± 0,05 → 1,02 ± 0,06 | +0,17 ± 0,28 |

- Sam plaster 0,10–0,25 przy WP70 jest niemierzalny w binach pT (N_S: 86 ± 63 i 59 ± 20 w 3–5 GeV; 3–30 GeV, 427929: Δμ = −76 ± 140 keV).
  Klasyfikator zostawia w plasterku 3–4 % sygnału (428770: 2803 → 2899). Jedyne komórki plastra z pomiarem: Δm bez klasyfikatora,
  428770 8–30 GeV: +183 ± 138 keV, k = 0,64 ± 0,23; 427929 3–30 GeV: −68 ± 266 keV.
- Bez klasyfikatora (oba runy, ≥ 0,10): Δm +61 ± 90 (5–8), −1 ± 39 keV (8–30), k = 1,17 ± 0,15 i 1,00 ± 0,06; m(Kπ) −1,01 ± 1,18 MeV, k = 1,05 ± 0,10 (8–30).
  3–5 GeV bez klasyfikatora nadal niemierzalne (k na granicy albo Z < 3 w co najmniej jednym cięciu).
- Runy przy ≥ 0,25 w 3–5 GeV różnią się o 2,5σ (−211 ± 63 wobec +21 ± 69 keV; χ² = 6,2/1). Średnia −105 ± 47 nie oddaje tego rozrzutu;
  ze skalowaniem błędu (√χ²) to −105 ± 116 keV.
- MC w samym MC (fd, sygnał z prawdy): dodanie plastra nie rusza Δm (δμ_MC ≤ 2 keV, σ(≥0,10)/σ(≥0,25) = 0,99–1,01). W m(Kπ) plaster leży niżej
  o 1,1–1,4 MeV bez klasyfikatora (δμ_MC = −0,33 ± 0,06 MeV dla ≥ 0,10; przy WP70 −0,02 ± 0,03) i jest o 3–7 % węższy.
- Wariant MC mb_nd (sygnał z prawdy względem fd, to samo cięcie): Δm −26 ± 32 keV (WP70, 3–30, ≥ 0,25), −19 ± 31 (≥ 0,10), szerokość 0,98 ± 0,04;
  w binach |Δμ| ≤ 74 keV przy błędach 34–75. m(Kπ): **+1,36 ± 0,87 MeV** (WP70, 3–30), +1,23 ± 0,66 (bez klasyfikatora, ≥ 0,10), σ_mb/σ_fd = 1,01–1,04.
  Czyli dane − mb dla m(Kπ) ≈ −1,3 ± 1,1 MeV. Różnica mb − fd ma 1,6–1,9σ i nie zależy od progu L_xy.

### 7.3 Odpowiedzi na pytania z wersji wstępnej
1. **Przesunięcie Δm o −127 ± 53 keV w 3–5 GeV po dodaniu plastra nie potwierdza się i nie ma go w MC.**
   - MC fd (prawda): δμ_MC = −0 ± 3 keV (plaster − reszta: −1 ± 19 keV). Selekcja fizycznie nie przesuwa piku.
   - 427929 z kompletem plików: δμ = −51 ± 35 keV (było −127 ± 53 na 226 plikach); 428770 z WP70: **+50 ± 40 keV** (przeciwny znak).
     Średnia obu runów po odjęciu MC: **−6 ± 26 keV**; 3–30 GeV: +2 ± 12 keV.
   - MC mb_nd dopasowane jak dane (z tłem, statystyka jak w danych): δμ = +76 ± 39 keV w 3–5 GeV przy zerowym przesunięciu prawdziwym.
     Takie odchylenia daje więc sama procedura (dopasowanie tła po dodaniu plastra o innym S/B). Wynik wstępny był fluktuacją tego typu,
     a nie niezgodnością dane/MC. Uwaga z §3 („część efektu z dopasowania tła”) była trafna.
   - Trend „dane coraz niżej przy małym L_xy” z plasterków powyżej 0,25 mm (§3) nie przedłuża się do 0,10: Δμ(≥ 0,10) = −89 ± 49 wobec −105 ± 47 keV.
2. **R w plasterku 0,10–0,25 jest większe także w MC — ale tylko w minbias, nie w forced-decay.**

   | selekcja, 3–30 GeV | 428770: R plaster / R ≥ 0,25 | 427929 | MC mb_nd | MC fd_nd |
   |---|---|---|---|---|
   | WP70 | 1,88 ± 0,34 / 1,04 ± 0,04 (×1,80 ± 0,32) | 1,89 ± 0,30 / 1,18 ± 0,06 (×1,60 ± 0,27) | 2,53 ± 0,68 / 1,01 ± 0,05 (×2,50 ± 0,70) | 1,14 ± 0,02 / 0,885 ± 0,005 (×1,29 ± 0,02) |
   | bez klasyfikatora | 2,15 ± 0,16 / 1,05 ± 0,05 (×2,05 ± 0,18) | 1,66 ± 0,13 / 1,41 ± 0,10 (×1,18 ± 0,14) | 1,95 ± 0,20 / 1,17 ± 0,11 (×1,66 ± 0,24) | 0,81 ± 0,01 / 0,735 ± 0,004 (×1,10 ± 0,01) |
   | score > 0,387 (jak §4) | ×2,27 ± 0,20 | ×1,47 ± 0,18 (było 1,47 ± 0,31) | ×2,05 ± 0,34 | ×1,21 ± 0,01 |

   Podwójny stosunek dane/mb: 0,67 ± 0,20 (WP70), 0,90 ± 0,14 (bez klasyfikatora), 0,89 ± 0,16 (score > 0,387) — zgodne z 1.
   Dane/fd: 1,31 ± 0,16, 1,36 ± 0,10, 1,51 ± 0,11 — fd wzrostu nie odtwarza. Wzrost R w plasterku to cecha próbki minimum-bias
   (w oknie Δm pod pikiem jest tam więcej par, w których dobra jest tylko jedna hipoteza), a nie niezgodność dane/MC.

### 7.4 Odbicie — dane/MC przy tym samym cięciu (`v2/out/tables_v2.md` §4–6)
| selekcja | bin | ≥ 0,25: 428770 / 427929 względem mb | **≥ 0,10: 428770 / 427929 względem mb** | plaster względem mb | ≥ 0,10 względem fd | χ² kształtu ≥ 0,10 vs fd / vs mb (ndf 24): 428770; 427929 |
|---|---|---|---|---|---|---|
| WP70 | 3–30 | 1,03 ± 0,06 / 1,17 ± 0,08 | **0,97 ± 0,07 / 1,11 ± 0,08** | 0,74 ± 0,24 / 0,75 ± 0,23 | 1,22 ± 0,05 / 1,40 ± 0,07 | 34 / 13; 27 / 21 |
| WP70 | 3–5 | 1,16 ± 0,12 / 1,19 ± 0,13 | **1,07 ± 0,11 / 1,09 ± 0,11** | 0,76 ± 0,22 / 0,76 ± 0,22 | 1,51 ± 0,09 / 1,54 ± 0,10 | 29 / 19; 27 / 26 |
| WP70 | 5–8 | 1,08 ± 0,11 / 1,18 ± 0,12 | **1,10 ± 0,11 / 1,17 ± 0,12** | n.m. | 1,12 ± 0,06 / 1,20 ± 0,07 | 22 / 21; 18 / 14 |
| WP70 | 8–30 | 0,96 ± 0,11 / 0,84 ± 0,12 | **0,96 ± 0,11 / 0,88 ± 0,13** | n.m. | 1,03 ± 0,04 / 0,93 ± 0,09 | 32 / —; — |
| bez klasyfikatora | 3–30 | 0,89 ± 0,10 / 1,21 ± 0,14 | **0,95 ± 0,08 / 1,02 ± 0,10** | 1,10 ± 0,14 / 0,85 ± 0,11 | 1,87 ± 0,07 / 2,02 ± 0,12 | 59 / 20; 24 / 18 |
| bez klasyfikatora | 3–5 | 0,88 ± 0,14 / 1,17 ± 0,21 | **0,97 ± 0,11 / 1,01 ± 0,12** | 1,11 ± 0,16 / 0,89 ± 0,14 | 2,10 ± 0,11 / 2,17 ± 0,15 | 45 / 15; 27 / 17 |
| bez klasyfikatora | 5–8 | 0,91 ± 0,14 / 1,14 ± 0,18 | **0,90 ± 0,16 / 0,80 ± 0,14** | n.m. | 1,26 ± 0,14 / 1,12 ± 0,13 | 29 / 34; 24 / 31 |
| bez klasyfikatora | 8–30 | 1,11 ± 0,14 / 0,98 ± 0,16 | **1,16 ± 0,17 / 1,02 ± 0,17** | n.m. | 1,09 ± 0,06 / 0,96 ± 0,09 | 28 / —; — |

- Błąd stosunków względem mb jest zdominowany przez statystykę mb_nd (5–13 % na bin) i jest wspólny dla obu runów.
- Normalizacja względem mb: wszystkie 16 komórek ≥ 0,10 mieści się w 1,5σ od 1; największe odchylenie przy WP70 to 1,17 ± 0,12.
- Normalizacja względem fd: dane są wyżej o 20–50 % przy WP70 poniżej 8 GeV i ~2 razy bez klasyfikatora; zgodność tylko w 8–30 GeV.
- Kształt: względem mb χ² = 13–26 na 24 (zgodny) poza 5–8 GeV bez klasyfikatora (31–34); względem fd w 428770 34–59 (p od 0,08 do < 10⁻⁴).
  Plaster względem ≥ 0,25 w tych samych danych: 21 i 14 (WP70), 22 i 21 (bez klasyfikatora) — próg nie zmienia kształtu w danych.
  W fd przy WP70 kształt plastra różni się od ≥ 0,25 (86/24), bez klasyfikatora nie (21/24); kształt mb względem fd w 3–30 GeV: 23–27/24 (w binie 5–8 GeV 30–51 przy 200–370 odbiciach w mb).

### 7.5 Co się zmieniło względem wersji wstępnej
| wielkość | wstępnie (MC z ≥ 0,25, 427929 z 226 plików) | teraz |
|---|---|---|
| δμ(Δm), 3–5 GeV, po dodaniu plastra | −127 ± 53 keV (tylko 427929) | −51 ± 35 (427929), +50 ± 40 (428770); razem − MC: −6 ± 26 keV |
| δμ(Δm), 3–30 GeV | −75 ± 32 keV | −31 ± 21, +16 ± 14; razem − MC: +2 ± 12 keV |
| Δμ(Δm) 427929 przy ≥ 0,25, 3–5 / 3–30 GeV | +60 ± 102 / +42 ± 64 | +21 ± 69 / +21 ± 44 |
| k(≥0,10)/k(≥0,25), Δm, 3–30 | 0,96 ± 0,05 (427929) | 0,95 ± 0,03 (427929), 1,01 ± 0,02 (428770) |
| δμ m(Kπ), 3–30, po dodaniu plastra | −1,03 ± 0,92 MeV (427929) | −0,19 ± 0,50, +0,31 ± 0,34; razem − MC: +0,17 ± 0,28 MeV |
| wpływ zamiany MC ≥ 0,25 → MC przy tym samym cięciu | założone 0 | Δm: ≤ 2 keV; m(Kπ): ≤ 0,03 MeV przy WP70, 0,33 MeV bez klasyfikatora |
| R plaster / R ≥ 0,25 | 1,47 ± 0,31 w danych, MC nieznane | dane 1,5–2,3; mb 1,7–2,5; fd 1,1–1,3 |
| R dane/MC przy 0,10 | brak | względem mb 0,97 ± 0,07 i 1,11 ± 0,08 (WP70, 3–30) |
| 428770 przy 0,10 | bez score, pT(D*), mierzalne tylko 8–30 | WP70, pT(D⁰), wszystkie biny Δm |

Założenie wersji wstępnej („MC nie zmienia się poniżej 0,25”) sprawdziło się dla Δm i dla m(Kπ) przy WP70; nie sprawdza się dla m(Kπ) bez
klasyfikatora (−0,33 MeV) i dla R.

### 7.6 Propozycje rozmiarów systematyk — wersja końcowa dla 3–30 GeV
| systematyka | wstępnie (§5) | **teraz** | uzasadnienie liczbowe |
|---|---|---|---|
| PEAK_POSITION, m(Kπ) D⁰ | ±0,8 MeV | **±0,8 MeV** | dane − fd przy ≥ 0,10: +0,00 ± 0,68 MeV (biny −0,60 ± 1,20, +0,12 ± 1,34, −0,52 ± 0,93); efekt progu po odjęciu MC +0,17 ± 0,28. Liczba ta sama, ale teraz z pomiaru (0,68 ⊕ 0,28 ⊕ 0,33 z zależności μ_MC od klasyfikatora w plasterku), nie z ekstrapolacji. Szablon MC musi być z tego samego cięcia L_xy. Jeśli kształt bierze z minbias zamiast forced-decay: mb − fd = +1,3 ± 0,9 MeV, wtedy ±1,3 MeV. |
| PEAK_POSITION, Δm 3–5 GeV | ±0,3 MeV albo wolne μ | **±0,2 MeV** (albo wolne μ) | oba runy ≥ 0,10: −89 ± 49 keV; najgorszy run −161 ± 73; rozrzut runów przy ≥ 0,25 ±116 keV. Dodatkowe −127 keV z progu odpada (−6 ± 26). |
| PEAK_POSITION, Δm 5–8 GeV | ±0,15 MeV | **±0,15 MeV** | −71 ± 36 keV; najgorszy run −108 ± 46; efekt progu +1 ± 8 |
| PEAK_POSITION, Δm 8–30 GeV | ±0,10 MeV | **±0,10 MeV** | +19 ± 30 keV; efekt progu +1 ± 2; bez klasyfikatora −1 ± 39 |
| PEAK_WIDTH, Δm | ±10 % | **±10 %** | k przy ≥ 0,10: 1,01 ± 0,08 / 1,04 ± 0,05 / 1,08 ± 0,03; zmiana z progiem 1,01 ± 0,02 i 0,95 ± 0,03; 427929 w 8–30 GeV ma k = 1,24 ± 0,09 (2,7σ) — dlatego nie zwężać |
| PEAK_WIDTH, m(Kπ) | ±10 % | **±10 %** | 0,96 ± 0,11 / 0,98 ± 0,11 / 0,92 ± 0,06; 3–30: 1,02 ± 0,06; MC: plaster węższy o 4–7 %, całość o ≤ 1,2 % |
| PEAK_POSITION / WIDTH, m(Kππ) D⁺ | ±0,5 MeV / ±10 % (±6 %) | bez zmian | selekcja D⁺ ta sama w v1.1 i v2.0; nie mierzono ponownie |
| Reflection, normalizacja | ±10 % (≥ 0,25), ±20 % do czasu MC | **±15 % przy R z MC minbias** (to samo cięcie L_xy) | dane/mb przy ≥ 0,10: 0,97 ± 0,07 i 1,11 ± 0,08; biny 0,88–1,17 przy błędach 0,11–0,13; plaster 0,75 ± 0,23 (WP70), 1,10 ± 0,14 i 0,85 ± 0,11 (bez klasyfikatora). Dokładność ogranicza statystyka mb_nd. **R z forced-decay nie używać**: dane/fd = 1,2–1,5 przy WP70, ~2 bez klasyfikatora; jeśli inaczej się nie da, skala ×1,3 ± 0,2 (WP70, < 8 GeV) |
| Reflection, kształt | wariant mb zamiast fd | **wariant: kształt mb zamiast fd** (bez zmian) | względem mb χ²/24 = 13–26; względem fd 34–59 w 428770; próg L_xy kształtu w danych nie zmienia (14–22/24) |

**Nadal ekstrapolacja:**
- **CharmNet.** Nasz klasyfikator (trenowany przy L_xy ≥ 0,25, L_xy jest cechą) zostawia w plasterku 3–4 % sygnału, więc liczby WP70 przy ≥ 0,10 to w 96 %
  ta sama próbka co ≥ 0,25. Selekcja, która plaster zachowuje, leży między WP70 a „bez klasyfikatora”. Bez klasyfikatora zmierzone są: Δm w 5–8 i 8–30 GeV,
  m(Kπ) w 8–30 GeV i R we wszystkich binach; **Δm i m(Kπ) w 3–5 GeV bez klasyfikatora — nie**. Tam propozycje ±0,2 MeV i ±0,8 MeV opierają się na WP70
  i na tym, że MC z prawdy nie widzi w plasterku przesunięcia Δm (−1 ± 19 keV przy WP70, +4 ± 10 bez) — to argument z MC, nie pomiar w danych.
  m(Kπ) w plasterku jest w MC o 1,1–1,4 MeV niżej bez klasyfikatora; danych, które by to potwierdziły, nie ma (plaster niemierzalny).
- **2–3 GeV.** Bez zmian: brak kandydatów w danych i w MC (tylko produkcja `dpt3to150`). Wszystkie liczby dotyczą 3–30 GeV; trend z pT (Δm: −89 → −71 → +19 keV;
  R dane/fd rośnie ku niskim pT) sugeruje, że w 2–3 GeV będzie gorzej niż w 3–5, ale to nie jest pomiar.
- **D⁺.** Nie mierzone przy v2.0 (selekcja identyczna); liczby z PEAK_pos.

### 7.7 Zastrzeżenia
1. Rozbieżność runów: Δm przy ≥ 0,25 w 3–5 GeV (2,5σ) oraz R bez klasyfikatora (≥ 0,25, 3–30 GeV: 1,05 ± 0,05 wobec 1,41 ± 0,10, 3,3σ; R plaster/R ≥ 0,25:
   2,05 ± 0,18 wobec 1,18 ± 0,14). Runy mają inne widmo pT i tryger (428770 zdominowany L1_J12, 427929 MBTS), a metoda odejmowania sidebandów przy S/B
   bez klasyfikatora ma własną systematykę, której nie oszacowałem. Stosunki R bez klasyfikatora traktować z dokładnością ~20 %.
2. 427929: wszystkie zdarzenia, bez GRL i bez wymogu trygera (kolumn nie ma w wyciągu).
3. mb_nd: 78 plików, 2780 prawdziwych D⁰ z tagiem. Wariant „mb” ma błędy 30–75 keV (Δm), 0,7–1,2 MeV (m(Kπ)) i 5–13 % (R); w plasterku przy WP70 27 %.
   Różnica m(Kπ) mb − fd = +1,3 ± 0,9 MeV nie jest rozstrzygnięta.
4. Kształt MC w plasterku przy WP70 w 8–30 GeV pochodzi z 22 zdarzeń (fit niezbieżny) — komórka oznaczona n.m.; nie wpływa na ≥ 0,10.
5. Błędy Δμ w tabelach to Hesse ⊕ błąd μ_MC; błędy różnic i stosunków po dodaniu plastra oraz błędy R — bootstrap po zdarzeniach (200 replik,
   połowa przedziału 16–84 %). Rozrzut od modelu tła (PEAK_pos: 11–34 keV, 0,3–0,8 MeV) nadal nie jest dodany; test na pseudo-danych mb (+76 ± 39 keV
   przy zerowym efekcie) pokazuje, że w 3–5 GeV jest on tego rzędu co błąd statystyczny.
6. Sygnał mb z prawdy w m(Kπ) dopasowany swobodnie (`v2/mb_truth.py`, n ≥ 300); fit z tłem Czebyszewa na czystym sygnale był zdegenerowany
   i tych wierszy (`rows_mb_nd_truth.csv`, obs = mKpi) nie używać.
7. D⁰ tylko przez tag D*; średnie „oba runy” nie uwzględniają wspólnego błędu MC (mały dla fd, dominujący dla mb).

### Pliki v2
`v2/run_v2.sh`, `v2/check_inputs.py`, `v2/fit_v2.py`, `v2/refl_v2.py`, `v2/mb_truth.py`, `v2/summarize_v2.py`;
`v2/out/tables_v2.md`, `summary_v2.json`, `check_inputs.{txt,json}`, `mc_fd.csv`, `mc_fd_boot.csv`, `rows_{data_428770,data_427929,mb_nd,mb_nd_truth}.csv`,
`boot_{data_428770,data_427929,mb_nd}.csv`, `mb_truth_free.csv`, `refl_rows.csv`, `refl_boot.csv`, `refl_shape.csv`; `v2/logs/`.
