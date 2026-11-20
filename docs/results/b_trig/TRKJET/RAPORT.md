# TRKJET — czy deficyt izolacji D⁰ w trygerze dżetowym to „D⁰ w dżecie”? (projekt B)

Wygenerowane: `analyze.py`, wejście `runs/lxplus`, 2026-11-17 05:37. Liczby: `out/trkjet.json`.

## 0. Werdykt

**H ODRZUCONA — przeważenie w zmiennych dżetowych (bez koła) nie domyka J12only lepiej niż pT×lead_pt×nSelBL.**

| kryterium | co sprawdza | liczby | wynik |
|---|---|---|---|
| K1 (O5) | ułamek sygnału: dżet spoza stożka (kol. O) pT > 5 GeV w ΔR < 0,8; J12only > 9:MB o > 3σ | J12only 0,169 ± 0,016, 9:MB 0,055 ± 0,014, MC 0,052 ± 0,001; różnica 5,4σ | TAK |
| K1 (A10 [koło]) | ułamek sygnału: dżet zawierający D⁰ (kol. A) pT > 10 GeV; J12only > 9:MB o > 3σ | J12only 0,611 ± 0,017, 9:MB 0,162 ± 0,016, MC 0,184 ± 0,002; różnica 19,2σ | TAK |
| K2 (O5) | SF_iso osobno w/poza dżetem zgodny z 1 (|pull| < 2) w 9:MB i J12only | w|9:MB: 0,43 ± 0,26 (słaby); w|8:J12only: 0,51 ± 0,11; poza|9:MB: 1,00 ± 0,04; poza|8:J12only: 0,47 ± 0,03 | NIE (miarodajnych 4/4) |
| K2 (O10) | SF_iso osobno w/poza dżetem zgodny z 1 (|pull| < 2) w 9:MB i J12only | w|9:MB: 1,47 ± 2,68 (niemiarodajny); w|8:J12only: 0,50 ± 0,14; poza|9:MB: 0,98 ± 0,04; poza|8:J12only: 0,46 ± 0,03 | NIE (miarodajnych 3/4) |
| K3 | przeważenie MC bez koła domyka 9:MB i J12only jednocześnie (|SF−1| < 2σ) | najlepszy: „pT×Σ pT półkuli przeciwnej”: 9:MB 0,988 ± 0,034, J12only 0,704 ± 0,036, χ²/ndf 80,7/4; odniesienie pT×lead_pt×nSelBL: J12only 0,761 ± 0,040, χ² 48,0/4 | NIE |
| K3 sufit [koło] | to samo ze zmienną z torów stożka | „pT×pT dżetu z D0 (kol. A) [koło]”: 9:MB 0,977 ± 0,023, J12only 0,859 ± 0,035, χ² 25,1/4 | (informacyjnie) |
| K4 | ε_iso(dane) w funkcji ΔR/pT dżetu kol. O: 9:MB vs J12only, p > 0,01 | o_drn5: χ²/ndf 156,6/8 (p = 8.24e-30); o_drn3: χ²/ndf 180,7/8 (p = 7.26e-35); o_ptmax08: χ²/ndf 190,2/4 (p = 4.8e-40) | NIE (testów z ndf ≥ 3: 3) |

Jak czytać: H wymaga K3 i K4 (przy tej samej topologii dżetowej MC opisuje izolację i dane różnych klas są zgodne). K1 sam nie wystarcza (tryger dżetowy z definicji wybiera zdarzenia z dżetem). Jeżeli K3 = NIE przy K1 = TAK: D⁰ rzeczywiście leży częściej przy dżecie, ale to nie tłumaczy deficytu — różnica siedzi w samym stożku (fragmentacja wokół D⁰ w danych ≠ MC przy tym samym otoczeniu).

## 1. Metoda i zastrzeżenia

- Pary D⁰←D* i score H wzięte z `b_trig` (`bt_<próbka>.npz`), bez ponownego liczenia. Selekcja: |m(Kπ) − m₀| < 40 MeV, score_H > 0.6, pT 3–30 GeV, komplet córek wycięty. 427929: MBTS_1 ∧ GRL. Cięcie izolacji: ptcone40 < 2.0 GeV (gałąź ntupla).
- Dżety: anti-kT (własna implementacja, schemat E), R = 0.4, tory pT > 0.5 GeV, |η| < 2.5, masa pionu. **fastjet nie jest dostępny** w `q_env` (LCG_110_swan) → własna implementacja w numba (`tj_kernel.py`, test jednostkowy w pliku).
- Kolekcje: **A** — D⁰ jako jedna pseudocząstka zamiast córek (dżet z D⁰, z); **B** — tory bez córek D⁰ i wolnego pionu; **O** — jak B, ale tylko tory z ΔR(tor, D⁰) > 0,4. Zmienne A i B zawierają tory stożka izolacji → związek z ptcone40 z konstrukcji (**[koło]**); test hipotezy opiera się na O, półkuli przeciwnej i pierścieniu 0,4–1,0.
- ε i ułamki: odejmowanie pasm Δm (SBS, `ISO_origin/iso_common.sbs_w`); fit jednoczesny pass/fail (`charmtag.dstar_refit.fit_passfail`, HESSE) przy N_sig(SBS) ≥ 150. Oznaczenia: „słaby” N_sig < 150, „niemiarodajny” N_sig < 50 albo N_sig/σ < 3.
- MC: fd_nd, truth ∈ {0,1}, Δm 143,5–147,5 MeV, ta sama selekcja; N = 41236, ε_MC = 0,664. mb_nd tylko w (a).
- Przeważenie: SF = Σw·pass / Σw·ε_MC(komórka) (`s2_universal_rw.sf_rw`), komórki z n_MC < 25 cofają się do mniejszej liczby wymiarów (kolumna „cof.” = udział takiego sygnału). **Zastrzeżenie:** fd_nd ma mało D⁰ przy twardych dżetach — duże „cof.” oznacza, że schemat nie jest testem.

**Zgodność własnego ptcone z gałęzią `ptcone40`** (Σ pT wszystkich torów drzewa `Tracks` w ΔR < 0,4 bez córek D⁰, wolny pion wliczony):

| próbka | par | komplet córek | identyczne (< 1 MeV) | |Δ| < 0,2 GeV | własne ≤ gałąź | ta sama decyzja (< 2 GeV) |
|---|---|---|---|---|---|---|
| data_427929 | 402470 | 0,9902 | 0,649 | 0,895 | 1,0000 | 0,9916 |
| data_428770 | 1233580 | 0,9917 | 0,719 | 0,955 | 1,0000 | 0,9950 |
| fd_nd | 178721 | 0,9952 | 0,826 | 0,985 | 1,0000 | 0,9947 |
| mb_nd | 207021 | 0,9922 | 0,775 | 0,978 | 1,0000 | 0,9930 |

Różnica jest jednostronna (własne ≤ gałąź, kwanty 0,1–0,2 GeV): gałąź liczy tory od pT > 0,1 GeV, drzewo `Tracks` ma pT > 0,2 GeV. Gałąź zawiera wolny pion D* (wersja bez pionu zgadza się w < 1 %).

## 2. Odtworzenie punktu wyjścia — ε_iso i SF bez wag per klasa

| klasa | N_sig (SBS) | ε SBS | ε fit | ε_MC | SF | jakość |
|---|---|---|---|---|---|---|
| 9:MB | 957 ± 46 | 0,655 ± 0,024 | 0,649 ± 0,028 | 0,664 | 0,98 ± 0,04 (fit) | ok |
| 9:MB&!j20&!J12 | 717 ± 39 | 0,752 ± 0,026 | 0,721 ± 0,029 | 0,664 | 1,09 ± 0,04 (fit) | ok |
| 9:MB&j20&!J12 | 210 ± 22 | 0,383 ± 0,051 | 0,440 ± 0,063 | 0,664 | 0,66 ± 0,09 (fit) | ok |
| 8:all | 2796 ± 83 | 0,324 ± 0,013 | 0,338 ± 0,016 | 0,664 | 0,51 ± 0,02 (fit) | ok |
| 8:J12only | 2194 ± 73 | 0,293 ± 0,015 | 0,302 ± 0,018 | 0,664 | 0,46 ± 0,03 (fit) | ok |
| 8:j20&!J12 | 294 ± 27 | 0,396 ± 0,044 | 0,412 ± 0,052 | 0,664 | 0,62 ± 0,08 (fit) | ok |
| 8:MBTSbased | 369 ± 30 | 0,404 ± 0,040 | 0,415 ± 0,045 | 0,664 | 0,63 ± 0,07 (fit) | ok |
| 8:other | 238 ± 23 | 0,480 ± 0,048 | 0,532 ± 0,060 | 0,664 | 0,80 ± 0,09 (fit) | ok |

Oczekiwane z ISO_origin (pełna statystyka): 9:MB 0,98 ± 0,04; 8:J12only 0,46 ± 0,03; 8:MBTSbased 0,62 ± 0,07; 8:other 0,80 ± 0,09.

## 3. (a) Ułamek sygnału D⁰ w dżecie / przy dżecie

| własność | 9:MB | 9:MB&!j20&!J12 | 9:MB&j20&!J12 | 8:all | 8:J12only | 8:j20&!J12 | 8:MBTSbased | 8:other | MC fd_nd | MC mb_nd |
|---|---|---|---|---|---|---|---|---|---|---|
| **A5** [koło]: dżet zawierający D0 (kol. A) pT > 5 GeV i ΔR(D0, oś) < 0,4 | 0,602 ± 0,024 | 0,519 ± 0,027 | 0,843 ± 0,049 | 0,842 ± 0,012 | 0,864 ± 0,013 | 0,764 ± 0,042 | 0,732 ± 0,038 | 0,817 ± 0,046 | 0,639 ± 0,002 | 0,619 ± 0,015 |
| **A10** [koło]: jw., pT > 10 GeV | 0,162 ± 0,016 | 0,070 ± 0,013 | 0,379 ± 0,048 | 0,562 ± 0,015 | 0,611 ± 0,017 | 0,398 ± 0,043 | 0,399 ± 0,039 | 0,366 ± 0,044 | 0,184 ± 0,002 | 0,192 ± 0,012 |
| **A20** [koło]: jw., pT > 20 GeV | 0,026 ± 0,007 | 0,000 ± 0,000 | 0,044 ± 0,021 | 0,293 ± 0,013 | 0,354 ± 0,015 | 0,058 ± 0,020 | 0,088 ± 0,022 | 0,042 ± 0,017 | 0,027 ± 0,001 | 0,022 ± 0,005 |
| **z<0.5** [koło]: z = pT(D0)/pT(dżetu z D0) < 0,5 | 0,061 ± 0,017 | 0,028 ± 0,016 | 0,131 ± 0,046 | 0,194 ± 0,015 | 0,221 ± 0,018 | 0,083 ± 0,042 | 0,078 ± 0,039 | 0,128 ± 0,041 | 0,065 ± 0,001 | 0,062 ± 0,008 |
| **B2** [koło]: ΔR(D0, dżet otoczenia kol. B, pT > 2 GeV) < 0,4 | 0,231 ± 0,022 | 0,152 ± 0,023 | 0,437 ± 0,052 | 0,534 ± 0,015 | 0,573 ± 0,016 | 0,394 ± 0,046 | 0,397 ± 0,042 | 0,394 ± 0,048 | 0,237 ± 0,002 | 0,210 ± 0,013 |
| **B3** [koło]: jw., pT > 3 GeV | 0,148 ± 0,019 | 0,073 ± 0,019 | 0,339 ± 0,050 | 0,433 ± 0,015 | 0,481 ± 0,017 | 0,288 ± 0,044 | 0,279 ± 0,040 | 0,235 ± 0,045 | 0,153 ± 0,002 | 0,142 ± 0,011 |
| **B5** [koło]: jw., pT > 5 GeV | 0,078 ± 0,015 | 0,018 ± 0,012 | 0,214 ± 0,044 | 0,322 ± 0,015 | 0,378 ± 0,017 | 0,135 ± 0,038 | 0,131 ± 0,036 | 0,099 ± 0,037 | 0,074 ± 0,001 | 0,079 ± 0,009 |
| **B8** [koło]: jw., pT > 8 GeV | 0,027 ± 0,010 | 0,003 ± 0,006 | 0,061 ± 0,031 | 0,213 ± 0,014 | 0,261 ± 0,016 | 0,025 ± 0,029 | 0,039 ± 0,028 | 0,044 ± 0,027 | 0,032 ± 0,001 | 0,029 ± 0,005 |
| **O3**: dżet spoza stożka (kol. O) pT > 3 GeV w ΔR < 0,8 | 0,142 ± 0,020 | 0,102 ± 0,021 | 0,266 ± 0,051 | 0,313 ± 0,015 | 0,324 ± 0,017 | 0,305 ± 0,045 | 0,290 ± 0,041 | 0,244 ± 0,046 | 0,148 ± 0,002 | 0,167 ± 0,012 |
| **O5**: jw., pT > 5 GeV | 0,055 ± 0,014 | 0,038 ± 0,013 | 0,136 ± 0,040 | 0,160 ± 0,013 | 0,169 ± 0,016 | 0,150 ± 0,036 | 0,147 ± 0,033 | 0,102 ± 0,035 | 0,052 ± 0,001 | 0,064 ± 0,008 |
| **O10**: jw., pT > 10 GeV | 0,003 ± 0,006 | -0,002 ± 0,002 | 0,011 ± 0,019 | 0,060 ± 0,009 | 0,069 ± 0,011 | 0,038 ± 0,020 | 0,042 ± 0,019 | 0,006 ± 0,015 | 0,009 ± 0,000 | 0,012 ± 0,003 |
| **away10**: Σ pT torów w półkuli przeciwnej > 10 GeV | 0,761 ± 0,017 | 0,703 ± 0,022 | 0,943 ± 0,019 | 0,957 ± 0,005 | 0,966 ± 0,005 | 0,967 ± 0,015 | 0,957 ± 0,014 | 0,874 ± 0,026 | 0,771 ± 0,002 | 0,760 ± 0,014 |

\* N_sig < 150 albo niemiarodajny. Ułamki są dla sygnału (wagi SBS), nie dla wszystkich kandydatów.

## 4. (b) SF izolacji osobno „w dżecie” i „poza dżetem”

SF = ε_dane / ε_MC w tej samej kategorii (MC bez wag). Dla definicji [koło] ε „w dżecie” jest małe z konstrukcji — patrz tylko na „poza”.

| podział | strona | klasa | N_sig | ε dane | metoda | ε MC (N) | SF | jakość |
|---|---|---|---|---|---|---|---|---|
| A10 [koło] | w | 9:MB | 155 ± 16 | 0,162 ± 0,034 | fit | 0,186 (7577) | 0,87 ± 0,18 | ok |
| A10 [koło] | w | 9:MB&!j20&!J12 | 51 ± 9 | 0,194 ± 0,061 | SBS | 0,186 (7577) | 1,04 ± 0,33 | słaby |
| A10 [koło] | w | 9:MB&j20&!J12 | 80 ± 12 | 0,126 ± 0,040 | SBS | 0,186 (7577) | 0,68 ± 0,22 | słaby |
| A10 [koło] | w | 8:all | 1570 ± 57 | 0,100 ± 0,009 | fit | 0,186 (7577) | 0,54 ± 0,05 | ok |
| A10 [koło] | w | 8:J12only | 1339 ± 53 | 0,081 ± 0,009 | fit | 0,186 (7577) | 0,43 ± 0,05 | ok |
| A10 [koło] | w | 8:j20&!J12 | 117 ± 15 | 0,241 ± 0,045 | SBS | 0,186 (7577) | 1,29 ± 0,24 | słaby |
| A10 [koło] | w | 8:MBTSbased | 147 ± 17 | 0,206 ± 0,038 | SBS | 0,186 (7577) | 1,11 ± 0,21 | słaby |
| A10 [koło] | w | 8:other | 87 ± 12 | 0,183 ± 0,044 | SBS | 0,186 (7577) | 0,98 ± 0,24 | słaby |
| A10 [koło] | poza | 9:MB | 803 ± 42 | 0,751 ± 0,031 | fit | 0,771 (33659) | 0,97 ± 0,04 | ok |
| A10 [koło] | poza | 9:MB&!j20&!J12 | 667 ± 38 | 0,766 ± 0,030 | fit | 0,771 (33659) | 0,99 ± 0,04 | ok |
| A10 [koło] | poza | 9:MB&j20&!J12 | 131 ± 19 | 0,538 ± 0,071 | SBS | 0,771 (33659) | 0,70 ± 0,09 | słaby |
| A10 [koło] | poza | 8:all | 1225 ± 59 | 0,628 ± 0,031 | fit | 0,771 (33659) | 0,81 ± 0,04 | ok |
| A10 [koło] | poza | 8:J12only | 854 ± 50 | 0,626 ± 0,038 | fit | 0,771 (33659) | 0,81 ± 0,05 | ok |
| A10 [koło] | poza | 8:j20&!J12 | 177 ± 22 | 0,546 ± 0,077 | fit | 0,771 (33659) | 0,71 ± 0,10 | ok |
| A10 [koło] | poza | 8:MBTSbased | 222 ± 25 | 0,567 ± 0,067 | fit | 0,771 (33659) | 0,74 ± 0,09 | ok |
| A10 [koło] | poza | 8:other | 151 ± 20 | 0,717 ± 0,082 | fit | 0,771 (33659) | 0,93 ± 0,11 | ok |
| B3 [koło] | w | 9:MB | 141 ± 21 | -0,029 ± 0,015 | SBS | 0,014 (6309) | -2,01 ± 1,08 | słaby |
| B3 [koło] | w | 9:MB&!j20&!J12 | 53 ± 14 | -0,049 ± 0,034 | SBS | 0,014 (6309) | -3,46 ± 2,44 | słaby |
| B3 [koło] | w | 9:MB&j20&!J12 | 71 ± 13 | -0,024 ± 0,018 | SBS | 0,014 (6309) | -1,65 ± 1,31 | słaby |
| B3 [koło] | w | 8:all | 1212 ± 58 | 0,001 ± 0,003 | fit | 0,014 (6309) | 0,04 ± 0,24 | ok |
| B3 [koło] | w | 8:J12only | 1055 ± 54 | 0,002 ± 0,003 | fit | 0,014 (6309) | 0,14 ± 0,23 | ok |
| B3 [koło] | w | 8:j20&!J12 | 85 ± 16 | -0,007 ± 0,013 | SBS | 0,014 (6309) | -0,50 ± 0,89 | słaby |
| B3 [koło] | w | 8:MBTSbased | 103 ± 18 | -0,007 ± 0,010 | SBS | 0,014 (6309) | -0,50 ± 0,74 | słaby |
| B3 [koło] | w | 8:other | 56 ± 13 | -0,000 ± 0,000 | SBS | 0,014 (6309) | -0,00 ± 0,00 | słaby |
| B3 [koło] | poza | 9:MB | 815 ± 41 | 0,777 ± 0,027 | fit | 0,781 (34927) | 0,99 ± 0,03 | ok |
| B3 [koło] | poza | 9:MB&!j20&!J12 | 665 ± 36 | 0,805 ± 0,028 | fit | 0,781 (34927) | 1,03 ± 0,04 | ok |
| B3 [koło] | poza | 9:MB&j20&!J12 | 139 ± 18 | 0,587 ± 0,061 | SBS | 0,781 (34927) | 0,75 ± 0,08 | słaby |
| B3 [koło] | poza | 8:all | 1589 ± 58 | 0,593 ± 0,022 | fit | 0,781 (34927) | 0,76 ± 0,03 | ok |
| B3 [koło] | poza | 8:J12only | 1143 ± 50 | 0,568 ± 0,027 | fit | 0,781 (34927) | 0,73 ± 0,03 | ok |
| B3 [koło] | poza | 8:j20&!J12 | 209 ± 21 | 0,622 ± 0,060 | fit | 0,781 (34927) | 0,80 ± 0,08 | ok |
| B3 [koło] | poza | 8:MBTSbased | 266 ± 24 | 0,607 ± 0,051 | fit | 0,781 (34927) | 0,78 ± 0,07 | ok |
| B3 [koło] | poza | 8:other | 182 ± 19 | 0,737 ± 0,065 | fit | 0,781 (34927) | 0,94 ± 0,08 | ok |
| O5 | w | 9:MB | 53 ± 14 | 0,203 ± 0,124 | SBS | 0,475 (2164) | 0,43 ± 0,26 | słaby |
| O5 | w | 9:MB&!j20&!J12 | 26 ± 9 | 0,388 ± 0,175 | SBS | 0,475 (2164) | 0,82 ± 0,37 | niemiarodajny |
| O5 | w | 9:MB&j20&!J12 | 30 ± 9 | 0,102 ± 0,148 | SBS | 0,475 (2164) | 0,21 ± 0,31 | niemiarodajny |
| O5 | w | 8:all | 448 ± 42 | 0,269 ± 0,049 | fit | 0,475 (2164) | 0,57 ± 0,10 | ok |
| O5 | w | 8:J12only | 370 ± 39 | 0,241 ± 0,052 | fit | 0,475 (2164) | 0,51 ± 0,11 | ok |
| O5 | w | 8:j20&!J12 | 45 ± 12 | 0,412 ± 0,123 | SBS | 0,475 (2164) | 0,87 ± 0,26 | niemiarodajny |
| O5 | w | 8:MBTSbased | 55 ± 14 | 0,411 ± 0,117 | SBS | 0,475 (2164) | 0,87 ± 0,25 | słaby |
| O5 | w | 8:other | 24 ± 9 | 0,196 ± 0,173 | SBS | 0,475 (2164) | 0,41 ± 0,36 | niemiarodajny |
| O5 | poza | 9:MB | 903 ± 43 | 0,674 ± 0,028 | fit | 0,674 (39072) | 1,00 ± 0,04 | ok |
| O5 | poza | 9:MB&!j20&!J12 | 691 ± 38 | 0,739 ± 0,028 | fit | 0,674 (39072) | 1,10 ± 0,04 | ok |
| O5 | poza | 9:MB&j20&!J12 | 181 ± 20 | 0,494 ± 0,065 | fit | 0,674 (39072) | 0,73 ± 0,10 | ok |
| O5 | poza | 8:all | 2349 ± 71 | 0,351 ± 0,017 | fit | 0,674 (39072) | 0,52 ± 0,03 | ok |
| O5 | poza | 8:J12only | 1825 ± 62 | 0,314 ± 0,019 | fit | 0,674 (39072) | 0,47 ± 0,03 | ok |
| O5 | poza | 8:j20&!J12 | 250 ± 24 | 0,416 ± 0,053 | fit | 0,674 (39072) | 0,62 ± 0,08 | ok |
| O5 | poza | 8:MBTSbased | 315 ± 27 | 0,419 ± 0,045 | fit | 0,674 (39072) | 0,62 ± 0,07 | ok |
| O5 | poza | 8:other | 213 ± 21 | 0,539 ± 0,063 | fit | 0,674 (39072) | 0,80 ± 0,09 | ok |
| O10 | w | 9:MB | 3 ± 6 | 0,598 ± 1,092 | SBS | 0,407 (361) | 1,47 ± 2,68 | niemiarodajny |
| O10 | w | 9:MB&!j20&!J12 | -2 ± 2 | – | SBS | 0,407 (361) | – | niemiarodajny |
| O10 | w | 9:MB&j20&!J12 | 2 ± 4 | 2,221 ± 3,010 | SBS | 0,407 (361) | 5,45 ± 7,40 | niemiarodajny |
| O10 | w | 8:all | 168 ± 27 | 0,210 ± 0,056 | fit | 0,407 (361) | 0,52 ± 0,14 | ok |
| O10 | w | 8:J12only | 153 ± 26 | 0,205 ± 0,056 | fit | 0,407 (361) | 0,50 ± 0,14 | ok |
| O10 | w | 8:j20&!J12 | 12 ± 6 | 0,604 ± 0,265 | SBS | 0,407 (361) | 1,48 ± 0,66 | niemiarodajny |
| O10 | w | 8:MBTSbased | 16 ± 7 | 0,375 ± 0,209 | SBS | 0,407 (361) | 0,92 ± 0,52 | niemiarodajny |
| O10 | w | 8:other | 0 ± 4 | -7,827 ± 111,119 | SBS | 0,407 (361) | -19,22 ± 272,89 | niemiarodajny |
| O10 | poza | 9:MB | 952 ± 45 | 0,650 ± 0,028 | fit | 0,666 (40875) | 0,98 ± 0,04 | ok |
| O10 | poza | 9:MB&!j20&!J12 | 719 ± 39 | 0,722 ± 0,029 | fit | 0,666 (40875) | 1,08 ± 0,04 | ok |
| O10 | poza | 9:MB&j20&!J12 | 206 ± 22 | 0,424 ± 0,062 | fit | 0,666 (40875) | 0,64 ± 0,09 | ok |
| O10 | poza | 8:all | 2631 ± 78 | 0,343 ± 0,017 | fit | 0,666 (40875) | 0,52 ± 0,03 | ok |
| O10 | poza | 8:J12only | 2044 ± 69 | 0,305 ± 0,019 | fit | 0,666 (40875) | 0,46 ± 0,03 | ok |
| O10 | poza | 8:j20&!J12 | 282 ± 26 | 0,407 ± 0,052 | fit | 0,666 (40875) | 0,61 ± 0,08 | ok |
| O10 | poza | 8:MBTSbased | 354 ± 29 | 0,419 ± 0,046 | fit | 0,666 (40875) | 0,63 ± 0,07 | ok |
| O10 | poza | 8:other | 237 ± 23 | 0,532 ± 0,059 | fit | 0,666 (40875) | 0,80 ± 0,09 | ok |

## 5. (c) SF po przeważeniu MC w zmiennych dżetowych

| schemat | 9:MB | 9:MB&!j20&!J12 | 9:MB&j20&!J12 | 8:all | 8:J12only | 8:j20&!J12 | 8:MBTSbased | 8:other | χ²/ndf vs 1 (4 klasy rozłączne) | cof. J12only |
|---|---|---|---|---|---|---|---|---|---|---|
| brak | 0,986 ± 0,036 | 1,133 ± 0,040 | 0,577 ± 0,077 | 0,488 ± 0,020 | 0,441 ± 0,022 | 0,596 ± 0,067 | 0,609 ± 0,060 | 0,724 ± 0,073 | 703,5/4 | 0,00 |
| pT | 0,977 ± 0,036 | 1,077 ± 0,039 | 0,629 ± 0,084 | 0,632 ± 0,027 | 0,600 ± 0,031 | 0,664 ± 0,076 | 0,695 ± 0,070 | 0,776 ± 0,078 | 197,0/4 | 0,00 |
| pT×lead_pt×nSelBL (ISO_origin) | 0,985 ± 0,034 | 1,053 ± 0,035 | 0,712 ± 0,094 | 0,774 ± 0,033 | 0,761 ± 0,040 | 0,733 ± 0,082 | 0,786 ± 0,077 | 0,837 ± 0,081 | 48,0/4 | 0,00 |
| pT×O: pT dżetu w ΔR<0,8 | 0,976 ± 0,035 | 1,066 ± 0,037 | 0,653 ± 0,087 | 0,641 ± 0,026 | 0,606 ± 0,030 | 0,698 ± 0,078 | 0,718 ± 0,070 | 0,788 ± 0,077 | 196,7/4 | 0,03 |
| pT×O: (ΔR, pT) najbliższego dżetu pT>3 | 0,974 ± 0,035 | 1,059 ± 0,037 | 0,660 ± 0,088 | 0,661 ± 0,027 | 0,627 ± 0,031 | 0,701 ± 0,079 | 0,728 ± 0,072 | 0,808 ± 0,079 | 162,5/4 | 0,10 |
| pT×O: (ΔR, pT) najbliższego dżetu pT>5 | 0,978 ± 0,035 | 1,062 ± 0,037 | 0,668 ± 0,089 | 0,677 ± 0,028 | 0,647 ± 0,033 | 0,708 ± 0,080 | 0,737 ± 0,073 | 0,804 ± 0,080 | 134,4/4 | 0,09 |
| pT×Σ pT półkuli przeciwnej | 0,988 ± 0,034 | 1,062 ± 0,036 | 0,698 ± 0,093 | 0,729 ± 0,031 | 0,704 ± 0,036 | 0,741 ± 0,083 | 0,777 ± 0,076 | 0,820 ± 0,080 | 80,7/4 | 0,00 |
| pT×pierścień 0,4–1,0 (własny) | 0,983 ± 0,034 | 1,059 ± 0,036 | 0,691 ± 0,091 | 0,678 ± 0,028 | 0,647 ± 0,032 | 0,719 ± 0,080 | 0,747 ± 0,073 | 0,802 ± 0,077 | 138,4/4 | 0,01 |
| pT×O pT(ΔR<0,8)×nSelBL | 0,986 ± 0,034 | 1,061 ± 0,035 | 0,703 ± 0,093 | 0,706 ± 0,029 | 0,676 ± 0,034 | 0,740 ± 0,082 | 0,771 ± 0,075 | 0,816 ± 0,079 | 106,8/4 | 0,06 |
| pT×O pT(ΔR<0,8)×pierścień | 0,981 ± 0,034 | 1,062 ± 0,036 | 0,679 ± 0,089 | 0,674 ± 0,028 | 0,640 ± 0,032 | 0,711 ± 0,079 | 0,742 ± 0,072 | 0,809 ± 0,078 | 146,1/4 | 0,01 |
| pT×B: (ΔR, pT) najbliższego dżetu pT>3 [koło] | 0,991 ± 0,028 | 1,019 ± 0,029 | 0,840 ± 0,088 | 0,831 ± 0,026 | 0,831 ± 0,031 | 0,806 ± 0,073 | 0,809 ± 0,063 | 0,858 ± 0,067 | 43,1/4 | 0,10 |
| pT×pT dżetu z D0 (kol. A) [koło] | 0,977 ± 0,023 | 0,998 ± 0,023 | 0,876 ± 0,078 | 0,866 ± 0,028 | 0,859 ± 0,035 | 0,844 ± 0,070 | 0,835 ± 0,061 | 0,962 ± 0,067 | 25,1/4 | 0,01 |

\* N_sig < 150. Klasy rozłączne do χ²: 9:MB, 8:J12only, 8:MBTSbased, 8:other.

Reszta w binach pT(D⁰) (ISO_origin: J12only 0,97 / 0,72 / 0,50 po pT×lead_pt×nSelBL):

| schemat | klasa | 3–5 | 5–8 | 8–30 |
|---|---|---|---|---|
| brak | 9:MB | 0,94 ± 0,05 | 1,08 ± 0,05 | 0,89 ± 0,12 * |
| brak | 8:J12only | 0,73 ± 0,06 | 0,53 ± 0,04 | 0,30 ± 0,03 |
| brak | 8:MBTSbased | 0,69 ± 0,11 | 0,59 ± 0,10 * | 0,74 ± 0,11 * |
| pT | 9:MB | 0,94 ± 0,05 | 1,08 ± 0,05 | 0,91 ± 0,11 * |
| pT | 8:J12only | 0,74 ± 0,06 | 0,54 ± 0,04 | 0,43 ± 0,04 |
| pT | 8:MBTSbased | 0,69 ± 0,11 | 0,61 ± 0,11 * | 0,83 ± 0,12 * |
| pT×lead_pt×nSelBL (ISO_origin) | 9:MB | 0,95 ± 0,05 | 1,08 ± 0,05 | 0,91 ± 0,11 * |
| pT×lead_pt×nSelBL (ISO_origin) | 8:J12only | 0,97 ± 0,08 | 0,72 ± 0,06 | 0,50 ± 0,05 |
| pT×lead_pt×nSelBL (ISO_origin) | 8:MBTSbased | 0,81 ± 0,13 | 0,66 ± 0,11 * | 0,89 ± 0,13 * |
| pT×O: pT dżetu w ΔR<0,8 | 9:MB | 0,93 ± 0,05 | 1,08 ± 0,05 | 0,92 ± 0,11 * |
| pT×O: pT dżetu w ΔR<0,8 | 8:J12only | 0,78 ± 0,06 | 0,58 ± 0,05 | 0,39 ± 0,03 |
| pT×O: pT dżetu w ΔR<0,8 | 8:MBTSbased | 0,74 ± 0,12 | 0,61 ± 0,10 * | 0,81 ± 0,12 * |
| pT×O: (ΔR, pT) najbliższego dżetu pT>3 | 9:MB | 0,93 ± 0,05 | 1,08 ± 0,05 | 0,92 ± 0,12 * |
| pT×O: (ΔR, pT) najbliższego dżetu pT>3 | 8:J12only | 0,80 ± 0,06 | 0,61 ± 0,05 | 0,40 ± 0,04 |
| pT×O: (ΔR, pT) najbliższego dżetu pT>3 | 8:MBTSbased | 0,76 ± 0,12 | 0,61 ± 0,10 * | 0,81 ± 0,12 * |
| pT×O: (ΔR, pT) najbliższego dżetu pT>5 | 9:MB | 0,94 ± 0,05 | 1,08 ± 0,05 | 0,92 ± 0,12 * |
| pT×O: (ΔR, pT) najbliższego dżetu pT>5 | 8:J12only | 0,81 ± 0,07 | 0,64 ± 0,05 | 0,42 ± 0,04 |
| pT×O: (ΔR, pT) najbliższego dżetu pT>5 | 8:MBTSbased | 0,75 ± 0,12 | 0,63 ± 0,11 * | 0,86 ± 0,12 * |
| pT×Σ pT półkuli przeciwnej | 9:MB | 0,95 ± 0,05 | 1,09 ± 0,05 | 0,92 ± 0,11 * |
| pT×Σ pT półkuli przeciwnej | 8:J12only | 0,88 ± 0,07 | 0,66 ± 0,05 | 0,48 ± 0,04 |
| pT×Σ pT półkuli przeciwnej | 8:MBTSbased | 0,80 ± 0,13 | 0,65 ± 0,11 * | 0,90 ± 0,13 * |
| pT×pierścień 0,4–1,0 (własny) | 9:MB | 0,94 ± 0,05 | 1,09 ± 0,05 | 0,92 ± 0,11 * |
| pT×pierścień 0,4–1,0 (własny) | 8:J12only | 0,83 ± 0,07 | 0,61 ± 0,05 | 0,42 ± 0,04 |
| pT×pierścień 0,4–1,0 (własny) | 8:MBTSbased | 0,77 ± 0,12 | 0,62 ± 0,11 * | 0,86 ± 0,12 * |
| pT×O pT(ΔR<0,8)×nSelBL | 9:MB | 0,94 ± 0,05 | 1,09 ± 0,05 | 0,93 ± 0,11 * |
| pT×O pT(ΔR<0,8)×nSelBL | 8:J12only | 0,87 ± 0,07 | 0,65 ± 0,05 | 0,43 ± 0,04 |
| pT×O pT(ΔR<0,8)×nSelBL | 8:MBTSbased | 0,82 ± 0,13 | 0,62 ± 0,11 * | 0,87 ± 0,12 * |
| pT×O pT(ΔR<0,8)×pierścień | 9:MB | 0,94 ± 0,05 | 1,09 ± 0,05 | 0,92 ± 0,11 * |
| pT×O pT(ΔR<0,8)×pierścień | 8:J12only | 0,82 ± 0,07 | 0,61 ± 0,05 | 0,41 ± 0,04 |
| pT×O pT(ΔR<0,8)×pierścień | 8:MBTSbased | 0,77 ± 0,12 | 0,62 ± 0,10 * | 0,84 ± 0,12 * |
| pT×B: (ΔR, pT) najbliższego dżetu pT>3 [koło] | 9:MB | 0,95 ± 0,04 | 1,09 ± 0,04 | 0,97 ± 0,10 * |
| pT×B: (ΔR, pT) najbliższego dżetu pT>3 [koło] | 8:J12only | 0,92 ± 0,05 | 0,77 ± 0,05 | 0,70 ± 0,05 |
| pT×B: (ΔR, pT) najbliższego dżetu pT>3 [koło] | 8:MBTSbased | 0,79 ± 0,10 | 0,81 ± 0,11 * | 0,86 ± 0,11 * |
| pT×pT dżetu z D0 (kol. A) [koło] | 9:MB | 0,95 ± 0,03 | 1,06 ± 0,03 | 0,89 ± 0,08 * |
| pT×pT dżetu z D0 (kol. A) [koło] | 8:J12only | 0,90 ± 0,06 | 0,79 ± 0,05 | 0,85 ± 0,06 |
| pT×pT dżetu z D0 (kol. A) [koło] | 8:MBTSbased | 0,75 ± 0,08 | 0,84 ± 0,10 * | 1,10 ± 0,13 * |

## 6. (d) ε_iso w danych w funkcji ΔR do dżetu — uniwersalność

**ΔR(D0, najbliższy dżet kol. O pT > 5)** (`o_drn5`)

| próbka | 0–0.6 | 0.6–0.8 | 0.8–1 | 1–1.4 | 1.4–2 | 2–3.2 | 3.2–50 | brak dżetu |
|---|---|---|---|---|---|---|---|---|
| 9:MB | 0,06 ± 0,16 (26) | 0,36 ± 0,18 (28) | 0,68 ± 0,19 (25) | 0,57 ± 0,10 (64) | 0,40 ± 0,11 (68) | 0,54 ± 0,05 (204) | 0,57 ± 0,07 (91) | 0,83 ± 0,03 (450) |
| 8:J12only | 0,21 ± 0,05 (190) | 0,26 ± 0,06 (181) | 0,34 ± 0,05 (204) | 0,35 ± 0,04 (384) | 0,35 ± 0,04 (376) | 0,28 ± 0,02 (601) | 0,18 ± 0,05 (138) | 0,25 ± 0,05 (123) |
| 8:MBTSbased | 0,23 ± 0,15 (27) | 0,58 ± 0,19 (27) | 0,37 ± 0,13 (32) | 0,08 ± 0,21 (27) | 0,37 ± 0,09 (72) | 0,46 ± 0,07 (109) | – | 0,68 ± 0,09 (62) |
| 8:other | – | – | – | – | 0,54 ± 0,13 (35) | 0,33 ± 0,10 (56) | – | 0,70 ± 0,07 (69) |
| MC fd_nd | 0,43 ± 0,02 (946) | 0,51 ± 0,01 (1218) | 0,56 ± 0,01 (1192) | 0,59 ± 0,01 (2523) | 0,61 ± 0,01 (4041) | 0,59 ± 0,01 (8599) | 0,58 ± 0,01 (2823) | 0,75 ± 0,00 (19894) |

χ²/ndf (biny z N_sig > 25 w obu): 9:MB vs 8:J12only: 156,6/8 (p = 8.2e-30); 9:MB vs 8:MBTSbased: 11,0/7 (p = 0.14); 8:J12only vs 8:MBTSbased: 28,6/7 (p = 0.00017); 8:J12only vs 8:other: 28,4/3 (p = 3e-06); 9:MB vs MC: 17,4/8 (p = 0.026); 8:J12only vs MC: 486,2/8 (p = 6.4e-100); 8:MBTSbased vs MC: 22,3/7 (p = 0.0022); 8:other vs MC: 8,2/3 (p = 0.042).

**ΔR(D0, najbliższy dżet kol. O pT > 3)** (`o_drn3`)

| próbka | 0–0.6 | 0.6–0.8 | 0.8–1 | 1–1.4 | 1.4–2 | 2–3.2 | 3.2–50 | brak dżetu |
|---|---|---|---|---|---|---|---|---|
| 9:MB | 0,30 ± 0,12 (58) | 0,41 ± 0,10 (78) | 0,55 ± 0,10 (75) | 0,64 ± 0,07 (127) | 0,63 ± 0,06 (151) | 0,69 ± 0,05 (195) | 0,78 ± 0,07 (71) | 0,84 ± 0,03 (200) |
| 8:J12only | 0,25 ± 0,04 (388) | 0,30 ± 0,04 (326) | 0,39 ± 0,04 (295) | 0,29 ± 0,03 (455) | 0,34 ± 0,03 (351) | 0,21 ± 0,03 (301) | 0,18 ± 0,07 (46) | 0,31 ± 0,09 (33) |
| 8:MBTSbased | 0,11 ± 0,12 (45) | 0,41 ± 0,10 (63) | 0,29 ± 0,11 (52) | 0,37 ± 0,12 (52) | 0,51 ± 0,09 (62) | 0,61 ± 0,08 (72) | – | – |
| 8:other | – | 0,30 ± 0,13 (36) | – | 0,50 ± 0,11 (44) | 0,47 ± 0,12 (39) | 0,48 ± 0,15 (26) | – | 0,80 ± 0,09 (31) |
| MC fd_nd | 0,47 ± 0,01 (2792) | 0,55 ± 0,01 (3326) | 0,60 ± 0,01 (2931) | 0,64 ± 0,01 (5472) | 0,66 ± 0,01 (6640) | 0,66 ± 0,00 (9447) | 0,69 ± 0,01 (2346) | 0,81 ± 0,00 (8282) |

χ²/ndf (biny z N_sig > 25 w obu): 9:MB vs 8:J12only: 180,7/8 (p = 7.3e-35); 9:MB vs 8:MBTSbased: 10,5/6 (p = 0.11); 8:J12only vs 8:MBTSbased: 28,9/6 (p = 6.2e-05); 8:J12only vs 8:other: 21,6/5 (p = 0.00061); 9:MB vs MC: 7,1/8 (p = 0.52); 8:J12only vs MC: 594,7/8 (p = 3.2e-123); 8:MBTSbased vs MC: 27,2/6 (p = 0.00013); 8:other vs MC: 10,0/5 (p = 0.076).

**ΔR(D0, najbliższy dżet kol. B pT > 5) [koło dla ΔR < 0,4]** (`b_drn5`)

| próbka | 0–0.2 | 0.2–0.4 | 0.4–0.6 | 0.6–0.8 | 0.8–1 | 1–1.4 | 1.4–2 | 2–3.2 | 3.2–50 | brak dżetu |
|---|---|---|---|---|---|---|---|---|---|---|
| 9:MB | 0,00 ± 0,00 (43) | 0,00 ± 0,00 (32) | 0,16 ± 0,12 (41) | – | – | 0,63 ± 0,10 (59) | 0,51 ± 0,13 (50) | 0,64 ± 0,05 (173) | 0,59 ± 0,07 (84) | 0,85 ± 0,03 (436) |
| 8:J12only | 0,00 ± 0,00 (512) | -0,00 ± 0,00 (317) | 0,35 ± 0,08 (122) | 0,43 ± 0,08 (114) | 0,47 ± 0,06 (149) | 0,52 ± 0,04 (260) | 0,51 ± 0,04 (257) | 0,49 ± 0,03 (334) | 0,37 ± 0,08 (68) | 0,42 ± 0,07 (73) |
| 8:MBTSbased | – | -0,01 ± 0,00 (29) | 0,41 ± 0,17 (26) | – | 0,39 ± 0,13 (31) | – | 0,49 ± 0,10 (54) | 0,53 ± 0,07 (94) | – | 0,65 ± 0,07 (65) |
| 8:other | – | – | – | – | – | – | 0,48 ± 0,11 (41) | 0,44 ± 0,11 (44) | – | 0,82 ± 0,07 (59) |
| MC fd_nd | 0,00 ± 0,00 (1564) | 0,00 ± 0,00 (1494) | 0,40 ± 0,01 (1150) | 0,62 ± 0,02 (1002) | 0,64 ± 0,01 (1049) | 0,66 ± 0,01 (2254) | 0,68 ± 0,01 (3622) | 0,67 ± 0,01 (7563) | 0,66 ± 0,01 (2472) | 0,79 ± 0,00 (19066) |

χ²/ndf (biny z N_sig > 25 w obu): 9:MB vs 8:J12only: 45,3/6 (p = 4.1e-08); 9:MB vs 8:MBTSbased: 9,4/4 (p = 0.052); 8:J12only vs 8:MBTSbased: 6,9/6 (p = 0.33); 8:J12only vs 8:other: 16,4/3 (p = 0.00093); 9:MB vs MC: 13,1/6 (p = 0.042); 8:J12only vs MC: 107,0/9 (p = 6e-19); 8:MBTSbased vs MC: 21,5/6 (p = 0.0015); 8:other vs MC: 8,3/3 (p = 0.039).

**ΔR(D0, dżet wiodący kol. A)** (`a_dr_j1`)

| próbka | 0–0.1 | 0.1–0.4 | 0.4–0.8 | 0.8–1.4 | 1.4–2.4 | 2.4–2.9 | ≥ 2.9 |
|---|---|---|---|---|---|---|---|
| 9:MB | 0,70 ± 0,03 (536) | -0,06 ± 0,13 (35) | – | 0,69 ± 0,13 (43) | 0,66 ± 0,07 (107) | 0,77 ± 0,10 (65) | 0,63 ± 0,06 (146) |
| 8:J12only | 0,10 ± 0,01 (717) | 0,02 ± 0,01 (174) | 0,58 ± 0,20 (45) | 0,41 ± 0,08 (109) | 0,57 ± 0,04 (358) | 0,41 ± 0,04 (311) | 0,34 ± 0,03 (488) |
| 8:MBTSbased | 0,37 ± 0,04 (154) | – | – | – | 0,53 ± 0,14 (45) | 0,49 ± 0,12 (43) | 0,40 ± 0,11 (61) |
| 8:other | 0,46 ± 0,05 (118) | – | – | – | 0,54 ± 0,13 (34) | – | 0,68 ± 0,18 (28) |
| MC fd_nd | 0,65 ± 0,00 (23129) | 0,13 ± 0,01 (2504) | 0,67 ± 0,02 (807) | 0,78 ± 0,01 (1684) | 0,81 ± 0,01 (4467) | 0,75 ± 0,01 (3253) | 0,75 ± 0,01 (5392) |

χ²/ndf (biny z N_sig > 25 w obu): 9:MB vs 8:J12only: 440,6/6 (p = 5.2e-92); 9:MB vs 8:MBTSbased: 49,3/4 (p = 5e-10); 8:J12only vs 8:MBTSbased: 33,1/4 (p = 1.1e-06); 8:J12only vs 8:other: 45,8/3 (p = 6.2e-10); 9:MB vs MC: 13,9/6 (p = 0.031); 8:J12only vs MC: 2184,0/7 (p = 0); 8:MBTSbased vs MC: 60,0/4 (p = 3e-12); 8:other vs MC: 17,0/3 (p = 0.00072).

**pT dżetu kol. O w ΔR < 0,8 [GeV]** (`o_ptmax08`)

| próbka | 0–1e-06 | 1e-06–3 | 3–5 | 5–8 | 8–12 | 12–20 | ≥ 20 |
|---|---|---|---|---|---|---|---|
| 9:MB | 0,81 ± 0,03 (238) | 0,66 ± 0,03 (584) | 0,46 ± 0,10 (82) | 0,10 ± 0,21 (28) | – | – | – |
| 8:J12only | 0,38 ± 0,04 (265) | 0,28 ± 0,02 (1219) | 0,32 ± 0,04 (343) | 0,19 ± 0,06 (172) | 0,23 ± 0,08 (97) | 0,22 ± 0,13 (47) | 0,35 ± 0,11 (62) |
| 8:MBTSbased | 0,63 ± 0,08 (55) | 0,41 ± 0,05 (208) | 0,16 ± 0,11 (55) | 0,38 ± 0,14 (34) | – | – | – |
| 8:other | 0,65 ± 0,09 (45) | 0,50 ± 0,06 (136) | 0,41 ± 0,14 (35) | – | – | – | – |
| MC fd_nd | 0,77 ± 0,00 (10278) | 0,66 ± 0,00 (24840) | 0,54 ± 0,01 (3954) | 0,50 ± 0,01 (1507) | 0,44 ± 0,02 (424) | 0,37 ± 0,04 (168) | 0,40 ± 0,06 (65) |

χ²/ndf (biny z N_sig > 25 w obu): 9:MB vs 8:J12only: 190,2/4 (p = 4.8e-40); 9:MB vs 8:MBTSbased: 26,4/4 (p = 2.6e-05); 8:J12only vs 8:MBTSbased: 16,4/4 (p = 0.0025); 8:J12only vs 8:other: 19,6/3 (p = 0.0002); 9:MB vs MC: 6,2/4 (p = 0.19); 8:J12only vs MC: 582,5/7 (p = 1.4e-121); 8:MBTSbased vs MC: 36,8/4 (p = 2e-07); 8:other vs MC: 8,8/3 (p = 0.032).

**Σ pT półkuli przeciwnej [GeV]** (`h_s_away`)

| próbka | 0–4 | 4–7 | 7–10 | 10–15 | 15–22 | 22–35 | ≥ 35 |
|---|---|---|---|---|---|---|---|
| 9:MB | 0,94 ± 0,04 (49) | 0,85 ± 0,05 (76) | 0,86 ± 0,04 (104) | 0,81 ± 0,06 (121) | 0,65 ± 0,05 (213) | 0,58 ± 0,06 (204) | 0,38 ± 0,06 (193) |
| 8:J12only | – | 0,25 ± 0,10 (27) | 0,26 ± 0,08 (37) | 0,23 ± 0,06 (76) | 0,30 ± 0,05 (159) | 0,33 ± 0,03 (447) | 0,28 ± 0,02 (1441) |
| 8:MBTSbased | – | – | – | – | 0,49 ± 0,09 (56) | 0,48 ± 0,08 (100) | 0,25 ± 0,06 (169) |
| 8:other | – | – | – | 0,47 ± 0,10 (32) | 0,44 ± 0,10 (43) | 0,39 ± 0,11 (56) | 0,40 ± 0,09 (77) |
| MC fd_nd | 0,86 ± 0,01 (2414) | 0,81 ± 0,01 (3308) | 0,79 ± 0,01 (3737) | 0,73 ± 0,01 (6628) | 0,68 ± 0,01 (8349) | 0,59 ± 0,00 (10613) | 0,46 ± 0,01 (6187) |

χ²/ndf (biny z N_sig > 25 w obu): 9:MB vs 8:J12only: 161,6/6 (p = 2.7e-32); 9:MB vs 8:MBTSbased: 6,1/3 (p = 0.11); 8:J12only vs 8:MBTSbased: 6,8/3 (p = 0.08); 8:J12only vs 8:other: 7,1/4 (p = 0.13); 9:MB vs MC: 10,5/7 (p = 0.16); 8:J12only vs MC: 344,3/6 (p = 2.6e-71); 8:MBTSbased vs MC: 18,8/3 (p = 0.0003); 8:other vs MC: 15,8/4 (p = 0.0033).

W nawiasach N_sig (SBS; MC: liczba kandydatów). „brak dżetu” = w zdarzeniu nie ma dżetu kolekcji powyżej progu.

## 7. Zastrzeżenia

- Dżety torowe (tylko naładowane, pT > 0,5 GeV) ≠ dżety kalorymetryczne trygera L1_J12/HLT_j20; pT dżetu torowego ≈ 0,6 pT pełnego, z dużym rozrzutem.
- Kolekcja O obcina dżet, w którym leży D⁰ (tory w ΔR < 0,4 usunięte): „dżet O w ΔR < 0,8” to reszta dżetu wokół D⁰ albo dżet sąsiedni.
- MC nie ma trygera; w fd_nd twarde dżety są rzadkie → przy dużych pT dżetu komórki MC puste (kolumna „cof.”). Domknięcie przy dużym „cof.” nie jest dowodem.
- SBS zakłada kształt tła Δm niezależny od zmiennej wewnątrz grupy (bin pT × pass/fail); przy podziale w/poza dżetem kontrola = fit (kolumna „metoda”).
- Pary fałszywe w MC (D⁰ i D* z różnych wpisów o tym samym eventNumber; ułamek w tabeli wejść `same_entry`) są odrzucone.
