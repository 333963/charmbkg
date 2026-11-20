# DLZ_dplus — D⁺: H + Dalitz + cechy kinematyczne bez dE/dx (raport, 16.11.2026)

**Krótko.** Brakujący bieg FEAT (D⁺ `H+ALL_noiso_nopid` = dalitz + ptratio + maxdR + absimp, 8 cech) policzony.
Zysk Zproxy WP70 vs H: **MC +24 / +17 / +20 %** (3–5 / 5–8 / 8–30 GeV), **transfer ×1,41 / 1,24 / 1,22**, test data +40 / +29 / +22 %;
C2ST tła +0,006 (≤ 0,01), brak rzeźbienia (R bliżej 1 niż H). Względem H+dalitz: +12 / +5 / +10 % (MC), +18 / +6 / +4 % (transfer).
Wynik wcześniejszego „ALL_noiso” (+153 % w teście data) był w całości artefaktem dE/dx. Zestaw zapisany w prerejestracji
`PROCEDURA.md` jako **H+DLZ**, z H+dalitz jako jedynym wariantem zapasowym. [zdanie o oczekiwanym zysku w parach D–D (projekt A) pominięte w kopii w repo charmbkg, 20.11.2026]
Żadnych liczb z sygnału D⁺ w danych ani z par nie liczono.

> **Uwaga do kopii w repo `charmbkg` (20.11.2026).** To kopia raportu z `charm_work/DLZ_dplus/`. Pominięto dwa fragmenty dotyczące
> projektu A (pary D–D): jedno zdanie w akapicie „Krótko” i treść §7; miejsca są oznaczone nawiasem kwadratowym. Pliku `PROCEDURA.md`
> (prerejestracja; zawiera wynik par) nie ma w repo. Reszta raportu bez zmian.

## 0. Ustawienia i przebieg
- Skrypty FEAT bez edycji, uruchamiane przez `scripts/dlz_run.py` (dodaje warianty do `GROUPS['Dplus']`, przekierowuje `BASE` na `DLZ_dplus/`;
  `frames` = dowiązanie do `FEAT_charmnet/frames`; score H, H_seed2, H+dalitz, H+ptratio, H+maxdR, H+absimp, H+ALL_noiso skopiowane z FEAT jako cache
  → te same repliki bootstrapu, liczby H identyczne z FEAT). Łańcuch: `scripts/run_dlz.sh`. Nowe warianty: `H+ALL_noiso_nopid`, `H+dalitz_ptratio`.
- Pierwsza próba (17:27) przerwana w trakcie foldu 0 (logi: `logs/attempt1_przerwany/`); ponowiona 18:07 (setsid), komplet, bez błędów.
- Protokół jak FEAT: XGB H, 5 foldów po `evkey` (blok zdarzenia), cap 300k, 200 replik bootstrapu Poissona po zdarzeniach; szum treningu 1,5 %, próg zysku 4,5 %.

## 1. Zproxy na MC (test mc), ΔZ/Z_H ± bootstrap sparowany — `out/Dplus_mc/zproxy.csv`
| wariant | WP70 3–5 | WP70 5–8 | WP70 8–30 | WP50 3–5 | WP80 3–5 |
|---|---|---|---|---|---|
| H (Z) | 6,68 | 5,57 | 5,16 | — | — |
| H_seed2 | −0,4 ± 0,6 | −0,7 ± 1,2 | +2,5 ± 3,2 | | |
| +dalitz | +10,3 ± 0,9 | +11,4 ± 1,7 | +9,2 ± 4,7 | +12,8 | +11,0 |
| +dalitz_ptratio | +17,2 ± 1,0 | +11,9 ± 1,9 | +17,4 ± 5,2 | +18,9 | +16,4 |
| **+ALL_noiso_nopid (H+DLZ)** | **+23,6 ± 1,1** | **+16,9 ± 2,1** | **+20,1 ± 5,4** | +28,6 | +23,2 |
| +ALL_noiso (z pid, FEAT) | +25,5 ± 1,2 | +16,4 ± 2,2 | +23,0 ± 5,8 | +27,0 | +24,9 |

- nopid ≈ ALL_noiso (różnica ≤ 2 % w 3–5 / 5–8): dE/dx CharmNet nic nie dodaje — zgodnie z FEAT (+pid sam +0,7 %).
- nopid vs dalitz: 1,121 / 1,049 / 1,100 (> 4,5 % we wszystkich binach; 5–8 granicznie).
- Ważność cech (gain, fold 0): dMinPtRatio 0,25, dMinPt 0,10, dMaxPtRatio 0,07, mKπ₁,₂/m_D 0,02 + 0,02, maxdR 0,013, absImpact 0,013.
  Największy wkład mają cechy pT córek (miękki ślad) → ryzyko dla ε_S w danych (wydajność śladów niskiego pT); stąd SF/sPlot w binach dMinPt w PROCEDURA §3.

## 2. Przeniesienie na dane (`transfer.py`, model MC) — `out/transfer_Dplus.csv`, WP70
| wariant | ε_B SB danych 3–5 / 5–8 / 8–30 | Z_transfer/Z_H |
|---|---|---|
| H | 0,0069 / 0,0073 / 0,0153 | 1 |
| H_seed2 | | 0,99 / 1,00 / 1,03 |
| +dalitz | 0,0048 / 0,0053 / 0,0111 | 1,20 / 1,17 / 1,18 |
| +dalitz_ptratio | 0,0041 / 0,0052 / 0,0103 | 1,29 / 1,18 / 1,21 |
| **+ALL_noiso_nopid** | 0,0034 / 0,0047 / 0,0102 | **1,41 / 1,24 / 1,22** (± 0,01–0,02) |
| +ALL_noiso (z pid) | 0,0035 / 0,0048 / 0,0102 | 1,40 / 1,23 / 1,22 |
ε_S^dane (korekta izolacji) = 0,707 / 0,701 / 0,704 ≈ ε_S MC — korekta izolacji nie działa na ten model (brak cech izolacji).

## 3. Test na danych (tło trening = dalekie SB danych) — `out/Dplus_data/zproxy.csv`, WP70
| wariant | 3–5 | 5–8 | 8–30 |
|---|---|---|---|
| H_seed2 | 0,0 ± 0,6 | −0,2 ± 0,9 | +0,2 ± 1,2 |
| +dalitz | +20,3 ± 1,1 | +15,9 ± 1,4 | +16,1 ± 2,0 |
| **+ALL_noiso_nopid** | **+40,3 ± 1,6** | **+28,5 ± 1,7** | **+22,3 ± 2,3** |
| +ALL_noiso (z pid) | +153 ± 4 | +150 ± 8 | +159 ± 11 |
Bez dE/dx test data zgadza się z transferem (1,40 / 1,29 / 1,22 vs 1,41 / 1,24 / 1,22) → nowe cechy nie uczą „MC vs dane”.
+150 % ALL_noiso to w całości dE/dx (C2ST 0,915).

## 4. Modelowanie tła (C2ST SB danych vs MC mb, ten sam podzbiór 300k + 300k co FEAT, 3 ziarna) — `out/modelling_Dplus.json`
| zestaw | C2ST (ziarna 0/1/2) | średnia | Δ vs H |
|---|---|---|---|
| H | 0,5894 / 0,5891 / 0,5895 | 0,5893 | — |
| +dalitz | 0,5899 / 0,5892 / 0,5905 | 0,5899 | +0,0005 |
| +dalitz_ptratio | 0,5943 / 0,5934 / 0,5942 | 0,5940 | +0,0046 |
| +ALL_noiso_nopid | 0,5948 / 0,5944 / 0,5958 | 0,5950 | **+0,0056** |
Spełnia kryterium ≤ 0,01; przyrost pochodzi głównie od stosunków pT (+0,004, jak FEAT). Rozrzut ziaren ±0,0005.

## 5. Rzeźbienie (`sculpt.py`) — `out/Dplus_{mc,data}/sculpt.json`
| model | wariant | R WP70 MC truth≥8 (SR) | R WP70 dane (pasmo) | nachylenie ε_B [/100 MeV] |
|---|---|---|---|---|
| test mc | H | 0,890 ± 0,013 | 0,946 ± 0,015 | −0,23 |
| test mc | +dalitz | 0,953 ± 0,015 | 0,970 ± 0,018 | −0,26 |
| test mc | +dalitz_ptratio | 0,944 ± 0,016 | 0,973 ± 0,019 | −0,26 |
| test mc | **+ALL_noiso_nopid** | **0,969 ± 0,018** | **0,998 ± 0,021** | −0,28 |
| test data | H | 0,893 ± 0,013 | 0,927 ± 0,015 | −0,25 |
| test data | **+ALL_noiso_nopid** | 0,972 ± 0,017 | 1,003 ± 0,021 | −0,30 |
Brak rzeźbienia; R bliżej 1 niż w H. Nachylenie ε_B(m) nieco większe (−0,28 vs −0,23 /100 MeV; Dalitz zawiera m_D) — gładkie, liniowe, bez piku;
bramka G3 dopuszcza ±0,10 od H (tu Δ = −0,05) i wymaga, by model tła fitu opisał ε_B(m) w SB.

## 6. Werdykt (kryterium FEAT §0)
| zestaw | zysk MC (3–5 / 5–8) | transfer | ΔC2ST | rzeźbienie | werdykt |
|---|---|---|---|---|---|
| H+dalitz | +10 / +11 % | 1,20 / 1,17 / 1,18 | +0,001 | nie | GO (zapas) |
| H+dalitz_ptratio | +17 / +12 % | 1,29 / 1,18 / 1,21 | +0,005 | nie | GO (pośredni, nie wybrany) |
| **H+ALL_noiso_nopid** | **+24 / +17 %** | **1,41 / 1,24 / 1,22** | **+0,006** | **nie** | **GO — wybrany (H+DLZ)** |
Wybór między dalitz a pełnym zestawem oparty wyłącznie na MC/transfer/C2ST; koszt pełnego zestawu: C2ST +0,005 więcej i większa zależność od
modelowania miękkich śladów (walidacja §3 PROCEDURA).

## 7. Oczekiwany zysk w parach (szczegóły PROCEDURA §6)
[treść §7 (liczby dotyczące par D–D, projekt A) pominięta w kopii w repo charmbkg, 20.11.2026; oryginał: `charm_work/DLZ_dplus/RAPORT.md`]

## 8. Pliki
- `PROCEDURA.md` — prerejestracja H+DLZ (zestaw, trening, WP, SF, bramki G1–G7, kolejność, raportowanie).
- `scripts/dlz_run.py`, `scripts/run_dlz.sh`; logi `logs/` (pierwsza przerwana próba `logs/attempt1_przerwany/`).
- `out/Dplus_mc/` (zproxy.csv, sculpt.json, scores_*), `out/Dplus_data/`, `out/transfer_Dplus.csv`, `out/modelling_Dplus.json`.
- FEAT_charmnet: nic nie nadpisane (tylko odczyt ramek i skopiowane cache score).

## 9. Wynik na danych (pełne h5 v2.0, oba runy; etap `dlz_data`, 17.11.2026 17:21) — NIE PRZYJĘTE wg kryteriów zapisanych z góry
Kryteria: `DATA_2026-11-17/KRYTERIA.md` (zapisane przed obejrzeniem piku). Tabele: `DATA_2026-11-17/tables.md`, werdykt `verdict.json`. WP70, biny pT 3–5 / 5–8 / 8–30 GeV.

| kryterium | 427929 | 428770 | wynik |
|---|---|---|---|
| K1 zysk σ(N)/N | +20,7 ± 2,9 / +9,1 ± 3,0 / +4,2 ± 2,9 % (średnio 11 %) | +32,9 ± 3,5 / +10,9 ± 3,1 / +6,4 ± 1,4 % (średnio 17 %) | PASS |
| K2 rzeźbienie | R 0,93 → 1,03 (3–5), 0,99 → 1,11 ± 0,04 (5–8); nachylenie −0,24 → −0,30 /100 MeV; χ² SB < 1,5 | R 0,93 → 1,03, 0,95 → 1,04; nachylenie −0,25 → −0,31 | FAIL (R w 5–8 i 8–30 dalej od 1 niż H) |
| K3 ΔSF_rel | 0,73 ± 0,07 / 0,88 ± 0,07 / 0,86 ± 0,06 | 0,53 ± 0,06 / 0,74 ± 0,09 / 0,88 ± 0,03 | FAIL |
| K4 ΔC2ST sygnału | −0,032 ± 0,006 | −0,028 ± 0,003 | PASS |
| K5 pik | szerokość DLZ/H 1,10 ± 0,04 (5–8) | 1,13 ± 0,05 (5–8), 1,04 ± 0,02 (8–30) | FAIL (ostrzeżenie) |
| K6 foldy | 0 naruszeń; AUC trening − OOF 0,006 | | PASS |
| K7 ziarna | ΔN_S < 1σ | ΔN_S < 1σ | PASS |

Interpretacja K3: ta sama wielkość policzona na MC minbias (gdzie prawdziwa wartość to 1) daje dla H+DLZ 0,59 ± 0,07 (3–5), 0,81 ± 0,07 (5–8), 0,93 ± 0,09 (8–30), a dla H 1,20 ± 0,14 / 0,92 ± 0,08 / 0,98 ± 0,10. Miara więc nie domyka się dla H+DLZ w 3–5 i 5–8 GeV: fit odniesienia przy WP90 zawyża N (N_DLZ/N_H przy WP90 = 2,3 w 427929 i 3,9 w 428770 w binie 3–5; χ²/ndf 1,8). Deficyt ΔSF_rel w danych ma tę samą wielkość co niedomknięcie na MC, więc nie jest dowodem złego modelowania cech, ale kalibracji ε_S dla H+DLZ tym fitem wykazać się nie da. W 8–30 GeV (domknięcie 0,93 ± 0,09) 428770 daje 0,88 ± 0,03 — nierozstrzygnięte.

Werdykt: zysk statystyczny potwierdzony na danych (K1), sygnał modelowany nie gorzej niż w H (K4), ale K2 i K3 nie przechodzą → H zostaje nominałem dla D⁺. Wariant zapasowy H+dalitz ma ten sam problem (ΔSF_rel 0,54 ± 0,11 / 0,92 ± 0,07 / 0,86 ± 0,06 w 427929). Warunek powrotu: pomiar SF, który domyka się na MC przy luźnym WP (odniesienie bez cięcia na score albo fit jednoczesny w kategoriach score), i model tła opisujący R > 1 w 5–8 GeV.

### 9.1 Obserwacja po fakcie (nie kryterium): ΔSF_rel podzielone przez domknięcie na MC
ΔSF_rel(dane)/domknięcie(MC mb), WP70, 3–5 / 5–8 / 8–30 GeV: 427929 1,24 ± 0,19 / 1,09 ± 0,13 / 0,92 ± 0,11; 428770 0,90 ± 0,15 / 0,91 ± 0,14 / 0,95 ± 0,10. Wszystko zgodne z 1 w granicach 10–20 %, czyli deficyt SF w danych jest tej samej wielkości co artefakt fitu. Zastrzeżenie: domknięcie liczone na MC minbias o innym S/B niż dane, więc korekty nie wolno stosować jako kalibracji; to tylko wskazanie, że cechy nie są źle modelowane na poziomie > 20 %. D⁺ nie ma tagu, więc odniesienia bez cięcia na score nie da się dopasować (3–5 GeV: 14,5 mln kandydatów na kilka tysięcy sygnału) — ograniczenie jest zasadnicze, nie techniczne.

### 9.2 Skan cięcia i obciążenie fitu (R1/R2 z RECON_lit_B3; `SCAN_2026-11-17/`, 17.11 23:00) — warunek przyjęcia nadal niespełniony
Kryteria zapisane z góry (`SCAN_2026-11-17/KRYTERIA.md`). Liczone z histogramów po cięciu i 200 replik etapu `dlz_data` (danych per kandydat lokalnie nie ma): WP 50/70/80/90 (bez WP60), sygnał pseudo-danych parametryczny, tło exp(pol4) z pasm bocznych, 100 pseudo-eksperymentów na komórkę; w 3–5 GeV kształt ustalony na MC.
- R = (N/ε_MC)_DLZ/(N/ε_MC)_H przy WP70, surowe → po korekcie obciążenia: 427929 1,08 → 1,02 ± 0,04 / 1,10 → 1,06 ± 0,04 / 1,00 → 0,95 ± 0,04; 428770 1,25 → 1,11 ± 0,07 / 1,17 → 1,01 ± 0,05 / 1,03 → 0,98 ± 0,02. Wszystko w 2σ od 1; średnio po WP i binach 0,990 ± 0,015 (427929) i 1,010 ± 0,012 (428770). → ε_S(H+DLZ)/ε_S(H) = 1 w granicach 4–7 % per komórka; deficyt K3 to w większości artefakt fitu.
- Obciążenie fitu exp(pol2) przy WP70, H / H+DLZ: 427929 1,01/1,07, 0,99/1,03, 1,00/1,05; 428770 1,01/1,14, 0,93/1,07, 1,00/1,05. H+DLZ zawyża N_S o 7–14 % w 3–5 GeV i 3–7 % wyżej (tło po cięciu DLZ ma krzywiznę, której exp(pol2) nie opisuje; χ²/ndf w pasmach 2–3,6 przy WP90) → **część zysku K1 to większy N_S, nie mniejszy błąd**.
- Werdykt wg kryteriów: 427929 nierozstrzygnięte (test czuły tylko w 5–8 GeV), 428770 FAIL (płaskość 8–30 GeV, gdzie samo H nie jest płaskie, χ² 13,5/3). Kryterium czułości na bezwzględnym obciążeniu okazało się kruche (post hoc: odjęcie ogonów sygnału z pasm zmienia klasyfikację komórek, R_corr bez zmian).
- Zastrzeżenia: R2 nie widzi składowych piko-podobnych (D_s→KKπ, D⁺→πππ, resztkowy D*), więc obciążenie jest dolnym ograniczeniem; 428770 jetowy przy MC minbias.
- Warunek powrotu: tło exp(pol3/4) dla obu klasyfikatorów, powtórzone K1 i skan; w 3–5 GeV nie luźniej niż WP70. H zostaje nominałem.

### 9.3 K1b: zysk z elastyczniejszym tłem (z histogramów; `SCAN_2026-11-17/K1b/`, 17.11 23:40) — niepotwierdzony wg kryterium z góry
Zysk w definicji K1 odtworzony dokładnie dla pol2/pol2. F-test: H+DLZ w 3–5 GeV wymaga pol3 w obu runach (p = 4·10⁻⁵, 8·10⁻¹¹), H zostaje przy pol2; N_S(H+DLZ) spada o 5,5 % i 8,6 %, N_DLZ/N_H z 1,08/1,25 do 1,02/1,14.
Zysk g = [σ/N]_H/[σ/N]_DLZ − 1 przy WP70 (uwaga: inna definicja niż K1, g = g_K1/(1−g_K1)), rząd z F-testu: 427929 +14,6 ± 8,2 / +10,0 ± 7,1 / +4,4 ± 8,1 %, średnio +9,7 ± 4,5 %; 428770 +31,1 ± 10,3 / +28,9 ± 8,9 (sztuczne: H pol4, DLZ pol2; przy równych rzędach +12,3 lub +7,9) / +5,8 ± 5,1 %, średnio +21,9 ± 4,9 % (≈ +16 % przy równych rzędach). Pull width 0,87–1,18 (błąd fitu pokrywa rozrzut). Obciążenie > 5 % w 3 komórkach, ale słabo wyznaczone (tło prawdy exp(pol6), wariant pasm zmienia je o 3–10 %).
Zastrzeżenia: błędy zawyżone ok. 2× (bootstrap po binach, niesparowany). → rozstrzygnięcie w K1c na danych per kandydat (pobrane z EOS do `SCAN_2026-11-17/cand/`, 1,3 GB).

### 9.4 K1c: wersja na danych per kandydat (`SCAN_2026-11-17/K1c/`, 18.11 00:30) — NIEROZSTRZYGNIĘTE; H zostaje
Odtworzenie `res_RUN` dokładne (24 fity, ΔN = 0). 200 sparowanych replik po zdarzeniach, 200 pseudo-eksperymentów z prawdziwymi kandydatami MC. Zysk WP70 w definicji K1 [%], 3–5 / 5–8 / 8–30 GeV / średnio:
| run | pol2/pol2 | F-test | pol3/pol3 | pol4/pol4 |
|---|---|---|---|---|
| 427929 | 20,6±3,3 / 9,1±2,9 / 4,2±2,8 / 11,3±1,7 | 12,7±3,9 / 9,1±2,9 / 4,2±2,8 / **8,7±1,9** | 19,2±3,6 / 8,3±3,3 / 4,4±3,6 / **10,7±2,0** | 17,6±4,1 / 6,0±4,1 / 3,0±4,3 / 8,9±2,4 |
| 428770 | 32,6±3,4 / 10,9±3,1 / 6,4±1,4 / 16,7±1,6 | 23,5±4,1 / 22,3±3,2 (H pol4, DLZ pol2) / 5,4±2,6 / **17,1±1,9** | 29,7±4,0 / 9,3±3,5 / 5,8±1,7 / **14,9±1,9** | 26,1±4,6 / 7,2±4,3 / 5,4±2,6 / 12,9±2,3 |
- Zysk realny, prawie cały z 3–5 GeV; w 427929 na progu 10 % (zależnie od reguły rzędu tła), w 428770 13–17 %. Rząd wybierany osobno dla klasyfikatorów przesuwa zysk o 7–12 pkt → jeśli kiedykolwiek, to wspólny rząd.
- Błąd fitu wiarygodny dla pol2–pol4 (pull width 0,89–1,19; σ_bootstrap/σ_fit 0,89–1,08).
- R = (N/ε)_DLZ/(N/ε)_H przy pol3/pol3, bez korekty: 427929 WP70 1,03 / 1,08 / 0,96 (± 0,04–0,05), średnio po WP 0,99–1,02 ± 0,025 → deficyt K3 był artefaktem exp(pol2). Nadmiar przy WP80–90 w 3–5 GeV (1,18–1,21); w 428770 3–5 GeV R rośnie z WP także dla pol4 (1,04 → 1,63) — otwarte.
- Obciążenia fitu tą metodą nie da się wyznaczyć na 5 %: wynik zależy od wariantu pasm (b różni się o 0,03–0,26), a kontrola w przerwie bez sygnału myli się w tle o −2,5…+8,6 % (V1), ±15 % (V2). „Kryterium obciążenia niespełnione” = „nie wykazano”.
- Próg zysku nie był ślepy (wartości centralne znane z K1b).
- Warunek powrotu: narzędzie do obciążenia niezależne od interpolacji pasm (MC minbias z pełnym tłem przy S/B danych albo kontrola na D⁰ z tagiem D*), wspólny rząd pol3 jako reguła zapisana z góry.

### 9.5 K1d: obciążenie fitu z MC minimum-bias z pełnym tłem (`SCAN_2026-11-17/K1d/`, 18.11 16:10) — obciążenie „w normie” NIE wykazane; linia H+DLZ zamknięta, H zostaje
Kryteria zapisane z góry (`KRYTERIA_K1d.md`). b = N_fit/N_true, WP70, wspólne tło exp(pol3), H / H+DLZ, 3–5 / 5–8 / 8–30 GeV (bootstrap po zdarzeniach MC, 200 replik):
- MC wprost: 1,04 ± 0,05 / 1,09 ± 0,04; 1,11 ± 0,07 / 1,12 ± 0,06; 0,98 ± 0,07 / 1,01 ± 0,05. Δ(DLZ − H): +0,05 ± 0,04, +0,01 ± 0,05, +0,03 ± 0,04.
- MC przeważone do S/B i kształtu pasm 427929: praktycznie to samo (wagi 0,86–1,16). Do 428770: wagi 3,5–16, fit niestabilny (5–8 GeV: 1,51 ± 0,22 / 1,52 ± 0,18) — ekstrapolacja, nie pomiar.
- Kryterium |b − 1| < 5 %: niespełnione dla H+DLZ w 3–5 i 5–8 GeV (ok. 2σ od 1); H przechodzi tylko dzięki większemu błędowi przy tej samej wartości centralnej w 5–8 GeV. Różnica DLZ − H zgodna z 0 przy WP70. Żadna komórka nie ma mocy na próg 5 % (σ_b 4–8 %).
- „Fałszywy sygnał” to kombinatoryka, nie refleksje (3–5 GeV: +5,1 % H / +8,5 % H+DLZ z jednej klasy prawdy); częściowo zrekonstruowany powab netto ≤ 1 %. Post hoc: w 3–5 GeV Δ rośnie z WP (+0,069 ± 0,027 przy WP50, +0,159 ± 0,047 przy WP80, 3,4σ) → H+DLZ rzeźbi kombinatorykę pod pikiem przy luźnych WP.
- Zysk prawdziwy na MC (σ/N_true, pol3): +12,2 ± 2,1 % — zgodny z danymi 427929 (+10,7 ± 2,0 %).
- Wniosek: obciążenie obu klasyfikatorów przy pol3 to +4…+12 % w 3–8 GeV, podobne dla H i H+DLZ przy WP70; progu 5 % tą statystyką MC nie da się wykazać. Rozstrzygnięcie wymagałoby ≥ 4× większego MC minbias. Uwaga uboczna dla nominału: także H ma b = 1,11 ± 0,07 w 5–8 GeV przy pol3.
