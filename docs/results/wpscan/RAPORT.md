# WPSCAN — skan punktu pracy drzewa H i sieci z prawdziwym fitem piku (projekt B)

19.11.2026, ~21:05 CET. Kryteria: `KRYTERIA.md` (zapisane przed histogramami i fitami pełnej statystyki; §4 = doprecyzowanie pulli po teście dymnym mechaniki).
**Cały wynik jest post hoc względem etapów 4 i 5a** (znane były σ_N/N przy WP70 i proxy z pasm). Bez AFP, tylko pojedyncze D, okno piku inkluzywnego otwarte.

Wykonanie: jeden job condora (etap ORCH `wpscan`, klaster 12793588, 175 s, 1 rdzeń) → `runs/lxplus/wps_hist.npz` (masa 1 MeV × 19 przedziałów score, progi z MC fd_nd
dla ε = 0,10…0,95, 30 replik bootstrapu po zdarzeniach) z gotowych plików `app/d150` etapu 4 (komplet: fd_nd 24, mb_nd 78, 427929 458, 428770 97 plików). Fity lokalnie
(`wps_fit.py`, 648 zadań, ~15 tys. fitów, 10 min): model `dn4_ana.fit_cats` bez zmian. Odczyt: `wps_report.py` (T1–T4, wg kryteriów), `wps_extra.py` (T5–T6, post hoc).
Pełne skany per punkt: `wps_tables.md`; liczby: `wps_fit.json`, `wps_summary.json`.
Liczba porównań: 12 komórek × 2 klasyfikatory × 17 punktów vs WP70 = 408 stosunków (skorelowanych); decyzji P1 12 (drzewo) + 12 (sieć, opisowo), P3 12. Bez poprawki na wielokrotność.

## Wnioski

1. **Drzewo H — gdzie zaostrzenie realnie zmniejsza błąd uzysku.** Wynik odporny (T5: sama σ_HESSE/N przy ustalonych ε, bez wyboru optimum):
   - 3–5 GeV: błąd × 0,68 / 0,52 (D⁰ 427929, ε = 0,50 / 0,35), × 0,75 / 0,60 (D⁰ 428770), × 0,66 / 0,80 (D⁺ 427929), × 0,74 / 0,64 (D⁺ 428770);
   - 5–8 GeV: D⁰ × 0,83 / 0,85 (427929) i × 0,77 / 0,70 (428770); D⁺ × 1,01 / 1,11 (427929 — brak zysku) i × 0,80 / 0,88 (428770);
   - 8–30 GeV: × 0,91–1,34 — WP70 jest w optimum albo zaostrzenie szkodzi.
   Proxy z pasm (Z/Z_70 = 2,0 i 1,5 dla D⁰ w 427929) zawyżało: realnie odpowiada mu 1/0,52 ≈ 1,9 (3–5 GeV, tylko przy ε = 0,35) i 1/0,85 ≈ 1,2 (5–8 GeV).
2. **Odczyt wg kryteriów (T1, F = stat ⊕ tło ±1, wygładzone):** P1 „opłaca się” w 7/12 komórek drzewa: D⁰ 427929 5–8 (ε_opt 0,35; g = 0,65 ± 0,10) i 8–30 (0,60; 0,83 ± 0,07),
   D⁰ 428770 3–5 (0,30; 0,38 ± 0,11), 5–8 (0,45; 0,52 ± 0,08), 8–30 (0,40; 0,56 ± 0,11), D⁺ 427929 8–30 (0,80; 0,84 ± 0,07), D⁺ 428770 3–5 (0,20; 0,71 ± 0,06).
   Niepotwierdzone: D⁰ 427929 3–5 (g = 0,80 ± 0,11, warunek 2σ), D⁺ 427929 5–8 (0,80 ± 0,15); zjedzone przez dryf: D⁺ 427929 3–5, D⁺ 428770 8–30; płasko: D⁺ 428770 5–8.
   **Zastrzeżenie do P1 przy 8–30 GeV:** zysk formalny pochodzi tam z systematyki tła przy WP70 (pojedyncza duża różnica rzędu ±1: 0,154 / 0,088 / 0,077 wobec stat 0,10 / 0,05 / 0,085),
   a stosunek samych błędów statystycznych jest ≥ 0,85–1,16 — tych trzech komórek NIE traktuję jako wskazania do zmiany WP. Systematyka tła (max z dwóch różnic) skacze między
   sąsiednimi punktami o czynnik kilka; wygładzanie 3-punktowe tego nie usuwa.
3. **N/ε_MC dryfuje w dół przy zaostrzaniu (T2, T5, T6) — to jest główny koszt.** Drzewo: średnio po binach d(0,50) = −5 ± 6 % / −11 ± 4 % (D⁰ 427929 / 428770), −4,5 ± 2,8 % / −5,3 ± 2,5 % (D⁺);
   d(0,35) = −12,5 ± 7,7 % / −9,0 ± 4,0 % (D⁰), −8,5 ± 3,2 % / −12,2 ± 2,3 % (D⁺). Na prawdzie MC (mb_nd wobec ε z fd_nd) drzewo nie dryfuje (0 ± 1–2 %; największe −3,6 ± 2,3 %),
   więc to różnica dane/MC w kształcie score, nie skład próbki MC: **SF drzewa maleje z zaostrzaniem o ~5 % (ε = 0,5) do ~10 % (ε = 0,35)**, mocniej w 428770 przy 8–30 GeV (−13 ± 3 %).
   Wg kryterium P2 dryf jest istotny w optimum (> 2σ) w 3/12 komórek drzewa (D⁺ 427929 3–5; D⁺ 428770 5–8 i 8–30) i tam przewyższa zysk (g_eff > 1); w komórkach D⁰ jest tego samego rzędu
   (−7…−21 %), ale błąd dryfu (6–16 %) nie pozwala go ani potwierdzić, ani wykluczyć. Zysk statystyczny przy 3–5 GeV (błąd 12 % → 6–8 %) jest tej samej wielkości co dryf (5–20 %):
   **ciasny punkt pracy wymaga własnego pomiaru SF przy tym ε** (SF z tagu D* przy WP70 się nie przenosi); bez niego zmiana WP zamienia błąd statystyczny na systematykę ε.
4. **Sieć w optimum vs drzewo w optimum (T3):** D⁺ 427929: ρ = 0,66 ± 0,09 / 0,91 ± 0,08 / 1,00 ± 0,04 (3–5 / 5–8 / 8–30 GeV; przy WP70 było 0,72 / 0,84 / 0,90 — G_fit etapu 4 = 1,11 / 1,03 ⇔ 0,90 / 0,97;
   tu dodatkowo systematyka tła i wygładzanie); D⁺ 428770: 0,71 ± 0,04 / 0,54 ± 0,05 / 0,72 ± 0,05. D⁰ 427929: — / 1,00 ± 0,09 / 0,81 ± 0,12; D⁰ 428770: — / 1,17 ± 0,10 / 0,70 ± 0,06; D⁰ 3–5 GeV nierozstrzygnięte
   (sieć przy WP70 nie domyka się na MC: −0,36 ± 0,06). Czyli: po zoptymalizowaniu punktu pracy obu klasyfikatorów **sieć w czystym minimum-bias daje realny zysk tylko dla D⁺ przy 3–5 GeV
   (~34 %; sama stat. HESSE przy ustalonym ε = 0,5: ~18 %)**, przy 5–8 i 8–30 GeV 0–9 % (nieistotne), dla D⁰ nic przy 5–8 i 19 ± 12 % przy 8–30; w 428770 zysk 28–46 % dla D⁺ i 30 % dla D⁰ 8–30.
   Ale N/ε sieci dryfuje mocniej niż drzewa: D⁰ już na prawdzie MC (−6…−13 % przy ε = 0,35: ε z fd_nd nie przenosi się na mb_nd — efekt wolnego pionu D*), na danych −16 ± 6 % (427929) i −24 ± 4 % (428770);
   D⁺ w 427929 stabilne (−2 ± 3 %), w 428770 −20 ± 2 %. Zysk sieci dotyczy więc precyzji N, nie N/ε.
5. **Adekwatność modelu (T4):** χ² i pulle na danych w porządku we wszystkich 432 fitach (żaden punkt nie odpadł przez p(χ²) < 10⁻³ ani pull ≥ 4,5; najniższe p ≈ 0,003 dla sieci D⁺ 8–30; poza tym 5 punktów na granicy dmu/sc i 1 fit nieważny).
   Odbicie K↔π przy ciasnych cięciach nie psuje fitu D⁰ (χ²/ndf 26–63/44 w całym skanie). Punkty odpadają przez domknięcie na MC (pseudo-dane mb_nd, błąd 4–12 %, sąsiednie punkty skorelowane):
   D⁰ sieć 3–5 GeV przy ε ≥ 0,70 (−0,30…−0,36), D⁰ drzewo 5–8 przy 0,60–0,65 (−0,15…−0,19), D⁺ 3–5 (−0,10…−0,12 przy ε 0,25–0,50), D⁺ 5–8 przy ε ≥ 0,65 (+0,12…+0,29).
   Skala szerokości sc skacze między sąsiednimi punktami przy małym S/B (0,83–1,16), a N nie jest monotoniczne w ε (podzbiory zagnieżdżone!) — fit przy luźnym WP i niskim pT jest słabo
   określony; bootstrap daje tam błąd do 1,8 × HESSE (np. D⁰ 427929 3–5 WP70: 530 → 970). Modelu nie poprawiałem.

**Rekomendacja:** nominału nie zmieniać dziś. Kandydat do sprawdzenia: drzewo H przy ε ≈ 0,4–0,5 w binie 3–5 GeV (oba gatunki) i dla D⁰ przy 5–8 GeV — pod warunkiem pomiaru SF przy tym ε
(tag D* dla D⁰; dla D⁺ brak metody — tylko dryf z tego skanu). Przy 8–30 GeV zostaje WP70. Optymalizacja punktu pracy daje przy niskim pT więcej (25–50 % w błędzie stat.) niż sieć przy WP70
(3–11 % dla D⁺ w 427929), ale kosztuje systematykę ε tej samej wielkości.

## Zastrzeżenia
- Post hoc; optimum wybierane na tych samych danych (17 punktów) — g z T1 zaniżone efektem wyboru; dlatego wnioski opieram na T5 (ustalone ε).
- Błąd g i ρ obejmuje tylko część statystyczną (30 replik → sam błąd rms ±13 %); systematyka tła jest zaszumiona i dominuje F w wielu komórkach WP70.
- Jedno ziarno sieci; sieć D⁰ nie ma kalibracji ε na próbce inkluzywnej (etap 4).
- Dryf d mierzy zmianę SF względem WP70, nie SF bezwzględny; zakłada poprawny model fitu wzdłuż skanu (na MC N_fit/N_prawda odchyla się w części punktów o 10–20 %, z błędem 4–12 %).
- 428770 = „all” (mieszanka klas triggera; etap 4 pokazał SF zależny od klasy przy 8–30 GeV) — klasy J12only / notJ12only są w histogramach, nie fitowane.
- Bin 2–3 GeV poza zakresem (węzeł `ana23` etapu 4 wciąż w kolejce).

## Nierozstrzygnięte
- D⁰ 3–5 GeV, sieć (oba runy): WP70 nie domyka się na MC → brak porównania sieć/drzewo w tym binie.
- Czy dryf N/ε drzewa dla D⁰ jest rzeczywisty per bin (błędy 6–16 %): rozstrzygnie SF z tagu D* przy ε = 0,4–0,5.
- D⁰ 427929 3–5: zysk stat. duży (× 0,5–0,7), ale formalnie niepotwierdzony (g = 0,80 ± 0,11).

## Tabele wg kryteriów

### T1. Optimum punktu pracy i zysk względem WP70 (F = √(stat² + tło²)/N; wygładzone 3 punkty)

| gatunek | próbka | pT | klasyf. | ε_opt | F(WP70) [stat / tło] | F(ε_opt) [stat / tło] | g = F_opt/F_70 | g surowe | g sam stat | sąsiedzi F/F70 | werdykt P1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D0 | 427929 | 3-5 | drzewo | 0,20 | 0,246 [0,214 / 0,170] | 0,196 [0,072 / 0,153] | 0,798 ± 0,109 | 0,621 | 0,338 | 0,54 / 0,81 | zysk niepotwierdzony (2σ) |
| D0 | 427929 | 3-5 | sieć | — | — | — | — | — | — | — | NIEROZSTRZYGNIĘTE: WP70 nieadekwatny |
| D0 | 427929 | 5-8 | drzewo | 0,35 | 0,137 [0,113 / 0,038] | 0,089 [0,081 / 0,027] | 0,646 ± 0,097 | 0,719 | 0,722 | 0,80 / 0,72 | OPŁACA SIĘ |
| D0 | 427929 | 5-8 | sieć | 0,40 | 0,138 [0,092 / 0,136] | 0,089 [0,081 / 0,015] | 0,643 ± 0,120 | 0,503 | 0,887 | 0,58 / 0,54 | OPŁACA SIĘ |
| D0 | 427929 | 8-30 | drzewo | 0,60 | 0,164 [0,103 / 0,154] | 0,137 [0,109 / 0,067] | 0,833 ± 0,071 | 0,689 | 1,055 | 0,78 / 0,74 | OPŁACA SIĘ |
| D0 | 427929 | 8-30 | sieć | 0,50 | 0,132 [0,154 / 0,041] | 0,111 [0,105 / 0,012] | 0,841 ± 0,109 | 0,659 | 0,677 | 0,74 / 0,69 | zysk niepotwierdzony (2σ) |
| D0 | 428770 | 3-5 | drzewo | 0,30 | 0,258 [0,135 / 0,221] | 0,097 [0,096 / 0,010] | 0,375 ± 0,111 | 0,371 | 0,707 | 0,39 / 0,36 | OPŁACA SIĘ |
| D0 | 428770 | 3-5 | sieć | — | — | — | — | — | — | — | NIEROZSTRZYGNIĘTE: WP70 nieadekwatny |
| D0 | 428770 | 5-8 | drzewo | 0,45 | 0,133 [0,083 / 0,094] | 0,069 [0,058 / 0,028] | 0,520 ± 0,078 | 0,512 | 0,695 | 0,63 / 0,51 | OPŁACA SIĘ |
| D0 | 428770 | 5-8 | sieć | 0,30 | 0,113 [0,109 / 0,066] | 0,081 [0,080 / 0,008] | 0,722 ± 0,149 | 0,631 | 0,734 | 0,60 / 0,68 | zysk niepotwierdzony (dryf zjada zysk) |
| D0 | 428770 | 8-30 | drzewo | 0,40 | 0,094 [0,046 / 0,088] | 0,053 [0,053 / 0,016] | 0,563 ± 0,106 | 0,556 | 1,160 | 0,49 / 0,56 | OPŁACA SIĘ |
| D0 | 428770 | 8-30 | sieć | 0,60 | 0,042 [0,044 / 0,007] | 0,037 [0,038 / 0,005] | 0,887 ± 0,060 | 0,859 | 0,862 | 0,80 / 0,86 | zysk niepotwierdzony (2σ) |
| Dplus | 427929 | 3-5 | drzewo | 0,45 | 0,108 [0,077 / 0,031] | 0,069 [0,070 / 0,011] | 0,637 ± 0,065 | 0,853 | 0,908 | 0,82 / — | zysk niepotwierdzony (dryf zjada zysk) |
| Dplus | 427929 | 3-5 | sieć | 0,50 | 0,078 [0,071 / 0,031] | 0,046 [0,042 / 0,010] | 0,588 ± 0,117 | 0,562 | 0,598 | — / 0,62 | zysk niepotwierdzony (dryf zjada zysk) |
| Dplus | 427929 | 5-8 | drzewo | 0,50 | 0,060 [0,046 / 0,029] | 0,048 [0,046 / 0,013] | 0,795 ± 0,147 | 0,878 | 1,006 | 0,90 / 0,88 | zysk niepotwierdzony (2σ) |
| Dplus | 427929 | 5-8 | sieć | 0,45 | 0,050 [0,047 / 0,005] | 0,044 [0,043 / 0,013] | 0,865 ± 0,045 | 0,959 | 0,922 | 0,93 / 0,89 | OPŁACA SIĘ |
| Dplus | 427929 | 8-30 | drzewo | 0,80 | 0,093 [0,085 / 0,077] | 0,078 [0,072 / 0,011] | 0,844 ± 0,066 | 0,637 | 0,848 | 0,75 / 0,66 | OPŁACA SIĘ |
| Dplus | 427929 | 8-30 | sieć | 0,80 | 0,083 [0,081 / 0,011] | 0,078 [0,074 / 0,020] | 0,933 ± 0,021 | 0,940 | 0,916 | 1,01 / 0,92 | mały zysk, bez rekomendacji |
| Dplus | 428770 | 3-5 | drzewo | 0,20 | 0,084 [0,077 / 0,025] | 0,059 [0,065 / 0,010] | 0,707 ± 0,063 | 0,813 | 0,847 | 0,78 / 0,65 | OPŁACA SIĘ |
| Dplus | 428770 | 3-5 | sieć | 0,50 | 0,071 [0,060 / 0,027] | 0,042 [0,040 / 0,007] | 0,588 ± 0,076 | 0,608 | 0,656 | — / 0,65 | zysk niepotwierdzony (dryf zjada zysk) |
| Dplus | 428770 | 5-8 | drzewo | 0,30 | 0,074 [0,057 / 0,047] | 0,073 [0,056 / 0,034] | 0,978 ± 0,073 | 0,885 | 0,985 | 1,08 / 0,98 | WP70 w optimum (płasko) |
| Dplus | 428770 | 5-8 | sieć | 0,45 | 0,064 [0,054 / 0,047] | 0,039 [0,036 / 0,004] | 0,611 ± 0,053 | 0,509 | 0,671 | 0,56 / 0,57 | zysk niepotwierdzony (dryf zjada zysk) |
| Dplus | 428770 | 8-30 | drzewo | 0,40 | 0,048 [0,033 / 0,018] | 0,035 [0,031 / 0,019] | 0,733 ± 0,045 | 0,970 | 0,950 | 0,91 / 0,93 | zysk niepotwierdzony (dryf zjada zysk) |
| Dplus | 428770 | 8-30 | sieć | 0,45 | 0,031 [0,026 / 0,034] | 0,025 [0,024 / 0,009] | 0,830 ± 0,075 | 0,596 | 0,924 | 0,61 / 0,58 | zysk niepotwierdzony (dryf zjada zysk) |

### T2. Stabilność N/ε_MC wzdłuż skanu (d = [N/ε]/[N/ε]_70 − 1; punkty adekwatne, ε ∈ [0,2; 0,95])

| gatunek | próbka | pT | klasyf. | rms d | χ²/ndf (d = 0) | największe odchylenie (ε: d ± σ) | d(ε_opt) | istotny w optimum | g_eff (z dryfem) |
|---|---|---|---|---|---|---|---|---|---|
| D0 | 427929 | 3-5 | drzewo | 0,138 | 6,8/14 | 0,20: -0.208 ± 0.158 | -0.208 ± 0.158 | nie | 0,798 |
| D0 | 427929 | 5-8 | drzewo | 0,081 | 7,0/13 | 0,85: -0.237 ± 0.106 | -0.075 ± 0.125 | nie | 0,646 |
| D0 | 427929 | 5-8 | sieć | 0,177 | 23,9/14 | 0,75: +0.207 ± 0.093 | +0.145 ± 0.150 | nie | 0,643 |
| D0 | 427929 | 8-30 | drzewo | 0,116 | 20,3/15 | 0,40: -0.200 ± 0.092 | -0.008 ± 0.063 | nie | 0,833 |
| D0 | 427929 | 8-30 | sieć | 0,167 | 54,4/15 | 0,35: -0.255 ± 0.076 | -0.083 ± 0.071 | nie | 0,841 |
| D0 | 428770 | 3-5 | drzewo | 0,141 | 7,6/13 | 0,85: +0.355 ± 0.274 | -0.146 ± 0.143 | nie | 0,375 |
| D0 | 428770 | 5-8 | drzewo | 0,107 | 20,4/13 | 0,50: -0.166 ± 0.069 | -0.134 ± 0.070 | nie | 0,520 |
| D0 | 428770 | 5-8 | sieć | 0,178 | 150,3/15 | 0,25: -0.247 ± 0.042 | -0.244 ± 0.074 | TAK | 2,284 |
| D0 | 428770 | 8-30 | drzewo | 0,114 | 86,4/15 | 0,55: -0.195 ± 0.045 | -0.092 ± 0.050 | nie | 0,563 |
| D0 | 428770 | 8-30 | sieć | 0,205 | 514,2/15 | 0,20: -0.427 ± 0.035 | -0.060 ± 0.041 | nie | 0,887 |
| Dplus | 427929 | 3-5 | drzewo | 0,151 | 64,6/13 | 0,20: -0.217 ± 0.064 | -0.152 ± 0.071 | TAK | 1,543 |
| Dplus | 427929 | 3-5 | sieć | 0,140 | 47,6/11 | 0,20: -0.216 ± 0.057 | -0.154 ± 0.057 | TAK | 2,065 |
| Dplus | 427929 | 5-8 | drzewo | 0,059 | 17,1/11 | 0,95: -0.139 ± 0.064 | -0.056 ± 0.037 | nie | 0,795 |
| Dplus | 427929 | 5-8 | sieć | 0,016 | 2,6/12 | 0,75: -0.022 ± 0.021 | +0.023 ± 0.046 | nie | 0,865 |
| Dplus | 427929 | 8-30 | drzewo | 0,058 | 11,9/15 | 0,35: -0.093 ± 0.054 | -0.048 ± 0.061 | nie | 0,844 |
| Dplus | 427929 | 8-30 | sieć | 0,050 | 13,0/15 | 0,40: +0.104 ± 0.066 | +0.002 ± 0.039 | nie | 0,933 |
| Dplus | 428770 | 3-5 | drzewo | 0,104 | 9,9/13 | 0,75: -0.092 ± 0.054 | +0.044 ± 0.088 | nie | 0,707 |
| Dplus | 428770 | 3-5 | sieć | 0,138 | 65,3/10 | 0,20: -0.249 ± 0.050 | -0.140 ± 0.043 | TAK | 2,060 |
| Dplus | 428770 | 5-8 | drzewo | 0,116 | 43,5/11 | 0,30: -0.176 ± 0.057 | -0.176 ± 0.057 | TAK | 2,567 |
| Dplus | 428770 | 5-8 | sieć | 0,138 | 71,7/12 | 0,25: -0.201 ± 0.052 | -0.127 ± 0.045 | TAK | 2,079 |
| Dplus | 428770 | 8-30 | drzewo | 0,107 | 172,3/15 | 0,30: -0.198 ± 0.029 | -0.105 ± 0.031 | TAK | 2,289 |
| Dplus | 428770 | 8-30 | sieć | 0,180 | 732,7/15 | 0,20: -0.307 ± 0.023 | -0.137 ± 0.021 | TAK | 4,528 |

### T3. Sieć w swoim optimum vs drzewo w swoim optimum (ρ = F_min(sieć)/F_min(drzewo))

| gatunek | próbka | pT | ε_opt drzewo / sieć | F_min drzewo | F_min sieć | ρ (optimum) | ρ przy WP70 | ρ sam stat (optimum) | werdykt P3 |
|---|---|---|---|---|---|---|---|---|---|
| D0 | 427929 | 3-5 | — | — | — | — | — | — | NIEROZSTRZYGNIĘTE |
| D0 | 427929 | 5-8 | 0,35 / 0,40 | 0,089 | 0,089 | 1,002 ± 0,085 | 1,007 | 1,001 | sieć nie lepsza |
| D0 | 427929 | 8-30 | 0,60 / 0,50 | 0,137 | 0,111 | 0,813 ± 0,119 | 0,805 | 0,960 | zysk mały |
| D0 | 428770 | 3-5 | — | — | — | — | — | — | NIEROZSTRZYGNIĘTE |
| D0 | 428770 | 5-8 | 0,45 / 0,30 | 0,069 | 0,081 | 1,172 ± 0,099 | 0,845 | 1,384 | sieć nie lepsza |
| D0 | 428770 | 8-30 | 0,40 / 0,60 | 0,053 | 0,037 | 0,699 ± 0,056 | 0,444 | 0,710 | sieć realnie lepsza |
| Dplus | 427929 | 3-5 | 0,45 / 0,50 | 0,069 | 0,046 | 0,662 ± 0,088 | 0,718 | 0,611 | sieć realnie lepsza |
| Dplus | 427929 | 5-8 | 0,50 / 0,45 | 0,048 | 0,044 | 0,910 ± 0,078 | 0,837 | 0,941 | zysk mały |
| Dplus | 427929 | 8-30 | 0,80 / 0,80 | 0,078 | 0,078 | 0,996 ± 0,035 | 0,901 | 1,024 | zysk mały |
| Dplus | 428770 | 3-5 | 0,20 / 0,50 | 0,059 | 0,042 | 0,706 ± 0,035 | 0,848 | 0,610 | sieć realnie lepsza |
| Dplus | 428770 | 5-8 | 0,30 / 0,45 | 0,073 | 0,039 | 0,537 ± 0,045 | 0,859 | 0,647 | sieć realnie lepsza |
| Dplus | 428770 | 8-30 | 0,40 / 0,45 | 0,035 | 0,025 | 0,719 ± 0,048 | 0,635 | 0,764 | sieć realnie lepsza |

### T4. Adekwatność fitu — punkty odrzucone w ε ∈ [0,2; 0,95]

| gatunek | próbka | pT | klasyf. | adekwatnych / 16 | odrzucone (ε: powód) |
|---|---|---|---|---|---|
| D0 | 427929 | 3-5 | drzewo | 15 | 0.95: domknięcie MC +1.06 ± 0.39 |
| D0 | 427929 | 3-5 | sieć | 11 | 0.70: domknięcie MC -0.36 ± 0.06; 0.75: domknięcie MC -0.33 ± 0.08; 0.85: domknięcie MC -0.30 ± 0.10; 0.90: domknięcie MC -0.30 ± 0.10; 0.95: granica dmu/sc |
| D0 | 427929 | 5-8 | drzewo | 14 | 0.60: domknięcie MC -0.15 ± 0.06; 0.65: domknięcie MC -0.19 ± 0.06 |
| D0 | 427929 | 5-8 | sieć | 15 | 0.85: granica dmu/sc |
| D0 | 427929 | 8-30 | drzewo | 16 | — |
| D0 | 427929 | 8-30 | sieć | 16 | — |
| D0 | 428770 | 3-5 | drzewo | 14 | 0.90: granica dmu/sc; 0.95: granica dmu/sc, domknięcie MC +1.06 ± 0.39 |
| D0 | 428770 | 3-5 | sieć | 12 | 0.70: domknięcie MC -0.36 ± 0.06; 0.75: domknięcie MC -0.33 ± 0.08; 0.85: domknięcie MC -0.30 ± 0.10; 0.90: domknięcie MC -0.30 ± 0.10 |
| D0 | 428770 | 5-8 | drzewo | 14 | 0.60: domknięcie MC -0.15 ± 0.06; 0.65: domknięcie MC -0.19 ± 0.06 |
| D0 | 428770 | 5-8 | sieć | 16 | — |
| D0 | 428770 | 8-30 | drzewo | 16 | — |
| D0 | 428770 | 8-30 | sieć | 16 | — |
| Dplus | 427929 | 3-5 | drzewo | 14 | 0.30: domknięcie MC -0.12 ± 0.04; 0.50: domknięcie MC -0.12 ± 0.05 |
| Dplus | 427929 | 3-5 | sieć | 12 | 0.25: domknięcie MC -0.10 ± 0.04; 0.35: domknięcie MC -0.11 ± 0.03; 0.40: domknięcie MC -0.10 ± 0.03; 0.45: domknięcie MC -0.10 ± 0.03 |
| Dplus | 427929 | 5-8 | drzewo | 12 | 0.65: domknięcie MC +0.13 ± 0.05; 0.80: domknięcie MC +0.17 ± 0.05; 0.85: domknięcie MC +0.21 ± 0.06; 0.90: domknięcie MC +0.21 ± 0.07 |
| Dplus | 427929 | 5-8 | sieć | 13 | 0.80: domknięcie MC +0.12 ± 0.04; 0.85: domknięcie MC +0.15 ± 0.05; 0.95: domknięcie MC +0.29 ± 0.08 |
| Dplus | 427929 | 8-30 | drzewo | 16 | — |
| Dplus | 427929 | 8-30 | sieć | 16 | — |
| Dplus | 428770 | 3-5 | drzewo | 14 | 0.30: domknięcie MC -0.12 ± 0.04; 0.50: domknięcie MC -0.12 ± 0.05 |
| Dplus | 428770 | 3-5 | sieć | 11 | 0.25: domknięcie MC -0.10 ± 0.04; 0.35: domknięcie MC -0.11 ± 0.03; 0.40: domknięcie MC -0.10 ± 0.03; 0.45: domknięcie MC -0.10 ± 0.03; 0.80: fit nieważny |
| Dplus | 428770 | 5-8 | drzewo | 12 | 0.65: domknięcie MC +0.13 ± 0.05; 0.80: domknięcie MC +0.17 ± 0.05; 0.85: domknięcie MC +0.21 ± 0.06; 0.90: domknięcie MC +0.21 ± 0.07 |
| Dplus | 428770 | 5-8 | sieć | 13 | 0.80: domknięcie MC +0.12 ± 0.04; 0.85: domknięcie MC +0.15 ± 0.05; 0.95: domknięcie MC +0.29 ± 0.08 |
| Dplus | 428770 | 8-30 | drzewo | 16 | — |
| Dplus | 428770 | 8-30 | sieć | 16 | — |

## Diagnostyki post hoc (poza kryteriami)

### T5 (post hoc). Ustalone punkty ε = 0,50 i 0,35 względem WP70 — bez wyboru optimum

| gatunek | próbka | pT | klasyf. | σ_HESSE/N: WP70 | stosunek (0,50) | stosunek (0,35) | d dane (0,50) | d dane (0,35) | d prawda MC mb_nd (0,50) | (0,35) |
|---|---|---|---|---|---|---|---|---|---|---|
| D0 | 427929 | 3-5 | drzewo | 0,117 | 0,68 | 0,52 | -0,056 ± 0,179 | -0,062 ± 0,220 | +0,000 ± 0,011 | -0,003 ± 0,017 |
| D0 | 427929 | 3-5 | sieć | 0,114 | 0,60 | 0,54 | -0,133 ± 0,106 | -0,168 ± 0,097 | -0,035 ± 0,011 | -0,064 ± 0,017 |
| D0 | 427929 | 5-8 | drzewo | 0,096 | 0,83 | 0,85 | +0,024 ± 0,138 | -0,075 ± 0,125 | -0,000 ± 0,014 | +0,018 ± 0,023 |
| D0 | 427929 | 5-8 | sieć | 0,084 | 0,95 | 0,99 | +0,098 ± 0,119 | +0,293 ± 0,159 | -0,025 ± 0,015 | -0,085 ± 0,024 |
| D0 | 427929 | 8-30 | drzewo | 0,089 | 1,05 | 1,34 | -0,072 ± 0,075 | -0,179 ± 0,110 | -0,023 ± 0,022 | -0,016 ± 0,035 |
| D0 | 427929 | 8-30 | sieć | 0,154 | 0,60 | 0,80 | -0,083 ± 0,071 | -0,255 ± 0,076 | -0,062 ± 0,024 | -0,135 ± 0,036 |
| D0 | 428770 | 3-5 | drzewo | 0,124 | 0,75 | 0,60 | -0,044 ± 0,156 | -0,128 ± 0,143 | +0,000 ± 0,011 | -0,003 ± 0,017 |
| D0 | 428770 | 3-5 | sieć | 0,103 | 0,90 | 0,62 | -0,064 ± 0,138 | -0,123 ± 0,132 | -0,035 ± 0,011 | -0,064 ± 0,017 |
| D0 | 428770 | 5-8 | drzewo | 0,082 | 0,77 | 0,70 | -0,166 ± 0,069 | -0,029 ± 0,085 | -0,000 ± 0,014 | +0,018 ± 0,023 |
| D0 | 428770 | 5-8 | sieć | 0,109 | 0,48 | 0,54 | -0,130 ± 0,068 | -0,224 ± 0,062 | -0,025 ± 0,015 | -0,085 ± 0,024 |
| D0 | 428770 | 8-30 | drzewo | 0,043 | 0,91 | 0,93 | -0,080 ± 0,060 | -0,105 ± 0,047 | -0,023 ± 0,022 | -0,016 ± 0,035 |
| D0 | 428770 | 8-30 | sieć | 0,035 | 1,01 | 1,11 | -0,115 ± 0,039 | -0,267 ± 0,045 | -0,062 ± 0,024 | -0,135 ± 0,036 |
| Dplus | 427929 | 3-5 | drzewo | 0,073 | 0,66 | 0,80 | -0,064 ± 0,067 | -0,212 ± 0,068 | +0,003 ± 0,011 | +0,019 ± 0,018 |
| Dplus | 427929 | 3-5 | sieć | 0,048 | 0,82 | 0,81 | -0,154 ± 0,057 | -0,156 ± 0,065 | +0,023 ± 0,011 | +0,030 ± 0,018 |
| Dplus | 427929 | 5-8 | drzewo | 0,046 | 1,01 | 1,11 | -0,056 ± 0,037 | -0,011 ± 0,049 | -0,011 ± 0,015 | -0,036 ± 0,023 |
| Dplus | 427929 | 5-8 | sieć | 0,041 | 1,00 | 1,07 | +0,008 ± 0,044 | +0,027 ± 0,045 | -0,007 ± 0,015 | +0,004 ± 0,024 |
| Dplus | 427929 | 8-30 | drzewo | 0,058 | 1,05 | 1,20 | -0,010 ± 0,053 | -0,093 ± 0,054 | +0,010 ± 0,022 | -0,018 ± 0,035 |
| Dplus | 427929 | 8-30 | sieć | 0,055 | 1,01 | 1,17 | +0,074 ± 0,048 | +0,037 ± 0,071 | -0,034 ± 0,022 | -0,024 ± 0,034 |
| Dplus | 428770 | 3-5 | drzewo | 0,077 | 0,74 | 0,64 | +0,031 ± 0,084 | +0,027 ± 0,084 | +0,003 ± 0,011 | +0,019 ± 0,018 |
| Dplus | 428770 | 3-5 | sieć | 0,050 | 0,79 | 0,86 | -0,140 ± 0,043 | -0,219 ± 0,044 | +0,023 ± 0,011 | +0,030 ± 0,018 |
| Dplus | 428770 | 5-8 | drzewo | 0,053 | 0,80 | 0,88 | -0,057 ± 0,042 | -0,153 ± 0,052 | -0,011 ± 0,015 | -0,036 ± 0,023 |
| Dplus | 428770 | 5-8 | sieć | 0,037 | 0,95 | 1,06 | -0,096 ± 0,043 | -0,161 ± 0,056 | -0,007 ± 0,015 | +0,004 ± 0,024 |
| Dplus | 428770 | 8-30 | drzewo | 0,024 | 1,03 | 1,18 | -0,063 ± 0,033 | -0,129 ± 0,027 | +0,010 ± 0,022 | -0,018 ± 0,035 |
| Dplus | 428770 | 8-30 | sieć | 0,021 | 1,09 | 1,29 | -0,093 ± 0,018 | -0,198 ± 0,021 | -0,034 ± 0,022 | -0,024 ± 0,034 |

### T6 (post hoc). Średni ważony dryf d po 3 binach pT

| gatunek | klasyf. | próbka | d(0,50) | d(0,35) |
|---|---|---|---|---|
| D0 | drzewo | 427929 | -0,051 ± 0,062 | -0,125 ± 0,077 |
| D0 | drzewo | 428770 | -0,111 ± 0,043 | -0,090 ± 0,040 |
| D0 | sieć | 427929 | -0,060 ± 0,053 | -0,158 ± 0,056 |
| D0 | sieć | 428770 | -0,116 ± 0,033 | -0,243 ± 0,035 |
| Dplus | drzewo | 427929 | -0,045 ± 0,028 | -0,085 ± 0,032 |
| Dplus | drzewo | 428770 | -0,053 ± 0,025 | -0,122 ± 0,023 |
| Dplus | sieć | 427929 | -0,009 ± 0,028 | -0,017 ± 0,033 |
| Dplus | sieć | 428770 | -0,100 ± 0,016 | -0,198 ± 0,018 |

## Korekta (19.11, po WPTIGHT)
„Dryf N/ε_MC” z §2–§3 nie jest zależnością SF dane/MC od ε. Ten sam fit pass/fail na MC (prawda = 1,00) daje dla D⁺ 0,949 ± 0,023 / 0,935 ± 0,027 / 0,931 ± 0,028 (ε = 0,50 / 0,40 / 0,35), dane 0,949 ± 0,025 (427929) i 0,951 ± 0,016 (428770); stosunek podwójny dane/MC = 1,00 ± 0,04. Jest to obciążenie fitu (zawyżony sygnał w kategorii luźnej, małe S/B), w granicach ±4–5 %. Zysk statystyczny ciasnego punktu istnieje tylko przy SF znanym z niezależnej próbki: z SF_rel z tych samych danych wynik jest algebraicznie równy uzyskowi przy WP70. Źródło: `WPTIGHT_2026-11-19/RAPORT.md`.
