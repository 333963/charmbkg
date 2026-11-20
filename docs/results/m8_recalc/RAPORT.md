# M8_przeliczenie — wielkości z noty przeliczone nowym nominalnym modelem fitu D⁰ („M8”) (projekt B)

20.11.2026, ~00:15 CET. **Całość post hoc** (wszystkie stare liczby i uzyski M8 przy WP70 były znane). M8 wzięty z FITBIAS bez zmian
(`FITBIAS_2026-11-19/skrypty/fb_fit.py`, `fb_var.py` importowane, nie kopiowane): ograniczenie Gaussa (dmu 0,8 MeV; skala 10 %), szablon DMM + CM z mb_nd (T, WP80 danego binu)
z κ(ε) × 1,364 na danych, exp(Czebyszew 2), zakres 1,70–2,05; σ = max(HESSE, bootstrap po zdarzeniach). Test odrzucenia szablonu przez dane: p(χ²) < 10⁻³ (próg adekwatności z WPSCAN §2).
Bez AFP, pojedyncze D, bez jobów. Jedyne pobranie: pliki score per kandydat MC z EOS (`dnet4/app/d150`, fd_nd 24 + mb_nd 78 plików, 241 MB → `in/app_d150/`).
Kontrola mechaniki: M0 odtwarza WPSCAN T5 (0,117; 0,68 / 0,52 …), FITBIAS R1 (4105 ± 463 …), etap 4 (R_fit 1,1155 / 0,9469 / 0,8797; HESSE co do cyfry) i D23 (448,2 / 860,4 / 501,3 / 391,1).
Nie dostrajano niczego; tam, gdzie przepis M8 nie istnieje (bin 2–3 GeV, sieć, kategorie sieć × drzewo), użyto jawnie oznaczonych rozszerzeń/wariantów.

## Odpowiedź w skrócie

| # | wielkość | stare (M0 / etap 4) | nowe (M8) | wniosek |
|---|---|---|---|---|
| 1 | σ_N/N względem WP70, D⁰ 3–5 GeV, ε = 0,50 / 0,35 (HESSE) | 0,68 / 0,52 (427929); 0,75 / 0,60 (428770) | 0,76 ± 0,06 / 0,61 ± 0,06; 0,70 ± 0,07 / 0,58 ± 0,06 (Δ ≤ 1,4σ) | **utrzymany**: zysk 24–42 % (z max(HESSE, boot): 21–52 %) |
| 1 | to samo, 5–8 GeV | 0,83 / 0,85 (427929); 0,77 / 0,70 (428770) | 1,10 ± 0,07 / 0,99 ± 0,08 (+3,8σ / +1,8σ); 0,79 ± 0,08 / 0,88 ± 0,05 (+0,2σ / +2,1σ) | **zmieniony**: w 427929 zysku nie ma; w 428770 12–21 % (max: 6–12 %) |
| 1 | to samo, 8–30 GeV | 1,05 / 1,34; 0,91 / 0,93 | 1,07 / 1,15; 428770: M8 odrzucony w całym skanie (p ≤ 10⁻⁶) → M0 | utrzymany: WP70 w optimum |
| 2 | D⁰ 2–3 GeV: N (drzewo 427 / sieć 427 / drzewo 428 / sieć 428) | 448 ± 54 / 860 ± 95 / 501 ± 66 / 391 ± 55 | M8 bez szablonu: 452 ± 38 / 867 ± 84 / 459 ± 51 / 438 ± 54; z szablonem 3–5 GeV (wariant): 472 ± 41 / 886 ± 88 / 482 ± 54 / 449 ± 50 | N bez zmian (|Δ| ≤ 0,8σ), błąd × 0,78–0,87 (wobec S4 z tym samym bootstrapem) |
| 2 | Z_LR | 8,6 / 9,1 / 7,8 / 6,8 | 12,1 / 11,8 / 9,4 / 10,0 (wariant z szablonem: 12,1 / 11,6 / 9,4 / 9,9) | rośnie (szerszy zakres), nadal zależne od rzędu tła |
| 2 | y(D⁰)/y(D⁺) wobec oczekiwanego 0,82–1,03 | 1,31 ± 0,16 / 1,23 ± 0,15 / 1,53 ± 0,26 / 1,28 ± 0,22 | 1,32 ± 0,14 / 1,29 ± 0,14 / 1,58 ± 0,24 / 1,33 ± 0,21 (z szablonem: 1,38 / 1,32 / 1,63 / 1,39) | **nadwyżka nie maleje** (× 1,4–1,7 wobec środka oczekiwania; było × 1,3–1,65; z szablonem lekko rośnie) |
| 3 | R_fit − R_prawda (D⁰, WP70, ziarno 1), 3–5 / 5–8 / 8–30 | +0,158 ± 0,154 / +0,014 ± 0,055 / −0,025 ± 0,068 → bramka 0/3 | +0,188 ± 0,120 / +0,050 ± 0,050 / −0,017 ± 0,047 → bramka 1/3 (8–30) | **werdykt bez zmian** (A niespełnione); 8–30 GeV przechodzi na granicy (ziarna: TAK / TAK / NIE przy σ 0,047 / 0,047 / 0,052) |
| 3 | G_fit na MC z fitu łącznego (ziarno 1) | 1,50 / 1,25 / 1,05 | 1,54 / 1,29 / 1,04 | bez zmian |
| 4 | XSEC `WYNIKI.md` | D⁰ ze starego modelu (388 / 148 µb; D⁺/D⁰ 0,36 / 0,45 / 0,38) | uzgodnione z `tabele.md`: 376 ± 31 ± 57 / 132 ± 11 ± 24 µb; 0,397 / 0,466 / 0,431 | poprawione; w `tabele.md` został nieaktualny nawias i kolumny „po poprawce domknięcia” |

## 1. Skan punktu pracy D⁰ (WPSCAN / WPTIGHT) — `tabele/wpscan_m8.md`

- **Adekwatność M8 wzdłuż skanu (W0).** 427929: bez odrzuceń w 3 binach (najniższe p = 0,038 przy 8–30 WP70); 2 punkty na granicy dmu/sc w 5–8 (ε = 0,85; 0,55). 428770 3–5 i 5–8: M8 przyjęty dla ε ≤ 0,85
  (przy 3–5 GeV ε = 0,95 i 0,90 odrzucone: p = 6·10⁻¹⁸ i 10⁻⁴ — poza zakresem pytania). **428770 8–30: dane odrzucają szablon w 12 z 18 punktów (wszystkie ε = 0,30–0,85; p = 10⁻¹³…5·10⁻⁴; M0: 0,4–0,5) → zostaje M0** (jak w FITBIAS przy WP70).
- **3–5 GeV — wniosek utrzymany.** Stosunek σ_HESSE/N do WP70: 427929 0,68 ± 0,12 / 0,52 ± 0,09 → 0,76 ± 0,06 / 0,61 ± 0,06; 428770 0,75 ± 0,09 / 0,60 ± 0,07 → 0,70 ± 0,07 / 0,58 ± 0,06 (różnice sparowane +0,8σ / +1,4σ i −0,6σ / −0,3σ).
  Z max(HESSE, bootstrap): 0,61 / 0,48 i 0,79 / 0,64 (stare 0,45 / 0,33 — zawyżony zysk, bo WP70 w M0 miał bootstrap 1,8 × HESSE — i 0,72 / 0,67). Zysk 24–42 % (HESSE), 21–52 % (max).
  Błąd bezwzględny w ciasnym punkcie też maleje: 427929 σ/N = 0,069 / 0,054 (M0: 0,096 / 0,071), 428770 0,078 / 0,064 (0,097 / 0,090).
- **5–8 GeV — wniosek zmieniony.** 427929: 0,83 ± 0,08 / 0,85 ± 0,10 → **1,10 ± 0,07 / 0,99 ± 0,08** (+3,8σ / +1,8σ sparowane) — M8 poprawił WP70 (0,113 → 0,082), a ciasnych punktów nie (0,090 / 0,081); zysku nie ma.
  428770: 0,77 ± 0,09 / 0,70 ± 0,07 → 0,79 ± 0,08 / 0,88 ± 0,05 (max: 0,94 / 0,88) — zysk 12–21 % (HESSE) albo 6–12 % (max). Zamiast „15–30 %” jest „0 % w 427929, 6–21 % w 428770”.
- **8–30 GeV — bez zmiany:** 427929 1,07 ± 0,08 / 1,15 ± 0,35; 428770 (M0) 0,91 / 0,93 w HESSE, 1,43 / 0,88 z bootstrapem. WP70 w optimum.
- **Dryf N/ε_MC (W2):** 427929: d(0,50) −0,051 ± 0,062 → +0,005 ± 0,047, d(0,35) −0,125 ± 0,077 → −0,068 ± 0,050; 428770 (nominał: M8, w 8–30 M0): −0,111 ± 0,043 → −0,063 ± 0,040 i −0,090 ± 0,040 → −0,126 ± 0,032.
  W 427929 z M8 żadne d per komórka nie przekracza 1σ; w 428770 odchylenie przenosi się między punktami (5–8 GeV: d(0,50) −0,166 ± 0,069 → −0,050 ± 0,066, d(0,35) −0,029 ± 0,085 → −0,159 ± 0,049). Bramka domknięcia SF_rel na MC nadal niezaliczona (FITBIAS: σ 0,127 / 0,057 / 0,053) — bez wniosku o SF.
- **Pass/fail wewnątrz WP70 (W3; algebra WPTIGHT):** stosunek błędów [E łączny]/[N_70] HESSE dla D⁰: M0 0,80–1,10 → M8 0,92–1,09 w 5 komórkach z M8 (428770 8–30, M0: 0,90–0,96). Wniosek „SF_rel z tych samych danych zwraca błąd WP70” utrzymany i ostrzejszy
  (pozorne 0,80–0,82 w 427929 3–5 znika: 0,92–1,05).
- Domknięcie mb_nd wzdłuż skanu (W4, z FITBIAS): M8 przy 5–8 i 8–30 GeV w granicach ±0,08 dla wszystkich ε; 3–5 GeV WP70 −0,23 ± 0,10 zostaje, przy ε ≤ 0,50: −0,02…−0,06 ± 0,06.

## 2. Pik D⁰ w 2–3 GeV — `tabele/d23_m8.md`

- **Szablon DMM + CM dla 2–3 GeV nie istnieje:** mb_nd nie ma tego binu (T: tylko 3–5 / 5–8 / 8–30; `mb_nd_d150` też). Nominał tutaj = **M8 bez składowej DMM + CM** (ograniczenie Gaussa + zakres 1,70–2,05; „M8n”).
  Warianty z jawnym zastrzeżeniem: szablon i κ z 3–5 GeV — (a) wg przepisu M8, κ_T(ε ≤ 0,5) = 1,154 × 1,364; (b) κ z prawdy mb_nd 3–5 GeV dla tego samego klasyfikatora i ε (drzewo 1,06–1,10; sieć 0,59–0,93) × 1,364. (a) i (b) dają to samo N (≤ 0,6 %).
- **Uzyski** (punkty z werdyktu etapu 4; S4 = model etapu 4: kształt stały, 1,75–2,00; σ(S4) tutaj 49 / 103 / 64 / 62 z innym ziarnem bootstrapu niż w etapie 4: 54 / 95 / 66 / 55):
  427929 drzewo WP10: 448 → 452 ± 38 (wariant z szablonem 472 ± 41; ten fit formalnie nieważny wg kryterium EDM, N stabilne); 427929 sieć WP20: 860 → 867 ± 84 (886 ± 88); 428770 drzewo WP10: 501 → 459 ± 51 (482 ± 54); 428770 sieć WP10: 391 → 438 ± 54 (449 ± 50).
  Zmiany ≤ 0,8σ starego błędu; błąd × 0,78–0,87. p(χ²) = 0,06–0,78 (427929 drzewo 0,06 we wszystkich modelach, także starym).
- **Istotność** Z_LR: 8,6 / 9,1 / 7,8 / 6,8 → 12,1 / 11,8 / 9,4 / 10,0. Wzrost pochodzi z szerszego zakresu (lepiej określone tło), nie z kształtu sygnału. Zależność od rzędu tła zostaje: D23 dla zakresu 1,70–2,05 z exp(pol4) podaje 7,9 / 8,1 / 7,0 / 6,4
  (kształt stały) — tego M8 nie zmienia, bo rząd tła w M8 jest ustalony (2).
- **N/ε_MC:** stabilne między WP30 / WP20 / WP10 w M8 (np. 427929 drzewo 4044 ± 470 / 4770 ± 395 / 4693 ± 394); ⟨y⟩ z trzech rozłącznych wycinków ε < 0,3: 4551 ± 350 / 4428 ± 353 (427929 drzewo / sieć), 4822 ± 473 / 4051 ± 434 (428770).
- **Nadwyżka względem D⁺ nie maleje:** y(D⁰)/y(D⁺) = 1,31 ± 0,16 → 1,32 ± 0,14 (427929 drzewo), 1,23 ± 0,15 → 1,29 ± 0,14 (sieć), 1,53 ± 0,26 → 1,58 ± 0,24 i 1,28 ± 0,22 → 1,33 ± 0,21 (428770); wobec środka oczekiwania (0,82–1,03 → 0,925) to × 1,33–1,65 → × 1,39–1,71.
  Wariant z szablonem 3–5 GeV **zwiększa** N o 2–5 % (1,38 / 1,32 / 1,63 / 1,39): składowa DMM + CM leży głównie poniżej piku i w modelu bez niej jest wchłaniana przez tło, nie przez sygnał. Hipoteza „nadwyżka = niemodelowana składowa pikująca” nie znajduje potwierdzenia.

## 3. Sieć vs drzewo dla D⁰ (etap 4): R i G_fit — `tabele/R_m8.md`, `tabele/gfit_m8.md`

- **Histogramów łącznych (sieć × drzewo) lokalnie nie było** (`wps_hist.npz` ma klasyfikatory osobno; `D23_.../in/*_d150.npz` tylko 3–5 GeV). Kategorie b / n / t zbudowano z plików score per kandydat MC pobranych z EOS (bez joba). **Danych nie pobierano i nie otwierano** (R na danych dla D⁰ jest zablokowane bramką etapu 4).
- M8 dla kategorii łącznych jest **rozszerzeniem**: κ_c per kategoria z prawdy mb_nd (truth CM {2, 3} + DMM {6, 7} do sygnału; b / n / t = 1,14 / 1,10 / 1,89, 1,01 / 0,92 / 1,70, 0,89 / 0,93 / 0,95); kształt szablonu, ograniczenie i zakres jak w M8.
  Kontrola przybliżenia FITBIAS: κ z prawdy dla właściwego drzewa H przy WP70 = 1,23 / 1,11 / 0,90 wobec 1,22 / 1,09 / 0,92 z interpolacji T użytej w M8 — zgodne.
- **Domknięcie (WP70, ziarno 1):** R_fit − R_prawda = +0,188 ± 0,120 / +0,050 ± 0,050 / −0,017 ± 0,047 (M0: +0,158 ± 0,138 [etap 4: ± 0,154] / +0,014 ± 0,056 / −0,025 ± 0,063 [± 0,068]). σ maleje o 13–25 %; liczby FITBIAS dla pass/fail zagnieżdżonego (0,127 / 0,057 / 0,053) potwierdzone dla sieć × drzewo (0,120 / 0,050 / 0,047).
  Bramka (σ ≤ 0,05 ∧ |Δ| ≤ 0,04 ∧ ≤ 2,5σ): 3–5 NIE (σ); **5–8 NIE** — σ na progu (0,050 / 0,049 / 0,051 dla ziaren 1–3), ale odchylenie +0,050 / +0,063 / +0,088 przekracza 0,04 (w M0: +0,014 / +0,037 / +0,046; ta sama realizacja mb_nd, ziarna skorelowane);
  **8–30 TAK dla ziaren 1 i 2 (σ 0,047), NIE dla ziarna 3 (0,052)** — na granicy; rms ze 100 replik ma błąd ±7 %, więc „zaliczone” nie jest odporne. Bez składowej DMM + CM (M8n) to samo: zysk pochodzi z zakresu i ograniczenia kształtu.
- **Werdykt etapu 4 bez zmian:** A (domknięcie MC, ziarno 1) niespełnione w 2 z 3 binów → sieć D⁰ nadal bez kalibracji na próbce inkluzywnej. Formalnie M8 otwierałby pomiar R na danych tylko dla 8–30 GeV i tylko w 427929 (w 428770 dane odrzucają M8 przy 8–30 GeV, a z M0 bramka nie przechodzi) — nie wykonano (post hoc, wynik graniczny; wymaga pobrania 405 MB danych).
- **G_fit:** na MC z fitu łącznego 1,50 / 1,25 / 1,05 → 1,54 / 1,29 / 1,04. Na danych z fitów jednej kategorii (`wps_hist.npz`; sieć z κ z prawdy przy jej WP70 — rozszerzenie; σ = max(HESSE, 30 replik), błąd G rzędu 15–20 %):
  427929 1,39 / 1,23 / 0,67 → 1,27 / 0,99 / 1,04; 428770 1,03 / 0,76 / 1,05 → 1,11 / 1,35 / (M8 odrzucony: sieć p = 4·10⁻¹⁸, drzewo 10⁻⁷ → zostaje M0: 1,05). Próg 1,15 spełniony w 1 z 3 binów na run; liczby M0 były zdominowane przez niestabilność bootstrapu (np. 0,67).
  Domknięcie uzysku sieci przy WP70 na mb_nd: −0,36 ± 0,06 / −0,08 ± 0,06 / −0,11 ± 0,09 → −0,15 ± 0,08 / +0,03 ± 0,05 / +0,06 ± 0,07 (3–5 GeV przestaje być „nieadekwatny” wg progu 2,5σ WPSCAN).

## 4. XSEC: spójność `WYNIKI.md` z `out/tabele.md`

- `WYNIKI.md` (22:20) był w całości sprzed podmiany M8 (`tabele.md` 22:25). Poprawione liczby D⁰ i pochodne w §1, §2, §3, §4, §5 (uzysk na zdarzenie), §6 (poz. 2), §7; kopia `WYNIKI.md.bak_przed_m8`. D⁺ bez zmian.
- **Niespójność w `out/tabele.md` (plik generowany — nie ruszany, do poprawy w `xsec.py`):** hak `fitbias_override.json` podmienia N, stat, `fit_model_rel`, `fit_bias_rel`, ale nie `closure_b`. Dlatego (i) nawias w T1 przy obciążeniu fitu D⁰ (−27 ± 12 / −16 ± 7 / −14 ± 9) jest ze starego modelu
  (M8: −23 ± 10 / −4 ± 7 / +6 ± 8; sama systematyka 25,3 / 8,3 / 9,4 % jest poprawna), (ii) kolumny „po poprawce domknięcia MC” w T1 i T2 dzielą uzysk M8 przez domknięcie M0. Poprawnie [przeliczone ręcznie, w `WYNIKI.md` oznaczone W]:
  σ × L 828 / 182 / 58 (zamiast 864 / 209 / 71), D⁺/D⁰ 0,290 ± 0,063 / 0,403 ± 0,056 / 0,458 ± 0,095 (zamiast 0,278 / 0,351 / 0,375). Nominał (bez poprawki) nie jest dotknięty.
- Skutki M8 w odczycie: D⁰ 5–8 GeV 376 ± 31 ± 57 µb jest teraz 1,7σ poniżej CMS (było 1,1σ) i 1,6σ poniżej Pythii; podwójny stosunek runów D⁺/D⁰ przy 5–8 GeV odchyla się o 2,1σ (było 1,3σ), D⁰ wobec D* przy 8–30 GeV o 2,1σ (tam 427929 ma M8, a 428770 M0).

## Co się utrzymuje / zmienia

- Utrzymane: zysk z zaostrzenia przy 3–5 GeV (24–42 %); WP70 w optimum przy 8–30 GeV; brak zysku z SF_rel z tych samych danych; pik D⁰ w 2–3 GeV (N, pozycja, szerokość) i jego nadwyżka N/ε wobec D⁺; brak domknięcia R dla D⁰ i werdykt etapu 4; nominalny WP70.
- Zmienione: zysk przy 5–8 GeV (427929: brak; 428770: 6–21 %); Z_LR piku 2–3 GeV 9,4–12,1 (z zastrzeżeniem rzędu tła); domknięcie R przy 8–30 GeV na granicy bramki (2 z 3 ziaren); dryf D⁰ w 427929 zgodny z zerem przy ε = 0,50; liczby D⁰ w XSEC `WYNIKI.md`.

## Zastrzeżenia i braki

- Post hoc; M8 sam nie ma poprawnego pokrycia (FITBIAS: HESSE zaniżony ~1,3×); 30 replik bootstrapu w skanie (błąd rms ±13 %), 100 w R i w punktach 2–3 GeV.
- Rozszerzenia poza przepis M8: (i) 2–3 GeV — brak szablonu, nominał bez składowej DMM + CM; (ii) κ dla sieci i dla kategorii sieć × drzewo z prawdy mb_nd (w próbce, na której liczone jest domknięcie); kształt szablonu wspólny (T, WP80).
- Bin 2–3 GeV: brak MC tła → brak domknięcia M8 w tym binie; S/B w mb_nd 3–5 zawyżone (BR × 1,364).
- **Braki:** (1) R dla D⁰ na danych — nie otwarte (bramka; dane `app/d150` 1,5 GB na EOS, niepobrane); (2) systematyka tła (rząd ±1) w figurze zasługi skanu z M8 — nie przeliczona (wnioski noty opierają się na samej części statystycznej, T5);
  (3) Sec. 6.6 (testy składowych z MC) i Sec. 6.11 (przewidywane uzyski tagowane) — poza zakresem, nadal stary model; (4) `xsec.py` nie podmienia `closure_b` — `tabele.md` wymaga regeneracji po poprawce haka (nie wykonane: plik generowany, poza zakresem);
  (5) rząd tła exp(pol4) z M8 dla piku 2–3 GeV — nie liczony (M8 ustalony), cytowane liczby D23.

## Pliki

`RAPORT.md`; `NOTA_teksty_EN.md` (gotowe akapity i tabele do Sec. 5.11, 5.12, 6.12, 5.13); `skrypty/`: `m8_wpscan.py`, `m8_wpscan_rep.py` (pkt 1), `m8_d23.py` (pkt 2), `m8_R.py`, `m8_gfit.py` (pkt 3), `m8_tabs.py`;
`tabele/`: `wpscan_m8.md`, `d23_m8.md`, `R_m8.md`, `gfit_m8.md`; `out/`: json + logi; `in/app_d150/` (MC z EOS, 241 MB).
Zmienione poza katalogiem: `XSEC_2026-11-19/WYNIKI.md` (+ `.bak_przed_m8`); notatka i wiersz indeksu w `charm-bkg-notes`. Nie ruszano: noty, podsumowania, `ORCH`, `Q_lxplus`, `TAGGED_*`, `INCAL_*`, `tabele.md`, `xsec.py`.
