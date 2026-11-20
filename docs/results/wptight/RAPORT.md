# WPTIGHT — skalibrowany ciasny punkt pracy drzewa H: SF_rel z fitu pass/fail inkluzywnego piku i bilans błędu uzysku (projekt B)

19.11.2026, ~21:15 CET. Kryteria: `KRYTERIA.md` (zapisane przed wszystkimi fitami, także przed domknięciem MC). **Całość post hoc względem etapów 4 / 5a / WPSCAN.**
Praca lokalna, bez nowych jobów: wejście = `WPSCAN_2026-11-19/runs/lxplus/wps_hist.npz` (masa × 19 przedziałów score, 30 replik bootstrapu po zdarzeniach).
Tylko drzewo `xgb.H.sb.n.s1`, bez AFP, pojedyncze D. Fity: `wpt_fit.py` (model `dn4_ana.fit_cats` bez zmian; MC 18 zadań 2 min, dane 72 zadania 5 min, ~9 tys. fitów),
odczyt: `wpt_report.py` → `wpt_tables.md` (M1, D1–D6), liczby `wpt_summary.json`, surowe `wpt_mc.json`, `wpt_data.json`.
Liczba porównań: 36 werdyktów (12 komórek × 3 skorelowane ε), bez poprawki na wielokrotność; SF_rel per klasa triggera 36 liczb opisowo.
Kontrola mechaniki: fit jednej kategorii WP70 i dryf d odtwarzają WPSCAN co do cyfry (np. D⁺ 427929 3–5: stat 0,077 / tło 0,031; D⁰ 427929 3–5: d(0,50) = −0,056 ± 0,179).

## Odpowiedź

**Nie: ciasny punkt pracy z poprawką SF_rel zmierzoną na inkluzywnym piku nie daje odpornie mniejszego błędu całkowitego niż WP70.** Trzy powody, każdy wystarczający:

1. **Algebra (zapisana w kryteriach przed fitami).** N_ciasny/(f_MC·SF_rel) z SF_rel z tych samych danych to dokładnie N_A + N_B, czyli uzysk WP70 z fitu dwóch kategorii.
   Zysk statystyczny zaostrzenia z WPSCAN (25–48 %) istniał tylko przy założeniu, że ułamek f = ε(ciasny | WP70) jest znany; gdy f jest wolnym parametrem, błąd wraca:
   stosunek błędów HESSE [E łączny]/[N_70] = **0,81–0,90 tylko w D⁺ 427929 3–5 GeV i 0,87–0,91 w D⁺ 428770 5–8 GeV; w pozostałych 10 komórkach 0,90–1,08** (średnio ≈ 0,98).
   Korelacja N_A ↔ f jest dodatnia (0,2–0,6), więc suma naiwna σ(N_A) ⊕ σ(SF_rel) zawyża błąd o 10–25 % (D2) — ale nawet poprawny błąd nie schodzi poniżej WP70 o więcej niż ~10–19 %.
2. **Werdykty wg kryteriów (D1):** z 36 — D⁰ 18 × NIEROZSTRZYGNIĘTE (domknięcie), D⁺: 1 × „opłaca się”, 2 × „zysk niepotwierdzony”, 3 × „< 10 %”, 8 × „szkodzi”, 4 × nierozstrzygnięte (G2).
   Jedyny pozytywny (D⁺ 428770 5–8, ε = 0,40: g = 0,70 ± 0,06) ma sąsiadów ε = 0,50 i 0,35 z g = 1,40 i 1,30 — to szum systematyki tła (max z dwóch różnic: 0,086 / 0,022 / 0,082),
   nie wynik; warunek odporności (≥ 2 z 3 ε) niespełniony w żadnej komórce.
3. **Domknięcie na MC pokazuje, że SF_rel z fitu jest obciążony tak samo jak „dryf” na danych** (niżej) — poprawka SF_rel wnosiłaby własną systematykę ~5 %.

## Wyniki

### Domknięcie na MC (M1) — bramka G1 (|SF − 1| ≤ 2,5σ), G2 (σ ≤ 0,05)
| | 3–5 GeV (ε = 0,50 / 0,40 / 0,35) | 5–8 GeV | 8–30 GeV |
|---|---|---|---|
| D⁰ SF_rel^MC | 1,12 ± 0,20 / 1,31 ± 0,17 / 1,33 ± 0,18 | 1,05 ± 0,06 / 1,02 ± 0,08 / 1,04 ± 0,08 | 0,97 ± 0,07 / 1,08 ± 0,11 / 1,15 ± 0,13 |
| D⁰ bramka | NIE (G2) × 3 | NIE (G2) × 3 | NIE (G2) × 3 |
| D⁺ SF_rel^MC | 0,917 ± 0,041 / 0,960 ± 0,048 / 0,980 ± 0,053 | 0,966 ± 0,035 / 0,904 ± 0,043 / 0,906 ± 0,040 | 0,959 ± 0,044 / 0,948 ± 0,048 / 0,926 ± 0,060 |
| D⁺ bramka | TAK / TAK / NIE (G2) | TAK / TAK / TAK | TAK / TAK / NIE (G2) |
- **D⁰ nie przechodzi nigdzie** (σ = 0,17–0,20 / 0,06–0,08 / 0,07–0,13). Metoda R etapu 4 miała 0,15 / 0,06 / 0,07 — **zagnieżdżone cięcia tego samego klasyfikatora nie są lepsze** niż kategorie sieć/drzewo.
  Powód ten sam: kategoria luźna ma S/B w oknie 0,01–0,02 (3–5 GeV) i odbicie K↔π; dodatkowo uzysk D⁰ 3–5 GeV na MC nie domyka się w żadnym wariancie (E/N_prawda = 0,77 ± 0,11–0,14, N_70/N_prawda = 0,74 ± 0,16).
- **D⁺ przechodzi w 7/9 punktów, ale „bez mocy” i z jednym znakiem:** wszystkie 9 wartości < 1 (pulle −0,4…−2,4). Średnia po trzech binach (post hoc, biny niezależne):
  **SF_rel^MC = 0,949 ± 0,023 / 0,935 ± 0,027 / 0,931 ± 0,028** (ε = 0,50 / 0,40 / 0,35; −2,2 / −2,4 / −2,4σ; trzy ε skorelowane), przy prawdzie z liczenia 0,96–1,02 (odstępstwo tylko 5–8 GeV: 0,966 ± 0,020).
  Fit zawyża sygnał w kategorii luźnej (to samo widać w WPSCAN T4: D⁺ 5–8 przy ε ≥ 0,65 +0,12…+0,29, i w LOWPT §1: N_fit/N_true 0,99 → 1,17 ku luźnym WP). Bramka 2,5σ per komórka tego nie odrzuca.

### Dane: SF_rel ± stat ± tło i bilans błędu (D1; pełne tabele w `wpt_tables.md`)
| gatunek | run | pT | SF_rel (0,50) | SF_rel (0,40) | SF_rel (0,35) | σ_tot ciasny (0,50 / 0,40 / 0,35) | σ_tot WP70 | g ± σ_g (0,50 / 0,40 / 0,35) | werdykt |
|---|---|---|---|---|---|---|---|---|---|
| D⁺ | 427929 | 3–5 | 0,951 ± 0,050 ± 0,016 | 0,904 ± 0,054 ± 0,027 | 0,878 ± 0,053 ± 0,020 | 0,072 / 0,074 / 0,085 | 0,083 | 0,87 ± 0,10 / 0,90 ± 0,12 / 1,03 ± 0,15 | niepotwierdzony × 2, nierozstrz. (G2) |
| D⁺ | 428770 | 3–5 | 1,001 ± 0,065 ± 0,028 | 1,029 ± 0,076 ± 0,036 | 1,036 ± 0,085 ± 0,033 | 0,087 / 0,081 / 0,082 | 0,081 | 1,08 ± 0,10 / 1,00 ± 0,09 / 1,02 ± 0,07 | szkodzi × 2, nierozstrz. |
| D⁺ | 427929 | 5–8 | 0,941 ± 0,033 ± 0,019 | 0,925 ± 0,043 ± 0,021 | 0,969 ± 0,043 ± 0,018 | 0,050 / 0,061 / 0,055 | 0,054 | 0,93 ± 0,04 / 1,12 ± 0,13 / 1,02 ± 0,12 | < 10 %, szkodzi × 2 |
| D⁺ | 428770 | 5–8 | 0,969 ± 0,037 ± 0,004 | 0,943 ± 0,041 ± 0,014 | 0,916 ± 0,043 ± 0,019 | 0,103 / 0,052 / 0,096 | 0,074 | 1,40 ± 0,11 / 0,70 ± 0,06 / 1,30 ± 0,13 | szkodzi, OPŁACA SIĘ (nieodporne), szkodzi |
| D⁺ | 427929 | 8–30 | 0,956 ± 0,041 ± 0,006 | 0,935 ± 0,051 ± 0,005 | 0,948 ± 0,057 ± 0,006 | 0,107 / 0,106 / 0,081 | 0,115 | 0,93 ± 0,09 / 0,93 ± 0,08 / 0,71 ± 0,07 | < 10 % × 2, nierozstrz. |
| D⁺ | 428770 | 8–30 | 0,942 ± 0,018 ± 0,005 | 0,904 ± 0,022 ± 0,001 | 0,879 ± 0,025 ± 0,006 | 0,040 / 0,039 / 0,048 | 0,038 | 1,05 ± 0,11 / 1,03 ± 0,12 / 1,26 ± 0,13 | szkodzi × 2, nierozstrz. |
| D⁰ | 427929 | 3–5 | 0,94 ± 0,09 ± 0,05 | 1,04 ± 0,13 ± 0,11 | 1,02 ± 0,13 ± 0,09 | 0,217 / 0,172 / 0,182 | 0,273 | 0,79 ± 0,26 / 0,63 ± 0,23 / 0,67 ± 0,24 | NIEROZSTRZ. (domknięcie) — opisowo |
| D⁰ | 428770 | 3–5 | 0,92 ± 0,10 ± 0,05 | 0,90 ± 0,10 ± 0,08 | 0,89 ± 0,10 ± 0,11 | 0,244 / 0,199 / 0,174 | 0,259 | 0,94 ± 0,11 / 0,77 ± 0,10 / 0,67 ± 0,08 | j.w. |
| D⁰ | 427929 | 5–8 | 1,03 ± 0,09 ± 0,04 | 0,99 ± 0,10 ± 0,04 | 1,03 ± 0,11 ± 0,08 | 0,099 / 0,123 / 0,113 | 0,119 | 0,83 ± 0,15 / 1,03 ± 0,19 / 0,95 ± 0,29 | j.w. |
| D⁰ | 428770 | 5–8 | 0,93 ± 0,06 ± 0,03 | 0,97 ± 0,07 ± 0,03 | 0,94 ± 0,07 ± 0,01 | 0,105 / 0,090 / 0,186 | 0,126 | 0,84 ± 0,25 / 0,72 ± 0,21 / 1,48 ± 0,49 | j.w. |
| D⁰ | 427929 | 8–30 | 0,91 ± 0,08 ± 0,08 | 0,78 ± 0,08 ± 0,15 | 0,80 ± 0,08 ± 0,14 | 0,199 / 0,143 / 0,152 | 0,186 | 1,07 ± 0,11 / 0,77 ± 0,04 / 0,82 ± 0,05 | j.w. |
| D⁰ | 428770 | 8–30 | 0,946 ± 0,030 ± 0,029 | 0,924 ± 0,036 ± 0,021 | 0,922 ± 0,039 ± 0,015 | 0,075 / 0,060 / 0,060 | 0,100 | 0,75 ± 0,15 / 0,61 ± 0,10 / 0,60 ± 0,09 | j.w. |
- Kontrola Δ = E/N_70 − 1 zgodna z 0 wszędzie (|Δ| ≤ 1,8σ). Fity adekwatne we wszystkich 36 punktach (p(χ²) ≥ 0,016, bez granic dmu/sc).
- D⁰ (opisowo): formalne g = 0,6–0,8 pochodzi w 427929 3–5 z niestabilności fitu jednej kategorii przy WP70 (bootstrap 1,83 × HESSE; sam HESSE daje 0,93–0,96), a w 428770 8–30 i 3–5
  z systematyki tła WP70 (0,088 i 0,221, pojedyncza duża różnica rzędu ±1) — te same dwa artefakty co w WPSCAN, nie informacja z kategoryzacji.
- **„Dryf” to w granicach błędu obciążenie fitu, nie różnica dane/MC (post hoc, D3):** średnie po binach SF_rel dane D⁺ = 0,949 ± 0,025 / 0,923 ± 0,030 / 0,936 ± 0,031 (427929) i
  0,951 ± 0,016 / 0,918 ± 0,019 / 0,896 ± 0,022 (428770); stosunek podwójny dane/MC-fit = **1,00 ± 0,04 / 0,99 ± 0,04 / 1,01 ± 0,05** (427929) i 1,00 ± 0,03 / 0,98 ± 0,04 / 0,96 ± 0,04 (428770; MC wspólne dla obu runów).
  Per komórka D⁺ stosunek podwójny 0,90–1,09 z błędem 0,05–0,11. Wniosek WPSCAN §3 („SF drzewa maleje z zaostrzaniem o 5–10 %, bo prawda MC nie dryfuje”) **wymaga korekty**: prawda MC
  z liczenia nie dryfuje, ale fit na MC z pełnym tłem dryfuje o −5…−7 % (2,2–2,4σ), tyle samo co dane. Realna zmiana SF drzewa D⁺ między WP70 a ε = 0,35–0,50 jest zgodna z zerem na poziomie ±4–5 %.
  D⁰: stosunek podwójny 0,95 ± 0,07 … 0,84 ± 0,06 (428770, ε = 0,35: −2,7σ), ale MC-fit ma tam błędy 5–6 % po uśrednieniu i nie przechodzi bramki — bez wniosku.

### 428770: SF_rel per klasa triggera (D4, opisowo) — histogramy mają podział J12only / notJ12only
Brak istotnej zależności: |J12only − notJ12only| ≤ 1,4σ we wszystkich 18 punktach. D⁺ 8–30: 0,940 ± 0,020 vs 0,952 ± 0,048 (ε = 0,50), 0,874 ± 0,027 vs 0,895 ± 0,054 (0,35).
D⁰ 5–8: notJ12only niżej w trzech ε (0,85–0,86 ± 0,09 vs 0,97–1,04 ± 0,10; +1,0…+1,4σ, skorelowane, bramka niezaliczona). Zależności 3,3σ znanej dla sieci (etap 4, 8–30 GeV) drzewo nie pokazuje.
Drobniejszego podziału (MBTSbased / other) w histogramach nie ma — brak, bez nowego joba.

### Bilans bezwzględny N/ε (D5): σ_abs = σ_tot ⊕ σ[SF_70]
SF_rel kalibruje tylko ciasny względem WP70; SF_70 wchodzi do obu estymatorów jednakowo. Status SF_70: D⁰ ≥ 5 GeV — tag D* (WP95: 1,05 ± 0,11 i 0,91 ± 0,07 w 427929; 0,84 ± 0,10 i 1,01 ± 0,04 w 428770;
tabeli dla WP70 w bazie nie znalazłem, przyjmuję ten sam rząd 5–10 %); D⁰ 3–5 GeV i D⁺ — **nieznany** (LOWPT proponował ±10 % / ±15 % jako systematykę z testu stabilności).
- D⁺ 427929 3–5 (jedyna komórka z zyskiem ≥ 10 % przechodząca bramkę): 0,072 vs 0,083 → przy σ[SF_70] = 5 / 10 / 15 %: 0,088 vs 0,097 (g_abs 0,91), **0,123 vs 0,130 (0,95)**, 0,166 vs 0,171 (0,97).
- D⁺ pozostałe komórki: g_abs = 0,96–1,16 przy 10 %. D⁰ (opisowo): 3–5 GeV 0,20–0,24 vs 0,29 przy 10 % (g_abs 0,68–0,95), 5–8 GeV 0,84–1,31.
- **Gdy σ[SF_70] ≳ 10 %, ciasny punkt nie zmienia błędu całkowitego N/ε o więcej niż 5 % w żadnej komórce D⁺**; dominuje nieznany SF_70.
- Wariant „SF_rel z 428770 przeniesiony na 427929” (D6, dane niezależne): SF_rel obu runów zgodne (|Δ| ≤ 1,5σ), ale σ_x/σ_tot(WP70) = 1,10–1,48 dla D⁺ 3–5 i 5–8 — też nie pomaga.

## Rozstrzygnięcia względem wcześniejszych prac
- **LOWPT_calib — sprzeczności nie ma.** Tam: stosunek uzysków z niezależnych fitów względem najluźniejszego WP90 (S/B w oknie 0,028 przy 3–5 GeV) → ρ(70/90) ± 0,095 (427929), ± 0,164 (428770),
  czułość 2σ 21–30 %. Tu: odniesienie WP70 (S/B 0,19) i fit łączny kategorii rozłącznych → SF_rel(0,50 | 0,70) ± 0,053 (427929), ± 0,071 (428770) przy 3–5 GeV; dryf WPSCAN (dwa niezależne fity) ± 0,067 / 0,084;
  „3–8 %” z WPSCAN to średnie po trzech binach (T6). Różnica błędu ×2 pochodzi z punktu odniesienia, dalsze ~20 % z fitu łącznego. LOWPT sam podawał ρ(70/80) ± 0,057.
  Werdykt LOWPT („nieużyteczna jako kalibracja na 5 %”) zostaje i jest wzmocniony: precyzja 5 % jest osiągalna, ale **dokładność nie** — obciążenie fitu na MC −5…−8 % przy 3–5 GeV (−8,3 ± 4,1 % przy ε = 0,50),
  znane tylko z błędem 4–5 %; LOWPT §7 widział to samo dla pary 70/80 (+6,6 ± 5,1 %).
- **Etap 4, metoda R dla D⁰:** zagnieżdżone cięcia nie poprawiają (σ domknięcia 0,06–0,20 > 0,05 we wszystkich binach). D⁰ nadal bez kalibracji na próbce inkluzywnej; zostaje tag D*.
  Dla D⁺ nieco lepiej niż w etapie 4 (3–5 GeV: 0,041–0,048 przy ε = 0,50 / 0,40 wobec 0,054).
- **WPSCAN:** rekomendacja „kandydat ε ≈ 0,4–0,5 pod warunkiem pomiaru SF” — warunek spełnialny tylko przez **niezależny** SF (tag D* dla D⁰); SF z tych samych danych inkluzywnych zwraca błąd WP70.

## Zastrzeżenia
- Post hoc; 30 replik (błąd rms ±13 %); σ_tło = max z dwóch różnic rzędu — zaszumione i decydujące o połowie werdyktów (D⁺ 428770 5–8); szum ten nie jest w σ_g.
- σ_g: większy z dwóch błędów względnych (rms stosunku HESSE po replikach; jackknife stosunku rms) przeniesiony na g — implementacja ustalona przed odczytem tabel, po uruchomieniu fitów.
- Bramka G1 per komórka nie widzi wspólnego obciążenia −5 % (D⁺); średnie po binach i stosunek podwójny dane/MC są diagnostyką post hoc, nie kryterium.
- Domknięcie na MC: pseudo-dane mb_nd o statystyce ≈ jednego runu 427929 i składzie minimum-bias — dla 428770 (dżetowy) nie jest reprezentatywne; repliki MC po binach, nie po zdarzeniach.
- Szablony sygnału per kategoria z MC (wspólne dmu, sc); wariant jednego wspólnego szablonu zmienia SF_rel o ≤ 0,04 (D3). f_refl per kategoria stałe z MC (D⁰).
- SF_70 dla WP70 z tagu D* nie został tu policzony ani odnaleziony w bazie — scenariusze 5 / 10 / 15 %.
- Jedno ziarno drzewa; bin 2–3 GeV poza zakresem; klasy triggera tylko J12only / notJ12only.

## Rekomendacja
1. **Nominał zostaje WP70; linii „ciasny punkt + SF_rel z inkluzywnego piku” nie rozwijać** (wynik negatywny, wpis w bazie wiedzy).
2. Ciasny punkt ma sens tylko z SF mierzonym na **niezależnej** próbce: D⁰ ≥ 5 GeV — tag D* przy ε = 0,4–0,5 (do zrobienia osobno, wymaga joba); D⁰ 3–5 GeV i D⁺ — brak metody.
3. Błąd uzysku przy niskim pT ogranicza dziś nieznany SF_70 (D⁺, D⁰ 3–5 GeV), nie punkt pracy — tam warto inwestować najpierw.
4. W WPSCAN/RAPORT §3 i w notatce bazy dopisać korektę: dryf N/ε_MC drzewa D⁺ jest odtwarzany przez fit na MC (obciążenie kategorii luźnej), nie jest dowodem zależności SF od ε.

## Pliki
`KRYTERIA.md`, `RAPORT.md`, `wpt_fit.py`, `wpt_report.py`, `wpt_tables.md` (M1, D1–D6), `wpt_summary.json`, `wpt_mc.json`, `wpt_data.json`, logi `wpt_*.log`.
