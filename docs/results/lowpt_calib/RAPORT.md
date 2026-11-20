# LOWPT_calib — kalibracja ε_S klasyfikatora skanem cięcia na inkluzywnym D⁺ (i D⁰): studium wykonalności, 18.11.2026

Kryteria zapisane przed liczeniem: `KRYTERIA.md` (15:10). Tabele: `tables.md` (T1–T8), liczby: `verdict.json`, surowe: `out/scan.json`, `out/d0_qual.json`.
Klasyfikator H (`H_s12345`, XGB bez izolacji), dane inkluzywne bez AFP, runy osobno. Liczenie lokalne ~4 min (200 replik, bez redukcji).

## 0. Werdykt wg kryteriów: **NIEUŻYTECZNA** (jako kalibracja/pomiar SF w 3–5 GeV)
- **K_A niespełnione**: czułość 2σ na ρ(70/90) w 3–5 GeV, 427929 = **20,6 %** (próg 15 %).
- **K_B spełnione formalnie, bez mocy**: wstrzyknięte 0,90 → odzyskane 0,922 ± 0,120 (obciążenie +2,5 ± 13,4 %). Błąd testu jest 2,5× większy niż próg 5 %,
  więc „spełnione” znaczy tylko „nie wykryto obciążenia większego niż ~13 %”.
- Metoda nadaje się jako **test spójności / górne ograniczenie** na zależność SF od WP na poziomie 10–20 % (2σ), nie jako kalibracja na 5 %.

## 1. Co sprawdziłem w źródłach
- Pliki `cand/dlz_data_*.npz` i OOF MC: pT ∈ [3,00; 30) GeV — **binu 2–3 GeV nie ma** (ani w danych, ani w MC; H trenowany dla pT ≥ 3 GeV).
  Skąd wziąć (nie pobierałem): produkcja v2.0/v2.1 ma osobną część 2–3 GeV (pliki „d23”, GITLAB §2c/§2e); trzeba ekstrakcji kandydatów D⁺ z pT ≥ 2 GeV
  z ROOT v2.0 obu runów (EOS) + MC, oraz klasyfikatora wytrenowanego na 2–3 GeV (nasz H tam nie istnieje).
- Histogramy po cięciu H odtwarzają `K1c/out/prep.json` dokładnie (oba runy); kandydaci „luźni” w plikach = suma WP95 wszystkich wariantów, więc H ≥ WP90 jest kompletne.
- Potwierdzone: skan płaski w większości komórek; exp(pol2) zmienia ρ(70/90) względem pol3 o 0–7 % (najwięcej w 427929 3–8 GeV); obciążenie N_fit/N_true (MC, pol3) rośnie ku luźnym WP
  (5–8 GeV: 0,99 → 1,17 od WP50 do WP90, T4); 428770 nie jest płaski w 8–30 GeV.

## 2. Dane: SF względny ρ(a/b) = N/ε/N/ε, tło exp(pol3), bootstrap sparowany po zdarzeniach (200)
| run | pT | ρ(50/90) | ρ(70/90) | ρ(70/80) | nachylenie d ln y/dε | czułość 2σ: 70/90 · 50/90 · 70/80 |
|---|---|---|---|---|---|---|
| 427929 | 3–5 | 0,911 ± 0,100 | 0,916 ± 0,095 | 0,999 ± 0,057 | +0,13 ± 0,21 | **20,6 %** · 22,0 % · 11,5 % |
| 427929 | 5–8 | 1,184 ± 0,116 | 1,169 ± 0,100 | 0,909 ± 0,046 | −0,07 ± 0,21 | 17,2 % · 19,5 % · 10,2 % |
| 427929 | 8–30 | 0,916 ± 0,093 | 0,969 ± 0,087 | 1,060 ± 0,054 | +0,17 ± 0,23 | 18,0 % · 20,2 % · 10,2 % (udane repliki 68 %) |
| 428770 | 3–5 | 1,132 ± 0,173 | 1,098 ± 0,164 | 1,041 ± 0,077 | −0,25 ± 0,27 | 29,9 % · 30,6 % · 14,8 % |
| 428770 | 5–8 | 0,894 ± 0,127 | 0,918 ± 0,123 | 0,901 ± 0,056 | +0,30 ± 0,25 | 26,7 % · 28,4 % · 12,4 % |
| 428770 | 8–30 | 0,876 ± 0,048 | 0,939 ± 0,044 | 0,930 ± 0,025 | **+0,39 ± 0,12** | 9,4 % · 11,0 % · 5,3 % |
- Wszystkie ρ(x/90) w 3–5 GeV zgodne z 1 w 1σ. Odstępstwa > 2σ: 427929 5–8 GeV ρ(80/90) = 1,286 ± 0,103 (fluktuacja/obciążenie fitu przy WP90: y(WP90) niskie);
  428770 8–30 GeV nachylenie 3σ, ρ(50/80) = 0,867 ± 0,032 (znane: próbka jetowa, MC minbias jej nie opisuje).
- Błąd statystyczny MC ilorazu ε: 0,002 (pomijalny). pol4 zmienia ρ(70/90) o ≤ 0,02 w 3–5 GeV, do 0,04–0,06 w 5–8 i 8–30 GeV.
- 427929 8–30 GeV: 32 % replik odrzuconych (skala szerokości przy dolnej granicy 0,7) — komórka opisowa.
- Prosty wzór dla selekcji zagnieżdżonych (σ² = |σ_a² − σ_b²|, hesjany) zaniża błąd ρ(70/90) o 2–20 % (σ_boot/σ_wzór 1,02–1,20).

## 3. Kontrola z tagiem D* (D⁰; R3, fit kształtu Δm) — **inny gatunek, inny klasyfikator: tylko zgodność rzędu wielkości**
D* nie ma WP90, więc porównanie na parach względem WP80 (T6): 427929 — zgodne (|Δ| ≤ 1,0σ; np. 5–8 GeV 70/80: D⁺ 0,909 ± 0,046 vs D* 0,965 ± 0,034);
428770 — D⁺ systematycznie niżej o 2,2–2,8σ (8–30 GeV 70/80: 0,930 ± 0,025 vs 1,010 ± 0,015). Rząd wielkości odchyleń ten sam (≤ 10–15 %), ale w runie jetowym
skan D⁺ pokazuje spadek ku ciasnym WP, którego tag D* dla D⁰ nie widzi. Nie rozstrzygam, czy to różnica gatunków, czy obciążenie fitu D⁺.

## 4. Wstrzyknięcie na MC z pełnym tłem (mb, statystyka ≈ 427929; sygnał prawdy ze score ≥ WP70 z wagą 0,65–0,67 ⇒ prawdziwe ρ(70/90) = 0,900)
| pT | bez wstrzyknięcia: ρ_rec(70/90) (prawda 1) | wstrzyknięte: ρ_rec(70/90) (prawda 0,90) | obciążenie [%] | ρ(70/80): prawda → odzyskane | pull width (σ ze wzoru zagnieżdżonego) |
|---|---|---|---|---|---|
| 3–5 | 1,010 ± 0,093 | **0,922 ± 0,120** | +2,5 ± 13,4 | 0,942 → 1,039 ± 0,074 (+10,3 ± 7,8 %) | 1,12 |
| 5–8 | 0,946 ± 0,074 | 0,851 ± 0,094 | −5,5 ± 10,6 | 0,936 → 0,882 ± 0,059 (−5,8 ± 6,1 %) | 1,50 |
| 8–30 | 0,879 ± 0,074 | 0,757 ± 0,087 | −15,9 ± 9,6 | 0,938 → 0,942 ± 0,072 (+0,4 ± 7,5 %) | 10 (odporny 2,9) |
- W 3–5 GeV skan odzyskuje wstrzyknięcie w granicach błędu, ale błąd (12–13 %) nie pozwala potwierdzić progu 5 %.
- W 5–8 i 8–30 GeV ρ(x/90) jest ściągane w dół o 5–16 % także bez wstrzyknięcia (obciążenie fitu przy WP90: b = 1,17 ± 0,10 i 1,11 ± 0,09) —
  1,6–2,3σ, ten sam znak w obu wariantach (te same zdarzenia MC, więc to jeden test, nie dwa).
- Pull width > 1 wszędzie: błąd trzeba brać z bootstrapu; w 8–30 GeV wzór zagnieżdżony zawodzi (błędy obu WP prawie równe).

## 5. Ograniczenia metody
1. **Tylko zmiana względna.** Skan mierzy SF(a)/SF(b); wspólny czynnik (np. SF = 0,85 przy każdym WP) jest niewidoczny. Punkt odniesienia to najluźniejsze cięcie,
   gdzie S/B jest najmniejsze (427929 3–5 GeV: 0,028 przy WP90 vs 0,14 przy WP70) i błąd oraz obciążenie fitu największe — to ono ustala czułość.
2. **Model tła.** pol2 → pol3 przesuwa ρ(70/90) o 0,00–0,07; pol3 → pol4 o ≤ 0,02 (3–5 GeV) do 0,06. Obciążenie fitu zależne od WP (kombinatoryka pod pikiem)
   udaje złą ε; MC minbias wyznacza je tylko na 7–13 %.
3. **428770**: czułość 30 % w 3–5 GeV, a MC minbias nie opisuje próbki — korekty obciążenia z MC nie ma.
4. **2–3 GeV**: nie policzone (brak kandydatów). S/B spada z pT (WP70, 427929: 1,34 / 0,72 / 0,14 w 8–30 / 5–8 / 3–5 GeV), więc czułość będzie gorsza niż 20 % — domysł, nie pomiar.
5. **D⁰ inkluzywne** (lokalnie: `B_trig/runs/lxplus/bt_*.npz`, pola `all_*`, okno 1,75–2,00 GeV; `scripts/d0_qual.py`, bez bootstrapu, jakościowo, T7):
   - odbicie K↔π ma w MC rms ≈ 68 MeV w oknie 250 MeV i stosunek odbicie/sygnał zależny od WP (0,53 przy WP90 → 0,66 przy WP50; w ±40 MeV 0,19 → 0,26),
     czyli skan D⁰ wymaga, żeby MC opisywało **wydajność klasyfikatora dla odbicia** — dodatkowa, niekalibrowana wielkość;
   - błędy N: 15–60 % na komórkę, przy swobodnym R fit nieokreślony (R = 0…2,3 z błędem 0,2–4,8);
   - 3–5 GeV: ρ(70/90) = 0,66 ± 0,08 (427929) i 0,70 ± 0,13 (428770) niezależnie od R (0, R_MC, 1,3 R_MC) — skan nie jest płaski, a przyczyny (fit w wąskim oknie vs prawdziwa
     różnica ε) nie da się tymi danymi rozdzielić; D⁺ w tym samym binie daje 0,92 ± 0,10;
   - uwaga do sprawdzenia: wydajki D⁰ z tego fitu są podejrzanie małe względem liczby D⁰ z tagiem D* (427929 5–8 GeV: N/ε ≈ 550–770 ± 350 vs 462 z tagiem) — nie zweryfikowałem,
     czy `all_*` to pełna próbka. Znaczenie kodów prawdy 100/101 = odbicie przyjąłem z `btrig_extract.py` (`mark_swaps`).
   Wniosek: D⁰ inkluzywne nie nadaje się do skanu bez szerszego okna masy i osobnej walidacji odbicia; D⁰ kalibrować tagiem D*.

## 6. Rekomendacja dla analizy z CharmNet
- **Nie używać skanu jako kalibracji (SF ± błąd) w 3–5 GeV; używać jako systematyki ε_S „z testu stabilności”**: w 427929 skan D⁺ wyklucza zależność SF od WP
  większą niż ~20 % (2σ, para 70/90) lub ~12 % (para 70/80). Rozsądna niepewność do przypisania per bin (1σ, połowa czułości): **3–5 GeV ≈ 10 %** (427929; 15 % w 428770),
  5–8 GeV ≈ 9 % (13 %), 8–30 GeV ≈ 9 % (5 %) — w 5–8 i 8–30 GeV tag D* jest lepszy (4–12 %) i bezwzględny. **2–3 GeV: brak liczby; spodziewane > 10 %.**
- To jest niepewność na zmianę SF między WP, nie na samo SF. Wartość bezwzględną w 3–5 GeV trzeba przenieść z 5–8 GeV (tag D*) z dodatkową niepewnością ekstrapolacji w pT.
- Tło: **exp(pol3) wspólne dla wszystkich WP** (pol2 za sztywne; pol4 jako wariant). WP: progi CharmNet + co najmniej dwa luźniejsze i dwa ciaśniejsze; najluźniejszy punkt nie luźniej
  niż ε_S ≈ 80 % w 3–5 GeV (przy 90 % błąd odniesienia podwaja czułość). Błąd różnic wyłącznie z bootstrapu sparowanego po zdarzeniach.
- Najpierw 427929 (czysty MBTS, MC minbias go opisuje); 428770 tylko jako kontrola.
- Do powtórzenia z CharmNet: (i) score CharmNet per kandydat D⁺ w danych i MC (ONNX + YAML albo jego pliki), (ii) ε_S^MC(WP) i kształty sygnału per WP, (iii) fity + bootstrap,
  (iv) test wstrzyknięcia na MC z pełnym tłem przy S/B analizy z CharmNet, (v) bin 2–3 GeV z produkcji v2.0. CharmNet używa izolacji (ptcone), nasz H nie — S/B i korelacja score–tło będą inne,
  więc liczby stąd są wskazówką rzędu wielkości.

## 7. Post hoc (po obejrzeniu wyników, poza kryteriami)
- Para 70/80 zamiast 70/90 dałaby czułość 11,5 % < 15 % w 427929 3–5 GeV, ale na MC ρ(70/80) ma tam obciążenie +6,6 ± 5,1 % (bez wstrzyknięcia) i +10,3 ± 7,8 % (ze wstrzyknięciem) —
  próg 5 % nie byłby wykazany. Werdyktu nie zmieniam.
- K_B oceniłem względem σ_stat (rozrzut ρ_rec przy stałej ε); względem błędu obciążenia (prawda per replika) wynik ten sam.

## 8. Zastrzeżenia
H ≠ CharmNet; MC mb ma statystykę jednego runu 427929 (błąd obciążenia 7–13 %); wstrzyknięcie to skok wagi przy progu WP70, nie gładka deformacja score;
kształt sygnału w 3–5 GeV ustalony z MC (jak K1c); kontrola D* porównuje różne gatunki; D⁰ tylko jakościowo; liczby z literatury (RECON_lit_B3) niezweryfikowane i nieużyte.

## 9. Pliki
`KRYTERIA.md`, `RAPORT.md`, `tables.md`, `verdict.json`, `scripts/{scan.py, report.py, d0_qual.py}`, `out/{scan.json, d0_qual.json}`, `logs/`.
