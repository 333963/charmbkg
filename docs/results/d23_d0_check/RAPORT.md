# Weryfikacja piku D⁰→Kπ w binie pT 2–3 GeV (etap 4 DNET) — 19.11.2026

Niezależna próba obalenia wyniku `DNET/runs/lxplus4/dnet4_d23.json` (WERDYKT_etap4 §5). Dane inkluzywne, bez AFP.
Wszystko tutaj jest **post hoc** (po obejrzeniu wyniku etapu 4); progów i kryteriów nie zapisywano z góry.

## 0. Werdykt

**Pik jest prawdziwy** (sygnał D⁰ w obu runach i dla obu klasyfikatorów); hipoteza „pik w całości wyrzeźbiony z tła” jest odrzucona kilkoma
niezależnymi testami. **Nie jest wykluczona domieszka nie-sygnałowa rzędu ≤ 20–30 % uzysku** (testy 1a, „pik ∝ tło” i 5 mają taką czułość), a uzysk
N/ε pozostaje niekalibrowany. Istotność zależy od modelu tła: nominalnie Z_LR = 8,6 / 9,1 (427929 drzewo / sieć) i 7,8 / 6,8 (428770), ale przy
exp(pol4) na wąskim zakresie 1,75–2,00 spada do 5,0 / 5,0 i 3,6 / 3,2 — w 428770 kryterium „Z ≥ 4 dla sąsiednich rzędów” jest spełnione tylko dla ±1.

| test | wynik | odczyt |
|---|---|---|
| 0 odtworzenie | N, σ, Z etapu 4 odtworzone co do cyfry (448,2 ± 54; 860,4; 501,3; 391,1); druga, niezależna implementacja fitu: 448,7 / 860,4 / 501,8 / 391,2 | kod etapu 4 poprawny |
| 1a tło MC z prawdą | tylko mb_nd **3–5 GeV** (w 2–3 GeV MC tła nie ma): N_fałsz·f_okno/B_okno = +0,09 ± 0,09 (drzewo WP10), +0,02 ± 0,05 (drzewo WP20), +0,01 ± 0,07 (sieć WP20), −0,17 ± 0,11 (sieć WP10); pik w danych wymaga S/B = 0,68 / 0,33 | „cały pik z tła” odrzucone (~7σ dla drzewa WP10); w skali danych 427929 drzewo WP10: N_fałsz = +57 ± 57 z 448 |
| 1b zły znak | brak próbki bez sygnału (Kπ neutralne; par jednoimiennych w tabelach nie ma) | niesprawdzalne |
| 1c ε_B(m) w danych | ln ε_B liniowy w masie: nachylenie −0,55 ± 0,03 / −0,48 ± 0,02 na 100 MeV (drzewo), −0,21 / −0,27 (sieć); krzywizna zgodna z 0 (+0,03 ± 0,04; −0,03 ± 0,03; +0,03 ± 0,02; +0,05 ± 0,04) | brak „namiotu” wokół m₀; drzewo silnie zależy od masy, ale monotonicznie |
| 1d drzewo (tło z pasm) vs sieć (tło MC + DisCo) | kategorie rozłączne: R = y_sieć/y_drzewo = 0,87 ± 0,10 (427929), 0,84 ± 0,14 (428770), średnio 0,86 ± 0,08 | zgodne w 1,7σ; R < 1 jak dla D⁺ przy pT ≥ 5 GeV (0,87–0,94), więc to raczej SF sieci niż rzeźbienie drzewa |
| pik ∝ tło? | N427/N428 przy WP10 = 0,89 ± 0,15 (drzewo), 0,98 ± 0,19 (sieć); tło B427/B428 = 0,50 / 0,71; D⁺: 1,13 ± 0,16 | „pik skaluje się z tłem” odrzucone na 2,7σ / 1,5σ; składowa ∝ tło: +92 [0; +164] z 449 (drzewo) |
| 2 odbicie K↔π | model ma odbicie (kształt i f_refl z MC, ustalone): f_refl = 0,92 (drzewo WP10), 0,71 (sieć WP20); obie hipotezy trzymane; 1,40 / 1,25 kand. na zdarzenie; w oknie 75 % (drzewo) / 52–55 % (sieć) kandydatów ma bliźniaka, który też przechodzi; N przy f_refl = 0 / ×0,5 / ×1,5 / wolne: 453 / 454 / 437 / 444 (nom. 449) | N niewrażliwe (≤ 3 %); błąd z bootstrapu po zdarzeniach = HESSE (51,5 vs 49,1) |
| 3 parametry piku | w etapie 4 **ustalone z MC**; wolne (błąd z bootstrapu): Δμ = −2,6 ± 2,0; +2,8 ± 2,1; −2,3 ± 2,2; −2,0 ± 2,3 MeV; skala szerokości 1,00 ± 0,05; 0,98 ± 0,04; 0,97 ± 0,04; 1,00 ± 0,03; Gauss: σ_dane/σ_MC = 0,87 ± 0,20; 0,99 ± 0,22; 1,10 ± 0,26; 1,14 ± 0,27 | zgodne z MC i z binem 3–5 GeV (Δμ −1,7…+1,2 MeV, sc 0,92–1,07 dla WP50–20) |
| 4 stabilność | wycinki rozłączne ε: χ²/ndf = 8,1/4 (p 0,09), 3,4/4, 4,4/4, 0,5/4; ⟨y⟩ = 4416 ± 410 / 4266 ± 403 (427929 drzewo / sieć), 4867 ± 534 / 3954 ± 472 (428770) | stabilne; jedno odchylenie: 427929 drzewo, wycinek ε 0,3–0,2 daje y = 744 ± 1441 (p = 0,02 dla trzech ciasnych) |
| 5 wiarygodność | y(D⁰)/y(D⁺) = 1,30 ± 0,16 / 1,23 ± 0,15 (427929); oczekiwane 0,82–1,03 | rząd wielkości zgodny; D⁰ o ~1,3–1,5× za dużo (≈ 2σ stat), tak samo jak w 3–5 GeV |

## 1. Czym różniła się selekcja z 18.11 (D23 §5, nota Sec. 6.12: Z ≤ 3,1)

Tamten test D⁰ był **wyłącznie bez klasyfikatora** (histogramy pełnej statystyki, 34,7 M / 63,7 M kandydatów, L_xy ≥ 0,10 mm, tło pol3–5, bez odbicia
w modelu); punkty pracy H sprawdzano tylko dla D⁺. Etap 4: L_xy ≥ 0,25 mm, 15,7 M / 29,5 M kandydatów, cięcie ε_MC = 0,1–0,2 (tło × 1/1500–1/4000).
Spójność: tamte N = 13,4–15,2 tys. ± 4,9–6,9 tys. (427929) zawierają odbicie i luźniejszy L_xy; z obecnego y ≈ 4,4 tys. przy ε = 1 oczekuje się
Z ≈ 1–1,5. **Wyniki nie są sprzeczne** — tam po prostu nie było czułości. Przy luźnych cięciach obecnych modeli też nic nie widać (wycinki ε 0,5–0,3: S/B 0,004–0,05).

## 2. Szczegóły

**Dane wejściowe.** Etap 4 zapisał score per kandydat (`app/d23/*.npz`, 2,6 GB; `app/d150` dla pT ≥ 3). Zredukowane na lxplus (`scripts/reduce_lx.py`,
6 przebiegów interaktywnych, każdy < 2 min) do `in/*.npz`: histogramy (m × score) wszystkich kandydatów + wiersze po luźnym cięciu (WP50) + score
bliźniaka K↔π. Liczba kandydatów zgodna z JSON (15 698 388 / 29 495 643).

**Dlaczego podejrzenie było zasadne.** Zestaw cech H zawiera log pT, |η|, log p_K, log p_π,min, |cos θ*| — z nich masę da się odtworzyć, a drzewo
`…sb…` uczy się sygnał MC kontra pasma boczne danych. ε_B(m) drzewa rzeczywiście zmienia się 5,5× między 1,70 a 2,05 GeV (sieć 2×). Zależność jest
jednak monotoniczna (bez krzywizny), a tło MC z prawdą w 3–5 GeV po tym samym cięciu nie ma garbu (`fig/kontrola_tloMC_3-5_i_epsB_2-3.png`).
Dodatkowo drzewo trenowano na pT ≥ 3 GeV — wyuczona funkcja masy nie ekstrapolowałaby się poniżej 3 GeV na pik dokładnie w m₀ o szerokości MC.

**Test 1a — ograniczenie.** fd_nd w 2–3 GeV zawiera tylko sygnał (truth 0/1), mb_nd dla tego binu nie istnieje. Użyto mb_nd 3–5 GeV z progami WP z fd_nd
3–5 (liczby w tabeli; `out/s02.log`). Domknięcie z pełnym tłem w mb_nd 3–5: N_fit/N_prawda = 0,99 ± 0,08 (drzewo WP10), 0,96 ± 0,06 (WP20),
0,94 ± 0,07 (sieć WP20), 0,88 ± 0,09 (sieć WP10). Uwaga: S/B w mb_nd jest zawyżone (BR(D⁰→Kπ) ×1,364), więc domknięcie jest łatwiejsze niż w danych.
Kategoria truth 2–7 (inne rozpady powabu) nie daje piku: N_fałsz = −6 ± 10 (drzewo WP10).

**Dane 3–5 GeV przy ciasnym cięciu** (kontrola, czy ciasne cięcie „dodaje” pik): y nie rośnie przy zacieśnianiu — 427929 drzewo 5650 ± 392 (WP50) →
5250 ± 382 (WP10); sieć 5419 → 3918 ± 365; 428770 drzewo 8149 → 6767; sieć 6654 → 5446. Spadek sieci przy WP10 (−25 %) = ε_dane < ε_MC, nie rzeźbienie.
y(2–3)/y(3–5) = 0,77–0,88 (427929), 0,61–0,69 (428770), wobec D⁺ 0,79 ± 0,08 i 0,52 ± 0,07.

**Warianty fitu** (N; N/σ; Z_LR), 427929 drzewo WP10: exp1 431; 12,8 · exp3 446; 9,0 · exp4 417 ± 71; 5,9; 5,0 · pol2 433 · pol3 451 · zakres 1,70–2,05:
exp2 460 ± 38, exp3 450 ± 41, exp4 447 ± 52 (Z_LR 7,9) · zakres 1,78–1,96 exp2: 417 ± 65 (5,7). Pozostałe trzy punkty: `out/s03.log`. Rozrzut N ≤ 7 %
(427929), do −25 % przy exp4 w 428770 (378 ± 97, 302 ± 88). Zakres nominalny 1,75–2,00 jest wąski (±115/135 MeV); szerszy 1,70–2,05 daje mniejsze błędy i te same N.
Z etapu 4 dla rzędu −1 (12,8) jest zawyżone: to N/σ(HESSE) ze zbyt sztywnym tłem.

**Wielokrotność.** Dla D⁰ 2–3 obejrzano 10 fitów (2 runy × 2 klasyfikatory × WP wybrany na ślepo z pasm + sąsiedzi); pozycja piku ustalona.
Z_LR 8,57 → 8,3 globalnie; 9,09 → 8,8; 7,81 → 7,5; 6,81 → 6,5. Nieistotne. Oba runy łącznie (dane niezależne): drzewo 11,6.

**Parametry piku — uwaga do etapu 4.** Błędy HESSE dla dmu/sc w `dnet4_d23.json` (np. sc 0,957 ± 0,017) są zaniżone ~2× (bootstrap: ± 2,0 MeV, ± 0,05);
wiarygodność w (dmu, sc) ma lokalne minima (szablon z dystrybuanty). Profil po sc (427929 drzewo): Δ(−2lnL) = 2,9 / 0 / 3,1 dla sc = 0,8 / 1,0 / 1,2.
σ_eff sygnału MC: 22,4–23,5 MeV (2–3 GeV, WP10–20) vs 22,3–23,8 (3–5 GeV).

**Test 5 — założenia.** y(D⁺) = 3441 ± 260 / 3055 ± 363 (D23, H WP20). Stosunek wydajności preselekcji z fd_nd (te same 15 plików): n(D⁰)/n(D⁺) =
51 920 / 27 709 = 1,87; przy założeniu, że próbka wymusza rozpad każdego D (gen. D⁰/D⁺ = 1/0,521, Pythia A3) daje ε⁰/ε⁺ = 0,98 — **założenie
niesprawdzone**. Wtedy skorygowane D⁺/D⁰ = 0,32 ± 0,04 / 0,34 ± 0,04 (427929), 0,26–0,32 (428770) wobec 0,36 ± 0,15 w 3–5 GeV (XSEC) i 0,40–0,47
(świat). BR w mb_nd (×1,364) tu nie wchodzi (fd_nd, stosunki). ε przy WP10–20 poza zakresem treningu niekalibrowane dla obu gatunków.

## 3. Czego nie dało się sprawdzić
- Tła z prawdą w samym binie 2–3 GeV (brak MC minbias) — test 1a jest ekstrapolacją z 3–5 GeV.
- Próbki bez sygnału w danych (zły znak, pary jednoimienne) — nie ma w tabelach.
- Zależności od cech (L_xy, dE/dx, czas życia) — w plikach `app` są tylko m, pT, score; tabele cech były tymczasowe.
- Tagu D* w 2–3 GeV (D23: niewidoczny) → brak kalibracji ε; N/ε to nie przekrój.
- Normalizacji widma MC 2–3 / 3–5 (różne produkcje fd_nd: 15 i 24 pliki).

## 4. Pliki
`scripts/`: `reduce_lx.py` (lxplus), `vlib.py`, `s01_repro.py`, `s02_mc_sculpt.py`, `s03_data.py`, `s04_more.py`, `s05_plots.py`; `out/s0{2,3,4}.{log,json}`;
`fig/masa_D0_2-3_{427929,428770}_{drzewo,siec}_WP*.png`, `fig/kontrola_tloMC_3-5_i_epsB_2-3.png`; `in/*.npz` (50 MB).
Na lxplus pozostał katalog `/tmp/klyczak_d23ver` (węzeł lxplus9103, 7 plików npz + skrypt + log, ~50 MB) — do usunięcia. Import `dn4_ana.py` tylko do odczytu.
