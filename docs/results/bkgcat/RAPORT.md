# BKGCAT — rozkład tła MC na kategorie prawdy dla drzewa H (2026-11-20)

**Status: PRELIMINARY, w całości POST HOC względem FITBIAS** (pytania postawione po obejrzeniu wyników FITBIAS; ta sama realizacja mb_nd — lokalizuje efekty, nie potwierdza ich niezależnie). Tylko MC + inkluzywne piki danych; żadnych selekcji z tagiem AFP. FITBIAS, XSEC i nota bez zmian.

## 1. Co zrobiono

- Pakiet `charm_work/Q_lxplus/bkgcat` (`bkg_hist.py`, `bkg_merge.py`, `bkg_job.sh`, `bkg_common.sh`, `bkg.sub`): dla każdego z 78 plików mb_nd (Dpt3to150, v2.0) łączy pliki app drzewa H (`D0.xgb.H.sb.n.s1`, `Dplus.xgb.H.sb.n.s1`) z h5 i z ROOT-em (TruthDmesons) i wypełnia histogram masa (1 MeV) × przedział score (progi WPSCAN) × kategoria prawdy, per gatunek i bin pT.
- Przebieg: interaktywnie na lxplus (2 procesy, `/tmp/klyczak_bkgcat`), **nic nie wysłano na condora**. Wynik: `/eos/user/k/klyczak/charm_q/bkgcat/{parts/ (78), bkgcat.npz, bkg.log, DONE}`; kopia w `in/`.
- Kontrole: 19 220 986 zdarzeń, 8 320 523 wpisów TruthDmesons; dopasowanie app ↔ h5 ↔ ROOT 100 % (D⁰ 4 344 489, D⁺ 10 811 220 kandydatów; 16 duplikatów klucza w D⁰); łącze reco ↔ prawda 100 % dla Sig i DMM; **sumy po kategoriach zgodne co do zdarzenia z WPSCAN** we wszystkich komórkach. Wagi w mb_nd ≡ 1 → histogramy ważone = nieważone.
- Analiza lokalna: `skrypty/bc_br.py`, `bc_comp.py`, `bc_m8.py` (+ `bc_lib.py`); tabele w `tabele/`, liczby w `out/`.

Kategorie: Sig; odbicie K↔π (Sig & swap); DecayModeMismatch (DMM, per tryb prawdziwy); CharmMismatch (CM, per hadron); PartComb (PC, per rodzic i per NoB/SomeB/AllB — dwa równoległe podziały tej samej składowej); Comb: tory z jednego hadronu b / z lekkiego rezonansu / inne z b / bez powabu (OtherNoB).

## 2. Skład tła przy WP70 (okno sygnału; % tła = bez sygnału) — `tabele/sklad.md`

| składowa | D⁰ 3–5 | 5–8 | 8–30 | D⁺ 3–5 | 5–8 | 8–30 |
|---|---|---|---|---|---|---|
| Comb bez powabu (OtherNoB) | 83,4 | 60,4 | 38,0 | 51,9 | 19,9 | 7,6 |
| PartComb | 10,4 | 17,1 | 19,6 | 41,3 | 64,5 | 64,4 |
| odbicie K↔π | 2,4 | 10,8 | 21,4 | — | — | — |
| DMM | 1,3 | 5,3 | 11,0 | 0,4 | 1,6 | 2,7 |
| CM | 0,5 | 1,8 | 2,8 | 1,5 | 6,5 | 14,2 |
| Comb z b (trzy podklasy) | 2,0 | 4,6 | 7,2 | 4,9 | 7,5 | 11,0 |
| tło/sygnał w oknie | 14,8 | 3,2 | 1,5 | 5,7 | 1,4 | 0,7 |

- DMM D⁰ (okno, 3–5 GeV): K⁻π⁺π⁰ 325, π⁺π⁻π⁰ 121, półleptonowe „K⁻ + X” 79, K⁻K⁺ 33 (306 w pasmach — pik poniżej okna), π⁺π⁻ 24 (125 w pasmach — powyżej), DCS K⁺π⁻ 13. CM D⁰: D⁺ 198, D_s 29.
- CM D⁺: D_s 130, D⁰ 47, D*⁰ 40, D*⁺ 33. DMM D⁺ w oknie znikome (π⁺π⁺π⁻π⁰ 30); K⁻K⁺π⁺ i π⁺π⁺π⁻ leżą w pasmach (113 i 62).
- PartComb D⁺: rodzice rozłożeni szeroko (D⁺ 2291, D*⁰ 1412, D⁰ 1166, D*⁺ 1036, D_s 829, Λ_c 200; z b: AllB 1952).
- WP90, WP50, przedział [70, 50): `tabele/sklad.md`; kształty w masie (25 MeV): `tabele/ksztalt.md`.

## 3. Dołek D⁰ 3–5 GeV i nadwyżka D⁺: struktura czy fluktuacja — `tabele/spur.md`

Metoda: do histogramu samej składowej dodany sygnał Asimova (N_prawda × szablon fd), fit szablon sygnału + exp(Czebyszew), N_spur = N_fit − N_wstrzyknięte, w % N_prawda; błąd = rms 60 replik Poissona składowej. (Pierwsza wersja — fit samej składowej — opierała się o dolną granicę N w `fb_fit` i dawała fałszywie małe błędy HESSE; porzucona.) Rząd tła + 1 nie zmienia niczego (kolumna w tabeli).

**D⁰ 3–5 GeV, WP70:** reszta (wszystko − Sig − odbicie) **−22,1 ± 9,1 % (2,4σ)**; domknięcie M8 z szablonem wprost −22,0 ± 9,7 %. Rozkład:

| składowa | N_spur/N_prawda [%] | Z |
|---|---|---|
| Comb bez powabu | −16,7 ± 10,9 | −1,5 |
| PartComb | −1,3 ± 3,5 | −0,4 |
| DMM + CM | −3,1 ± 1,4 | −2,2 |
| Comb z b i reszta | ≈ −1 | — |

- Dołek siedzi w kombinatoryce **bez powabu** (83 % tła), nie w żadnej kategorii powabowej (PartComb −1,3 ± 3,5 %; podskładowe w `tabele/spur.md`).
- Znak nie jest stały w rozłącznych przedziałach score (głębokość względem liczby kandydatów składowej): [95, 90) +0,7 ± 0,5 %, [90, 70) +1,0 ± 0,4 %, [70, 50) −1,4 ± 0,8 %, WP50 −3,1 ± 1,8 % (χ² zgodności 10,9/3, p ≈ 0,012); stąd WP90 +21,6 ± 18,7 %, WP50 −11,3 ± 6,7 %. W 5–8 i 8–30 GeV wzór się nie powtarza (WP70: +0,4 ± 6,9 i −11,8 ± 6,7 %).
- Globalnie: 86 testów składowych niepikujących, |Z| > 2 w 4 (oczekiwane 3,9).
- **Wniosek: fluktuacja tła bez powabu na poziomie 2,4σ (lokalnie), bez mechanizmu fizycznego; nie jest to struktura konkretnej kategorii.** Zastrzeżenie: zmiana znaku z twardością cięcia (p ≈ 0,01) mogłaby też oznaczać rzeźbienie masy przez score — nierozstrzygalne na jednej realizacji mb_nd; rozstrzygnie niezależna próbka MC.
- Rzeczywista struktura: **DMM ma ujemny sygnał pozorny** — średnia po binach −4,6 ± 1,0 % (4,8σ; WP90 −3,6 ± 1,0, WP50 −3,3 ± 1,0). To efekt znany z FITBIAS (pasma zawyżone przez K⁻K⁺ i π⁺π⁻, okno zasilane K⁻π⁺π⁰) — obsługuje go szablon pikujący M8.

**D⁺, kategoria luźna [70, 50):** reszta +29,1 ± 15,1 % (3–5), średnia po binach +12,8 ± 5,4 % (2,4σ).

| składowa | 3–5 | 5–8 | 8–30 | średnia | Z |
|---|---|---|---|---|---|
| PartComb | +19,1 ± 10,7 | +7,3 ± 7,4 | +10,0 ± 8,4 | +10,7 ± 4,9 | +2,2 |
| CM | +4,6 ± 1,8 | +3,6 ± 2,3 | +5,4 ± 4,1 | +4,3 ± 1,3 | +3,2 |
| Comb bez powabu | +9,8 ± 12,2 | −0,6 ± 5,5 | −4,4 ± 3,3 | −2,8 ± 2,8 | −1,0 |
| DMM | −0,6 ± 0,7 | +1,3 ± 1,1 | −3,2 ± 2,1 | −0,3 ± 0,6 | −0,5 |

- **CM pikuje naprawdę** (głównie D_s; WP90 +2,5 ± 0,6 %, 4,4σ; WP70 +1,4 ± 0,6 %) — jest w szablonie M8.
- **PartComb: nadwyżka spójna co do znaku we wszystkich binach pT**, WP70 +4,0 ± 1,6 % (2,5σ), WP90 +5,0 ± 2,4 %, WP50 +0,8 ± 1,4 %. Nie wskazuje jednego rodzica (największe: rodzic D⁺ i Λ_c, po ~2σ). **Wskazanie 2–2,5σ, nie ustalenie**; jeśli prawdziwe, jest to składowa, której M8 nie modeluje (PartComb nie wchodzi do szablonu pikującego).
- Kombinatoryka bez powabu w D⁺ jest płaska (zgodna z zerem).

## 4. Szablon składowych pikujących wprost z drzewa H — `tabele/m8.md`, `out/peak_templates_H.npz`

- Plik: klucze `<sp>|<bin>|<kat>|{h5,h5_pdg,kap,kap_pdg,nsig,npk}` + `edges|<sp>`; kategorie WP95…WP10, B50, B40, B35, L90-70. Wczytanie do M8: `bc_m8.patch(k, sp, bn, elo, ehi, mode, r)` podmienia `peakT`/`kap` w słowniku z `fb_var.wcat` (mode `ref` | `dir` | `pdg`); FITBIAS bez zmian. Odtworzenie zapisanych uzysków FITBIAS M8 co do cyfry.
- κ wprost vs interpolowane (WP70): D⁰ 1,231 ± 0,018 / 1,108 ± 0,024 / 0,903 ± 0,033 wobec 1,215 / 1,089 / 0,922; D⁺ 0,527 / 0,534 / 0,507 wobec 0,534 / 0,546 / 0,514. Kształt szablonu T vs H zgodny (χ² 5–14 na 13). **Interpolacja κ była poprawna w granicach 1σ.**
- Domknięcie mb_nd: M8 dir ≈ M8 ref (różnice ≤ 2 pkt proc. przy WP70).

Uzyski D⁰ w danych, WP70, M8:

| run | bin | FITBIAS M8 | dir, r = 1 | Δ łącznie | w tym sam szablon (r = 1) | PDG − dir |
|---|---|---|---|---|---|---|
| 427929 | 3–5 | 4105 ± 366 | 4059 ± 357 | −1,1 % | −0,5 % | −0,1 % |
| 428770 | 3–5 | 6403 ± 639 | 6322 ± 613 | −1,3 % | −0,2 % | −0,3 % |
| 427929 | 5–8 | 1530 ± 126 | 1536 ± 124 | +0,4 % | −0,8 % | −0,2 % |
| 428770 | 5–8 | 4846 ± 302 | 4829 ± 302 | −0,3 % | −0,1 % | −0,2 % |
| 427929 | 8–30 | 681 ± 54 | 700 ± 54 | +2,8 % | −0,7 % | +0,1 % |
| 428770 | 8–30 | 6829 ± 236 | 6987 ± 238 | +2,3 % | −5,6 % | +0,2 % |

- Wszystkie zmiany poniżej błędu statystycznego (0,1–0,7σ). Przy 8–30 GeV dominują dwa przeciwne efekty: szablon wprost (−0,7 / −5,6 %) i zdjęcie skali κ × 1,364 (+3,5 / +8,4 %).
- 428770 8–30 GeV nadal odrzuca model (p = 3·10⁻⁵ wobec 1·10⁻⁷) — szablon wprost tego nie naprawia.
- D⁺: zmiany ≤ 0,7 %. Wyjątek: 427929 8–30 po przeważeniu PDG +5,3 % (661 wobec 627) — skok fitu, nie efekt wag (κ zmienia się o 2 %); traktować jako niestabilność.

## 5. BR w MC wobec PDG 2026 — `tabele/br.md`, `out/br.json`

PDG: listingi 2026 pobrane z otwartej strony PDG (numery Γ w tabeli). Lista córek w TruthDmesons zawiera tylko π±, K±, π⁰ — tryb = topologia widzialna + klasa brakującej masy (nic / radiacyjny / K⁰ / η / inne). **D_s i Λ_c nie są zapisywane w TruthDmesons — tablicy dla nich nie da się zbudować z tej produkcji.**

- Większość trybów w ±5 %: K⁻π⁺ 3,974 vs 3,936 % (1,010), K⁻π⁺π⁺ 9,404 vs 9,31 % (1,010), K⁻π⁺π⁰ 1,018, K⁻K⁺ 1,020, π⁺π⁻ 1,019, D*⁺→D⁰π⁺ 67,72 vs 67,7 %.
- Największe odchylenia (MC/PDG): **D⁺→K⁺K⁻π⁺π⁰ 4,53** (3,00 vs 0,662 %), D⁺→π⁺π⁺π⁻ 1,47, D⁰→K⁺K⁻π⁰ 1,41, D⁰→K⁺K⁻π⁺π⁻ 1,28, D⁰→K⁻π⁺2π⁰ 1,26, D⁺→2π⁺π⁻π⁰ 1,24, D⁺→K⁻2π⁺π⁰ 1,11; za mało: D⁰→K⁻π⁺η(γγ) 0,54, D⁺→K⁰2π⁺π⁻ 0,70.
- Przeważenie do PDG: κ D⁰ −1,4 … −1,7 %, D⁺ −1,3 … −1,8 %; uzyski w danych ≤ 0,3 % (pkt 4). **Dla szablonu pikującego bez znaczenia.**

### Znalezisko uboczne: „decayType” sygnału jest topologiczny

Prawdziwy `decayType` 10 (D⁰ „Kπ”) i 4 (D⁺ „Kππ”) oznacza „K∓π±[π±] i żadnych innych π±, K±, π⁰” — zawiera rozpady z niewidocznym K⁰, η→γγ i innymi neutralnymi:

| | ekskl. + radiacyjny | + K⁰ | + η | + inne | „cały decayType”/PDG |
|---|---|---|---|---|---|
| D⁰, decayType 10 | 73,62 ± 0,13 % | 11,82 % | 7,44 % | 7,11 % | 1,371 |
| D⁺, decayType 4 | 94,01 ± 0,07 % | 5,63 % | 0,00 % | 0,37 % | 1,075 |

1. **Czynniki „BR_MC/PDG” 1,364 (D⁰) i 1,068 (D⁺) używane jako skala κ w danych w FITBIAS są artefaktem tej definicji** — BR trybu sygnału w MC zgadza się z PDG w 1 %. Właściwa skala to r = 1 (zastosowana w pkt 4; efekt na uzyski ≤ 3,5 %, przy 428770 8–30 +8 %).
2. **XSEC (do weryfikacji w studium XSEC; nic tam nie zmieniono):** jeżeli mianownik C liczy N_gen po `decayType`, to zawiera rozpady nierekonstruowalne w oknie masy, a licznik (Sig w zakresie mas) jest w 100 % ekskluzywny (0 kandydatów z brakującą masą na 3653 / 1919 / 848). Wtedy σ_popr/σ_XSEC = **0,734 (D⁰), 0,939 (D⁺)**.

## 6. Tło MC w 2–3 GeV

W naszych plikach nie ma: mb_nd to wyłącznie produkcja Dpt3to150 (minimalne pT kandydata 3,000 GeV). Produkcja Dpt2to3 dla mb_nd istnieje według dokumentacji produkcji v2.0 (80 plików, 72,1 mln kandydatów, na NERSC), nie była na naszej liście zbiorów. Kontener rucio `user.malewick.mc23.minbiasND.CharmAnalysis.Dpt2to3GeV.v2.0.Oct28_EXT0` istnieje, ale dziś ma 0 plików (podobnie kontener Dpt3to150 mb_nd). **Potrzebna kopia z NERSC albo przywrócenie repliki** — pakiet bkgcat przyjmie ją po dodaniu progu etapu 4 DNET dla tego binu.

## 7. Zastrzeżenia

- Jedna realizacja mb_nd, ta sama co w FITBIAS: istotności z pkt 3 są lokalne i nie są niezależnym potwierdzeniem.
- Błędy sygnału pozornego z 60 (grupy) / 30 (podskładowe) replik — dokładność samego błędu ~10–13 %.
- Błędy uzysków w danych: HESSE (FITBIAS: rzeczywisty rozrzut ~1,3×).
- Porównania z PDG dla trybów z K⁰ zakładają K⁰_L = K⁰_S; klasa braku wymaga |η| < 2,5 (precyzja float).
- `bkg.sub` (wariant condorowy) napisany, nieuruchamiany.
- D_s, Λ_c: brak prawdy generatora w ntuplach — nie da się sprawdzić ich BR ani udziału w PartComb inaczej niż przez rodzica kandydata.

## 8. Pozostałości

- lxplus `/tmp/klyczak_bkgcat` (węzeł lxplus9103; 11 MB: parts/, bkgcat.npz, logi, probe1–4.py, partial.npz).
- AFS `~/private/charm_q/Q_lxplus/bkgcat` (5 plików pakietu).
- EOS `/eos/user/k/klyczak/charm_q/bkgcat/` (parts/ 78 plików, bkgcat.npz, bkg.log, DONE; ~10 MB).
- Lokalnie: katalog tymczasowy (listingi PDG pdf/txt, częściowe scalenie, skrypty próbne), `skrypty/__pycache__`.
