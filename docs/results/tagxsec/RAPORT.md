# TAGXSEC — fiducjalny przekrój czynny D⁰ + proton w AFP i stosunek σ_tag/σ_incl, run 427929 (projekt B), 20.11.2026

Kryteria zapisane przed liczeniem: `KRYTERIA.md` (00:10 CET); **dopisek §10 (~00:25) powstał po pierwszym przeliczeniu** — zmiana definicji mianownika C na prawdziwy tryb rozpadu po uwadze z BKGCAT / CVER.
**Wszystko jest post hoc względem odczytu TAGGED (19.11, 23:51)** — N_tag i jego zakres 30–55 były znane.
Żadnych nowych fitów pików w zdarzeniach z tagiem; N_tag, N_acc, fity per strona wzięte z `TAGGED/out/unblind.json`. Wynik lokalny; nota, TAGGED, XSEC, AFP_akceptancja, `pipeline.yaml` nietknięte; pętla ORCH nie restartowana; bez condora.
Tabele (generowane): `tabele/wyniki_excl.md` (mianownik C po prawdziwym trybie — obowiązuje dla ρ), `tabele/wyniki_xsec.md` (mianownik po `decayType`, jak XSEC). Liczby: `out/tagxsec_{excl,xsec}.json`. Akapit do noty: `akapit_do_noty.md`.

## 0. Wynik w skrócie

1. **Wynik główny — stosunek: R^fid = σ_tag^fid/σ_incl (D⁰, 3–30 GeV) = 0,99 −0,34 +0,35 (stat) ± 0,28 (syst) %** (nominał, N_tag = 54,7 metodą TAGGED); **wariant „okno − pasma” (N = 29,9): 0,54 ± 0,23 ± 0,09 %.**
   Różnica nominał → wariant (−46 %) to jawna, jednostronna systematyka modelu szerokich składowych; nie jest w ±0,28. R^fid nie zależy od definicji mianownika C (0,99 % w obu).
2. **σ_tag^fid(D⁰ + anty-D⁰; 3 ≤ pT < 30 GeV, |η| < 2,2; proton 0,03 ≤ ξ ≤ 0,10 po stronie A lub C, bez cięcia t) = 13,6 −4,7 +4,8 (stat) ± 4,6 (syst) µb** (nominał); **wariant „okno − pasma”: 7,4 ± 3,1 ± 1,9 µb.**
   C z prawdziwym trybem K∓π± w liczniku i mianowniku — zgodnie z werdyktem CVER (`CVER_2026-11-20/RAPORT.md`: σ = σ_XSEC × 0,734 ± 0,002 dla D⁰, × 0,939 ± 0,001 dla D⁺, niezależnie od pT). σ_incl(3–30) użyte w R^fid: 1375 µb.
   Własne zliczenie per D daje ten sam czynnik: 0,734 / 0,734 / 0,737 per bin w `mb_nd` (D⁺ 0,939); łączny czynnik na σ_tag^fid 0,730 (z małą zmianą ρ).
   Przypis: przy niepoprawionym C z XSEC (mianownik po `decayType`) byłoby 18,6 −6,4 +6,6 ± 6,5 µb (wariant 10,1 ± 4,3 ± 2,6 µb) — nieobowiązujące, tylko do śledzenia zmian (`tabele/wyniki_xsec.md`).
3. **C_D(SD)/C_D(ND)** (licznik i mianownik po prawdziwym trybie): przyjęte **ρ = 0,96 ± 0,12** (D⁰, 3–30 GeV). Z krotności (C_D(n) z `mb_nd` przeważone do nSelBL kandydatów tagowanych w danych) 0,962 ± 0,012; wprost z `mb_sd` 1,08 ± 0,09 (98 zrekonstruowanych D⁰), z protonem w oknie 1,17 ± 0,19 (25).
   Błąd przyjętego ρ to w całości niedomknięcie metody na MC SD (12 ± 10 %), czyli statystyka `mb_sd`. D⁺: wprost 0,79 ± 0,09, z krotności 0,98 (2,1σ — napięcie wg kryteriów, D⁺ opisowo); oba gatunki razem wprost 0,93 ± 0,07 — zgodne z 0,96.
4. **Zastrzeżenie zakresu:** 64 % mianownika C̄ i 71 % widma pochodzi z binu 3–5 GeV, gdzie SF klasyfikatora jest scenariuszem (1 ± 10 %), a domknięcie fitu inkluzywnego ma 25 %. σ_tag^fid w µb jest więc „pomiarem przy scenariuszu”; R^fid jest od tego w większości wolny.
   σ_tag dla ≥ 5 GeV nie podaję (N_tag(5–30) z TAGGED to artefakt granicy N_B → 0).
5. **Porównania (R, poziom cząstek):** Pythia miękka A3, 3–30 GeV: 0,35 % (obwiednia 0,14–1,40) → **zgodny** dla nominału (+1,4σ od centralnej) i wariantu (+0,8σ). Scenariusz twardy 0,66 % (0,58–1,41): zgodny (+0,7σ / −0,5σ).
   Literatura (Łuszczak–Maciuła–Szczurek, Siddikov–Schmidt) po przeliczeniu okna ξ: R ≈ 0,4–1,2 % — ten sam rząd. Bezwzględnie: Pythia 6,6 µb (2,6–26; z jej własnym σ_incl) wobec 13,6 (7,4) µb, LMS 1,1–2,5 µb.
   Przewidywany uzysk zrekonstruowany z Pythii przy poprawionym C: ~27 D⁰ (11–105) wobec 54,7 (29,9) zmierzonych; dawne ~20 liczone z C z XSEC było zaniżone × 0,73.
6. **D⁺ (opisowo, kontrolki niespełnione):** R^fid = 0,55 −0,20 +0,22 ± 0,20 % (wariant 0,60 ± 0,23 ± 0,11 %); σ_tag^fid = 4,0 −1,5 +1,6 ± 1,5 µb (wariant 4,3 ± 1,7 ± 0,9 µb; przypis: 4,2 przy C z XSEC). R^fid(D⁺) zgadza się z wariantem D⁰ i leży 1σ poniżej nominału D⁰.

## 1. Definicja i wzór

Poziom cząstek: D⁰ + anty-D⁰, prompt + nieprompt, 3 ≤ pT < 30 GeV, |η| < 2,2, w zdarzeniu z protonem 0,030 ≤ ξ ≤ 0,100 po stronie A lub C (dowolny proces, bez cięcia t, zdarzenie dwustronne liczone raz), √s = 13,6 TeV; liczone mezony.
Cięcie ξ_cal i weto wierzchołków pile-up są **poprawione wydajnością**, nie są częścią definicji (uzasadnienie: KRYTERIA §1).

σ_tag^fid = (N_tag − N_acc)/(C̄_tag · BR · L · ⟨A·ε_AFP⟩ · ε_trig · ε_ξcal · ε_npu), C̄_tag = Σ_b g_b · C_b^ND · SF_b · ρ_b; R^fid = σ_tag^fid/Σ_b σ_incl,b.

**Definicja C (po dopisku §10):** licznik = wiersze app z kategorią prawdy 0/1, `ok`, D⁰: `swap` ≠ 1, masa w zakresie szablonu, WP70 (jak XSEC; odtwarza 3653 / 1919 / 848 co do zliczenia).
Mianownik = `TruthDmesons` z `decayType` 10 (D⁰) / 4 (D⁺) **∧ klasa brakującej masy 0 lub 1** (m_D − m_wid < 30 MeV albo m_brak < 100 MeV, liczone z czteropędów córek; kod z `bkgcat/bkg_hist.py`), |η| < 2,2, bin pT prawdy — per D, w `mb_sd` i `mb_nd` osobno.
Licznika nie sprawdzałem osobno pod kątem trybu (BKGCAT: 0 kandydatów z brakującą masą na 3653 / 1919 / 848 w `mb_nd`; dla `mb_sd` przyjęte).

| czynnik | mianownik po prawdziwym trybie | mianownik po `decayType` (XSEC) | źródło |
|---|---|---|---|
| N_tag − N_acc | 54,7 − 0,20 (wariant 29,9 − 0,20) | to samo | TAGGED |
| C_b^ND (3–5 / 5–8 / 8–30) | 0,223 / 0,300 / 0,396 | 0,164 / 0,220 / 0,292 | `mb_nd`; mianowniki 16 366 / 6394 / 2143 wobec 22 306 / 8706 / 2909 |
| SF_b | 1 ± 0,10 (scen.) / 1,005 ± 0,126 / 0,965 ± 0,097 | to samo | XSEC |
| ρ_b (z krotności, dane); łącznie | 0,963 / 0,950 / 0,984; 0,962 ± 0,120 | 0,962 / 0,939 / 0,966; 0,956 ± 0,137 | §2 |
| g_b (MC SD z protonem w oknie; 2191 prawdziwych D⁰, wszystkie mody) | 0,713 / 0,227 / 0,060 | to samo | `mb_sd`; inkluzywne zmierzone: 0,73 / 0,20 / 0,07 |
| **C̄_tag** | **0,241** | 0,176 | udziały binów 64 / 27 / 9 % |
| BR, L | 3,951 %, 0,4648 nb⁻¹ | | XSEC, LUMI |
| ⟨A·ε_AFP⟩ z wagami stron | 0,922 (A: 0,919, C: 0,924; N_A = 17,5, N_C = 36,3 z istniejących fitów K7 TAGGED) | | różnica wobec średniej 0,9215 pomijalna |
| ε_trig, ε_ξcal, ε_npu | 1,000, 1,000, 0,985 | | AFP_akceptancja, nota §6.11, §6.7 |
| mianownik wzoru | 4,020 µb⁻¹ | 2,934 µb⁻¹ | |

Kontrola: f_true(TAGGED) 0,878 % / (0,922 · 0,985 · 0,962) = 1,00 % wobec R^fid = 0,99 % (reszta = widmo pT).

## 2. C_D(SD)/C_D(ND)

**Co zrobiono.** Score drzewa H dla `mb_sd` nie istniał (EOS `dnet4/app/d150`: tylko data, fd_nd, mb_nd). Zastosowany istniejącym kodem `dnet4` (`dn4_extract` → `dn2_extract` → `dn4_apply`, modele foldów `.xgb.`), 10 plików, interaktywnie (~3 min, 2 wątki).
Potem tabele zrekonstruowanego sygnału i prawdy z krotnością zdarzenia (`tx_mc.py`) oraz prawdziwy tryb per D z córek (`tx_excl.py`; `mb_sd` 10 plików, `mb_nd` 78 plików, 4 min).

| bin pT | C_ND | C_SD (wszystkie) | ρ₁ | C_SD (proton w oknie) |
|---|---|---|---|---|
| 3–5 | 3653/16366 = 0,223 ± 0,003 | 58/243 = 0,239 ± 0,027 | 1,07 ± 0,12 | 16/63 = 0,254 ± 0,055 |
| 5–8 | 1919/6394 = 0,300 ± 0,006 | 28/83 = 0,337 ± 0,052 | 1,12 ± 0,17 | 9/18 = 0,50 ± 0,12 |
| 8–30 | 848/2143 = 0,396 ± 0,011 | 12/30 = 0,400 ± 0,089 | 1,01 ± 0,23 | 0/3 |
| **3–30 (wagi g_b)** | | | **1,08 ± 0,09** | 1,17 ± 0,19 |

- **Ograniczenie statystyki:** `mb_sd` (1,72 M zdarzeń) ma 356 prawdziwych D⁰→K∓π± w obszarze fid. i 98 zrekonstruowanych przy WP70 (z protonem w oknie: 84 i 25). Oszacowanie wprost ma 9 % (16 % z protonem) — powyżej progu 30 zdarzeń z kryteriów, ale za mało na poprawkę per bin.
- **Oszacowanie 2 (krotność).** C_b(n) w `mb_nd` (nSelBL, przedziały o ≥ 100 zrekonstruowanych): zależność słaba poza najniższym przedziałem (3–5 GeV: ~30 % niżej dla n < 7, dalej bez trendu; X3).
  Kandydaci tagowani w danych (T2 ∧ MBTS ∧ WP70) mają medianę nSelBL 19 (okno) / 18 (pasma) wobec 28 dla sygnału `mb_nd` i 19 dla sygnału `mb_sd` — MC SD dobrze oddaje krotność zdarzeń tagowanych.
  ρ₂ = 0,962 ± 0,012 (okno, nominał), 0,966 (pasma), 0,967 (okno − pasma; część przedziałów ujemna → opisowo).
- **Domknięcie metody 2 na MC SD:** ρ₂ z krotności sygnału `mb_sd` = 0,960 (nSelBL) / 0,967 (nTracksPt500Eta25), wprost 1,08 ± 0,09 → niedomknięcie +12 ± 10 %. Przy tej statystyce nie da się odróżnić fluktuacji od efektu (np. inny rozkład η lub otoczenia D w SD).
  D⁺ ma niedomknięcie przeciwnego znaku (−17 ± 10 %: wprost 0,79 ± 0,09, z krotności 0,94); oba gatunki razem wprost 0,93 ± 0,07. Znaki przeciwne wskazują raczej na fluktuację niż wspólny efekt.
- **Przyjęte (reguła z kryteriów): ρ = 0,962 ± 0,120** (stat 0,012 ⊕ warianty 0,004 ⊕ niedomknięcie 0,119). ρ₁ wobec ρ₂(dane): 1,2σ (D⁰); D⁺ 2,1σ — formalnie „napięcie”, nominału nie zmieniam (D⁺ i tak opisowo).
- Definicja mianownika prawie nie zmienia ρ (0,956 → 0,962): udział prawdziwego trybu jest taki sam w SD (0,742 ± 0,020) i ND (0,734).
- Wymóg wierzchołka pierwotnego (filtr produkcji, poza C): 1,8 % prawdziwych D⁰ z protonem w `mb_sd` pochodzi ze zdarzeń o < 5 torach (0,4 % o ≤ 2) — w budżecie 1,8 %, bez poprawki.

## 3. Budżet (D⁰, mianownik po prawdziwym trybie; względny, %)

| pozycja | σ_tag^fid | R^fid | uwagi |
|---|---|---|---|
| statystyka N_tag | −34 / +36 | −34 / +36 | wariant: ±42 |
| kształt fitu tagowanego (TAGGED) | 23,2 | 23,2 | dla wariantu nie liczony |
| **model szerokich składowych (nominał → wariant)** | **−46** | **−46** | jednostronnie, osobno, poza kwadraturą |
| definicja mianownika C (prawdziwy tryb; werdykt CVER) | zastosowana (× 0,730 wobec C z XSEC); niepewność czynnika 0,3 % — pominięta | 0 | |
| obciążenie fitu (domknięcie inkluzywne M8) | 19,4 | skraca się (założenie) | 25 % w 3–5 GeV |
| ρ = C_D(tag)/C_D(ND) | 12,5 | 12,5 | statystyka `mb_sd` |
| SF klasyfikatora | 7,3 | 1,3 (reszta) | 3–5 GeV scenariusz |
| stat. N_incl / model tła N_incl | – | 8,4 / 5,5 | |
| A / ε_AFP / ustawienie-optyka | 3,1 / 2,1 / 3,4 | to samo | AFP razem 5,1 |
| nieprompt (z XSEC), stat. MC C_b, widmo pT (g_b) | 2,4 / 1,1 / < 0,1 | – / 0,2 / < 0,1 | |
| ε_wierzch, ε_npu, ε_ξcal, ε_trig, N_acc | 1,8 / 0,7 / 0,5 / 0,1 / 0,2 | 1,8 / 0,7 / 0,5 / 0,5 / 0,2 | |
| BR, L | 0,8 / 2,2 | skracają się | L zachowawczo 7 %: syst 34,8 % |
| **syst razem** | **34,1** (bez kształtu 25,0) | **28,8** (bez kształtu 17,0) | |

**Co się skraca w R^fid:** BR i L dokładnie; definicja mianownika C dokładnie (ten sam czynnik per bin w liczniku i mianowniku); C_b^ND, SF_b i obciążenie fitu inkluzywnego — do różnicy widm pT (g_b tagowane ≈ inkluzywne: reszta ≤ 1,3 %).
**Co się nie skraca:** ρ, A, ε_AFP, ustawienie/optyka, ε_ξcal, ε_npu, ε_wierzch(SD), N_acc, statystyka i model tła N_incl, kształt fitu tagowanego i model szerokich składowych.

## 4. Porównania

| | R^fid [%] | σ_tag^fid [µb] (poprawione C; po ukośniku przypis: C z XSEC albo inna normalizacja) | werdykt R (nominał / wariant) |
|---|---|---|---|
| **pomiar, nominał** | 0,99 −0,34 +0,35 ± 0,28; 68 %: [0,54; 1,44] | 13,6 −4,7 +4,8 ± 4,6 / 18,6 −6,4 +6,6 ± 6,5 | |
| **pomiar, wariant** | 0,54 ± 0,23 ± 0,09; 68 %: [0,29; 0,78] | 7,4 ± 3,1 ± 1,9 / 10,1 ± 4,3 ± 2,6 | |
| Pythia miękka A3 (SD_var; obwiednia zbioru A) | 0,35 (0,14–1,40) | 6,6 (2,6–26) przy σ_incl Pythii 1871 µb; 4,8 (1,9–19) przy σ_incl zmierzonym 1375 µb | zgodny (+1,4σ) / zgodny (+0,8σ) |
| Pythia twarda (osobny mechanizm) | 0,66 (0,58–1,41) | 12,3 (10,9–26) / 9,0 (8,0–19) | zgodny (+0,7σ) / zgodny (−0,5σ) |
| LMS 1412.3132 (IS + absorpcja, x_IP < 0,1), przeliczone | ≈ 0,5–1,2 | ≈ 1,1–2,5 | rząd wielkości, bez werdyktu |
| Siddikov–Schmidt 2008.12446 (dipol, cc̄, 13 TeV), przeliczone | ≈ 0,4–1,0 | – | rząd wielkości |

- **Pythia:** R per bin pT z `SD_var/out/results.json` (A3off: 0,375 / 0,324 / 0,268 %) zważone przekrojem inkluzywnym Pythii A3 per bin (1230 / 481 / 160 µb, XSEC T5) → 0,353 %; obwiednia × 0,40 … × 3,97 przeniesiona z R(2–30 GeV) — założenie, że czynniki nie zależą od pT.
  **`SD_var_2026-11-19/RAPORT_domkniecie.md` nie istniał w chwili składania (00:17 i ponownie ~00:40; job `charm_sdvar_close` biegł)** — użyto `RAPORT.md` / `results.json` z 19.11; po domknięciu obwiednię trzeba podmienić.
  Porównanie jest na poziomie cząstek, więc **bez czynnika detektorowego** i bez C — nie zależy od błędu mianownika.
  Przewidywania uzysku liczone jako Pythia × C_XSEC × BR były zaniżone o 0,734: z czynnikiem AFP 0,92 (zamiast 0,55) i poprawionym C przewidywane N(D⁰, WP70; 3–30 GeV) = 6,6 µb × 4,02 µb⁻¹ ≈ 27 (11–105 w obwiedni; ~20 przy C z XSEC, 12 przy starym czynniku 0,55); twarda ~49. Zmierzone 54,7 (wariant 29,9).
  Uwaga: przy mianowniku po prawdziwym trybie σ_incl(D⁰, 3–30) = 1375 µb, tj. dane/Pythia A3 ≈ 0,73 inkluzywnie (CVER: 276 i 97 µb w 5–8 i 8–30 GeV — zgodne z użytymi tu 276,2 i 97,1).
- **Wniosek:** R^fid leży w górnej połowie obwiedni miękkiej dyfrakcji Pythii i przy centralnej wartości scenariusza twardego; nie rozróżnia ich i żadnego nie wyklucza na 95 %. Wariant „okno − pasma” leży bliżej wartości centralnych.
- **Literatura (tylko z notatek bazy wiedzy; wg notatek liczby źródłowe „do sprawdzenia w PDF”):**
  LMS: σ_SD(D⁰ + anty-D⁰; |y| < 2,5, pT > 3,5 GeV, x_IP < 0,1, Pomeron + Reggeon 25 %, S_G = 0,05, 14 TeV) ≈ 3,56 µb; σ_SD/σ_incl ≈ 2 %. Siddikov–Schmidt: σ_SD/σ_incl(cc̄) = 1,59 % przy 13 TeV.
  Przeliczenie [założenia jawne]: pT ≥ 3 zamiast 3,5: × 1,36; |η| < 2,2 zamiast 2,5: × 0,89 (oba z prawdy `mb_sd`; y ≈ η); 14 → 13,6 TeV: × 0,98 [Z]; ułamek x_IP < 0,1 przypadający na 0,03–0,10: **0,26–0,6 [Z]**
  (dolna wartość = rozkład płaski w ln ξ od 10⁻³; górna = rosnący z ξ; rozkładu w ξ w notatkach nie ma) → σ ≈ 1,1–2,5 µb, R ≈ 0,5–1,2 % (LMS), 0,4–1,0 % (S–S, przy założeniu x_IP < 0,1 i R(D⁰) = R(cc̄)).
  Bezwzględnie LMS leży 3–12 razy poniżej pomiaru, ale rachunek jest kolinearny LO: 3,56 µb / 2 % = 178 µb inkluzywnie wobec ~1100–1550 µb z danych w tej samej akceptancji; miarodajny jest stosunek, a ten zgadza się co do rzędu.
  Obszar fiducjalny zawiera też proton wiodący ND (w Pythii 12–39 % tagów) i miękką dyfrakcję, nieobecne w rachunkach twardej SD.
- **Liczba pochodna (niezalecana):** ekstrapolacja do SD całkowitego × 3,9 (2,3–4,1) dałaby σ_SD(D⁰) ~ 50–70 µb (nominał) — tylko orientacyjnie, zależna od modelu.

## 5. Dla ≥ 5 GeV (składniki pod przyszły pomiar; bez σ_tag)

g_b(5–8, 8–30) = 0,79 / 0,21; C̄_tag(5–30) = 0,305 (prawdziwy tryb; 0,221 przy C z XSEC), SF zmierzony z tagu D* (12,5 / 10 %); Pythia R(5–30) = 0,31 % → ~6 (2–23) zrekonstruowanych D⁰ przy WP70 w 427929 przy normalizacji do zmierzonego σ_incl (niezależne od definicji C); przy bezwzględnym σ Pythii i poprawionym C ~10 (7 przy C z XSEC).
Potrzebny N_tag(5–30) modelem bez degeneracji szerokich składowych (propozycja TAGGED §7.1) — niewykonane.

## 6. Zastrzeżenia

1. **Post hoc** względem odczytu TAGGED; liczby N_tag znane przy pisaniu kryteriów. Istotność sygnału tagowanego (Z ≈ 3,0) bez zmian i bez systematyki modelu.
2. **Odstępstwo od kryteriów:** definicja mianownika C zmieniona po pierwszym przeliczeniu (KRYTERIA §10), z zewnętrznego powodu, potwierdzona potem werdyktem CVER; wynik przy starej definicji zachowany tylko jako przypis. R^fid i ρ zmieniły się o < 1 %; σ_tag^fid × 0,730.
3. **N_tag(D⁰) 30–55:** nominał i wariant obok siebie; wariant jest z konstrukcji zaniżony, nominał może być zawyżony. Dla wariantu systematyka kształtu nie była liczona, a błąd statystyczny nie zawiera niepewności ekstrapolacji tła.
4. **Bin 3–5 GeV = scenariusz** (SF = 1 ± 10 %, obciążenie fitu 25 % bez poprawki): dotyczy σ_tag^fid w µb; w R^fid prawie się skraca. W bazie wiedzy jest nowsza notatka o bezwzględnym SF przy WP70 (0,97 ± 0,03) — nie użyta (wejściem był XSEC); zmieniłaby σ o ≲ 3 %.
5. **ρ:** niedomknięcie 12 % to statystyka `mb_sd`, nie zmierzony efekt; zejście poniżej 5 % wymaga ~10 razy więcej MC SD z D (próbki wymuszonych rozpadów w SD w h5 v2.0 nie ma). SF klasyfikatora w zdarzeniach z protonem przyjęty jak inkluzywnie.
6. Krotność w danych wzięta z kandydatów w oknie (sygnał + tło), nie z czystego sygnału; warianty zmieniają ρ o ≤ 0,5 %.
7. **Obszar fiducjalny nie jest „σ_SD”:** zawiera proton wiodący ND i DD; A = 0,97 ± 0,03 pokrywa udział ND 17–39 %. Wymóg „dokładnie jeden proton” i cięcie ξ_cal są wydajnościami (ε_AFP; ε_ξcal = 1,000 +0/−0,005 — przeniesienie z kandydatów tła na sygnał jest założeniem noty).
8. ε_npu = 0,985 ± 0,007 z próbki inkluzywnej; przy małej krotności strata „drugi wierzchołek” może być mniejsza (do −0,9 % na σ, w ramach błędu).
9. Obwiednia Pythii pochodzi sprzed domknięcia SD_var i jest przeniesiona z 2–30 GeV; przeliczenia literatury zawierają założenie o rozkładzie ξ o rozpiętości × 2,3.
10. D⁺ tylko opisowo (K3 WP50 i K4 niespełnione; SF scenariusz we wszystkich binach; ρ(D⁺) wprost 0,79 ± 0,09, 2,1σ od wartości z krotności).
11. Okno sygnału D⁺ liczone warunkiem |m − m₀| < 0,045 daje 99 kandydatów wobec 100 w histogramie 5 MeV TAGGED (krawędź binu) — dotyczy tylko rozkładu krotności.

## 7. Co uruchomiono na lxplus, pozostałości

- `bin/lx -check`: połączenie było (lxplus9103). `condor_q` sprawdzone: `charm_dnet6`, `charm_sdvar_close`, później `charm_xsatlas` (inne etapy) — nietknięte; obliczenia bez condora.
- Interaktywnie w `/tmp/klyczak_tagxsec`: `tx_app.sh mb_sd` (~3 min), `tx_mc.py mb_sd` / `mb_nd` (50 s / 48 s), `tx_excl.py mb_sd` / `mb_nd` (25 s / 4 min, 2 procesy). Kod `dnet4` / `dnet2` czytany z `~/private/charm_q/Q_lxplus` bez zmian; nic nie zapisano do `charm_q/dnet4`.
- Wyniki na EOS: `/eos/user/k/klyczak/charm_q/tagxsec/` (`app/app_mb_sd_0_*.npz` × 10, `tx_mb_{sd,nd}.npz`, `ex_mb_{sd,nd}.npz`, logi). Lokalnie: `runs/lxplus/`.
- Wpadka własna: pętla oczekiwania `until … pgrep -f tx_excl.py` na lxplus dopasowała samą siebie i wisiała ~7 min (zatrzymana, proces usunięty `pkill` po własnym wzorcu; innych procesów nie dotyczyło).
- **Pozostałości (nie usuwałem):** lxplus9103 `/tmp/klyczak_tagxsec/` (~45 MB: `app/`, `wd/` 37 MB plików pośrednich, `gen_mb_nd.npz`, tabele, skrypty, logi); lokalnie `skrypty/lxplus/gen_mb_nd.npz` (1,4 MB), `skrypty/tx_ana.py.bak_przed_trybem`,
  `out/tagxsec.json` i `tabele/wyniki.md` (pierwsze przeliczenie przed zmianą definicji; treścią równe wersjom `_xsec`).

## 8. Pliki (`~/Documents/charm/charm_work/TAGXSEC_2026-11-20/`)

- `KRYTERIA.md` (§1–9 przed liczeniem, §10 dopisek po), `RAPORT.md`, `akapit_do_noty.md`.
- `skrypty/tx_ana.py [excl|xsec]` (lokalnie: ρ, złożenie, budżet, tabele), `skrypty/lxplus/{tx_app.sh, tx_mc.py, tx_excl.py, push.sh}`.
- `tabele/wyniki_excl.md`, `tabele/wyniki_xsec.md` (X1 C per bin, X2 ρ, X3 C(n), X4 wynik i budżet dla D⁰ i D⁺, nominał i wariant), `out/tagxsec_{excl,xsec}.json`, `runs/lxplus/` (tabele MC i logi).
