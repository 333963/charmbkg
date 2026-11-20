# KONSTR — tło „przez konstrukcję” w trzech słabych miejscach pomiaru (projekt B), 20.11.2026

Kryteria: `KRYTERIA.md` (00:57 CET, przed ekstrakcją i fitami; `out/KRYTERIA_czas.txt`). **Całość post hoc** (FITBIAS / M8 / BKGCAT, odczyt TAGGED 19.11, D23 były znane); nic nie zastępuje nominału.
Bez condora, bez zmian w `ORCH` i nocie.
Ekstrakcja: `skrypty/kx_extract.py`, interaktywnie na lxplus9103 `/tmp/klyczak_konstr` (2 procesy; 14 s + 144 s), MC lokalnie. Kontrola: liczby wierszy 427929 WP70 = `tg_incl` co do zliczenia (6 komórek).
Fity okna sygnału: A 16 (+ 8 diagnostycznych), B 13 (+ 5 w dodatku), C 56. Bez poprawki na wielokrotność.

## Werdykty

| | co zbudowano | test rzeźbienia / zamknięcia | przed → po | werdykt |
|---|---|---|---|---|
| **A** odbicie K↔π z danych | pomiar liczby wierszy odbicia z piku w `m_swap` bliźniaka (poza skrzyżowaniem pasm), kształt fd, normalizacja z danych | estymator „okno − pasma” rzeźbi (3–5 GeV: samo tło udaje 2725 ± 424 odbić); wersja z fitem per przedział: 8–30 GeV 0,89 ± 0,14 prawdy (kryterium 10 % nierozstrzygalne), 5–8 OK na MC, ale w danych χ² kształtu 30–40/8; 3–5 nie do użycia | 428770 8–30 WP70: M8 6829 ± 236, p = 1,1·10⁻⁷ → M8-R 6812 ± 235, p = 1,1·10⁻⁷; M0 7716 ± 332 (p 0,43) → M0-R 7716 ± 332 (0,43) | **brak zysku** — odbicie nie jest przyczyną odrzucenia M8 |
| **B** tło fitu tagowanego | kształt tła z próbki nietagowanej tej samej selekcji (205 967 wierszy D⁰), wolne tylko N_S, N_B | pasma tagowane = nietagowane bez wag (p 0,45; D⁺ 0,40); zamknięcie 3 × 200 pseudo-próbek: ⟨pull⟩ −0,07…−0,01, pokrycie 0,69–0,72 | N_tag 54,7 −18,7 +19,4 (Z 3,05) → 54,1 −18,3 +19,0 (Z 3,05); rozrzut modelu 30–55 → 44–54 | **zysk umiarkowany** (wariant / sprawdzian, nie nominał): centralna wartość i błąd bez zmian, rozrzut modelu ~2× mniejszy |
| **C** pik D⁰ 2–3 GeV | tło = rozkład masy z odwróconego cięcia × exp(c₁x) | pasma: 427929 drzewo NIE (p 0,004–0,011), pozostałe 3 OK; MC 3–5 GeV: sygnał pozorny ujemny, −2…−22 % tła w oknie (10 z 16 testów niezaliczonych) | Z_LR funkcyjnie 3,2–12,0 (rozstęp 5,8–7,0) → konstrukcyjnie 6,8–15,3 (rozstęp 4,9–6,8) | **brak zysku** — zależność przeszła z rzędu wielomianu na rząd transferu |

## A. Szablon odbicia K↔π z danych (`tabele/A_mc.md`, `A_mc2.md`, `A_dane_428770.md`, `A_dane_427929_pomiar.md`)

- **Konstrukcja.** Histogram 2D (m × m_swap) wierszy przechodzących WP70 własnym score. Odbicie = pik w m_swap przy m₀ w przedziałach m (25 MeV); 5 przedziałów skrzyżowania |m − m₀| < 50 MeV pominiętych; do 9 punktów dopasowana normalizacja A_R szablonu odbicia fd.
- **Test na mb_nd z prawdą.** E1 (okno − pasma liniowo, wg kryteriów): 3–5 GeV WP70 A/prawda = 1,92 ± 0,16, samo tło +6,4σ; 8–30: 0,86 ± 0,13, ale χ² wobec prawdy 22,6/9 (p = 0,007) → **niezaliczone**.
  **Odstępstwo 1** (po teście MC, przed danymi): estymator E2 — fit per przedział [kształt sygnału w m_swap + wielomian 2] w ±120 MeV; pierwsze wypróbowane ustawienie, bez wyboru spośród czterech sprawdzonych.
  E2, WP70: 3–5 GeV 1,42 ± 0,16, samo tło +3,2σ (NIE); 5–8: 1,09 ± 0,11 (χ² 16,5/9), przy czystości danych 1,12 ± 0,03 (NIE o 2 pkt); **8–30: 0,89 ± 0,14 (χ² 7,9/9), samo tło −3 ± 63, przy czystości danych 0,88 ± 0,05** — w 1σ, lecz próg 10 % nie do rozstrzygnięcia tą statystyką MC; przy WP50 / WP90 ten sam znak (0,79 / 0,68).
  Skala szerokości s_R niemierzalna (0,84–1,60 między testami) → s_R = 1: „szablon z danych” to kształt MC z normalizacją z danych.
- **Wg kryteriów szablon nie jest poświadczony w żadnej komórce.** **Odstępstwo 2:** mimo zapisu „fitów okna nie wykonuję” policzono M8-R / M0-R jako diagnostykę (odpowiedź na pytanie z planu); liczby orientacyjne.
- **Pomiar na danych (bez fitu okna):** A_R/N_S względem f_refl z MC — 428770 8–30: WP70 1,14 ± 0,09, WP50 1,02 ± 0,05; 5–8: 1,11 ± 0,21 / 0,98 ± 0,11 (χ² kształtu 40/8 i 17/8 — w danych 5–8 estymator już rzeźbi); 427929 WP70 (5–8 / 8–30): 1,15 ± 0,27 (χ² 30/8) / 1,16 ± 0,18 (7,9/8).
  Udział odbicia w danych zgodny z MC w ≤ 1,6σ; dokładność 5–9 % nie lepsza niż dotychczasowa z tagu D* (0,957 ± 0,041).
- **Fity 428770 8–30 GeV WP70** (HESSE; bootstrap 60): M0 7716 ± 332 (518), p = 0,43; M8 6829 ± 236 (301), p = 1,1·10⁻⁷; M8-R (f_refl 0,81 zamiast 0,71) 6812 ± 235 (327), p = 1,1·10⁻⁷; M0-R 7716 ± 332 (524), p = 0,43. Zmiana N ≤ 0,3 %, p bez zmian.
  WP50: M8 5377 (p 3,5·10⁻⁷) → M8-R 4684 (p 3·10⁻⁹; przeskok minimum, sc 1,10 → 0,97) — niestabilność samego M8. 5–8 GeV: zmiany ≤ 0,3 %.
- **Co odrzuca M8 (diagnostyka poza kryteriami, `out/ka_diag_428770.log`):** M8 bez składowej DMM + CM: 7504 ± 235, **p = 0,27**; M7 (kształt wolny, 1,70–2,05, bez DMM + CM) 7534 ± 237, p = 0,29; M8 z rzędem tła + 1 / + 2: p ~ 10⁻⁷; M8 z kształtem wolnym: p ~ 10⁻⁷.
  **Odrzucenie pochodzi w całości od szablonu DMM + CM** (κ i kształt z mb_nd dla próbki dżetowej), nie od odbicia, ograniczenia kształtu ani zakresu. Uzysk przesuwa się przy tym o −9 % (6829 wobec 7504).
- Wniosek praktyczny: dla 428770 8–30 GeV kandydat na nominał to „M8 bez DMM + CM” (7504 ± 235, p 0,27; błąd 0,71 × M0) — do decyzji w etapie fitów; bliźniaki nie dają szablonu DMM + CM.

## B. Tło fitu tagowanego z konstrukcji — POST HOC po odczycie TAGGED, lokalnie (`tabele/B_prep.md`, `B_wyniki.md`, `B_extra.md`)

- **Konstrukcja.** D⁰, 427929 MBTS, T2, 3–30 GeV, WP70. b̂(m) = składowa tła z fitu M8 rozkładu nietagowanego (MB ∧ GRL ∧ npu = 0) per bin pT, zsumowana; Ŝ jak w TAGGED; fit tagowany: wolne N_S, N_B.
- **Czy przeważać (tylko pasma, 516 wierszy D⁰ / 267 D⁺):** próbka tagowana ma dużo mniejszą krotność (mediana nSelBL 18 wobec 38), ale kształt masy pasm jest zgodny bez wag: D⁰ χ² 7,9/8 (p 0,45), nachylenie stosunku −0,28 ± 0,39 /GeV; z wagami nSelBL 6,9/8, nSelBL × pT 6,7/8. D⁺: 15,8/15 (p 0,40), nachylenie −0,57 ± 0,40.
  Nominał konstrukcji = bez wag (reguła „najprostszy z p ≥ 0,05”). Moc testu mała (błąd nachylenia 0,4 /GeV wobec nachyleń −0,08…−0,29 /GeV zmierzonych w TEMPLATE 17.11 dla składowej przypadkowej).
  Uwaga: przedziały nSelBL przestawiono na kwantyle pasm tagowanych po obejrzeniu samych pasm (pierwsze były na kwantylach próbki nietagowanej) — przed zamknięciem i fitem okna.
- **Zamknięcie, 200 pseudo-próbek na wiersz** (losowanie wierszy nietagowanych ∝ waga „prawdy”; pełna procedura na każdej):

| gatunek | wagi prawdy / procedury | N oczek. | konstrukcja: obciążenie; rms N̂; ⟨pull⟩; pokrycie | TAGGED (tło wolne): obciążenie; rms N̂; ⟨pull⟩; pokrycie |
|---|---|---|---|---|
| D⁰ | brak / brak | 23,1 | −0,9 ± 1,3; 18,2; −0,07; 0,69 | +0,5 ± 1,3; 18,9; +0,01; 0,67 |
| D⁰ | nSelBL × pT / brak | 43,0 | −0,8 ± 1,2; 17,3; −0,05; 0,72 | −1,1 ± 1,3; 17,8; −0,07; 0,70 |
| D⁰ | nSelBL × pT / nSelBL × pT | 43,0 | −0,1 ± 1,3; 18,3; −0,01; 0,69 | +0,1 ± 1,4; 19,5; −0,01; 0,70 |
| D⁺ | brak / brak | 19,2 | +0,4 ± 0,7; 10,2; +0,02; 0,66 | −1,2 ± 0,7; 10,3; −0,13; 0,71 |
| D⁺ | nSelBL × pT / brak | 36,6 | −0,4 ± 0,7; 10,3; −0,06; 0,70 | −2,1 ± 0,8; 10,8; −0,20; 0,69 |
| D⁺ | nSelBL × pT / nSelBL × pT | 36,6 | +0,9 ± 0,7; 10,5; +0,06; 0,72 | −1,2 ± 0,8; 10,9; −0,13; 0,73 |

  Wszystkie zaliczone. Zysk statystyczny mały: rms N̂ × 0,94–0,97 (D⁰), × 0,95–0,99 (D⁺); dla D⁺ znika obciążenie tła wolnego (−1,2…−2,1 → zgodne z 0; K4 z TAGGED).
  Efekt uboczny: próbka nietagowana przeważona do krotności tagowanej ma ~1,9 × większy udział sygnału (43 wobec 23 D⁰ na ~770 kandydatów) — przy tym samym ułamku D na kandydata co w zdarzeniach nietagowanych o tej krotności oczekiwane byłoby ~43 D⁰; zmierzone 54 leży 0,6σ wyżej.
- **Fity okna (D⁰):**

| wariant | N_tag | −/+ (68 %) | Z | tło w oknie |
|---|---|---|---|---|
| TAGGED (tło wolne, rząd z pasm) — odtworzenie | 54,7 | −18,7 +19,4 | 3,05 | 188 |
| **konstrukcja nominalna (M8, bez wag)** | **54,1** | −18,3 +19,0 | 3,05 | 189 |
| wagi nSelBL / nSelBL × pT | 54,6 / 54,4 | −18,2 +18,8 | 3,11 / 3,10 | 187 |
| zakres 1,75–2,00 | 46,0 | −19,7 +20,4 | 2,39 | 204 |
| (ii) DMM + CM w tle: M7 / M0 (1,75–2,00) | 46,6 / 44,3 | −15,6 +16,3 / −18,5 +19,0 | 3,10 / 2,45 | 217 / 222 |
| dodatek: wszystkie szerokie składowe w tle (M7 bez odbicia; N = wąski pik) | 52,8 | −18,1 +19,0 | 3,02 | – |
| okno − pasma (TAGGED §7) | 29,9 | – | – | 228 |

  Ułamek z mianownikiem tego samego modelu (`B_extra.md`): 0,87 −0,30 +0,31 % (M8), 0,78 (M7), 0,65 (M0), 0,81 (M7 bez odbicia), 0,69 (M0 bez odbicia).
  D⁺ (kontrolnie): 29,1 −10,8 +11,5 (Z 2,88) → 32,2 −10,9 +11,5 (Z 3,21); zakres 1,75–2,10: 26,9; okno − pasma 31,7.
- **Odczyt.** (1) Tło wolne i tło z konstrukcji dają to samo N_tag (54,7 / 54,1) i tę samą istotność — swoboda kształtu tła nie była źródłem rozrzutu. (2) Dolny kraniec „~30” nie jest odtwarzany przez żaden wariant z tłem o kształcie nietagowanym: nawet przy wszystkich szerokich składowych w tle wąski pik daje 53;
  „30” wynikało z ekstrapolacji samych pasm tagowanych (tło w oknie 228 wobec 187–222). (3) Rozrzut po wariantach modelu (ii): 7,5–9,8 zdarzenia (kryterium < 12 spełnione); pełny zakres wariantów 44–55; największa pozycja to zakres fitu (−8,1), statystyczna w naturze.
  (4) Systematyka kształtu liczona jak w TAGGED, ale z konstrukcją: zakres 8,1 ⊕ model składowych 7,5 ⊕ wagi 0,5 = 11,0 (było 12,6 i dodatkowo niepokryte „30”). **Werdykt wg kryteriów: zysk** (zamknięcie zaliczone, rozrzut (ii) < 12, σ_stat nie większe).
- Zastrzeżenia: post hoc po odczycie; założenie „proporcje szerokich składowych w zdarzeniach z tagiem jak inkluzywnie” zostaje (konstrukcja bada tylko jego wpływ: 54 ↔ 47); test kształtu pasm ma małą moc; f_true, przewidywania i TAGXSEC nieprzeliczane.

## C. Pik D⁰ 2–3 GeV: tło z odwróconego cięcia (`tabele/C_prep.md`, `C_wyniki.md`)

- **Konstrukcja.** B(m) = H_ctrl(m) · exp(c₁x [+ c₂x²]); H_ctrl: K1 = przedział ε 0,5–0,3 (14–60 × więcej kandydatów niż w punkcie pracy), K2 = odrzucone ε > 0,5 (1600–9000 ×). Trzeciej próbki z kryteriów („reszta poniżej WP10”) nie zrobiono — pokrywa się z K1 + K2.
  Zanieczyszczenie K1 sygnałem poprawione iteracyjnie (N/ε · Δε). Bliźniaki K↔π nie dają próbki tła; złego znaku nie ma.
- **(a) pasma danych:** 428770 drzewo i sieć, 427929 sieć: p = 0,03–1,0, rząd 2 niepotrzebny (2ΔlnL ≤ 1,4) → OK. **427929 drzewo WP10: p = 0,004–0,011 i 2ΔlnL = 5,4–5,8 w zakresie 1,75–2,00 → NIE** (ten histogram ma p ≈ 0,05 także z każdym tłem funkcyjnym).
- **(b) mb_nd 3–5 GeV, samo tło:** sygnał pozorny zawsze ujemny: drzewo WP10 −4 ± 13 … −27 ± 12 (−2…−10 % tła w oknie), drzewo WP20 do −84 ± 24 (−3,5σ), sieć −19 ± 10 … −67 ± 18 (−9…−22 % tła, do −3,8σ); 10 z 16 testów niezaliczonych (|Z| ≥ 2). Tło funkcyjne exp2 na tych samych histogramach: −22…+18, w 1,5σ.
  Transfer exp(c₁x) nie przenosi kształtu: cięcie wzbogaca składowe powabowe o innym kształcie (BKGCAT: DMM ma ujemny sygnał pozorny).
- **Okno (drzewo WP10 / sieć; `C_wyniki.md`):** rząd transferu 1 → N zmienia się o −14…+3 % (np. 427929 drzewo 449 → 431 i 460 → 404; 428770 drzewo 502 → 437–440 i 459 → 423–433), σ × 0,67–0,79, Z_LR 11–15; rząd transferu 2 → wyniki jak exp2 (N co do 2 %, σ i Z co do 0,1).
  Rozstęp Z_LR: funkcyjne 5,0–12,0 / 5,0–11,7 / 3,6–9,4 / 3,2–9,9 → konstrukcyjne 8,5–15,3 / 9,1–15,0 / 7,7–12,6 / 6,8–12,5 (427929 drzewo / sieć, 428770 drzewo / sieć). Rozstęp 4,9–6,8 wobec 5,8–7,0 — kryterium „≤ połowa” niespełnione.
- **Werdykt: brak zysku.** Istotność nadal zależy od jednego wyboru (rząd transferu zamiast rzędu wielomianu), a wariant sztywny (rząd 1) ma na MC obciążenie ujemne i w 427929 (drzewo) nie opisuje pasm. Uboczne: najmniejsze Z po wariantach konstrukcyjnych to 6,8–9,1 (funkcyjnie 3,2–5,0, wyłącznie z exp4 na wąskim zakresie) — pik istnieje w każdym wariancie; uzysk zależy od modelu tła na poziomie ±7–10 %.

## Odstępstwa od kryteriów

1. A: estymator E1 zastąpiony przez E2 po niezaliczeniu testu MC (przed danymi). 2. A: fity -R wykonane mimo braku poświadczenia (diagnostyka). 3. A: bootstrap 60 zamiast 100 i tylko w komórce docelowej (obciążenie maszyny przez równoległe przebiegi); dodatkowe fity diagnostyczne M7, M8 bez DMM + CM, rząd + 1 / + 2.
4. B: przedziały nSelBL zmienione po obejrzeniu pasm (przed zamknięciem); zamknięcie w trzech układach wag zamiast jednego; dodatek `B_extra.md` (f z mianownikiem modelu, wariant bez odbicia) po fitach okna. 5. C: dwie próbki kontrolne zamiast trzech; kryterium (c) zrealizowane jako poprawka iteracyjna.

## Pozostałości (nic nie usunięto)

- lxplus9103 `/tmp/klyczak_konstr/` (~19 MB: `kx_extract.py`, `run.sh`, `run.log`, `kx_data_427929.npz`, `kx_data_428770.npz`, `DONE`). Nic na EOS ani AFS; żadnych jobów condora.
- Lokalnie: `KONSTR_2026-11-20/in/` (19 MB), `skrypty/__pycache__` (jeśli powstał), `out/kb_prep_sb.json` (wynik pośredni).
- Proces `ka_data.py` z pierwszego, zbyt wolnego przebiegu zakończony sygnałem.

## Pliki

`KRYTERIA.md`, `RAPORT.md`, `NOTA_tekst_EN.md`; `skrypty/`: `kx_extract.py`, `ka_lib.py`, `ka_mc.py`, `ka_mc2.py`, `ka_data.py` (A), `kb_tag.py`, `kb_extra.py` (B), `kc_d23.py` (C);
`tabele/`: `A_mc.md`, `A_mc2.md`, `A_dane_428770.md`, `A_dane_427929_pomiar.md`, `B_prep.md`, `B_wyniki.md`, `B_extra.md`, `C_prep.md`, `C_wyniki.md`; `out/`: json + logi; `in/`: histogramy 2D i wiersze.
