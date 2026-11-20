# SIG_constr — „sygnał z konstrukcji”: tag D*⁺→D⁰π⁺ w danych zamiast sygnału MC (projekt B)

Stan: 17.11.2026 08:15. Kryteria: `KRYTERIA.md` (zapisane 17.11 07:17, przed liczbami). Tabele pełne: `out/tables.md` (T1–T8). Liczby: `out/*.json`.
Część B (embedding): pełny bieg `sigemb` pobrany i przeliczony 17.11 08:50 (§6.3–6.4; aktualizacja raportu 09:00).

## 0. Werdykty

| pytanie | werdykt | liczby decydujące |
|---|---|---|
| K0 statystyka | trening od zera **niewykonalny** dla H (n_eff 102–191 ≪ N*(5 %) = 3000); dla zadania z izolacją na granicy (n_eff 440 / 894 vs N*(5 %) = 1000); przeważanie wykonalne tylko pT-całkowo w H06 | T1, T2 |
| (a) zysk istotności, klasyfikator H (18 cech) | **NIE** — każdy wariant z danych przegrywa albo remisuje | 9:MB: przeważenie 0,962 ± 0,025, od zera 0,899 ± 0,030, dostrojenie 0,980 ± 0,017; 8:J12only: 0,972 ± 0,013 / 0,930 ± 0,016 / 0,972 ± 0,013 |
| (a) zysk istotności, 18 cech + izolacja | **TAK w 8:J12only, NIE w 9:MB** | J12only WP80: przeważenie 1,055 ± 0,015, od zera 1,074 ± 0,020, dostrojenie 1,093 ± 0,018; 9:MB: 0,995 ± 0,016 / 0,957 ± 0,020 / 1,002 ± 0,013 |
| (b) mniejsza systematyka ε_S | **NIE** dla izolacji (wagi bez koła nie domykają binów pT w J12only; przenoszenie między klasami i na inkluzywne D⁰ zawodzi); dla H formalnie domyka, ale przeważenie nie zmienia ε (≤ 2 %) — brak zysku | J12only ptcone40 < 2: bez koła 0,896 ± 0,088 (3–30), 0,705 ± 0,115 (5–8, pull +2,6), 0,623 ± 0,125 (8–30, pull +3,0); z kołem 0,992 ± 0,131 |
| (c) kalibracja izolacji dla sieci typu CharmNet w próbce dżetowej | **WARUNKOWO**: c1 TAK, c2 TAK, c3 NIE dla cięcia bezwzględnego przy pT 8–30 GeV (−10 %), TAK dla 3–8 GeV i dla izolacji względnej | ε_S modelu z ISO z tagu: błąd wzgl. 3,0 % (J12only), 3,7 % (9:MB); przenoszenie ptcone40 < 2: +0,7 / −1,2 / −10,0 % |
| (B) embedding | pełny bieg: kierunek losowy/prostopadły **nie odtwarza** ε tagu (Δ = 0,422 ± 0,016 w J12only); wariant bez koła WO: Δ = 0,210 ± 0,016; WB domyka (0,022 ± 0,021), ale jest w ~80 % kołowy → aktywność skorelowana z D⁰; embedding mierzy składową nieskorelowaną, nie jest kalibracją | §6 |

Rekomendacja (§7): nie trenować ani nie przeważać H na sygnale z tagu; izolację w próbce dżetowej brać jako ε z danych per strumień (pomiar wprost z tagu D*), w zmiennej względnej ptcone40/pT; sieć z izolacją na wejściu opłaca się tylko w strumieniu dżetowym i tylko uczona/dostrojona na danych tego strumienia.

## 1. Metoda (skrót; definicje w `KRYTERIA.md`)

- Próbka: pary D⁰←D* z `bt_*.npz`, |m(Kπ) − m₀| < 40 MeV, pT 3–30 GeV. Klasy: 9:MB (427929, MBTS_1 ∧ GRL), 8:J12only, 8:j20&!J12, 8:MBTSbased, 8:other, 8:all (428770). Obszary: PRE (bez cięcia score), H06 (score_H > 0,6).
- Sygnał z danych: wagi odejmowania pasm Δm (SBS, `iso_common.sbs_w`) — metoda główna. sWeights (własny fit Δm) podane tylko przy czystości okna ≥ 0,25; przy niższej fit jest niestabilny (N_s → 0), a przy czystości 0,3–0,5 zaniża N_sig o 20–50 % względem SBS. Zgodność n_eff obu metod dopiero przy czystości ≥ 0,8 (237 vs 240).
- Modele (xgboost, 5 foldów po zdarzeniach, wszystko OOF): MC (sygnał MC vs pasma Δm danych), RW (sygnał MC przeważony w(x) = p_dane/p_MC, klasyfikator stosunku z ujemnymi wagami SBS), DATA (od zera: okno Δm vs pasma, CWoLa), STACK (dostrojenie: DATA z marginesem modelu MC), H (sam score). Zadanie H: PRE, pT 5–30, 18 cech. Zadanie I: H06, pT 3–30, 18 cech + ptcone40/100 (+ względne) + score.
- Z_rel = ε_S(wariant)/ε_S(model MC) przy tym samym ε_B (z pasm Δm, WP modelu MC), ε_S z danych przez SBS; błąd: bootstrap sparowany po zdarzeniach (100 replik).
- Wyciek masy/Δm: R² regresji m(Kπ) i Δm z wyjścia modeli ≤ 0,002 (próg 0,02) — brak.
- Domknięcie: ε_S(MC przeważone, OOF) vs ε_S(dane) z fitu Δm pass/fail (SBS gdy N_sig < 150); błąd = stat danych ⊕ stat MC ⊕ rozrzut wag (jackknife po 5 modelach foldów).

## 2. n_eff (K0) — `fig1_neff.png`, T1

| klasa | obszar | 3–5 GeV | 5–8 GeV | 8–30 GeV | 3–30 GeV |
|---|---|---|---|---|---|
| 9:MB | PRE | N 638 ± 106, n_eff 36 | 356 ± 33, 117 | 141 ± 15, 94 | 1135 ± 112, 102 |
| 9:MB | H06 | 538 ± 40, 184 | 302 ± 19, 240 | 116 ± 11, 107 | 956 ± 46, 440 |
| 8:J12only | PRE | 870 ± 180, 23 | 604 ± 83, 54 | 1277 ± 56, 516 | 2752 ± 205, 179 |
| 8:J12only | H06 | 629 ± 56, 126 | 497 ± 31, 249 | 1067 ± 35, 911 | 2193 ± 73, 894 |
| 8:j20&!J12 | H06 | 129 ± 22, 33 | 74 ± 11, 45 | 93 ± 10, 87 | 296 ± 27, 122 |
| 8:MBTSbased | H06 | 152 ± 25, 37 | 94 ± 13, 56 | 124 ± 12, 116 | 370 ± 30, 149 |
| 8:all | H06 | 888 ± 64, 191 | 664 ± 35, 352 | 1246 ± 38, 1074 | 2797 ± 83, 1146 |

Sygnał MC dla porównania: PRE 48 041, H06 39 865 (3–5: 21 243, 5–8: 13 009, 8–30: 5613).
Bez cięcia score czystość okna w 3–5 GeV to 3–6 % i n_eff spada do 23–36 — w tym binie sygnał z tagu jest „niemierzalny” wg K0 (< 100). Próg 300 dla przeważania spełniają tylko: 9:MB H06 3–30 (440), 8:J12only H06 3–30 i 8–30 (894, 911), 8:J12only PRE 8–30 (516), 8:all.

## 3. Krzywa uczenia — `fig2_krzywa_uczenia.png`, T2

Z_rel (ε_S przy ε_B z WP80 modelu pełnego) w funkcji liczby zdarzeń sygnału MC w treningu:

| N | 100 | 300 | 1000 | 3000 | 10 000 | N*(2 / 5 / 10 %) |
|---|---|---|---|---|---|---|
| zadanie H (18 cech, PRE) | 0,81 | 0,88 | 0,945 | 0,970 | 0,992 | 10 000 / 3000 / 1000 |
| zadanie I (+ izolacja, H06) | 0,82 | 0,94 | 0,976 | 0,994 | – | 3000 / 1000 / 300 |

Dostępne n_eff (102–191 dla H; 440–894 dla I) leży 1–2 rzędy poniżej tego, co potrzebne, by trening od zera nie tracił. Obserwowana strata treningu od zera na danych (−7…−10 % dla H, −4 % dla I w 9:MB) zgadza się z odczytem krzywej przy tych n_eff — czyli strata jest statystyczna, a nie wynika z niedopasowania MC.
Dopisek wg K(a): wariant od zera wyrównałby model MC dla H przy n_eff ≈ 10⁴ (≈ 50–100× obecnej próbki D*-tag bez cięcia score).

## 4. Warianty treningu na danych — `fig3_warianty.png`, T3

Z_rel względem modelu MC, ε_B z WP80 modelu MC:

| zadanie | klasa | RW (MC przeważone) | DATA (od zera) | STACK (dostrojenie) | sam score H |
|---|---|---|---|---|---|
| H | 9:MB | 0,962 ± 0,025 | 0,899 ± 0,030 | 0,980 ± 0,017 | 0,982 ± 0,030 |
| H | 8:J12only | 0,972 ± 0,013 | 0,930 ± 0,016 | 0,972 ± 0,013 | 0,939 ± 0,015 |
| I | 9:MB | 0,995 ± 0,016 | 0,957 ± 0,020 | 1,002 ± 0,013 | 0,778 ± 0,026 |
| I | 8:J12only | 1,055 ± 0,015 | 1,074 ± 0,020 | 1,093 ± 0,018 | 0,693 ± 0,023 |

8:J12only, zadanie I, inne punkty pracy: WP95 — 1,013 / 0,997 / 1,028; WP50 — 1,109 ± 0,025 / 1,100 ± 0,035 / 1,125 ± 0,030.

- Klasyfikator H: MC wygrywa wszędzie. Zgodne z B2/B3 (C2ST 0,512; NO-GO retreningu) — tu potwierdzone dodatkowo per klasa triggera i z dostrojeniem.
- Izolacja jako cecha daje dużo (sam score traci 22 % w 9:MB i 31 % w J12only względem modelu z izolacją), ale w J12only model MC ma złą izolację: jego ε_S na danych to 0,671 ± 0,020 przy 0,883 w MC (WP80). Warianty z danych odzyskują 5–12 % Z. W 9:MB izolacja w MC jest dobra (ε_S dane 0,750 ± 0,028 vs MC 0,791) i zysku nie ma.
- Zastrzeżenie: zysk w J12only jest względem modelu MC z izolacją, w tej samej klasie; warianty z danych są uczone na tej klasie (OOF po zdarzeniach). Błąd Z_rel nie zawiera rozrzutu samego treningu.

## 5. Domknięcie przeważenia sygnału (K(b)) — `fig4_domkniecie.png`, T4–T6

### 5.1 Izolacja (H06), ε_dane/ε_MC,przeważone dla ptcone40 < 2 GeV

| klasa, pT | ε_dane | wagi pT | EV (bez koła) | EV+F18 (bez koła) | ISO (z kołem) |
|---|---|---|---|---|---|
| 9:MB 3–30 | 0,649 ± 0,029 | 0,96 ± 0,04 | 0,96 ± 0,04 | 0,96 ± 0,04 | 0,96 ± 0,08 |
| 8:J12only 3–30 | 0,304 ± 0,018 | 0,619 ± 0,039 (pull +9,9) | 0,896 ± 0,088 (+1,2) | 0,867 ± 0,118 (+1,1) | 0,992 ± 0,131 (+0,1) |
| 8:J12only 3–5 | 0,568 ± 0,060 | 0,753 ± 0,080 (+3,1) | 0,844 ± 0,116 (+1,4) | – | 0,929 ± 0,205 (+0,3) |
| 8:J12only 5–8 | 0,355 ± 0,036 | 0,576 ± 0,058 (+7,3) | 0,705 ± 0,115 (+2,6) | – | 0,890 ± 0,173 (+0,6) |
| 8:J12only 8–30 | 0,119 ± 0,010 | 0,432 ± 0,047 (+12,1) | 0,623 ± 0,125 (+3,0) | – | 0,994 ± 0,103 (+0,1) |
| 8:j20&!J12 3–30 | 0,409 ± 0,052 | 0,683 ± 0,086 (+3,7) | 0,753 ± 0,122 (+2,0) | 0,754 ± 0,113 (+2,2) | 0,890 ± 0,242 (+0,5) |

- Wagi **bez koła** (cechy zdarzenia, pierścień ptcone100 − ptcone40, 18 cech) zdejmują około 2/3 deficytu pT-całkowo, ale nie domykają binów 5–8 i 8–30 GeV w J12only (pull +2,6 i +3,0) ani 8:j20&!J12 (+2,0…+2,2). Kryterium TAK niespełnione. Zgodne z TRKJET2 („deficytu nie zdejmuje zmienna zdarzenia”).
- Wagi **z kołem** domykają wszędzie (|pull| ≤ 0,6), z konstrukcji; ich błąd (10–24 %, zdominowany przez rozrzut wag) jest WIĘKSZY niż błąd pomiaru ε wprost z tagu (6–10 %) i większy niż obecna systematyka 8–9 % na bin. Przeważanie nic tu nie dodaje do pomiaru bezpośredniego.
- Przenoszenie wag między klasami (T6) zawodzi: 9:MB → 8:J12only 0,609 ± 0,050 (pull +7,7; z kołem 0,616 ± 0,120, +3,2); 8:J12only → 9:MB 1,233 ± 0,077 (−3,0). Działa tylko między strumieniami dżetowymi (8:J12only → 8:j20&!J12: 0,911 ± 0,127; z kołem 0,989 ± 0,166). Wagi są własnością strumienia.

### 5.2 Score H (PRE), punkty WP95/80/50

Wszystkie 48 punktów (9:MB i 8:all × 5–8, 8–30 GeV × 4 zestawy wag) mają |pull| ≤ 1,6. Ale: (i) domykają już same wagi pT; (ii) przeważenie zmienia ε_MC o ≤ 0,02 (wyjątek: wagi F18 w 8–30 GeV, +0,03…+0,06, w stronę od danych); (iii) czułość testu to 10–13 % w 5–8 GeV i 3–5 % w 8–30 GeV (8:all), czyli statystyka tagu. Punkt 8:all 5–8 WP95 zostaje na 0,84 ± 0,10 przy każdym zestawie wag — przeważenie nie rusza znanego SF_H = 0,83 ± 0,10.
Wniosek: dla H metoda formalnie domyka, lecz nie zmniejsza systematyki — ta jest ograniczona statystyką tagu, a nie modelem sygnału.

### 5.3 Przenoszenie D⁰←D* → inkluzywne D⁰ (MC fd_nd, prawda) — `fig5_przenoszenie.png`, T7

[ε(D⁰←D*, przeważone w pT) − ε(inkluzywne)]/ε:

| cięcie | 3–5 GeV | 5–8 GeV | 8–30 GeV |
|---|---|---|---|
| ptcone40 < 1 GeV | −1,9 % | −7,1 % | −32 % |
| ptcone40 < 2 GeV | +0,7 % (1,6σ) | −1,2 % (−1,7σ) | **−10,0 % (−5,9σ)** |
| ptcone40 < 4 GeV | +0,2 % | +0,6 % | −2,4 % |
| ptcone40/pT < 0,5 | +0,6 % | +0,3 % | −0,9 % |
| score H, WP80 | −1,1 % | −1,1 % | −1,0 % |

Wolny pion siedzi w stożku D⁰ i przy wysokim pT sam zjada próg bezwzględny. Wyjęcie π_s ze stożka pogarsza sprawę (+6…+18 %: D⁰ z D* są wtedy bardziej izolowane niż inkluzywne). Cięcie względne przenosi się w 1 %. Score H: −1 % (WP80), −3 % (WP50) w 3–8 GeV.

## 6. Embedding (część B)

### 6.1 Czy TRKJET2 już odpowiada?
Częściowo. TRKJET2-T1 (werdykt: „NIE piedestał”) mierzy stożki prostopadłe/losowe w zdarzeniach z kandydatem D. Własny odczyt z `tj2_*.npz` (`s6_perp.py`, T8), ε(Σ pT < 2 GeV), H06, 3–30 GeV:

| klasa | własny stożek (bez π_s) | prostopadły | losowy | Δ_corr = prostop. − własny | Δ_corr w MC |
|---|---|---|---|---|---|
| 9:MB | 0,723 ± 0,023 (MC 0,736) | 0,901 ± 0,017 (MC 0,891) | 0,872 ± 0,018 | 0,178 ± 0,029 | 0,155 |
| 8:J12only | 0,387 ± 0,016 (MC 0,594) | 0,799 ± 0,014 (MC 0,884) | 0,808 ± 0,014 | 0,412 ± 0,021 | 0,290 |
| 8:j20&!J12 | 0,488 ± 0,045 (MC 0,679) | 0,731 ± 0,039 | 0,791 ± 0,038 | 0,243 ± 0,060 | 0,211 |

Rozkład deficytu w J12only (0,594 − 0,387 = 0,207): piedestał 0,085, aktywność skorelowana z D⁰ ponad MC 0,122. Rośnie z pT: Δ_corr 0,155 / 0,354 / 0,592 w binach (MC 0,098 / 0,192 / 0,448).
TRKJET2 nie daje: nośników bez wymagania kandydata D (nieobciążonych selekcją), kierunku „jak prawdziwe D⁰ względem dżetu torowego”, kinematyki D⁰ z MC. Stąd paczka — zgodnie z K(B).

### 6.2 Paczka i test
- `Q_lxplus/sigemb/` (wzorzec trkjet2; jądro `trkjet2/tj2_kernel.py` bez zmian): dawcy = 39 865 prawdziwych D⁰ z MC (`emb_donors.npz`), nośniki = zdarzenia ROOT v2.0 (data_428770 co 20., data_427929 co 40., mb_nd co 40.), warianty kierunku: losowy (rnd), prostopadły do dżetu wiodącego (perpJ), przy dżecie torowym z rozkładu ΔR prawdziwych D⁰ (nearJ). Wyniki na EOS `charm_q/sigemb`, pobranie do `SIG_constr/runs/lxplus/`. Etap `ORCH/stages/sigemb.sh`.
- Test przez `bin/lx` w /tmp (1–3 pliki, usunięte): 428770 — 8448 nośników ze 100 592 wpisów w 20 s; 427929 — 3503 ze 136 304.
- Analiza `s7_embed.py`: przeważenie rnd ∪ nearJ do rozkładów sygnału z tagu: WO (zmienne dżetów POZA stożkiem — bez koła), WB (+ odległość do najbliższego dżetu), WA (+ pT tego dżetu — koło).

### 6.3 Wynik pełnego biegu `sigemb` (`out/s7_embed.md`, `fig6_embedding.png`)

Nośników: 9:MB 596 222, 8:J12only 219 607, 8:j20&!J12 57 071, MC mb_nd 480 468. H06, ε(Σ pT < 2 GeV w ΔR < 0,4).

| klasa, pT | ε_tag | rnd | perpJ | WO (bez koła) | WB | Δ(rnd) | Δ(WO) | Δ(WB) |
|---|---|---|---|---|---|---|---|---|
| 9:MB 3–30 | 0,723 ± 0,023 | 0,962 | 0,972 | 0,794 ± 0,003 | 0,750 ± 0,020 | 0,239 ± 0,023 | 0,071 ± 0,023 | 0,028 ± 0,030 |
| 8:J12only 3–30 | 0,387 ± 0,016 | 0,809 | 0,856 | 0,597 ± 0,004 | 0,408 ± 0,014 | 0,422 ± 0,016 | 0,210 ± 0,016 | 0,022 ± 0,021 |
| 8:J12only 3–5 | 0,619 ± 0,049 | 0,809 | 0,856 | 0,608 ± 0,010 | 0,633 ± 0,043 | 0,189 ± 0,049 | −0,012 ± 0,050 | 0,014 ± 0,065 |
| 8:J12only 5–8 | 0,428 ± 0,030 | 0,809 | 0,856 | 0,582 ± 0,007 | 0,467 ± 0,028 | 0,380 ± 0,030 | 0,154 ± 0,031 | 0,038 ± 0,041 |
| 8:J12only 8–30 | 0,229 ± 0,013 | 0,809 | 0,856 | 0,597 ± 0,003 | 0,248 ± 0,012 | 0,579 ± 0,013 | 0,367 ± 0,014 | 0,019 ± 0,018 |
| 8:j20&!J12 3–30 | 0,488 ± 0,045 | 0,840 | 0,875 | 0,630 ± 0,005 | 0,582 ± 0,040 | 0,352 ± 0,045 | 0,142 ± 0,046 | 0,094 ± 0,060 |
| MC mb_nd ← fd_nd 3–30 | 0,730 ± 0,002 | 0,963 | 0,976 | 0,793 ± 0,001 | 0,758 ± 0,002 | 0,233 ± 0,002 | 0,063 ± 0,002 | 0,028 ± 0,003 |

- Δ(rnd) > 3σ we wszystkich klasach → aktywność w stożku jest skorelowana z D⁰ (dżet c), nie jest tłem zdarzenia. Kontrola ε(rnd) ≥ ε_tag spełniona.
- Kontrola 9:MB dane vs MC (ta sama procedura): rnd 0,962 / 0,963, WO 0,794 ± 0,003 / 0,793 ± 0,001, WB 0,750 ± 0,020 / 0,758 ± 0,002 — zgodne w 2σ. Procedura jest poprawna; w MC znane jest też obciążenie własne: Δ(WO) = 0,063, Δ(WB) = 0,028 ± 0,003 (9σ od zera — WB też nie domyka ściśle).
- Wariant bez koła (WO) nie domyka: 3,1σ w 9:MB, 13σ w J12only (w binie 8–30 GeV Δ = 0,37) → wg K(B) embedding **nie jest kalibracją**; mierzy składową nieskorelowaną.
- Nośniki nieobciążone mają ε(rnd) 0,96 w 9:MB wobec 0,87–0,90 w zdarzeniach z kandydatem D (TRKJET2) — zdarzenia z D są o ~0,07 bardziej zajęte już poza stożkiem.

### 6.4 Czy WB wolno nazwać kalibracją? (`s7b_wb_circ.py`, `out/s7b_wb_circ.json`)

WB przeważa embedding do rozkładu (pT, ΔR) najbliższego dżetu torowego pT > 2 GeV, który dla prawdziwego D⁰ jest zbudowany z jego towarzyszy. Dżet o pT > 2 GeV w ΔR < 0,4 („IN”) ma swoje tory w stożku, więc sam przesądza o wyniku cięcia Σ pT < 2 GeV.

| | 8:J12only | 9:MB | MC fd_nd |
|---|---|---|---|
| udział sygnału z dżetem IN (f_IN) | 0,564 | 0,236 | 0,236 |
| f_IN dla embeddingu losowego | 0,143 | 0,029 | – |
| ε_tag przy IN | 0,039 ± 0,007 | 0,046 ± 0,028 | 0,062 |
| ε_tag przy OUT | 0,855 ± 0,020 | 0,926 ± 0,016 | 0,936 |
| ε^embed przy OUT (WB) | 0,895 | 0,958 | – |
| Δ(OUT) = embed − tag | +0,040 ± 0,020 | +0,032 ± 0,016 | – |
| część wariancji pass/fail ustalona przez komórkę WB | 0,79 | 0,78 | – |
| udział sygnału w komórkach z ε < 0,1 lub > 0,9 | 0,80 | 0,83 | – |

- WB jest kołowe w ~80 %: ε ≈ (1 − f_IN)·ε_OUT, a f_IN pochodzi z samego tagu (J12only: 0,436 × 0,855 + 0,564 × 0,039 = 0,395 wobec zmierzonego 0,387). Zgodność Δ(WB) = 0,022 ± 0,021 jest więc w większości tożsamością, nie testem.
- Treść niekołowa to wyłącznie część OUT: tam embedding jest o 0,03–0,04 za wysoki (2σ w każdej klasie, ten sam znak co w MC) — resztkowa aktywność skorelowana poniżej progu dżetu 2 GeV.
- Co WB jednak pokazuje: deficyt J12only to geometria D⁰–dżet. f_IN rośnie z 0,236 (9:MB i MC — identycznie) do 0,564 w J12only; przy f_IN z MC J12only miałoby ε = 0,66. Z 0,34 różnicy do MC około 0,27 (80 %) to wzrost f_IN, reszta (0,07) to niższe ε_OUT (piedestał 0,936 → 0,855).
- Werdykt wg `KRYTERIA.md`: kryterium kalibracji dotyczy wariantu bez koła (WO) i jest niespełnione. WB **nie jest kalibracją** — dla inkluzywnego D⁰ wymagałoby f_IN z danych, czyli znów tagu D* (z tym samym problemem przenoszenia przy 8–30 GeV, §5.3). WB jest **dowodem mechanizmu**: MC opisuje stożek przy zadanej geometrii (ε_IN, ε_OUT zgodne w ~0,03–0,08), a nie opisuje, jak często D⁰ siedzi w dżecie w próbce z trygerem dżetowym. Praktyczny wniosek: jedyną wielkością do wzięcia z danych per strumień jest f_IN (jedna liczba na bin pT), co jest równoważne pomiarowi ε wprost.

## 7. Wnioski i rekomendacja (część C)

1. **(a)** „Sygnał z konstrukcji” nie zwiększa istotności klasyfikatora H: sygnał MC w 18 cechach jest wystarczająco dobry, a próbka tagu 50–100× za mała. Zysk +5…+12 % (> 3σ) istnieje wyłącznie dla modelu z izolacją w strumieniu dżetowym — i jest w istocie naprawą złej izolacji w MC, nie nową informacją.
2. **(b)** Przeważanie sygnału MC do tagu nie zmniejsza systematyki ε_S. Dla izolacji wagi bez koła nie domykają, wagi z kołem są gorsze od pomiaru wprost, a wagi nie przenoszą się między strumieniami. Dla H nie ma czego poprawiać w granicach statystyki tagu.
3. **(c)** Kalibracja izolacji z tagu D* jest używalna pod dwoma warunkami: (i) osobno per strumień triggera, (ii) zmienna względna ptcone40/pT albo zakres pT < 8 GeV; dla cięcia bezwzględnego w 8–30 GeV trzeba dodać systematykę 10 % (albo poprawkę z MC). Niesprawdzone: przenoszenie D*-tag → inkluzywne dla samego wyjścia sieci z izolacją (potrzebny model na MC fd_nd; ~1 h lokalnie).
4. **(B)** Embedding nie zastąpi tagu: daje piedestał na zdarzeniach nieobciążonych i rozdziela deficyt J12only na piedestał (~0,085) i aktywność skorelowaną (~0,12 ponad MC). Pełny bieg (§6.4): deficyt to geometria D⁰–dżet — udział D⁰ z dżetem torowym > 2 GeV w stożku rośnie z 0,24 (9:MB, MC) do 0,56 (J12only) i tłumaczy ~80 % różnicy; wariant WB to pokazuje, ale jest kołowy, więc to materiał do opisu mechanizmu i systematyki, nie kalibracja.

Rekomendacja: zostać przy modelu H z MC i SF z tagu; ε_iso brać z danych per strumień w zmiennej względnej; jeżeli sieć typu CharmNet ma używać izolacji w próbce 428770 — dostrojenie na tagu D* tego strumienia (wariant STACK), z ε_S mierzonym OOF na tym samym tagu (błąd 3 % pT-całkowo). Koszt: trening lokalny minuty; pełny embedding ≤ 6 h zegarowych na lxplus.

## 8. Co wymaga jobów lxplus

Etap `sigemb` (wykonany, wyniki pobrane 17.11 08:41):
```
  - name: sigemb
    needs: [trkjet2]
    max_hours: 6
```
Analiza po pobraniu wykonana: `s7_embed.py` → `out/s7_embed.{json,md}`, `s7b_wb_circ.py` → `out/s7b_wb_circ.json`, `s8_figs.py` (fig6 z pełnego biegu). Dalszych jobów nie potrzeba.

## 9. Zastrzeżenia

- Lokalny komputer był przeciążony: xgboost w ustawieniach zredukowanych (eta 0,25, ≤ 150 drzew, max_bin 64, tło ≤ 40 tys., 100 replik bootstrapu). Krzywa uczenia i Z_rel są wewnętrznie spójne (te same ustawienia), wartości bezwzględne AUC mogą być o ~1 % niższe niż przy pełnym treningu.
- Błąd wag w domknięciu to jackknife po 5 modelach foldów (przybliżenie; dla wag z kołem dominuje).
- ε_B liczone z górnego pasma Δm; tło pod pikiem może mieć nieco inny skład.
- sWeights niewiarygodne poniżej czystości ~0,5; wszystkie wyniki na SBS.
- Rozdział §5.1: 9:MB → brak deficytu, więc domknięcie tam nie testuje metody, tylko brak szkody.

## 10. Pliki

Skrypty: `sc_common.py`, `sc_ml.py`, `s1_neff.py`, `s2_learn.py`, `s3_variants.py`, `s4_closure.py`, `s5_transfer.py`, `s6_perp.py`, `s7_embed.py`, `s7b_wb_circ.py`, `s8_figs.py`, `s9_tables.py`, `run_chain2.sh`.
Rysunki: `fig1_neff.png`, `fig2_krzywa_uczenia.png`, `fig3_warianty.png`, `fig4_domkniecie.png`, `fig5_przenoszenie.png`, `fig6_embedding.png`.
Paczka: `Q_lxplus/sigemb/`, etap `ORCH/stages/sigemb.sh`. Test: `runs/emb_test/`.
