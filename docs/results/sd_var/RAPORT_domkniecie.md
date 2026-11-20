# SD_var — domknięcie: brakujące części CRqcd i przeliczone przewidywania liczby tagowanych D (projekt B, 19/20.11.2026)

Poziom cząstek (Pythia 8.317) + uzyski INKLUZYWNE z danych. **Żadnego fitu ani liczby z okna sygnału w zdarzeniach z tagiem.** Noty, `ORCH/pipeline.yaml`, stanu etapu `sdvar` i pętli nie zmieniano.
Tabele: `out/pred_tables.md`, `out/pred.json` (skrypt `predict.py`), `out/tables.md`, `out/results.json` (`analyze_var.py`); kopia stanu sprzed domknięcia: `out.bak_przed_domknieciem/`, `analyze_var.py.bak_przed_domknieciem`.

## 0. Oświadczenie o niezależności przewidywań

- Wyników pomiaru tagowanego (`TAGGED_2026-11-19/out`, `tabele`, `fig`, `runs`, `RAPORT.md`) **nie otwierałem**; przeczytałem tylko `TAGGED_2026-11-19/KRYTERIA.md` (definicje; plik zawiera też przeliczenie dawnych przewidywań na ułamki — to nie wynik).
- Metoda i pierwsze liczby przewidywań (`predict.py`, bieg próbny) powstały **19.11.2026 o 23:49 CET**; ostateczny bieg 20.11.2026 o 02:01 CET po dosłaniu części CRqcd (zmiana obwiedni < 0,3 %).
- Zmiany po 23:49 wynikają wyłącznie z trzech wejść spoza regionu tagowanego: (i) dosłane części CRqcd; (ii) decyzja z 20.11: pomiar A·ε_AFP = 0,92 ± 0,05 (`AFP_akceptancja_2026-11-20/PLAN.md`) — podane obie wersje; (iii) decyzja z 20.11: wycofanie „BR w MC ≠ PDG” (`CVER_2026-11-20/RAPORT.md`) — usunięta skala κ × 1,364 w uzysku D⁰.
- **Co zobaczyłem mimochodem (jawnie):** ok. 00:05 CET, szukając notatki do aktualizacji w bazie wiedzy, przeczytałem notatkę XSEC2 (studium przekroju inkluzywnego, nie TAGGED): uzyski D w zdarzeniach z bitem `HLT_noalg_L1AFP_A_OR_C` (R1, bez wymogu tagu offline): D⁰ 45 ± 30 / 9 ± 10 / 5 ± 4, D⁺ 12 ± 14 / 0 ± 6 / 5 ± 6, oraz P(L1_AFP | zdarzenie z D) = 0,41 ± 0,22 % (428770) i ≈ 0,5 ± 0,3 % (427929). To liczby z danych pokrewne przewidywaniu dla R1. **Nie użyłem ich i niczego do nich nie dopasowałem** (rachunek R1 był już zapisany; po tym zmieniłem w nim tylko czynnik ε_AFP 0,94/0,95 z pkt. (ii), −1 %).
- **20.11 ok. 02:20 CET — po zamknięciu liczb (bieg ostateczny 02:01) i po napisaniu §1–§5 tego raportu** — przy wyszukiwaniu (`grep`) notatki SD_var w bazie wiedzy wyświetliła mi się **nazwa pliku** notatki studium tagowanego (TAGXSEC) zawierająca w tytule wynik: „R_fid 0,99 procent (wariant 0,54), sigma 13,6 µb, rho 0,96 ± 0,12”. Pliku nie otworzyłem. Żadna liczba przewidywań nie została po tym zmieniona (znaczniki czasu `out/pred.json`, `out/results.json`: 02:01).
  Plan `AFP_akceptancja` odnotowuje, że jego autor widział tytuł notatki TAGGED z wynikiem; ja do 02:20 żadnego tytułu z wynikiem tagowanym nie widziałem (lista `git status` z 00:05 zawierała tylko notatki INCAL i XSEC2); o 02:20 — jak wyżej.

## 1. Części CRqcd: co wysłano, stan, przyczyny

**Kolejka przed wysyłką (19.11, 23:43):** żadnych jobów `charm_sdvar` — cztery wiszące ponowienia z RAPORT §7.1 zostały zdjęte 19.11 o 14:38 razem z DAG-iem (`sdvar.dag.dagman.out`: „Received SIGUSR1 … Aborting DAG”, rescue001). W kolejce były tylko `charm_dnet6` i `charm_tagged` — nietknięte.

**Wysłane (19.11, 23:45; własny DAG `Q_lxplus/sd_var/closure/sdclose.dag`, batch `charm_sdvar_close`, 9 węzłów, 1 rdzeń / 2 GB, `+MaxRuntime` 1,5 h (dv) i 4 h (nv), `RETRY 0`, submit w czystym env, `-maxjobs 10`):**

| węzeł | ziarno | wynik | CPU |
|---|---|---|---|
| `g_dv_CRqcd_101` (1 M zd.) | 1230101 | OK | 8,7 min |
| `g_nv_CRqcd_101…106` (6 × 1 M zd.; 4 brakujące + 2 zapasowe) | 1260101–1260106 | **6/6 OK**, żaden nie zawisł | 1,3–2,2 h |
| `d_dv_CRqcd_8` (diagnostyka: stare ziarno 1230008 pod gdb) | 1230008 | segfault odtworzony, ślad stosu | 6,5 min |
| `d_nv_CRqcd_8` (diagnostyka: stare ziarno 1260008 pod gdb, SIGINT po 3 h 10 min) | 1260008 | zawieszenie odtworzone (3 h 08 min CPU bez końca), ślad stosu po SIGINT | 3,1 h |

Paczka `sd_var` użyta bez zmian (`sd_job.sh`, `bin/pythia_sd`); nowe pliki tylko w `sd_var/closure/` (`sdclose.sub`, `sdclose.dag`, `sd_diag.sh`). Bez węzła merge: scalenie dwóch modeli zrobione ręcznie tym samym `sd_merge.py` (kilka sekund) do **nowego** katalogu EOS `charm_q/sd_var/merged_close/` → lokalnie `runs/lxplus/merged_close/`; dotychczasowe `merged/` i `DONE` nietknięte.

**Stan końcowy:** `dv_CRqcd` **10/10 części** (10,0 M zd., 7152 D z tagiem), `nv_CRqcd` **22 części** (22,0 M zd., 3644 D z tagiem; plan 20 — użyto wszystkich, które się skończyły, reguła przyjęta przed wysyłką; dobór „po ukończeniu” jest ten sam co w pierwszym biegu). Jedna próba wystarczyła.

**Przyczyny (z logów condora i z diagnostyki):**
- **Segfault `g_dv_CRqcd_8`: błąd Pythii 8.317 w rekoneksji QCD (`ColourReconnection:mode = 1`), nie pamięć i nie węzeł.** Ślad stosu (gdb, ziarno 1230008): `ColourReconnection::getJunctionIndices` ← `mDip` ← `makeAllPseudoParticles` ← `nextNew` ← `BeamRemnants::add` ← `PartonLevel::next` ← `Pythia::next`. Stos ma 8 ramek (nie przepełnienie stosu przez rekurencję), RSS 6 MB / zużycie wg condora 15 MB przy limicie 3 GB; trzy próby na trzech różnych węzłach + czwarta pod gdb padają tak samo po 5–7 min CPU → jedno konkretne zdarzenie wyznaczone przez ziarno (odczyt nieprawidłowego dipola/złącza). Częstość: 1 zdarzenie na 1,1·10⁷ dyfrakcyjnych z CRqcd.
- **Zawieszenia `nv_CRqcd` (ziarna 1260008, 12, 13, 17):** z logu condora — joby zużywały CPU do końca (0,99 rdzenia), pamięć stała 12–13 MB, zdjęte przez `SYSTEM_PERIODIC_REMOVE` po 8 h (limit `workday`); ponowienia tych samych ziaren znów liczyły > 20 min bez końca. To pętla obliczeniowa w generatorze, nie pamięć, nie I/O i nie tempo węzła. **Diagnostyka (ziarno 1260008 pod gdb, trzecie uruchomienie tego ziarna): znów nie kończy; po 3 h 08 min CPU (zwykły job: 1,3–2,2 h) ostatnie zapisane zdarzenie z D ma numer 736 482 z 10⁶, a proces stoi w `JunctionSplitting::splitJunGluons` (→ `StringZ::zLund` → `Rndm::flat`) ← `JunctionSplitting::checkColours` ← `BeamRemnants::addNew` ← `BeamRemnants::add` ← `PartonLevel::next`; RSS 11 MB.** Czyli nieskończona (albo praktycznie nieskończona) pętla losowania przy rozdzielaniu gluonów między złączami w nowym modelu resztek wiązki (`BeamRemnants:remnantMode = 1`, wymagany przez CR mode 1) — deterministyczna dla ziarna. Zastrzeżenie: to jedna próbka stosu jednego ziarna; pozostałych trzech ziaren nie powtarzałem.
  Częstość: 4 ziarna na 26 jobów ND z CRqcd (2,6·10⁷ zdarzeń), żadne z 6 nowych; w próbkach bez CRqcd (`nv_pT0lo/hi`, 40 jobów) zero. Konfiguracja wspólna dla obu awarii: `crqcd.cmnd` (`ColourReconnection:mode = 1`, `BeamRemnants:remnantMode = 1`).
- Wpływ na R: zdarzenia gubione przez awarie to ≤ 1 na 5·10⁶; pierwszy bieg tracił przez nie całe joby (losowo względem fizyki), więc R nie jest obciążone. Potwierdza to zgodność przed/po: CRqcd 0,512 ± 0,008 → 0,512 ± 0,007 %, NDCRqcd 0,561 ± 0,008 → 0,563 ± 0,008 %.

**Skutek dla obwiedni:** D⁰ × 0,40 … × 3,97 → **× 0,40 … × 3,98** (R = 0,384 ± 0,005 %, zakres 0,154–1,529 %); D⁺ × 0,38 … × 3,79; D* × 0,39 … × 3,98. Grupa L (ND CRqcd) × 1,57 → × 1,58. Kryteria i zbiory bez zmian; `analyze_var.py` dostał tylko opcję `--over` (nadpisanie modeli z `merged_close`).

## 2. Jak liczone są przewidywania i czy wchodzi w nie BR albo C z MC

Oba scenariusze mają postać **N_tag = Σ_bin N_incl(dane) × ułamek**:
- **M (generator):** ułamek = R_bin(Pythia, poziom cząstek, gatunek D bez wymogu kanału rozpadu) × A·ε_AFP × ε_ξcal;
- **S (pasma boczne):** ułamek = udział kandydatów pasm bocznych w regionie T2 (dane; liczby kandydatów) × (1 − f_acc), z przypisanym czynnikiem 0,5–2.

| droga | czy zawiera C_XSEC × BR_PDG albo uzysk MC | skutek wycofania „BR w MC ≠ PDG” (CVER: czynniki 0,734 / 0,939) |
|---|---|---|
| M, wszystkie gatunki, 427929 i R1 | **nie** — R jest stosunkiem dla całego gatunku (BR się skraca), N_incl to uzysk z danych | **bez zmian** |
| S, wszystkie gatunki, 427929 i R1 | **nie** — same dane | **bez zmian** |
| przeniesienie 427929 → R1 | nie (stosunki świetlności, wydajności L1 AFP, weta, ξ_cal) | bez zmian |
| uzysk inkluzywny D⁰ modelem M8 | pośrednio: składowa pikująca (inne mody D⁰ + inne hadrony powabne) była skalowana κ × 1,364 | **skala usunięta** (κ_dane = κ_MC): N_incl D⁰ WP70 6316 → 6330 (per bin −0,7 / +1,2 / +3,5 %), WP80 8042 → 8041; przewidywania +0,1 % / −0,2 % |
| uzysk inkluzywny D⁺ (M0) | nie (brak składowej pikującej) | bez zmian |
| dawne liczby noty (N_incl z fitu C0, κ z MC względem sygnału prawdy) | nie przez BR; sygnał prawdy w MC to tryb ekskluzywny (CVER pkt 48) | bez zmian |

Żadne przewidywanie liczby tagowanych D nie jest liczone drogą „Pythia/FONLL × C_XSEC × BR_PDG × L”, więc **żadne nie jest zaniżone o 1/0,734 ani 1/0,939 i żadnej poprawki BR nie stosuję**. (Zaniżone o te czynniki są tylko porównania przekrojów w XSEC — poza tym raportem.)

## 3. Czynnik detektorowy: stary 0,55 i nowy 0,92

**Stary (TAGREG §3, 18.11):** A·ε_AFP = 0,55 (zakres 0,35–0,80) — **założenie**, złożone z: tory NEAR × FAR 0,90–0,95 (literatura), dopasowanie/rekonstrukcja protonu > 0,9 (literatura), akceptancja geometryczna w t / pT protonu **przyjęta** 0,6–0,9 (brak źródła). 0,92 × 0,9 × 0,67 ≈ 0,55. Osobno mnożone ε_ξcal = 0,99.
Czego stary czynnik **nie** zawierał: weta `nPileupVertices == 0`, wymogu „dokładnie jeden proton” (jawnie), wydajności trygera, tagów przypadkowych.

**Nowy (`AFP_akceptancja_2026-11-20/PLAN.md` §0, §2, §3, §5):** A = 0,97 ± 0,03 (fast-sim z optyką 2022 + krawędzie z danych; okno 0,03–0,10 bez cięcia t) × ε_AFP = 0,95 ± 0,02 (tag-and-probe near/far; **zawiera** tor w obu stacjach i wymóg „dokładnie jeden proton”) = **0,92 ± 0,05**.

| składnik | w starym 0,55 | w nowym 0,92 | w rachunku tutaj |
|---|---|---|---|
| akceptancja geometryczna (ξ, t) | tak (0,6–0,9, założona) | tak (0,97) | raz |
| wydajność torów near/far i rekonstrukcji protonu | tak (≈ 0,83) | tak (w 0,95) | raz |
| wymóg „dokładnie jeden proton” (showery, drugi tor) | niejawnie | **tak** (w 0,95) | raz — nie dodaję osobno |
| cięcie ξ_cal | nie (osobno 0,99) | nie | osobno: stara wersja 0,99; nowa 1,00 dla 427929 (pomiar TAGREG 1,000 +0/−0,005); dla R1 stosunek 0,974 w czynniku k |
| weto npu = 0 | nie | nie | **nie mnożę**: skraca się w ułamku N_tag/N_incl, gdy oba mają weto; uzyski niżej są dla N_incl bez weta → przy wecie w obu × 0,985 (427929). Dla R1 stosunek wet 0,941/0,985 jest w k |
| tryger (MBTS przy tagu 0,9998) | nie | nie | pominięty (≤ 0,1 %); dla R1 wydajność L1 AFP w k |
| tagi przypadkowe | nie | nie | osobna kolumna (dodają się do prawdziwych) |

Nowy/stary = 0,92 × 1,00 / (0,55 × 0,99) = **1,690** (sama A·ε: × 1,673; ε_ξcal 0,99 → 1,00: × 1,010). Scenariusz S nie zależy od tego czynnika.
Uwaga: nowa A dotyczy także protonów wiodących ND (A = 0,86 dla ND jest już w 0,97); niepewność ± 5 % jest pomijalna wobec obwiedni modeli.

## 4. Tabela końcowa przewidywań (prawdziwe tagi; pT 3–30 GeV; klasyfikator H)

Składniki wspólne: R nominalne (oficjalne A3) w binach 3–5 / 5–8 / 8–30 GeV: D⁰ 0,375 / 0,324 / 0,268 %, D⁺ 0,376 / 0,306 / 0,268 %, D* 0,376 / 0,318 / 0,292 %; obwiednia modeli D⁰ × 0,40 … × 3,98, D⁺ × 0,38 … × 3,79, D* × 0,39 … × 3,98;
ułamek w pasmach bocznych P_SB: D⁰ 0,365 / 0,354 / 0,470 %, D⁺ 0,233 / 0,232 / 0,324 %, D* 0,229 / 0,256 / 0,210 %; f_acc = 0,0087 (427929), 0,038 (R1; zakres 0,038–0,086).
Uzysk inkluzywny 427929 (nasza kopia, L = 0,4648 nb⁻¹, bez weta npu): D⁰ **model M8 bez skali κ**: WP70 4077 / 1548 / 705 = 6330 ± 496, WP80 5508 / 1690 / 843 = 8041 ± 753; D⁺ model M0: WP70 5885 ± 353, WP80 6583 ± 516; D* z tagu (Δm): 885 ± 43 (bez zmian; selekcja WP80).

### 4.1 Ułamek tagowany (niezależny od punktu pracy) [%]

| gatunek | S: pasma boczne (× 0,5–2) | M: generator, stary czynnik 0,55 × 0,99 (obwiednia modeli) | M: generator, nowy czynnik 0,92 × 1,00 |
|---|---|---|---|
| D⁰ | 0,370 (0,185–0,74) | 0,191 (0,077–0,76) | **0,322 (0,129–1,28)** |
| D⁺ | 0,241 (0,121–0,48) | 0,187 (0,072–0,71) | **0,317 (0,121–1,20)** |
| D* | 0,233 (0,116–0,47) | 0,190 (0,073–0,76) | **0,320 (0,124–1,28)** |

### 4.2 Liczba tagowanych D

| gatunek | run / region | WP | N_incl | S: pasma (× 0,5–2) | M stary czynnik | **M nowy czynnik** | przypadkowe | droga liczenia; wpływ wycofania BR |
|---|---|---|---|---|---|---|---|---|
| D⁰ | 427929, MBTS | 70 | 6330 | 23,5 (11,7–47) | 12,1 (4,8–48) | **20,4 (8,2–81)** | 0,21 | ułamek × uzysk z danych; bez zmian (κ: +0,1 %) |
| D⁰ | 427929, MBTS | 80 | 8041 | 29,8 (14,9–60) | 15,5 (6,2–61) | **26,1 (10,5–104)** | 0,26 | jw. (κ: −0,2 %) |
| D⁺ | 427929, MBTS | 70 | 5885 | 14,2 (7,1–28) | 11,0 (4,2–42) | **18,6 (7,1–71)** | 0,12 | jw.; bez zmian |
| D⁺ | 427929, MBTS | 80 | 6583 | 15,8 (7,9–32) | 12,3 (4,7–47) | **20,8 (8,0–79)** | 0,14 | jw.; bez zmian |
| D* | 427929, MBTS | (80) | 885 | 2,1 (1,0–4) | 1,7 (0,6–7) | **2,8 (1,1–11)** | 0,02 | jw.; bez zmian |
| D⁰ | 428770, R1 (0,4687 nb⁻¹) | 70 | — | 14,9 (7,4–30) | 7,7 (3,1–31) | **13,0 (5,2–52)** | 0,5–0,6 | 427929 × k; bez zmian |
| D⁰ | 428770, R1 | 80 | — | 18,9 (9,4–38) | 9,8 (3,9–39) | **16,6 (6,6–66)** | 0,7 | jw. |
| D⁺ | 428770, R1 | 70 | — | 9,0 (4,5–18) | 7,0 (2,7–27) | **11,8 (4,5–45)** | 0,4–0,5 | jw. |
| D⁺ | 428770, R1 | 80 | — | 10,1 (5,0–20) | 7,8 (3,0–30) | **13,2 (5,0–50)** | 0,4–0,5 | jw. |
| D* | 428770, R1 | (80) | — | 1,3 (0,7–3) | 1,1 (0,4–4) | **1,8 (0,7–7)** | 0,05–0,07 | jw. |

- Przedziały: S — czynnik 0,5–2 (założenie „sygnał tagowany tak często jak tło”); M — obwiednia modeli dyfrakcji (rozrzut modeli, nie przedział ufności). Z dawnym zakresem A·ε 0,35–0,80 dołożonym w kwadraturze (ln) stara wersja M dla D⁰ WP80 to 5,6–65; w nowej wersji ± 5 % nie zmienia zakresu.
- Nie wliczone w przedziały: błąd N_incl (D⁰ ± 8–9 % stat.; domknięcie fitu M8 na MC 25 / 8 / 9 % w binach, ≈ 19 % w sumie; D⁺ ± 6–8 %), C_D(dyfrakcja)/C_D(ND) (wydajność klasyfikatora w zdarzeniach z przerwą — niesprawdzona).
- **R1 = 427929 × k, k = 0,635 (zakres 0,60–0,77):** L 0,4687/0,4648 = 1,008 × wydajność L1 AFP dla tagu offline 0,684 (zdarzenia T2 z kandydatem D⁰ w 427929 z bitem L1 AFP: 7139/10438; strona A ≈ 0,91, C ≈ 0,46) × weto npu 0,941/0,985 = 0,955 × ξ_cal 0,974/1,000 × ε_AFP 0,94/0,95.
  Dół zakresu: weto jak w klasie „inne trygery” (0,88); góra: stosunek liczby zdarzeń T2 z kandydatem D⁰ (8267/10438) × (1 − f_acc). Założenia: ta sama wydajność L1 AFP w obu runach; zdarzenia z D obecne w drzewie R1 (drzewo ma 44 % zdarzeń tego łańcucha — wymóg wierzchołka, dla zdarzeń z D ≈ 1).
  Tagi przypadkowe w R1: 4–9 % prawdziwych. Uzysk inkluzywny w R1 nie istnieje (tryger wymaga protonu) — tylko N_tag i N_tag/L.
- **R2 (`j20_L1AFP_A_OR_C` ∧ ¬R1, 27,65 nb⁻¹) — opisowo.** Liczby z modelu nie podaję: R dotyczy próbki minimum-bias, a wydajność `j20` dla zdarzeń z D i przerwą nie jest znana (w R1 bit `j20` ma 13 % zdarzeń z D⁰). Kierunek: R maleje z pT (0,375 → 0,268 %; w scenariuszu twardym × 4,9 mniej przy 8–30 niż przy 2–3 GeV), a wymóg dżetu wybiera zdarzenia o dużej aktywności, więc ułamek tagowany powinien być mniejszy niż w minimum-bias.
  Kotwica z pasm dla całego 428770 (N_incl „all” WP70: D⁰ 6336 / 4834 / 7716, D⁺ 4118 / 3451 / 6458; P_SB 0,93 / 0,89 / 0,68 % i 0,60 / 0,66 / 0,55 %; f_acc 0,048): D⁰ ≈ 147, D⁺ ≈ 79; R2 ma 77 % zdarzeń tagowanych z kandydatem D⁰ → **rząd 110 D⁰ i 60 D⁺ (× 0,5–2)**, z 5–9 % tagów przypadkowych; orientacyjnie, bo skład trygerowy sygnału i pasm w skimie nie musi być ten sam.

### 4.3 Co się zmieniło względem noty (Sec. 6.11; WP80, 427929)

| gatunek | nota: S | teraz: S | nota: M | teraz: M, stary czynnik | teraz: M, nowy czynnik | przyczyny |
|---|---|---|---|---|---|---|
| D⁰ | 35 (18–71) | 29,8 (14,9–60) | 18,6 (7,4–74) | 15,5 (6,2–61) | 26,1 (10,5–104) | N_incl 9596 (fit C0) → 8041 (M8): −16 %; A·ε 0,55 → 0,92: × 1,69 |
| D⁺ | 18 (9–36) | 15,8 (7,9–32) | 31 (5–120) | 12,3 (4,7–47) | 20,8 (8,0–79) | R 0,75 → 0,384 % z R(pT): 31 → 14,4; N_incl 7599 (C0) → 6583 (M0, WPSCAN): −13 %; A·ε × 1,69 |
| D* | 2 (1–4) | 2,1 (1,0–4) | 4 (0,6–14) | 1,7 (0,6–7) | 2,8 (1,1–11) | R 0,75 → 0,384 %; A·ε × 1,69; N_incl bez zmian |

Sprawdzenie ciągłości: dawne N_incl (C0) z nowym R i starym czynnikiem dają D⁰ 18,6 (liczba z noty), D⁺ 14,4, D* 1,7.
Po pomiarze A·ε scenariusze M i S są dla D⁰ zgodne co do wartości centralnej (26 wobec 30); dla D⁺ M jest powyżej S (21 wobec 16), bo P_SB(D⁺) < P_SB(D⁰), a model daje R(D⁺)/R(D⁰) = 0,98–1,05.
Przy WP70 (punkt pracy pomiaru tagowanego): D⁰ S 23,5 (11,7–47), M 20,4 (8,2–81); D⁺ S 14,2 (7,1–28), M 18,6 (7,1–71).

## 5. Tekst do noty (Sec. 6.11) — po angielsku, gotowy do wklejenia

Zastępuje akapit „Expected number of tagged D in 427929 …” i ostatnie zdania akapitu „Model uncertainty …” o brakujących częściach.

> **Expected number of tagged D (a prediction; updated 20.11).** Both scenarios have the form N_tag = Σ over pT bins of N_incl × fraction, with N_incl the inclusive yield fitted in data (our copy of run 427929, classifier H, pT 3–30 GeV; D⁰ with fit model M8 without any branching-fraction scale of the peaking component, D⁺ with the nominal model M0). In the sideband scenario the fraction is the share of sideband candidates in the tagged region times (1 − f_acc), with an assigned factor 0.5–2. In the generator scenario it is R(pT) × A·ε_AFP × ε_ξcal, with R from the official A3 configuration and the envelope of the diffraction-model variation (D⁰ × 0.40 to × 3.98, D⁺ × 0.38 to × 3.79, D* × 0.39 to × 3.98). The detector factor assumed on 18.11, A·ε_AFP = 0.55 (0.35–0.80; track and proton reconstruction about 0.83 from the literature times an assumed geometric acceptance of 0.6–0.9) with ε_ξcal = 0.99, is replaced by the measured A·ε_AFP = 0.92 ± 0.05 (fast-simulation acceptance 0.97 ± 0.03 in 0.03 ≤ ξ ≤ 0.10 times the near/far tag-and-probe efficiency 0.95 ± 0.02, which includes the requirement of exactly one proton) with ε_ξcal = 1.00 in 427929; the generator scenario rises by a factor 1.69. The pile-up-vertex veto is not part of the fraction: it cancels in N_tag/N_incl when applied to both, and the yields below refer to N_incl without the veto (× 0.985 with it). No prediction contains a Monte Carlo efficiency or a branching fraction, so the withdrawn statement "branching fractions in the MC differ from PDG" does not change them; removing the corresponding scale from the D⁰ fit changes the D⁰ yield by +0.2 % (WP70) and 0.0 % (WP80).
>
> | species | run, region | working point | N_incl | sideband scenario (× 0.5–2) | generator, detector factor 0.55 × 0.99 | generator, measured factor 0.92 × 1.00 | accidental |
> |---|---|---|---|---|---|---|---|
> | D⁰ | 427929, MBTS | WP80 | 8041 ± 753 | 30 (15–60) | 15.5 (6.2–61) | **26 (10–104)** | 0.3 |
> | D⁰ | 427929, MBTS | WP70 | 6330 ± 496 | 23 (12–47) | 12.1 (4.8–48) | **20 (8–81)** | 0.2 |
> | D⁺ | 427929, MBTS | WP80 | 6583 ± 516 | 16 (8–32) | 12.3 (4.7–47) | **21 (8–79)** | 0.14 |
> | D⁺ | 427929, MBTS | WP70 | 5885 ± 353 | 14 (7–28) | 11.0 (4.2–42) | **19 (7–71)** | 0.12 |
> | D* | 427929, MBTS | (WP80, D* tag) | 885 ± 43 | 2.1 (1.0–4) | 1.7 (0.6–7) | **2.8 (1.1–11)** | 0.02 |
> | D⁰ | 428770, R1 | WP70 | – | 15 (7–30) | 7.7 (3.1–31) | **13 (5–52)** | 0.5–0.6 |
> | D⁺ | 428770, R1 | WP70 | – | 9 (5–18) | 7.0 (2.7–27) | **12 (5–45)** | 0.4–0.5 |
> | D* | 428770, R1 | (WP80, D* tag) | – | 1.3 (0.7–3) | 1.1 (0.4–4) | **1.8 (0.7–7)** | 0.05–0.07 |
>
> The corresponding tagged-to-inclusive fractions, which do not depend on the working point, are 0.37 % (D⁰), 0.24 % (D⁺) and 0.23 % (D*) in the sideband scenario and 0.32 % (0.13–1.28 %), 0.32 % (0.12–1.20 %) and 0.32 % (0.12–1.28 %) in the generator scenario with the measured detector factor (0.19 % with the old one). With respect to the numbers of 18–19.11 at WP80 (D⁰ 35 or 18.6, D⁺ 18 or 31, D* 2 or 4): the inclusive yields are lower (D⁰ 9596 → 8041 with model M8, D⁺ 7599 → 6583), the model values for D⁺ and D* now use R = 0.384 % with its pT dependence instead of 0.75 %, and the detector factor is the measured one. The ranges are a spread of assumptions and models, not confidence intervals; they do not include the uncertainty of N_incl (about 20 % for D⁰ from the fit closure, 6–8 % for D⁺) or a possible difference of the classifier efficiency between diffractive and non-diffractive events. For region R1 of run 428770 (0.4687 nb⁻¹, `HLT_noalg_L1AFP_A_OR_C`) the prediction is that of 427929 times 0.635 (range 0.60–0.77): luminosity ratio 1.008, L1 AFP efficiency for an offline tag 0.684 (0.91 on side A, 0.46 on side C), ratio of the pile-up-veto efficiencies 0.955, of the ξ_cal-cut efficiencies 0.974 and of the AFP efficiencies 0.989; there is no inclusive yield in R1, because the trigger requires a proton. For R2 (27.65 nb⁻¹, jet trigger on L1 AFP) no model number is given: the efficiency of the jet requirement for diffractive events with a D is not known, and R falls with pT; the sideband anchor on the whole run gives of the order of 110 D⁰ and 60 D⁺ (× 0.5–2) with 5–9 % accidental tags. A differential measurement in 427929 is not feasible; one pT-integrated point is.
>
> **Missing parts of the colour-reconnection variation (closed on 20.11).** The parts missing on 19.11 were produced with new seeds: the diffractive sample now has 10 of 10 parts and the non-diffractive one 22 (20 planned). The envelope is unchanged: R(D⁰) = 0.384 ± 0.005 %, 0.154–1.53 % (× 0.40 to × 3.98). The crash is a segmentation fault inside the QCD colour-reconnection code of Pythia 8.317 (`ColourReconnection::getJunctionIndices`), reproducible for one seed, about one event in 10⁷; the hangs are a non-terminating loop in the junction handling of the new beam-remnant model (`JunctionSplitting::splitJunGluons`, one stack sample of one of the four seeds), four seeds in 2.6·10⁷ non-diffractive events. Both are properties of the generator configuration (`ColourReconnection:mode = 1`, `BeamRemnants:remnantMode = 1`), not of memory or of the batch nodes; the lost events are too few to bias R (0.512 % before and after for the reconnection scenario).
> <!-- src: charm_work/SD_var_2026-11-19/RAPORT_domkniecie.md §1–§4; out/pred_tables.md; out/tables.md §2; charm_work/AFP_akceptancja_2026-11-20/PLAN.md §0, §5; charm_work/CVER_2026-11-20/RAPORT.md pkt 4–5 -->

## 6. Pozostałości i zastrzeżenia

1. Na lxplus (nieusunięte): `~/private/charm_q/Q_lxplus/sd_var/closure/` (DAG, logi dagmana, ~30 kB), logi `logs/sd_var/close_*` i `sdclose.log`; EOS `charm_q/sd_var/merged_close/` (3 pliki, 0,3 MB) oraz nowe części `dv_CRqcd/part_1230101.*`, `nv_CRqcd/part_12601{01..06}.*`. Z wcześniejszego przebiegu nadal leżą trzy puste `core.*` w `sd_var/` i `lxplus932:/tmp/klyczak_sdvar`.
2. `EOS …/sd_var/merged/{dv,nv}_CRqcd.json` to nadal stan 9/10 i 16/20 (nie nadpisywałem, żeby nie ruszać wyniku etapu `sdvar`); aktualne są pliki w `merged_close/`. Ponowne uruchomienie kroku `merge` paczki scali wszystkie części (10 i 22).
3. Uzyski D⁰ w `XSEC_2026-11-19/fitbias_override.json` i w `TAGGED KRYTERIA.md` (4105 / 1530 / 681) mają jeszcze skalę κ × 1,364; tu użyłem `FITBIAS out/var_data_k1.json` (κ_MC: 4077 / 1548 / 705). BKGCAT podaje dla „r = 1 + szablon wprost” nieco inne przesunięcia (−1,1 / +0,4 / +2,8 %) — różnica ≤ 1 % uzysku, bez znaczenia dla przewidywań.
4. D*: N_incl = 885 pochodzi z selekcji z klasyfikatorem (WP80); pomiar tagowany używa D* bez klasyfikatora — porównywać ułamek (0,23 % / 0,32 %), nie liczbę.
5. Scenariusz M nie zawiera: Reggeonu, innego generatora, C_D(SD)/C_D(ND), zależności A od t poza fast-sim. Skala σ_SD × 0,64 (ALFA) i podwójne liczenie CRqcd (grupy X i L) — jak w RAPORT §4, §6.
6. Przeniesienie do R1 zakłada tę samą wydajność L1 AFP w obu runach (pomiar `AFP_akceptancja`: A 0,904 / 0,915, C 0,434 / 0,458 — zgodne) i nie poprawia na inną krotność przy μ = 0,05.
7. Otwarte decyzje: wklejenie akapitu z §5 do noty; czy nominałem w nocie ma być WP70 czy WP80; czy nadpisać `merged/` na EOS plikami z `merged_close/`.
