# XSEC_final — końcowy rachunek σ(D⁰), σ(D⁺), run 427929 (projekt B, 20.11.2026)

**Status: PRELIMINARY, POST HOC.** Definicje: `DEFINICJA.md` (zapisane przed fitami). Tabele: `tabele/final.md` (generowane przez `xsec_final.py` z `out/*.json`). Liczby: `out/final.json`.
Zastępuje: `XSEC_2026-11-19/out/tabele.md` T1 / T1b / T2 i `WYNIKI.md` §1, §2, §7, §8.3, §8.5; `CVER` §5; `SFABS_wynik` „Przekrój”; `XSEC_atlas` §C. Żaden z tych katalogów nie został zmieniony.
Tylko piki inkluzywne, bez tagu protonu. Praca lokalna (bez lxplus).

## 0. Wynik

L = 0,4648 ± 0,0103 nb⁻¹. D + anty-D, prompt + nieprompt, |η(D)| < 2,2. σ = N / (C_ekskl × SF × BR_PDG × L).

| | bin [GeV] | N | C_ekskl [%] | SF | **σ [µb] ± stat ± syst ± L** | razem | dane/Pythia A3 | dane prompt/FONLL centr. (pasmo) | dane/CMS 13 TeV | dane/ALICE 13 TeV* |
|---|---|---|---|---|---|---|---|---|---|---|
| D⁰ | 3–5 | 4059 ± 468 | 22,32 | 0,967 ± 0,033 | **1028 ± 119 ± 274 ± 23** | 29 % | 0,84 ± 0,24 | 1,50 (w paśmie) | — | 0,95 ± 0,29 |
| D⁰ | 5–8 | 1536 ± 161 | 30,01 | 1,005 ± 0,126 | **278 ± 29 ± 47 ± 6** | 20 % | 0,58 ± 0,12 (−3,6σ) | 1,13 (w paśmie) | **0,53 ± 0,12 (−3,1σ)** | 0,82 ± 0,18 (−1,0σ) |
| D⁰ | 8–30 | 700 ± 70 | 39,57 | 0,965 ± 0,097 | **100 ± 10 ± 21 ± 2** | 23 % | 0,63 ± 0,14 (−2,6σ) | 1,11 (w paśmie) | 0,75 ± 0,20 (−1,2σ) | 0,86 ± 0,22 |
| D⁺ | 3–5 | 3547 ± 276 | 16,11 | 0,949 ± 0,048 | **536 ± 42 ± 73 ± 12** | 16 % | 0,84 ± 0,13 (−1,2σ) | 1,99 (+0,5σ nad pasmem) | — | 1,09 ± 0,20 |
| D⁺ | 5–8 | 1651 ± 76 | 23,06 | 0,993 ± 0,052 | **167 ± 8 ± 28 ± 4** | 17 % | 0,66 ± 0,11 (−3,0σ) | 1,69 (+0,4σ nad pasmem) | 0,88 ± 0,21 (−0,5σ) | 1,06 ± 0,21 |
| D⁺ | 8–30 | 687 ± 59 | 30,83 | 0,990 ± 0,051 | **52,0 ± 4,5 ± 7,4 ± 1,2** | 17 % | 0,61 ± 0,10 (−3,8σ) | 1,43 (+0,3σ nad pasmem) | 0,86 ± 0,19 (−0,7σ) | 1,07 ± 0,21 |

\* sprawdzian dodany po fakcie (§5.3), zależny od kształtu rapidity Pythii.

| bin | **D⁺/D⁰ ± stat ± syst** | Pythia A3 | ALICE | LHCb | e⁺e⁻ | CMS | **D*/D⁰ ± stat ± syst** | prawda mb_nd | CMS |
|---|---|---|---|---|---|---|---|---|---|
| 3–5 | **0,52 ± 0,07 ± 0,15** | 0,521 (0,0σ) | 0,458 (+0,4σ) | 0,405 (+0,7σ) | 0,395 (+0,8σ) | — | — | — | — |
| 5–8 | **0,60 ± 0,07 ± 0,13** | 0,528 (+0,5σ) | 0,470 (+0,9σ) | (+1,3σ) | (+1,4σ) | 0,37 ± 0,07 (+1,5σ) | **0,48 ± 0,08 ± 0,11** | 0,416 (+0,5σ) | 0,39 ± 0,06 (+0,6σ) |
| 8–30 | **0,52 ± 0,07 ± 0,12** | 0,534 (−0,1σ) | 0,430 (+0,6σ) | (+0,8σ) | (+0,9σ) | 0,46 ± 0,07 (+0,3σ) | **0,61 ± 0,09 ± 0,13** | 0,437 (+1,0σ) | 0,42 ± 0,06 (+1,1σ) |

σ(D*±) bez zmian (XSEC_atlas): 134 ± 18 ± 23 ± 3 µb (5–8), 60,6 ± 7,1 ± 6,8 ± 1,3 µb (8–30); wobec CMS 0,66 ± 0,17 (−1,8σ) / 1,08 ± 0,23; wobec ATLAS 7 TeV × K_FONLL 0,76 ± 0,19 / 0,90 ± 0,18; wobec Pythii 0,67 ± 0,15 / 0,87 ± 0,14.
D⁺ wobec ATLAS 7 TeV × K_FONLL (1,84 / 2,08): 0,84 ± 0,18 (−0,8σ) / 0,82 ± 0,17 (−1,0σ).

Względem XSEC z 19.11: σ(D⁰) × 0,753 / 0,740 / 0,760; σ(D⁺) × 0,990 / 0,946 / 0,948 (mianownik C, BR bez DCS, uzysk r = 1, SF z SFABS). Względem CVER (276 / 97 µb): +0,8 % / +3,3 %.

## 1. Rozstrzygnięcie rozbieżności κ (`skrypty/kappa_fit.py`, `out/kappa_{A,B,altmin}.json`)

**Obie liczby były policzone poprawnie, ale mierzą co innego; +6 % XSEC_atlas przy 5–8 GeV to przeskok fitu między dwoma minimami, nie efekt κ.**

| | 3–5 | 5–8 | 8–30 |
|---|---|---|---|
| FITBIAS (szablon T interpolowany, κ × 1,364) — baza BKGCAT | 4105 | 1530 | 681 |
| ten sam szablon, r = 1 | 4077 | 1548 | 705 |
| szablon z drzewa H, κ × 1,364 | 4080 | 1510 | 679 |
| **szablon z drzewa H, r = 1 (nominał)** | **4059** | **1536** | **700** |
| sam efekt skali κ (1,364 → 1), histogramy WPSCAN | −0,5 % | **+1,7 %** | **+3,1 %** |
| sam efekt szablonu (wprost − interpolowany, r = 1) | −0,5 % | −0,8 % | −0,7 % |
| łącznie wobec FITBIAS (= liczba BKGCAT) | −1,1 % | +0,4 % | +2,8 % |
| XSEC_atlas: N(κ × 1)/N(κ × 1,364), własne tabele, \|η_reco\| < 2,2, start domyślny | −0,4 % | **+6,1 %** | +2,6 % |
| to samo, fit κ × 1,364 startowany z rozwiązania κ × 1 (to samo minimum) | −0,5 % | **+1,4 %** | +2,6 % |
| to samo bez cięcia η (= histogram WPSCAN co do ≤ 3 wpisów w binie) | −0,5 % | +1,8 % | +3,1 % |

- **(a) baza:** BKGCAT podaje zmianę łączną (skala + szablon) względem zapisu FITBIAS; XSEC_atlas — samą skalę κ we własnym ficie. Efekt samej skali w obu rachunkach jest ten sam: +1,4…+1,8 % (5–8), +2,6…+3,1 % (8–30), −0,5 % (3–5).
- **(d) drugie minimum:** w próbce XSEC_atlas (5–8 GeV, z cięciem η) fit M8 ma dwa minima w skali szerokości sygnału: sc ≈ 1,00 (N = 1507 przy κ × 1) i sc ≈ 0,945 (N = 1443 przy κ × 1; 1421 przy κ × 1,364), różniące się o Δ(−2 ln L) = 0,1–0,3. Start domyślny trafia przy κ × 1 w pierwsze, przy κ × 1,273 i 1,364 w drugie → pozorne „+6 %”. Skan po κ pokazuje skok między r = 1,182 a 1,273 (1499 → 1426).
- **(b) próbka, (c) szablon:** bez znaczenia — tabele XSEC_atlas bez cięcia η odtwarzają histogram WPSCAN (19 078 wobec 19 085 kandydatów) i liczby BKGCAT (1539 / 1513 wobec 1536 / 1510); szablon z tabel XSEC_atlas i z `bkgcat.npz` dają to samo (κ 1,108 wobec 1,108).
- **Skutek uboczny (POST HOC, poza DEFINICJĄ):** to samo podwójne minimum istnieje w fitcie nominalnym na histogramach WPSCAN: 5–8 GeV N = 1468 (sc 0,953, Δ(−2 ln L) = −0,06 względem nominału 1536, sc 1,001), 3–5 GeV N = 4245 (sc 1,019, +0,14 względem 4059, sc 0,971); 8–30 GeV jedno minimum.
  Nominał zostaje zgodnie z DEFINICJĄ (start domyślny, jak FITBIAS/BKGCAT); różnicę wpisałem jako osobny wiersz budżetu „drugie minimum” (4,6 / 4,4 / 0 %) — **odstępstwo od DEFINICJI w stronę ostrożną**. Przyczyna leży w modelu M8 (szablon PCHIP × skala szerokości z ograniczeniem 10 %), którego nie wolno mi było dostrajać.

**Uzyski końcowe D⁰ (427929, WP70, M8, szablon z drzewa H, r = 1):** 4059 ± 468 (bootstrap 468; 1,3 × HESSE 463), 1536 ± 161 (1,3 × HESSE; bootstrap 111), 700 ± 70 (1,3 × HESSE; bootstrap 54). p(χ²) = 0,24 / 0,22 / 0,065.
Systematyka wariantów: szablon wprost vs interpolowany 0,5 / 0,8 / 0,7 %; κ ± 5,8 % → 0,07 / 0,19 / 0,44 % (± 20 %: 0,4 / 0,6 / 1,5 %); model (M4, M9) 7,5 / 1,5 / 12,8 %; domknięcie mb_nd −22,0 ± 9,7 / −4,6 ± 6,5 / +7,2 ± 7,8 %.

## 2. C ekskluzywne

Mianowniki `XSEC_atlas/out/truth_ngen.json` i CVER (tab. §2, „ekskl. bez γ + z γ”) **zgadzają się co do sztuki** (asercja w skrypcie): D⁰ 16 366 / 6394 / 2143, D⁺ 20 141 / 7923 / 2734.
C_ekskl = 22,32 / 30,01 / 39,57 % (D⁰), 16,11 / 23,06 / 30,83 % (D⁺); C_ekskl/C_XSEC = 1,363 / 1,362 / 1,357 i 1,065.
BR: D⁰ 3,936 % (bez DCS — mianownik bez DCS, DCS w MC siedzi w składowej pikującej; CVER używał 3,951 %, stąd +0,4 % względem jego 276 µb), D⁺ 9,31 %.
Uwaga: C_ekskl D⁰ w XSEC_atlas (21,37 / 28,57 / 37,61 %) jest niższe, bo tam licznik ma cięcie |η_reco| < 2,2; σ z tamtej wersji: 288 / 96 µb (WP70), 286 / 107 µb (cięcia) — zgodne z 278 / 100 µb.

## 3. SF

Nominał wg DEFINICJI §4. D⁺ w wariantach (σ [µb] ± stat ± syst ± L):

| bin | nominał (SFABS) | SF = 1 ± 0,10 | SFABS ± 0,08 |
|---|---|---|---|
| 3–5 | 536 ± 42 ± 73 ± 12 | 509 ± 40 ± 82 ± 11 | 536 ± 42 ± 82 ± 12 |
| 5–8 | 166,7 ± 7,7 ± 27,6 ± 3,7 | 165,5 ± 7,6 ± 30,8 ± 3,7 | 166,7 ± 7,7 ± 29,4 ± 3,7 |
| 8–30 | 52,0 ± 4,5 ± 7,4 ± 1,2 | 51,5 ± 4,5 ± 8,6 ± 1,2 | 52,0 ± 4,5 ± 8,1 ± 1,2 |

Sprawdzian cięciowy XSEC_atlas (1/R, 5–30 GeV, niezależny od SFABS, nie wchodzi do nominału): **D⁺ 0,97 ± 0,05 ± 0,05**, D⁰ 0,95 ± 0,06 ± 0,08 — zgodne z SF nominalnymi (0,99 ± 0,05; 1,005 ± 0,126 i 0,965 ± 0,097). Wybór wariantu SF zmienia σ(D⁺) o ≤ 5 % przy 3–5 GeV i ≤ 1 % wyżej.

## 4. Budżet (pełna tabela: `tabele/final.md` F2)

| [%] | D⁰ 3–5 | D⁰ 5–8 | D⁰ 8–30 | D⁺ 3–5 | D⁺ 5–8 | D⁺ 8–30 |
|---|---|---|---|---|---|---|
| stat | 11,5 | 10,5 | 10,1 | 7,8 | 4,6 | 8,7 |
| fit: model / obciążenie / drugie minimum | 7,5 / 24,0 / 4,6 | 1,5 / 8,0 / 4,4 | 12,8 / 10,6 / 0 | 3,9 / 7,8 / — | 2,9 / 12,4 / — | 7,7 / 5,6 / — |
| pikująca: szablon / κ | 0,5 / 0,1 | 0,8 / 0,2 | 0,7 / 0,4 | — | — | — |
| SF | 3,4 | 12,5 | 10,1 | 5,0 | 5,2 | 5,2 |
| śledzenie / d₀ + materiał / wierzchołek / generator (nowe) | 5,2 / 1,5 / — / 2,3 | 5,2 / 1,5 / — / 2,1 | 5,2 / 1,5 / — / 2,2 | 7,8 / 1,5 / 3,0 / 2,5 | 7,8 / 1,5 / 3,0 / 2,3 | 7,8 / 1,5 / 3,0 / 2,2 |
| stat MC / BR / nieprompt / odbicie / zdarzenia | 1,5 / 0,8 / 3,1 / 0,1 / 0,2 | 1,9 / 0,8 / 1,7 / 0,1 / 0,2 | 2,7 / 0,8 / 0,4 / 0,1 / 0,2 | 1,6 / 1,5 / 1,6 / — / 0,2 | 2,1 / 1,5 / 0,6 / — / 0,2 | 2,9 / 1,5 / 0,6 / — / 0,2 |
| **razem syst** | **26,7** | **16,8** | **20,5** | **13,7** | **16,5** | **14,3** |
| L | 2,2 | 2,2 | 2,2 | 2,2 | 2,2 | 2,2 |

Podwójne liczenie — rozstrzygnięcia: (1) wiersz dE/dx z XSEC_atlas pominięty (SF mierzony na danych go zawiera); (2) „κ bez skali BR” usunięte z wiersza „model fitu” (to teraz nominał), w zamian szablon wprost/interpolowany i κ ± 5,8 % z własnego refitu (zastępuje wiersz XSEC_atlas 0,9 / 0,4 %, który był liczony względem skoku minimum);
(3) śledzenie nie jest w SF (SF = klasyfikator przy danej rekonstrukcji) — dodane w całości; (4) d₀/materiał 1,5 % zostawione w całości, choć część klasyfikatorowa jest w SF (zachowawczo).
Błąd stat. D⁰ jest większy niż w XSEC (8 → 10 % przy ≥ 5 GeV), bo HESSE × 1,3 zamiast HESSE. D⁺ zostaje przy max(HESSE, bootstrap); z × 1,3 byłoby 9,4 / 5,9 / 8,7 % (informacyjnie).

**Korelacje w D⁺/D⁰:** skracają się L i wybór zdarzeń; śledzenie → 1 tor (2,6 %), d₀/materiał → 1,0 %, generator → 1,5 %, wierzchołek 3,0 % zostaje; SF: 3–5 GeV SF⁺/SF⁰ = 0,981 ± 0,042 (wspólny SF_ref skrócony), ≥ 5 GeV prawie nieskorelowane (waga komórki 427929 w SF_ref 0,06–0,07) → 12,9 / 11,3 %.
Budżet stosunku: stat 13,9 / 11,4 / 13,3 %; obciążenie fitu 25,3 / 14,8 / 12,0 %; model 8,4 / 3,3 / 14,9 %; razem 31 / 24 / 27 %. Stosunek ogranicza fit D⁰ (3–5) i SF⁰ z tagu D* (≥ 5 GeV).

## 5. Porównania

### 5.1 CMS 13 TeV (arXiv:2107.01476) — lista kontrolna definicji (praca otwarta 20.11: arXiv HTML v2; HEPData niedostępne — 403)

| czynnik | CMS | my | liczba przeliczenia | źródło |
|---|---|---|---|---|
| cząstka + antycząstka | **suma** (podpisy tabel „D⁰ + D̄⁰”, „D⁺ + D⁻”; równ. 2: „including the charge-conjugate final state”) | suma | 1 | CMS §6, tab. dσ/dpT; **nie ma czynnika 2** — potwierdza to też rząd wielkości względem ALICE (§5.3) |
| zmienna i zakres | \|η\| < 2,1 (**η, nie y**: „pseudorapidity is used instead of rapidity to facilitate the comparison with the ATLAS measurement”) | \|η\| < 2,2 | × 0,960 / 0,961 (D⁰), 0,959 / 0,964 (D⁺) z prawdy mb_nd; płasko 0,955 | CMS §1; `truthD_mb.parquet` |
| prompt | tylko prompt; nieprompt odjęty **z symulacji** (Pythia 8 CUETP8M1), 5–17 % zależnie od binu i mezonu; wartości per bin tylko na rys. 3 | prompt + nieprompt | × (1 − f_np) = 0,893 / 0,874 (D⁰), 0,908 / 0,889 (D⁺); f_np z mb_nd, niepewność ± 50 % f_np (5–7 % na wyniku) | CMS §4.4; XSEC `mc.f_np_all`. FONLL dałby f_np ≈ 12 / 16 % → stosunek niżej o 1,5 / 4 % |
| tabela | **dσ/dpT [µb/GeV]** (równ. 2 dzieli przez ΔpT) | σ w binie | Σ wartość × Δ: 5–6, 6–7, 7–8 (Δ = 1) → 432 ± 22 ± 39; 8–12 … 24–40 → 107,8 ± 2,9 ± 10,3 | CMS tab. dσ/dpT; suma XSEC odtworzona co do cyfry |
| zakres pT | 8–40 | 8–30 | − 0,5 µb (0,5 %; 30–40 GeV przy widmie ~pT⁻⁵) | własne |
| energia | 13 TeV | 13,6 TeV | ÷ 1,036 / 1,0455 (FONLL prompt, \|η\| < 2,2) | XSEC2 `xsec2/fonll` |
| świetlność CMS | 2,5 %, w kolumnie syst | osobno 2,2 % | — | CMS tab. syst. |
| BR | PDG 2018 (wartości nie cytowane; błąd 0,8 % D⁰) | PDG 2026 3,936 % | zmiana PDG 2018 → 2026 < 1 % (niesprawdzone w źródle) | CMS §5 |
| łącznie | | | **× 0,827 / 0,803 (D⁰), 0,840 / 0,820 (D⁺)** | |

**Wynik sprawdzenia: nie ma błędu definicji.** Każdy czynnik jest bliski 1 i udokumentowany; czynnik 2 wykluczony, frakcja nieprompt zmienia stosunek o kilka procent. D⁰ 5–8 GeV: 230 ± 48 µb w definicji CMS wobec 432 ± 45 → **0,53 ± 0,12 (−3,1σ)** — napięcie zostaje (było 0,53 ± 0,11).
D⁰ 8–30: 0,75 ± 0,20 (−1,2σ); D⁺: 0,88 ± 0,21, 0,86 ± 0,19; D*: 0,66 ± 0,17 (−1,8σ), 1,08 ± 0,23.
Nie zweryfikowane wprost: wartości z HEPData (odczyt z tabeli HTML; kolejność kolumn sprawdzona na błędach stat. względem tabeli uzysków), frakcje nieprompt CMS per bin (rys. 3).

### 5.2 FONLL, Pythia A3, ATLAS 7 TeV

- **FONLL 13,6 TeV (prompt, 2 f σ; pasmo skal ⊕ mas ⊕ PDF):** D⁰ 221 (138–344) i 79 (60–107) µb wobec danych prompt 248 ± 52 i 87 ± 21 → **w paśmie**, 1,13 / 1,11 wartości centralnej; 3–5 GeV: 929 ± 275 wobec 619 (266–1115). D⁺: 151 ± 27 wobec 89 (55–140), 46 ± 8 wobec 32 (24–44) → tuż nad górną krawędzią (+0,4 / +0,3σ); 3–5: 491 ± 81 wobec 247 (104–447), +0,5σ.
- **Pythia 8 A3 (wszystkie D gatunku w obszarze, × σ_ND/N_zd ÷ f_ND — niezależne od `decayType`):** dane/Pythia 0,84 / 0,58 / 0,63 (D⁰), 0,84 / 0,66 / 0,61 (D⁺), 0,67 / 0,87 (D*); istotność −2,6…−3,8σ przy ≥ 5 GeV (błąd tylko danych). Oba gatunki jednakowo — Pythia A3 zawyża o ~1,6 przy ≥ 5 GeV (CMS: „dane poniżej Pythii, powyżej FONLL”).
- **ATLAS 7 TeV × K_FONLL:** D⁺ 0,84 ± 0,18 / 0,82 ± 0,17; D* 0,76 ± 0,19 / 0,90 ± 0,18 — poniżej, w granicach 1,2σ.

### 5.3 Sprawdzian wobec ALICE 13 TeV (dodany po fakcie; nie było go w DEFINICJI)

Dane prompt przeliczone do |y| < 0,5 kształtem rapidity Pythii (f_y = 0,246 wobec 0,227 przy płaskim η), ÷ 2 (średnia cząstka/antycząstka), ÷ K_E: **D⁰ 0,95 ± 0,29 / 0,82 ± 0,18 / 0,86 ± 0,22; D⁺ 1,09 ± 0,20 / 1,06 ± 0,21 / 1,07 ± 0,21** — zgodne.
To samo przeliczenie zastosowane do CMS daje CMS/ALICE = **1,54 ± 0,20 (D⁰ 5–8)**, 1,15 ± 0,15 (D⁰ 8–40), 1,20 ± 0,21 i 1,25 ± 0,19 (D⁺).
Odczyt: napięcie D⁰ 5–8 GeV jest w dużej części napięciem CMS–ALICE; nasze D⁰ leży między nimi, bliżej ALICE. Zastrzeżenia: wartości ALICE to sumy binów z `XSEC WYNIKI.md` §5 (niezweryfikowane tutaj), błąd ALICE przyjęty 8 %, kształt y z Pythii.

### 5.4 Pozostałe napięcia

1. **D⁰ 5–8 GeV wobec CMS: −3,1σ** (nie błąd definicji; CMS w tym binie 1,5 × ALICE).
2. **Dane/Pythia A3 ≈ 0,6 przy ≥ 5 GeV dla D⁰, D⁺ i D*** (−2,6…−3,8σ) — wspólne dla gatunków, więc raczej generator niż efekt jednego kanału.
3. **Stosunki wewnętrzne ciągną w jedną stronę:** D⁺/D⁰ = 0,60 przy 5–8 GeV (+1,3σ wobec LHCb, +1,5σ wobec CMS), D*/D⁰ = 0,48 / 0,61 (+0,5 / +1,0σ wobec prawdy mb_nd). Każde z osobna nieistotne; razem sugerują σ(D⁰) niższe o 15–30 % względem D⁺ i D*. Selekcja cięciowa bez klasyfikatora daje ten sam poziom D⁰ (XSEC_atlas), więc kandydaci to model fitu M8 (podwójne minimum, obciążenie 8–24 %) i preselekcja D⁰, nie klasyfikator.
4. **D⁰ 3–5 GeV:** domknięcie M8 na MC −22 ± 10 % (fluktuacja tła wg BKGCAT) — w budżecie 24 %, bez poprawki; z poprawką σ = 1317 µb.
5. D⁰ 8–30: p(χ²) = 0,065, model M4/M9 daje +11…+13 %.

## 6. Zastrzeżenia

- Wszystko post hoc; uzyski i poziom napięcia z CMS były znane przed zapisaniem DEFINICJI.
- Odstępstwa od DEFINICJI: (i) wiersz „drugie minimum” (§1) dodany po fitach; (ii) sprawdzian ALICE (§5.3) dodany po fakcie; (iii) K dla D⁺ (ATLAS 7 TeV) z samego prompt — zgodnie z DEFINICJĄ, ale niepewność K tylko z min/min.
- Wiersze śledzenia, d₀/materiału, wierzchołka i krotności są z literatury ATLAS, nie z pomiaru na tych danych.
- SF D⁺ i D⁰ 3–5 GeV opiera się w 81–97 % na przeniesieniu z D⁰ ≥ 5 GeV (SFABS; warunkowe); SF_ref zdominowany przez run 428770.
- Kalibracja świetlności 2,2 % to założenie (zachowawczo 7 %).
- σ(D*) i jej budżet przepisane z XSEC_atlas bez przeliczania; korelacja statystyczna D* z próbką D⁰ pominięta (błąd D*/D⁰ zawyżony).
- Stosunki CMS (D⁺/D⁰, D*/D⁰) to moja arytmetyka z tabeli (syst: suma kwadratowa minus wspólna świetlność) — CMS takich stosunków nie podaje.
- Mianownik `mode == 1` w `GEN_mcstat`/`gen_norm` (projekt A) — nie ruszane (uwaga CVER §5 pkt 6 pozostaje otwarta).

## 7. Pliki

`DEFINICJA.md`, `RAPORT.md`, `xsec_final.py`, `skrypty/kappa_fit.py` (A: histogramy WPSCAN; B: tabele XSEC_atlas; C: minima), `skrypty/xsec_orig_kopia_2026-11-19.py` (kopia oryginału, nieużywana w rachunku),
`out/{final.json, kappa_A.json, kappa_B.json, kappa_altmin.json, kappa_*.log}`, `tabele/final.md`, `PODMIANY.md`, `NOTA_6.13_EN.md`.
Baza wiedzy: `charm-bkg-notes/wyniki/Przekrój inkluzywny D0 i D+ 427929 końcowy - 278 i 100 µb, D+ 167 i 52 µb, rozbieżność kappa to drugie minimum fitu (XSEC_final).md` + wiersz w `00_INDEX.md`.
Pozostałości: brak poza tym katalogiem (PDF pracy CMS pobrany do katalogu tymczasowego, nieużyty).
