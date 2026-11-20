# AFP-MIX v2: tło przypadkowego tagu AFP z mieszania protonów — opis metody

Kod: `charmbkg/afp_mix.py` (metoda, bramki, procedura), `charmbkg/afp_mix_toy.py` (toy, scenariusze naruszające założenia,
tabela mocy), testy `tests/test_afp_mix.py`; sterowanie na plikach parquet: `analysis/afp_mix/`. Dokument opisuje algorytm. **Wyniki na danych nie są częścią repozytorium.**

## 1. Model i konwencje

Zdarzenie z tagiem = prawdziwa dyfrakcja ⊕ niezależny proton (pile-up bez wierzchołka, halo). Obiekt tła: kandydaci D biorcy
+ zawartość AFP dawcy z tej samej warstwy. Konwencje charmpp (gałąź minimum-bias): tag = dokładnie jeden proton po stronie,
w oknie 0,030 ≤ ξ ≤ 0,100 (`tag_def='nside'`; toy: `'inwindow'`); region A∨C; cięcie ξ_cal ≤ ξ1_AFP + 0,02 na każdej otagowanej stronie.
N_tag = N_diff·(1 − f̄) + N_acc; f_acc = N_acc/N_tag.

**Zakres ważności (W2):** tylko zdarzenia wybrane trygerami nie-AFP. Zdarzenia wyłącznie z łańcuchów L1AFP są wykluczone
u biorców i dawców i metoda ich nie opisuje.

## 2. Estymator A (nominał)

- **Pula:** dawcy z wetem ξ_cal > ξ_hi + δ + margines = 0,14 po stronie (własny proton nie może przejść cięcia).
  Region A∨C = suma mieszań z **pul jednostronnych** (tag przypadkowy po obu stronach naraz, O(p²), pominięty). Pula dwustronna
  została odrzucona (inny skład i aktywność niż u biorców).
- **Warstwy:** (run, LB) × klasa nSelBL (10, 15, 20, 25, 30, 40, 50, 80) × „z D/bez D”; hierarchia (run, LB, klasa) → (run, μ, klasa)
  → (run, klasa) → (run), pełna pula na każdym poziomie (`n_min_pool` = 50). Częstość protonu przypadkowego zależy od aktywności
  i od obecności D (W1); pula ma inną aktywność niż biorcy z małym ξ_cal. Warstwowanie usuwa to z konstrukcji, bez zawężania puli
  (pełna statystyka, W3). Biorcy (z D) dostają więc dawców z D tej samej klasy aktywności.
- **Overlay „union” na jawnych parach** (biorca, otagowany dawca): para liczy się, gdy region przechodzi z dawcą, a bez niego nie.
  Zastępuje iteracyjną poprawkę „własnego okna” z v1. Biorca wykluczony z własnej puli.
- **Błąd stat.:** bootstrap Poissona po otagowanych dawcach, wagi wspólne dla wszystkich wariantów → błąd różnicy wariantów.
- **Ekstrapolacja pula → obszar użycia:** profil R = obs/pred w przedziałach ξ_cal puli (0,14; 0,2; 0,3; 0,5; 1), dopasowanie liniowe
  w ln ξ_cal, κ = R w ⟨ln ξ_cal⟩ obszaru użycia (ważone p_acc). **Bez korekty**; systematyka √((κ−1)² + σ_κ²). Uzasadnienie:
  model liniowy w ln ξ_cal jest najprostszym, który dane potrafią ograniczyć (4 przedziały, ramię ≈ 1,5 jednostki ln); korekta
  wymagałaby wiary w model poza zakresem, więc pełne odchylenie wchodzi do błędu.

## 3. Drugie estymatory (reguła „połowa różnicy”, CMS-TOTEM 2002.12146)

- **B(i), styl ALFA (1911.00453):** absolutne P(tag) z próbki quasi-minbias o pełnej aktywności (łańcuchy L1MBTS, weto ξ_cal po obu
  stronach, nSelBL ≥ 20) w kwantylowych przedziałach μ, bez klas aktywności i bez klasy D. Próbka quasi-minbias
  zastępuje tryger losowy, gdy ten jest niedostępny; gdy łańcuchy MBTS stanowią większość próbki, brane są tylko zdarzenia bez D (rozłączne z pulą A).
- **B(ii), styl CMS-TOTEM („przesuwane weto”):** N_B2 = N_A × N_obs/N_pred w obszarze kontrolnym: tag po stronie, ξ_cal ≤ 0,12,
  Δ = ξ_cal − ξ_AFP > δ + 0,025 (czysto przypadkowe z konstrukcji, najbliżej obszaru użycia).
- Do błędu wchodzi ½·max|A − B|. Zgodność metod = główne bramki G7, G8.

## 4. Bramki

Werdykt PASS = żadna bramka werdyktu nie ma FAIL (kryteria: |z| < 2,5 albo p ≥ 0,01; N/A przy < 5 oczekiwanych).
G4 z v1 usunięta (D-test jednostronny bez mocy); G3 i G6 z v1 → I3, I6, tylko informacyjne. G1 i G5: iloraz wiarygodności
Poissona (odporny na zera, łączenie przedziałów do ≥ 5 oczekiwanych), p dla G5 z pseudoeksperymentów z replikami bootstrapu.

| bramka | rola | co testuje | czego NIE testuje |
|---|---|---|---|
| G0_inputs | werdykt | pokrycie warstw (≥ 99 % biorców w komórce zachowującej klasę aktywności: LB×klasa → μ×klasa → run×klasa; 0 bez puli), brak zdarzeń wyłącznie z trygera AFP | poprawności fizycznej; zakresu poza trygerami nie-AFP (W2) — tam metoda nie obowiązuje |
| G1_pool_homogeneity | werdykt | czy po warstwowaniu częstość tagu w puli nie zależy od ξ_cal dawcy POWYŻEJ progu (Poisson LR obs vs pred w przedziałach ξ_cal) | obszaru poniżej progu (obszar użycia); zależności, które znoszą się w przedziale |
| G2_disjoint_closure | werdykt | domknięcie na zbiorach ROZŁĄCZNYCH w obszarze anty-ξ_cal: tagi zdarzeń puli z najliczniejszej klasy trygera vs przewidywanie z dawców pozostałych klas (te same warstwy LB × aktywność × D) | obszaru użycia (ξ_cal małe); zależności od ξ_cal; przy jednej klasie trygera N/A (v1: G2 była w ~90 % tożsamością) |
| I2_hasD_vs_noD | informacyjna | (informacyjna) skład puli W1: tagi zdarzeń z D w obszarze anty vs przewidywanie z dawców BEZ D, per klasa aktywności — uzasadnia klasę „z D/bez D” w warstwach | niczego o estymacie nominalnej (dawcy bez D nie są w niej używani) — nie wchodzi do werdyktu |
| G5_cr_shape | werdykt | kształt (Δ, ξ_AFP) tagów w obszarze kontrolnym Δ > δ + 0,025 u biorców SPOZA puli (ξ_cal ≤ 0,12) vs mieszanie; Poisson LR, p z pseudoeksperymentów z bootstrapem dawców | normalizacji (to G7); obszaru Δ ≤ δ; przy < ~50 zdarzeniach nie ma mocy |
| G7_methods_A_vs_Bii | werdykt | GŁÓWNA: częstość przypadkowych najbliżej obszaru użycia — N_obs/N_pred w obszarze kontrolnym (ξ_cal 0,075–0,12, Δ > 0,045; „przesuwane weto”) = stosunek metod B(ii)/A | obszaru ξ_cal < 0,075 (tam prawdziwe tagi są nieodróżnialne); przeciek własnych protonów przez rozdzielczość ξ_cal zawyża R (skan marginesu) |
| G8_methods_A_vs_Bi | werdykt | GŁÓWNA: zgodność z drugą metodą — absolutne P(tag) z próbki quasi-minbias o pełnej aktywności (weto ξ_cal po obu stronach, nSelBL ≥ 20) w przedziałach μ, bez warstw aktywności (styl ALFA: pełna niezależność) | obszaru użycia (obie pule mają ξ_cal > próg); przy tej samej próbce trygerowej pule częściowo się pokrywają |
| G9_opposite_side | werdykt | obszar „dwuprotonowy”: zdarzenia z prawdziwym tagiem (po cięciu) po stronie przeciwnej i ξ_cal > próg po stronie badanej — częstość tagu po stronie badanej vs mieszanie (czy zdarzenia dyfrakcyjne mają tę samą częstość przypadkowych) | zależności od ξ_cal po stronie badanej; CD z protonem poza oknem nie szkodzi, ale CD z rozdzielczością ξ_cal zawyża obs |
| I3_mu_slope | informacyjna | (informacyjna) nachylenie f_tag(μ) dane vs mieszanie | niczego przy σ_b(dane) ≫ b_mix — nie wchodzi do werdyktu |
| I6_eta_asymmetry | informacyjna | (informacyjna) f_acc z asymetrii η* kandydatów przy założonym A_diff | A_diff z Pythii nie pasuje do kandydatów inkluzywnych — nie wchodzi do werdyktu |

Rozważone i odrzucone jako „wolne od prawdziwych tagów”: sidebandy masy (kombinatoryka w zdarzeniach dyfrakcyjnych niesie prawdziwy
proton — tylko informacyjnie pik/SB), pasmo tuż pod progiem puli z Δ ≤ δ (prawdziwe tagi). Wolne z konstrukcji są: Δ > δ + margines
(G5, G7) i strona przeciwna do prawdziwego tagu przy ξ_cal > próg (G9).

## 5. Systematyka

σ(N_acc)/N_acc = stat ⊕ pula ⊕ ekstrapolacja ⊕ ½Δ metod. „Pula” = max po wariantach (margines 0 / 0,05 / 0,10, LB ↔ μ, zgrubne klasy
aktywności, + klasa trygera) z max(|Δ| − σ_Δ, 0). Ekstrapolacja i ½Δ(B(ii)) dotyczą tego samego pytania — suma kwadratowa jest
konserwatywna. Warianty „dawcy z D i bez D” oraz „tylko LB” (≈ v1) są raportowane informacyjnie.

## 6. Toy: scenariusze naruszające założenia i moc bramek

`python -m charmbkg.afp_mix_toy --power --seeds 6` (2 runy × 60 LB × 2500 zdarzeń; ułamek seedów z FAIL per bramka; pokrycie = ułamek
seedów z |A − prawda| < 2σ_tot deklarowanego).

| scenariusz | A/prawda | tylko LB (≈ v1)/prawda | B(i)/prawda | B(ii)/prawda | ⟨σ_tot⟩ dekl. | pokrycie 2σ | G0 | G1 | G2 | G5 | G7 | G8 | G9 | werdykt FAIL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| nominal | 0.95 ± 0.04 | 0.94 | 0.95 | 0.85 | 0.20 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| rate_vs_nsel | 0.99 ± 0.05 | 0.92 | 0.89 | 1.12 | 0.21 | 1.00 | 0.00 | 0.00 | 0.17 | 0.00 | 0.00 | 0.00 | 0.00 | 0.17 |
| rate_vs_trig | 1.07 ± 0.03 | 1.11 | 0.90 | 1.18 | 0.19 | 1.00 | 0.00 | 0.00 | 1.00 | 0.00 | 0.00 | 0.33 | 0.00 | 1.00 |
| pool_noD_x2 | 0.97 ± 0.04 | 1.59 | 1.88 | 0.90 | 0.51 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | 0.00 | 1.00 |
| xical_step_0.12_x1.5 | 0.71 ± 0.02 | 0.70 | 0.70 | 1.01 | 0.29 | 0.83 | 0.00 | 0.00 | 0.00 | 0.00 | 0.17 | 0.00 | 0.00 | 0.17 |
| xical_step_0.12_x2 | 0.49 ± 0.01 | 0.50 | 0.49 | 0.92 | 0.49 | 0.33 | 0.00 | 0.00 | 0.17 | 0.00 | 1.00 | 0.00 | 0.00 | 1.00 |
| xical_step_0.07_x2 | 0.64 ± 0.02 | 0.61 | 0.60 | 0.65 | 0.22 | 0.33 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| xical_pow_0.3 | 0.55 ± 0.03 | 0.52 | 0.53 | 0.76 | 0.66 | 1.00 | 0.00 | 0.83 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.83 |
| halo_xi_tilt_low | 0.69 ± 0.02 | 0.70 | 0.71 | 1.08 | 0.34 | 1.00 | 0.00 | 0.00 | 0.17 | 0.00 | 0.67 | 0.00 | 0.00 | 0.67 |
| xical_smear_0.5 | 1.12 ± 0.08 | 1.15 | 1.19 | 1.87 | 0.44 | 1.00 | 0.00 | 0.00 | 0.00 | 0.17 | 0.83 | 0.00 | 0.00 | 0.83 |
| independence_broken | 0.94 ± 0.05 | 0.70 | 0.59 | 0.99 | 0.26 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | 0.00 | 1.00 |
| afp_trigger | 1.02 ± 0.02 | 0.68 | 0.54 | 1.12 | 0.27 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | 0.00 | 1.00 |

Wnioski:
- zależność od aktywności i od obecności D (W1: `rate_vs_nsel`, `pool_noD_x2`, `independence_broken`, `afp_trigger`) — A domyka
  z konstrukcji (v1-like: 0,7–1,6); B(i) odstaje i G8 daje FAIL, czyli werdykt ostrzega, choć A jest poprawne;
- zależność od ξ_cal powyżej progu (`xical_pow`) → G1; skok tuż pod progiem (0,12) → G7 dopiero przy ×2 (×1,5: moc 17 %);
- **martwe pole: skok częstości poniżej ξ_cal ≈ 0,075 (`xical_step_0.07_x2`) — bias −36 %, żadna bramka nie wykrywa, pokrycie 33 %.**
  Tam prawdziwe tagi są nieodróżnialne od przypadkowych; ograniczenie nieusuwalne na tych danych;
- nominał: 0 FAIL na 6 seedach tu; w 40 seedach jednorunowych 1 fałszywy FAIL (G1); w 26 biegach formatu parquet jeden G2 z z = 3,1.

## 7. Ograniczenia

- Zakres: tylko trygery nie-AFP (W2). Dla zdarzeń z L1AFP potrzebna osobna metoda.
- Obszar ξ_cal < 0,075 nietestowalny (martwe pole z §6); ekstrapolacja liniowa w ln ξ_cal jest założeniem, nie pomiarem.
- Bez trygera losowego B(i) nie jest w pełni niezależne (ta sama próbka trygerowa, weto ξ_cal jak w A).
- Kandydaci bez score klasyfikatora (luźna selekcja): f_acc dla próbki po cięciu na score trzeba liczyć osobno.
- I6 (η*): A_diff z generatora nie stosuje się do kandydatów inkluzywnych — bramka tylko informacyjna.

## 8. Literatura

- CMS-TOTEM, arXiv:2002.12146 (dijety SD 8 TeV) — dwie metody tła nakładek: protony z zero-bias oraz „metoda HERA” (ξ_TOTEM z przypadków
  z ξ_CMS > 0,12, normalizacja w obszarze zdominowanym przez tło); różnica metod jako systematyka. Wzór dla puli A, B(ii) i reguły połowy różnicy.
- ATLAS ALFA, arXiv:1911.00453 (inkluzywna SD 8 TeV) — „overlay” z próbki kontrolnej o pełnej aktywności, absolutne prawdopodobieństwo × N. Wzór dla B(i).
- ATLAS AFP, arXiv:2009.14537 (ℓℓ + p) — losowe parowanie ξ_AFP z próbki kontrolnej, normalizacja na zdarzeniach niedopasowanych; zmiana definicji puli jako systematyka.
- ATLAS AFP, arXiv:2304.10953 (ALP → γγ) — podmiana danych AFP z innych zdarzeń w przedziałach chwilowej świetlności (warstwy LB/μ); test zanieczyszczenia puli sygnałem.
