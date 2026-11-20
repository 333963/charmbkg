# TEMPLATE — szablon tła przypadkowego tagu AFP (projekt B, AFP-MIX v2, WSTĘPNE)

Stan: 17.11.2026. Wynik i wnioski: `RAPORT.md`. Tutaj: jak odtworzyć i schemat plików.

Metoda bez zmian: `~/Projects/charm/charm-bkg/charmbkg/afp_mix.py` v2 (estymator A, `mix_region`), systematyki względne z
`../runs/lxplus/out_v2/result_v2.json`. Selekcja zdarzeń i kandydatów w `common.py` to kopia logiki `../run_data.py` v2
(N_tag i N_acc zgadzają się z `out_v2` co do zdarzenia; kontrola w `make_template.py`).

## Odtworzenie

    cd ~/Documents/charm/charm_work/AFP_mix/TEMPLATE
    PY=~/Documents/charm/.venv_h5/bin/python
    $PY make_template.py     # ~3 min, ~10 GB RAM → template.json, template.csv, template_shapes.csv
    $PY jump_tests.py        # ~4 min             → jump_tests.json   (ZADANIE 2, testy (i)–(iii))
    $PY make_report.py       # sekundy            → jump_summary.json, RAPORT.md

Wątki BLAS ograniczone do 4 (zmienne `*_NUM_THREADS`). `explore1.py` — jednorazowy przegląd wejścia (rozkład ξ protonów, zakres pT).

## Zakres ważności

- Runy 427929 i 428770 osobno (nie łączyć). h5 v2.0, kandydaci **bez score** (luźna selekcja slim: D⁺ m 1,72–2,02; D⁰ m 1,75–2,00,
  Lxy > 0,25; D* Δm 0,139–0,160; pT ≥ 3 GeV). Po selekcji score wszystko trzeba przeliczyć.
- Tylko zdarzenia z trygerów nie-AFP: GRL ∧ ≥ 1 zapisany łańcuch bez L1AFP ∧ nPileupVertices == 0. Dla zdarzeń z trygera L1AFP
  (w 428770 to 95,6 % tagów) szablon NIE obowiązuje.
- Tag: dokładnie jeden proton po stronie, 0,030 ≤ ξ_AFP ≤ 0,100; cięcie ξ_cal ≤ ξ_AFP + 0,02 na każdej otagowanej stronie; region A, C albo A∨C.
- Przedziały pT: 3–5, 5–8, 8–30, 3–30 GeV. Przedziału 2–3 GeV nie ma w danych. `3-150_v2check` to wiersz kontrolny (porównanie z `RAPORT_DANE_v2.md`).
- To są liczności tagów i frakcje tła. Plik nie zawiera uzysków sygnału, fitów masy, przekrojów ani histogramów masy próbki tagowanej.

## `template.csv` — wiersz = run × region × gatunek × przedział pT

Liczności bezwymiarowe (kandydaci), masy w GeV, pT w GeV.

| pole | znaczenie |
|---|---|
| `run` | numer runu |
| `region` | `A`, `C`, `AorC` (A∨C; zdarzenia z tagiem po obu stronach zostają) |
| `species` | `Dplus`, `D0`, `Dstar` |
| `pt_bin`, `pt_lo`, `pt_hi` | przedział pT kandydata [lo, hi) |
| `mass_var`, `win_lo`, `win_hi` | zmienna okna (`Dmeson_mass` albo `Dmeson_delta_mass`) i okno piku [lo, hi] |
| `n_cand_all_win` | wszyscy kandydaci w oknie (bez wymogu tagu) |
| `n_tag` | kandydaci w oknie w zdarzeniach z tagiem w regionie po cięciu ξ_cal (sygnał + kombinatoryka, bez fitu) |
| `n_acc` | przewidywana liczba kandydatów w oknie z tagiem przypadkowym = Σ p_acc po kandydatach |
| `n_acc_err_stat` | bootstrap Poissona otagowanych dawców (B = 500) |
| `n_acc_err_pool`, `_extrap`, `_half`, `_syst` | `n_acc` × błąd względny runu z v2: warianty puli, ekstrapolacja w ξ_cal, ½ różnicy metod A–B; `_syst` = suma kwadratowa trzech |
| `f_acc` | `n_acc / n_tag` |
| `f_acc_err_stat` | `f_acc_err_stat_donors` ⊕ `f_acc_err_stat_ntag` (bootstrap dawców ⊕ Poisson N_tag z wagami = liczba kandydatów w zdarzeniu) |
| `f_acc_err_pool`, `_extrap`, `_half`, `_syst`, `_tot` | jak dla `n_acc`; `_tot` = stat ⊕ syst |
| `mean_p_acc_all_win` | `n_acc / n_cand_all_win` — średnie prawdopodobieństwo tagu przypadkowego na kandydata nietagowanego |
| `shape_chi2`, `shape_ndf`, `shape_p` | test kształtu masy: składowa przypadkowa vs wszyscy kandydaci bez wymogu tagu; kowariancja bootstrap dawców ⊕ próbkowanie biorców bez podwójnego liczenia |
| `shape_chi2_cons`, `shape_p_cons` | wariant zachowawczy (bootstrap + pełne próbkowanie; wariancja zawyżona) |
| `shape_chi2_samp`, `shape_p_samp` | tylko próbkowanie biorców przy zrealizowanych wagach (H0 ścisłe: masa niezależna od wag) |
| `shape_chi2_null_perm` | kalibracja `_samp`: to samo po przetasowaniu p_acc między biorcami (oczekiwane ≈ ndf) |
| `shape_n_acc_fullrange`, `shape_n_all_fullrange` | sumy histogramów w pełnym zakresie masy |
| `peakfrac_ratio` ± `_err` | (udział okna piku w h_acc)/(udział okna piku w próbce nietagowanej); biny w całości w oknie |
| `shape_slope_per_GeV` ± `shape_slope_err` | liniowe nachylenie stosunku kształtów acc/nietagowane względem masy [1/GeV] |
| `shape_boot_over_sampling_var`, `shape_overlap_frac`, `shape_cov_min_eig_rel` | diagnostyka kowariancji (ślady i najmniejsza istotna wartość własna) |

Systematyki są wspólne dla runu (wyznaczone w v2 dla regionu A∨C inkluzywnie) i przypisane każdemu wierszowi.
Proponowana zmiana składnika „ekstrapolacja” (RAPORT.md §2.4) jest w `jump_summary.json` → `proposal`, **nie** w tych kolumnach.

## `template_shapes.csv` — wiersz = bin masy

| pole | znaczenie |
|---|---|
| `run`, `region`, `species`, `pt_bin`, `mass_var` | jak wyżej |
| `bin_lo`, `bin_hi` | brzegi binu [GeV]; D⁺ 1,72–2,02 i D⁰ 1,75–2,00 po 10 MeV, D* Δm 0,139–0,160 po 1 MeV. Masa zapisana z krokiem 0,1 MeV; wartość równa brzegowi należy do binu wyższego, ostatni bin zawiera górny kraniec |
| `w_acc` | składowa przypadkowa: Σ p_acc po kandydatach w binie (oczekiwana liczba kandydatów z tagiem przypadkowym) |
| `w_acc_err_stat` | bootstrap dawców, błąd bezwzględny (zawiera wspólną normalizację; biny silnie skorelowane) |
| `n_all` | wszyscy kandydaci w binie bez wymogu tagu |
| `ratio_acc_over_all` ± `ratio_err` | stosunek kształtów unormowanych do 1 w pełnym zakresie; błąd z diagonali kowariancji testu kształtu |

Użycie: tło przypadkowe w próbce tagowanej = histogram `w_acc` (kształt i normalizacja), z błędem normalizacji z `template.csv`.

## `template.json`

`samples[]`: `rows` (= `template.csv`), `shapes` (= `template_shapes.csv`, listy per kategoria), `check_vs_v2`, `syst_rel_run` (pula, ekstrapolacja,
½Δ metod, κ), `usage_profile` (rozkład wag par „acc” w ξ_cal strony tagu × strony przeciwnej, ⟨nSelBL⟩ — wejście do ZADANIA 2).

## `jump_tests.json`, `jump_summary.json` (ZADANIE 2)

- `jump_tests.json`: surowe liczności. `test_i` — tagi obserwowane i przewidywane w puli (ξ_cal,s > 0,14) w przedziałach ξ_cal strony przeciwnej,
  podpróbki `hasD`/`noD`/`all`/`*_npo0`/`*_npo1` (bez/z protonem po stronie przeciwnej), tablice 2D, zależność od nSelBL/μ/LB;
  `mirror_all`, `mirror_npo0`, `nominal_npo0`, `act_fine` — warianty mieszania względem nominału (błąd ze wspólnego bootstrapu);
  `test_ii` — protony „niemożliwe”: tablice obs/pred o osiach [wycinek ξ_p, ξ_cal strony tagu, margines, ξ_cal strony przeciwnej, hasD].
- `jump_summary.json`: R = obs/pred z błędami, test symetrii stron, `proposal` (proponowana systematyka i f_acc).
