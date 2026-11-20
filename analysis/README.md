# analysis/ — analysis scripts outside the package

Two kinds of content:

1. **Driver scripts** that use `charmbkg`: `afp_mix/`, `tagtrig/` (described below, in Polish).
2. **Archived analysis steps** (from 1.2.0; older working directories added in 1.3.0): the scripts of each study of 4–20.11.2026
   exactly as they were run, one directory per topic with a short `README.md`; their reports, tables and numbers are in
   [`docs/results/<topic>/`](../docs/results/README.md). These scripts keep the paths of the original working area, read inputs that
   are not in the repository and are not expected to run from here; `python -m compileall -q analysis` passes. They are outside the
   scope of `ruff` and of the default `pytest` run. Which script produced which headline number: table „Wynik → skrypt → raport” in
   the top-level [overview](../PRZEGLAD.md). Files left out on purpose (code derived from the charmpy GitLab repository, code of
   project A, snapshots of packages, binaries) are named in the `README.md` of the topic concerned.

| group | topics |
|---|---|
| network and classifier inputs | `dnet` (`lxplus/dnet1` … `dnet6`), `feat_charmnet`, `dlz_dplus`, `dlz_d0`, `cn_clone`, `sig_constr` |
| working point and fit | `wpscan`, `wptight`, `fitbias`, `fitbias_dplus`, `m8_recalc`, `bkgcat`, `catfit`, `fit_scorebins`, `lowpt_calib`, `peak_pos`, `peak_lxy010`, `konstr`, `d23_plan`, `d23_d0_check` |
| background under the peak | `b1_salad_v1`, `salad_427929`, `est_spread`, `mc_templates`, `dstws`, `fmm_data` |
| efficiency | `sfabs`, `ev_dstar`, `iso_pis`, `iso_sf`, `iso_ped`, `veto_pv` |
| trigger and data sample | `b_trig` (+ `TRKJET`, `TRKJET2`, `AFPSTREAM`, `ISO_origin`, `BKG_427929`), `tagreg`, `tagtrig/lxplus` |
| normalisation and cross-section | `lumi`, `cver`, `xsec_atlas`, `xsec_final` (`v2/` = final), `reg428` |
| proton-tagged measurement | `tagged`, `tagxsec`, `afp_acceptance`, `afp_mix` (+ `AFPTRIG`, `FORCED_2026-11-18`, `TEMPLATE`, `lxplus/afp_ext`), `sd_pred`, `sd_var` |
| figures | `presentation` |

Topic with a report only (no scripts): `br_frag` — see [`docs/results/`](../docs/results/README.md).

## afp_mix/ — tło przypadkowego tagu AFP (metoda: [docs/afp_mix_method.md](../docs/afp_mix_method.md))

Kolejność:

1. `afp_extract.py --sample S --h5dir DIR --part K --nparts N --out-prefix WD/S.pK` — z plików h5 (datasety `evtable`, `events`)
   robi odchudzoną tabelę zdarzeń (AFP, ξ_cal, aktywność, bity trygera) i tabelę kandydatów D (luźne okna masy) w parquet; jedna część
   = co N-ty plik próbki.
2. `afp_merge.py --out OUT --wd WD --samples "a b" --nparts N [--nfiles "a:n b:m"]` — scala części, sprawdza komplet plików i sumy wierszy
   → `<próbka>_{events,cands}.parquet`, `<próbka>.meta.json`.
3. `run_data.py --in DIR --out DIR` — selekcja zdarzeń (GRL, tryger nie-AFP, bez wierzchołków pile-up), klasy trygera, potem
   `charmbkg.afp_mix.run_procedure` (estymator A, drugie estymatory B(i)/B(ii), bramki, rozbicie błędu) i tabele per gatunek / pT / μ / LB;
   raport i json w `--out`. `run_data.py --make-toy DIR` zapisuje toy (`charmbkg.afp_mix_toy`) w tym samym formacie parquet.

Test (toy w formacie z kroku 2, domknięcie z prawdą, zgodność z `afp_mix.mix_region`, wykluczenie trygera AFP):
`python -m pytest -q analysis/afp_mix/test_run_data.py` (wymaga `pyarrow` — od 1.3.0 w `.[dev]`; test jest też krokiem CI).

Dodane 20.11.2026 (zarchiwizowane, nie uruchamiają się z repo — wejścia w obszarze roboczym): `lxplus/afp_ext/` (paczka HTCondor
kroków 1–2), `AFPTRIG/` (metoda kształtu dla strumienia z trygerem AFP), `FORCED_2026-11-18/` (wykluczenie tagów wymuszonych),
`TEMPLATE/` (szablon N_tag / N_acc / f_acc). Raporty: [`docs/results/afp_mix/`](../docs/results/afp_mix/).

## tagtrig/ — bity wszystkich łańcuchów trygera dla zdarzeń z tagiem AFP

1. `tt_bits.py --sample S --rootdir DIR --part K --nparts N --out-prefix WD/S.pK` — z drzewa `CharmAnalysis` (ROOT) zapisuje dla zdarzeń
   z tagiem AFP spakowane bity łańcuchów HLT i pozycji L1, liczniki kandydatów oraz, per blok świetlności, liczniki zdarzeń i preskale.
2. `tt_merge.py --out OUT --wd WD --samples S --nparts N [--nfiles S:n]` — scala części → `bits_<run>.parquet`, `ps_perLB_<run>.csv`,
   `lb_summary_<run>.csv`, `chains_<run>.json`, `summary.json` (same liczniki).

Test: `python analysis/tagtrig/test_tt.py` (wymaga `uproot`, `awkward`, `pyarrow`).

Dodane 20.11.2026: `lxplus/tagtrig/` — paczka HTCondor tych kroków oraz `tt_cool.py`, `tt_lumi.sh` (preskale z COOL i świetlność
per blok; wymagają środowiska ATLAS na lxplus). Wyniki: [`docs/results/tagreg/`](../docs/results/tagreg/) (`runs/lxplus/`).

## Uwagi

- Testy z tego katalogu nie wchodzą do domyślnego `pytest` (`testpaths = tests`). W CI uruchamiany jest tylko
  `analysis/afp_mix/test_run_data.py`; `tagtrig/test_tt.py` i `b_trig/tests/` mają dodatkowe zależności.
- `PYTHONPATH` nie jest potrzebny, gdy pakiet jest zainstalowany (`pip install -e .`); `run_data.py` dodaje korzeń repo sam
  (albo katalog ze zmiennej `CHARMBKG`).
