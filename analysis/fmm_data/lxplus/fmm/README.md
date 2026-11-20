# FMM — pary bliźniacze D0 K↔π z h5 v2.0 → f_MM pod pikiem (projekt B, 17.11.2026)

Otwarty punkt METHODS_NOTE §7 / B1: „f_MM in data can be measured directly from `hypo_swap` — not yet done”.
**Ustalenie:** `Dmeson_hypo_swap` to flaga prawdy MC (w danych −1), a MisMatch* (truth 8–11) to pary „jeden tor z hadronu powabnego +
tor przypadkowy”, nie zamiany K↔π. Obserwabla bliźniacza (m_swap / bliźniak w oknie) ma AUC ≈ 0,52 dla MisMatch* vs Other
(pełne MC v1.2, `HYP_d0/runs/lxplus/mb_nd.parquet`: bliźniak w oknie 12,6 % vs 10,7 %) → **nie identyfikuje** MisMatch*.
Paczka dostarcza więc wejście do **ograniczenia** (przedziału) Δf_MM, nie pomiaru — patrz `FMM_data/analyze.py`, `FMM_data/RAPORT.md`.

## Co robi
`fmm_extract.py` (logika = `B_trig/btrig_extract.py`, ten sam model `runs/lxplus_v11/models/D0_H`, charmtag z `ev_dstar/src`):
D0 (decayType 10), L_xy ≥ 0,25 mm, pT 3–30, pełne okno v2.0 1,70–2,05 GeV → npz: m, m_swap (analitycznie), bliźniak (n_grp, tw_m, tw_score),
score H + fold, wszystkie cechy H (float16), dE/dx i kinematyka torów, prawda MC (truth, hypo_swap, pdgId/rodzice córek), bity triggera,
GRL, nSelBL, lead_pt, μ. Próbki: data_428770 (30 % zdarzeń, hash), data_427929 (100 %), mb_nd (100 %), fd_nd (tylko truth ≤ 7).
DAG `charm_fmm`: 10 + 10 + 8 + 3 węzłów (longlunch, 2 CPU) → `merge` → `/eos/user/k/klyczak/charm_q/fmm/{fmm_<próbka>.npz, fmm_summary.json, DONE}`.

## Uruchomienie (tylko przez ORCH — `ORCH/stages/fmm.sh`; wpis do `pipeline.yaml` dopisywany osobno)
```
  - name: fmm         # FMM: pary bliźniacze D0 z h5 v2.0 → charm_q/fmm → FMM_data/runs/lxplus; analiza FMM_data/analyze.py (projekt B)
    needs: [h5_v20]
    max_hours: 6
```
Ręcznie na lxplus: `./submit_fmm.sh check | submit | status` (submit w czystym env, `condor_q` w check, FORCE=1 nadpisuje DONE).
Zmienne: `FMM_FRAC`, `FMM_GROUP`, `FMM_NFILES`, `FMM_OUT`, `FMM_LXY` (fmm_common.sh).

## Test (17.11, lxplus /tmp, posprzątane)
| próbka | pliki | czas | RAM | wiersze | npz |
|---|---|---|---|---|---|
| mb_nd | 2 | 18 s | 0,7 GB | 110 467 | 10,4 MB |
| fd_nd | 2 | 18 s | 0,9 GB | 39 997 | 3,8 MB |
| data_428770 (frac 0,3) | 1 | 29 s | 1,2 GB | 76 930 | 6,9 MB |
| data_427929 | 4 | 9 s | 0,4 GB | 42 277 | 3,9 MB |
Ekstrapolacja: węzeł ≤ ~5 min CPU (+ kolejka), wynik ≈ 0,67 + 0,45 + 0,41 + 0,05 ≈ **1,6 GB** do pobrania. `tw_m` = `m_swap` do < 0,4 MeV.
`submit_fmm.sh check` OK, `make_dag_fmm.sh` 31 + 1 węzłów, `fmm_job.sh` DRYRUN OK, `fmm_merge.py` OK, `bash -n` OK.
Pliki testowe lokalnie: `FMM_data/runs/test/` (+ `out/RAPORT_TEST.md` — tylko test kodu).
Uwaga: `make_dag_fmm.sh` tworzy `Q_lxplus/logs/fmm/` i `fmm/lists/` (AFS, małe); niczego nie kasuje.

## Analiza lokalna
`~/Documents/charm/.venv_h5/bin/python charm_work/FMM_data/analyze.py` (domyślnie `--in FMM_data/runs/lxplus --out FMM_data`) → `RAPORT.md`, `results.json`.
