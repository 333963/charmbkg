# MCTMPL — histogramy do testu szablonów tła z MC (zadanie #7 planu prac; projekt B, 17.11.2026)

Cel: sprawdzić na danych szablony tła, które w ficie TRExFitter głównego łańcucha mają kształt i normalizację z MC bez systematyk (`Reflection`, `CharmMismatch`,
`DecayModeMismatch`, `PartComb`; swobodny tylko `Comb`). Opis, kryteria i struktura wyniku: `charm_work/MC_templates/RAPORT.md`, `KRYTERIA.md`.

## Co robi
`mctmpl_hist.py` (charmtag z `ev_dstar/src`, modele `runs/lxplus_v11/models/{D0_H,Dplus_H}`): D⁰ (L_xy ≥ 0,25 mm) i D⁺ (L_xy ≥ 0,5 mm), pT 3–30, |η| < 2,2,
okna 1,70–2,05 / 1,70–2,20 GeV → **tylko histogramy** (npz ~0,1–0,3 MB na węzeł), osie [pT, WP H (95/80/50), klasa trygera | kategoria prawdy, …]:
masa × komórka aktywności (nSelBL × pT toru wiodącego), |impact| i ct w SR/SB, tag D* (Δm × masa poprawnej hipotezy / bliźniaka), D⁺: Δm_min, m(KK), m(KKπ) w oknie φ,
m(pKπ), m(Kππ) w oknie D*/φ. Dane: GRL; **zdarzenia regionu sygnału głównego łańcucha (tag AFP ∧ nPileupVertices = 0 ∧ ξ_cal ≤ ξ_AFP + 0,02) są odrzucane i nieliczone.**
Progi WP: `wp_thresholds.json` (kwantyle score sygnału prawdy z 4 plików fd_nd, 17.11; odtworzenie: `mctmpl_hist.py --thr-out`).
DAG `charm_mctmpl`: 14 + 10 + 8 + 3 węzłów (longlunch, 2 CPU) → `merge` (suma + `mctmpl_analyze.py`) →
`/eos/user/k/klyczak/charm_q/mctmpl/{mctmpl_<próbka>.npz, mctmpl_summary.json, ana/WYNIKI.md, ana/results.json, DONE}`.

## Uruchomienie (tylko przez ORCH — `ORCH/stages/mctmpl.sh`; wpis do `pipeline.yaml` dopisywany osobno)
```
  - name: mctmpl      # MCTMPL: histogramy do testu szablonów tła MC (zad. #7) z h5 v2.0 → charm_q/mctmpl → MC_templates/runs/lxplus; analiza MC_templates/analyze.py (projekt B)
    needs: [h5_v20]
    max_hours: 6
```
Ręcznie na lxplus: `./submit_mctmpl.sh check | submit | status` (submit w czystym env, `condor_q` w check, FORCE=1 nadpisuje DONE).
Zmienne: `MT_GROUP`, `MT_NFILES`, `MT_OUT`, `MT_WP` (mctmpl_common.sh).

## Test (17.11, lxplus /tmp, posprzątane)
| próbka | pliki | czas | RAM | kandydaci D⁰ / D⁺ | npz |
|---|---|---|---|---|---|
| mb_nd | 2 | 16 s | 1,0 GB | 104 712 / 277 646 | 0,16 MB |
| fd_nd | 1 | 15 s | 1,2 GB | 104 783 / 268 504 | 0,24 MB |
| data_428770 | 1 | 30 s | 2,3 GB | 244 953 / 709 980 | 0,27 MB |
| data_427929 | 4 | 10 s | 0,5 GB | 40 020 / 111 902 | 0,07 MB |
Ekstrapolacja: węzeł ≤ ~5 min CPU (+ kolejka), wynik < 5 MB. `submit_mctmpl.sh check` OK, `make_dag_mctmpl.sh` 35 + 1 węzłów, `mctmpl_job.sh` DRYRUN OK,
`mctmpl_merge.py` OK, `mctmpl_analyze.py` w env LCG OK (iminuit 2.32; 9 s na plikach testowych), `bash -n` OK.
Pliki testowe lokalnie: `MC_templates/runs/test/` (`out/`, `out_minn/` — tylko test kodu, liczby bez znaczenia).
Uwaga: `make_dag_mctmpl.sh` tworzy `Q_lxplus/logs/mctmpl/` (już istnieje, pusty) i `mctmpl/lists/` (AFS, małe); niczego nie kasuje.

## Analiza lokalna
`~/Documents/charm/.venv_h5/bin/python charm_work/MC_templates/analyze.py` (domyślnie `--in MC_templates/runs/lxplus --out MC_templates/out`) → `out/WYNIKI.md`,
`out/results.json`. To ten sam kod co w węźle merge; gotowy wynik z condora: `MC_templates/runs/lxplus/ana/WYNIKI.md`.
