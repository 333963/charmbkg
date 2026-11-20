# AFP_EXT — slim evtable + kandydaci D z h5 v2.0 dla AFP-MIX (projekt B, 17.11.2026; przygotowane, NIEWYSŁANE)

Wejście (tylko odczyt): `/eos/user/k/klyczak/charm_h5_v20/parts/{data_428770,data_427929}/*_none.h5` (datasety `events`, `evtable`; schemat `../h5_v20/README.md`,
nazwy sprawdzone odczytem 17.11). Wyjście: `/eos/user/k/klyczak/charm_q/afp_ext/` — `<próbka>_events.parquet`, `<próbka>_cands.parquet`,
`<próbka>.meta.json`, `afp_summary.json`, `DONE`; części w `parts/` (ORCH ich nie pobiera).

- `afp_extract.py` (numpy/h5py/pyarrow): na plik h5 — `evtable` → zdarzenia slim (event_number, lumiBlock, actual_mu, passGRL_MinBiasAFP,
  nPileupVertices, nSelBL, xi_cal_A/C, afp_nProtons/nInWindow/xi1/xi2/tag_A/C, 11 × `trigPassed_*`, 5 × `trigTBP_*`, n_cand, n_cand_sel);
  `events` → kandydaci slim (event_number, Dmeson_species 0/1/3, Dmeson_mass, Dmeson_delta_mass, Dmeson_pt, Dmeson_eta, Dmeson_Lxy).
  Selekcja kandydatów (score nie ma w h5 → luźne okna): D⁺ m ∈ [1,72; 2,02]; D⁰ m ∈ [1,75; 2,00] i L_xy ≥ 0,25 mm; D* Δm ∈ [0,139; 0,160].
  Kwantyzacja: m, Δm 0,1 MeV; pT 10 MeV; η, L_xy 0,001. `run_number` i preskale `trigPS_*` (stałe w LB) → meta json (`trigps`).
  Kontrole: komplet kolumn, jeden run na plik, event_number unikalny, każdy kandydat ma zdarzenie w evtable.
- `afp_merge.py` (finish): strumieniowe scalenie części; każdy plik h5 dokładnie raz, liczba plików = `AFP_NFILES` (97 / 458), sumy wierszy.
- Tryger: wszystkie 11 łańcuchów `TRIG_STORE` jest w obu runach, w tym **`HLT_noalg_L1AFP_A_OR_C`** i `HLT_mb_sptrk_pt4_L1AFP_A_OR_C`
  → zdarzenia wyłącznie z trygera AFP można wykluczyć (G0). Łańcuchów spoza listy nie ma w h5: w 428770 ~17 % zdarzeń nie ma żadnego
  zapisanego łańcucha (analiza je odrzuca w trybie `nonafp`).

## Uruchomienie (lxplus, `~/private/charm_q/Q_lxplus/afp_ext`; normalnie ORCH etap `afp_ext`)
`./submit_afp.sh check|submit|status|fetch` — batch `charm_afp`; DAG: 2 próbki × `AFP_NPARTS`=8 węzłów (1 CPU, 4 GB, longlunch, RETRY 2) → `afp_finish`.
`check`: liczba `*_none.h5` = oczekiwana (nie wymaga `parts/DONE` etapu h5_v20), `condor_q`, brak `DONE`, kwota AFS < 85 %.
Submit w czystym env (`env -u LD_LIBRARY_PATH -u PYTHONPATH`, `condor_submit_dag -f`). `DRYRUN=1 ./submit_afp.sh submit` = DAG + komendy.
Scratch `mktemp` + trap (plik h5 kopiowany xrdcp pojedynczo, fallback FUSE), wyniki tylko EOS, logi `Q_lxplus/logs/afp_ext`.
Zmienne: `AFP_OUT`, `H5V20`, `AFP_SAMPLES`, `AFP_NPARTS`, `AFP_NFILES`, `AFP_DISK`, `AFP_MAX_FILES` (test).

## Test (17.11, lxplus /tmp, skasowane)
2 pliki 428770 + 20 plików 427929: ekstrakcja 6–20 s / plik 428770 (84–95 k zdarzeń, 0,53–0,59 M kandydatów slim z 1,14 M), 1,4 s / plik 427929;
merge OK (i poprawnie odmawia przy brakującej części); `check`, `make_dag_afp.sh` (16 węzłów + finish), `DRYRUN=1 afp_job.sh` OK.
Szacunek pełnego biegu: 428770 ≈ 8 M zdarzeń → ~180 MB + kandydaci ~55 M → ~400 MB; 427929 ≈ 29 M zdarzeń → ~440 MB + kandydaci ~16 M → ~130 MB;
**razem ≈ 1,2 GB** (< 3 GB). Czas: węzeł 428770 (12–13 plików) 2–5 min, 427929 (57–58 plików) ~2 min, finish ~3 min; ~0,5 h z kolejką.

## Po biegu (lokalnie)
`~/Documents/charm/.venv_h5/bin/python ~/Documents/charm/charm_work/AFP_mix/run_data.py --save-p` → `AFP_mix/runs/lxplus/out/RAPORT_DANE.md` (WSTĘPNE).
