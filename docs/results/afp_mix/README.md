# AFP_mix — tło przypadkowego tagu AFP z mieszania protonów na danych (projekt B, WSTĘPNE)

**v2 (17.11.2026).** Metoda, bramki, toy: repo `~/Projects/charm/charm-bkg` (`charmbkg/afp_mix.py`, `docs/afp_mixing.md`). Tutaj: bieg na danych.
v1: `*.bak_v1` obok plików (README, `run_data.py`, `test_run_data.py`) oraz wynik v1 w `runs/lxplus/out/` (nietknięty).
Podstawa zmian: `charm_work/REV_afp_mix/RAPORT.md`, `charm_work/RECON_lit_B/RAPORT.md`.

## Bieg
Wejście: `runs/lxplus/<próbka>_{events,cands}.parquet` (z lxplus `Q_lxplus/afp_ext`, już pobrane).

    ~/Documents/charm/.venv_h5/bin/python ~/Documents/charm/charm_work/AFP_mix/run_data.py --in runs/lxplus --out runs/lxplus/out_v2

→ `runs/lxplus/out_v2/RAPORT_DANE_v2.md`, `result_v2.json`, `run_data_v2.log`. Czas ~15 min (oba runy), RAM szczytowo ~15 GB.
Opcje: `--samples`, `--save-p` (p_tag/p_acc per zdarzenie), `--tag-def`, `--trigger-mode`, `--n-boot`, `--n-pseudo`.

## Co liczy `run_data.py` v2
- Selekcja: GRL ∧ ≥ 1 zapisany łańcuch bez L1AFP ∧ `nPileupVertices == 0`. **Zakres: tylko trygery nie-AFP (W2).**
- Estymator A: `afp_mix.run_procedure` — warstwy (LB × nSelBL × „z D/bez D”), pule jednostronne ξ_cal > 0,14, overlay union na parach.
- Klasy trygera (3 = L1J12, 2 = j20, 1 = sptrk_pt4_L1MBTS_2, 0 = pozostałe L1MBTS): G2 i wariant puli.
- B(i): łańcuchy L1MBTS (quasi-minbias; `HLT_noalg_L1RD0_FILLED` ma 0 zdarzeń w obu runach), weto po obu stronach, nSelBL ≥ 20;
  gdy łańcuchy MBTS to ≥ 50 % próbki (427929) — tylko zdarzenia bez D. B(ii): obszar kontrolny ξ_cal ≤ 0,12, Δ > 0,045.
- Wynik: f_acc z rozbiciem błędu (stat, pula, ekstrapolacja, ½Δ metod), bramki v2 z opisem „co testuje / czego NIE”, tabele per region,
  μ, LB, gatunek × pT (pik i sidebandy; bez fitów), porównanie z v1.
- **Żadnych fitów masy, uzysków sygnału ani przekrojów** — tylko frakcje tła przypadkowego i liczności tagów.

## Wynik v2 (WSTĘPNE)
| run | f_acc (A∨C) | werdykt | v1 |
|---|---|---|---|
| 427929 | 0,0087 ± 0,0044 | PASS (G9 N/A) | 0,0087 |
| 428770 | 0,0477 ± 0,0135 | PASS (G9 N/A) | 0,0459 |

## Test
`~/Documents/charm/.venv_h5/bin/python -m pytest -q -p no:cacheprovider charm_work/AFP_mix/test_run_data.py` (5 testów): toy wyeksportowany
do formatu afp_ext (`run_data.py --make-toy DIR`) → domknięcie z prawdą, zgodność z `afp_mix.mix_region`, wykluczenie trygera AFP.
