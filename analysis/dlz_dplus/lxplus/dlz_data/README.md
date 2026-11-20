# DLZ_DATA — H+DLZ dla D⁺ na danych + przesiew cech D⁰ (17.11.2026, projekt B)

Tylko pojedyncze D (D⁺→Kππ, D⁰→Kπ). Paczka nie czyta ani nie liczy par D–D (projekt A, ślepy).
Prerejestracja: `charm_work/DLZ_dplus/PROCEDURA.md`; kryteria: `charm_work/DLZ_dplus/DATA_2026-11-17/KRYTERIA.md` (zapisane przed fitem piku).

## DAG (`make_dag_dlz.sh` → `dlz.dag`, batch `charm_dlz`)
| węzeł | co | zasoby |
|---|---|---|
| `frm_fd_*` (3), `frm_mb_*` (6), `frm_sb_*` (10) | ramki z h5 v2.0 (L_xy ≥ 0,25 mm, weto D*): sygnał fd truth 0/1; mb wszystko; SB danych 428770 \|m−m0\| > 70 MeV, 6 % zdarzeń | 1 cpu, 8 GB |
| `train` | H / H+DLZ / H+dalitz, seed 12345 (+ H, H+DLZ seed 23456); 5 foldów po zdarzeniach; XGB H; cap tła 3·10⁶; progi WP z sygnału OOF | 8 cpu, 32 GB, workday |
| `ext_data_428770_*` (20), `ext_data_427929_*` (10) | `apply.py`: score modelu foldu (hasz (run,event)); wiersze „luźne” (≥ WP95 dowolnego modelu) + histogramy masy wszystkich kandydatów | 2 cpu, 10 GB |
| `merge` | `dlz_<próbka>.npz`, `dlz_summary.json` | 16 GB |
| `ana` | `mc_checks.py` (R_sculpt MC truth ≥ 8, closure fitu na mb, AUC trening/OOF) + `analyze.py <run>`: fity masy (2 Gaussy z MC + exp(pol2)), bootstrap po zdarzeniach (200), zysk σ(N)/N, ΔSF_rel, pik, rzeźbienie w SB, sPlot (ε warunkowe, cechy 1D, C2ST sygnału), połowy danych | 4 cpu, 16 GB |
| `d0` | h5 v1.2 (fd wszystkie, mb 40, dane 6 plików) → ramki FEAT + nowe cechy (`d0/d0_feats.py`) → `study` (mc: 12 wariantów; data: 5), `transfer`, `sculpt`, `d0_sculpt` (fałszywy sygnał, regresja masy), `modelling` (ΔC2ST) | 8 cpu, 32 GB, workday |
| `fin` | `DONE` | |

Wyjście: `/eos/user/k/klyczak/charm_q/dlz_data/{parts,train,ana,d0,d0_scores,dlz_*.npz,DONE}`. Scratch: `mktemp` + `trap`. AFS: tylko logi `Q_lxplus/logs/dlz_data`.
Etap ORCH: `ORCH/stages/dlz_data.sh` (pobiera `ana/`, `d0/`, `train/thr.json`, `train/models/*/meta.json`, `dlz_summary.json` → `DLZ_dplus/DATA_2026-11-17/runs/lxplus/`).

## Użycie (lxplus; wysyłka przez ORCH)
```
./submit_dlz.sh check      # pliki, h5, środowisko, kolejka, AFS
./submit_dlz.sh submit     # check + DAG + condor_submit_dag (czyste env)
./submit_dlz.sh status
```
Zmienne: `dlz_common.sh` (`DLZ_CAP`, `DLZ_SBFRAC`, `DLZ_SETS`, `DLZ_NBOOT`, `DLZ_D0=0` wyłącza węzeł D⁰, `DLZ_FAST=1` test 40 drzew).

## Pliki
`dlzlib.py`, `apply.py`, `analyze.py`, `fitlib.py`, `mc_checks.py` = kopie z `DLZ_dplus/DATA_2026-11-17/scripts/` (źródło prawdy tam; w kopii `dlzlib.py`
zmieniona tylko domyślna ścieżka `CHARMBKG_SRC`). `d0/` = skrypty `FEAT_charmnet/scripts` (kopie; zmieniona tylko ścieżka charmbkg) + `DLZ_d0/scripts`.
`src/charmbkg` = migawka `~/Projects/charm/charm-bkg/charmbkg` (bez `bkggen`; `src/SNAPSHOT.txt`).

## Test (17.11)
- lokalnie: `dlz_job.sh frm/train/ext/merge/ana` z `DLZ_NOENV=1 DLZ_FAST=1` na 2+2+1+6 plikach;
- lxplus `/tmp/klyczak_dlzt` (bez condora, bez EOS/AFS): to samo na plikach z `/eos/.../charm_h5_v20/parts` + węzeł `d0` w trybie `D0_MODE=screen`.
