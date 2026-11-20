#!/bin/bash
# DLZ_d0: przesiew (tryb screen: lr 0,15, cap 150k, 1 wątek) wszystkich grup; potem potwierdzenie najlepszej w pełnym protokole H (run_d0.sh).
cd /home/k/Documents/charm/charm_work/DLZ_d0
P=/home/k/Documents/charm/.venv_h5/bin/python
R=scripts/d0_run.py
export D0_MODE=screen D0_NJ=${D0_NJ:-1} OMP_NUM_THREADS=${D0_NJ:-1}
set -x
$P $R study --test mc --variants H,H_seed2,H+d0prod,H+ip,H+time,H+NOKIN,H+kin2,H+kinall,H+pidbb,H+NOKIN_pid,H+ALLnew > logs/screen_study_mc.log 2>&1
$P $R transfer > logs/screen_transfer.log 2>&1
$P $R sculpt mc > logs/screen_sculpt_mc.log 2>&1
$P scripts/d0_sculpt.py > logs/screen_d0_sculpt.log 2>&1
$P $R modelling > logs/screen_modelling.log 2>&1
echo SCREEN_DONE
