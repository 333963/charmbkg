#!/bin/bash
# DLZ_d0: łańcuch (3 wątki). Uruchamiać po zwolnieniu rdzeni.
cd /home/k/Documents/charm/charm_work/DLZ_d0
P=/home/k/Documents/charm/.venv_h5/bin/python
R=scripts/d0_run.py
export OMP_NUM_THREADS=3
set -x
$P $R study --test mc --variants H,H_seed2,H+ALL_noiso,H+ptratio,H+d0prod,H+ip,H+time,H+kin2,H+kinall,H+pidbb,H+NOKIN,H+NOKIN_pid,H+ALLnew > logs/study_mc.log 2>&1
$P $R transfer > logs/transfer.log 2>&1
$P $R sculpt mc > logs/sculpt_mc.log 2>&1
$P scripts/d0_sculpt.py > logs/d0_sculpt.log 2>&1
$P $R modelling > logs/modelling.log 2>&1
echo D0_DONE
