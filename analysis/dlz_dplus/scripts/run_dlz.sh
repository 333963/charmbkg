#!/bin/bash
cd /home/k/Documents/charm/charm_work/DLZ_dplus
P=/home/k/Documents/charm/.venv_h5/bin/python
R=scripts/dlz_run.py
set -x
$P $R study --test mc --variants H,H_seed2,H+dalitz,H+ptratio,H+maxdR,H+absimp,H+ALL_noiso,H+ALL_noiso_nopid,H+dalitz_ptratio > logs/study_mc.log 2>&1
$P $R transfer > logs/transfer.log 2>&1
$P $R sculpt mc > logs/sculpt_mc.log 2>&1
$P $R modelling > logs/modelling.log 2>&1
$P $R study --test data --variants H,H_seed2,H+ALL_noiso,H+dalitz,H+ALL_noiso_nopid > logs/study_data.log 2>&1
$P $R sculpt data > logs/sculpt_data.log 2>&1
echo DLZ_DONE
