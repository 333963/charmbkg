#!/bin/bash
cd /work/b1d0/work
export CHARMTAG_SRC=/work/b1d0/in/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run() { python3 /work/b1d0/work/b1_study_d0.py "$@" > "/work/b1d0/work/logs/$(echo "$@" | tr ' -' '__').log" 2>&1; }
( run --mode smooth --kind hgb; run --mode real_t27 --kind hgb --mc_truth_min 2; run --mode cat_local --kind hgb --seed 1 ) &
( run --mode swap_local --kind hgb --swap_factor 0.5; run --mode swap_plus --kind hgb --mc_truth_min 2; run --mode smooth --kind hgb --strength 1.6 --seed 1 ) &
wait
echo FIX_DONE > /work/b1d0/work/logs/FIX_DONE
