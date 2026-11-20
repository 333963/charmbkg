#!/bin/bash
cd /work/b1d0/work
export CHARMTAG_SRC=/work/b1d0/in/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run() { python3 /work/b1d0/work/b1_study_d0.py "$@" > "/work/b1d0/work/logs/$(echo "$@" | tr ' -' '__').log" 2>&1; }
while [ ! -f /work/b1d0/work/logs/ALL_DONE ]; do sleep 10; done
( run --mode smooth --kind hgb --strength 1.6; run --mode smooth --kind sk_mlp --strength 1.6; run --mode cat_local+smooth --kind hgb --strength 1.6 ) &
( run --mode smooth --kind hgb --strength 1.6 --leak 0.02; run --mode swap_local --kind hgb --swap_factor 0.5 --seed 1; run --mode smooth --kind hgb --seed 1 ) &
wait
echo EXTRA_DONE > /work/b1d0/work/logs/EXTRA_DONE
