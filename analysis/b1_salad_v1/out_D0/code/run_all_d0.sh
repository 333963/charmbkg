#!/bin/bash
cd /work/b1d0/work
export CHARMTAG_SRC=/work/b1d0/in/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run() { python3 /work/b1d0/work/b1_study_d0.py "$@" > "/work/b1d0/work/logs/$(echo "$@" | tr ' -' '__').log" 2>&1; }
( run --mode smooth --kind hgb; run --mode cat_local --kind hgb; run --mode swap_local --kind hgb --swap_factor 1.5; run --mode swap_local --kind hgb --swap_factor 0.5; run --mode cat_step --kind hgb; run --mode cuts_local --kind hgb; run --mode smooth --kind hgb --leak 0.02; run --mode smooth --kind sk_mlp; run --mode cat_local_MM --kind hgb ) &
( while pgrep -f b1_data_d0.py >/dev/null; do sleep 5; done; run --mode cat_global --kind hgb; run --mode swap_plus --kind hgb; run --mode real_t27 --kind hgb; run --mode cat_local+smooth --kind hgb; run --mode smooth --kind hgb --leak 0.05; run --mode cat_local --kind sk_mlp; run --mode cat_step --kind sk_mlp; run --mode swap_local --kind sk_mlp --swap_factor 1.5 ) &
wait
echo ALL_DONE > /work/b1d0/work/logs/ALL_DONE
