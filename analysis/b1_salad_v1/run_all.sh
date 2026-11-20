#!/bin/bash
cd /work/b1/work
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run() { python3 b1_study.py "$@" > "logs/$(echo "$@" | tr ' -' '__').log" 2>&1; }
( run --mode smooth --kind hgb; run --mode cat_local --kind hgb; run --mode cat_step --kind hgb; run --mode cuts_local --kind hgb; run --mode smooth --kind hgb --leak 0.02; run --mode smooth --kind sk_mlp ) &
( run --mode cat_global --kind hgb; run --mode cat_local+smooth --kind hgb; run --mode smooth --kind hgb --leak 0.05; run --mode cat_local --kind sk_mlp; run --mode cat_step --kind sk_mlp ) &
wait
echo ALL_DONE
