#!/bin/sh
# SIG_constr: kolejka lokalna nr 2 (1 rdzeń)
P=/home/k/Documents/charm/.venv_h5/bin/python
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
$P s5_transfer.py > out/s5_transfer.log 2>&1
$P s6_perp.py > out/s6_perp.log 2>&1
$P s3_variants.py 8:J12only > out/s3_8J12.log 2>&1
$P s1_neff.py > out/s1_neff.log 2>&1
$P s6_perp.py > out/s6_perp.log 2>&1
