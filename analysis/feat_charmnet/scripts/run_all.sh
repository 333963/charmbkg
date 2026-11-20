#!/bin/bash
# FEAT_charmnet: pełny bieg lokalny (sekwencyjnie, 8 wątków XGB)
cd /home/k/Documents/charm/charm_work/FEAT_charmnet
P=/home/k/Documents/charm/.venv_h5/bin/python
set -x
$P scripts/study.py --species D0 --test mc --cap 300000 --nboot 200
$P scripts/study.py --species Dplus --test mc --cap 300000 --nboot 200
$P scripts/study.py --species D0 --test data --cap 300000 --nboot 200 --variants H,H_seed2,H+iso2,H+ALL,H+ALL_noiso
$P scripts/study.py --species Dplus --test data --cap 300000 --nboot 200 --variants H,H_seed2,H+iso2,H+ALL,H+ALL_noiso
echo RUN_ALL_DONE
