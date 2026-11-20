#!/bin/bash
# Generuje cn.dag: ext_<próbka>_<g> (CN_GROUP plików h5_v20 na węzeł) → train. NIE wysyła.
# Listy w lists/<znacznik czasu>_g<CN_GROUP>/ (nowy katalog przy każdym wywołaniu → bez kasowania i bez starych list w DAG-u).
set -euo pipefail
CN=$(cd "$(dirname "$0")" && pwd); export CN
source "$CN/cn_common.sh"
: "${CN_LISTS:=$CN/lists/$(date -u +%Y%m%dT%H%M%S)_g$CN_GROUP}"
mkdir -p "$Q/logs/cnclone" 2>/dev/null || true; mkdir -p "$CN_LISTS"
S=$CN/cn.sub
PENV="CN_OUT=$CN_OUT H5V20=$H5V20 CN_SAMPLES=$(echo $CN_SAMPLES | tr ' ' ,) CN_FRAC_DATA=$CN_FRAC_DATA CN_FRAC_MB=$CN_FRAC_MB CN_LXY=$CN_LXY CN_FOLDS=$CN_FOLDS CN_EPOCHS=$CN_EPOCHS CN_CPUS=$CN_CPUS CN_MAX_SIG=$CN_MAX_SIG CN_MAX_BKG=$CN_MAX_BKG CN_FD_BKG=$CN_FD_BKG"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${CN_DAG:-$CN/cn.dag}
EXT=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $CN_SAMPLES; do
    ls "$H5V20/$s/"*_"$CN_PRESEL".h5 2>/dev/null | grep -v '\.part$' > "$CN_LISTS/$s.all" || true
    [ -s "$CN_LISTS/$s.all" ] || { echo "STOP: brak plików h5_v20 dla $s" >&2; exit 1; }
    split -d -a 3 -l "$CN_GROUP" "$CN_LISTS/$s.all" "$CN_LISTS/${s}_"
    for l in "$CN_LISTS/${s}_"[0-9][0-9][0-9]; do
      g=$((10#${l##*_})); n="ext_${s}_$g"; EXT+=("$n")
      mem=4GB; case $s in data_428770) mem=8GB ;; esac
      echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"ext $s $l $g\" cpus=\"1\" mem=\"$mem\" disk=\"3GB\" flavour=\"longlunch\" $V"; echo "RETRY $n 2"
    done
  done
  echo "JOB train $S"; echo "VARS train node=\"train\" args=\"train\" cpus=\"$CN_CPUS\" mem=\"16GB\" disk=\"8GB\" flavour=\"workday\" $V"; echo "RETRY train 1"
  echo "PARENT ${EXT[*]} CHILD train"
} > "$DAG"
echo "cn.dag ($DAG): ekstrakcja ${#EXT[@]} węzłów + train → $CN_OUT"
