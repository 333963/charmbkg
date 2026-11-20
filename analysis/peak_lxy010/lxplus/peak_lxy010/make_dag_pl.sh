#!/bin/bash
# Generuje pl.dag: jeden węzeł ekstrakcji na próbkę (równolegle) → finish. NIE wysyła.
set -euo pipefail
PL=$(cd "$(dirname "$0")" && pwd); export PL
source "$PL/pl_common.sh"
mkdir -p "$Q/logs/peaklxy"
S=$PL/pl.sub
# środowisko węzłów (pl.sub: getenv=False, environment="$(penv)"). SAMPLES ma spacje → w penv jako lista po przecinku
# PL_SAMPLES_CSV; pl_common.sh odtwarza z niej SAMPLES (finish sprawdza komplet tego samego zestawu).
PENV="PL_OUT=$PL_OUT H5V20=$H5V20 NOMINAL=$NOMINAL DATA_FRAC=$DATA_FRAC MAX_GB=$MAX_GB SCORE_MODE=$SCORE_MODE"
[ -n "${PL_EXTRA:-}" ] && PENV+=" PL_EXTRA=$PL_EXTRA"
PENV+=" PL_SAMPLES_CSV=$(echo $SAMPLES | tr ' ' ',')"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
{
  echo "# wygenerowane $(date -u +%FT%TZ): PL_OUT=$PL_OUT H5V20=$H5V20 NOMINAL=$NOMINAL SAMPLES=$SAMPLES DATA_FRAC=$DATA_FRAC SCORE_MODE=$SCORE_MODE PL_EXTRA=${PL_EXTRA:-}"
  for s in $SAMPLES; do
    mem=6GB; fl=workday
    case $s in data_*) mem=12GB ;; esac
    echo "JOB pl_$s $S"; echo "VARS pl_$s step=\"$s\" cpus=\"2\" mem=\"$mem\" flavour=\"$fl\" $V"; echo "RETRY pl_$s 1"
  done
  echo "JOB pl_finish $S"; echo "VARS pl_finish step=\"finish\" cpus=\"1\" mem=\"4GB\" flavour=\"espresso\" $V"
  echo "PARENT $(for s in $SAMPLES; do printf 'pl_%s ' "$s"; done)CHILD pl_finish"
} > "$PL/pl.dag"
echo "pl.dag: $(echo "$SAMPLES" | wc -w) węzły ekstrakcji + finish → $PL_OUT"
