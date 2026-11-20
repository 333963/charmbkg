#!/bin/bash
# Generuje tj.dag: ext_<próbka>_<k> (k < TJ_NPARTS[próbka]) → merge. NIE wysyła.
set -euo pipefail
TJ=$(cd "$(dirname "$0")" && pwd); export TJ
source "$TJ/tj_common.sh"
mkdir -p "$Q/logs/trkjet" 2>/dev/null || true
S=$TJ/tj.sub
PENV="TJ_OUT=$TJ_OUT H5V20=$H5V20 RAWV20=$RAWV20 BT_IN=$BT_IN TJ_SAMPLES=$(echo $TJ_SAMPLES | tr ' ' ,) TJ_NPARTS=$(echo $TJ_NPARTS | tr ' ' ,)"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
DAG=${TJ_DAG:-$TJ/tj.dag}
V="penv=\"$PENV\""
nodes=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $TJ_SAMPLES; do
    np=$(tj_np "$s")
    for ((k = 0; k < np; k++)); do
      j=ext_${s}_$k; nodes+=("$j")
      echo "JOB $j $S"; echo "VARS $j node=\"$j\" args=\"ext $s $k\" cpus=\"1\" mem=\"4GB\" disk=\"$TJ_DISK\" flavour=\"longlunch\" $V"; echo "RETRY $j 2"
    done
  done
  echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"8GB\" disk=\"$TJ_DISK\" flavour=\"longlunch\" $V"; echo "RETRY merge 1"
  echo "PARENT ${nodes[*]} CHILD merge"
} > "$DAG"
echo "tj.dag ($DAG): ${#nodes[@]} węzłów + merge → $TJ_OUT"
