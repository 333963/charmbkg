#!/bin/bash
# Generuje tj2.dag: ext_<próbka>_<k> (k < TJ2_NPARTS[próbka]) → merge. NIE wysyła.
set -euo pipefail
TJ2=$(cd "$(dirname "$0")" && pwd); export TJ2
source "$TJ2/tj2_common.sh"
mkdir -p "$Q/logs/trkjet2" 2>/dev/null || true
S=$TJ2/tj2.sub
PENV="TJ2_OUT=$TJ2_OUT H5V20=$H5V20 RAWV20=$RAWV20 BT_IN=$BT_IN TJ2_SAMPLES=$(echo $TJ2_SAMPLES | tr ' ' ,) TJ2_NPARTS=$(echo $TJ2_NPARTS | tr ' ' ,)"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
DAG=${TJ2_DAG:-$TJ2/tj2.dag}
V="penv=\"$PENV\""
nodes=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $TJ2_SAMPLES; do
    np=$(tj2_np "$s")
    for ((k = 0; k < np; k++)); do
      j=ext_${s}_$k; nodes+=("$j")
      echo "JOB $j $S"; echo "VARS $j node=\"$j\" args=\"ext $s $k\" cpus=\"1\" mem=\"4GB\" disk=\"$TJ2_DISK\" flavour=\"longlunch\" $V"; echo "RETRY $j 2"
    done
  done
  echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"8GB\" disk=\"$TJ2_DISK\" flavour=\"longlunch\" $V"; echo "RETRY merge 1"
  echo "PARENT ${nodes[*]} CHILD merge"
} > "$DAG"
echo "tj2.dag ($DAG): ${#nodes[@]} węzłów + merge → $TJ2_OUT"
