#!/bin/bash
# Generuje dstws.dag: ext_<próbka>_<k> (k < DW_NPARTS[próbka]) → merge. NIE wysyła.
set -euo pipefail
DW=$(cd "$(dirname "$0")" && pwd); export DW
source "$DW/dstws_common.sh"
mkdir -p "$Q/logs/dstws" 2>/dev/null || true
S=$DW/dstws.sub
PENV="DW_OUT=$DW_OUT H5V20=$H5V20 RAWV20=$RAWV20 NOMINAL=$NOMINAL DW_SRC=$DW_SRC DW_LXY=$DW_LXY DW_SAMPLES=$(echo $DW_SAMPLES | tr ' ' ,) DW_NPARTS=$(echo $DW_NPARTS | tr ' ' ,)"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
DAG=${DW_DAG:-$DW/dstws.dag}
V="penv=\"$PENV\""
nodes=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $DW_SAMPLES; do
    np=$(dw_np "$s")
    for ((k = 0; k < np; k++)); do
      j=ext_${s}_$k; nodes+=("$j")
      echo "JOB $j $S"; echo "VARS $j node=\"$j\" args=\"ext $s $k\" cpus=\"1\" mem=\"$DW_MEM\" disk=\"$DW_DISK\" flavour=\"longlunch\" $V"; echo "RETRY $j 2"
    done
  done
  echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"12GB\" disk=\"6GB\" flavour=\"longlunch\" $V"; echo "RETRY merge 1"
  echo "PARENT ${nodes[*]} CHILD merge"
} > "$DAG"
echo "dstws.dag ($DAG): ${#nodes[@]} węzłów + merge → $DW_OUT"
