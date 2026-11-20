#!/bin/bash
# Generuje emb.dag: ext_<próbka>_<k> (k < EMB_NPARTS[próbka]) → merge. NIE wysyła.
set -euo pipefail
EMB=$(cd "$(dirname "$0")" && pwd); export EMB
source "$EMB/emb_common.sh"
mkdir -p "$Q/logs/sigemb" 2>/dev/null || true
S=$EMB/emb.sub
PENV="EMB_OUT=$EMB_OUT RAWV20=$RAWV20 EMB_SAMPLES=$(echo $EMB_SAMPLES | tr ' ' ,) EMB_NPARTS=$(echo $EMB_NPARTS | tr ' ' ,) EMB_PRESCALE=$(echo $EMB_PRESCALE | tr ' ' ,)"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
DAG=${EMB_DAG:-$EMB/emb.dag}
V="penv=\"$PENV\""
nodes=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $EMB_SAMPLES; do
    np=$(emb_np "$s")
    for ((k = 0; k < np; k++)); do
      j=ext_${s}_$k; nodes+=("$j")
      echo "JOB $j $S"; echo "VARS $j node=\"$j\" args=\"ext $s $k\" cpus=\"1\" mem=\"4GB\" disk=\"$EMB_DISK\" flavour=\"longlunch\" $V"; echo "RETRY $j 2"
    done
  done
  echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"8GB\" disk=\"$EMB_DISK\" flavour=\"longlunch\" $V"; echo "RETRY merge 1"
  echo "PARENT ${nodes[*]} CHILD merge"
} > "$DAG"
echo "emb.dag ($DAG): ${#nodes[@]} węzłów + merge → $EMB_OUT"
