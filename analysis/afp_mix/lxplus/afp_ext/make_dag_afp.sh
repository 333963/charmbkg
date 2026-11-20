#!/bin/bash
# Generuje afp.dag: węzeł na (próbka, część) → finish (scalenie). NIE wysyła.
set -euo pipefail
AX=$(cd "$(dirname "$0")" && pwd); export AX
source "$AX/afp_common.sh"
mkdir -p "$Q/logs/afp_ext"
S=$AX/afp.sub
PENV="AFP_OUT=$AFP_OUT H5V20=$H5V20 AFP_SAMPLES=$(echo $AFP_SAMPLES | tr ' ' ,) AFP_NPARTS=$AFP_NPARTS AFP_NFILES=$(echo $AFP_NFILES | tr ' ' ,)"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
DAG=${AFP_DAG:-$AX/afp.dag}
V="penv=\"$PENV\""
nodes=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $AFP_SAMPLES; do
    for ((k = 0; k < AFP_NPARTS; k++)); do
      j=afp_${s}_$k; nodes+=("$j")
      echo "JOB $j $S"; echo "VARS $j step=\"$s\" part=\"$k\" cpus=\"1\" mem=\"4GB\" disk=\"$AFP_DISK\" flavour=\"longlunch\" $V"; echo "RETRY $j 2"
    done
  done
  echo "JOB afp_finish $S"; echo "VARS afp_finish step=\"finish\" part=\"0\" cpus=\"1\" mem=\"4GB\" disk=\"$AFP_DISK\" flavour=\"longlunch\" $V"; echo "RETRY afp_finish 1"
  echo "PARENT ${nodes[*]} CHILD afp_finish"
} > "$DAG"
echo "afp.dag ($DAG): ${#nodes[@]} węzłów + finish → $AFP_OUT"
