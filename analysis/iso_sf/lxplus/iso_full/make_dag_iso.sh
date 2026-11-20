#!/bin/bash
# Generuje iso.dag: węzeł na (próbka, część) → finish (scalenie). NIE wysyła.
set -euo pipefail
IX=$(cd "$(dirname "$0")" && pwd); export IX
source "$IX/iso_common.sh"
mkdir -p "$Q/logs/iso_full"
S=$IX/iso.sub
PENV="ISO_OUT=$ISO_OUT H5V12=$H5V12 TAG_DIR=$TAG_DIR ISO_SAMPLES=$(echo $ISO_SAMPLES | tr ' ' ,) ISO_NPARTS_DATA=$ISO_NPARTS_DATA"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
DAG=${ISO_DAG:-$IX/iso.dag}
V="penv=\"$PENV\""
nodes=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $ISO_SAMPLES; do
    n=$(ix_nparts "$s")
    for ((k = 0; k < n; k++)); do
      j=iso_${s}_$k; nodes+=("$j")
      echo "JOB $j $S"; echo "VARS $j step=\"$s\" part=\"$k\" cpus=\"1\" mem=\"4GB\" disk=\"$ISO_DISK\" flavour=\"workday\" $V"; echo "RETRY $j 1"
    done
  done
  echo "JOB iso_finish $S"; echo "VARS iso_finish step=\"finish\" part=\"0\" cpus=\"1\" mem=\"4GB\" disk=\"$ISO_DISK\" flavour=\"longlunch\" $V"
  echo "PARENT ${nodes[*]} CHILD iso_finish"
} > "$DAG"
echo "iso.dag ($DAG): ${#nodes[@]} węzłów + finish → $ISO_OUT"
