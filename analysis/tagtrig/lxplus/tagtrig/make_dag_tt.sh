#!/bin/bash
# Generuje tt.dag: (próbka × część) węzłów bits + lumi (niekrytyczny: SCRIPT POST /bin/true → zawsze „sukces”) → finish. NIE wysyła.
set -euo pipefail
TT=$(cd "$(dirname "$0")" && pwd); export TT
source "$TT/tt_common.sh"
mkdir -p "$Q/logs/tagtrig"
S=$TT/tt.sub
PENV="TT_OUT=$TT_OUT RAW20=$RAW20 PTPROD=$PTPROD TT_SAMPLES=$(echo $TT_SAMPLES | tr ' ' ,) TT_NPARTS=$TT_NPARTS TT_NFILES=$(echo $TT_NFILES | tr ' ' ,)"
PENV+=" TT_STEP=$TT_STEP TT_MAX_MB=$TT_MAX_MB TT_LUMI=$TT_LUMI TT_ASETUP=$TT_ASETUP TT_LUMI_TIMEOUT=$TT_LUMI_TIMEOUT"
[ -n "$TT_LUMITAG" ] && PENV+=" TT_LUMITAG=$TT_LUMITAG"
[ -n "$TT_GRL_XML" ] && PENV+=" TT_GRL_XML=$TT_GRL_XML"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
for v in "$TT_OUT" "$RAW20" "$TT_ASETUP" "$TT_LUMITAG" "$TT_GRL_XML" "$TT_STEP"; do case "$v" in *" "*) echo "STOP: spacja w zmiennej biegu ($v)"; exit 1 ;; esac; done
DAG=${TT_DAG:-$TT/tt.dag}
V="penv=\"$PENV\""
nodes=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $TT_SAMPLES; do
    fl=longlunch; case $s in data_428770) fl=workday ;; esac
    for ((k = 0; k < TT_NPARTS; k++)); do
      j=tt_${s}_$k; nodes+=("$j")
      echo "JOB $j $S"; echo "VARS $j step=\"$s\" part=\"$k\" cpus=\"1\" mem=\"4GB\" disk=\"$TT_DISK\" flavour=\"$fl\" $V"; echo "RETRY $j 2"
    done
  done
  if [ "$TT_LUMI" = 1 ]; then
    nodes+=(tt_lumi)
    echo "JOB tt_lumi $S"; echo "VARS tt_lumi step=\"lumi\" part=\"0\" cpus=\"1\" mem=\"4GB\" disk=\"2GB\" flavour=\"longlunch\" $V"
    echo "SCRIPT POST tt_lumi /bin/true"
  fi
  echo "JOB tt_finish $S"; echo "VARS tt_finish step=\"finish\" part=\"0\" cpus=\"1\" mem=\"4GB\" disk=\"2GB\" flavour=\"longlunch\" $V"; echo "RETRY tt_finish 1"
  echo "PARENT ${nodes[*]} CHILD tt_finish"
} > "$DAG"
echo "tt.dag ($DAG): ${#nodes[@]} węzłów (bits$([ "$TT_LUMI" = 1 ] && echo ' + lumi')) + finish → $TT_OUT"
