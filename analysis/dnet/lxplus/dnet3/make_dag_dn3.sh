#!/bin/bash
# Generuje dn3.dag: [smoke] → train_<blok> (GPU: dn3_gpu.sub gdy DN3_GPU=1; inaczej CPU 8 rdzeni, flavour tomorrow) → ana. NIE wysyła.
set -euo pipefail
DN3=$(cd "$(dirname "$0")" && pwd); export DN3
source "$DN3/dn3_common.sh"
mkdir -p "$Q/logs/dnet3" 2>/dev/null || true
S=$DN3/dn3.sub; TS=$S
PENV="DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT DN3_OUT=$DN3_OUT DN3_BLOCKS=${DN3_BLOCKS// /+} DN3_GPU=$DN3_GPU DN3_EPOCHS=$DN3_EPOCHS DN3_CAP_SIG=$DN3_CAP_SIG DN3_CAP_BKG=$DN3_CAP_BKG DN3_CPUS=$DN3_CPUS DN_FOLDS=$DN_FOLDS DN_LAM=$DN_LAM DN_NBOOT=$DN_NBOOT DN_FAST=$DN_FAST DN_REDO=${DN_REDO:-0}"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${DN3_DAG:-$DN3/dn3.dag}
TR=()
tfl=tomorrow; tcpu=$DN3_CPUS
[ "$DN3_GPU" = 1 ] && { tfl=workday; tcpu=4; TS=$DN3/dn3_gpu.sub; }
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  if [ "$DN3_SMOKE" = 1 ]; then
    echo "JOB smoke $S"; echo "VARS smoke node=\"smoke\" args=\"smoke\" cpus=\"2\" mem=\"6GB\" disk=\"2GB\" flavour=\"microcentury\" $V"; echo "RETRY smoke 1"
  fi
  for b in $DN3_BLOCKS; do
    n="train_$b"; TR+=("$n")
    echo "JOB $n $TS"; echo "VARS $n node=\"$n\" args=\"train $b\" cpus=\"$tcpu\" mem=\"24GB\" disk=\"8GB\" flavour=\"$tfl\" $V"; echo "RETRY $n 1"
  done
  [ "$DN3_SMOKE" = 1 ] && echo "PARENT smoke CHILD ${TR[*]}"
  echo "JOB ana $S"; echo "VARS ana node=\"ana\" args=\"ana\" cpus=\"8\" mem=\"16GB\" disk=\"4GB\" flavour=\"workday\" $V"; echo "RETRY ana 1"
  echo "PARENT ${TR[*]} CHILD ana"
} > "$DAG"
echo "dn3.dag ($DAG): smoke $DN3_SMOKE + trening ${#TR[@]} (GPU=$DN3_GPU) + ana → $DN3_OUT"
