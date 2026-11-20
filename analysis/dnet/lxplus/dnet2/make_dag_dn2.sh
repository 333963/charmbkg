#!/bin/bash
# Generuje dn2.dag: ext_<próbka>_<k> (24) → train_<blok> (9; GPU gdy DN2_GPU=1) → ana. NIE wysyła. DN2_SKIP_EXT=1: bez ekstrakcji.
set -euo pipefail
DN2=$(cd "$(dirname "$0")" && pwd); export DN2
source "$DN2/dn2_common.sh"
mkdir -p "$Q/logs/dnet2" 2>/dev/null || true
S=$DN2/dn2.sub
PENV="DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT RAWV20=$RAWV20 DN2_NPARTS=${DN2_NPARTS// /+} DN2_FRAC=${DN2_FRAC// /+} DN2_BLOCKS=${DN2_BLOCKS// /+} DN2_GPU=$DN2_GPU DN2_EPOCHS=$DN2_EPOCHS DN2_CAP_SIG=$DN2_CAP_SIG DN2_CAP_BKG=$DN2_CAP_BKG DN2_CPUS=$DN2_CPUS DN_FOLDS=$DN_FOLDS DN_LAM=$DN_LAM DN_NBOOT=$DN_NBOOT DN_FAST=$DN_FAST"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${DN2_DAG:-$DN2/dn2.dag}
EXT=(); TR=()
tfl=tomorrow; tcpu=$DN2_CPUS; [ "$DN2_GPU" = 1 ] && { tfl=workday; tcpu=4; }
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  if [ "${DN2_SKIP_EXT:-0}" != 1 ]; then
    for s in $DN_SAMPLES; do
      np_=$(dn2_kv "$DN2_NPARTS" "$s" 4)
      for ((k = 0; k < np_; k++)); do
        n="ext_${s}_$k"; EXT+=("$n")
        echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"ext $s $k $np_\" cpus=\"1\" gpus=\"0\" mem=\"12GB\" disk=\"6GB\" flavour=\"workday\" $V"; echo "RETRY $n 2"
      done
    done
  fi
  for b in $DN2_BLOCKS; do
    n="train_$b"; TR+=("$n")
    echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"train $b\" cpus=\"$tcpu\" gpus=\"$DN2_GPU\" mem=\"24GB\" disk=\"8GB\" flavour=\"$tfl\" $V"; echo "RETRY $n 1"
  done
  [ "${DN2_SKIP_EXT:-0}" = 1 ] || echo "PARENT ${EXT[*]} CHILD ${TR[*]}"
  echo "JOB ana $S"; echo "VARS ana node=\"ana\" args=\"ana\" cpus=\"8\" gpus=\"0\" mem=\"16GB\" disk=\"4GB\" flavour=\"workday\" $V"; echo "RETRY ana 1"
  echo "PARENT ${TR[*]} CHILD ana"
} > "$DAG"
echo "dn2.dag ($DAG): ekstrakcja ${#EXT[@]} + trening ${#TR[@]} (GPU=$DN2_GPU) + ana → $DN2_OUT"
