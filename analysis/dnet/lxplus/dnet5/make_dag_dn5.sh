#!/bin/bash
# Generuje dn5.dag (NIE wysyła):  fold_<blok>_<k> (k = 0…K−1, równolegle) → merge_<blok> ; wszystkie merge → ana → DONE ; eqbench (miękki, bez dzieci).
# Wszystkie fold/merge są miękkie (błąd → znacznik fail/ na EOS, DAG idzie dalej; ana liczy to, co jest). Sonda GPU jest poza DAG (submit_dn5.sh).
set -euo pipefail
DN5=$(cd "$(dirname "$0")" && pwd); export DN5
source "$DN5/dn5_common.sh"
mkdir -p "$Q/logs/dnet5" 2>/dev/null || true
S=$DN5/dn5.sub
PENV="DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT DN4_OUT=$DN4_OUT DN5_OUT=$DN5_OUT DN5_EPOCHS=$DN5_EPOCHS DN5_CAP_SIG=$DN5_CAP_SIG DN5_CAP_BKG=$DN5_CAP_BKG DN5_CPUS=$DN5_CPUS DN5_ANA_CPUS=$DN5_ANA_CPUS DN_FOLDS=$DN_FOLDS DN_LAM=$DN_LAM DN_NBOOT=$DN_NBOOT DN_FAST=$DN_FAST DN_REDO=${DN_REDO:-0} DN5_KEEP_B64=$(printf '%s' "$DN5_KEEP" | base64 -w0)"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
echo "$PENV" > "$DN5/dn5.penv"
V="penv=\"$PENV\""
DAG=${DN5_DAG:-$DN5/dn5.dag}
node() {  # node <nazwa> <args> <cpus> <mem> <disk> <flavour> <soft> <retry>
  echo "JOB $1 $S"; echo "VARS $1 node=\"$1\" args=\"$2\" cpus=\"$3\" mem=\"$4\" disk=\"$5\" flavour=\"$6\" soft=\"$7\" $V"; echo "RETRY $1 $8"
}
MG=(); NF=0
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for b in $DN5_BLOCKS; do
    if [ "$DN5_TRAIN" = gpu ]; then
      echo "JOB gtrain_$b $DN5/dn5_gtrain.sub"; echo "VARS gtrain_$b node=\"gtrain_$b\" args=\"gtrain $b\" wait=\"$DN5_GPU_TRAIN_WAIT\" soft=\"1\" $V"; echo "RETRY gtrain_$b 0"
      MG+=("gtrain_$b"); continue
    fi
    F=()
    for ((k = 0; k < DN_FOLDS; k++)); do x="fold_${b}_$k"; F+=("$x"); NF=$((NF + 1)); node "$x" "fold $b $k" "$DN5_CPUS" "$DN5_FOLD_MEM" 6GB workday 1 2; done
    node "merge_$b" "merge $b" "$DN5_CPUS" "$DN5_MERGE_MEM" 8GB workday 1 1
    echo "PARENT ${F[*]} CHILD merge_$b"
    MG+=("merge_$b")
  done
  node ana "ana" "$DN5_ANA_CPUS" 12GB 8GB workday 0 1
  echo "PARENT ${MG[*]} CHILD ana"
  [ "$DN5_EQBENCH" = 1 ] && [ ! -s "$DN5_OUT/test/bench.json" ] && node eqbench "eqbench" 8 8GB 4GB longlunch 1 0
} > "$DAG"
echo "dn5.dag ($DAG): tryb $DN5_TRAIN, foldy CPU $NF + łączenie ${#MG[@]} + ana + eqbench $DN5_EQBENCH = $(grep -c '^JOB' "$DAG") węzłów → $DN5_OUT (rdzeni/fold $DN5_CPUS, maxjobs $DN5_MAXJOBS)"
