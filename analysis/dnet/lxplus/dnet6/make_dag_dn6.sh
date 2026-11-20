#!/bin/bash
# Generuje dn6.dag (NIE wysyła). Tylko bloki, których tag NIE jest gotowy na EOS (meta + sc + models/meta) — gotowe są pomijane (brak podwójnych wyników).
#   DN6_TRAIN=gpu:  gtrain_<blok> (1 job GPU = 5 foldów + łączenie; sam się usuwa po DN6_GPU_WAIT s w kolejce)
#   DN6_TRAIN=cpu:  fold_<blok>_<k> (k = 0…K−1, równolegle; 4 rdzenie) → merge_<blok>
#   wszystkie treningi → ana → DONE ;  ana → app_<próbka>_<k> (miękkie, bramka P1) → anar (miękki) → DONE_R     [DN6_APP=1]
set -euo pipefail
DN6=$(cd "$(dirname "$0")" && pwd); export DN6
source "$DN6/dn6_common.sh"
mkdir -p "$Q/logs/dnet6" 2>/dev/null || true
S=$DN6/dn6.sub
PENV="DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT DN4_OUT=$DN4_OUT DN6_OUT=$DN6_OUT H5V20=$H5V20 RAWV20=$RAWV20 DN4_MR_D0=$DN4_MR_D0 DN4_NBOOT=$DN4_NBOOT DN6_DM_MAX=$DN6_DM_MAX DN6_SIGN=$DN6_SIGN DN6_EPOCHS=$DN6_EPOCHS DN6_CAP_SIG=$DN6_CAP_SIG DN6_CAP_BKG=$DN6_CAP_BKG DN6_CPUS=$DN6_CPUS DN6_BLOCKS=${DN6_BLOCKS// /+} DN_FOLDS=$DN_FOLDS DN_LAM=$DN_LAM DN_LXY=$DN_LXY DN_NBOOT=$DN_NBOOT DN_FAST=$DN_FAST DN_REDO=${DN_REDO:-0} DN6_APP_NOGATE=${DN6_APP_NOGATE:-0} DN6_APP_REF=$DN6_APP_REF DN6_APP_TAG_B64=$(printf '%s' "$DN6_APP_TAG" | base64 -w0) DN6_KEEP_B64=$(printf '%s' "$DN6_KEEP" | base64 -w0)"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
echo "$PENV" > "$DN6/dn6.penv"
V="penv=\"$PENV\""
DAG=${DN6_DAG:-$DN6/dn6.dag}
node() {  # node <nazwa> <args> <cpus> <mem> <disk> <flavour> <soft> <retry>
  echo "JOB $1 $S"; echo "VARS $1 node=\"$1\" args=\"$2\" cpus=\"$3\" mem=\"$4\" disk=\"$5\" flavour=\"$6\" soft=\"$7\" $V"; echo "RETRY $1 $8"
}
MG=(); NF=0; NG=0; SKIP=(); APP=()
{
  echo "# wygenerowane $(date -u +%FT%TZ) tryb $DN6_TRAIN: $PENV"
  for b in $DN6_BLOCKS; do
    tag=$(dn6_tag "$b"); [ -n "$tag" ] || { echo "nieznany blok $b" >&2; exit 2; }
    if dn6_ready "$tag" && [ "${DN_REDO:-0}" != 1 ]; then SKIP+=("$b"); continue; fi
    if [ "$DN6_TRAIN" = gpu ]; then
      echo "JOB gtrain_$b $DN6/dn6_gtrain.sub"; echo "VARS gtrain_$b node=\"gtrain_$b\" args=\"gtrain $b\" wait=\"$DN6_GPU_WAIT\" $V"; echo "RETRY gtrain_$b 0"
      MG+=("gtrain_$b"); NG=$((NG + 1)); continue
    fi
    F=()
    for ((k = 0; k < DN_FOLDS; k++)); do x="fold_${b}_$k"; F+=("$x"); NF=$((NF + 1)); node "$x" "fold $b $k" "$DN6_CPUS" "$DN6_FOLD_MEM" 6GB workday 0 2; done
    node "merge_$b" "merge $b" "$DN6_CPUS" "$DN6_MERGE_MEM" 8GB longlunch 0 1
    echo "PARENT ${F[*]} CHILD merge_$b"
    MG+=("merge_$b")
  done
  node ana "ana" "$DN6_CPUS" "$DN6_ANA_MEM" 8GB workday 0 1
  [ ${#MG[@]} -gt 0 ] && echo "PARENT ${MG[*]} CHILD ana"
  if [ "$DN6_APP" = 1 ]; then
    for kv in $DN6_APP150; do
      s=${kv%%:*}; n=${kv##*:}
      for ((k = 0; k < n; k++)); do x="app_${s}_$k"; APP+=("$x"); node "$x" "app $s $k $n" "$DN6_CPUS" 16GB 10GB workday 1 2; done
    done
    node anar "anar" "$DN6_CPUS" 16GB 8GB workday 1 1
    echo "PARENT ana CHILD ${APP[*]}"
    echo "PARENT ${APP[*]} CHILD anar"
  fi
} > "$DAG"
echo "dn6.dag ($DAG): tryb $DN6_TRAIN — GPU $NG, foldy CPU $NF, łączenia $(( ${#MG[@]} - NG )), ana, app ${#APP[@]}$([ "$DN6_APP" = 1 ] && echo ' + anar') = $(grep -c '^JOB' "$DAG") węzłów; tagi już gotowe (pominięte): ${SKIP[*]:-brak} → $DN6_OUT"
