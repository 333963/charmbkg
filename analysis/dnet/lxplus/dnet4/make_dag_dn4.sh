#!/bin/bash
# Generuje dn4.dag (NIE wysyła):
#   [smoke] → extdp_<próbka>_<k> (tory otoczenia D+) ─→ train_dps{1,2,3}
#           → train_{xgbD0,xgbDp,d0s1,d0s2,d0s3,d0sb}
#   wszystkie train → app_d150_<próbka>_<k> (dane obu runów + mb_nd + fd_nd, 100 % zdarzeń) → ana (domknięcie MC → R na danych; G, SF D*) → DONE
#   train_{xgbD0,xgbDp,d0s1,dps1} → app_d23_<próbka>_<k> (pT 2–3 GeV; węzły miękkie) ; app_d23* + ana → ana23 → DONE23
# Węzły miękkie (soft=1: błąd → znacznik fail/ na EOS, DAG idzie dalej): bloki z DN4_SOFT, wszystkie app_d23, ana23.
set -euo pipefail
DN4=$(cd "$(dirname "$0")" && pwd); export DN4
source "$DN4/dn4_common.sh"
mkdir -p "$Q/logs/dnet4" 2>/dev/null || true
S=$DN4/dn4.sub
PENV="DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT DN4_OUT=$DN4_OUT H5V20=$H5V20 RAWV20=$RAWV20 H5D23=$H5D23 RAWD23=$RAWD23 DN4_MR_D0=$DN4_MR_D0 DN4_EPOCHS=$DN4_EPOCHS DN4_CAP_SIG=$DN4_CAP_SIG DN4_CAP_BKG=$DN4_CAP_BKG DN4_CPUS=$DN4_CPUS DN4_APP_CPUS=$DN4_APP_CPUS DN4_NBOOT=$DN4_NBOOT DN_FOLDS=$DN_FOLDS DN_LAM=$DN_LAM DN_LXY=$DN_LXY DN_NBOOT=$DN_NBOOT DN_FAST=$DN_FAST DN_REDO=${DN_REDO:-0}"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${DN4_DAG:-$DN4/dn4.dag}
node() {  # node <nazwa> <args> <cpus> <mem> <disk> <flavour> <soft> <retry>
  echo "JOB $1 $S"; echo "VARS $1 node=\"$1\" args=\"$2\" cpus=\"$3\" mem=\"$4\" disk=\"$5\" flavour=\"$6\" soft=\"$7\" $V"; echo "RETRY $1 $8"
}
issoft() { local x; for x in $DN4_SOFT; do [ "$x" = "$1" ] && { echo 1; return; }; done; echo 0; }
EXT=(); TR=(); TRDP=(); TR23=(); A150=(); A23=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  [ "$DN4_SMOKE" = 1 ] && node smoke "smoke" 2 6GB 3GB longlunch 0 1
  need_dp=0; for b in $DN4_BLOCKS; do case $b in dps*) need_dp=1 ;; esac; done
  if [ "$need_dp" = 1 ] && [ "${DN4_SKIP_EXTDP:-0}" != 1 ]; then
    for kv in $DN4_EXTDP; do
      s=${kv%%:*}; n=${kv##*:}
      for ((k = 0; k < n; k++)); do x="extdp_${s}_$k"; EXT+=("$x"); node "$x" "extdp $s $k $n" 1 12GB 6GB workday 0 2; done
    done
  fi
  for b in $DN4_BLOCKS; do
    x="train_$b"
    case $b in
      xgb*) node "$x" "train $b" 4 20GB 6GB workday "$(issoft "$b")" 1 ;;
      *)    node "$x" "train $b" "$DN4_CPUS" 24GB 8GB tomorrow "$(issoft "$b")" 1 ;;
    esac
    TR+=("$x")
    case $b in dps*) TRDP+=("$x") ;; esac
    case $b in xgbD0|xgbDp|d0s1|dps1) TR23+=("$x") ;; esac
  done
  for kv in $DN4_APP150; do
    s=${kv%%:*}; n=${kv##*:}
    for ((k = 0; k < n; k++)); do x="app_d150_${s}_$k"; A150+=("$x"); node "$x" "app d150 $s $k $n" "$DN4_APP_CPUS" 16GB 10GB workday 0 2; done
  done
  for kv in $DN4_APP23; do
    s=${kv%%:*}; n=${kv##*:}
    for ((k = 0; k < n; k++)); do x="app_d23_${s}_$k"; A23+=("$x"); node "$x" "app d23 $s $k $n" "$DN4_APP_CPUS" 20GB 12GB workday 1 2; done
  done
  node ana "ana" 8 32GB 8GB workday 0 1
  node pis "pis" "$DN4_APP_CPUS" 16GB 4GB longlunch 1 1          # diagnostyka wolnego pionu D* (MC, inferencja; miękki)
  [ ${#A23[@]} -gt 0 ] && node ana23 "ana23" 4 32GB 8GB workday 1 1
  if [ "$DN4_SMOKE" = 1 ]; then
    first=("${TR[@]}"); [ ${#EXT[@]} -gt 0 ] && first+=("${EXT[@]}")
    echo "PARENT smoke CHILD ${first[*]}"
  fi
  [ ${#EXT[@]} -gt 0 ] && [ ${#TRDP[@]} -gt 0 ] && echo "PARENT ${EXT[*]} CHILD ${TRDP[*]}"
  echo "PARENT ${TR[*]} CHILD ${A150[*]} pis"
  echo "PARENT ${A150[*]} CHILD ana"
  if [ ${#A23[@]} -gt 0 ]; then
    echo "PARENT ${TR23[*]} CHILD ${A23[*]}"
    echo "PARENT ${A23[*]} ana CHILD ana23"
  fi
} > "$DAG"
echo "dn4.dag ($DAG): smoke $DN4_SMOKE + extdp ${#EXT[@]} + trening ${#TR[@]} + app150 ${#A150[@]} + app23 ${#A23[@]} + ana$([ ${#A23[@]} -gt 0 ] && echo ' + ana23') = $(grep -c '^JOB' "$DAG") węzłów → $DN4_OUT"
