#!/bin/bash
# Generuje dn.dag: ext_<próbka>_<g> (24) → merge → train_<blok> (10) → ana. NIE wysyła.
# DN_SKIP_EXT=1: bez ekstrakcji i merge (dane już w $DN_OUT/data) — np. dołożenie bloków.
set -euo pipefail
DN=$(cd "$(dirname "$0")" && pwd); export DN
source "$DN/dn_common.sh"
: "${DN_LISTS:=$DN/lists/$(date -u +%Y%m%dT%H%M%S)}"
mkdir -p "$Q/logs/dnet1" 2>/dev/null || true; mkdir -p "$DN_LISTS"
S=$DN/dn.sub
PENV="DN_OUT=$DN_OUT H5V20=$H5V20 DN_SAMPLES=${DN_SAMPLES// /+} DN_GROUPS=${DN_GROUPS// /+} DN_BLOCKS=${DN_BLOCKS// /+} DN_FRAC_D0=$DN_FRAC_D0 DN_FRAC_DP=$DN_FRAC_DP DN_FRAC_DP_MB=$DN_FRAC_DP_MB DN_LXY=$DN_LXY DN_FOLDS=$DN_FOLDS DN_EPOCHS=$DN_EPOCHS DN_CPUS=$DN_CPUS DN_CAP_SIG=$DN_CAP_SIG DN_CAP_BKG=$DN_CAP_BKG DN_LAM=$DN_LAM DN_LAM_SCAN=$DN_LAM_SCAN DN_NBOOT=$DN_NBOOT DN_FAST=$DN_FAST"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${DN_DAG:-$DN/dn.dag}
EXT=(); TR=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  if [ "${DN_SKIP_EXT:-0}" != 1 ]; then
    for s in $DN_SAMPLES; do
      ls "$H5V20/$s/"*_"$DN_PRESEL".h5 2>/dev/null | grep -v '\.part$' > "$DN_LISTS/$s.all" || true
      [ -s "$DN_LISTS/$s.all" ] || { echo "STOP: brak plików h5_v20 dla $s" >&2; exit 1; }
      split -d -a 3 -l "$(dn_group "$s")" "$DN_LISTS/$s.all" "$DN_LISTS/${s}_"
      for l in "$DN_LISTS/${s}_"[0-9][0-9][0-9]; do
        g=$((10#${l##*_})); n="ext_${s}_$g"; EXT+=("$n")
        mem=6GB; case $s in data_428770|mb_nd) mem=10GB ;; esac
        echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"ext $s $l $g\" cpus=\"1\" mem=\"$mem\" disk=\"4GB\" flavour=\"longlunch\" $V"; echo "RETRY $n 2"
      done
    done
    echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"16GB\" disk=\"12GB\" flavour=\"longlunch\" $V"; echo "RETRY merge 1"
    echo "PARENT ${EXT[*]} CHILD merge"
  fi
  for b in $DN_BLOCKS; do
    n="train_$b"; TR+=("$n")
    echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"train $b\" cpus=\"$DN_CPUS\" mem=\"16GB\" disk=\"6GB\" flavour=\"workday\" $V"; echo "RETRY $n 1"
  done
  [ "${DN_SKIP_EXT:-0}" = 1 ] || echo "PARENT merge CHILD ${TR[*]}"
  echo "JOB ana $S"; echo "VARS ana node=\"ana\" args=\"ana\" cpus=\"8\" mem=\"16GB\" disk=\"4GB\" flavour=\"workday\" $V"; echo "RETRY ana 1"   # 8 rdzeni: ~650 fitów Δm z profilem
  echo "PARENT ${TR[*]} CHILD ana"
} > "$DAG"
echo "dn.dag ($DAG): ekstrakcja ${#EXT[@]} + merge + trening ${#TR[@]} + ana → $DN_OUT"
