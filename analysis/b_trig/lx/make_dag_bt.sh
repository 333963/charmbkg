#!/bin/bash
# Generuje bt.dag: trg_<próbka> (ścieżka A, ROOT) + ext_<próbka>_<g> (ścieżka B, h5_v20, BT_GROUP plików) → merge. NIE wysyła.
# BT_PATH=A|B|AB (domyślnie AB; A = tylko tabele zdarzeń z ROOT — nie potrzebuje h5_v20).
set -euo pipefail
BT=$(cd "$(dirname "$0")" && pwd); export BT
source "$BT/bt_common.sh"
: "${BT_PATH:=AB}"
: "${BT_LISTS:=$BT/lists/g$BT_GROUP}"                        # osobny katalog na BT_GROUP (bez kasowania starych list)
mkdir -p "$Q/logs/b_trig" 2>/dev/null || true; mkdir -p "$BT_LISTS"
S=$BT/bt.sub
PENV="BT_OUT=$BT_OUT H5V20=$H5V20 RAWV20=$RAWV20 NOMINAL=$NOMINAL BT_FRAC_DATA=$BT_FRAC_DATA BT_FRAC_MB=$BT_FRAC_MB BT_LXY=$BT_LXY"
case "$PENV" in *\'*|*\"*|*[[:space:]][[:space:]]*) echo "STOP: cudzysłów w zmiennych"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${BT_DAG:-$BT/bt.dag}
EXT=(); TRG=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): BT_PATH=$BT_PATH $PENV"
  if [[ $BT_PATH == *A* ]]; then
    for s in $BT_TRG; do
      TRG+=("trg_$s")
      echo "JOB trg_$s $S"; echo "VARS trg_$s node=\"trg_$s\" args=\"trg $s\" cpus=\"1\" mem=\"8GB\" disk=\"2GB\" flavour=\"workday\" $V"; echo "RETRY trg_$s 1"
    done
  fi
  if [[ $BT_PATH == *B* ]]; then
    for s in $BT_SAMPLES; do
      ls "$H5V20/$s/"*_"$BT_PRESEL".h5 2>/dev/null | grep -v '\.part$' > "$BT_LISTS/$s.all" || true
      [ -s "$BT_LISTS/$s.all" ] || { echo "STOP: brak plików h5_v20 dla $s" >&2; exit 1; }
      split -d -a 3 -l "$BT_GROUP" "$BT_LISTS/$s.all" "$BT_LISTS/${s}_"
      for l in "$BT_LISTS/${s}_"[0-9][0-9][0-9]; do
        g=$((10#${l##*_})); n="ext_${s}_$g"; EXT+=("$n")
        mem=6GB; case $s in data_*) mem=10GB ;; esac
        echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"ext $s $l $g\" cpus=\"2\" mem=\"$mem\" disk=\"4GB\" flavour=\"longlunch\" $V"; echo "RETRY $n 1"
      done
    done
    echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"16GB\" disk=\"8GB\" flavour=\"longlunch\" $V"
    echo "PARENT ${EXT[*]} ${TRG[*]} CHILD merge"            # DONE (merge) dopiero po trg_* → fetch ma komplet
  fi
} > "$DAG"
echo "bt.dag ($DAG): BT_PATH=$BT_PATH, ekstrakcja ${#EXT[@]} węzłów → $BT_OUT"
