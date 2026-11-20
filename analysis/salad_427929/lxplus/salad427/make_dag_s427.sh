#!/bin/bash
# Generuje s427.dag: ext_<próbka>_<g> → study_<gatunek> (D0, Dplus równolegle) → fin. NIE wysyła.
# Listy w lists/<znacznik czasu>/ (nowy katalog przy każdym wywołaniu → bez kasowania i bez starych list w DAG-u).
set -euo pipefail
S4=$(cd "$(dirname "$0")" && pwd); export S4
source "$S4/s427_common.sh"
: "${S4_LISTS:=$S4/lists/$(date -u +%Y%m%dT%H%M%S)}"
mkdir -p "$Q/logs/salad427" 2>/dev/null || true; mkdir -p "$S4_LISTS"
S=$S4/s427.sub
PENV="S4_OUT=$S4_OUT H5V20=$H5V20 S4_SAMPLES=$(echo $S4_SAMPLES | tr ' ' ,) S4_SPECIES=$(echo $S4_SPECIES | tr ' ' ,) S4_FRAC_J12=$S4_FRAC_J12 S4_LXY=$S4_LXY S4_TASKS=${S4_TASKS//,/+} S4_CPUS=$S4_CPUS S4_FAST=$S4_FAST"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${S4_DAG:-$S4/s427.dag}
EXT=(); ST=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $S4_SAMPLES; do
    ls "$H5V20/$s/"*_"$S4_PRESEL".h5 2>/dev/null | grep -v '\.part$' > "$S4_LISTS/$s.all" || true
    [ -s "$S4_LISTS/$s.all" ] || { echo "STOP: brak plików h5_v20 dla $s" >&2; exit 1; }
    split -d -a 3 -l "$(s4_kv "$S4_GROUP" "$s")" "$S4_LISTS/$s.all" "$S4_LISTS/${s}_"
    for l in "$S4_LISTS/${s}_"[0-9][0-9][0-9]; do
      g=$((10#${l##*_})); n="ext_${s}_$g"; EXT+=("$n")
      mem=4GB; case $s in data_428770) mem=8GB ;; esac
      echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"ext $s $l $g\" cpus=\"1\" mem=\"$mem\" disk=\"3GB\" flavour=\"longlunch\" $V"; echo "RETRY $n 2"
    done
  done
  for sp in $S4_SPECIES; do
    n="study_$sp"; ST+=("$n")
    echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"study $sp\" cpus=\"$S4_CPUS\" mem=\"24GB\" disk=\"8GB\" flavour=\"tomorrow\" $V"; echo "RETRY $n 1"
  done
  echo "JOB fin $S"; echo "VARS fin node=\"fin\" args=\"fin\" cpus=\"1\" mem=\"1GB\" disk=\"1GB\" flavour=\"espresso\" $V"
  echo "PARENT ${EXT[*]} CHILD ${ST[*]}"
  echo "PARENT ${ST[*]} CHILD fin"
} > "$DAG"
echo "s427.dag ($DAG): ekstrakcja ${#EXT[@]} węzłów + study ${#ST[@]} + fin → $S4_OUT"
