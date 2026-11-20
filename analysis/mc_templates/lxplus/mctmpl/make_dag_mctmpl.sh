#!/bin/bash
# Generuje mctmpl.dag: ext_<próbka>_<g> (MT_GROUP plików h5_v20 na węzeł) → merge. NIE wysyła. Niczego nie kasuje
# (listy w osobnym katalogu na (próbka, rozmiar grupy, liczba plików) — nadpisywane tą samą treścią).
set -euo pipefail
MT=$(cd "$(dirname "$0")" && pwd); export MT
source "$MT/mctmpl_common.sh"
: "${MT_LISTS:=$MT/lists}"
mkdir -p "$Q/logs/mctmpl" 2>/dev/null || true
S=$MT/mctmpl.sub
c() { echo "$1" | tr ' ' ,; }
PENV="MT_OUT=$MT_OUT H5V20=$H5V20 NOMINAL=$NOMINAL MT_WP=$MT_WP MT_SAMPLES=$(c "$MT_SAMPLES") MT_NFILES=$(c "$MT_NFILES") MT_GROUP=$(c "$MT_GROUP")"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${MT_DAG:-$MT/mctmpl.dag}
EXT=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $MT_SAMPLES; do
    g=$(mt_kv "$MT_GROUP" "$s" 10); nf=$(mt_kv "$MT_NFILES" "$s"); ng=$(mt_ngroups "$s")
    D="$MT_LISTS/${s}_g${g}_n${nf}"; mkdir -p "$D"
    ls "$H5V20/$s/"*_"$MT_PRESEL".h5 2>/dev/null | grep -v '\.part$' | sort > "$D/all" || true
    [ "$(wc -l < "$D/all")" = "$nf" ] || { echo "STOP: $s — plików $(wc -l < "$D/all") ≠ MT_NFILES $nf" >&2; exit 1; }
    split -d -a 3 -l "$g" "$D/all" "$D/g"
    for ((k = 0; k < ng; k++)); do
      l=$(printf '%s/g%03d' "$D" "$k"); [ -s "$l" ] || { echo "STOP: brak listy $l" >&2; exit 1; }
      n="ext_${s}_$k"; EXT+=("$n")
      mem=6GB; case $s in data_428770) mem=12GB ;; esac
      echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"ext $s $l $k\" cpus=\"2\" mem=\"$mem\" disk=\"4GB\" flavour=\"longlunch\" $V"; echo "RETRY $n 2"
    done
  done
  echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"4GB\" disk=\"4GB\" flavour=\"microcentury\" $V"; echo "RETRY merge 1"
  echo "PARENT ${EXT[*]} CHILD merge"
} > "$DAG"
echo "mctmpl.dag ($DAG): ${#EXT[@]} węzłów ekstrakcji + merge → $MT_OUT"
