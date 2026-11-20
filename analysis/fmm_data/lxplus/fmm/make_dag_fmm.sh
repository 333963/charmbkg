#!/bin/bash
# Generuje fmm.dag: ext_<próbka>_<g> (FMM_GROUP plików h5_v20 na węzeł) → merge. NIE wysyła. Niczego nie kasuje
# (listy w osobnym katalogu na (próbka, rozmiar grupy, liczba plików) — nadpisywane tą samą treścią).
set -euo pipefail
FM=$(cd "$(dirname "$0")" && pwd); export FM
source "$FM/fmm_common.sh"
: "${FMM_LISTS:=$FM/lists}"
mkdir -p "$Q/logs/fmm" 2>/dev/null || true
S=$FM/fmm.sub
c() { echo "$1" | tr ' ' ,; }
PENV="FMM_OUT=$FMM_OUT H5V20=$H5V20 NOMINAL=$NOMINAL FMM_LXY=$FMM_LXY FMM_SAMPLES=$(c "$FMM_SAMPLES") FMM_NFILES=$(c "$FMM_NFILES") FMM_GROUP=$(c "$FMM_GROUP") FMM_FRAC=$(c "$FMM_FRAC")"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${FMM_DAG:-$FM/fmm.dag}
EXT=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for s in $FMM_SAMPLES; do
    g=$(fm_kv "$FMM_GROUP" "$s" 10); nf=$(fm_kv "$FMM_NFILES" "$s"); ng=$(fm_ngroups "$s")
    D="$FMM_LISTS/${s}_g${g}_n${nf}"; mkdir -p "$D"
    ls "$H5V20/$s/"*_"$FMM_PRESEL".h5 2>/dev/null | grep -v '\.part$' | sort > "$D/all" || true
    [ "$(wc -l < "$D/all")" = "$nf" ] || { echo "STOP: $s — plików $(wc -l < "$D/all") ≠ FMM_NFILES $nf" >&2; exit 1; }
    split -d -a 3 -l "$g" "$D/all" "$D/g"
    for ((k = 0; k < ng; k++)); do
      l=$(printf '%s/g%03d' "$D" "$k"); [ -s "$l" ] || { echo "STOP: brak listy $l" >&2; exit 1; }
      n="ext_${s}_$k"; EXT+=("$n")
      mem=6GB; case $s in data_428770) mem=10GB ;; esac
      echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"ext $s $l $k\" cpus=\"2\" mem=\"$mem\" disk=\"4GB\" flavour=\"longlunch\" $V"; echo "RETRY $n 2"
    done
  done
  echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"16GB\" disk=\"10GB\" flavour=\"longlunch\" $V"; echo "RETRY merge 1"
  echo "PARENT ${EXT[*]} CHILD merge"
} > "$DAG"
echo "fmm.dag ($DAG): ${#EXT[@]} węzłów ekstrakcji + merge → $FMM_OUT"
