#!/bin/bash
# Generuje dlz.dag: frm_{fd,mb,sb}_<g> → train → ext_<próbka>_<g> → merge → ana ; d0 (niezależny) ; fin po (ana, d0). NIE wysyła.
set -euo pipefail
DZ=$(cd "$(dirname "$0")" && pwd); export DZ
source "$DZ/dlz_common.sh"
: "${DLZ_LISTS:=$DZ/lists}"
mkdir -p "$Q/logs/dlz_data" 2>/dev/null || true
S=$DZ/dlz.sub
c() { echo "$1" | tr ' ' ,; }
PENV="DLZ_OUT=$DLZ_OUT H5V20=$H5V20 H5V12=$H5V12 DLZ_LXY=$DLZ_LXY DLZ_SBFRAC=$DLZ_SBFRAC DLZ_CAP=$DLZ_CAP DLZ_NJ=$DLZ_NJ DLZ_NBOOT=$DLZ_NBOOT DLZ_FAST=$DLZ_FAST DLZ_D0=$DLZ_D0 D0_FD=$D0_FD D0_MB=$D0_MB D0_DATA=$D0_DATA D0_VARIANTS=$D0_VARIANTS D0_VARIANTS_DATA=$D0_VARIANTS_DATA DLZ_DATA=$(c "$DLZ_DATA") DLZ_NFILES=$(c "$DLZ_NFILES") DLZ_GROUP=$(c "$DLZ_GROUP") DLZ_SETS=$(c "$DLZ_SETS")"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych"; exit 1 ;; esac
V="penv=\"$PENV\""
DAG=${DLZ_DAG:-$DZ/dlz.dag}
FRM=(); EXT=()
lists() {  # <próbka> <klucz grupy> → katalog list
  local s=$1 k=$2 g nf D
  g=$(dz_kv "$DLZ_GROUP" "$k" 10); nf=$(dz_kv "$DLZ_NFILES" "$s")
  D="$DLZ_LISTS/${s}_${k}_g${g}_n${nf}"; mkdir -p "$D"
  ls "$H5V20/$s/"*_"$DLZ_PRESEL".h5 2>/dev/null | grep -v '\.part$' | sort > "$D/all" || true
  [ "$(wc -l < "$D/all")" = "$nf" ] || { echo "STOP: $s — plików $(wc -l < "$D/all") ≠ DLZ_NFILES $nf" >&2; exit 1; }
  split -d -a 3 -l "$g" "$D/all" "$D/g"; echo "$D"
}
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for t in fd mb sb; do
    case $t in fd) s=fd_nd; k=fd_nd ;; mb) s=mb_nd; k=mb_nd ;; sb) s=data_428770; k=sb ;; esac
    D=$(lists "$s" "$k"); ng=$(dz_ngroups "$s" "$k")
    for ((i = 0; i < ng; i++)); do
      l=$(printf '%s/g%03d' "$D" "$i"); [ -s "$l" ] || { echo "STOP: brak listy $l" >&2; exit 1; }
      n="frm_${t}_$i"; FRM+=("$n")
      echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"frm $t $l $i\" cpus=\"1\" mem=\"8GB\" disk=\"4GB\" flavour=\"longlunch\" $V"; echo "RETRY $n 2"
    done
  done
  echo "JOB train $S"; echo "VARS train node=\"train\" args=\"train\" cpus=\"$DLZ_NJ\" mem=\"32GB\" disk=\"10GB\" flavour=\"tomorrow\" $V"; echo "RETRY train 1"
  echo "PARENT ${FRM[*]} CHILD train"
  for s in $DLZ_DATA; do
    D=$(lists "$s" "$s"); ng=$(dz_ngroups "$s")
    for ((i = 0; i < ng; i++)); do
      l=$(printf '%s/g%03d' "$D" "$i"); [ -s "$l" ] || { echo "STOP: brak listy $l" >&2; exit 1; }
      n="ext_${s}_$i"; EXT+=("$n")
      echo "JOB $n $S"; echo "VARS $n node=\"$n\" args=\"ext $s $l $i\" cpus=\"2\" mem=\"10GB\" disk=\"6GB\" flavour=\"longlunch\" $V"; echo "RETRY $n 2"
    done
  done
  echo "PARENT train CHILD ${EXT[*]}"
  echo "JOB merge $S"; echo "VARS merge node=\"merge\" args=\"merge\" cpus=\"1\" mem=\"16GB\" disk=\"10GB\" flavour=\"longlunch\" $V"; echo "RETRY merge 1"
  echo "PARENT ${EXT[*]} CHILD merge"
  echo "JOB ana $S"; echo "VARS ana node=\"ana\" args=\"ana\" cpus=\"4\" mem=\"16GB\" disk=\"10GB\" flavour=\"workday\" $V"; echo "RETRY ana 1"
  echo "PARENT merge CHILD ana"
  echo "JOB fin $S"; echo "VARS fin node=\"fin\" args=\"fin\" cpus=\"1\" mem=\"1GB\" disk=\"1GB\" flavour=\"espresso\" $V"
  if [ "$DLZ_D0" = 1 ]; then
    echo "JOB d0 $S"; echo "VARS d0 node=\"d0\" args=\"d0\" cpus=\"$DLZ_NJ\" mem=\"32GB\" disk=\"20GB\" flavour=\"workday\" $V"; echo "RETRY d0 1"
    echo "PARENT ana d0 CHILD fin"
  else
    echo "PARENT ana CHILD fin"
  fi
} > "$DAG"
echo "dlz.dag ($DAG): ${#FRM[@]} frm + train + ${#EXT[@]} ext + merge + ana + d0($DLZ_D0) + fin → $DLZ_OUT"
