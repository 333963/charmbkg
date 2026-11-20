#!/bin/bash
# Generuje xa.dag (NIE wysyla): ext_<probka>_<k> -> fin.
set -euo pipefail
XA=$(cd "$(dirname "$0")" && pwd)
Q=$(cd "$XA/.." && pwd)
mkdir -p "$Q/logs/xsatlas"
NODES=${XA_NODES:-data_427929:8 mb_nd:6 fd_nd:3}
E=()
{
  echo "# wygenerowane $(date -u +%FT%TZ)"
  for x in $NODES; do
    IFS=: read -r s n <<< "$x"
    for ((k = 0; k < n; k++)); do
      nd="ext_${s}_$k"; E+=("$nd")
      echo "JOB $nd $XA/xa.sub"; echo "VARS $nd node=\"$nd\" args=\"ext $s $k $n\" mem=\"4GB\" flavour=\"longlunch\""; echo "RETRY $nd 2"
    done
  done
  echo "JOB fin $XA/xa.sub"; echo "VARS fin node=\"fin\" args=\"fin\" mem=\"1GB\" flavour=\"espresso\""; echo "RETRY fin 1"
  echo "PARENT ${E[*]} CHILD fin"
} > "$XA/xa.dag"
echo "xa.dag: $(grep -c '^JOB' "$XA/xa.dag") wezlow"
