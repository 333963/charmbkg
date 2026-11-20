#!/bin/bash
# Generuje cd.dag: ext_<próbka>_<g> (pliki *_none.h5 próbki dzielone równo na ≤ CD_NGMAX części) → finish. NIE wysyła.
# Test: CD_DAG=, CD_LISTS=, CD_LOGDIR= przekierowują pliki poza paczkę.
set -euo pipefail
CD=$(cd "$(dirname "$0")" && pwd); export CD
source "$CD/cd_common.sh"
: "${CD_LISTS:=$CD/lists}"
DAG=${CD_DAG:-$CD/cd.dag}
mkdir -p "$CD_LOGDIR" "$CD_LISTS"
S=$CD/cd.sub
NG=""; BODY=""; EXT=()
for s in $SAMPLES; do
  D="$CD_LISTS/$s"; mkdir -p "$D"; rm -f "$D"/g[0-9][0-9][0-9]
  cd_files "$s" > "$D/all"
  n=$(wc -l < "$D/all"); [ "$n" -gt 0 ] || { echo "STOP: brak $H5D23/$s/*_none.h5" >&2; exit 1; }
  gm=$(cd_kv "$CD_NGMAX" "$s" 4); per=$(( (n + gm - 1) / gm ))
  split -d -a 3 -l "$per" "$D/all" "$D/g"
  ng=$(ls "$D"/g[0-9][0-9][0-9] | wc -l); NG="$NG $s:$ng"
  for ((i = 0; i < ng; i++)); do
    l=$(printf '%s/g%03d' "$D" "$i"); node="ext_${s}_$i"; EXT+=("$node")
    BODY+="JOB $node $S"$'\n'"VARS $node node=\"$node\" args=\"ext $s $l $i\" cpus=\"$CD_NJ\" mem=\"8GB\" disk=\"6GB\" flavour=\"workday\" __V__"$'\n'"RETRY $node 2"$'\n'
  done
  echo "  $s: $n plików → $ng części po ≤ $per"
done
PENV="CD_OUT=$CD_OUT H5D23=$H5D23 DP_MODELS=$DP_MODELS DP_THR=$DP_THR CD_PT_LO=$CD_PT_LO CD_PT_HI=$CD_PT_HI CD_LXY=$CD_LXY CD_SC_FULL=$CD_SC_FULL CD_WIN=$CD_WIN CD_PS_WIN=$CD_PS_WIN CD_PS_OUT=$CD_PS_OUT CD_MAX_GB=$CD_MAX_GB CD_HARD_GB=$CD_HARD_GB CD_NJ=$CD_NJ"
PENV+=" CD_SAMPLES_CSV=$(echo $SAMPLES | tr ' ' ',') CD_NG_CSV=$(echo $NG | tr ' ' ',')"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  printf '%s' "${BODY//__V__/$V}"
  echo "JOB finish $S"; echo "VARS finish node=\"finish\" args=\"finish\" cpus=\"1\" mem=\"6GB\" disk=\"6GB\" flavour=\"longlunch\" $V"; echo "RETRY finish 1"
  echo "PARENT ${EXT[*]} CHILD finish"
} > "$DAG"
echo "cd.dag ($DAG): ${#EXT[@]} węzłów ekstrakcji + finish → $CD_OUT"
