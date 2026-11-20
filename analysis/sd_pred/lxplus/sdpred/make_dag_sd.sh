#!/bin/bash
# make_dag_sd.sh — sd.dag: gen × (Σ njobs z models.txt) → merge. NIE wysyła. Ziarno = 10000·indeks_modelu + numer joba (stałe, niezależne od SD_ONLY).
set -euo pipefail
SD=$(cd "$(dirname "$0")" && pwd); export SD
source "$SD/sd_common.sh"
mkdir -p "$Q/logs/sdpred"
S=$SD/sd.sub; D=${SD_DAG:-$SD/sd.dag}
PENV="SD_OUT=$SD_OUT FORCE=${FORCE:-0}"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
nodes=(); n=0
{
  echo "# wygenerowane $(date -u +%FT%TZ): SD_OUT=$SD_OUT SD_ONLY=$SD_ONLY SD_NEV=$SD_NEV SD_NJOBS=$SD_NJOBS"
  while read -r idx model nj nev cm; do
    for c in ${cm//,/ }; do [ -f "$SD/models/$c.cmnd" ] || { echo "STOP: brak models/$c.cmnd ($model)" >&2; exit 1; }; done
    fl=longlunch; case $model in nd_*) fl=workday ;; esac   # ND 1 M ≈ 42 min; reszta ≤ 15 min
    for j in $(seq 1 "$nj"); do
      seed=$((idx * 10000 + j)); nd=g_${model}_$j
      echo "JOB $nd $S"; echo "VARS $nd step=\"gen\" model=\"$model\" seed=\"$seed\" nev=\"$nev\" cm=\"$cm\" flavour=\"$fl\" penv=\"$PENV\""
      nodes+=("$nd"); n=$((n+1))
    done
  done < <(sd_models)
  echo "JOB merge $S"; echo "VARS merge step=\"merge\" model=\"-\" seed=\"0\" nev=\"0\" cm=\"-\" flavour=\"longlunch\" penv=\"$PENV\""
  echo "PARENT ${nodes[*]} CHILD merge"
  echo "RETRY ALL_NODES 2"
} > "$D"
echo "sd.dag: $n × gen + merge → $SD_OUT"
