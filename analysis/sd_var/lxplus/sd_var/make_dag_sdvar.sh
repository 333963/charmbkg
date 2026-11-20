#!/bin/bash
# make_dag_sdvar.sh — sdvar.dag: gen × (Σ njobs z models.txt) + xsec × (modele dv_/ref_) → merge. NIE wysyła.
# Ziarno = 1000000 + 10000·indeks_modelu + numer joba (stałe, niezależne od SD_ONLY/SD_GROUPS; rozłączne z sdpred); xsec: numer joba 0.
set -euo pipefail
SD=$(cd "$(dirname "$0")" && pwd); export SD
source "$SD/sd_common.sh"
mkdir -p "$SD_LOGS"
S=$SD/sdvar.sub; D=${SD_DAG:-$SD/sdvar.dag}
PENV="SD_OUT=$SD_OUT FORCE=${FORCE:-0}"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
nodes=(); n=0; nx=0
{
  echo "# wygenerowane $(date -u +%FT%TZ): SD_OUT=$SD_OUT SD_ONLY=$SD_ONLY SD_GROUPS=$SD_GROUPS SD_NEV=$SD_NEV SD_NJOBS=$SD_NJOBS"
  while read -r idx model nj nev cm; do
    for c in ${cm//,/ }; do [ -f "$SD/models/$c.cmnd" ] || { echo "STOP: brak models/$c.cmnd ($model)" >&2; exit 1; }; done
    fl=longlunch; case $model in nv_*) fl=workday ;; esac   # ND 1 M ≈ 42 min; reszta ≤ 20 min
    if [ "${SD_NEV:-}" != "" ] && [ "$SD_NEV" -le 100000 ]; then fl=espresso; fi
    for j in $(seq 1 "$nj"); do
      seed=$((1000000 + idx * 10000 + j)); nd=g_${model}_$j
      echo "JOB $nd $S"; echo "VARS $nd step=\"gen\" model=\"$model\" seed=\"$seed\" nev=\"$nev\" cm=\"$cm\" flavour=\"$fl\" penv=\"$PENV\""
      nodes+=("$nd"); n=$((n+1))
    done
    if sd_needs_xsec "$model"; then
      nd=x_${model}
      echo "JOB $nd $S"; echo "VARS $nd step=\"xsec\" model=\"$model\" seed=\"$((1000000 + idx * 10000))\" nev=\"$SD_XSEC_NEV\" cm=\"$cm\" flavour=\"espresso\" penv=\"$PENV\""
      nodes+=("$nd"); nx=$((nx+1))
    fi
  done < <(sd_models)
  echo "JOB merge $S"; echo "VARS merge step=\"merge\" model=\"-\" seed=\"0\" nev=\"0\" cm=\"-\" flavour=\"longlunch\" penv=\"$PENV\""
  echo "PARENT ${nodes[*]} CHILD merge"
  echo "RETRY ALL_NODES 2"
} > "$D"
echo "sdvar.dag: $n × gen + $nx × xsec + merge → $SD_OUT"
