#!/bin/bash
# Węzeł DAG DNET etap 1. Użycie:
#   dn_job.sh ext <próbka> <lista> <grupa>   dn_extract na h5_v20 → $DN_OUT/parts/dn_<próbka>_<grupa>.npz
#   dn_job.sh merge                          części → $DN_OUT/data/dn_<próbka>.npz
#   dn_job.sh train <blok>                   konfiguracje bloku (dn_feat.BLOCKS) → $DN_OUT/{sc,meta,onnx}/<tag>.*
#   dn_job.sh ana                            analyze.py compute → $DN_OUT/dnet_results.json + DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
DN=$(cd "$(dirname "$0")" && pwd); export DN
for v in DN_SAMPLES DN_GROUPS DN_BLOCKS; do [ -n "${!v:-}" ] && printf -v "$v" '%s' "${!v//+/ }"; done   # z penv („+”) → spacje
source "$DN/dn_common.sh"
step=$1; s=${2:-}; lst=${3:-}; grp=${4:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/dn_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */dn_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }   # tylko własny katalog mktemp
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] DNET1 $step $s $grp host $(hostname) WD=$WD DN_OUT=$DN_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then dn_env; mkdir -p "$DN_OUT"/{parts,data,sc,meta,onnx} 2>/dev/null || true; fi
case $step in
  ext)
    mapfile -t FILES < "$lst"
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    fd0=1.0; fdp=1.0
    case $s in data_*) fd0=$DN_FRAC_D0; fdp=$DN_FRAC_DP ;; mb_*) fdp=$DN_FRAC_DP_MB ;; esac
    run python3 -u "$DN/dn_extract.py" --tag "$(dn_tag "$s")" --frac-d0 "$fd0" --frac-dp "$fdp" --lxy-min "$DN_LXY" \
        --out "$WD/dn_${s}_${grp}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || dn_put "$WD/dn_${s}_${grp}.npz" "$DN_OUT/parts/dn_${s}_${grp}.npz" ;;
  merge)
    run python3 -u "$DN/dn_train.py" merge --parts "$DN_OUT/parts" --out "$WD/data"
    if [ "${DRYRUN:-0}" != 1 ]; then for f in "$WD"/data/dn_*.npz; do dn_put "$f" "$DN_OUT/data/$(basename "$f")"; done; fi ;;
  train)
    export DN_DONE_DIR=$DN_OUT
    run python3 -u "$DN/dn_train.py" train --data "$DN_OUT/data" --out "$WD/out" --block "$s"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for d in sc onnx meta; do for f in "$WD/out/$d"/*; do [ -f "$f" ] && dn_put "$f" "$DN_OUT/$d/$(basename "$f")"; done; done   # meta na końcu = znacznik „gotowe”
    fi ;;
  ana)
    fast=(); [ "$DN_FAST" = 1 ] && fast=(--fast)
    run python3 -u "$DN/analyze.py" compute --data "$DN_OUT/data" --sc "$DN_OUT" --out "$WD/dnet_results.json" --nproc "${DN_ANA_CPUS:-8}" --nboot "$DN_NBOOT" "${fast[@]}"
    if [ "${DRYRUN:-0}" != 1 ]; then
      dn_put "$WD/dnet_results.json" "$DN_OUT/dnet_results.json"
      date -u +%FT%TZ > "$WD/DONE"; dn_put "$WD/DONE" "$DN_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK DNET1 $step $s $grp"
