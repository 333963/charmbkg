#!/bin/bash
# Węzeł DAG DNET etap 2. Użycie:
#   dn2_job.sh ext <próbka> <część> <nczęści>   dn2_extract (Tracks z ROOT v2.0) → $DN2_OUT/parts/dn2_<próbka>_<część>.npz
#   dn2_job.sh train <blok>                      dn2_train → $DN2_OUT/{sc,meta,onnx}
#   dn2_job.sh ana                               analyze.py compute (etap 2 + odniesienia z etapu 1, wiersze wspólne) → dnet2_results.json + DONE
set -eo pipefail
DN2=$(cd "$(dirname "$0")" && pwd); export DN2
for v in DN_SAMPLES DN_GROUPS DN_BLOCKS DN2_NPARTS DN2_FRAC DN2_BLOCKS; do [ -n "${!v:-}" ] && printf -v "$v" '%s' "${!v//+/ }"; done
source "$DN2/dn2_common.sh"
step=$1; s=${2:-}; part=${3:-0}; np_=${4:-1}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/dn2_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */dn2_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] DNET2 $step $s $part host $(hostname) WD=$WD DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then
  [ "$step" = train ] && export DN2_USE_GPU=$DN2_GPU
  dn2_env; mkdir -p "$DN2_OUT"/{parts,sc,meta,onnx} 2>/dev/null || true
fi
case $step in
  ext)
    mx=(); [ -n "${DN2_MAXFILES:-}" ] && mx=(--max-files "$DN2_MAXFILES")
    run python3 -u "$DN2/dn2_extract.py" --sample "$s" --data "$DN_OUT/data" --rawdir "$RAWV20/$s" --part "$part" --nparts "$np_" \
        --frac-a "$(dn2_kv "$DN2_FRAC" "$s" 1.0)" --out "$WD/dn2_${s}_${part}.npz" "${mx[@]}"
    [ "${DRYRUN:-0}" = 1 ] || dn_put "$WD/dn2_${s}_${part}.npz" "$DN2_OUT/parts/dn2_${s}_${part}.npz" ;;
  train)
    export DN_DONE_DIR=$DN2_OUT
    run python3 -u "$DN2/dn2_train.py" train --data "$DN_OUT/data" --env "$DN2_OUT/parts" --out "$WD/out" --block "$s" --epochs "$DN2_EPOCHS" \
        --cap-sig "$DN2_CAP_SIG" --cap-bkg "$DN2_CAP_BKG" --threads "$DN2_CPUS"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for d in sc onnx meta; do for f in "$WD/out/$d"/*; do [ -f "$f" ] && dn_put "$f" "$DN2_OUT/$d/$(basename "$f")"; done; done
    fi ;;
  ana)
    fast=(); [ "$DN_FAST" = 1 ] && fast=(--fast)
    run python3 -u "$DN1/analyze.py" compute --data "$DN_OUT/data" --sc "$DN_OUT,$DN2_OUT" --keep '^D0\.(ds|tr)\.' --out "$WD/dnet2_results.json" \
        --nproc "${DN_ANA_CPUS:-8}" --nboot "$DN_NBOOT" "${fast[@]}"
    if [ "${DRYRUN:-0}" != 1 ]; then
      dn_put "$WD/dnet2_results.json" "$DN2_OUT/dnet2_results.json"
      date -u +%FT%TZ > "$WD/DONE"; dn_put "$WD/DONE" "$DN2_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK DNET2 $step $s $part"
