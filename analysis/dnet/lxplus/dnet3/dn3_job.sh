#!/bin/bash
# Węzeł DAG DNET etap 3. Użycie:
#   dn3_job.sh smoke            dn3_train.py selftest --need-torch (dane syntetyczne, CPU): tokeny, fit3 ≡ fit przy λ_act = 0, wszystkie bloki 1 epoka
#   dn3_job.sh train <blok>     dn3_train → $DN3_OUT/{sc,meta}
#   dn3_job.sh ana              analyze.py compute (etap 3 + odniesienia z etapów 1 i 2, wiersze wspólne) → dnet3_results.json + DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
DN3=$(cd "$(dirname "$0")" && pwd); export DN3
for v in DN_SAMPLES DN_GROUPS DN_BLOCKS DN2_NPARTS DN2_FRAC DN2_BLOCKS DN3_BLOCKS; do [ -n "${!v:-}" ] && printf -v "$v" '%s' "${!v//+/ }"; done
source "$DN3/dn3_common.sh"
# konfiguracje do analizy: odniesienia z etapu 2 (tr same córki, tr e04 — ziarno 1) + wszystko z etapu 3 (analyze.py dokłada XGB H i sieci etapu 1)
: "${DN3_KEEP:=^D0\.tr\.H\+(dau|e04)\.mc\.d[0-9]+\.s1\$|^D0\.tr\.H\+e04[a-z0-9]+\.mc\.|^D0\.tr\.H\+e04\.mc\.d[0-9]+(a[SB]|rep)}"
step=$1; s=${2:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/dn3_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */dn3_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }   # tylko własny katalog mktemp
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] DNET3 $step $s host $(hostname) WD=$WD DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT DN3_OUT=$DN3_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then
  [ "$step" = train ] && export DN3_USE_GPU=$DN3_GPU
  dn3_env; mkdir -p "$DN3_OUT"/{sc,meta} 2>/dev/null || true
fi
case $step in
  smoke)
    run python3 -u "$DN3/dn3_train.py" selftest --need-torch --dir "$WD/selftest" ;;
  train)
    export DN_DONE_DIR=$DN3_OUT
    run python3 -u "$DN3/dn3_train.py" train --data "$DN_OUT/data" --env "$DN2_OUT/parts" --out "$WD/out" --block "$s" --epochs "$DN3_EPOCHS" \
        --cap-sig "$DN3_CAP_SIG" --cap-bkg "$DN3_CAP_BKG" --threads "$DN3_CPUS"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for d in sc meta; do for f in "$WD/out/$d"/*; do [ -f "$f" ] && dn_put "$f" "$DN3_OUT/$d/$(basename "$f")"; done; done
    fi ;;
  ana)
    fast=(); [ "$DN_FAST" = 1 ] && fast=(--fast)
    run python3 -u "$DN1/analyze.py" compute --data "$DN_OUT/data" --sc "$DN_OUT,$DN2_OUT,$DN3_OUT" --keep "$DN3_KEEP" --out "$WD/dnet3_results.json" \
        --nproc "${DN_ANA_CPUS:-8}" --nboot "$DN_NBOOT" "${fast[@]}"
    if [ "${DRYRUN:-0}" != 1 ]; then
      dn_put "$WD/dnet3_results.json" "$DN3_OUT/dnet3_results.json"
      date -u +%FT%TZ > "$WD/DONE"; dn_put "$WD/DONE" "$DN3_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK DNET3 $step $s"
