#!/bin/bash
# Węzeł DAG CNCLONE. Użycie:
#   cn_job.sh ext <próbka> <lista> <grupa>   cn_extract na h5_v20 → $CN_OUT/parts/cn_<próbka>_<grupa>.npz
#   cn_job.sh train                          sklejenie + trening 2 klonów (K foldów + model pełny) + score → $CN_OUT/cn_*.npz, cn_summary.json, DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
CN=$(cd "$(dirname "$0")" && pwd); export CN
[ -n "${CN_SAMPLES:-}" ] && CN_SAMPLES=${CN_SAMPLES//,/ }   # z penv (przecinki) → spacje
source "$CN/cn_common.sh"
step=$1; s=${2:-}; lst=${3:-}; grp=${4:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/cn_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */cn_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }   # tylko własny katalog mktemp
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] CNCLONE $step $s $grp host $(hostname) WD=$WD CN_OUT=$CN_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then cn_env; mkdir -p "$CN_OUT/parts" 2>/dev/null || true; fi
case $step in
  ext)
    mapfile -t FILES < "$lst"
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    fb=(); [ "$CN_FD_BKG" = 1 ] && fb=(--fd-bkg)
    run python3 -u "$CN/cn_extract.py" --tag "$(cn_tag "$s")" --all-frac "$(cn_frac "$s")" --lxy-min "$CN_LXY" "${fb[@]}" \
        --out "$WD/cn_${s}_${grp}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || cn_put "$WD/cn_${s}_${grp}.npz" "$CN_OUT/parts/cn_${s}_${grp}.npz" ;;
  train)
    mkdir -p "$WD/out"
    run python3 -u "$CN/cn_train.py" --parts "$CN_OUT/parts" --out "$WD/out" --need "$CN_SAMPLES" --folds "$CN_FOLDS" --epochs "$CN_EPOCHS" \
        --threads "$CN_CPUS" --max-sig "$CN_MAX_SIG" --max-bkg "$CN_MAX_BKG"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/out/cn_*.npz "$WD"/out/cn_summary.json; do cn_put "$f" "$CN_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; cn_put "$WD/DONE" "$CN_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK CNCLONE $step $s $grp"
