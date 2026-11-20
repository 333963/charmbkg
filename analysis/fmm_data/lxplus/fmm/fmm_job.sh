#!/bin/bash
# Węzeł DAG FMM. Użycie:
#   fmm_job.sh ext <próbka> <lista> <grupa>   fmm_extract na h5_v20 → $FMM_OUT/parts/fmm_<próbka>_<grupa>.npz
#   fmm_job.sh merge                          sklejenie → $FMM_OUT/fmm_<próbka>.npz + fmm_summary.json + DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
FM=$(cd "$(dirname "$0")" && pwd); export FM
for v in FMM_SAMPLES FMM_NFILES FMM_GROUP FMM_FRAC; do [ -n "${!v:-}" ] && export "$v=${!v//,/ }"; done   # z penv (przecinki) → spacje
source "$FM/fmm_common.sh"
step=$1; s=${2:-}; lst=${3:-}; grp=${4:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/fmm_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */fmm_*) rm -rf "${WD:?}" ;; esac; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] FMM $step $s $grp host $(hostname) WD=$WD FMM_OUT=$FMM_OUT"
[ "${DRYRUN:-0}" = 1 ] || fm_env
case $step in
  ext)
    mapfile -t FILES < "$lst"
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    X=(); [ "$(fm_tag "$s")" = fd ] && X=(--truth-max 7)
    run python3 -u "$FM/fmm_extract.py" --models "$NOMINAL/models/D0_H" --tag "$(fm_tag "$s")" --frac "$(fm_kv "$FMM_FRAC" "$s" 1.0)" \
        --lxy-min "$FMM_LXY" "${X[@]}" --out "$WD/fmm_${s}_${grp}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || fm_put "$WD/fmm_${s}_${grp}.npz" "$FMM_OUT/parts/fmm_${s}_${grp}.npz" ;;
  merge)
    NG=""; for x in $FMM_SAMPLES; do NG="$NG $x:$(fm_ngroups "$x")"; done
    run python3 -u "$FM/fmm_merge.py" --parts "$FMM_OUT/parts" --out "$WD" --samples "$FMM_SAMPLES" --ngroups "$NG" --nfiles "$FMM_NFILES"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/fmm_*.npz "$WD"/fmm_summary.json; do fm_put "$f" "$FMM_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; fm_put "$WD/DONE" "$FMM_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK FMM $step $s $grp"
