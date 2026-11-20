#!/bin/bash
# Węzeł CAND_D23. Użycie:
#   cd_job.sh ext <próbka> <lista plików h5> <g>   → $CD_OUT/parts/{dp_<s>_<g>.parquet, hist_<s>_<g>.npz, part_<s>_<g>.json}
#   cd_job.sh finish                               → $CD_OUT/{dp_<s>.parquet, hist_<s>.npz, summary.json, DONE}
# Scratch mktemp + trap (kasowany); wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
CD=$(cd "$(dirname "$0")" && pwd); export CD
source "$CD/cd_common.sh"
step=$1; s=${2:-}; lst=${3:-}; grp=${4:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/cd23_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */cd23_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] CAND_D23 $step $s $grp host $(hostname) WD=$WD CD_OUT=$CD_OUT H5D23=$H5D23 pT=[$CD_PT_LO,$CD_PT_HI)"
[ "${DRYRUN:-0}" = 1 ] || cd_env
case $step in
  ext)
    mapfile -t FILES < "$lst"; [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    # shellcheck disable=SC2046
    run python3 -u "$CD/cd_extract.py" --sample "$s" --tag "$grp" --out-dir "$WD/o" $(cd_xargs) "${FILES[@]}"
    if [ "${DRYRUN:-0}" != 1 ]; then
      [ -f "$WD/o/dp_${s}_${grp}.parquet" ] && cd_put "$WD/o/dp_${s}_${grp}.parquet" "$CD_OUT/parts/dp_${s}_${grp}.parquet"
      cd_put "$WD/o/hist_${s}_${grp}.npz" "$CD_OUT/parts/hist_${s}_${grp}.npz"
      cd_put "$WD/o/part_${s}_${grp}.json" "$CD_OUT/parts/part_${s}_${grp}.json"      # ostatni = znacznik kompletu części
    fi ;;
  finish)
    [ "${DRYRUN:-0}" = 1 ] || mkdir -p "$WD/o"
    run python3 -u "$CD/cd_merge.py" --parts "$CD_OUT/parts" --out "$WD/o" --samples "$SAMPLES" --ngroups "${CD_NG:-}" \
        --max-gb "$CD_MAX_GB" --hard-gb "$CD_HARD_GB" --thr "$DP_THR"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/o/dp_*.parquet "$WD"/o/hist_*.npz "$WD/o/summary.json"; do cd_put "$f" "$CD_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; cd_put "$WD/DONE" "$CD_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK CAND_D23 $step $s $grp"
