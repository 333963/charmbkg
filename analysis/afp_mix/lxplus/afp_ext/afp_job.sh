#!/bin/bash
# Węzeł AFP_EXT. Użycie: afp_job.sh <próbka> <część> | afp_job.sh finish. Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi.
set -eo pipefail
AX=$(cd "$(dirname "$0")" && pwd); export AX
[ -n "${AFP_SAMPLES:-}" ] && AFP_SAMPLES=${AFP_SAMPLES//,/ }   # z penv (przecinki) → spacje
[ -n "${AFP_NFILES:-}" ] && AFP_NFILES=${AFP_NFILES//,/ }
source "$AX/afp_common.sh"
step=$1; part=${2:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/afp_${step}_${part}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] AFP_EXT $step/$part host $(hostname) WD=$WD AFP_OUT=$AFP_OUT H5V20=$H5V20"
if [ "${DRYRUN:-0}" != 1 ]; then ax_env; mkdir -p "$AFP_OUT/parts" 2>/dev/null || true; fi
case $step in
  finish)
    run python3 -u "$AX/afp_merge.py" --out "$AFP_OUT" --wd "$WD" --samples "$AFP_SAMPLES" --nparts "$AFP_NPARTS" --nfiles "$AFP_NFILES"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/*_events.parquet "$WD"/*_cands.parquet "$WD"/*.meta.json "$WD"/afp_summary.json; do ax_put "$f" "$AFP_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; ax_put "$WD/DONE" "$AFP_OUT/DONE"
    fi ;;
  data_*)
    o=$step.p$part
    run python3 -u "$AX/afp_extract.py" --sample "$step" --h5dir "$H5V20/$step" --part "$part" --nparts "$AFP_NPARTS" \
        --out-prefix "$WD/$o" --wd "$WD" ${AFP_MAX_FILES:+--max-files "$AFP_MAX_FILES"}
    if [ "${DRYRUN:-0}" != 1 ]; then
      for x in events.parquet cands.parquet meta.json; do ax_put "$WD/$o.$x" "$AFP_OUT/parts/$o.$x"; done
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK AFP_EXT $step/$part"
