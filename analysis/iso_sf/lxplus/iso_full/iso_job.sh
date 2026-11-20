#!/bin/bash
# Węzeł ISO-FULL. Użycie: iso_job.sh <próbka> <część> | iso_job.sh finish. Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi.
set -eo pipefail
IX=$(cd "$(dirname "$0")" && pwd); export IX
[ -n "${ISO_SAMPLES:-}" ] && ISO_SAMPLES=${ISO_SAMPLES//,/ }   # z penv (przecinki) → spacje
source "$IX/iso_common.sh"
step=$1; part=${2:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/iso_${step}_${part}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] ISO $step/$part host $(hostname) WD=$WD ISO_OUT=$ISO_OUT TAG_DIR=$TAG_DIR"
if [ "${DRYRUN:-0}" != 1 ]; then ix_env; mkdir -p "$ISO_OUT/parts"; fi
case $step in
  finish)
    run python3 -u "$IX/iso_merge.py" --out "$ISO_OUT" --tag-dir "$TAG_DIR" --wd "$WD" --samples "$ISO_SAMPLES" --nparts-data "$ISO_NPARTS_DATA"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/*_tag.parquet "$WD"/*_tag.meta.json "$WD"/iso_summary.json; do ix_put "$f" "$ISO_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; ix_put "$WD/DONE" "$ISO_OUT/DONE"
    fi ;;
  data_428770|fd_nd|mb_nd)
    n=$(ix_nparts "$step"); o=${step}_tag.p${part}
    run python3 -u "$IX/iso_extract.py" --sample "$step" --tag "$TAG_DIR/$step.parquet" --h5dir "$H5V12/$step" \
        --part "$part" --nparts "$n" --out "$WD/$o.parquet"
    if [ "${DRYRUN:-0}" != 1 ]; then
      ix_put "$WD/$o.parquet" "$ISO_OUT/parts/$o.parquet"; ix_put "$WD/$o.meta.json" "$ISO_OUT/parts/$o.meta.json"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK ISO $step/$part"
