#!/bin/bash
# Węzeł DSTWS. Użycie: dstws_job.sh ext <próbka> <część> | dstws_job.sh merge. Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → komendy.
set -eo pipefail
DW=$(cd "$(dirname "$0")" && pwd); export DW
[ -n "${DW_SAMPLES:-}" ] && DW_SAMPLES=${DW_SAMPLES//,/ }   # z penv (przecinki) → spacje
[ -n "${DW_NPARTS:-}" ] && DW_NPARTS=${DW_NPARTS//,/ }
source "$DW/dstws_common.sh"
step=$1; s=${2:-}; part=${3:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/dw_${step}_${s}_${part}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD MPLCONFIGDIR=$WD/mpl
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] DSTWS $step $s $part host $(hostname) WD=$WD DW_OUT=$DW_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then dw_env; mkdir -p "$DW_OUT/parts" 2>/dev/null || true; fi
case $step in
  ext)
    o=dw_${s}_p$part.npz
    run python3 -u "$DW/dstws_extract.py" --sample "$s" --models "$NOMINAL/models/D0_H" --h5dir "$H5V20/$s" --rawdir "$RAWV20/$s" \
        --lxy-min "$DW_LXY" --part "$part" --nparts "$(dw_np "$s")" --out "$WD/$o" ${DW_MAX_FILES:+--max-files "$DW_MAX_FILES"}
    [ "${DRYRUN:-0}" = 1 ] || dw_put "$WD/$o" "$DW_OUT/parts/$o" ;;
  merge)
    run python3 -u "$DW/dstws_merge.py" --parts "$DW_OUT/parts" --out "$WD" --samples "$DW_SAMPLES" --nparts "$DW_NPARTS"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/dw_*.npz "$WD"/dw_summary.json; do dw_put "$f" "$DW_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; dw_put "$WD/DONE" "$DW_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK DSTWS $step $s $part"
