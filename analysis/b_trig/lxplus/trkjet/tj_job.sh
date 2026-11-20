#!/bin/bash
# Węzeł TRKJET. Użycie: tj_job.sh ext <próbka> <część> | tj_job.sh merge. Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → komendy.
set -eo pipefail
TJ=$(cd "$(dirname "$0")" && pwd); export TJ
[ -n "${TJ_SAMPLES:-}" ] && TJ_SAMPLES=${TJ_SAMPLES//,/ }   # z penv (przecinki) → spacje
[ -n "${TJ_NPARTS:-}" ] && TJ_NPARTS=${TJ_NPARTS//,/ }
source "$TJ/tj_common.sh"
step=$1; s=${2:-}; part=${3:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/tj_${step}_${s}_${part}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] TRKJET $step $s $part host $(hostname) WD=$WD TJ_OUT=$TJ_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then tj_env; mkdir -p "$TJ_OUT/parts" 2>/dev/null || true; fi
case $step in
  ext)
    bt=$BT_IN/bt_$s.npz
    [ "${DRYRUN:-0}" = 1 ] || bt=$(fetch_input "$BT_IN/bt_$s.npz" "$WD")      # kopia na scratch (xrdcp), inaczej FUSE
    o=tj_${s}_p$part.npz
    run python3 -u "$TJ/tj_extract.py" --sample "$s" --bt "$bt" --h5dir "$H5V20/$s" --rawdir "$RAWV20/$s" \
        --part "$part" --nparts "$(tj_np "$s")" --out "$WD/$o" ${TJ_MAX_FILES:+--max-files "$TJ_MAX_FILES"}
    [ "${DRYRUN:-0}" = 1 ] || tj_put "$WD/$o" "$TJ_OUT/parts/$o" ;;
  merge)
    run python3 -u "$TJ/tj_merge.py" --parts "$TJ_OUT/parts" --out "$WD" --samples "$TJ_SAMPLES" --nparts "$TJ_NPARTS" --bt "$BT_IN"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/tj_*.npz "$WD"/tj_summary.json; do tj_put "$f" "$TJ_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; tj_put "$WD/DONE" "$TJ_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK TRKJET $step $s $part"
