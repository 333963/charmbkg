#!/bin/bash
# Węzeł TRKJET2. Użycie: tj2_job.sh ext <próbka> <część> | tj2_job.sh merge. Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → komendy.
set -eo pipefail
TJ2=$(cd "$(dirname "$0")" && pwd); export TJ2
[ -n "${TJ2_SAMPLES:-}" ] && TJ2_SAMPLES=${TJ2_SAMPLES//,/ }   # z penv (przecinki) → spacje
[ -n "${TJ2_NPARTS:-}" ] && TJ2_NPARTS=${TJ2_NPARTS//,/ }
source "$TJ2/tj2_common.sh"
step=$1; s=${2:-}; part=${3:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/tj2_${step}_${s}_${part}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] TRKJET2 $step $s $part host $(hostname) WD=$WD TJ2_OUT=$TJ2_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then tj2_env; mkdir -p "$TJ2_OUT/parts" 2>/dev/null || true; fi
case $step in
  ext)
    bt=$BT_IN/bt_$s.npz
    [ "${DRYRUN:-0}" = 1 ] || bt=$(fetch_input "$BT_IN/bt_$s.npz" "$WD")      # kopia na scratch (xrdcp), inaczej FUSE
    o=tj2_${s}_p$part.npz
    run python3 -u "$TJ2/tj2_extract.py" --sample "$s" --bt "$bt" --h5dir "$H5V20/$s" --rawdir "$RAWV20/$s" \
        --part "$part" --nparts "$(tj2_np "$s")" --out "$WD/$o" ${TJ2_MAX_FILES:+--max-files "$TJ2_MAX_FILES"}
    [ "${DRYRUN:-0}" = 1 ] || tj2_put "$WD/$o" "$TJ2_OUT/parts/$o" ;;
  merge)
    run python3 -u "$TJ2/tj2_merge.py" --parts "$TJ2_OUT/parts" --out "$WD" --samples "$TJ2_SAMPLES" --nparts "$TJ2_NPARTS" --bt "$BT_IN"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/tj2_*.npz "$WD"/tj2_summary.json; do tj2_put "$f" "$TJ2_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; tj2_put "$WD/DONE" "$TJ2_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK TRKJET2 $step $s $part"
