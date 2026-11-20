#!/bin/bash
# Węzeł SIGEMB. Użycie: emb_job.sh ext <próbka> <część> | emb_job.sh merge. Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → komendy.
set -eo pipefail
EMB=$(cd "$(dirname "$0")" && pwd); export EMB
[ -n "${EMB_SAMPLES:-}" ] && EMB_SAMPLES=${EMB_SAMPLES//,/ }
[ -n "${EMB_NPARTS:-}" ] && EMB_NPARTS=${EMB_NPARTS//,/ }
[ -n "${EMB_PRESCALE:-}" ] && EMB_PRESCALE=${EMB_PRESCALE//,/ }
source "$EMB/emb_common.sh"
step=$1; s=${2:-}; part=${3:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/emb_${step}_${s}_${part}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] SIGEMB $step $s $part host $(hostname) WD=$WD EMB_OUT=$EMB_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then emb_env; mkdir -p "$EMB_OUT/parts" 2>/dev/null || true; fi
case $step in
  ext)
    o=emb_${s}_p$part.npz
    run python3 -u "$EMB/emb_extract.py" --sample "$s" --rawdir "$RAWV20/$s" --donors "$EMB/emb_donors.npz" --kernel-dir "$Q/trkjet2" \
        --part "$part" --nparts "$(emb_np "$s")" --prescale "$(emb_ps "$s")" --out "$WD/$o" ${EMB_MAX_FILES:+--max-files "$EMB_MAX_FILES"}
    [ "${DRYRUN:-0}" = 1 ] || emb_put "$WD/$o" "$EMB_OUT/parts/$o" ;;
  merge)
    run python3 -u "$EMB/emb_merge.py" --parts "$EMB_OUT/parts" --out "$WD" --samples "$EMB_SAMPLES" --nparts "$EMB_NPARTS"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/emb_*.npz "$WD"/emb_summary.json; do emb_put "$f" "$EMB_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; emb_put "$WD/DONE" "$EMB_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK SIGEMB $step $s $part"
