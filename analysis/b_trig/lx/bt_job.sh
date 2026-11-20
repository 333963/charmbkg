#!/bin/bash
# Węzeł DAG B_trig. Użycie:
#   bt_job.sh trg <próbka>                 ścieżka A: tabela zdarzeń z ROOT v2.0 → $BT_OUT/trg_<próbka>.npz
#   bt_job.sh ext <próbka> <lista> <grupa> ścieżka B: btrig_extract na h5_v20 → $BT_OUT/parts/bt_<próbka>_<grupa>.npz
#   bt_job.sh merge                        sklejenie części → $BT_OUT/bt_<próbka>.npz + bt_summary.json + DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
BT=$(cd "$(dirname "$0")" && pwd); export BT
source "$BT/bt_common.sh"
step=$1; s=${2:-}; lst=${3:-}; grp=${4:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/bt_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] BT $step $s $grp host $(hostname) WD=$WD BT_OUT=$BT_OUT"
[ "${DRYRUN:-0}" = 1 ] || bt_env
case $step in
  trg)
    mapfile -t FILES < <(ls "$RAWV20/$s"/*dpt3to150*.root)
    run python3 -u "$BT/trg_root_extract.py" --out "$WD/trg_${s}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || bt_put "$WD/trg_${s}.npz" "$BT_OUT/trg_${s}.npz" ;;
  ext)
    mapfile -t FILES < "$lst"
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    run python3 -u "$BT/btrig_extract.py" --models "$NOMINAL/models/D0_H" --tag "$(bt_tag "$s")" --all-frac "$(bt_frac "$s")" \
        --lxy-min "$BT_LXY" --out "$WD/bt_${s}_${grp}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || bt_put "$WD/bt_${s}_${grp}.npz" "$BT_OUT/parts/bt_${s}_${grp}.npz" ;;
  merge)
    run python3 -u "$BT/bt_merge.py" --parts "$BT_OUT/parts" --out "$WD" --need "$BT_SAMPLES"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/bt_*.npz "$WD"/bt_summary.json; do bt_put "$f" "$BT_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; bt_put "$WD/DONE" "$BT_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK BT $step $s $grp"
