#!/bin/bash
# Węzeł / przebieg BKGCAT: bkg_job.sh <k> <n>  — k-ta z n grup plików mb_nd (wg numeru NNNNNN) → $BKG_OUT/parts/bkg_<NNNNNN>.npz (gotowe pomijane).
#                          bkg_job.sh merge    — suma parts → $BKG_OUT/{bkgcat.npz, DONE}.
# Wejścia czytane wprost z EOS (fuse). Scratch: $BKG_WD (interaktywnie, np. /tmp/klyczak_bkgcat — nie kasowany) albo mktemp + trap (condor). AFS tylko logi.
# BKG_TMAX [s] — po tym czasie nie zaczyna kolejnego pliku (przebiegi interaktywne ≤ 10 min). Bez AFP, tylko MC.
set -eo pipefail
BKG=$(cd "$(dirname "$0")" && pwd); export BKG
source "$BKG/bkg_common.sh"
set -u
if [ -n "${BKG_WD:-}" ]; then WD=$BKG_WD; mkdir -p "$WD"; else WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/bkg_XXXXXX"); trap 'cd /; rm -rf "$WD"' EXIT; fi
cd "$WD"; mkdir -p parts
q_env >/dev/null 2>&1; export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
T0=$(date +%s)
if [ "${1:-}" = merge ]; then
  for f in "$BKG_OUT"/parts/bkg_*.npz; do [ -s "parts/$(basename "$f")" ] || cp "$f" parts/; done
  python3 -u "$BKG/bkg_merge.py" --parts parts --thr "$BKG_THR" --out bkgcat.npz
  put_output bkgcat.npz "$BKG_OUT/bkgcat.npz"; date -u +%FT%TZ > DONE; put_output DONE "$BKG_OUT/DONE"; echo "OK merge"; exit 0
fi
K=$1; N=$2; i=0
for h in $(ls "$BKG_H5"/mb_nd_*_none.h5 | sort); do
  id=$(basename "$h" | sed 's/mb_nd_\([0-9]*\)_none.h5/\1/'); i=$((i + 1))
  [ $(( (i - 1) % N )) -eq "$K" ] || continue
  o="bkg_$id.npz"
  [ -s "$BKG_OUT/parts/$o" ] && { echo "$o już jest — pomijam"; continue; }
  [ $(( $(date +%s) - T0 )) -gt "${BKG_TMAX:-100000}" ] && { echo "BKG_TMAX — przerwano przed $id"; exit 3; }
  python3 -u "$BKG/bkg_hist.py" --fid "$id" --app "$BKG_APP" --h5 "$BKG_H5" --raw "$BKG_RAW" --thr "$BKG_THR" --out "parts/$o"
  put_output "parts/$o" "$BKG_OUT/parts/$o"
done
echo "[$(date -u +%FT%TZ)] OK grupa $K/$N"
