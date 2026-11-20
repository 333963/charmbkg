#!/bin/bash
# Węzeł TAGGED: wejścia czytane wprost z EOS (fuse), wynik → $TG_OUT. Scratch mktemp + trap; AFS tylko logi.
set -eo pipefail
TG=$(cd "$(dirname "$0")" && pwd); export TG
source "$TG/tg_common.sh"
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/tg_XXXXXX"); trap 'cd /; rm -rf "$WD"' EXIT; cd "$WD"
echo "[$(date -u +%FT%TZ)] TAGGED host $(hostname) TG_APP=$TG_APP TG_OUT=$TG_OUT"
q_env >/dev/null 2>&1; export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
python3 -u "$TG/tg_extract.py" --app "$TG_APP" --h5 "$TG_H5" --wps "$TG_WPS" --out "$WD/o" 2>&1 | tee "$WD/tg.log"
mkdir -p "$TG_OUT" 2>/dev/null || true
for f in "$WD"/o/*.npz "$WD/tg.log"; do put_output "$f" "$TG_OUT/$(basename "$f")"; done
date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$TG_OUT/DONE"
echo "[$(date -u +%FT%TZ)] OK TAGGED"
