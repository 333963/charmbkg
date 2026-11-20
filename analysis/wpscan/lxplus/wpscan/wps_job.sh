#!/bin/bash
# Węzeł WPSCAN: pliki app czytane wprost z EOS (fuse), wynik → $WPS_OUT/{wps_hist.npz, wps.log, DONE}. Scratch mktemp + trap; AFS tylko logi.
set -eo pipefail
WPS=$(cd "$(dirname "$0")" && pwd); export WPS
source "$WPS/wps_common.sh"
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/wps_XXXXXX"); trap 'cd /; rm -rf "$WD"' EXIT; cd "$WD"
echo "[$(date -u +%FT%TZ)] WPSCAN host $(hostname) WPS_APP=$WPS_APP WPS_OUT=$WPS_OUT"
q_env >/dev/null 2>&1; export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
python3 -u "$WPS/wps_hist.py" --app "$WPS_APP" --out "$WD/wps_hist.npz" --nrep "$WPS_NREP" 2>&1 | tee "$WD/wps.log"
mkdir -p "$WPS_OUT" 2>/dev/null || true
put_output "$WD/wps_hist.npz" "$WPS_OUT/wps_hist.npz"
put_output "$WD/wps.log" "$WPS_OUT/wps.log"
date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$WPS_OUT/DONE"
echo "[$(date -u +%FT%TZ)] OK WPSCAN"
