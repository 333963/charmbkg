#!/bin/bash
# Węzeł DAG MCTMPL. Użycie:
#   mctmpl_job.sh ext <próbka> <lista> <grupa>   mctmpl_hist na h5_v20 → $MT_OUT/parts/mctmpl_<próbka>_<grupa>.npz
#   mctmpl_job.sh merge                          suma histogramów → $MT_OUT/mctmpl_<próbka>.npz + mctmpl_summary.json + DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
MT=$(cd "$(dirname "$0")" && pwd); export MT
for v in MT_SAMPLES MT_NFILES MT_GROUP; do [ -n "${!v:-}" ] && export "$v=${!v//,/ }"; done   # z penv (przecinki) → spacje
source "$MT/mctmpl_common.sh"
step=$1; s=${2:-}; lst=${3:-}; grp=${4:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/mctmpl_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */mctmpl_*) rm -rf "${WD:?}" ;; esac; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] MCTMPL $step $s $grp host $(hostname) WD=$WD MT_OUT=$MT_OUT"
[ "${DRYRUN:-0}" = 1 ] || mt_env
case $step in
  ext)
    mapfile -t FILES < "$lst"
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    run python3 -u "$MT/mctmpl_hist.py" --models "$NOMINAL/models" --wp "$MT_WP" --tag "$(mt_tag "$s")" --run "$(mt_run "$s")" \
        --out "$WD/mctmpl_${s}_${grp}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || mt_put "$WD/mctmpl_${s}_${grp}.npz" "$MT_OUT/parts/mctmpl_${s}_${grp}.npz" ;;
  merge)
    NG=""; for x in $MT_SAMPLES; do NG="$NG $x:$(mt_ngroups "$x")"; done
    run python3 -u "$MT/mctmpl_merge.py" --parts "$MT_OUT/parts" --out "$WD" --samples "$MT_SAMPLES" --ngroups "$NG" --nfiles "$MT_NFILES"
    # analiza (lekka, ~kilka min) od razu na condorze → $MT_OUT/ana/{WYNIKI.md,results.json}; błąd analizy nie blokuje histogramów
    run python3 -u "$MT/mctmpl_analyze.py" --in "$WD" --out "$WD/ana" || echo "UWAGA: mctmpl_analyze.py nie powiodło się (histogramy są; analiza lokalnie)"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/mctmpl_*.npz "$WD"/mctmpl_summary.json; do mt_put "$f" "$MT_OUT/$(basename "$f")"; done
      for f in "$WD"/ana/WYNIKI.md "$WD"/ana/results.json; do [ -s "$f" ] && mt_put "$f" "$MT_OUT/ana/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; mt_put "$WD/DONE" "$MT_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK MCTMPL $step $s $grp"
