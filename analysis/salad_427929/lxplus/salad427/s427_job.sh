#!/bin/bash
# Węzeł DAG SALAD427. Użycie:
#   s427_job.sh ext <próbka> <lista> <grupa>   s427_extract na h5_v20 → $S4_OUT/parts/s427_<próbka>_<grupa>.npz
#   s427_job.sh study <gatunek>                s427_study (tagger H427, SALAD k-fold per klasa, transfer, B1) → $S4_OUT/s427_<gatunek>*.json|npz
#   s427_job.sh fin                            sprawdza komplet json → $S4_OUT/DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
S4=$(cd "$(dirname "$0")" && pwd); export S4
[ -n "${S4_SAMPLES:-}" ] && S4_SAMPLES=${S4_SAMPLES//,/ }     # z penv (przecinki) → spacje
[ -n "${S4_SPECIES:-}" ] && S4_SPECIES=${S4_SPECIES//,/ }
[ -n "${S4_TASKS:-}" ] && S4_TASKS=${S4_TASKS//+/,}
source "$S4/s427_common.sh"
step=$1; s=${2:-}; lst=${3:-}; grp=${4:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/s427_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */s427_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }   # tylko własny katalog mktemp
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] SALAD427 $step $s $grp host $(hostname) WD=$WD S4_OUT=$S4_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then s4_env; mkdir -p "$S4_OUT/parts" 2>/dev/null || true; fi
case $step in
  ext)
    mapfile -t FILES < "$lst"
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    run python3 -u "$S4/s427_extract.py" --tag "$(s4_tag "$s")" --sample "$s" --lxy-min "$S4_LXY" --frac-j12 "$S4_FRAC_J12" \
        --out "$WD/s427_${s}_${grp}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || s4_put "$WD/s427_${s}_${grp}.npz" "$S4_OUT/parts/s427_${s}_${grp}.npz" ;;
  study)
    mkdir -p "$WD/out"
    run python3 -u "$S4/s427_study.py" --species "$s" --parts "$S4_OUT/parts" --out "$WD/out" --tasks "$S4_TASKS" --threads "$S4_CPUS"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/out/s427_"$s"*.json "$WD"/out/s427_"$s"_w.npz; do [ -f "$f" ] && s4_put "$f" "$S4_OUT/$(basename "$f")"; done
    fi ;;
  fin)
    if [ "${DRYRUN:-0}" != 1 ]; then
      for sp in $S4_SPECIES; do [ -s "$S4_OUT/s427_$sp.json" ] || { echo "brak $S4_OUT/s427_$sp.json"; exit 3; }; done
      date -u +%FT%TZ > "$WD/DONE"; s4_put "$WD/DONE" "$S4_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK SALAD427 $step $s $grp"
