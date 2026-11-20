#!/bin/bash
# Węzeł DAG DLZ_DATA. Użycie:
#   dlz_job.sh frm <fd|mb|sb> <lista> <g>     ramki treningowe z h5 v2.0 → $DLZ_OUT/parts/frm_<tag>_<g>.npz
#   dlz_job.sh train                          trening H / H+DLZ / H+dalitz (PROCEDURA §2) → $DLZ_OUT/train/{models,oof,thr.json}
#   dlz_job.sh ext <próbka> <lista> <g>       apply na danych → $DLZ_OUT/parts/dlz_<próbka>_<g>.npz
#   dlz_job.sh merge                          → $DLZ_OUT/dlz_<próbka>.npz + dlz_summary.json
#   dlz_job.sh ana                            fity masy/bootstrap/sPlot per run + kontrole MC → $DLZ_OUT/ana/
#   dlz_job.sh d0                             przesiew nowych cech D⁰ (MC + SB danych, h5 v1.2) → $DLZ_OUT/d0/
#   dlz_job.sh fin                            DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
DZ=$(cd "$(dirname "$0")" && pwd); export DZ
for v in DLZ_DATA DLZ_NFILES DLZ_GROUP DLZ_SETS; do [ -n "${!v:-}" ] && export "$v=${!v//,/ }"; done   # z penv (przecinki) → spacje
source "$DZ/dlz_common.sh"
step=$1; a2=${2:-}; lst=${3:-}; grp=${4:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/dlz_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */dlz_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
putdir() { local f; ( cd "$1" && find . -type f ) | while read -r f; do dz_put "$1/${f#./}" "$2/${f#./}"; done; }
echo "[$(date -u +%FT%TZ)] DLZ_DATA $step $a2 $grp host $(hostname) WD=$WD DLZ_OUT=$DLZ_OUT"
[ "${DRYRUN:-0}" = 1 ] || dz_env
case $step in
  frm)
    mapfile -t FILES < "$lst"; [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    X=(); [ "$a2" = sb ] && X=(--frac "$DLZ_SBFRAC")
    run python3 -u "$DZ/dlz_frames.py" --tag "$a2" "${X[@]}" --lxy-min "$DLZ_LXY" --out "$WD/frm_${a2}_${grp}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || dz_put "$WD/frm_${a2}_${grp}.npz" "$DLZ_OUT/parts/frm_${a2}_${grp}.npz" ;;
  train)
    run python3 -u "$DZ/dlz_train.py" --parts "$DLZ_OUT/parts" --out "$WD/train"
    [ "${DRYRUN:-0}" = 1 ] || putdir "$WD/train" "$DLZ_OUT/train" ;;
  ext)
    mapfile -t FILES < "$lst"; [ ${#FILES[@]} -gt 0 ] || { echo "pusta lista $lst"; exit 2; }
    run python3 -u "$DZ/apply.py" --models "$DLZ_OUT/train/models" --thr "$DLZ_OUT/train/thr.json" --nthread 2 --lxy-min "$DLZ_LXY" \
        --out "$WD/dlz_${a2}_${grp}.npz" "${FILES[@]}"
    [ "${DRYRUN:-0}" = 1 ] || dz_put "$WD/dlz_${a2}_${grp}.npz" "$DLZ_OUT/parts/dlz_${a2}_${grp}.npz" ;;
  merge)
    NG=""; for x in $DLZ_DATA; do NG="$NG $x:$(dz_ngroups "$x")"; done
    run python3 -u "$DZ/dlz_merge.py" --parts "$DLZ_OUT/parts" --out "$WD" --samples "$DLZ_DATA" --ngroups "$NG"
    if [ "${DRYRUN:-0}" != 1 ]; then for f in "$WD"/dlz_*.npz "$WD"/dlz_summary.json; do dz_put "$f" "$DLZ_OUT/$(basename "$f")"; done; fi ;;
  ana)
    export DLZ_BASE=$WD/base
    if [ "${DRYRUN:-0}" != 1 ]; then
      mkdir -p "$DLZ_BASE/slim" "$DLZ_BASE/out"
      cp -r "$DLZ_OUT/train/oof" "$DLZ_OUT/train/thr.json" "$DLZ_BASE/"; cp -r "$DLZ_OUT/train/models" "$DLZ_BASE/"
      for x in $DLZ_DATA; do cp "$DLZ_OUT/dlz_$x.npz" "$DLZ_BASE/slim/data_${x#data_}.npz"; done
    fi
    run python3 -u "$DZ/mc_checks.py"
    for x in $DLZ_DATA; do run python3 -u "$DZ/analyze.py" "${x#data_}" --nboot "$DLZ_NBOOT"; done
    [ "${DRYRUN:-0}" = 1 ] || putdir "$DLZ_BASE/out" "$DLZ_OUT/ana" ;;
  d0)
    export D0_BASE=$WD/d0 FEAT_SCRIPTS=$DZ/d0 D0_FRAMES_IN=$WD/fr D0_NJ=$DLZ_NJ D0_MODE=${D0_MODE:-full} FEAT_DIR=$WD/nofeat DSTAR_ISO_FIT=$DZ/d0/dstar_iso_fit.json
    P="python3 -u"; R=$DZ/d0/d0_run.py
    if [ "${DRYRUN:-0}" != 1 ]; then
      mkdir -p "$WD/fr" "$D0_BASE/out" "$WD/in/fd_nd" "$WD/in/mb_nd" "$WD/in/data_428770"; cp "$DZ/d0/dstar_iso_fit.json" "$D0_BASE/out/"
      ls "$H5V12/fd_nd/"*_none.h5 | sort | { if [ "$D0_FD" = all ]; then cat; else head -n "$D0_FD"; fi; } | xargs -I{} ln -s {} "$WD/in/fd_nd/"
      ls "$H5V12/mb_nd/"*_none.h5 | sort | head -n "$D0_MB" | xargs -I{} ln -s {} "$WD/in/mb_nd/"
      nd=$(ls "$H5V12/data_428770/"*_none.h5 | wc -l); st=$(( nd / D0_DATA )); [ "$st" -ge 1 ] || st=1
      ls "$H5V12/data_428770/"*_none.h5 | sort | awk -v s="$st" -v n="$D0_DATA" '(NR-1)%s==0 && c<n {print; c++}' | xargs -I{} ln -s {} "$WD/in/data_428770/"
      ls "$WD"/in/* | head -50
    fi
    for s in fd_nd mb_nd data_428770; do run $P "$DZ/d0/extract.py" "$WD/in/$s" "$WD/fr/${s}_D0.pkl" D0; done
    run $P "$DZ/d0/d0_feats.py"
    run $P "$R" study --test mc --variants "$D0_VARIANTS"
    run $P "$R" transfer
    run $P "$R" sculpt mc
    run $P "$DZ/d0/d0_sculpt.py"
    run $P "$R" modelling
    run $P "$R" study --test data --variants "$D0_VARIANTS_DATA"
    run $P "$R" sculpt data
    if [ "${DRYRUN:-0}" != 1 ]; then
      DO=$D0_BASE/out; [ "$D0_MODE" = screen ] && DO=$D0_BASE/screen/out
      ( cd "$DO" && find . -type f \( -name '*.json' -o -name '*.csv' -o -name '*.log' \) ) | while read -r f; do dz_put "$DO/${f#./}" "$DLZ_OUT/d0/${f#./}"; done
      ( cd "$DO" && find . -type f -name 'scores_*.npz' ) | while read -r f; do dz_put "$DO/${f#./}" "$DLZ_OUT/d0_scores/${f#./}"; done
    fi ;;
  fin)
    if [ "${DRYRUN:-0}" != 1 ]; then date -u +%FT%TZ > "$WD/DONE"; dz_put "$WD/DONE" "$DLZ_OUT/DONE"; fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK DLZ_DATA $step $a2 $grp"
