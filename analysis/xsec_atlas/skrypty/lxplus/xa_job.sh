#!/bin/bash
# Wezel DAG XSEC_atlas.  xa_job.sh ext <probka> <k> <n>  -> $XA_OUT/parts/xa_<probka>_<k>.npz ;  xa_job.sh fin -> $XA_OUT/DONE
# Srodowisko i sciezki z ../incal/inc_common.sh (tylko source; bez zmian). Scratch mktemp + trap; wyjscia tylko EOS; AFS tylko logi.
set -eo pipefail
XA=$(cd "$(dirname "$0")" && pwd)
INC=$(cd "$XA/../incal" && pwd); export INC
source "$INC/inc_common.sh"
: "${XA_OUT:=$EOS_OUT/charm_q/xsatlas}"
: "${XA_NODES:=data_427929:8 mb_nd:6 fd_nd:3}"
XA_NODES=${XA_NODES//+/ }
step=$1; shift
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/xa_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */xa_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] XSATLAS $step $* host $(hostname) WD=$WD XA_OUT=$XA_OUT"
case $step in
  ext)
    s=$1; k=$2; n=$3; o="xa_${s}_${k}.npz"
    inc_env
    export OMP_NUM_THREADS=1
    python3 -u "$XA/xa_extract.py" --sample "$s" --k "$k" --n "$n" --h5dir "$H5V20" --models "$INC_MODELS" --max-files "${XA_MAXFILES:-0}" --out "$WD/$o"
    put_output "$WD/$o" "${XA_DEST:-$XA_OUT/parts}/$o" ;;
  fin)
    miss=0
    for x in $XA_NODES; do
      IFS=: read -r s n <<< "$x"
      for ((k = 0; k < n; k++)); do [ -s "$XA_OUT/parts/xa_${s}_${k}.npz" ] || { echo "BRAK xa_${s}_${k}.npz"; miss=$((miss + 1)); }; done
    done
    [ "$miss" = 0 ] || { echo "BLAD: brak $miss czesci - bez DONE"; exit 1; }
    date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$XA_OUT/DONE" ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK XSATLAS $step $*"
