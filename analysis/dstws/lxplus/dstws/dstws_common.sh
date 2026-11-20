# dstws_common.sh — paczka DSTWS (17.11; projekt B: tło D* przez konstrukcję — D0 × tor z drzewa Tracks: RS/WS/ROT/MIX; source, nie uruchamiać).
# Wymaga ../common.sh. Na lxplus: ~/private/charm_q/Q_lxplus/dstws.
DW=${DW:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$DW/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                    # EOS_OUT, RUNS, q_env, put_output
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                # wejście (odczyt): <próbka>/<próbka>_<NNNNNN>_none.h5 (kandydaci D0/D*, cechy H)
: "${RAWV20:=$EOS_OUT/charm_raw_v20}"                    # wejście (odczyt): ROOT v2.0, drzewo Tracks; tylko *dpt3to150*
: "${NOMINAL:=$RUNS}"                                    # models/D0_H (odczyt) — ten sam model co EV/B_trig
: "${DW_SRC:=$Q/ev_dstar/src}"                           # charmtag 2.0.2 (jak B_trig)
: "${DW_OUT:=$EOS_OUT/charm_q/dstws}"                    # wyjście: parts/dw_<próbka>_p<k>.npz → dw_<próbka>.npz, dw_summary.json, DONE
: "${DW_SAMPLES:=data_428770 data_427929 fd_nd mb_nd}"
: "${DW_NPARTS:=data_428770:12 data_427929:6 fd_nd:3 mb_nd:4}"   # węzłów na próbkę (~8 plików 428770 ≈ 10 min / 76 plików 427929 ≈ 5 min / 8 fd / 20 mb na węzeł; test 17.11)
: "${DW_LXY:=0.25}"
: "${DW_DISK:=3GB}"; : "${DW_MEM:=8GB}"
export DW Q H5V20 RAWV20 NOMINAL DW_SRC DW_OUT DW_SAMPLES DW_NPARTS DW_LXY DW_DISK DW_MEM
dw_np() { local kv; for kv in $DW_NPARTS; do [ "${kv%%:*}" = "$1" ] && { echo "${kv##*:}"; return; }; done; echo 4; }
dw_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2" && rm -f "$1"; fi; }
dw_env() {
  q_env
  export PYTHONPATH="$DW:$DW_SRC:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
  python3 -c 'import numpy, pandas, h5py, uproot, awkward, numba, yaml, charmtag.features, charmtag.model, charmtag.samples' \
    || { echo "BŁĄD: brak numpy/h5py/uproot/awkward/numba/charmtag w q_env"; return 1; }
}
