# tj2_common.sh — paczka TRKJET2 (17.11; projekt B: dżety torowe anty-kT i otoczenie D0 dla par D0←D* z b_trig; source, nie uruchamiać). Wymaga ../common.sh.
TJ2=${TJ2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$TJ2/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                    # EOS_OUT, q_env, put_output, fetch_input
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                # wejście (odczyt): <próbka>/<próbka>_<NNNNNN>_none.h5 (wiersz → event_entry, córki)
: "${RAWV20:=$EOS_OUT/charm_raw_v20}"                    # wejście (odczyt): ROOT v2.0, drzewo Tracks; tylko *dpt3to150*
: "${BT_IN:=$EOS_OUT/charm_q/b_trig}"                    # wejście (odczyt): bt_<próbka>.npz z etapu b_trig (pary + score H)
: "${TJ2_OUT:=$EOS_OUT/charm_q/trkjet2}"                   # wyjście: parts/tj2_<próbka>_p<k>.npz → tj2_<próbka>.npz, tj2_summary.json, DONE
: "${TJ2_SAMPLES:=data_428770 data_427929 fd_nd mb_nd}"
: "${TJ2_NPARTS:=data_428770:12 data_427929:12 fd_nd:3 mb_nd:5}"   # węzłów na próbkę (~8 plików ROOT 428770 / ~38 plików 427929 na węzeł)
: "${TJ2_DISK:=3GB}"
export TJ2 Q H5V20 RAWV20 BT_IN TJ2_OUT TJ2_SAMPLES TJ2_NPARTS TJ2_DISK
tj2_np() { local kv; for kv in $TJ2_NPARTS; do [ "${kv%%:*}" = "$1" ] && { echo "${kv##*:}"; return; }; done; echo 4; }
tj2_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2" && rm -f "$1"; fi; }
tj2_env() {
  q_env
  export PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
  python3 -c 'import numpy, h5py, uproot, awkward, numba' || { echo "BŁĄD: brak numpy/h5py/uproot/awkward/numba w q_env"; return 1; }
}
