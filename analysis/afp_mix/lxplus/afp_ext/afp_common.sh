# afp_common.sh — paczka AFP_EXT (17.11; projekt B: slim evtable + kandydaci D z h5 v2.0 dla AFP-MIX; source, nie uruchamiać). Wymaga ../common.sh.
AX=${AX:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$AX/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                    # EOS_OUT, q_env, put_output
: "${H5V20:=/eos/user/k/klyczak/charm_h5_v20/parts}"     # wejście (tylko odczyt): <próbka>/<próbka>_<NNNNNN>_none.h5 (events + evtable)
: "${AFP_OUT:=$EOS_OUT/charm_q/afp_ext}"                 # wyjście: <próbka>_events.parquet, <próbka>_cands.parquet, .meta.json, afp_summary.json, DONE; części w parts/
: "${AFP_SAMPLES:=data_428770 data_427929}"
: "${AFP_NPARTS:=8}"                                     # węzłów na próbkę
: "${AFP_NFILES:=data_428770:97 data_427929:458}"        # oczekiwane pliki *_none.h5 (h5_v20/samples_v20.txt, 17.11)
: "${AFP_DISK:=4GB}"
export AX Q H5V20 AFP_OUT AFP_SAMPLES AFP_NPARTS AFP_NFILES AFP_DISK
ax_nexp() { local kv; for kv in $AFP_NFILES; do [ "${kv%%:*}" = "$1" ] && { echo "${kv##*:}"; return; }; done; echo 0; }
ax_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2" && rm -f "$1"; fi; }
ax_env() {
  q_env
  export PYTHONDONTWRITEBYTECODE=1
  python3 -c 'import numpy, h5py, pyarrow.parquet' || { echo "BŁĄD: brak numpy/h5py/pyarrow w q_env"; return 1; }
}
