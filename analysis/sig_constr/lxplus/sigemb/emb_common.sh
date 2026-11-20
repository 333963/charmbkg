# emb_common.sh — paczka SIGEMB (17.11; projekt B, SIG_constr: embedding kinematyki prawdziwego D0 z MC w zdarzenia danych; source, nie uruchamiać).
# Wymaga ../common.sh oraz ../trkjet2/tj2_kernel.py (jądro TRKJET2 używane bez zmian: anty-kT, stożki, półkule — nie duplikujemy).
EMB=${EMB:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$EMB/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                    # EOS_OUT, q_env, put_output
: "${RAWV20:=$EOS_OUT/charm_raw_v20}"                    # wejście (odczyt): ROOT v2.0, drzewa Tracks + CharmAnalysis; tylko *dpt3to150*
: "${EMB_OUT:=$EOS_OUT/charm_q/sigemb}"                  # wyjście: parts/emb_<próbka>_p<k>.npz → emb_<próbka>.npz, emb_summary.json, DONE
: "${EMB_SAMPLES:=data_428770 data_427929 mb_nd}"
: "${EMB_NPARTS:=data_428770:12 data_427929:12 mb_nd:5}"
: "${EMB_PRESCALE:=data_428770:20 data_427929:40 mb_nd:40}"  # co który nośnik (hasz eventNumber); cel ~0,5 M nośników na próbkę (test 17.11: 100 k / 45 k / 246 k wpisów na plik)
: "${EMB_DISK:=3GB}"
export EMB Q RAWV20 EMB_OUT EMB_SAMPLES EMB_NPARTS EMB_PRESCALE EMB_DISK
emb_kv() { local kv; for kv in $1; do [ "${kv%%:*}" = "$2" ] && { echo "${kv##*:}"; return; }; done; echo "$3"; }
emb_np() { emb_kv "$EMB_NPARTS" "$1" 4; }
emb_ps() { emb_kv "$EMB_PRESCALE" "$1" 1; }
emb_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2" && rm -f "$1"; fi; }
emb_env() {
  q_env
  export PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
  python3 -c 'import numpy, uproot, awkward, numba' || { echo "BŁĄD: brak numpy/uproot/awkward/numba w q_env"; return 1; }
}
