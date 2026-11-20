# mctmpl_common.sh — paczka MCTMPL (17.11.2026; projekt B: test szablonów tła z MC, zadanie #7 planu prac; source, nie uruchamiać).
# Na lxplus: ~/private/charm_q/Q_lxplus/mctmpl (= lokalnie charm_work/Q_lxplus/mctmpl; kopiuje etap ORCH stages/mctmpl.sh). Wymaga ../common.sh.
MT=${MT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$MT/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                        # EOS_OUT, RUNS, q_env, put_output
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                     # wejście (tylko odczyt): <próbka>/<próbka>_NNNNNN_none.h5
: "${NOMINAL:=$RUNS}"                                         # models/{D0_H,Dplus_H} (tylko odczyt)
: "${MT_SRC:=$Q/ev_dstar/src}"                                # charmtag 2.0.2 (jak B_trig/FMM)
: "${MT_OUT:=$EOS_OUT/charm_q/mctmpl}"                        # parts/mctmpl_<próbka>_<g>.npz, mctmpl_<próbka>.npz, mctmpl_summary.json, DONE
: "${MT_WP:=$MT/wp_thresholds.json}"                          # progi WP95/80/50 (zamrożone w paczce; mctmpl_hist.py --thr-out)
: "${MT_SAMPLES:=data_428770 data_427929 mb_nd fd_nd}"
: "${MT_NFILES:=data_428770:97 data_427929:458 mb_nd:78 fd_nd:24}"   # oczekiwane *_none.h5 (h5_v20, 17.11)
: "${MT_GROUP:=data_428770:7 data_427929:46 mb_nd:10 fd_nd:8}"       # plików h5 na węzeł
: "${MT_PRESEL:=none}"
export MT Q H5V20 NOMINAL MT_SRC MT_OUT MT_WP MT_SAMPLES MT_NFILES MT_GROUP MT_PRESEL
mt_kv() { local kv; for kv in $1; do [ "${kv%%:*}" = "$2" ] && { echo "${kv##*:}"; return; }; done; echo "${3:-0}"; }
mt_tag() { case $1 in data_*) echo data ;; fd_*) echo fd ;; mb_*) echo mb ;; *) echo "$1" ;; esac; }
mt_run() { case $1 in data_*) echo "${1#data_}" ;; *) echo 0 ;; esac; }
mt_ngroups() { local n g; n=$(mt_kv "$MT_NFILES" "$1"); g=$(mt_kv "$MT_GROUP" "$1" 10); echo $(( (n + g - 1) / g )); }
mt_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2"; fi; }
mt_env() {
  q_env
  export PYTHONPATH="$MT:$MT_SRC:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1
  python3 -c "import charmtag, charmtag.features, charmtag.model, h5py, yaml; print('charmtag:', charmtag.__file__)" \
    || { echo "BŁĄD: import charmtag/h5py/yaml"; return 1; }
}
