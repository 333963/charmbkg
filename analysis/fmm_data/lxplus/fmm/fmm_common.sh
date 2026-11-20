# fmm_common.sh — paczka FMM (17.11.2026; projekt B: pary bliźniacze D0 K↔π → f_MM pod pikiem; source, nie uruchamiać). Wymaga ../common.sh.
# Na lxplus: ~/private/charm_q/Q_lxplus/fmm (= lokalnie charm_work/Q_lxplus/fmm; kopiuje etap ORCH stages/fmm.sh).
FM=${FM:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$FM/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                        # EOS_OUT, RUNS, q_env, put_output
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                     # wejście (tylko odczyt): <próbka>/<próbka>_NNNNNN_none.h5
: "${NOMINAL:=$RUNS}"                                         # models/D0_H (tylko odczyt) — ten sam model co B_trig/EV
: "${FMM_SRC:=$Q/ev_dstar/src}"                               # charmtag 2.0.2 (jak B_trig)
: "${FMM_OUT:=$EOS_OUT/charm_q/fmm}"                          # parts/fmm_<próbka>_<g>.npz, fmm_<próbka>.npz, fmm_summary.json, DONE
: "${FMM_SAMPLES:=data_428770 data_427929 mb_nd fd_nd}"
: "${FMM_NFILES:=data_428770:97 data_427929:458 mb_nd:78 fd_nd:24}"   # oczekiwane *_none.h5 (h5_v20, 17.11)
: "${FMM_GROUP:=data_428770:10 data_427929:46 mb_nd:10 fd_nd:8}"      # plików h5 na węzeł
: "${FMM_FRAC:=data_428770:0.3 data_427929:1.0 mb_nd:1.0 fd_nd:1.0}"   # podpróbka zdarzeń (428770: 30 % wystarcza — systematyka dominuje)
: "${FMM_LXY:=0.25}"                                          # L_xy ≥ 0,25 mm = nominał v1.2 (h5_v20 ma 0,10)
: "${FMM_PRESEL:=none}"
export FM Q H5V20 NOMINAL FMM_SRC FMM_OUT FMM_SAMPLES FMM_NFILES FMM_GROUP FMM_FRAC FMM_LXY FMM_PRESEL
fm_kv() { local kv; for kv in $1; do [ "${kv%%:*}" = "$2" ] && { echo "${kv##*:}"; return; }; done; echo "${3:-0}"; }
fm_tag() { case $1 in data_*) echo data ;; fd_*) echo fd ;; mb_*) echo mb ;; *) echo "$1" ;; esac; }
fm_ngroups() { local n g; n=$(fm_kv "$FMM_NFILES" "$1"); g=$(fm_kv "$FMM_GROUP" "$1" 10); echo $(( (n + g - 1) / g )); }
fm_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2"; fi; }
fm_env() {
  q_env
  export PYTHONPATH="$FM:$FMM_SRC:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1
  python3 -c "import charmtag, charmtag.features, charmtag.model, h5py, yaml; print('charmtag:', charmtag.__file__)" \
    || { echo "BŁĄD: import charmtag/h5py/yaml"; return 1; }
}
