# s427_common.sh — paczka SALAD427 (17.11.2026; projekt B: SALAD v1.0 na runie 427929 (MBTS) i w klasach trygera 428770; source).
# Na lxplus: ~/private/charm_q/Q_lxplus/salad427 (= lokalnie charm_work/Q_lxplus/salad427; kopiuje etap ORCH stages/salad427.sh).
# Wymaga ../common.sh. TYLKO pojedyncze D (D0, D+), tylko pasma boczne i tło. Par D-D paczka nie czyta.
S4=${S4:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$S4/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                        # EOS_OUT, q_env, put_output
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                     # wejście (odczyt): <próbka>/<próbka>_NNNNNN_none.h5
: "${S4_OUT:=$EOS_OUT/charm_q/salad427}"                      # parts/s427_<próbka>_<g>.npz, s427_<gatunek>{,_adv}.json, s427_<gatunek>_w.npz, DONE
: "${S4_SAMPLES:=data_427929 data_428770 mb_nd fd_nd}"
: "${S4_GROUP:=data_428770:5 data_427929:46 mb_nd:13 fd_nd:8}" # plików h5 na węzeł ekstrakcji
: "${S4_FRAC_J12:=0.25}"                                      # podpróbka zdarzeń klasy J12only w 428770 (reszta klas 100 %)
: "${S4_LXY:=0.25}"                                           # L_xy >= 0,25 mm = nominał (h5_v20 ma 0,10)
: "${S4_SPECIES:=D0 Dplus}"
: "${S4_TASKS:=main,adv}"
: "${S4_CPUS:=8}"
: "${S4_FAST:=0}"                                             # 1 = test (2 foldy, 30 drzew, 20 replik)
: "${S4_PRESEL:=none}"
export S4 Q H5V20 S4_OUT S4_SAMPLES S4_GROUP S4_FRAC_J12 S4_LXY S4_SPECIES S4_TASKS S4_CPUS S4_FAST S4_PRESEL
s4_kv() { local kv; for kv in $1; do [ "${kv%%:*}" = "$2" ] && { echo "${kv##*:}"; return; }; done; echo "${3:-10}"; }
s4_tag() { case $1 in data_*) echo data ;; fd_*) echo fd ;; mb_*) echo mb ;; *) echo "$1" ;; esac; }
s4_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else mkdir -p "$(dirname "$2")" && cp -f "$1" "$2"; fi; }
s4_env() {
  [ "${S4_NOENV:-0}" = 1 ] || q_env
  export PYTHONPATH="$S4:$S4/src:${PYTHONPATH:-}" CHARMBKG_SRC="$S4/src" PYTHONDONTWRITEBYTECODE=1 S427_FAST=$S4_FAST \
         OMP_NUM_THREADS=$S4_CPUS MKL_NUM_THREADS=$S4_CPUS OPENBLAS_NUM_THREADS=$S4_CPUS
  python3 -c "import charmbkg, h5py, yaml, sklearn, pandas, scipy; from charmbkg.bkggen import salad_v1 as S; print('charmbkg:', charmbkg.__file__, 'sklearn', sklearn.__version__, 'salad yaml', S.load_cfg()['version'])" \
    || { echo "BŁĄD: import charmbkg/h5py/sklearn"; return 1; }
}
