# bt_common.sh — paczka B_trig (projekt B: odporność SF/ε_B/izolacji na trigger i run; source, nie uruchamiać). Wymaga ../common.sh.
# Na lxplus: ~/private/charm_q/Q_lxplus/b_trig (= lokalnie charm_work/B_trig/lx + *.py z charm_work/B_trig, kopiowane przez etap ORCH).
BT=${BT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$BT/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                        # EOS_OUT, RUNS, q_env, put_output
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                     # h5_v20 (etap h5_v20): <próbka>/<próbka>_NNNNNN_none.h5
: "${RAWV20:=$EOS_OUT/charm_raw_v20}"                         # ROOT v2.0 (etap rucio_v20), tylko *dpt3to150*
: "${NOMINAL:=$RUNS}"                                         # models/D0_H (tylko odczyt) — ten sam model co EV/B2
: "${BT_SRC:=$Q/ev_dstar/src}"                                # charmtag 2.0.2 (jak EV)
: "${BT_OUT:=$EOS_OUT/charm_q/b_trig}"                        # parts/, bt_<próbka>.npz, trg_<próbka>.npz, DONE
: "${BT_SAMPLES:=data_428770 data_427929 fd_nd mb_nd}"        # ścieżka B (h5_v20)
: "${BT_TRG:=data_428770 data_427929}"                        # ścieżka A (ROOT, tabela zdarzeń) — nie czeka na h5_v20
: "${BT_GROUP:=10}"                                           # plików h5 na węzeł ekstrakcji
: "${BT_FRAC_DATA:=0.3}"; : "${BT_FRAC_MB:=0.3}"               # podpróbka zdarzeń dla all_* (ε_B); pary zawsze 100 %
: "${BT_LXY:=0.25}"                                           # L_xy(D0, D*) ≥ 0,25 mm → zgodność z v1.2 (h5_v20 ma 0,10)
: "${BT_PRESEL:=none}"
export BT Q H5V20 RAWV20 NOMINAL BT_SRC BT_OUT BT_SAMPLES BT_TRG BT_GROUP BT_FRAC_DATA BT_FRAC_MB BT_LXY BT_PRESEL
bt_tag() { case $1 in data_*) echo data ;; fd_*) echo fd ;; mb_*) echo mb ;; *) echo "$1" ;; esac; }
bt_frac() { case $1 in data_*) echo "$BT_FRAC_DATA" ;; mb_*) echo "$BT_FRAC_MB" ;; *) echo 1.0 ;; esac; }
bt_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2" && rm -f "$1"; fi; }
bt_env() {
  q_env
  export PYTHONPATH="$BT:$BT_SRC:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1
  python3 -c "import charmtag, charmtag.kinematics, charmtag.model, uproot, h5py, btrig_common; print('charmtag:', charmtag.__file__)" \
    || { echo "BŁĄD: import charmtag/uproot/h5py/btrig_common"; return 1; }
}
