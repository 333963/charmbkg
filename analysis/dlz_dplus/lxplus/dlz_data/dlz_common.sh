# dlz_common.sh — paczka DLZ_DATA (17.11.2026; projekt B: H+DLZ dla D⁺ na danych wg DLZ_dplus/PROCEDURA.md + przesiew nowych cech D⁰; source).
# Na lxplus: ~/private/charm_q/Q_lxplus/dlz_data (= lokalnie charm_work/Q_lxplus/dlz_data; kopiuje etap ORCH stages/dlz_data.sh). Wymaga ../common.sh.
# TYLKO pojedyncze D (D⁺, D⁰). Par D–D paczka nie czyta ani nie liczy.
DZ=${DZ:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$DZ/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                        # EOS_OUT, q_env, put_output
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                     # wejście D⁺ (odczyt): <próbka>/<próbka>_NNNNNN_none.h5
: "${H5V12:=$EOS_OUT/charm_h5_v12/parts}"                     # wejście D⁰ (odczyt; jak FEAT_charmnet)
: "${DLZ_OUT:=$EOS_OUT/charm_q/dlz_data}"                     # parts/, train/{models,oof,thr.json}, dlz_<próbka>.npz, ana/, d0/, DONE
: "${DLZ_DATA:=data_428770 data_427929}"                      # próbki danych do apply
: "${DLZ_NFILES:=data_428770:97 data_427929:458 mb_nd:78 fd_nd:24}"   # oczekiwane *_none.h5 (h5_v20, 17.11)
: "${DLZ_GROUP:=data_428770:5 data_427929:46 mb_nd:13 fd_nd:8 sb:10}" # plików h5 na węzeł (sb = data_428770 dla tła treningowego)
: "${DLZ_SBFRAC:=0.06}"                                       # podpróbka zdarzeń SB 428770 do treningu (~5·10⁶ kand. → cap 3·10⁶)
: "${DLZ_LXY:=0.25}"                                          # L_xy ≥ 0,25 mm = nominał v1.2 (h5_v20 ma 0,10)
: "${DLZ_CAP:=3000000}"; : "${DLZ_NJ:=8}"; : "${DLZ_NBOOT:=200}"
: "${DLZ_SETS:=H:12345 HDLZ:12345 Hdal:12345 H:23456 HDLZ:23456}"
: "${DLZ_FAST:=0}"                                            # 1 = test (40 drzew)
: "${DLZ_D0:=1}"                                              # 1 = węzeł d0 (przesiew cech D⁰, MC + SB danych)
: "${D0_FD:=all}"; : "${D0_MB:=40}"; : "${D0_DATA:=6}"        # liczba plików h5 v1.2 dla D⁰ (fd: wszystkie; mb: pierwsze 40; dane: 6 równomiernie)
: "${D0_VARIANTS:=H,H_seed2,H+ALL_noiso,H+d0prod,H+ip,H+time,H+NOKIN,H+kin2,H+kinall,H+pidbb,H+NOKIN_pid,H+ALLnew}"
: "${D0_VARIANTS_DATA:=H,H_seed2,H+NOKIN,H+kinall,H+ALLnew}"
: "${DLZ_PRESEL:=none}"
export DZ Q H5V20 H5V12 DLZ_OUT DLZ_DATA DLZ_NFILES DLZ_GROUP DLZ_SBFRAC DLZ_LXY DLZ_CAP DLZ_NJ DLZ_NBOOT DLZ_SETS DLZ_FAST DLZ_D0 D0_FD D0_MB D0_DATA D0_VARIANTS D0_VARIANTS_DATA DLZ_PRESEL
dz_kv() { local kv; for kv in $1; do [ "${kv%%:*}" = "$2" ] && { echo "${kv##*:}"; return; }; done; echo "${3:-0}"; }
dz_ngroups() { local n g; n=$(dz_kv "$DLZ_NFILES" "$1"); g=$(dz_kv "$DLZ_GROUP" "${2:-$1}" 10); echo $(( (n + g - 1) / g )); }
dz_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else mkdir -p "$(dirname "$2")" && cp -f "$1" "$2"; fi; }
dz_env() {
  [ "${DLZ_NOENV:-0}" = 1 ] || q_env
  export PYTHONPATH="$DZ:$DZ/src:${PYTHONPATH:-}" CHARMBKG_SRC="$DZ/src" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=$DLZ_NJ MKL_NUM_THREADS=$DLZ_NJ OPENBLAS_NUM_THREADS=$DLZ_NJ
  python3 -c "import charmbkg, h5py, yaml, xgboost, sklearn, pandas, iminuit; print('charmbkg:', charmbkg.__file__, 'xgb', xgboost.__version__, 'iminuit', iminuit.__version__)" \
    || { echo "BŁĄD: import charmbkg/h5py/xgboost/iminuit"; return 1; }
}
