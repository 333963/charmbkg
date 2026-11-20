# dn_common.sh — paczka DNET etap 1 (17.11.2026, projekt B): XGB H vs sieci (MLP + dekorelacja od masy) dla D0→Kπ i D+→Kππ. source, nie uruchamiać.
# Na lxplus: ~/private/charm_q/Q_lxplus/dnet1 (= lokalnie charm_work/Q_lxplus/dnet1). Wymaga ../common.sh. Tylko pojedyncze D (bez par D–D).
DN=${DN:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$DN/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                    # EOS_OUT, q_env, put_output
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                # wejście (odczyt): <próbka>/<próbka>_NNNNNN_none.h5
: "${DN_OUT:=$EOS_OUT/charm_q/dnet}"                     # wyjście: parts/, data/dn_<próbka>.npz, sc/, meta/, onnx/, dnet_results.json, DONE
: "${DN_SAMPLES:=data_428770 data_427929 fd_nd mb_nd}"
: "${DN_GROUPS:=data_428770:10 data_427929:60 mb_nd:20 fd_nd:12}"   # plików h5 na węzeł ekstrakcji → 10 + 8 + 4 + 2 = 24 węzły
: "${DN_FRAC_D0:=0.10}"; : "${DN_FRAC_DP:=0.03}"        # podpróbka ZDARZEŃ danych w tabelach inkluzywnych (pary D0←D* zawsze 100 %)
: "${DN_FRAC_DP_MB:=0.3}"                                # podpróbka zdarzeń mb dla D+ (D0: 1,0)
: "${DN_LXY:=0.25}"                                      # L_xy ≥ [mm] (jak B_trig / dlz_data na v2.0)
: "${DN_FOLDS:=5}"; : "${DN_EPOCHS:=30}"; : "${DN_CPUS:=4}"
: "${DN_CAP_SIG:=400000}"; : "${DN_CAP_BKG:=800000}"
: "${DN_LAM:=5,25}"; : "${DN_LAM_SCAN:=1,5,25,100}"    # λ DisCo: nominalne (3 seedy) i skan
: "${DN_BLOCKS:=xgbA xgbB mlpB1 mlpB2 mlpC1 mlpC2 mlpC3 mlpC4 dpA dpC}"
: "${DN_NBOOT:=200}"; : "${DN_FAST:=0}"
: "${DN_PRESEL:=none}"
export DN Q H5V20 DN_OUT DN_SAMPLES DN_GROUPS DN_FRAC_D0 DN_FRAC_DP DN_FRAC_DP_MB DN_LXY DN_FOLDS DN_EPOCHS DN_CPUS DN_CAP_SIG DN_CAP_BKG DN_LAM \
       DN_LAM_SCAN DN_BLOCKS DN_NBOOT DN_FAST DN_PRESEL
dn_tag() { case $1 in data_*) echo data ;; fd_*) echo fd ;; mb_*) echo mb ;; *) echo "$1" ;; esac; }
dn_group() { local x; for x in $DN_GROUPS; do [ "${x%%:*}" = "$1" ] && { echo "${x##*:}"; return; }; done; echo 20; }
dn_put() { if declare -F put_output >/dev/null && [ "${DN_NOXRD:-0}" != 1 ]; then put_output "$1" "$2"; else mkdir -p "$(dirname "$2")" && cp -f "$1" "$2"; fi; }
dn_env() {
  q_env
  export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=$DN_CPUS MKL_NUM_THREADS=$DN_CPUS
  python3 -c 'import numpy, h5py, scipy, pandas, yaml, xgboost, torch, iminuit; print("torch", torch.__version__, "xgboost", xgboost.__version__)' \
    || { echo "BŁĄD: brak modułów w q_env"; return 1; }
}
