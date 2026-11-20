# cn_common.sh — paczka CNCLONE (17.11; projekt B: klon CharmNet D0 z izolacją ptcone40/100 i bez → SF z tagu D*; source, nie uruchamiać).
# Na lxplus: ~/private/charm_q/Q_lxplus/cnclone (= lokalnie charm_work/Q_lxplus/cnclone). Wymaga ../common.sh. TO KLON, nie CharmNet.
CN=${CN:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$CN/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                    # EOS_OUT, q_env, put_output
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                # wejście (odczyt): <próbka>/<próbka>_NNNNNN_none.h5
: "${CN_OUT:=$EOS_OUT/charm_q/cnclone}"                  # wyjście: parts/cn_<próbka>_<g>.npz → cn_<próbka>.npz, cn_models.npz, cn_summary.json, DONE
: "${CN_SAMPLES:=data_428770 data_427929 fd_nd mb_nd}"
: "${CN_GROUP:=10}"                                      # plików h5 na węzeł ekstrakcji
: "${CN_FRAC_DATA:=0.3}"; : "${CN_FRAC_MB:=1.0}"        # podpróbka zdarzeń tabel inkluzywnych a_* (pary zawsze 100 %)
: "${CN_LXY:=0.10}"                                      # L_xy(D0), L_xy(D*) ≥ [mm]; 0,10 = podłoga plików = selekcja treningu CharmNet (B_trig: 0,25)
: "${CN_FOLDS:=5}"; : "${CN_EPOCHS:=20}"; : "${CN_CPUS:=4}"
: "${CN_MAX_SIG:=1500000}"; : "${CN_MAX_BKG:=3000000}"
: "${CN_FD_BKG:=0}"                                      # 1 = tło treningu także z fd_nd (przepis: cała reszta MC)
: "${CN_PRESEL:=none}"
export CN Q H5V20 CN_OUT CN_SAMPLES CN_GROUP CN_FRAC_DATA CN_FRAC_MB CN_LXY CN_FOLDS CN_EPOCHS CN_CPUS CN_MAX_SIG CN_MAX_BKG CN_FD_BKG CN_PRESEL
cn_tag() { case $1 in data_*) echo data ;; fd_*) echo fd ;; mb_*) echo mb ;; *) echo "$1" ;; esac; }
cn_frac() { case $1 in data_*) echo "$CN_FRAC_DATA" ;; mb_*) echo "$CN_FRAC_MB" ;; *) echo 1.0 ;; esac; }
cn_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2"; fi; }
cn_env() {
  q_env
  export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=$CN_CPUS MKL_NUM_THREADS=$CN_CPUS
  python3 -c 'import numpy, h5py, scipy' || { echo "BŁĄD: brak numpy/h5py/scipy w q_env"; return 1; }
  python3 -c 'import torch; print("torch", torch.__version__)' 2>/dev/null || echo "UWAGA: brak torch → zamiennik sklearn HGB (flaga substitute w cn_summary.json)"
}
