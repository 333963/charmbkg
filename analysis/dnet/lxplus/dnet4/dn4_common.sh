# dn4_common.sh — paczka DNET etap 4 (19.11.2026, projekt B): własna sieć jako kandydat na nominał. source, nie uruchamiać.
#   (1) metoda R: względny SF sieć/drzewo z inkluzywnego piku D0/D+ na PEŁNEJ statystyce danych (modele foldów zapisane → węzły app),
#   (2) ziarna 1–3 D0.tr.H+e04.mc.d5, (3) sieć D+ (3 ziarna), (4) tło treningowe z pasm danych obu runów (sb2mc), (5) bin 2–3 GeV (dpt2to3).
# Na lxplus: ~/private/charm_q/Q_lxplus/dnet4. Wymaga ../dnet2 i ../dnet1 (kod wspólny, bez zmian) oraz wyników na EOS:
#   $DN_OUT/data (tabele etapu 1), $DN2_OUT/parts (tory otoczenia D0), $DN_OUT/{sc,meta}. Kryteria: DNET/KRYTERIA_etap4.md. Bez GPU.
DN4=${DN4:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$DN4/.." && pwd)}
DN2=${DN2:-$Q/dnet2}
# shellcheck disable=SC1091
source "$DN2/dn2_common.sh"                              # DN1, DN_OUT, DN2_OUT, H5V20, RAWV20, DN_FOLDS, DN_LAM, dn_put, dn2_env, dn2_kv, q_env
: "${DN4_OUT:=$EOS_OUT/charm_q/dnet4}"                   # wyjście: parts_dp/, sc/, meta/, models/, app/{d150,d23}/, *.json, DONE, DONE23, fail/
: "${H5D23:=$EOS_OUT/charm_h5_v20_d23/parts}"            # h5 v2.0 produkcji dpt2to3 (etap h5_d23)
: "${RAWD23:=$EOS_OUT/charm_raw_v20_d23}"                # ROOT v2.0 dpt2to3 (drzewo Tracks)
: "${DN4_BLOCKS:=xgbD0 xgbDp d0s1 d0s2 d0s3 dps1 dps2 dps3 d0sb}"
: "${DN4_SOFT:=d0sb dps2 dps3}"                          # bloki „miękkie”: niepowodzenie nie zatrzymuje DAG (znacznik fail/<węzeł> na EOS)
: "${DN4_EXTDP:=data_428770:10 data_427929:8 mb_nd:8 fd_nd:2}"        # węzły ekstrakcji torów otoczenia D+ (tabele etapu 1, 100 % wierszy)
: "${DN4_APP150:=data_428770:20 data_427929:20 mb_nd:8 fd_nd:2}"      # węzły app pT 3–30 GeV (h5 → cechy → tory → score; 100 % zdarzeń)
: "${DN4_APP23:=data_428770:32 data_427929:30 fd_nd:3}"               # węzły app pT 2–3 GeV; pusty = bez binu 2–3
: "${DN4_TAGS23:=\\.xgb\\.|\\.mc\\.d[0-9]+\\.s1\$}"       # modele stosowane w binie 2–3: drzewa odniesienia + sieci ziarna 1
: "${DN4_MR_D0:=1.70,2.05}"                              # zakres mas D0 w węzłach app (fit nominalny 1,75–2,00; szerszy = wariant zakresu)
: "${DN4_EPOCHS:=20}"; : "${DN4_CAP_SIG:=300000}"; : "${DN4_CAP_BKG:=500000}"; : "${DN4_CPUS:=8}"; : "${DN4_APP_CPUS:=4}"
: "${DN4_NBOOT:=100}"; : "${DN4_SMOKE:=1}"
export DN4 DN2 DN1 DN Q DN4_OUT H5D23 RAWD23 DN4_BLOCKS DN4_SOFT DN4_EXTDP DN4_APP150 DN4_APP23 DN4_TAGS23 DN4_MR_D0 DN4_EPOCHS DN4_CAP_SIG DN4_CAP_BKG \
       DN4_CPUS DN4_APP_CPUS DN4_NBOOT DN4_SMOKE
dn4_env() {   # CPU: q_env (LCG_110_swan + venv) przez dn2_env bez GPU
  DN2_USE_GPU=0; DN2_CPUS=${1:-$DN4_CPUS}; export DN2_USE_GPU DN2_CPUS
  dn2_env
  python3 -c 'import iminuit, xgboost, h5py, scipy; print("iminuit", iminuit.__version__, "xgboost", xgboost.__version__)' || { echo "BŁĄD: środowisko (iminuit/xgboost/h5py)"; return 1; }
}
dn4_h5dir()  { case $1 in d23) echo "$H5D23" ;; *) echo "$H5V20" ;; esac; }
dn4_rawdir() { case $1 in d23) echo "$RAWD23" ;; *) echo "$RAWV20" ;; esac; }
dn4_files()  { ls "$(dn4_h5dir "$1")/$2/"*_none.h5 2>/dev/null | sort; }     # dn4_files <tryb> <próbka>
