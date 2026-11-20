# cd_common.sh — paczka CAND_D23 (18.11.2026, projekt B): kandydaci D⁺ per kandydat (score H modeli foldów dlz_data, poza zakresem
# treningu) + liczniki/histogramy D⁰ i D* z h5 v2.0 produkcji dpt2to3 (D pT 2–3 GeV). source, nie uruchamiać. Wymaga ../common.sh (Q).
# Dane inkluzywne — żadnej selekcji AFP. TYLKO pojedyncze D; par D–D paczka nie czyta.
CD=${CD:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$CD/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                          # EOS_OUT, q_env, put_output
: "${H5D23:=$EOS_OUT/charm_h5_v20_d23/parts}"                   # wejście (etap h5_d23): <próbka>/<próbka>_NNNNNN_none.h5, DONE
: "${CD_OUT:=$EOS_OUT/charm_q/cand_d23}"                        # wyjście: dp_<s>.parquet, hist_<s>.npz, summary.json, DONE; części w parts/
: "${DP_MODELS:=$EOS_OUT/charm_q/dlz_data/train/models/H_s12345}"   # modele foldów D⁺ H (fold{0..4}.json, meta.json) — te same co dlz_data/SCAN
: "${DP_THR:=$EOS_OUT/charm_q/dlz_data/train/thr.json}"         # progi WP binu 3–5 GeV (tylko odniesienie w summary.json)
: "${SAMPLES:=data_428770 data_427929 fd_nd}"
: "${CD_NGMAX:=data_428770:12 data_427929:10 fd_nd:3}"          # maks. liczba węzłów ekstrakcji na próbkę (pliki dzielone równo)
: "${CD_NFILES:=data_428770:96 data_427929:450 fd_nd:15}"       # oczekiwane *_none.h5 (tylko ostrzeżenie)
: "${CD_PT_LO:=2}"; : "${CD_PT_HI:=3}"                          # okno pT [GeV) — obejście twardego pT ≥ 3: src/charmbkg/samples.py pt_range()
: "${CD_LXY:=0.25}"                                             # jak dlz_data (produkcja v2.0 ma dla D⁺ L_xy ≥ 0,5 → cięcie nieaktywne)
: "${CD_SC_FULL:=0.93}"; : "${CD_WIN:=0.06}"; : "${CD_PS_WIN:=0.05}"; : "${CD_PS_OUT:=0.01}"   # preskala wierszy (cd_extract.py)
: "${CD_MAX_GB:=1.8}"; : "${CD_HARD_GB:=2.0}"                   # budżet wyjścia (finish przerzedza po ps_u; > HARD = STOP bez DONE)
: "${CD_NJ:=4}"
: "${CD_LOGDIR:=$Q/logs/cand_d23}"
[ -n "${CD_SAMPLES_CSV:-}" ] && SAMPLES=${CD_SAMPLES_CSV//,/ }   # węzły condora (bez spacji w environment)
[ -n "${CD_NG_CSV:-}" ] && CD_NG=${CD_NG_CSV//,/ }               # faktyczna liczba części per próbka (z make_dag_cd.sh)
export CD Q H5D23 CD_OUT DP_MODELS DP_THR SAMPLES CD_NGMAX CD_NFILES CD_PT_LO CD_PT_HI CD_LXY CD_SC_FULL CD_WIN CD_PS_WIN CD_PS_OUT CD_MAX_GB CD_HARD_GB CD_NJ CD_LOGDIR
cd_kv() { local kv; for kv in $1; do [ "${kv%%:*}" = "$2" ] && { echo "${kv##*:}"; return; }; done; echo "${3:-0}"; }
cd_files() { ls "$H5D23/$1/"*_none.h5 2>/dev/null | grep -v '\.part$' | sort; }
cd_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else mkdir -p "$(dirname "$2")" && cp -f "$1" "$2"; fi; }
cd_xargs() { echo "--dp-models $DP_MODELS --pt-lo $CD_PT_LO --pt-hi $CD_PT_HI --lxy-min $CD_LXY --sc-full $CD_SC_FULL --win $CD_WIN --ps-win $CD_PS_WIN --ps-out $CD_PS_OUT --nthread $CD_NJ"; }
cd_env() {   # LCG + venv (q_env; CD_NOENV=1 = test lokalny) i NASZ src (kopia charmbkg z pt_range)
  [ "${CD_NOENV:-0}" = 1 ] || q_env
  export PYTHONPATH="$CD:$CD/src:${PYTHONPATH:-}" CHARMBKG_SRC="$CD/src" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=$CD_NJ
  python3 - <<'PY' || { echo "BŁĄD: charmbkg nie z $CD/src (bez pt_range) albo brak h5py/xgboost/pyarrow/yaml"; return 1; }
import os, charmbkg, h5py, xgboost, pyarrow, yaml  # noqa: F401
from charmbkg.samples import pt_range
assert charmbkg.__file__.startswith(os.environ["CD"] + "/src"), charmbkg.__file__
assert pt_range() == (3.0, 30.0) or os.environ.get("CHARMBKG_PT_RANGE")
PY
}
