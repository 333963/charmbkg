# pl_common.sh — paczka PEAK_LXY010 (18.11.2026, projekt B): D0 z tagiem D* (obie hipotezy K/π) z h5 **v2.0** (L_xy ≥ 0,10 mm)
# → parquet per próbka, do powtórki PEAK_pos / SWAP_d0 przy selekcji głównego łańcucha (L_xy ≥ 0,10 mm). Wzorzec: Q_lxplus/hyp_extract (etap hyp_dstar).
# source, nie uruchamiać. Wymaga ../common.sh (Q).
PL=${PL:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$PL/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                   # EOS_OUT, RUNS, q_env, put_output
: "${H5V20:=/eos/user/k/klyczak/charm_h5_v20/parts}"     # wejście: {fd_nd,mb_nd,data_428770,data_427929}/*_none.h5 (etap h5_v20; D0/D* z L_xy ≥ 0,10)
: "${NOMINAL:=${RUNS:-$EOS_OUT/charm_q/runs/lxplus_v11}}" # nominalny bieg: models/D0_H (fold*.json, meta.json) — tylko odczyt
: "${PL_OUT:=$EOS_OUT/charm_q/peak_lxy010}"             # wyjścia (EOS): <próbka>.parquet + .meta.json, extract_summary.json, thresholds.json, DONE
: "${SAMPLES:=fd_nd mb_nd data_428770 data_427929}"
[ -n "${PL_SAMPLES_CSV:-}" ] && SAMPLES=${PL_SAMPLES_CSV//,/ }   # węzły condora: lista z pl.dag (bez spacji w environment)
: "${DATA_FRAC:=1}"                                      # pełna statystyka (wiersze z tagiem D* to ~1–4 % okna)
: "${PL_EXTRA:=--dstar-only}"                            # tylko D0 z tagiem D*; pusty = wszystkie D0 w oknie (duże! ustaw wtedy DATA_FRAC/MAX_GB)
: "${MAX_GB:=0}"                                         # limit parquet [GB]; 0 = bez limitu (z --dstar-only nie działa i nie jest potrzebny)
: "${SCORE_MODE:=model}"                                 # v2.0 ma inne wiersze niż scores/D0_H (v1.1) → wynik z modeli foldów (fold po (run, event) jak w apply)
export PL Q H5V20 NOMINAL PL_OUT SAMPLES DATA_FRAC PL_EXTRA MAX_GB SCORE_MODE
pl_args() {   # argumenty ekstrakcji dla próbki $1
  local s=$1 frac=1
  case $s in data_*) frac=$DATA_FRAC ;; esac
  echo "--sample $s --h5 $H5V20/$s/*_none.h5 --models $NOMINAL/models/D0_H --score-mode $SCORE_MODE --frac $frac --max-gb $MAX_GB${PL_EXTRA:+ $PL_EXTRA}"
}
pl_put() {    # src dst: put_output z common.sh (gdy jest), inaczej cp + rm
  if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2" && rm -f "$1"; fi
}
pl_env() {    # LCG + venv (q_env) i NASZ src (snapshot Q_lxplus/hyp_extract/src = AA_merge + HYP.diff + tag D*) przed zainstalowanym charmtag
  q_env
  export PYTHONPATH="$PL/src:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1
  python3 - <<'PY' || { echo "BŁĄD: charmtag nie z $PL/src, bez hyp_extract/dstar_tag albo brak pyarrow"; return 1; }
import os, charmtag, pyarrow  # noqa: F401
from charmtag import hyp, hyp_extract, dstar_tag  # noqa: F401
assert charmtag.__file__.startswith(os.environ["PL"] + "/src"), charmtag.__file__
PY
}
