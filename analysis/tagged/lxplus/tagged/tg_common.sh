# tg_common.sh — paczka TAGGED (19.11.2026, projekt B): masa × score × tag AFP per kandydat (app etapu 4 + evtable h5 v2.0). source, nie uruchamiać.
TG=${TG:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$TG/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                      # EOS_OUT, q_env, put_output
: "${TG_APP:=$EOS_OUT/charm_q/dnet4/app/d150}"             # wejście (tylko odczyt)
: "${TG_H5:=$EOS_OUT/charm_h5_v20/parts}"                  # wejście (tylko odczyt): <próbka>/<próbka>_<NNNNNN>_none.h5 (evtable)
: "${TG_WPS:=$EOS_OUT/charm_q/wpscan/wps_hist.npz}"        # progi WP drzewa H
: "${TG_OUT:=$EOS_OUT/charm_q/tagged}"                     # wynik: tg_cand_<run>.npz, tg_incl_<run>.npz, tg.log, DONE
export TG Q TG_APP TG_H5 TG_WPS TG_OUT
