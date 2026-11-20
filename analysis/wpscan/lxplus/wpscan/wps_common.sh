# wps_common.sh — paczka WPSCAN (19.11.2026, projekt B): histogramy masa × score z plików app etapu 4 dla skanu punktu pracy. source, nie uruchamiać.
WPS=${WPS:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$WPS/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                      # EOS_OUT, q_env, put_output
: "${WPS_APP:=$EOS_OUT/charm_q/dnet4/app/d150}"            # wejście (tylko odczyt): app_<próbka>_*.npz etapu 4
: "${WPS_OUT:=$EOS_OUT/charm_q/wpscan}"                    # wynik: wps_hist.npz, wps.log, DONE
: "${WPS_NREP:=30}"
export WPS Q WPS_APP WPS_OUT WPS_NREP
