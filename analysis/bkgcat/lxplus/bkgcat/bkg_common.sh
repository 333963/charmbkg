# bkg_common.sh — paczka BKGCAT (20.11.2026, projekt B): histogramy masa × score drzewa H × kategoria prawdy dla mb_nd + tablica rozpadów z prawdy generatora.
# source, nie uruchamiać. Wzorzec: ../wpscan. POST HOC względem FITBIAS. Bez AFP, tylko pojedyncze D, tylko MC.
BKG=${BKG:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$BKG/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                      # EOS_OUT, q_env, put_output
: "${BKG_APP:=$EOS_OUT/charm_q/dnet4/app/d150}"            # wejście (odczyt): app_mb_nd_<k>_<NNNNNN>.npz (score per kandydat, etap 4)
: "${BKG_H5:=$EOS_OUT/charm_h5_v20/parts/mb_nd}"           # wejście (odczyt): mb_nd_<NNNNNN>_none.h5 (gałęzie prawdy kandydata)
: "${BKG_RAW:=$EOS_OUT/charm_raw_v20/mb_nd}"               # wejście (odczyt): ROOT v2.0 (TruthDmesons, DMesons_truthPt)
: "${BKG_THR:=$EOS_OUT/charm_q/wpscan/wps_hist.npz}"       # progi punktów pracy = WPSCAN (bez przeliczania)
: "${BKG_OUT:=$EOS_OUT/charm_q/bkgcat}"                    # wynik: parts/bkg_<NNNNNN>.npz, bkgcat.npz, bkg.log, DONE
: "${BKG_NGRP:=8}"
export BKG Q BKG_APP BKG_H5 BKG_RAW BKG_THR BKG_OUT BKG_NGRP
