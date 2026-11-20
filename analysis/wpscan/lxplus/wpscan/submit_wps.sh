#!/bin/bash
# submit_wps.sh check|submit|status — WPSCAN (1 job). Submit w czystym env (bez LCG). Stary wynik odkładany, nie kasowany.
set -uo pipefail
WPS=$(cd "$(dirname "$0")" && pwd); export WPS
source "$WPS/wps_common.sh"
mkdir -p "$Q/logs/wpscan"
check() { local ok=0 n
  for s in fd_nd mb_nd data_427929 data_428770; do n=$(ls "$WPS_APP"/app_${s}_*.npz 2>/dev/null | wc -l); echo "app $s: $n"; [ "$n" -gt 0 ] || { echo "STOP: brak app $s"; ok=1; }; done
  condor_q -batch 2>/dev/null | grep -q charm_wpscan && { echo "STOP: charm_wpscan w kolejce"; ok=1; }
  echo "AFS: $(fs lq ~ | awk 'NR==2{print $4}')"; [ $ok = 0 ] && echo "OK — można wysłać"; return $ok; }
case ${1:-check} in
  check) check ;;
  submit) check || exit 1
    if [ -e "$WPS_OUT/DONE" ]; then old="$WPS_OUT/old_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$old"; for f in DONE wps_hist.npz wps.log; do [ -e "$WPS_OUT/$f" ] && mv "$WPS_OUT/$f" "$old/"; done; fi
    cd "$WPS" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit wps.sub ;;
  status) condor_q -batch 2>/dev/null | grep charm_wpscan; ls -la "$WPS_OUT" 2>/dev/null ;;
esac
