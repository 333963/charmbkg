#!/bin/bash
# submit_tg.sh check|submit|status — TAGGED (1 job). Submit w czystym env (bez LCG). Stary wynik odkładany, nie kasowany.
set -uo pipefail
TG=$(cd "$(dirname "$0")" && pwd); export TG
source "$TG/tg_common.sh"
mkdir -p "$Q/logs/tagged"
check() { local ok=0 n
  for s in data_427929 data_428770; do n=$(ls "$TG_APP"/app_${s}_*.npz 2>/dev/null | wc -l); echo "app $s: $n"; [ "$n" -gt 0 ] || { echo "STOP: brak app $s"; ok=1; }
    n=$(ls "$TG_H5/$s"/${s}_*_none.h5 2>/dev/null | wc -l); echo "h5 $s: $n"; [ "$n" -gt 0 ] || { echo "STOP: brak h5 $s"; ok=1; }; done
  [ -f "$TG_WPS" ] || { echo "STOP: brak $TG_WPS"; ok=1; }
  condor_q -batch 2>/dev/null | grep -q charm_tagged && { echo "STOP: charm_tagged w kolejce"; ok=1; }
  echo "AFS: $(fs lq ~ | awk 'NR==2{print $4}')"; [ $ok = 0 ] && echo "OK — można wysłać"; return $ok; }
case ${1:-check} in
  check) check ;;
  submit) check || exit 1
    if [ -e "$TG_OUT/DONE" ]; then old="$TG_OUT/old_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$old"; for f in "$TG_OUT"/DONE "$TG_OUT"/*.npz "$TG_OUT"/tg.log; do [ -e "$f" ] && mv "$f" "$old/"; done; fi
    cd "$TG" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit tg.sub ;;
  status) condor_q -batch 2>/dev/null | grep charm_tagged; ls -la "$TG_OUT" 2>/dev/null ;;
esac
