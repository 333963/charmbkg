#!/bin/bash
# submit_tj2.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/trkjet2)
#  check  — wejścia: bt_<próbka>.npz (b_trig DONE), h5 v2.0 i ROOT v2.0 (dpt3to150) per próbka; nic w kolejce (charm_trkjet2); brak $TJ2_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_tj2.sh → condor_submit_dag tj2.dag (batch charm_trkjet2) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
TJ2=$(cd "$(dirname "$0")" && pwd); export TJ2
source "$TJ2/tj2_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s q n r
  for f in tj2_extract.py tj2_kernel.py tj2_merge.py; do [ -f "$TJ2/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -s "$BT_IN/DONE" ] || { echo "STOP: brak $BT_IN/DONE (etap b_trig)"; ok=1; }
  for s in $TJ2_SAMPLES; do
    n=$(ls "$H5V20/$s" 2>/dev/null | grep -c '_none\.h5$'); r=$(ls "$RAWV20/$s" 2>/dev/null | grep -c 'dpt3to150.*\.root$')
    echo "  $s: bt $( [ -s "$BT_IN/bt_$s.npz" ] && echo jest || echo BRAK ), h5 $n, ROOT $r, węzłów $(tj2_np "$s")"
    [ -s "$BT_IN/bt_$s.npz" ] && [ "$n" -gt 0 ] && [ "$r" -ge "$n" ] || { echo "STOP: $s — brak wejścia (bt/h5/ROOT; ROOT < h5?)"; ok=1; }
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_trkjet2|tj2\.dag' && { echo "STOP: TRKJET2 w kolejce"; ok=1; }
  [ -s "$TJ2_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: TRKJET2 już policzone ($(cat "$TJ2_OUT/DONE")); FORCE=1 nadpisuje albo inny TJ2_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  TJ2_OUT=$TJ2_OUT BT_IN=$BT_IN TJ2_NPARTS=$TJ2_NPARTS"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_tj2.sh submit)" ;;
  submit) [ -z "${TJ2_DAG:-}" ] || { echo "STOP: TJ2_DAG tylko do testów make_dag_tj2.sh"; exit 1; }
          check || exit 1
          "$TJ2/make_dag_tj2.sh" || exit 1
          chmod +x "$TJ2"/*.sh "$TJ2"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then for s in $TJ2_SAMPLES; do DRYRUN=1 "$TJ2/tj2_job.sh" ext "$s" 0; done; DRYRUN=1 "$TJ2/tj2_job.sh" merge; exit 0; fi
          ( tj2_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$TJ2" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_trkjet2 tj2.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_trkjet2|Total for klyczak'; ls -la "$TJ2_OUT" 2>/dev/null | head -20; ls "$TJ2_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$TJ2_OUT" -maxdepth 1 \( -name 'tj2_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
