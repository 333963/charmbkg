#!/bin/bash
# submit_tj.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/trkjet)
#  check  — wejścia: bt_<próbka>.npz (b_trig DONE), h5 v2.0 i ROOT v2.0 (dpt3to150) per próbka; nic w kolejce (charm_trkjet); brak $TJ_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_tj.sh → condor_submit_dag tj.dag (batch charm_trkjet) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
TJ=$(cd "$(dirname "$0")" && pwd); export TJ
source "$TJ/tj_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s q n r
  for f in tj_extract.py tj_kernel.py tj_merge.py; do [ -f "$TJ/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -s "$BT_IN/DONE" ] || { echo "STOP: brak $BT_IN/DONE (etap b_trig)"; ok=1; }
  for s in $TJ_SAMPLES; do
    n=$(ls "$H5V20/$s" 2>/dev/null | grep -c '_none\.h5$'); r=$(ls "$RAWV20/$s" 2>/dev/null | grep -c 'dpt3to150.*\.root$')
    echo "  $s: bt $( [ -s "$BT_IN/bt_$s.npz" ] && echo jest || echo BRAK ), h5 $n, ROOT $r, węzłów $(tj_np "$s")"
    [ -s "$BT_IN/bt_$s.npz" ] && [ "$n" -gt 0 ] && [ "$r" -ge "$n" ] || { echo "STOP: $s — brak wejścia (bt/h5/ROOT; ROOT < h5?)"; ok=1; }
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_trkjet|tj\.dag' && { echo "STOP: TRKJET w kolejce"; ok=1; }
  [ -s "$TJ_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: TRKJET już policzone ($(cat "$TJ_OUT/DONE")); FORCE=1 nadpisuje albo inny TJ_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  TJ_OUT=$TJ_OUT BT_IN=$BT_IN TJ_NPARTS=$TJ_NPARTS"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_tj.sh submit)" ;;
  submit) [ -z "${TJ_DAG:-}" ] || { echo "STOP: TJ_DAG tylko do testów make_dag_tj.sh"; exit 1; }
          check || exit 1
          "$TJ/make_dag_tj.sh" || exit 1
          chmod +x "$TJ"/*.sh "$TJ"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then for s in $TJ_SAMPLES; do DRYRUN=1 "$TJ/tj_job.sh" ext "$s" 0; done; DRYRUN=1 "$TJ/tj_job.sh" merge; exit 0; fi
          ( tj_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$TJ" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_trkjet tj.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_trkjet|Total for klyczak'; ls -la "$TJ_OUT" 2>/dev/null | head -20; ls "$TJ_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$TJ_OUT" -maxdepth 1 \( -name 'tj_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
