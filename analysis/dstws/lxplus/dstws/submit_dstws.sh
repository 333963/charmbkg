#!/bin/bash
# submit_dstws.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/dstws)
#  check  — wejścia: h5 v2.0 i ROOT v2.0 (dpt3to150) per próbka, model D0_H, charmtag (ev_dstar/src); nic w kolejce (charm_dstws);
#           brak $DW_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_dstws.sh → condor_submit_dag dstws.dag (batch charm_dstws) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
DW=$(cd "$(dirname "$0")" && pwd); export DW
source "$DW/dstws_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s q n r
  for f in dstws_extract.py dstws_kernel.py dstws_merge.py dstws_job.sh dstws.sub; do [ -f "$DW/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -s "$NOMINAL/models/D0_H/config_used.yaml" ] || { echo "STOP: brak modelu $NOMINAL/models/D0_H"; ok=1; }
  [ -f "$DW_SRC/charmtag/features.py" ] || { echo "STOP: brak charmtag w $DW_SRC"; ok=1; }
  for s in $DW_SAMPLES; do
    n=$(ls "$H5V20/$s" 2>/dev/null | grep -c '_none\.h5$'); r=$(ls "$RAWV20/$s" 2>/dev/null | grep -c 'dpt3to150.*\.root$')
    echo "  $s: h5 $n, ROOT $r, węzłów $(dw_np "$s")"
    [ "$n" -gt 0 ] && [ "$r" -ge "$n" ] || { echo "STOP: $s — brak wejścia (h5/ROOT; ROOT < h5?)"; ok=1; }
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_dstws|dstws\.dag' && { echo "STOP: DSTWS w kolejce"; ok=1; }
  [ -s "$DW_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: DSTWS już policzone ($(cat "$DW_OUT/DONE")); FORCE=1 nadpisuje albo inny DW_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  DW_OUT=$DW_OUT DW_NPARTS=$DW_NPARTS"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_dstws.sh submit)" ;;
  submit) [ -z "${DW_DAG:-}" ] || { echo "STOP: DW_DAG tylko do testów make_dag_dstws.sh"; exit 1; }
          check || exit 1
          "$DW/make_dag_dstws.sh" || exit 1
          chmod +x "$DW"/*.sh "$DW"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then for s in $DW_SAMPLES; do DRYRUN=1 "$DW/dstws_job.sh" ext "$s" 0; done; DRYRUN=1 "$DW/dstws_job.sh" merge; exit 0; fi
          ( dw_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$DW" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_dstws dstws.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_dstws|Total for klyczak'; ls -la "$DW_OUT" 2>/dev/null | head -20; ls "$DW_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$DW_OUT" -maxdepth 1 \( -name 'dw_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
