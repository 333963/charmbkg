#!/bin/bash
# submit_dn4.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/dnet4)
#  check  — pliki paczki, ../dnet2 i ../dnet1; wejścia na EOS: tabele etapu 1, części z torami D0 etapu 2 (komplet), h5/ROOT v2.0 (dpt3to150 i dpt2to3),
#           nazwy bloków; nic w kolejce (charm_dnet4); brak $DN4_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_dn4.sh → condor_submit_dag dn4.dag (batch charm_dnet4) w CZYSTYM env; DRYRUN=1 → tylko DAG
#           FORCE=1 przy istniejącym DONE: stare DONE/DONE23 → *.prev_<epoch> (nic nie jest kasowane); gotowe modele/części są pomijane (DN_REDO=1 liczy od nowa)
set -uo pipefail
DN4=$(cd "$(dirname "$0")" && pwd); export DN4
source "$DN4/dn4_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n q f b x kv
  for f in dn4_train.py dn4_lib.py dn4_apply.py dn4_extract.py dn4_ana.py dn4_pis.py dn4_job.sh dn4.sub make_dag_dn4.sh; do [ -f "$DN4/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for f in dn2_train.py dn2_models.py dn2_extract.py dn2_common.sh; do [ -f "$DN2/$f" ] || { echo "STOP: brak ../dnet2/$f"; ok=1; }; done
  for f in dn_feat.py dn_train.py dn_extract.py dn_models.py analyze.py btrig_common.py dn_common.sh; do [ -f "$DN1/$f" ] || { echo "STOP: brak ../dnet1/$f"; ok=1; }; done
  for b in $DN4_BLOCKS; do grep -q "\"$b\": \[" "$DN4/dn4_train.py" || { echo "STOP: nieznany blok $b (dn4_train.py BLOCKS)"; ok=1; }; done
  for s in $DN_SAMPLES; do
    [ -s "$DN_OUT/data/dn_$s.npz" ] || { echo "STOP: brak $DN_OUT/data/dn_$s.npz (etap dnet1)"; ok=1; }
    x=$(dn2_kv "$DN2_NPARTS" "$s" 4)
    n=$(ls "$DN2_OUT/parts/dn2_${s}_"*.npz 2>/dev/null | wc -l); echo "  tory otoczenia D0 $s: $n części (oczekiwane $x)"
    [ "$n" -eq "$x" ] || { echo "STOP: niepełne części etapu 2 dla $s"; ok=1; }
  done
  for kv in $DN4_APP150; do
    s=${kv%%:*}; n=$(dn4_files d150 "$s" | wc -l); x=$(ls "$RAWV20/$s"/*.root 2>/dev/null | wc -l); echo "  d150 $s: h5 $n, ROOT $x, węzłów ${kv##*:}"
    [ "$n" -gt 0 ] && [ "$x" -gt 0 ] || { echo "STOP: brak h5/ROOT v2.0 dla $s"; ok=1; }
  done
  for kv in $DN4_APP23; do
    s=${kv%%:*}; n=$(dn4_files d23 "$s" | wc -l); x=$(ls "$RAWD23/$s"/*.root 2>/dev/null | wc -l); echo "  d23 $s: h5 $n, ROOT $x, węzłów ${kv##*:}"
    [ "$n" -gt 0 ] && [ "$x" -gt 0 ] || { echo "STOP: brak h5/ROOT dpt2to3 dla $s (DN4_APP23='' wyłącza bin 2–3)"; ok=1; }
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_dnet4|dn4\.dag' && { echo "STOP: DNET4 w kolejce"; ok=1; }
  [ -s "$DN4_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$DN4_OUT/DONE")); FORCE=1 albo inny DN4_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  DN4_OUT=$DN4_OUT smoke=$DN4_SMOKE epok $DN4_EPOCHS bloki: $DN4_BLOCKS (miękkie: $DN4_SOFT)"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_dn4.sh submit)" ;;
  submit) [ -z "${DN4_DAG:-}" ] || { echo "STOP: DN4_DAG tylko do testów"; exit 1; }
          check || exit 1
          "$DN4/make_dag_dn4.sh" || exit 1
          chmod +x "$DN4"/*.sh "$DN4"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -c '^JOB' "$DN4/dn4.dag"; grep -h '^PARENT' "$DN4/dn4.dag" | cut -c1-200; exit 0; fi
          if [ "${FORCE:-0}" = 1 ]; then for f in DONE DONE23; do [ -e "$DN4_OUT/$f" ] && { mv "$DN4_OUT/$f" "$DN4_OUT/$f.prev_$(date +%s)" || exit 1; }; done; fi
          cd "$DN4" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_dnet4 dn4.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_dnet4|Total for klyczak'; ls "$DN4_OUT" 2>/dev/null | head -20
          echo "sc $(ls "$DN4_OUT/sc" 2>/dev/null | wc -l) meta $(ls "$DN4_OUT/meta" 2>/dev/null | wc -l) parts_dp $(ls "$DN4_OUT/parts_dp" 2>/dev/null | wc -l) app150 $(ls "$DN4_OUT/app/d150" 2>/dev/null | wc -l) app23 $(ls "$DN4_OUT/app/d23" 2>/dev/null | wc -l) fail $(ls "$DN4_OUT/fail" 2>/dev/null | tr '\n' ' ')" ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
