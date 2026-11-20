#!/bin/bash
# submit_dn2.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/dnet2)
#  check  — pliki paczki i ../dnet1, wyniki etapu 1 ($DN_OUT/data/dn_<próbka>.npz, meta odniesienia), ROOT v2.0 per próbka, pula GPU (informacyjnie),
#           nic w kolejce (charm_dnet2), brak $DN2_OUT/DONE, AFS < 85 %
#  submit — check + make_dag_dn2.sh → condor_submit_dag dn2.dag (batch charm_dnet2) w CZYSTYM env; DRYRUN=1 → tylko DAG
set -uo pipefail
DN2=$(cd "$(dirname "$0")" && pwd); export DN2
source "$DN2/dn2_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n q f
  for f in dn2_extract.py dn2_train.py dn2_models.py dn2_job.sh dn2.sub; do [ -f "$DN2/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for f in dn_feat.py dn_train.py dn_models.py analyze.py; do [ -f "$DN1/$f" ] || { echo "STOP: brak ../dnet1/$f"; ok=1; }; done
  for s in $DN_SAMPLES; do
    [ -s "$DN_OUT/data/dn_$s.npz" ] || { echo "STOP: brak $DN_OUT/data/dn_$s.npz (etap dnet1)"; ok=1; }
    n=$(ls "$RAWV20/$s/"*.root 2>/dev/null | wc -l); echo "  ROOT v2.0 $s: $n plików"; [ "$n" -gt 0 ] || { echo "STOP: brak ROOT v2.0 $s"; ok=1; }
  done
  [ -s "$DN_OUT/meta/D0.xgb.H.sb.n.s1.json" ] || { echo "STOP: brak score odniesienia XGB H z etapu 1"; ok=1; }
  if [ "$DN2_GPU" = 1 ]; then
    echo "  GPU w puli: $(condor_status -constraint 'TotalGpus > 0' -af State 2>/dev/null | sort | uniq -c | tr '\n' ' ')"
    ls -d /cvmfs/sft.cern.ch/lcg/views/LCG_1[01][0-9]*_cuda/x86_64-el9-gcc*-opt >/dev/null 2>&1 || { echo "STOP: brak widoku LCG z CUDA (DN2_GPU=0 → CPU)"; ok=1; }
  fi
  condor_q -batch 2>/dev/null | grep -qE 'charm_dnet2|dn2\.dag' && { echo "STOP: DNET2 w kolejce"; ok=1; }
  [ -s "$DN2_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$DN2_OUT/DONE")); FORCE=1 albo inny DN2_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  DN2_OUT=$DN2_OUT GPU=$DN2_GPU epok $DN2_EPOCHS bloki: $DN2_BLOCKS"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_dn2.sh submit)" ;;
  submit) [ -z "${DN2_DAG:-}" ] || { echo "STOP: DN2_DAG tylko do testów"; exit 1; }
          check || exit 1
          "$DN2/make_dag_dn2.sh" || exit 1
          chmod +x "$DN2"/*.sh "$DN2"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -c '^JOB' "$DN2/dn2.dag"; grep -h '^VARS train' "$DN2/dn2.dag" | head -1; exit 0; fi
          cd "$DN2" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_dnet2 dn2.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_dnet2|Total for klyczak'; ls "$DN2_OUT" 2>/dev/null | head; echo "parts $(ls "$DN2_OUT/parts" 2>/dev/null | wc -l) meta $(ls "$DN2_OUT/meta" 2>/dev/null | wc -l)" ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
