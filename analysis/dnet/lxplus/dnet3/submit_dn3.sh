#!/bin/bash
# submit_dn3.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/dnet3)
#  check  — pliki paczki, ../dnet2 i ../dnet1; wejścia na EOS: tabele etapu 1, części z torami etapu 2 (komplet, klucz pc40), score odniesień
#           (XGB H z etapu 1, tr.H+e04 z etapu 2); nazwy bloków; pula GPU (informacyjnie); nic w kolejce (charm_dnet3); brak $DN3_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_dn3.sh → condor_submit_dag dn3.dag (batch charm_dnet3) w CZYSTYM env; DRYRUN=1 → tylko DAG
#           FORCE=1 przy istniejącym DONE: stary DONE → DONE.prev_<epoch> (nic nie jest kasowane); gotowe konfiguracje są pomijane (DN_REDO=1 liczy od nowa)
set -uo pipefail
DN3=$(cd "$(dirname "$0")" && pwd); export DN3
source "$DN3/dn3_common.sh"
cmd=${1:-check}
REF2="D0.tr.H+e04.mc.d${DN_LAM%%,*}.s1"
check() {
  local ok=0 s n q f b x exp=0 p
  for f in dn3_train.py dn3_models.py dn3_job.sh dn3.sub dn3_gpu.sub make_dag_dn3.sh; do [ -f "$DN3/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for f in dn2_train.py dn2_models.py dn2_common.sh; do [ -f "$DN2/$f" ] || { echo "STOP: brak ../dnet2/$f"; ok=1; }; done
  for f in dn_feat.py dn_train.py dn_models.py analyze.py btrig_common.py dn_common.sh; do [ -f "$DN1/$f" ] || { echo "STOP: brak ../dnet1/$f"; ok=1; }; done
  for b in $DN3_BLOCKS; do grep -q "^    \"$b\": \[" "$DN3/dn3_train.py" || { echo "STOP: nieznany blok $b (dn3_train.py BLOCKS)"; ok=1; }; done
  for s in $DN_SAMPLES; do
    [ -s "$DN_OUT/data/dn_$s.npz" ] || { echo "STOP: brak $DN_OUT/data/dn_$s.npz (etap dnet1)"; ok=1; }
    x=$(dn2_kv "$DN2_NPARTS" "$s" 4); exp=$((exp + x))
    n=$(ls "$DN2_OUT/parts/dn2_${s}_"*.npz 2>/dev/null | wc -l); echo "  tory otoczenia $s: $n części (oczekiwane $x)"
    [ "$n" -eq "$x" ] || { echo "STOP: niepełne części etapu 2 dla $s (inne wiersze niż w odniesieniu)"; ok=1; }
  done
  p=$(ls "$DN2_OUT/parts/dn2_fd_nd_"*.npz 2>/dev/null | head -1)
  if [ -n "$p" ] && command -v unzip >/dev/null 2>&1; then
    for x in idx E nE pc40; do unzip -l "$p" 2>/dev/null | grep -q "D0__a_$x.npy" || { echo "STOP: w $p brak D0__a_$x"; ok=1; }; done
  else
    echo "  (bez kontroli kluczy npz: brak unzip albo części fd_nd)"
  fi
  [ -s "$DN_OUT/meta/D0.xgb.H.sb.n.s1.json" ] && [ -s "$DN_OUT/sc/D0.xgb.H.sb.n.s1.npz" ] || { echo "STOP: brak score odniesienia XGB H z etapu 1"; ok=1; }
  [ -s "$DN2_OUT/meta/$REF2.json" ] && [ -s "$DN2_OUT/sc/$REF2.npz" ] || { echo "STOP: brak score odniesienia $REF2 z etapu 2 ($DN2_OUT/sc)"; ok=1; }
  if [ "$DN3_GPU" = 1 ]; then
    echo "  GPU w puli: $(condor_status -constraint 'TotalGpus > 0' -af State 2>/dev/null | sort | uniq -c | tr '\n' ' ')"
    ls -d /cvmfs/sft.cern.ch/lcg/views/LCG_1[01][0-9]*_cuda/x86_64-el9-gcc*-opt >/dev/null 2>&1 || { echo "STOP: brak widoku LCG z CUDA (DN3_GPU=0 → CPU)"; ok=1; }
  fi
  condor_q -batch 2>/dev/null | grep -qE 'charm_dnet3|dn3\.dag' && { echo "STOP: DNET3 w kolejce"; ok=1; }
  [ -s "$DN3_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$DN3_OUT/DONE")); FORCE=1 albo inny DN3_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  DN3_OUT=$DN3_OUT GPU=$DN3_GPU smoke=$DN3_SMOKE epok $DN3_EPOCHS bloki: $DN3_BLOCKS"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_dn3.sh submit)" ;;
  submit) [ -z "${DN3_DAG:-}" ] || { echo "STOP: DN3_DAG tylko do testów"; exit 1; }
          check || exit 1
          "$DN3/make_dag_dn3.sh" || exit 1
          chmod +x "$DN3"/*.sh "$DN3"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -c '^JOB' "$DN3/dn3.dag"; grep -h '^JOB train' "$DN3/dn3.dag" | head -1; exit 0; fi
          if [ "${FORCE:-0}" = 1 ] && [ -e "$DN3_OUT/DONE" ]; then mv "$DN3_OUT/DONE" "$DN3_OUT/DONE.prev_$(date +%s)" || exit 1; fi
          cd "$DN3" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_dnet3 dn3.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_dnet3|Total for klyczak'; ls "$DN3_OUT" 2>/dev/null | head; echo "sc $(ls "$DN3_OUT/sc" 2>/dev/null | wc -l) meta $(ls "$DN3_OUT/meta" 2>/dev/null | wc -l)" ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
