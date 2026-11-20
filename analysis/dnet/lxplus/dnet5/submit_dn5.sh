#!/bin/bash
# submit_dn5.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/dnet5)
#  check  — pliki paczki i ../dnet{4,2,1}; wejścia na EOS: tabele etapu 1, score drzewa H etapu 1, tory otoczenia D0 (etap 2, komplet), tory otoczenia D+
#           (parts_dp etapu 4, komplet = extdp skończony), score odniesienia etapu 2; nic w kolejce (charm_dnet5); brak $DN5_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_dn5.sh → condor_submit_dag -maxjobs $DN5_MAXJOBS (batch charm_dnet5) w CZYSTYM env; potem sonda GPU (batch charm_dnet5gpu,
#           osobny job, sam się usuwa po $DN5_GPU_WAIT s w kolejce). DRYRUN=1 → tylko DAG. FORCE=1 przy DONE: stare DONE → DONE.prev_<epoch>.
set -uo pipefail
DN5=$(cd "$(dirname "$0")" && pwd); export DN5
source "$DN5/dn5_common.sh"
cmd=${1:-check}
CLEAN=(env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin)
check() {
  local ok=0 s n q f b x kv sp need_dp=0
  for f in dn5_train.py dn5_lib.py dn5_models.py dn5_ana.py dn5_job.sh dn5.sub dn5_gpu.sub dn5_gtrain.sub make_dag_dn5.sh; do [ -f "$DN5/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for f in dn4_train.py dn4_lib.py dn4_common.sh; do [ -f "$DN4/$f" ] || { echo "STOP: brak ../dnet4/$f"; ok=1; }; done
  for f in dn2_train.py dn2_models.py dn2_common.sh; do [ -f "$DN2/$f" ] || { echo "STOP: brak ../dnet2/$f"; ok=1; }; done
  for f in dn_feat.py dn_train.py dn_models.py analyze.py btrig_common.py dn_common.sh; do [ -f "$DN1/$f" ] || { echo "STOP: brak ../dnet1/$f"; ok=1; }; done
  for b in $DN5_BLOCKS; do grep -q "\"$b\": _v(" "$DN5/dn5_lib.py" || { echo "STOP: nieznany blok $b (dn5_lib.py BLOCKS)"; ok=1; }; case $b in eDp*) need_dp=1 ;; esac; done
  for s in $DN_SAMPLES; do
    [ -s "$DN_OUT/data/dn_$s.npz" ] || { echo "STOP: brak $DN_OUT/data/dn_$s.npz (etap dnet1)"; ok=1; }
    x=$(dn2_kv "$DN2_NPARTS" "$s" 4)
    n=$(ls "$DN2_OUT/parts/dn2_${s}_"*.npz 2>/dev/null | wc -l); echo "  tory otoczenia D0 $s: $n części (oczekiwane $x)"
    [ "$n" -eq "$x" ] || { echo "STOP: niepełne części etapu 2 dla $s"; ok=1; }
    if [ "$need_dp" = 1 ]; then
      x=$(dn2_kv "$DN4_EXTDP" "$s" 0)
      n=$(ls "$DN4_OUT/parts_dp/dn2_${s}_"*.npz 2>/dev/null | wc -l); echo "  tory otoczenia D+ $s: $n części (oczekiwane $x)"
      [ "$n" -eq "$x" ] || { echo "STOP: niepełne parts_dp etapu 4 dla $s (extdp jeszcze liczy?)"; ok=1; }
    fi
  done
  if [ "$need_dp" = 1 ] && condor_q -nobatch -constraint 'JobBatchName == "charm_dnet4"' -af Args 2>/dev/null | grep -q ' extdp '; then echo "STOP: węzły extdp etapu 4 jeszcze w kolejce"; ok=1; fi
  for sp in D0 Dplus; do [ -s "$DN_OUT/sc/$sp.xgb.H.sb.n.s1.npz" ] || { echo "STOP: brak score drzewa H etapu 1 ($sp)"; ok=1; }; done
  for s in 1 2 3; do [ -s "$DN2_OUT/sc/D0.tr.H+e04.mc.d${DN_LAM%%,*}.s$s.npz" ] || { echo "STOP: brak score odniesienia etapu 2 (ziarno $s)"; ok=1; }; done
  condor_q -batch 2>/dev/null | grep -qE 'charm_dnet5( |$)|dn5\.dag' && { echo "STOP: DNET5 w kolejce"; ok=1; }
  [ -s "$DN5_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$DN5_OUT/DONE")); FORCE=1 albo inny DN5_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  DN5_OUT=$DN5_OUT epok $DN5_EPOCHS rdzeni/fold $DN5_CPUS bloki: $DN5_BLOCKS  GPU: $DN5_GPU (czekanie ≤ ${DN5_GPU_WAIT}s)"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_dn5.sh submit)" ;;
  submit) [ -z "${DN5_DAG:-}" ] || { echo "STOP: DN5_DAG tylko do testów"; exit 1; }
          check || exit 1
          "$DN5/make_dag_dn5.sh" || exit 1
          chmod +x "$DN5"/*.sh "$DN5"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -c '^JOB' "$DN5/dn5.dag"; grep -h '^PARENT' "$DN5/dn5.dag" | cut -c1-200; exit 0; fi
          if [ "${FORCE:-0}" = 1 ] && [ -e "$DN5_OUT/DONE" ]; then mv "$DN5_OUT/DONE" "$DN5_OUT/DONE.prev_$(date +%s)" || exit 1; fi
          mkdir -p "$DN5_OUT/gpu" "$Q/logs/dnet5" 2>/dev/null || true
          cd "$DN5" && "${CLEAN[@]}" condor_submit_dag -f -maxjobs "$DN5_MAXJOBS" -batch-name charm_dnet5 dn5.dag; rc=$?
          if [ $rc -eq 0 ] && [ "$DN5_GPU" = 1 ] && [ ! -s "$DN5_OUT/gpu/probe.json" ] && ! condor_q -batch 2>/dev/null | grep -q 'charm_dnet5gpu'; then
            date +%s > "$DN5_OUT/gpu/submitted"
            "${CLEAN[@]}" condor_submit -batch-name charm_dnet5gpu -append "penv = $(cat "$DN5/dn5.penv") DN5_GPU_BLOCKS=${DN5_GPU_BLOCKS// /+}" -append "wait = $DN5_GPU_WAIT" dn5_gpu.sub \
              || echo "UWAGA: sonda GPU nie wysłana (nie blokuje etapu)"
          fi
          exit $rc ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_dnet5|Total for klyczak'; ls "$DN5_OUT" 2>/dev/null | head -20
          echo "foldy $(ls "$DN5_OUT"/folds/*/fold*.json 2>/dev/null | wc -l) sc $(ls "$DN5_OUT/sc" 2>/dev/null | wc -l) meta $(ls "$DN5_OUT/meta" 2>/dev/null | wc -l) gpu $(ls "$DN5_OUT/gpu" 2>/dev/null | tr '\n' ' ') test $(ls "$DN5_OUT/test" 2>/dev/null | tr '\n' ' ') fail $(ls "$DN5_OUT/fail" 2>/dev/null | tr '\n' ' ')" ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
