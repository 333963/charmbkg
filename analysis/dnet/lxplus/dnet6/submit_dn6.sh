#!/bin/bash
# submit_dn6.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/dnet6)
#  check  — pliki paczki i ../dnet{5,4,2,1}; wejścia na EOS (tabele etapu 1, tory otoczenia D0 etapu 2 — komplet, score drzewa H etapu 1, score pełnej sieci
#           etapu 2 ziarna 1–3; DN6_APP=1: h5/ROOT v2.0 i drzewo odniesienia etapu 4); nic w kolejce (charm_dnet6); brak $DN6_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_dn6.sh → condor_submit_dag -f -maxjobs $DN6_MAXJOBS (batch charm_dnet6) w CZYSTYM env. Tryb: DN6_TRAIN=gpu (domyślnie) | cpu.
#           Tagi gotowe na EOS są pomijane → ponowny submit w trybie cpu po nieudanym GPU liczy tylko brakujące ziarna. DRYRUN=1 → tylko DAG.
set -uo pipefail
DN6=$(cd "$(dirname "$0")" && pwd); export DN6
source "$DN6/dn6_common.sh"
cmd=${1:-check}
CLEAN=(env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin)
check() {
  local ok=0 s n q f b x kv
  for f in dn6_train.py dn6_mask.py dn6_apply.py dn6_ana.py dn6_study.py dn6_job.sh dn6.sub dn6_gtrain.sub make_dag_dn6.sh dn6_switch.sh; do [ -f "$DN6/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for f in dn5_train.py dn5_lib.py dn5_models.py dn5_common.sh; do [ -f "$DN5/$f" ] || { echo "STOP: brak ../dnet5/$f"; ok=1; }; done
  for f in dn4_train.py dn4_lib.py dn4_apply.py dn4_extract.py dn4_ana.py dn4_common.sh; do [ -f "$DN4/$f" ] || { echo "STOP: brak ../dnet4/$f"; ok=1; }; done
  for f in dn2_train.py dn2_models.py dn2_extract.py dn2_common.sh; do [ -f "$DN2/$f" ] || { echo "STOP: brak ../dnet2/$f"; ok=1; }; done
  for f in dn_feat.py dn_train.py dn_models.py analyze.py btrig_common.py dn_common.sh; do [ -f "$DN1/$f" ] || { echo "STOP: brak ../dnet1/$f"; ok=1; }; done
  for b in $DN6_BLOCKS; do [ -n "$(dn6_tag "$b")" ] || { echo "STOP: nieznany blok $b"; ok=1; }; done
  case $DN6_TRAIN in gpu|cpu) ;; *) echo "STOP: DN6_TRAIN=$DN6_TRAIN (gpu|cpu)"; ok=1 ;; esac
  for s in $DN_SAMPLES; do
    [ -s "$DN_OUT/data/dn_$s.npz" ] || { echo "STOP: brak $DN_OUT/data/dn_$s.npz (etap dnet1)"; ok=1; }
    x=$(dn2_kv "$DN2_NPARTS" "$s" 4)
    n=$(ls "$DN2_OUT/parts/dn2_${s}_"*.npz 2>/dev/null | wc -l); echo "  tory otoczenia D0 $s: $n części (oczekiwane $x)"
    [ "$n" -eq "$x" ] || { echo "STOP: niepełne części etapu 2 dla $s"; ok=1; }
  done
  [ -s "$DN_OUT/sc/D0.xgb.H.sb.n.s1.npz" ] || { echo "STOP: brak score drzewa H etapu 1"; ok=1; }
  for s in 1 2 3; do [ -s "$DN2_OUT/sc/D0.tr.H+e04.mc.d${DN_LAM%%,*}.s$s.npz" ] || { echo "STOP: brak score pełnej sieci etapu 2 (ziarno $s)"; ok=1; }; done
  if [ "$DN6_APP" = 1 ]; then
    [ -s "$DN4_OUT/models/$DN6_APP_REF/meta.json" ] || { echo "STOP: brak drzewa odniesienia $DN4_OUT/models/$DN6_APP_REF (DN6_APP=0 wyłącza app)"; ok=1; }
    for kv in $DN6_APP150; do
      s=${kv%%:*}; n=$(dn4_files d150 "$s" | wc -l); echo "  app $s: h5 $n, węzłów ${kv##*:}"
      [ "$n" -gt 0 ] || { echo "STOP: brak h5 v2.0 dla $s"; ok=1; }
    done
  fi
  n=0; for b in $DN6_BLOCKS; do if dn6_ready "$(dn6_tag "$b")"; then n=$((n + 1)); fi; done
  condor_q -batch 2>/dev/null | grep -qE 'charm_dnet6( |$)|dn6\.dag' && { echo "STOP: DNET6 w kolejce"; ok=1; }
  [ -s "$DN6_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$DN6_OUT/DONE")); FORCE=1 albo inny DN6_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  DN6_OUT=$DN6_OUT tryb $DN6_TRAIN (GPU czeka ≤ ${DN6_GPU_WAIT}s) bloki: $DN6_BLOCKS (gotowych tagów: $n) maska dm<$DN6_DM_MAX znak $DN6_SIGN app: $DN6_APP"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_dn6.sh submit)" ;;
  submit) [ -z "${DN6_DAG:-}" ] || { echo "STOP: DN6_DAG tylko do testów"; exit 1; }
          check || exit 1
          "$DN6/make_dag_dn6.sh" || exit 1
          chmod +x "$DN6"/*.sh "$DN6"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -c '^JOB' "$DN6/dn6.dag"; grep -h '^PARENT' "$DN6/dn6.dag" | cut -c1-200; exit 0; fi
          if [ "${FORCE:-0}" = 1 ] && [ -e "$DN6_OUT/DONE" ]; then mv "$DN6_OUT/DONE" "$DN6_OUT/DONE.prev_$(date +%s)" || exit 1; fi
          mkdir -p "$DN6_OUT" "$Q/logs/dnet6" 2>/dev/null || true
          cd "$DN6" && "${CLEAN[@]}" condor_submit_dag -f -maxjobs "$DN6_MAXJOBS" -batch-name charm_dnet6 dn6.dag; rc=$?
          [ $rc -eq 0 ] && echo "$DN6_TRAIN $(date +%s)" >> "$DN6/dn6.submits"
          exit $rc ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_dnet6|Total for klyczak'
          condor_q -nobatch -constraint 'JobBatchName == "charm_dnet6" && JobUniverse != 7' -af:j JobStatus RequestGPUs QDate Args 2>/dev/null | head -20
          ls "$DN6_OUT" 2>/dev/null | head -20
          echo "foldy $(ls "$DN6_OUT"/folds/*/fold*.json 2>/dev/null | wc -l) sc $(ls "$DN6_OUT/sc" 2>/dev/null | wc -l) meta $(ls "$DN6_OUT/meta" 2>/dev/null | wc -l) app $(ls "$DN6_OUT/app/d150" 2>/dev/null | wc -l) fail $(ls "$DN6_OUT/fail" 2>/dev/null | tr '\n' ' ')" ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
