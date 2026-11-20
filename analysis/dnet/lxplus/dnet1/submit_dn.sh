#!/bin/bash
# submit_dn.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/dnet1)
#  check  — pliki paczki, h5 v2.0 (pliki per próbka + schemat), moduły (torch, xgboost, onnx) w q_env, nic w kolejce (charm_dnet1), brak $DN_OUT/DONE, AFS < 85 %
#  submit — check + make_dag_dn.sh → condor_submit_dag dn.dag (batch charm_dnet1) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
DN=$(cd "$(dirname "$0")" && pwd); export DN
source "$DN/dn_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n q f
  for f in dn_feat.py dn_extract.py dn_train.py dn_models.py analyze.py btrig_common.py dn_job.sh dn.sub src/charmbkg/features.py src/charmtag/dstar_refit.py ext/bethe_bloch_charmpy.py; do
    [ -f "$DN/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for s in $DN_SAMPLES; do n=$(ls "$H5V20/$s/"*_"$DN_PRESEL".h5 2>/dev/null | grep -vc '\.part$'); echo "  h5_v20 $s: $n plików"; [ "$n" -gt 0 ] || { echo "STOP: brak h5_v20 $s"; ok=1; }; done
  f=$(ls "$H5V20/data_428770/"*_none.h5 2>/dev/null | head -1)
  [ -n "$f" ] && ( dn_env >/dev/null && cd "$DN" && python3 - "$f" <<'PY' ) || { echo "STOP: schemat h5_v20 / środowisko"; ok=1; }
import sys, h5py
import dn_feat as DF
names = h5py.File(sys.argv[1], "r")["events"].dtype.names
miss = [c for c in DF.EXTRA + ["Dmeson_delta_mass", "Dmeson_truth_category"] if c not in names]
import torch, onnx, onnxruntime
print("  schemat h5_v20:", "OK" if not miss else f"BRAK {miss}", "| torch", torch.__version__, "onnx", onnx.__version__, "| konfiguracji", len(DF.all_configs()))
sys.exit(1 if miss else 0)
PY
  condor_q -batch 2>/dev/null | grep -qE 'charm_dnet1|dn\.dag' && { echo "STOP: DNET1 w kolejce"; ok=1; }
  [ -s "$DN_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$DN_OUT/DONE")); FORCE=1 nadpisuje albo inny DN_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  DN_OUT=$DN_OUT L_xy≥$DN_LXY frac D0/D+=$DN_FRAC_D0/$DN_FRAC_DP foldów $DN_FOLDS epok $DN_EPOCHS λ=$DN_LAM bloki: $DN_BLOCKS"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_dn.sh submit)" ;;
  submit) [ -z "${DN_DAG:-}" ] || { echo "STOP: DN_DAG tylko do testów make_dag_dn.sh"; exit 1; }
          check || exit 1
          "$DN/make_dag_dn.sh" || exit 1
          chmod +x "$DN"/*.sh "$DN"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -c '^JOB' "$DN/dn.dag"; grep -h '^VARS' "$DN/dn.dag" | head -2; exit 0; fi
          # czyste env: LCG psuje condor_dagman (GLIBCXX_3.4.31)
          cd "$DN" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_dnet1 dn.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_dnet1|Total for klyczak'; ls "$DN_OUT" 2>/dev/null | head -20; echo "parts $(ls "$DN_OUT/parts" 2>/dev/null | wc -l) meta $(ls "$DN_OUT/meta" 2>/dev/null | wc -l)" ;;
  fetch)  find "$DN_OUT" -maxdepth 2 \( -name '*.json' -o -name '*.onnx' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
