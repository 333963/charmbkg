#!/bin/bash
# submit_cn.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/cnclone)
#  check  — pliki paczki, h5 v2.0 (DONE + pliki per próbka + schemat kolumn), torch w q_env, nic w kolejce (charm_cnclone), brak $CN_OUT/DONE, AFS < 85 %
#  submit — check + make_dag_cn.sh → condor_submit_dag cn.dag (batch charm_cnclone) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
#  Uwaga: stare części w $CN_OUT/parts z innym CN_GROUP nie są kasowane — przy zmianie CN_GROUP użyć nowego CN_OUT.
set -uo pipefail
CN=$(cd "$(dirname "$0")" && pwd); export CN
source "$CN/cn_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n q f
  for f in cn_feat.py cn_extract.py cn_train.py cn_job.sh cn.sub; do [ -f "$CN/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -f "$H5V20/DONE" ] || [ "${CN_NODONE:-0}" = 1 ] || { echo "STOP: brak $H5V20/DONE (konwersja h5_v20)"; ok=1; }
  for s in $CN_SAMPLES; do n=$(ls "$H5V20/$s/"*_"$CN_PRESEL".h5 2>/dev/null | grep -vc '\.part$'); echo "  h5_v20 $s: $n plików"; [ "$n" -gt 0 ] || { echo "STOP: brak h5_v20 $s"; ok=1; }; done
  f=$(ls "$H5V20/data_428770/"*_none.h5 2>/dev/null | head -1)
  [ -n "$f" ] && ( cn_env >/dev/null && cd "$CN" && python3 - "$f" <<'PY' ) || { echo "STOP: schemat h5_v20 / środowisko"; ok=1; }
import sys, h5py
import cn_feat as CF
names = h5py.File(sys.argv[1], "r")["events"].dtype.names
miss = [c for c in CF.RAW_COLS + list(CF.EVENT_COLS.values()) + ["trigPassed_" + c for c in CF.TRIG_CHAINS] if c not in names]
try:
    import torch; be = "torch " + torch.__version__
except Exception:
    be = "BRAK torch → zamiennik HGB"
print("  schemat h5_v20:", "OK" if not miss else f"BRAK {miss}", "| backend:", be); sys.exit(1 if miss else 0)
PY
  n=$(ls "$CN_OUT/parts" 2>/dev/null | wc -l)
  [ "$n" -gt 0 ] && [ "${FORCE:-0}" != 1 ] && echo "UWAGA: w $CN_OUT/parts jest już $n części (zostaną nadpisane te o tych samych nazwach)"
  condor_q -batch 2>/dev/null | grep -qE 'charm_cnclone|cn\.dag' && { echo "STOP: CNCLONE w kolejce"; ok=1; }
  [ -s "$CN_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$CN_OUT/DONE")); FORCE=1 nadpisuje albo inny CN_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  CN_OUT=$CN_OUT CN_LXY=$CN_LXY FRAC data/mb=$CN_FRAC_DATA/$CN_FRAC_MB foldów $CN_FOLDS epok $CN_EPOCHS"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_cn.sh submit)" ;;
  submit) [ -z "${CN_DAG:-}" ] || { echo "STOP: CN_DAG tylko do testów make_dag_cn.sh"; exit 1; }
          check || exit 1
          "$CN/make_dag_cn.sh" || exit 1
          chmod +x "$CN"/*.sh "$CN"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -h '^VARS' "$CN/cn.dag" | head -3; DRYRUN=1 "$CN/cn_job.sh" train; exit 0; fi
          # czyste env: LCG psuje condor_dagman (GLIBCXX_3.4.31)
          cd "$CN" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_cnclone cn.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_cnclone|Total for klyczak'; ls -la "$CN_OUT" 2>/dev/null | head -20; ls "$CN_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$CN_OUT" -maxdepth 1 \( -name 'cn_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
