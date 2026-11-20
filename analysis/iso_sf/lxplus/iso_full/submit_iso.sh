#!/bin/bash
# submit_iso.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/iso_full)
#  check  — wejścia: T4 $TAG_DIR/<próbka>.parquet+meta, h5 v1.2 z meta obecne; nic w kolejce (charm_iso); brak $ISO_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_iso.sh → condor_submit_dag iso.dag (batch charm_iso) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
IX=$(cd "$(dirname "$0")" && pwd); export IX
source "$IX/iso_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s q
  [ -f "$IX/iso_extract.py" ] && [ -f "$IX/iso_merge.py" ] || { echo "STOP: brak iso_extract.py/iso_merge.py"; ok=1; }
  for s in $ISO_SAMPLES; do
    [ -s "$TAG_DIR/$s.parquet" ] && [ -s "$TAG_DIR/$s.meta.json" ] || { echo "STOP: brak $TAG_DIR/$s.parquet|.meta.json"; ok=1; continue; }
    python3 - "$TAG_DIR/$s.meta.json" "$H5V12/$s" <<'PY' || ok=1
import json, os, sys
m = json.load(open(sys.argv[1])); miss = [f for f in m["files"] if not os.path.exists(os.path.join(sys.argv[2], f))]
print(f"  {m['sample']}: {len(m['files'])} plików w meta T4, brak w h5 v1.2: {len(miss)} {miss[:3]}, wierszy T4 {m['n_rows_out']}")
sys.exit(1 if miss else 0)
PY
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_iso|iso.dag' && { echo "STOP: ISO w kolejce"; ok=1; }
  [ -s "$ISO_OUT/DONE" ] && { echo "STOP: ISO już policzone ($(cat "$ISO_OUT/DONE")) — inny ISO_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  ISO_OUT=$ISO_OUT TAG_DIR=$TAG_DIR ISO_SAMPLES=$ISO_SAMPLES ISO_NPARTS_DATA=$ISO_NPARTS_DATA"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_iso.sh submit)" ;;
  submit) [ -z "${ISO_DAG:-}" ] || { echo "STOP: ISO_DAG tylko do testów make_dag_iso.sh"; exit 1; }
          check || exit 1
          "$IX/make_dag_iso.sh" || exit 1
          chmod +x "$IX"/*.sh "$IX"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then for s in $ISO_SAMPLES; do DRYRUN=1 "$IX/iso_job.sh" "$s" 0; done; DRYRUN=1 "$IX/iso_job.sh" finish; exit 0; fi
          ( ix_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$IX" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_iso iso.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_iso|Total for klyczak'; ls -la "$ISO_OUT" "$ISO_OUT/parts" 2>/dev/null | head -30 ;;
  fetch)  find "$ISO_OUT" -maxdepth 1 \( -name '*.parquet' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
