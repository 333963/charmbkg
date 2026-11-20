#!/bin/bash
# submit_s427.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/salad427)
#  check  — pliki paczki, h5 v2.0 (DONE + pliki per próbka + kolumny), import charmbkg/sklearn w q_env, nic w kolejce (charm_salad427),
#           brak $S4_OUT/DONE, AFS < 85 %
#  submit — check + make_dag_s427.sh → condor_submit_dag s427.dag (batch charm_salad427) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
#  Uwaga: części w $S4_OUT/parts o tych samych nazwach są nadpisywane; przy zmianie S4_GROUP użyć nowego S4_OUT (study czyta wszystkie części).
set -uo pipefail
S4=$(cd "$(dirname "$0")" && pwd); export S4
source "$S4/s427_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n q f
  for f in s427_lib.py s427_extract.py s427_study.py s427_job.sh s427.sub src/charmbkg/bkggen/salad_v1.py src/charmbkg/data/salad_v1.yaml; do
    [ -f "$S4/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -f "$H5V20/DONE" ] || [ "${S4_NODONE:-0}" = 1 ] || { echo "STOP: brak $H5V20/DONE (konwersja h5_v20)"; ok=1; }
  for s in $S4_SAMPLES; do n=$(ls "$H5V20/$s/"*_"$S4_PRESEL".h5 2>/dev/null | grep -vc '\.part$'); echo "  h5_v20 $s: $n plików"; [ "$n" -gt 0 ] || { echo "STOP: brak h5_v20 $s"; ok=1; }; done
  f=$(ls "$H5V20/data_427929/"*_none.h5 2>/dev/null | head -1)
  [ -n "$f" ] && ( s4_env && cd "$S4" && python3 - "$f" <<'PY' ) || { echo "STOP: schemat h5_v20 / środowisko"; ok=1; }
import sys, h5py
import s427_lib as SL
names = h5py.File(sys.argv[1], "r")["events"].dtype.names
miss = [c for c in [SL.UID_COL, SL.GRL_COL, "Dmeson_hypo_swap"] + ["trigPassed_" + c for c in SL.CHAINS] if c not in names]
print("  schemat h5_v20:", "OK" if not miss else f"BRAK {miss}"); sys.exit(1 if miss else 0)
PY
  n=$(ls "$S4_OUT/parts" 2>/dev/null | wc -l)
  [ "$n" -gt 0 ] && echo "UWAGA: w $S4_OUT/parts jest już $n części (te o tych samych nazwach zostaną nadpisane)"
  condor_q -batch 2>/dev/null | grep -qE 'charm_salad427|s427\.dag' && { echo "STOP: SALAD427 w kolejce"; ok=1; }
  [ -s "$S4_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$S4_OUT/DONE")); FORCE=1 nadpisuje albo inny S4_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  S4_OUT=$S4_OUT LXY=$S4_LXY FRAC_J12=$S4_FRAC_J12 gatunki: $S4_SPECIES zadania: $S4_TASKS FAST=$S4_FAST"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_s427.sh submit)" ;;
  submit) [ -z "${S4_DAG:-}" ] || { echo "STOP: S4_DAG tylko do testów make_dag_s427.sh"; exit 1; }
          check || exit 1
          "$S4/make_dag_s427.sh" || exit 1
          chmod +x "$S4"/*.sh "$S4"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -h '^VARS' "$S4/s427.dag" | head -3; DRYRUN=1 "$S4/s427_job.sh" study D0; exit 0; fi
          # czyste env: LCG psuje condor_dagman (GLIBCXX_3.4.31)
          cd "$S4" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_salad427 s427.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_salad427|Total for klyczak'; ls -la "$S4_OUT" 2>/dev/null | head -20; ls "$S4_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$S4_OUT" -maxdepth 1 \( -name 's427_*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
