#!/bin/bash
# submit_dlz.sh check|submit|status|fetch   (lxplus, ~/private/charm_q/Q_lxplus/dlz_data)
#  check  — pliki paczki, komplet h5 v2.0 (= DLZ_NFILES), h5 v1.2 dla D⁰, środowisko (charmbkg, xgboost, iminuit), schemat, nic w kolejce
#           (charm_dlz), brak $DLZ_OUT/DONE (FORCE=1 nadpisuje), AFS < 85 %
#  submit — check + make_dag_dlz.sh → condor_submit_dag dlz.dag (batch charm_dlz) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
DZ=$(cd "$(dirname "$0")" && pwd); export DZ
source "$DZ/dlz_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n e q f
  for f in dlzlib.py apply.py dlz_frames.py dlz_train.py dlz_merge.py analyze.py fitlib.py mc_checks.py dlz_job.sh dlz.sub d0/d0_run.py d0/study.py; do
    [ -f "$DZ/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -d "$DZ/src/charmbkg" ] || { echo "STOP: brak src/charmbkg"; ok=1; }
  for s in data_428770 data_427929 mb_nd fd_nd; do
    n=$(ls "$H5V20/$s" 2>/dev/null | grep -c "_${DLZ_PRESEL}\.h5$"); e=$(dz_kv "$DLZ_NFILES" "$s")
    echo "  h5_v20 $s: $n plików (oczekiwane $e)"
    [ "$n" -gt 0 ] && [ "$n" = "$e" ] || { echo "STOP: $s — niekompletne h5 v2.0 (zmiana liczby plików → DLZ_NFILES)"; ok=1; }
  done
  if [ "$DLZ_D0" = 1 ]; then for s in fd_nd mb_nd data_428770; do n=$(ls "$H5V12/$s" 2>/dev/null | grep -c "_none\.h5$"); echo "  h5_v12 $s: $n plików"; [ "$n" -gt 0 ] || { echo "STOP: brak h5 v1.2 $s"; ok=1; }; done; fi
  f=$(ls "$H5V20/data_427929/"*_none.h5 2>/dev/null | head -1)
  [ -n "$f" ] && ( dz_env >/dev/null && python3 - "$f" <<'PY' ) || { echo "STOP: schemat h5_v20 / środowisko"; ok=1; }
import sys, h5py
names = h5py.File(sys.argv[1], "r")["events"].dtype.names
need = ["Dmeson_mKpi1", "Dmeson_mKpi2", "Dmeson_maxdR", "Dmeson_impact", "Dmeson_Lxy", "passGRL_MinBiasAFP", "trigPassed_HLT_mb_mbts_L1MBTS_1", "trigPassed_HLT_noalg_L1J12", "nSelBL"]
miss = [c for c in need if c not in names]
print("  schemat h5_v20:", "OK" if not miss else f"BRAK {miss}"); sys.exit(1 if miss else 0)
PY
  condor_q -batch 2>/dev/null | grep -qE 'charm_dlz|dlz\.dag' && { echo "STOP: DLZ_DATA w kolejce"; ok=1; }
  [ -s "$DLZ_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$DLZ_OUT/DONE")); FORCE=1 nadpisuje albo inny DLZ_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  DLZ_OUT=$DLZ_OUT cap=$DLZ_CAP sbfrac=$DLZ_SBFRAC sets=[$DLZ_SETS] d0=$DLZ_D0 fast=$DLZ_FAST"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_dlz.sh submit)" ;;
  submit) [ -z "${DLZ_DAG:-}" ] || { echo "STOP: DLZ_DAG tylko do testów make_dag_dlz.sh"; exit 1; }
          check || exit 1
          "$DZ/make_dag_dlz.sh" || exit 1
          chmod +x "$DZ"/*.sh "$DZ"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -h '^VARS' "$DZ/dlz.dag" | head -3; DRYRUN=1 "$DZ/dlz_job.sh" merge; exit 0; fi
          cd "$DZ" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_dlz dlz.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_dlz|Total for klyczak'; ls -la "$DLZ_OUT" 2>/dev/null | head -20; ls "$DLZ_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$DLZ_OUT" -maxdepth 2 \( -name 'dlz_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
