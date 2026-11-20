#!/bin/bash
# submit_afp.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/afp_ext)
#  check  — wejścia: h5 v2.0 *_none.h5 w komplecie (liczba = AFP_NFILES) per próbka; nic w kolejce (charm_afp); brak $AFP_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_afp.sh → condor_submit_dag afp.dag (batch charm_afp) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
AX=$(cd "$(dirname "$0")" && pwd); export AX
source "$AX/afp_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s q n e
  [ -f "$AX/afp_extract.py" ] && [ -f "$AX/afp_merge.py" ] || { echo "STOP: brak afp_extract.py/afp_merge.py"; ok=1; }
  for s in $AFP_SAMPLES; do
    n=$(ls "$H5V20/$s" 2>/dev/null | grep -c '_none\.h5$'); e=$(ax_nexp "$s")
    echo "  $s: $n plików *_none.h5 (oczekiwane $e) w $H5V20/$s"
    [ "$n" -gt 0 ] && [ "$n" = "$e" ] || { echo "STOP: $s — niekompletne h5 v2.0 (etap h5_v20 trwa? zmiana liczby plików → AFP_NFILES)"; ok=1; }
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_afp|afp\.dag' && { echo "STOP: AFP_EXT w kolejce"; ok=1; }
  condor_q -batch 2>/dev/null | grep -qE 'charm_h5v20' && echo "UWAGA: charm_h5v20 jeszcze w kolejce (dane mogą być kompletne — decyduje liczba plików wyżej)"
  [ -s "$AFP_OUT/DONE" ] && { echo "STOP: AFP_EXT już policzone ($(cat "$AFP_OUT/DONE")) — inny AFP_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  AFP_OUT=$AFP_OUT H5V20=$H5V20 AFP_SAMPLES=$AFP_SAMPLES AFP_NPARTS=$AFP_NPARTS"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_afp.sh submit)" ;;
  submit) [ -z "${AFP_DAG:-}" ] || { echo "STOP: AFP_DAG tylko do testów make_dag_afp.sh"; exit 1; }
          check || exit 1
          "$AX/make_dag_afp.sh" || exit 1
          chmod +x "$AX"/*.sh "$AX"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then for s in $AFP_SAMPLES; do DRYRUN=1 "$AX/afp_job.sh" "$s" 0; done; DRYRUN=1 "$AX/afp_job.sh" finish; exit 0; fi
          ( ax_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$AX" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_afp afp.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_afp|Total for klyczak'; ls -la "$AFP_OUT" "$AFP_OUT/parts" 2>/dev/null | head -50 ;;
  fetch)  find "$AFP_OUT" -maxdepth 1 \( -name '*.parquet' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
