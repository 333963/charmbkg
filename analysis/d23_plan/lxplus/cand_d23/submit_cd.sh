#!/bin/bash
# submit_cd.sh check|submit|status|fetch   (lxplus, ~/private/charm_q/Q_lxplus/cand_d23; normalnie przez ORCH stages/candd23.sh)
#  check  — pliki paczki + src z pt_range, $H5D23/DONE (CD_ALLOW_NODONE=1 pomija), liczby plików h5 d23, modele foldów D⁺ na EOS,
#           nic w kolejce (charm_candd23), brak $CD_OUT/DONE (FORCE=1 nadpisuje), kwota AFS < 85 %, suchy przebieg ekstrakcji na
#           1 pliku (pierwszy plik pierwszej próbki danych i MC, 60 tys. wierszy, wynik w /tmp — kasowany)
#  submit — check + make_dag_cd.sh → condor_submit_dag -f cd.dag (batch charm_candd23) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
CD=$(cd "$(dirname "$0")" && pwd); export CD
source "$CD/cd_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n e q f k T done1=0 done2=0
  for f in cd_extract.py cd_merge.py cd_job.sh cd.sub make_dag_cd.sh src/charmbkg/samples.py; do [ -f "$CD/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  grep -q 'CHARMBKG_PT_RANGE' "$CD/src/charmbkg/samples.py" 2>/dev/null || { echo "STOP: src/charmbkg/samples.py bez pt_range (twarde pT ≥ 3!)"; ok=1; }
  if [ ! -s "$H5D23/DONE" ]; then
    if [ "${CD_ALLOW_NODONE:-0}" = 1 ]; then echo "UWAGA: brak $H5D23/DONE (CD_ALLOW_NODONE=1)"; else echo "STOP: brak $H5D23/DONE (etap h5_d23 nieskończony)"; ok=1; fi
  fi
  for s in $SAMPLES; do
    n=$(cd_files "$s" | wc -l); e=$(cd_kv "$CD_NFILES" "$s")
    echo "  h5 d23 $s: $n plików (oczekiwane $e)$([ "$n" = "$e" ] || echo '  ← UWAGA: inna liczba')"
    [ "$n" -gt 0 ] || { echo "STOP: brak $H5D23/$s/*_none.h5"; ok=1; }
  done
  for k in 0 1 2 3 4; do [ -s "$DP_MODELS/fold$k.json" ] || { echo "STOP: brak $DP_MODELS/fold$k.json"; ok=1; }; done
  [ -s "$DP_MODELS/meta.json" ] || echo "UWAGA: brak $DP_MODELS/meta.json (best_iteration z atrybutu modelu)"
  [ -s "$DP_THR" ] || echo "UWAGA: brak $DP_THR (summary bez progów odniesienia)"
  condor_q -batch 2>/dev/null | grep -qE 'charm_candd23|cd\.dag' && { echo "STOP: CAND_D23 w kolejce"; ok=1; }
  [ -s "$CD_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$CD_OUT/DONE")); FORCE=1 nadpisuje albo inny CD_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  if [ $ok = 0 ] && [ "${CD_NODRY:-0}" != 1 ]; then    # suchy przebieg na 1 pliku (dane + MC), w podpowłoce (LCG nie zostaje w env submitu)
    T=$(mktemp -d "${TMPDIR:-/tmp}/cd23_check_XXXXXX")
    for s in $SAMPLES; do
      case $s in data_*) [ $done1 = 1 ] && continue; done1=1 ;; *) [ $done2 = 1 ] && continue; done2=1 ;; esac
      f=$(cd_files "$s" | head -1)
      # shellcheck disable=SC2046
      ( cd_env >/dev/null && python3 -u "$CD/cd_extract.py" --sample "$s" --tag chk --out-dir "$T" $(cd_xargs) --max-rows 60000 "$f" 2>&1 | tail -2 \
        && python3 -c "import json,sys; j=json.load(open('$T/part_${s}_chk.json')); c=j['counts']; print('  suchy przebieg $s: D+', c['dp_sel'], 'zapisane', c['dp_kept'], 'pT', j['dp_pt_min'], j['dp_pt_max'], 'poza oknem', c['dp_pt_outside'], 'brak kolumn', j['missing_columns']); sys.exit(0 if c['dp_sel'] > 0 and j['dp_pt_min'] >= $CD_PT_LO and j['dp_pt_max'] < $CD_PT_HI else 1)" ) \
        || { echo "STOP: suchy przebieg $s ($f) nieudany (środowisko / schemat / modele / pT poza oknem)"; ok=1; }
    done
    case "$T" in */cd23_check_??????) rm -rf "${T:?}" ;; esac
  fi
  echo "kwota AFS: ${q:-?}%  CD_OUT=$CD_OUT H5D23=$H5D23 pT=[$CD_PT_LO,$CD_PT_HI) sc_full=$CD_SC_FULL ps_win=$CD_PS_WIN ps_out=$CD_PS_OUT max_gb=$CD_MAX_GB"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_cd.sh submit)" ;;
  submit) [ -z "${CD_DAG:-}" ] || { echo "STOP: CD_DAG tylko do testów make_dag_cd.sh"; exit 1; }
          check || exit 1
          "$CD/make_dag_cd.sh" || exit 1
          chmod +x "$CD"/*.sh "$CD"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -h '^VARS' "$CD/cd.dag" | head -3; DRYRUN=1 "$CD/cd_job.sh" finish; exit 0; fi
          # LCG view psuje condor_dagman (GLIBCXX_3.4.31) → submit w czystym env; -f: stare pliki DAG
          cd "$CD" && env -u LD_LIBRARY_PATH -u PYTHONPATH -u PYTHONHOME PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_candd23 cd.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_candd23|Total for klyczak'; ls -la "$CD_OUT" 2>/dev/null | head -20; echo "części: $(ls "$CD_OUT/parts" 2>/dev/null | grep -c '^part_')"
          [ -f "$CD_OUT/summary.json" ] && head -c 1500 "$CD_OUT/summary.json" ;;
  fetch)  find "$CD_OUT" -maxdepth 1 \( -name '*.parquet' -o -name '*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
