#!/bin/bash
# submit_pl.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/peak_lxy010)
#  check  — wejścia h5 v2.0 na EOS (pliki per próbka + DONE etapu h5_v20), modele D0_H (fold*.json, meta.json; brak = STOP, bo
#           score jest liczony z modeli), src z hyp_extract/dstar_tag, nic w kolejce (charm_peaklxy), brak $PL_OUT/DONE, kwota AFS < 85 %
#  submit — check + make_dag_pl.sh → condor_submit_dag pl.dag (batch charm_peaklxy); DRYRUN=1 → tylko DAG + komendy węzłów
#  status / fetch — kolejka + pliki $PL_OUT
# Zmienne: SAMPLES="fd_nd mb_nd data_428770 data_427929" PL_OUT= H5V20= NOMINAL= PL_EXTRA=--dstar-only DATA_FRAC=1 SCORE_MODE=model
set -uo pipefail
PL=$(cd "$(dirname "$0")" && pwd); export PL
source "$PL/pl_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n
  [ -f "$PL/src/charmtag/hyp_extract.py" ] && [ -f "$PL/src/charmtag/dstar_tag.py" ] || { echo "STOP: $PL/src bez charmtag/hyp_extract.py lub dstar_tag.py"; ok=1; }
  [ -s "$H5V20/DONE" ] || echo "UWAGA: brak $H5V20/DONE (etap h5_v20 niekompletny?) — sprawdź liczby plików niżej (oczekiwane: fd_nd 24, mb_nd 78, data_428770 97, data_427929 458)"
  for s in $SAMPLES; do
    n=$(ls "$H5V20/$s/"*_none.h5 2>/dev/null | grep -vc '\.part$')
    [ "$n" -gt 0 ] || { echo "STOP: brak $H5V20/$s/*_none.h5"; ok=1; }
    echo "  $s: $n plików h5 v2.0"
  done
  if [ "$SCORE_MODE" = model ]; then
    [ -s "$NOMINAL/models/D0_H/meta.json" ] && [ -s "$NOMINAL/models/D0_H/fold0.json" ] || { echo "STOP: brak $NOMINAL/models/D0_H/{meta.json,fold*.json} (SCORE_MODE=model)"; ok=1; }
  fi
  condor_q -batch 2>/dev/null | grep -qE 'charm_peaklxy|pl\.dag' && { echo "STOP: PEAK_LXY010 w kolejce"; ok=1; }
  [ -s "$PL_OUT/DONE" ] && { echo "STOP: już policzone ($(cat "$PL_OUT/DONE")) — inny PL_OUT albo usuń DONE"; ok=1; }
  local q; q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  PL_OUT=$PL_OUT H5V20=$H5V20 SCORE_MODE=$SCORE_MODE PL_EXTRA=${PL_EXTRA:-} DATA_FRAC=$DATA_FRAC"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_pl.sh submit)" ;;
  submit) check || exit 1
          "$PL/make_dag_pl.sh" || exit 1
          chmod +x "$PL"/*.sh
          if [ "${DRYRUN:-0}" = 1 ]; then for s in $SAMPLES finish; do DRYRUN=1 "$PL/pl_job.sh" "$s"; done; exit 0; fi
          ( pl_env ) || exit 1   # podpowłoka: LCG view psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$PL" && env -u LD_LIBRARY_PATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_peaklxy pl.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_peaklxy|Total for klyczak'; ls -la "$PL_OUT" 2>/dev/null | head -20
          [ -f "$PL_OUT/extract_summary.json" ] && head -c 2500 "$PL_OUT/extract_summary.json" ;;
  fetch)  find "$PL_OUT" -maxdepth 1 \( -name '*.parquet' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
