#!/bin/bash
# submit_tt.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/tagtrig)
#  check  — wejścia ROOT v2.0 na EOS (liczba *dpt3to150*.root = TT_NFILES per próbka); nic w kolejce (charm_tagtrig); brak $TT_OUT/DONE;
#           kwota AFS < 85 %; cvmfs Atheny (tylko UWAGA — węzeł lumi niekrytyczny); suchy przebieg: tt_bits.py na 1. pliku każdej próbki
#           (TT_CHECK_EVENTS zdarzeń, odczyt przez FUSE bez kopii) + tt_merge.py, w /tmp (kasowane), nic na EOS
#  submit — check + make_dag_tt.sh → condor_submit_dag -f tt.dag (batch charm_tagtrig) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy węzłów
#  status / fetch — kolejka + pliki $TT_OUT
# Zmienne: TT_OUT RAW20 TT_SAMPLES TT_NPARTS TT_NFILES TT_LUMI=0|1 TT_ASETUP TT_LUMITAG TT_GRL_XML TT_CHECK_EVENTS (0 = bez suchego przebiegu)
set -uo pipefail
TT=$(cd "$(dirname "$0")" && pwd); export TT
source "$TT/tt_common.sh"
cmd=${1:-check}
dry() {   # suchy przebieg: 1 plik × TT_CHECK_EVENTS zdarzeń per próbka → części → merge; wszystko w /tmp
  local tmp s rc=0
  tmp=$(mktemp -d /tmp/tt_check_XXXXXX) || return 1
  (
    trap 'rm -rf "$tmp"' EXIT
    tt_env >/dev/null 2>&1 || { echo "STOP: q_env bez numpy/awkward/uproot/pandas/pyarrow"; exit 1; }
    mkdir -p "$tmp/parts"
    for s in $TT_SAMPLES; do
      python3 "$TT/tt_bits.py" --sample "$s" --rootdir "$RAW20/$s" --ptprod "$PTPROD" --part 0 --nparts 1 --max-files 1 \
          --max-events "$TT_CHECK_EVENTS" --no-copy --step "$TT_STEP" --out-prefix "$tmp/parts/$s.p0" 2>&1 | tail -1 || exit 1
      [ -s "$tmp/parts/$s.p0.meta.json" ] || exit 1
    done
    python3 "$TT/tt_merge.py" --out "$tmp" --wd "$tmp/fin" --samples "$TT_SAMPLES" --nparts 1 --allow-incomplete --max-mb "$TT_MAX_MB" 2>&1 | tail -4 || exit 1
    [ -s "$tmp/fin/summary.json" ] || exit 1
    python3 - "$tmp/fin/summary.json" <<'PY' || exit 1
import json, sys
S = json.load(open(sys.argv[1]))
for s, r in S["samples"].items():
    assert r["n_chains"] >= 50 and not r["trig11_missing"], f"{s}: łańcuchów {r['n_chains']}, brak z 11: {r['trig11_missing']}"
    assert r["n_T0"] > 0, f"{s}: 0 zdarzeń tagowanych w suchym przebiegu"
    print(f"  suchy przebieg {s}: łańcuchów {r['n_chains']} + L1 {r['n_l1tbp']}, bez trigPS {len(r['ps_missing'])}, T0 {r['n_T0']} z {r['n_grl']} GRL")
PY
  ) || rc=1
  return $rc
}
check() {
  local ok=0 s q n e
  [ -f "$TT/tt_bits.py" ] && [ -f "$TT/tt_merge.py" ] && [ -f "$TT/tt_lumi.sh" ] || { echo "STOP: brak tt_bits.py/tt_merge.py/tt_lumi.sh"; ok=1; }
  for s in $TT_SAMPLES; do
    n=$(tt_nroot "$s"); e=$(tt_nexp "$s")
    echo "  $s: $n plików *${PTPROD}*.root (oczekiwane $e) w $RAW20/$s"
    [ "$n" -gt 0 ] && [ "$n" = "$e" ] || { echo "STOP: $s — brak albo inna liczba plików ROOT v2.0 (zmiana → TT_NFILES)"; ok=1; }
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_tagtrig|tt\.dag' && { echo "STOP: TAGTRIG w kolejce"; ok=1; }
  [ -s "$TT_OUT/DONE" ] && { echo "STOP: TAGTRIG już policzone ($(cat "$TT_OUT/DONE")) — inny TT_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}'); q=${q%.*}
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  if [ "$TT_LUMI" = 1 ]; then
    [ -f /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh ] || echo "UWAGA: brak cvmfs ATLASLocalRootBase — węzeł lumi zapisze tylko status"
    [ -n "$TT_LUMITAG" ] || echo "UWAGA: TT_LUMITAG pusty → bez iLumiCalc (świetlności per LB NIE będzie); tylko klucze/preskale z COOL"
  fi
  if [ $ok = 0 ] && [ "${TT_CHECK_EVENTS:-0}" -gt 0 ]; then dry || { echo "STOP: suchy przebieg (tt_bits.py/tt_merge.py na 1 pliku) nie przeszedł"; ok=1; }; fi
  echo "kwota AFS: ${q:-?}%  TT_OUT=$TT_OUT RAW20=$RAW20 TT_SAMPLES=$TT_SAMPLES TT_NPARTS=$TT_NPARTS TT_LUMI=$TT_LUMI"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_tt.sh submit)" ;;
  submit) [ -z "${TT_DAG:-}" ] || { echo "STOP: TT_DAG tylko do testów make_dag_tt.sh"; exit 1; }
          if [ "${DRYRUN:-0}" = 1 ]; then
            "$TT/make_dag_tt.sh" || exit 1
            for s in $TT_SAMPLES; do DRYRUN=1 "$TT/tt_job.sh" "$s" 0; done; [ "$TT_LUMI" = 1 ] && DRYRUN=1 "$TT/tt_job.sh" lumi; DRYRUN=1 "$TT/tt_job.sh" finish; exit 0
          fi
          check || exit 1
          "$TT/make_dag_tt.sh" || exit 1
          chmod +x "$TT"/*.sh "$TT"/*.py
          # czyste env: LCG przed condor_submit_dag psuje condor_dagman (GLIBCXX_3.4.31); q_env było tylko w podpowłoce (dry)
          cd "$TT" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_tagtrig tt.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_tagtrig|Total for klyczak'; ls -la "$TT_OUT" 2>/dev/null | head -30; echo "części: $(ls "$TT_OUT/parts" 2>/dev/null | grep -c 'meta\.json$')"
          [ -f "$TT_OUT/lumi_STATUS.txt" ] && tail -n 12 "$TT_OUT/lumi_STATUS.txt" ;;
  fetch)  find "$TT_OUT" -maxdepth 1 \( -name '*.parquet' -o -name '*.csv' -o -name '*.json' -o -name '*.txt' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
