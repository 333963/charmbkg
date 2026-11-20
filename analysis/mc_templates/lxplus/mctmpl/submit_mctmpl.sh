#!/bin/bash
# submit_mctmpl.sh check|submit|status|fetch   (lxplus, ~/private/charm_q/Q_lxplus/mctmpl)
#  check  — komplet h5 v2.0 (*_none.h5 = MT_NFILES), modele D0_H i Dplus_H, charmtag, progi WP, schemat, nic w kolejce
#           (charm_mctmpl), brak $MT_OUT/DONE (FORCE=1 nadpisuje), AFS < 85 %
#  submit — check + make_dag_mctmpl.sh → condor_submit_dag mctmpl.dag (batch charm_mctmpl) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
MT=$(cd "$(dirname "$0")" && pwd); export MT
source "$MT/mctmpl_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n e q f
  for f in mctmpl_hist.py mctmpl_merge.py mctmpl_analyze.py mctmpl_job.sh mctmpl.sub; do [ -f "$MT/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -s "$MT_WP" ] || { echo "STOP: brak progów $MT_WP"; ok=1; }
  [ -d "$MT_SRC/charmtag" ] || { echo "STOP: brak $MT_SRC/charmtag"; ok=1; }
  for s in D0_H Dplus_H; do [ -s "$NOMINAL/models/$s/meta.json" ] || { echo "STOP: brak $NOMINAL/models/$s"; ok=1; }; done
  for s in $MT_SAMPLES; do
    n=$(ls "$H5V20/$s" 2>/dev/null | grep -c "_${MT_PRESEL}\.h5$"); e=$(mt_kv "$MT_NFILES" "$s")
    echo "  $s: $n plików *_${MT_PRESEL}.h5 (oczekiwane $e), $(mt_ngroups "$s") węzłów"
    [ "$n" -gt 0 ] && [ "$n" = "$e" ] || { echo "STOP: $s — niekompletne h5 v2.0 (zmiana liczby plików → MT_NFILES)"; ok=1; }
  done
  f=$(ls "$H5V20/mb_nd/"*_none.h5 2>/dev/null | head -1)
  [ -n "$f" ] && ( mt_env >/dev/null && python3 - "$f" "$MT_WP" <<'PY' ) || { echo "STOP: schemat h5_v20 / progi"; ok=1; }
import sys, json, h5py
names = h5py.File(sys.argv[1], "r")["events"].dtype.names
need = ["Dmeson_hypo_swap", "Dmeson_truth_category", "Dmeson_truth_channel", "Dmeson_truthPdgId", "Dmeson_daughter0_trackId", "Dmeson_delta_mass",
        "Dmeson_mKpi1", "Dmeson_impact", "event_uid", "passGRL_MinBiasAFP", "nSelBL", "leadTrackPtBL", "nPileupVertices", "afp_tag_A", "afp_xi1_A",
        "xi_cal_A", "trigPassed_HLT_noalg_L1J12", "trigPassed_HLT_mb_mbts_L1MBTS_1", "trigPassed_HLT_noalg_L1AFP_A_OR_C"]
miss = [c for c in need if c not in names]
w = json.load(open(sys.argv[2]))
bad = [f"{s}:{p}:{k}" for s in ("D0", "Dplus") for p in ("3-5", "5-8", "8-30") for k in ("0.95", "0.80", "0.50") if not 0 < w[s][p][k] < 1]
print("  schemat h5_v20:", "OK" if not miss else f"BRAK {miss}", "| progi WP:", "OK" if not bad else f"ZŁE {bad}"); sys.exit(1 if miss or bad else 0)
PY
  condor_q -batch 2>/dev/null | grep -qE 'charm_mctmpl|mctmpl\.dag' && { echo "STOP: MCTMPL w kolejce"; ok=1; }
  [ -s "$MT_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$MT_OUT/DONE")); FORCE=1 nadpisuje"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  MT_OUT=$MT_OUT"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_mctmpl.sh submit)" ;;
  submit) [ -z "${MT_DAG:-}" ] || { echo "STOP: MT_DAG tylko do testów make_dag_mctmpl.sh"; exit 1; }
          check || exit 1
          "$MT/make_dag_mctmpl.sh" || exit 1
          chmod +x "$MT"/*.sh "$MT"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -h '^VARS' "$MT/mctmpl.dag" | head -3; DRYRUN=1 "$MT/mctmpl_job.sh" merge; exit 0; fi
          ( mt_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$MT" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_mctmpl mctmpl.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_mctmpl|Total for klyczak'; ls -la "$MT_OUT" 2>/dev/null | head -20; ls "$MT_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$MT_OUT" -maxdepth 1 \( -name 'mctmpl_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
