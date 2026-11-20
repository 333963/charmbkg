#!/bin/bash
# submit_fmm.sh check|submit|status|fetch   (lxplus, ~/private/charm_q/Q_lxplus/fmm)
#  check  — komplet h5 v2.0 (*_none.h5 = FMM_NFILES), model D0_H, charmtag, schemat (hypo_swap, trackId, trigPassed), nic w kolejce
#           (charm_fmm), brak $FMM_OUT/DONE (FORCE=1 nadpisuje), AFS < 85 %
#  submit — check + make_dag_fmm.sh → condor_submit_dag fmm.dag (batch charm_fmm) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
FM=$(cd "$(dirname "$0")" && pwd); export FM
source "$FM/fmm_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n e q f
  for f in fmm_extract.py fmm_merge.py fmm_job.sh fmm.sub; do [ -f "$FM/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -d "$FMM_SRC/charmtag" ] || { echo "STOP: brak $FMM_SRC/charmtag"; ok=1; }
  [ -s "$NOMINAL/models/D0_H/meta.json" ] || { echo "STOP: brak $NOMINAL/models/D0_H"; ok=1; }
  for s in $FMM_SAMPLES; do
    n=$(ls "$H5V20/$s" 2>/dev/null | grep -c "_${FMM_PRESEL}\.h5$"); e=$(fm_kv "$FMM_NFILES" "$s")
    echo "  $s: $n plików *_${FMM_PRESEL}.h5 (oczekiwane $e), $(fm_ngroups "$s") węzłów, frac $(fm_kv "$FMM_FRAC" "$s" 1.0)"
    [ "$n" -gt 0 ] && [ "$n" = "$e" ] || { echo "STOP: $s — niekompletne h5 v2.0 (zmiana liczby plików → FMM_NFILES)"; ok=1; }
  done
  f=$(ls "$H5V20/mb_nd/"*_none.h5 2>/dev/null | head -1)
  [ -n "$f" ] && ( fm_env >/dev/null && python3 - "$f" <<'PY' ) || { echo "STOP: schemat h5_v20"; ok=1; }
import sys, h5py
names = h5py.File(sys.argv[1], "r")["events"].dtype.names
need = ["Dmeson_hypo_swap", "Dmeson_truth_category", "Dmeson_daughter0_trackId", "Dmeson_daughter0_truthParentPdgId", "Dmeson_daughter0_pdgId",
        "event_uid", "passGRL_MinBiasAFP", "nSelBL", "leadTrackPtBL", "trigPassed_HLT_noalg_L1J12", "trigPassed_HLT_mb_mbts_L1MBTS_1"]
miss = [c for c in need if c not in names]
print("  schemat h5_v20:", "OK" if not miss else f"BRAK {miss}"); sys.exit(1 if miss else 0)
PY
  condor_q -batch 2>/dev/null | grep -qE 'charm_fmm|fmm\.dag' && { echo "STOP: FMM w kolejce"; ok=1; }
  [ -s "$FMM_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$FMM_OUT/DONE")); FORCE=1 nadpisuje"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  FMM_OUT=$FMM_OUT FMM_LXY=$FMM_LXY"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_fmm.sh submit)" ;;
  submit) [ -z "${FMM_DAG:-}" ] || { echo "STOP: FMM_DAG tylko do testów make_dag_fmm.sh"; exit 1; }
          check || exit 1
          "$FM/make_dag_fmm.sh" || exit 1
          chmod +x "$FM"/*.sh "$FM"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -h '^VARS' "$FM/fmm.dag" | head -3; DRYRUN=1 "$FM/fmm_job.sh" merge; exit 0; fi
          ( fm_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$FM" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_fmm fmm.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_fmm|Total for klyczak'; ls -la "$FMM_OUT" 2>/dev/null | head -20; ls "$FMM_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$FMM_OUT" -maxdepth 1 \( -name 'fmm_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
