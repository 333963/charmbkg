#!/bin/bash
# submit_bt.sh check|submit|status|fetch   (lxplus, ~/private/charm_q/Q_lxplus/b_trig). BT_PATH=A (tylko ROOT) | B | AB.
set -uo pipefail
BT=$(cd "$(dirname "$0")" && pwd); export BT
source "$BT/bt_common.sh"
: "${BT_PATH:=AB}"; export BT_PATH
cmd=${1:-check}
check() {
  local ok=0 s n q
  for f in btrig_common.py btrig_extract.py trg_root_extract.py bt_merge.py; do [ -f "$BT/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -d "$BT_SRC/charmtag" ] || { echo "STOP: brak $BT_SRC/charmtag"; ok=1; }
  [ -s "$NOMINAL/models/D0_H/meta.json" ] || { echo "STOP: brak $NOMINAL/models/D0_H"; ok=1; }
  if [[ $BT_PATH == *A* ]]; then
    for s in $BT_TRG; do n=$(ls "$RAWV20/$s/"*dpt3to150*.root 2>/dev/null | wc -l); echo "  ROOT $s: $n plików"; [ "$n" -gt 0 ] || { echo "STOP: brak ROOT $s"; ok=1; }; done
  fi
  if [[ $BT_PATH == *B* ]]; then
    [ -f "$H5V20/DONE" ] || [ "${BT_NODONE:-0}" = 1 ] || { echo "STOP: brak $H5V20/DONE (konwersja h5_v20 niezakończona)"; ok=1; }
    for s in $BT_SAMPLES; do n=$(ls "$H5V20/$s/"*_"$BT_PRESEL".h5 2>/dev/null | grep -vc '\.part$'); echo "  h5_v20 $s: $n plików"; [ "$n" -gt 0 ] || { echo "STOP: brak h5_v20 $s"; ok=1; }; done
    f=$(ls "$H5V20/data_428770/"*_none.h5 2>/dev/null | head -1)
    [ -n "$f" ] && ( bt_env >/dev/null && python3 - "$f" <<'PY' ) || { echo "STOP: schemat h5_v20"; ok=1; }
import sys, h5py
from btrig_common import find_trig_cols, find_col
names = h5py.File(sys.argv[1], "r")["events"].dtype.names
t = find_trig_cols(names); need = ["HLT_noalg_L1J12", "HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1"]
miss = [c for c in need if not (t.get(c) or {}).get("pass")] + [c for c in ("nSelBL", "leadTrackPtBL", "Dmeson_ptcone40", "Dmeson_Lxy", "passGRL_MinBiasAFP") if c not in names]
print("  schemat h5_v20:", "OK" if not miss else f"BRAK {miss}"); sys.exit(1 if miss else 0)
PY
  fi
  condor_q -batch 2>/dev/null | grep -qE 'charm_btrig' && { echo "STOP: B_trig w kolejce"; ok=1; }
  [ -s "$BT_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$BT_OUT/DONE")); FORCE=1 nadpisuje"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  BT_PATH=$BT_PATH BT_OUT=$BT_OUT BT_LXY=$BT_LXY FRAC data/mb=$BT_FRAC_DATA/$BT_FRAC_MB"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_bt.sh submit)" ;;
  submit) check || exit 1
          "$BT/make_dag_bt.sh" || exit 1
          chmod +x "$BT"/*.sh "$BT"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then grep -h '^VARS' "$BT/bt.dag" | head -3; exit 0; fi
          ( bt_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$BT" && env -u LD_LIBRARY_PATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_btrig bt.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_btrig|Total for klyczak'; ls -la "$BT_OUT" 2>/dev/null | head -20 ;;
  fetch)  find "$BT_OUT" -maxdepth 1 \( -name 'bt_*.npz' -o -name 'trg_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
