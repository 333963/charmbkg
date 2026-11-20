#!/bin/bash
# Węzeł DAG DNET etap 4. Użycie:
#   dn4_job.sh smoke                       dn4_train.py selftest --need-torch (dane syntetyczne: pule, tokeny, drzewo, sieci 1 epoka, zapis/odczyt
#                                          modeli, dn4_apply) + dn4_ana.py selftest (zabawki: domknięcie R, bin 2–3)
#   dn4_job.sh extdp <próbka> <k> <n>      tory otoczenia D+ dla tabel etapu 1 (dn2_extract --species Dplus) → $DN4_OUT/parts_dp
#   dn4_job.sh train <blok>                dn4_train → $DN4_OUT/{sc,meta,models}
#   dn4_job.sh app <tryb> <próbka> <k> <n> tryb d150|d23: k-ta z n grup plików h5 → dn4_extract (100 % zdarzeń) → dn2_extract (Tracks) → dn4_apply
#                                          → $DN4_OUT/app/<tryb>/app_<próbka>_<k>.npz
#   dn4_job.sh ana                         analyze.py compute (G, rzeźbienie, SF z tagu D* dla wszystkich ziaren; + odniesienie: score etapu 2)
#                                          + dn4_ana.py r (domknięcie MC → R na danych) → dnet4_results.json, dnet4_ref2_results.json, dnet4_R.json, DONE
#   dn4_job.sh ana23                       dn4_ana.py d23 → dnet4_d23.json, DONE23
# DN4_SOFTNODE=1 (VARS soft): niepowodzenie → znacznik $DN4_OUT/fail/<węzeł> i kod 0 (DAG idzie dalej bez tego wyniku).
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
DN4=$(cd "$(dirname "$0")" && pwd); export DN4
for v in DN_SAMPLES DN_GROUPS DN_BLOCKS DN2_NPARTS DN2_FRAC DN2_BLOCKS DN4_BLOCKS DN4_SOFT DN4_EXTDP DN4_APP150 DN4_APP23; do
  [ -n "${!v:-}" ] && printf -v "$v" '%s' "${!v//+/ }"
done
source "$DN4/dn4_common.sh"
step=$1; shift
set -u
NODE="${step}_$(echo "$*" | tr ' ' '_')"; NODE=${NODE%_}
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/dn4_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */dn4_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }   # tylko własny katalog mktemp
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
put() { [ "${DRYRUN:-0}" = 1 ] || dn_put "$1" "$2"; }
echo "[$(date -u +%FT%TZ)] DNET4 $step $* host $(hostname) WD=$WD DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT DN4_OUT=$DN4_OUT soft=${DN4_SOFTNODE:-0}"
body() {
  local s k n mode tag nf i f id o d
  if [ "${DRYRUN:-0}" != 1 ]; then
    DEV=cpu
    case "$step:${1:-}:${DN4_GPU:-0}" in
      train:xgb*:*) dn4_env "$DN4_CPUS" ;;
      train:*:1)    # 19.11: trening sieci na GPU (dn4_gpu.sub); brak CUDA na węźle = błąd (bez cichego zejścia na CPU przy 4 rdzeniach)
        DN4_CPUS=4; DN2_USE_GPU=1; DN2_CPUS=4; export DN2_USE_GPU DN2_CPUS; dn2_env
        python3 -c 'import torch, sys; ok = torch.cuda.is_available(); print("GPU:", torch.cuda.get_device_name(0) if ok else None); sys.exit(0 if ok else 3)' || { echo "BŁĄD: brak CUDA na węźle GPU"; return 3; }
        DEV=cuda ;;
      app:*|extdp:*|pis:*) dn4_env "$DN4_APP_CPUS" ;;
      *) dn4_env "$DN4_CPUS" ;;
    esac
    mkdir -p "$DN4_OUT"/{parts_dp,sc,meta,models,app/d150,app/d23,fail} 2>/dev/null || true
  fi
  case $step in
    smoke)
      run python3 -u "$DN4/dn4_train.py" selftest --need-torch --dir "$WD/selftest"
      run python3 -u "$DN4/dn4_ana.py" selftest --nboot 20 --dir "$WD/anatest" ;;
    extdp)
      s=$1; k=$2; n=$3
      run python3 -u "$DN2/dn2_extract.py" --sample "$s" --data "$DN_OUT/data" --rawdir "$RAWV20/$s" --part "$k" --nparts "$n" --species Dplus \
          --frac-a 1.0 --out "$WD/dn2_${s}_${k}.npz"
      put "$WD/dn2_${s}_${k}.npz" "$DN4_OUT/parts_dp/dn2_${s}_${k}.npz" ;;
    train)
      export DN_DONE_DIR=$DN4_OUT
      run python3 -u "$DN4/dn4_train.py" train --data "$DN_OUT/data" --env "$DN2_OUT/parts" --env-dp "$DN4_OUT/parts_dp" --out "$WD/out" --block "$1" \
          --epochs "$DN4_EPOCHS" --cap-sig "$DN4_CAP_SIG" --cap-bkg "$DN4_CAP_BKG" --threads "$DN4_CPUS" --device "${DEV:-cpu}"
      if [ "${DRYRUN:-0}" != 1 ]; then
        for d in "$WD"/out/models/*/; do [ -d "$d" ] || continue; for f in "$d"*; do put "$f" "$DN4_OUT/models/$(basename "$d")/$(basename "$f")"; done; done
        for d in sc meta; do for f in "$WD/out/$d"/*; do [ -f "$f" ] && put "$f" "$DN4_OUT/$d/$(basename "$f")"; done; done   # meta na końcu = „gotowe”
      fi ;;
    app)
      mode=$1; s=$2; k=$3; n=$4
      # Plik po pliku (pamięć ~6 GB na plik danych 428770; test 19.11): wyjście app_<próbka>_<k>_<NNNNNN>.npz na każdy plik h5; gotowe pliki pomijane.
      mapfile -t ALL < <(dn4_files "$mode" "$s")
      nf=${#ALL[@]}; [ "$nf" -gt 0 ] || { echo "brak h5 dla $mode $s"; return 2; }
      FILES=(); for ((i = k; i < nf; i += n)); do FILES+=("${ALL[$i]}"); done
      [ ${#FILES[@]} -gt 0 ] || { echo "grupa $k/$n pusta ($nf plików) — nic do zrobienia"; return 0; }
      echo "pliki h5 ($mode $s $k/$n): ${#FILES[@]} z $nf"
      case $s in data_*) tag=data ;; fd_*) tag=fd ;; *) tag=mb ;; esac
      if [ "$mode" = d23 ]; then export DN4_PT=2,3; TAGS=$DN4_TAGS23; else export DN4_PT=3,30; TAGS=""; fi
      export DN4_MR_D0
      for f in "${FILES[@]}"; do
        id=$(basename "$f" | sed -E 's/.*_([0-9]{6})_none\.h5$/\1/')
        case $id in [0-9][0-9][0-9][0-9][0-9][0-9]) ;; *) echo "nieoczekiwana nazwa h5: $f"; return 2 ;; esac
        o="app_${s}_${k}_${id}.npz"
        if [ -s "$DN4_OUT/app/$mode/$o" ] && [ "${DN_REDO:-0}" != 1 ]; then echo "$o już jest — pomijam"; continue; fi
        rm -rf "$WD/h5" "$WD/data" "$WD/parts"; mkdir -p "$WD/h5/$s" "$WD/data" "$WD/parts"; ln -s "$f" "$WD/h5/$s/"
        run python3 -u "$DN4/dn4_extract.py" --tag "$tag" --frac-d0 1.0 --frac-dp 1.0 --lxy-min "$DN_LXY" --out "$WD/data/dn_${s}.npz" "$f"
        run python3 -u "$DN2/dn2_extract.py" --sample "$s" --data "$WD/data" --rawdir "$(dn4_rawdir "$mode")/$s" --h5dir "$WD/h5" --part 0 --nparts 1 \
            --species D0,Dplus --frac-a 1.0 --out "$WD/parts/dn2_${s}_0.npz"
        run python3 -u "$DN4/dn4_apply.py" --sample "$s" --data "$WD/data" --env "$WD/parts" --models "$DN4_OUT/models" --tags "$TAGS" \
            --threads "$DN4_APP_CPUS" --out "$WD/$o"
        put "$WD/$o" "$DN4_OUT/app/$mode/$o"; rm -f "$WD/$o"
      done ;;
    ana)
      fast=(); [ "$DN_FAST" = 1 ] && fast=(--fast)
      run python3 -u "$DN4/dn4_ana.py" r --app "$DN4_OUT/app/d150" --out "$WD/dnet4_R.json" --nboot "$DN4_NBOOT"
      put "$WD/dnet4_R.json" "$DN4_OUT/dnet4_R.json"
      run python3 -u "$DN1/analyze.py" compute --data "$DN_OUT/data" --sc "$DN_OUT,$DN4_OUT" --keep '^(D0|Dplus)\.tr\.H\+e04\.' --sf-all \
          --out "$WD/dnet4_results.json" --nproc "${DN_ANA_CPUS:-8}" --nboot "$DN_NBOOT" "${fast[@]}"
      put "$WD/dnet4_results.json" "$DN4_OUT/dnet4_results.json"
      run python3 -u "$DN1/analyze.py" compute --data "$DN_OUT/data" --sc "$DN_OUT,$DN2_OUT" --keep '^D0\.tr\.H\+e04\.mc\.d[0-9]+\.s[123]$' --sf-all \
          --out "$WD/dnet4_ref2_results.json" --nproc "${DN_ANA_CPUS:-8}" --nboot "$DN_NBOOT" "${fast[@]}" \
          || echo "UWAGA: analiza odniesienia (score etapu 2) nieudana — dnet4_ref2_results.json brak"
      [ -s "$WD/dnet4_ref2_results.json" ] && put "$WD/dnet4_ref2_results.json" "$DN4_OUT/dnet4_ref2_results.json"
      if [ "${DRYRUN:-0}" != 1 ]; then date -u +%FT%TZ > "$WD/DONE"; put "$WD/DONE" "$DN4_OUT/DONE"; fi ;;
    pis)   # diagnostyka wolnego pionu D* (tylko MC fd_nd, inferencja) → dnet4_pis.json
      run python3 -u "$DN4/dn4_pis.py" --data "$DN_OUT/data" --env "$DN2_OUT/parts" --models "$DN4_OUT/models" --threads "$DN4_APP_CPUS" --out "$WD/dnet4_pis.json"
      put "$WD/dnet4_pis.json" "$DN4_OUT/dnet4_pis.json" ;;
    ana23)
      run python3 -u "$DN4/dn4_ana.py" d23 --app "$DN4_OUT/app/d23" --r "$DN4_OUT/dnet4_R.json" --out "$WD/dnet4_d23.json" --nboot "$DN4_NBOOT"
      put "$WD/dnet4_d23.json" "$DN4_OUT/dnet4_d23.json"
      if [ "${DRYRUN:-0}" != 1 ]; then date -u +%FT%TZ > "$WD/DONE23"; put "$WD/DONE23" "$DN4_OUT/DONE23"; fi ;;
    *) echo "nieznany krok $step"; return 2 ;;
  esac
  return 0
}
set +e
( set -eo pipefail; body "$@" )        # podpowłoka POZA listą ||/&& — inaczej bash ignoruje set -e w jej wnętrzu
rc=$?
if [ $rc -ne 0 ] && [ "${DN4_SOFTNODE:-0}" = 1 ]; then
  echo "[$(date -u +%FT%TZ)] MIĘKKI BŁĄD DNET4 $NODE rc=$rc — znacznik fail/$NODE, DAG idzie dalej"
  echo "$(date -u +%FT%TZ) rc=$rc host=$(hostname)" > "$WD/fail_$NODE"; dn_put "$WD/fail_$NODE" "$DN4_OUT/fail/$NODE" || true
  exit 0
fi
[ $rc -eq 0 ] && echo "[$(date -u +%FT%TZ)] OK DNET4 $NODE"
exit $rc
