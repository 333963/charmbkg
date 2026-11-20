#!/bin/bash
# Węzeł DAG DNET etap 5a. Użycie:
#   dn5_job.sh fold <blok> <k>     trening JEDNEGO foldu → $DN5_OUT/models/<tag>/fold<k>.pt + folds/<tag>/fold<k>.json (json na końcu = „gotowe”)
#   dn5_job.sh merge <blok>        modele foldów z EOS → score OOF wszystkich tabel etapu 1 → sc/<tag>.npz, models/<tag>/meta.json, meta/<tag>.json
#   dn5_job.sh gtrain <blok>       GPU: 5 foldów po kolei + łączenie w jednym jobie (te same pliki wyjściowe co fold × K + merge)
#   dn5_job.sh ana                 analyze.py compute (G, rzeźbienie, SF z tagu D*; format etapów 1–4) → dnet5_results.json;
#                                  dn5_ana.py (pary sparowane, kaskada H × E) → dnet5_cmp.json; DONE
#   dn5_job.sh eqbench             test równoważności stara/nowa ścieżka (dane syntetyczne) + czas epoki vs liczba wątków → test/{equiv,bench}.json
#   dn5_job.sh gpu                 sonda GPU (poza DAG): fold 0 bloków $DN5_GPU_BLOCKS na GPU → gpu/probe.json (czas w kolejce, czas epoki)
# DN5_SOFTNODE=1: niepowodzenie → znacznik $DN5_OUT/fail/<węzeł> i kod 0. Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi.
set -eo pipefail
DN5=$(cd "$(dirname "$0")" && pwd); export DN5
for v in DN_SAMPLES DN_GROUPS DN_BLOCKS DN2_NPARTS DN2_FRAC DN2_BLOCKS DN4_BLOCKS DN4_SOFT DN4_EXTDP DN4_APP150 DN4_APP23 DN5_BLOCKS DN5_GPU_BLOCKS; do
  [ -n "${!v:-}" ] && printf -v "$v" '%s' "${!v//+/ }"
done
[ -n "${DN5_KEEP_B64:-}" ] && DN5_KEEP=$(echo "$DN5_KEEP_B64" | base64 -d) && export DN5_KEEP
source "$DN5/dn5_common.sh"
step=$1; shift
set -u
NODE="${step}_$(echo "$*" | tr ' ' '_')"; NODE=${NODE%_}
T_START=$(date +%s)
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/dn5_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */dn5_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }   # tylko własny katalog mktemp
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
put() { [ "${DRYRUN:-0}" = 1 ] || dn_put "$1" "$2"; }
echo "[$(date -u +%FT%TZ)] DNET5 $step $* host $(hostname) WD=$WD DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT DN4_OUT=$DN4_OUT DN5_OUT=$DN5_OUT soft=${DN5_SOFTNODE:-0}"
COM=(--data "$DN_OUT/data" --env "$DN2_OUT/parts" --env-dp "$DN4_OUT/parts_dp" --hsc "$DN_OUT" --epochs "$DN5_EPOCHS" --cap-sig "$DN5_CAP_SIG" --cap-bkg "$DN5_CAP_BKG")
body() {
  local b k tag f t0 tsub
  if [ "${DRYRUN:-0}" != 1 ]; then
    case $step in
      gpu|gtrain) DN2_USE_GPU=1; DN2_CPUS=$DN5_CPUS; export DN2_USE_GPU DN2_CPUS; dn2_env ;;
      eqbench) dn4_env 8 ;;
      ana) dn4_env "$DN5_ANA_CPUS" ;;
      *) dn4_env "$DN5_CPUS" ;;
    esac
    mkdir -p "$DN5_OUT"/{models,folds,sc,meta,test,gpu,fail} 2>/dev/null || true
  fi
  case $step in
    fold)
      b=$1; k=$2; tag=$(dn5_tag "$b"); [ -n "$tag" ] || { echo "nieznany blok $b"; return 2; }
      if [ -s "$DN5_OUT/models/$tag/fold$k.pt" ] && [ -s "$DN5_OUT/folds/$tag/fold$k.json" ] && [ "${DN_REDO:-0}" != 1 ]; then echo "$tag fold $k już jest — pomijam"; return 0; fi
      run python3 -u "$DN5/dn5_train.py" fold "${COM[@]}" --out "$WD/out" --block "$b" --fold "$k" --threads "$DN5_CPUS" --device cpu
      put "$WD/out/models/$tag/fold$k.pt" "$DN5_OUT/models/$tag/fold$k.pt"
      put "$WD/out/folds/$tag/fold$k.json" "$DN5_OUT/folds/$tag/fold$k.json" ;;
    merge)
      b=$1; tag=$(dn5_tag "$b"); [ -n "$tag" ] || { echo "nieznany blok $b"; return 2; }
      if [ -s "$DN5_OUT/meta/$tag.json" ] && [ -s "$DN5_OUT/sc/$tag.npz" ] && [ "${DN_REDO:-0}" != 1 ]; then echo "$tag już złączony — pomijam"; return 0; fi
      for ((k = 0; k < DN_FOLDS; k++)); do [ -s "$DN5_OUT/models/$tag/fold$k.pt" ] && [ -s "$DN5_OUT/folds/$tag/fold$k.json" ] || { echo "brak foldu $k dla $tag"; return 3; }; done
      run python3 -u "$DN5/dn5_train.py" merge "${COM[@]}" --out "$WD/out" --block "$b" --models "$DN5_OUT/models" --folds-dir "$DN5_OUT/folds" --threads "$DN5_CPUS"
      put "$WD/out/sc/$tag.npz" "$DN5_OUT/sc/$tag.npz"
      put "$WD/out/models/$tag/meta.json" "$DN5_OUT/models/$tag/meta.json"
      put "$WD/out/meta/$tag.json" "$DN5_OUT/meta/$tag.json" ;;                      # meta na końcu = „gotowe”
    gtrain)   # 20.11-styl: cały wariant na GPU w jednym jobie — 5 foldów po kolei (każdy osobnym procesem → ta sama inicjalizacja co job foldu) + łączenie
      b=$1; tag=$(dn5_tag "$b"); [ -n "$tag" ] || { echo "nieznany blok $b"; return 2; }
      if [ -s "$DN5_OUT/meta/$tag.json" ] && [ -s "$DN5_OUT/sc/$tag.npz" ] && [ "${DN_REDO:-0}" != 1 ]; then echo "$tag już złączony — pomijam"; return 0; fi
      python3 -c 'import torch; t = torch.randn(64, 8, device="cuda"); print("CUDA OK", torch.cuda.get_device_name(0), float((t @ t.T).sum()))' \
        || { echo "BŁĄD: GPU nie działa na $(hostname) — węzeł gtrain nie liczy na CPU"; return 5; }
      mkdir -p "$WD/out/models/$tag" "$WD/out/folds/$tag"
      for ((k = 0; k < DN_FOLDS; k++)); do
        # fold policzony wcześniej na GPU (ponowienie węzła) → wziąć z EOS; foldy policzone na CPU są liczone od nowa (jeden wariant = jedno urządzenie)
        if [ -s "$DN5_OUT/models/$tag/fold$k.pt" ] && grep -q '"device": "cuda"' "$DN5_OUT/folds/$tag/fold$k.json" 2>/dev/null && [ "${DN_REDO:-0}" != 1 ]; then
          cp "$DN5_OUT/models/$tag/fold$k.pt" "$WD/out/models/$tag/"; cp "$DN5_OUT/folds/$tag/fold$k.json" "$WD/out/folds/$tag/"; echo "$tag fold $k (GPU) już jest"; continue
        fi
        run python3 -u "$DN5/dn5_train.py" fold "${COM[@]}" --out "$WD/out" --block "$b" --fold "$k" --threads "$DN5_CPUS" --device cuda
        grep -q '"device": "cuda"' "$WD/out/folds/$tag/fold$k.json" || { echo "BŁĄD: fold $k policzony nie na GPU"; return 5; }
        put "$WD/out/models/$tag/fold$k.pt" "$DN5_OUT/models/$tag/fold$k.pt"
        put "$WD/out/folds/$tag/fold$k.json" "$DN5_OUT/folds/$tag/fold$k.json"
      done
      run python3 -u "$DN5/dn5_train.py" merge "${COM[@]}" --out "$WD/out" --block "$b" --models "$WD/out/models" --folds-dir "$WD/out/folds" --threads "$DN5_CPUS" --device cuda
      put "$WD/out/sc/$tag.npz" "$DN5_OUT/sc/$tag.npz"
      put "$WD/out/models/$tag/meta.json" "$DN5_OUT/models/$tag/meta.json"
      put "$WD/out/meta/$tag.json" "$DN5_OUT/meta/$tag.json" ;;
    ana)
      fast=(); [ "$DN_FAST" = 1 ] && fast=(--fast)
      sc4=(); [ -d "$DN4_OUT/sc" ] && sc4=(--sc4 "$DN4_OUT")
      run python3 -u "$DN5/dn5_ana.py" --data "$DN_OUT/data" --sc "$DN_OUT,$DN2_OUT,$DN5_OUT" "${sc4[@]}" --out "$WD/dnet5_cmp.json" --nboot "$DN_NBOOT"
      put "$WD/dnet5_cmp.json" "$DN5_OUT/dnet5_cmp.json"
      run python3 -u "$DN1/analyze.py" compute --data "$DN_OUT/data" --sc "$DN_OUT,$DN2_OUT,$DN5_OUT" --keep "$DN5_KEEP" \
          --out "$WD/dnet5_results.json" --nproc "$DN5_ANA_CPUS" --nboot "$DN_NBOOT" "${fast[@]}"
      put "$WD/dnet5_results.json" "$DN5_OUT/dnet5_results.json"
      if [ "${DRYRUN:-0}" != 1 ]; then date -u +%FT%TZ > "$WD/DONE"; put "$WD/DONE" "$DN5_OUT/DONE"; fi ;;
    eqbench)
      run python3 -u "$DN5/dn5_train.py" equiv --dir "$WD/eq" --threads 4 --json "$WD/equiv.json" || echo "UWAGA: equiv rc=$?"
      [ -s "$WD/equiv.json" ] && put "$WD/equiv.json" "$DN5_OUT/test/equiv.json"
      run python3 -u "$DN5/dn5_train.py" bench --threads-list 1,2,4,8 --json "$WD/bench.json"
      put "$WD/bench.json" "$DN5_OUT/test/bench.json" ;;
    gpu)
      tsub=$(cat "$DN5_OUT/gpu/submitted" 2>/dev/null || echo 0)
      nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | tee "$WD/gpuname.txt" || true
      for b in $DN5_GPU_BLOCKS; do
        run python3 -u "$DN5/dn5_train.py" fold "${COM[@]}" --out "$WD/gpu" --block "$b" --fold 0 --threads "$DN5_CPUS" --device auto || echo "UWAGA: sonda $b rc=$?"
      done
      python3 - "$WD" "$tsub" "$T_START" <<'PY'
import glob, json, os, sys
wd, tsub, t0 = sys.argv[1], int(sys.argv[2] or 0), int(sys.argv[3])
r = {"t_submit": tsub, "t_start": t0, "queue_wait_s": (t0 - tsub) if tsub else None, "host": os.uname().nodename,
     "gpu": open(os.path.join(wd, "gpuname.txt")).read().strip() if os.path.isfile(os.path.join(wd, "gpuname.txt")) else None, "folds": {}}
for f in sorted(glob.glob(os.path.join(wd, "gpu", "folds", "*", "fold0.json"))):
    j = json.load(open(f))
    r["folds"][os.path.basename(os.path.dirname(f))] = {k: j.get(k) for k in ("device", "device_name", "sec", "sec_epoch", "sec_load", "sec_total", "epochs",
                                                                                    "n_pool", "n_train", "best_epoch", "threads", "nparams")}
json.dump(r, open(os.path.join(wd, "probe.json"), "w"), indent=1)
print("SONDA GPU:", json.dumps(r))
PY
      put "$WD/probe.json" "$DN5_OUT/gpu/probe.json" ;;
    *) echo "nieznany krok $step"; return 2 ;;
  esac
  return 0
}
set +e
( set -eo pipefail; body "$@" )        # podpowłoka POZA listą ||/&& — inaczej bash ignoruje set -e w jej wnętrzu
rc=$?
if [ $rc -ne 0 ] && [ "${DN5_SOFTNODE:-0}" = 1 ]; then
  echo "[$(date -u +%FT%TZ)] MIĘKKI BŁĄD DNET5 $NODE rc=$rc — znacznik fail/$NODE, DAG idzie dalej"
  echo "$(date -u +%FT%TZ) rc=$rc host=$(hostname)" > "$WD/fail_$NODE"; dn_put "$WD/fail_$NODE" "$DN5_OUT/fail/$NODE" || true
  exit 0
fi
[ $rc -eq 0 ] && echo "[$(date -u +%FT%TZ)] OK DNET5 $NODE ($(( $(date +%s) - T_START )) s)"
exit $rc
